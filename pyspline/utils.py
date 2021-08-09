# External modules
import numpy as np
from scipy import sparse
from . import libspline


class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a explicitly raised exception.
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| pySpline Error: "
        i = 17
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)
        Exception.__init__()


def writeTecplot1D(handle, name, data):
    """A Generic function to write a 1D data zone to a tecplot file.

    Parameters
    ----------
    handle : file handle
        Open file handle
    name : str
        Name of the zone to use
    data : array of size (N, ndim)
        1D array of data to write to file
    """
    nx = data.shape[0]
    ndim = data.shape[1]
    handle.write('Zone T="%s" I=%d\n' % (name, nx))
    handle.write("DATAPACKING=POINT\n")
    for i in range(nx):
        for idim in range(ndim):
            handle.write("%f " % (data[i, idim]))
        handle.write("\n")


def writeTecplot2D(handle, name, data):
    """A Generic function to write a 2D data zone to a tecplot file.

    Parameters
    ----------
    handle : file handle
        Open file handle
    name : str
        Name of the zone to use
    data : 2D np array of size (nx, ny, ndim)
        2D array of data to write to file
    """
    nx = data.shape[0]
    ny = data.shape[1]
    ndim = data.shape[2]
    handle.write('Zone T="%s" I=%d J=%d\n' % (name, nx, ny))
    handle.write("DATAPACKING=POINT\n")
    for j in range(ny):
        for i in range(nx):
            for idim in range(ndim):
                handle.write("%20.16g " % (data[i, j, idim]))
            handle.write("\n")


def writeTecplot3D(handle, name, data):
    """A Generic function to write a 3D data zone to a tecplot file.

    Parameters
    ----------
    handle : file handle
        Open file handle
    name : str
        Name of the zone to use
    data : 3D np array of size (nx, ny, nz, ndim)
        3D array of data to write to file
    """
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    ndim = data.shape[3]
    handle.write('Zone T="%s" I=%d J=%d K=%d\n' % (name, nx, ny, nz))
    handle.write("DATAPACKING=POINT\n")
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                for idim in range(ndim):
                    handle.write("%f " % (data[i, j, k, idim]))
                handle.write("\n")


def _writeHeader(f, ndim):
    """Write tecplot zone header depending on spatial dimension"""
    if ndim == 1:
        f.write('VARIABLES = "X"\n')
    elif ndim == 2:
        f.write('VARIABLES = "X", "Y"\n')
    else:
        f.write('VARIABLES = "X", "Y", "Z"\n')


def openTecplot(fileName, ndim):
    """A Generic function to open a Tecplot file to write spatial data.

    Parameters
    ----------
    fileName : str
        Tecplot filename. Should have a .dat extension.
    ndim : int
        Number of spatial dimensions. Must be 1, 2 or 3.

    Returns
    -------
    f : file handle
        Open file handle
    """
    f = open(fileName, "w")
    _writeHeader(f, ndim)

    return f


def closeTecplot(f):
    """Close Tecplot file opened with openTecplot()"""
    f.close()


def _assembleMatrix(data, indices, indptr, shape):
    """
    Generic assemble matrix function to create a CSR matrix

    Parameters
    ----------
    data : array
        Data values for matrix
    indices : int array
        CSR type indices
    indptr : int array
        Row pointer
    shape : tuple-like
        Actual shape of matrix

    Returns
    -------
    M : scipy csr sparse matrix
        The assembled matrix
    """
    M = sparse.csr_matrix((data, indices, indptr), shape)

    return M


def checkInput(inputVal, inputName, dataType, dataRank, dataShape=None):
    """This is a generic function to check the data type and sizes of
    inputs in functions where the user must supply proper
    values. Since Python does not do type checking on Inputs, this is
    required

    Parameters
    ----------
    input : int, float, or complex
        The input argument to check
    inputName : str
        The name of the variable, so it can be identified in an Error message
    dataType : str
        Numpy string dtype code
    dataRank : int
        Desired rank. 0 for scalar, 1 for vector 2 for matrix etc
    dataShape : tuple
        The required shape of the data.
        Scalar is denoted by ()
        Vector is denoted by (n, ) where n is the shape
        Matrix is denoted by (n, m) where n and m are rows/columns

    Returns
    -------
    output : various
        The input transformed to the correct type and guaranteed to
        the desired size and shape.
    """

    # Check the data rank:
    rank = np.ndim(inputVal)
    if rank != dataRank:
        raise Error("'%s' must have rank %d. Input was of rank %d." % (inputName, dataRank, rank))

    # Check the data type
    inputVal = np.array(inputVal)
    tmp = inputVal.astype(dataType)

    # Check if the values are the same:
    diff = (tmp - inputVal).flatten()
    if np.dot(diff, diff) > 10 * np.finfo(1.0).eps:
        raise Error("'%s' could not be safely cast to required type without losing information" % inputName)

    # Finally check that the shape is correct:
    if dataShape is not None:
        if tmp.shape != tuple(dataShape):
            raise Error(
                "'%s' was not the correct shape. Input was shape "
                "%s while requested dataShape was '%s'." % (inputName, repr(tmp.shape), repr(dataShape))
            )
    if dataRank == 0:
        return tmp.squeeze()
    else:
        return tmp


def trilinearVolume(*args):
    """This is a short-cut function to create a trilinear b-spline
    volume. It can be created with ``x`` **OR** with ``xmin`` and
    ``xmax``.

    Parameters
    ----------
    x : array of size (2, 2, 2, 3)
        Coordinates of the corners of the box.

    xmin : array of size (3)
        The extreme lower corner of the box
    xmax : array of size (3)
        The extreme upper corner of the box. In this case, by
        construction, the box will be coordinate axis aligned.
    """
    # Local modules
    from .pyVolume import Volume

    tu = [0, 0, 1, 1]
    tv = [0, 0, 1, 1]
    tw = [0, 0, 1, 1]
    ku = 2
    kv = 2
    kw = 2

    if len(args) == 1:
        return Volume(coef=args[0], tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    elif len(args) == 2:
        xmin = np.array(args[0]).astype("d")
        xmax = np.array(args[1]).astype("d")

        xLow = xmin[0]
        xHigh = xmax[0]
        yLow = xmin[1]
        yHigh = xmax[1]
        zLow = xmin[2]
        zHigh = xmax[2]

        coef = np.zeros((2, 2, 2, 3))
        coef[0, 0, 0, :] = [xLow, yLow, zLow]
        coef[1, 0, 0, :] = [xHigh, yLow, zLow]
        coef[0, 1, 0, :] = [xLow, yHigh, zLow]
        coef[1, 1, 0, :] = [xHigh, yHigh, zLow]
        coef[0, 0, 1, :] = [xLow, yLow, zHigh]
        coef[1, 0, 1, :] = [xHigh, yLow, zHigh]
        coef[0, 1, 1, :] = [xLow, yHigh, zHigh]
        coef[1, 1, 1, :] = [xHigh, yHigh, zHigh]

        return Volume(coef=coef, tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    else:
        raise Error("An unknown number of arguments was passed to trilinear  Volume")


def bilinearSurface(*args):
    """This is short-cut function to create a bilinear surface.

    Args can contain:

    1.  ``x`` array of size(4, 3).  The four corners of the array
        arranged in the coordinate direction orientation::

          2          3
          +----------+
          |          |
          |          |
          |          |
          +----------+
          0          1

    2. ``p1``, ``p2``, ``p3``, ``p4``. Individual corner points in CCW Ordering::

          3          2
          +----------+
          |          |
          |          |
          |          |
          +----------+
          0          1
    """
    # Local modules
    from .pySurface import Surface

    if len(args) == 1:
        # One argument passed in ... assume its X
        if len(args[0]) != 4:
            raise Error("A single argument passed to bilinear surface must contain 4 points and be of size (4, 3)")
        coef = np.zeros((2, 2, 3))
        coef[0, 0] = args[0][0]
        coef[1, 0] = args[0][1]
        coef[0, 1] = args[0][2]
        coef[1, 1] = args[0][3]
        return Surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)
    else:
        # Assume 4 arguments
        coef = np.zeros([2, 2, 3])
        coef[0, 0] = args[0]
        coef[1, 0] = args[1]
        coef[0, 1] = args[3]
        coef[1, 1] = args[2]
        return Surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)


def line(*args, **kwargs):
    """This is a short cut function to create a line as a pySpline Curve object.

    Args can contain: (pick one)

    1. ``X`` array of size(2, ndim). The two end points
    2. ``x1``, ``x2`` The two end points (each of size ndim)
    3. ``x```, dir=direction. A point and a displacement vector
    4. ``x1``, dir=direction, length=length. As 3. but with a specific length
    """
    # Local modules
    from .pyCurve import Curve

    if len(args) == 2:
        # Its a two-point type
        return Curve(coef=[args[0], args[1]], k=2, t=[0, 0, 1, 1])
    elif len(args) == 1:
        if len(args[0]) == 2:  # its X
            return Curve(coef=args[0], k=2, t=[0, 0, 1, 1])
        elif "dir" in kwargs:
            # We have point and direction
            if "length" in kwargs:
                x2 = args[0] + kwargs["dir"] / np.linalg.norm(kwargs["dir"]) * kwargs["length"]
            else:
                x2 = args[0] + kwargs["dir"]

            return Curve(coef=[args[0], x2], k=2, t=[0, 0, 1, 1])
        else:
            Error("Error: dir must be specified if only 1 argument is given")


def plane_line(pts_T, vec_T, p0, v1, v2):
    """Check a plane against multiple lines
    Python wrapper for libspline function

    Args:
        pts_T ([real]): [description]
        vec_T ([real]): [description]
        p0 ([real]): [description]
        v1 ([real]): [description]
        v2 ([real]): [description]

    Returns:
        [int]: [description]
        [real]: [description]
    """
    return libspline.plane_line(pts_T, vec_T, p0, v1, v2)


def tfi2d(e0, e1, e2, e3):
    """[summary]

    Args:
        e0 ([type]): [description]
        e1 ([type]): [description]
        e2 ([type]): [description]
        e3 ([type]): [description]

    Returns:
        [type]: [description]
    """
    return libspline.tfi2d(e0, e1, e2, e3)


def line_plane(pt, upVec, p0, v1, v2):
    """[summary]

    Args:
        pt ([type]): [description]
        upVec ([type]): [description]
        p0 ([type]): [description]
        v1 ([type]): [description]
        v2 ([type]): [description]

    Returns:
        [type]: [description]
        [type]: [description]
        [type]: [description]
    """
    return libspline.line_plane(pt, upVec, p0, v1, v2)
