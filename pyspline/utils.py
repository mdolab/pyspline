"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# ===========================================================================
# External Python modules
# ===========================================================================
import numpy
from scipy import sparse

# ===========================================================================
# Custom Python modules
# ===========================================================================
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
        1D aray of data to write to file
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
    data : 2D numpy array of sive (nx, ny, ndim)
        2D aray of data to write to file
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
    data : 3D numpy array of sive (nx, ny, nz, ndim)
        3D aray of data to write to file
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
    rank = numpy.ndim(inputVal)
    if rank != dataRank:
        raise Error("'%s' must have rank %d. Input was of rank %d." % (inputName, dataRank, rank))

    # Check the data type
    inputVal = numpy.array(inputVal)
    tmp = inputVal.astype(dataType)

    # Check if the values are the same:
    diff = (tmp - inputVal).flatten()
    if numpy.dot(diff, diff) > 10 * numpy.finfo(1.0).eps:
        raise Error("'%s' could not be safely cast to required type" "without losing information" % inputName)

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
