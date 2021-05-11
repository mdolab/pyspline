# Standard Python modules
import warnings

# External modules
import numpy
from scipy.sparse import linalg

# Local modules
from . import libspline
from .pyCurve import Curve
from .pySurface import Surface
from .utils import Error, _assembleMatrix, checkInput, closeTecplot, openTecplot, writeTecplot3D


class Volume(object):
    """
    Create an instance of a b-spline surface. There are two
    ways to initialize the class

    * **Creation**: Create an instance of the Volume class directly by
      supplying the required information: kwargs MUST contain the
      folloiwng information: ``ku, kv, kw, tu, tv, tw, coef``.

    * **LMS/Interpolation**: Create an instance of the Volume class by
      using an interpolating spline to given data points or a LMS
      spline. The following keyword argument information is required:

      1. ``ku`` and ``kv`` and ``kw`` Spline Orders

      2. ``X`` real arry size (Nu, Nv, Nw, nDim) of data to fit. **OR**
          1. ``x`` (3D) and ``y`` (3D)  ``z`` (3D) 3D volume interpolation/fitting

      3. ``u``, ``v``, ``w`` real array of size (Nu, Nv, Nw). Optional

      4. ``nCtlu``, ``nCtlv``, ``nCtlw``.  Specify number of control points. Only for
         LMS fitting.

    Parameters
    ----------
    ku : int
       Spline order in u
    kv : int
       Spline order in v
    kw : int
       Spline order in w
    nCtlu : int
       Number of control points in u
    nCtlv : int
       Number of control points in u
    nCtlw : int
       Number of control points in w
    coef : array, size (nCtlu, nCtl, nDim)
       b-spline coefficient array.
    tu : array, size(nCtlu + ku)
       knot array in u
    tv : array, size(nCtlv + kv)
       knot array in v
    tw : array, size(nCtlw + kw)
       knot array in w
    X : array, size (Nu, Nv, Nw, ndim)
       Full data array to fit
    x : array, size (Nu, Nv)
       Just x data to fit/interpolate
    y : array, size (Nu, Nv, Nw)
       Just y data to fit/interpolate
    z : array, size (Nu, Nv, Nw)
       Just z data to fit/interpolate
    u : array, size (Nu, Nv, Nw)
       Explict u parameters to use. Optional.
    v : array, size (Nu, Nv, Nw)
       Explict v parameters to use. Optional.
    w : array, size (Nu, Nv, Nw)
       Explict w parameters to use. Optional.
    nIter : int
       Number of Hoscheks parater corrections to run
    recompute : bool
       Specifies whether the actual fitting is completed.

    Notes
    -----
    The orientation of the nodes, edges and faces for the volumes is
    given below::

               NODES      |           EDGES         |           FACES
           6           7|             5          |
           #-----------#|       #------------#  |          #-----------#
          /           / |      /|           /|  |         /|          /|
         /           /  |     / |          / |  |        / |         / |
        /           /   |   6/  |        7/  |  |       /  |   1    /  |
       /           /    |   /   |10      /   |11|      /   |    ---------- 5
      /           /     |  /    |    4  /    |  |     /    |      /    |(back)
     #-----------#      | #------------#     |  |    #-----------#     |
     4           5      | |     |      |     |  |    |     |     |     |
                        | |     |      |     |  |    |     |     |     | <-3
           2           3| |     |   1  |     |  |2-> |     |     |     |
           #-----------#| |     #------|-----#  |    |     #-----|-----#
          /           / | |8   /       |9   /   |4 ----------    |    /
         /           /  | |   /        |   /    |    |   /       |   /
        /           /   | |  /2        |  /3    |    |  /      0 |  /
       /           /    | | /          | /      |    | /         | /
      /           /     | |/           |/       |    |/          |/
     #-----------#      | #------------#        |    #-----------#
     0           1      |         0             |
    """

    def __init__(self, recompute=True, **kwargs):
        self.faceSurfaces = [None, None, None, None, None, None]
        self.edgeCurves = [None, None, None, None, None, None, None, None, None, None, None, None]
        self.data = None
        self.udata = None
        self.vdata = None
        self.wdata = None
        if (
            "ku" in kwargs
            and "kv" in kwargs
            and "kw" in kwargs
            and "tu" in kwargs
            and "tv" in kwargs
            and "tw" in kwargs
            and "coef" in kwargs
        ):
            self.X = None
            self.u = None
            self.v = None
            self.w = None
            self.U = None
            self.V = None
            self.W = None
            self.interp = False

            self.ku = checkInput(kwargs["ku"], "ku", int, 0)
            self.kv = checkInput(kwargs["kv"], "kv", int, 0)
            self.kw = checkInput(kwargs["kw"], "kw", int, 0)
            self.coef = checkInput(kwargs["coef"], "coef", float, 4)
            self.nCtlu = self.coef.shape[0]
            self.nCtlv = self.coef.shape[1]
            self.nCtlw = self.coef.shape[2]
            self.nDim = self.coef.shape[3]
            self.tu = checkInput(kwargs["tu"], "tu", float, 1, (self.nCtlu + self.ku,))
            self.tv = checkInput(kwargs["tv"], "tv", float, 1, (self.nCtlv + self.kv,))
            self.tw = checkInput(kwargs["tw"], "tw", float, 1, (self.nCtlw + self.kw,))

            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.wmin = self.tw[0]
            self.wmax = self.tw[-1]
            self.origData = False
            self.setFaceSurfaces()
            self.setEdgeCurves()

        else:  # We have LMS/Interpolate
            # Do some checking on the number of control points
            assert (
                "ku" in kwargs
                and "kv" in kwargs
                and "kw" in kwargs
                and (
                    "X" in kwargs
                    or "x" in kwargs
                    or ("x" in kwargs and "y" in kwargs)
                    or ("x" in kwargs and "y" in kwargs and "z" in kwargs)
                )
            ), "Error: ku, kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate"

            if "X" in kwargs:
                self.X = numpy.array(kwargs["X"])
                if len(self.X.shape) == 1:
                    self.nDim = 1
                else:
                    self.nDim = self.X.shape[3]
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                x = checkInput(kwargs["x"], "x", float, 3)
                y = checkInput(kwargs["y"], "y", float, 3)
                z = checkInput(kwargs["z"], "z", float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[2], 3))
                self.X[:, :, :, 0] = x
                self.X[:, :, :, 1] = y
                self.X[:, :, :, 2] = z
                self.nDim = 3
            elif "x" in kwargs and "y" in kwargs:
                x = checkInput(kwargs["x"], "x", float, 3)
                y = checkInput(kwargs["y"], "y", float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[3], 3))
                self.X[:, :, :, 0] = x
                self.X[:, :, :, 1] = y
                self.nDim = 2
            elif "x" in kwargs:
                x = checkInput(kwargs["x"], "x", float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[3], 3))
                self.X[:, :, :, 0] = kwargs["x"]
                self.nDim = 1
            # enf if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.Nw = self.X.shape[2]
            self.ku = checkInput(kwargs["ku"], "ku", int, 0)
            self.kv = checkInput(kwargs["kv"], "kv", int, 0)
            self.kw = checkInput(kwargs["kw"], "kw", int, 0)

            if "nCtlu" in kwargs and "nCtlv" in kwargs and "nCtlw" in kwargs:
                self.nCtlu = checkInput(kwargs["nCtlu"], "nCtlu", int, 0)
                self.nCtlv = checkInput(kwargs["nCtlv"], "nCtlv", int, 0)
                self.nCtlw = checkInput(kwargs["nCtlw"], "nCtlw", int, 0)

                self.interp = False
            else:
                self.nCtlu = self.Nu
                self.nCtlv = self.Nv
                self.nCtlw = self.Nw
                self.interp = True

            self.origData = True

            # Sanity Check on Inputs
            if self.nCtlu >= self.Nu:
                self.nCtlu = self.Nu
            if self.nCtlv >= self.Nv:
                self.nCtlv = self.Nv
            if self.nCtlw >= self.Nw:
                self.nCtlw = self.Nw

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku:
                self.ku = self.Nu
            if self.Nv < self.kv:
                self.kv = self.Nv
            if self.Nw < self.kw:
                self.kw = self.Nw
            if self.nCtlu < self.ku:
                self.ku = self.nCtlu
            if self.nCtlv < self.kv:
                self.kv = self.nCtlv
            if self.nCtlw < self.kw:
                self.kw = self.nCtlw

            if "nIter" in kwargs:
                self.nIter = kwargs["nIter"]
            else:
                self.nIter = 1

            if "u" in kwargs and "v" in kwargs and "w" in kwargs:
                self.u = checkInput(kwargs["u"], "u", float, 1)
                self.v = checkInput(kwargs["v"], "v", float, 1)
                self.w = checkInput(kwargs["w"], "w", float, 1)
            else:
                if self.nDim == 3:
                    self.calcParameterization()
                else:
                    Error(
                        "Automatic parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise."
                    )

            self.umin = 0
            self.umax = 1
            self.vmin = 0
            self.vmax = 1
            self.wmin = 0
            self.wmax = 1
            self.calcKnots()
            self.setCoefSize()

            if recompute:
                self.recompute()
        # end if (Interpolation type)

    def recompute(self):
        """Recompute the volume if any driving data has been modified"""

        self.setCoefSize()

        vals, rowPtr, colInd = libspline.volume_jacobian_wrap(
            self.U,
            self.V,
            self.W,
            self.tu,
            self.tv,
            self.tw,
            self.ku,
            self.kv,
            self.kw,
            self.nCtlu,
            self.nCtlv,
            self.nCtlw,
        )

        N = _assembleMatrix(vals, colInd, rowPtr, (self.Nu * self.Nv * self.Nw, self.nCtlu * self.nCtlv * self.nCtlw))
        if self.interp:
            # Factorize once for efficiency
            solve = linalg.dsolve.factorized(N)
            for idim in range(self.nDim):
                self.coef[:, :, :, idim] = solve(self.X[:, :, :, idim].flatten()).reshape(
                    [self.nCtlu, self.nCtlv, self.nCtlw]
                )
        else:
            solve = linalg.dsolve.factorized(N.transpose() * N)
            for idim in range(self.nDim):
                rhs = N.transpose() * self.X[:, :, :, idim].flatten()
                self.coef[:, :, :, idim] = solve(rhs).reshape([self.nCtlu, self.nCtlv, self.nCtlw])

        self.setFaceSurfaces()
        self.setEdgeCurves()

    def setCoefSize(self):
        self.coef = numpy.zeros((self.nCtlu, self.nCtlv, self.nCtlw, self.nDim))

    def calcParameterization(self):
        """Compute distance based parametrization. Use the fortran
        function for this"""
        S, u, v, w = libspline.para3d(self.X.T)
        S = S.T
        self.u = u
        self.v = v
        self.w = w

        self.U = numpy.asarray(S[:, :, :, 0], order="c")
        self.V = numpy.asarray(S[:, :, :, 1], order="c")
        self.W = numpy.asarray(S[:, :, :, 2], order="c")

        return

    def calcKnots(self):
        """Determine the knots depending on if it is inerpolated or
        an LMS fit"""
        if self.interp:
            self.tu = libspline.knots_interp(self.u, numpy.array([], "d"), self.ku)
            self.tv = libspline.knots_interp(self.v, numpy.array([], "d"), self.kv)
            self.tw = libspline.knots_interp(self.w, numpy.array([], "d"), self.kw)
        else:
            self.tu = libspline.knots_lms(self.u, self.nCtlu, self.ku)
            self.tv = libspline.knots_lms(self.v, self.nCtlv, self.kv)
            self.tw = libspline.knots_lms(self.w, self.nCtlw, self.kw)

    def getValueCorner(self, corner):
        """Get the value of the volume spline on the corner.

        Parameters
        ----------
        corner : int
            Index of corner, 0<=corner<=7

        Returns
        -------
        value : float
            Volume spline evaluation at corner.
        """

        if corner not in range(0, 8):
            raise Error("Corner must be in range 0..7 inclusive8")

        if corner == 0:
            val = self.getValue(self.umin, self.vmin, self.wmin)
        elif corner == 1:
            val = self.getValue(self.umax, self.vmin, self.wmin)
        elif corner == 2:
            val = self.getValue(self.umin, self.vmax, self.wmin)
        elif corner == 3:
            val = self.getValue(self.umax, self.vmax, self.wmin)
        elif corner == 4:
            val = self.getValue(self.umin, self.vmin, self.wmax)
        elif corner == 5:
            val = self.getValue(self.umax, self.vmin, self.wmax)
        elif corner == 6:
            val = self.getValue(self.umin, self.vmax, self.wmax)
        elif corner == 7:
            val = self.getValue(self.umax, self.vmax, self.wmax)

        return val

    def getOrigValueCorner(self, corner):
        """Get the value of the original spline data on the corner if
        it exists

        Parameters
        ----------
        corner : int
            Index of corner, 0<=corner<=7

        Returns
        -------
        value : float
            Original data on corner.
        """
        if corner not in range(0, 8):
            raise Error("Corner must be in range 0..7 inclusive")

        if corner == 0:
            val = self.X[0, 0, 0]
        elif corner == 1:
            val = self.X[-1, 0, 0]
        elif corner == 2:
            val = self.X[0, -1, 0]
        elif corner == 3:
            val = self.X[-1, -1, 0]
        elif corner == 4:
            val = self.X[0, 0, -1]
        elif corner == 5:
            val = self.X[-1, 0, -1]
        elif corner == 6:
            val = self.X[0, -1, -1]
        elif corner == 7:
            val = self.X[-1, -1, -1]

        return val

    def getOrigValuesFace(self, face):
        """For a given face index, face, return the 4 corners and the
        values of the midpoints of the 4 edges on that face.

        Parameters
        ----------
        face : int
            Index of face, 0<=face<=5

        Returns
        -------
        coords : array of size (8, ndim)
            The first 4 entries are the corner, and the last 4 are the
            midpoints.
        """
        if face not in range(0, 6):
            raise Error("Face must be in range 0..5 inclusive")

        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu - 1) // 2, (self.Nu - 1) // 2]
        else:
            midu = [self.Nu // 2, self.Nu // 2 - 1]

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv - 1) // 2, (self.Nv - 1) // 2]
        else:
            midv = [self.Nv // 2, self.Nv // 2 - 1]

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw - 1) // 2, (self.Nw - 1) // 2]
        else:
            midw = [self.Nw // 2, self.Nw // 2 - 1]

        if face == 0:
            values = [
                self.X[0, 0, 0],
                self.X[-1, 0, 0],
                self.X[0, -1, 0],
                self.X[-1, -1, 0],
                0.5 * (self.X[midu[0], 0, 0] + self.X[midu[1], 0, 0]),
                0.5 * (self.X[midu[0], -1, 0] + self.X[midu[1], -1, 0]),
                0.5 * (self.X[0, midv[0], 0] + self.X[0, midv[1], 0]),
                0.5 * (self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0]),
            ]
        elif face == 1:
            values = [
                self.X[0, 0, -1],
                self.X[-1, 0, -1],
                self.X[0, -1, -1],
                self.X[-1, -1, -1],
                0.5 * (self.X[midu[0], 0, -1] + self.X[midu[1], 0, -1]),
                0.5 * (self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1]),
                0.5 * (self.X[0, midv[0], -1] + self.X[0, midv[1], -1]),
                0.5 * (self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1]),
            ]
        elif face == 2:
            values = [
                self.X[0, 0, 0],
                self.X[0, -1, 0],
                self.X[0, 0, -1],
                self.X[0, -1, -1],
                0.5 * (self.X[0, midv[0], 0] + self.X[0, midv[1], 0]),
                0.5 * (self.X[0, midv[0], -1] + self.X[0, midv[1], -1]),
                0.5 * (self.X[0, 0, midw[0]] + self.X[0, 0, midw[1]]),
                0.5 * (self.X[0, -1, midw[0]] + self.X[0, -1, midw[1]]),
            ]
        elif face == 3:
            values = [
                self.X[-1, 0, 0],
                self.X[-1, -1, 0],
                self.X[-1, 0, -1],
                self.X[-1, -1, -1],
                0.5 * (self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0]),
                0.5 * (self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1]),
                0.5 * (self.X[-1, 0, midw[0]] + self.X[-1, 0, midw[1]]),
                0.5 * (self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]]),
            ]
        elif face == 4:
            values = [
                self.X[0, 0, 0],
                self.X[-1, 0, 0],
                self.X[0, 0, -1],
                self.X[-1, 0, -1],
                0.5 * (self.X[midu[0], 0, 0] + self.X[midu[1], 0, 0]),
                0.5 * (self.X[midu[0], 0, -1] + self.X[midu[1], 0, -1]),
                0.5 * (self.X[0, 0, midw[0]] + self.X[0, 0, midw[1]]),
                0.5 * (self.X[-1, 0, midw[0]] + self.X[-1, 0, midw[1]]),
            ]
        elif face == 5:
            values = [
                self.X[0, -1, 0],
                self.X[-1, -1, 0],
                self.X[0, -1, -1],
                self.X[-1, -1, -1],
                0.5 * (self.X[midu[0], -1, 0] + self.X[midu[1], -1, 0]),
                0.5 * (self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1]),
                0.5 * (self.X[0, -1, midw[0]] + self.X[0, -1, midw[1]]),
                0.5 * (self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]]),
            ]

        return numpy.array(values)

    def getMidPointEdge(self, edge):
        """Get the midpoint of the edge using the original data.

        Parameters
        ----------
        edge : int
            Edge index. Must be 0<edge<11.

        Returns
        -------
        midpoint : array of length nDim
            Mid point of edge
        """
        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu - 1) // 2, (self.Nu - 1) // 2]
        else:
            midu = [self.Nu // 2, self.Nu // 2 - 1]

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv - 1) // 2, (self.Nv - 1) // 2]
        else:
            midv = [self.Nv // 2, self.Nv // 2 - 1]

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw - 1) // 2, (self.Nw - 1) // 2]
        else:
            midw = [self.Nw // 2, self.Nw // 2 - 1]

        if edge == 0:
            val = self.X[midu[0], 0, 0] + self.X[midu[1], 0, 0]
        elif edge == 1:
            val = self.X[midu[0], -1, 0] + self.X[midu[1], -1, 0]
        elif edge == 2:
            val = self.X[0, midv[0], 0] + self.X[0, midv[1], 0]
        elif edge == 3:
            val = self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0]
        elif edge == 4:
            val = self.X[midu[0], 0, -1] + self.X[midu[1], 0, -1]
        elif edge == 5:
            val = self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1]
        elif edge == 6:
            val = self.X[0, midv[0], -1] + self.X[0, midv[1], -1]
        elif edge == 7:
            val = self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1]
        elif edge == 8:
            val = self.X[0, 0, midw[0]] + self.X[0, 0, midw[1]]
        elif edge == 9:
            val = self.X[-1, 0, midw[0]] + self.X[-1, 0, midw[1]]
        elif edge == 10:
            val = self.X[0, -1, midw[0]] + self.X[0, -1, midw[1]]
        elif edge == 11:
            val = self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]]

        return val

    def getMidPointFace(self, face):
        """Get the midpoint of the face using the original data.

        Parameters
        ----------
        face : int
            Face index. Must be 0, 1, 2, 3, 4 or 5

        Returns
        -------
        midpoint : array of length nDim
            Mid point of face
        """
        if face not in range(0, 6):
            raise Error("Face must be in range 0..5 inclusive")
        if not self.origData:
            raise Error("No original data for this surface")

        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu - 1) // 2, (self.Nu - 1) // 2]
        else:
            midu = [self.Nu // 2, self.Nu // 2 - 1]

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv - 1) // 2, (self.Nv - 1) // 2]
        else:
            midv = [self.Nv // 2, self.Nv // 2 - 1]

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw - 1) // 2, (self.Nw - 1) // 2]
        else:
            midw = [self.Nw // 2, self.Nw // 2 - 1]

        if face == 0:
            val = 0.25 * (
                self.X[midu[0], midv[0], 0]
                + self.X[midu[1], midv[0], 0]
                + self.X[midu[0], midv[1], 0]
                + self.X[midu[1], midv[1], 0]
            )
        elif face == 1:
            val = 0.25 * (
                self.X[midu[0], midv[0], -1]
                + self.X[midu[1], midv[0], -1]
                + self.X[midu[0], midv[1], -1]
                + self.X[midu[1], midv[1], -1]
            )
        elif face == 2:
            val = 0.25 * (
                self.X[0, midv[0], midw[0]]
                + self.X[0, midv[1], midw[0]]
                + self.X[0, midv[0], midw[1]]
                + self.X[0, midv[1], midw[1]]
            )
        elif face == 3:
            val = 0.25 * (
                self.X[-1, midv[0], midw[0]]
                + self.X[-1, midv[1], midw[0]]
                + self.X[-1, midv[0], midw[1]]
                + self.X[-1, midv[1], midw[1]]
            )
        elif face == 4:
            val = 0.25 * (
                self.X[midu[0], 0, midw[0]]
                + self.X[midu[1], 0, midw[0]]
                + self.X[midu[0], 0, midw[1]]
                + self.X[midu[1], 0, midw[1]]
            )
        elif face == 5:
            val = 0.25 * (
                self.X[midu[0], -1, midw[0]]
                + self.X[midu[1], -1, midw[0]]
                + self.X[midu[0], -1, midw[1]]
                + self.X[midu[1], -1, midw[1]]
            )

        return val

    def setFaceSurfaces(self):
        """Create face spline objects for each of the faces"""

        self.faceSurfaces[0] = Surface(ku=self.ku, kv=self.kv, tu=self.tu, tv=self.tv, coef=self.coef[:, :, 0, :])
        self.faceSurfaces[1] = Surface(ku=self.ku, kv=self.kv, tu=self.tu, tv=self.tv, coef=self.coef[:, :, -1, :])
        self.faceSurfaces[2] = Surface(ku=self.ku, kv=self.kw, tu=self.tu, tv=self.tw, coef=self.coef[:, 0, :, :])
        self.faceSurfaces[3] = Surface(ku=self.ku, kv=self.kw, tu=self.tu, tv=self.tw, coef=self.coef[:, -1, :, :])
        self.faceSurfaces[4] = Surface(ku=self.kv, kv=self.kw, tu=self.tv, tv=self.tw, coef=self.coef[0, :, :, :])
        self.faceSurfaces[5] = Surface(ku=self.kv, kv=self.kw, tu=self.tv, tv=self.tw, coef=self.coef[-1, :, :, :])

    def setEdgeCurves(self):
        """Create edge spline objects for each edge"""

        self.edgeCurves[0] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, 0, 0])
        self.edgeCurves[1] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, -1, 0])
        self.edgeCurves[2] = Curve(k=self.kv, t=self.tv, coef=self.coef[0, :, 0])
        self.edgeCurves[3] = Curve(k=self.kv, t=self.tv, coef=self.coef[-1, :, 0])
        self.edgeCurves[4] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, 0, -1])
        self.edgeCurves[5] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, -1, -1])
        self.edgeCurves[6] = Curve(k=self.kv, t=self.tv, coef=self.coef[0, :, -1])
        self.edgeCurves[7] = Curve(k=self.kv, t=self.tv, coef=self.coef[-1, :, -1])
        self.edgeCurves[8] = Curve(k=self.kw, t=self.tw, coef=self.coef[0, 0, :])
        self.edgeCurves[9] = Curve(k=self.kw, t=self.tw, coef=self.coef[-1, 0, :])
        self.edgeCurves[10] = Curve(k=self.kw, t=self.tw, coef=self.coef[0, -1, :])
        self.edgeCurves[11] = Curve(k=self.kw, t=self.tw, coef=self.coef[-1, -1, :])

    def getBasisPt(self, u, v, w, vals, istart, colInd, lIndex):
        """This function should only be called from pyBlock The purpose
        is to compute the basis function for a u, v, w point and add
        it to pyBlocks's global dPt/dCoef
        matrix. vals, rowPtr, colInd is the CSR data and lIndex in
        the local -> global mapping for this volume"""

        return libspline.getbasisptvolume(
            u, v, w, self.tu, self.tv, self.tw, self.ku, self.kv, self.kw, vals, colInd, istart, lIndex.T
        )

    def __call__(self, u, v, w):
        """
        Equivalant to getValue()
        """
        return self.getValue(u, v, w)

    def getValue(self, u, v, w):
        """Get the value at the volume points(s) u, v, w. This is the
        main evaluation routine for the volume object.

        Parameters
        ----------
        u : scalar, vector or matrix or tensor of values
            u position
        v : scalar, vector or matrix or tensor of values
            v position
        w : scalar, vector or matrix or tensor of values
            w position

        Returns
        -------
        values : scalar, vector, matrix or tensor of values
           The spline evaluation at (u, v, w)
        """
        u = numpy.atleast_3d(u).T
        v = numpy.atleast_3d(v).T
        w = numpy.atleast_3d(w).T

        if not u.shape == v.shape == w.shape:
            raise Error("u and v must have the same shape")

        vals = libspline.eval_volume(u, v, w, self.tu, self.tv, self.tw, self.ku, self.kv, self.kw, self.coef.T)
        return vals.squeeze().T

    def getValueEdge(self, edge, s):
        """Get the value at the volume points(s) u, v, w

        Parameters
        ----------
        edge : int
             Index of edge. Must be between 0 and 11.
        s : float or array
            Parameter position(s) along edge to evaluate.

        Returns
        -------
        values : array
            Array of values evaluated along edge.
        """
        if edge == 0:
            u = s
            v = self.vmin
            w = self.wmin
        elif edge == 1:
            u = s
            v = self.vmax
            w = self.wmin
        elif edge == 2:
            u = self.umin
            v = s
            w = self.wmax
        elif edge == 3:
            u = self.umax
            v = s
            w = self.umin
        elif edge == 4:
            u = s
            v = self.vmin
            w = self.wmax
        elif edge == 5:
            u = s
            v = self.vmax
            w = self.wmax
        elif edge == 6:
            u = self.umin
            v = s
            w = self.wmax
        elif edge == 7:
            u = self.umax
            v = s
            w = self.wmax
        elif edge == 8:
            u = self.umin
            v = self.vmin
            w = s
        elif edge == 9:
            u = self.umax
            v = self.vmin
            w = s
        elif edge == 10:
            u = self.umin
            v = self.vmax
            w = s
        elif edge == 11:
            u = self.umax
            v = self.vmax
            w = s

        u = numpy.atleast_3d(u).T
        v = numpy.atleast_3d(v).T
        w = numpy.atleast_3d(w).T

        if not u.shape == v.shape == w.shape:
            raise Error("u, v, and w must have the same shape")

        vals = libspline.eval_volume(u, v, w, self.tu, self.tv, self.tw, self.ku, self.kv, self.kw, self.coef.T)
        return vals.squeeze().T

    def getBounds(self):
        """Determine the extents of the volume

        Returns
        -------
        xMin : array of length 3
            Lower corner of the bounding box
        xMax : array of length 3
            Upper corner of the bounding box
        """
        if self.nDim != 3:
            raise Error("getBounds is only defined for nDim = 3")

        cx = self.coef[:, :, :, 0].flatten()
        cy = self.coef[:, :, :, 1].flatten()
        cz = self.coef[:, :, :, 2].flatten()

        Xmin = numpy.zeros(self.nDim)
        Xmin[0] = min(cx)
        Xmin[1] = min(cy)
        Xmin[2] = min(cz)

        Xmax = numpy.zeros(self.nDim)
        Xmax[0] = max(cx)
        Xmax[1] = max(cy)
        Xmax[2] = max(cz)

        return Xmin, Xmax

    def projectPoint(self, x0, nIter=25, eps=1e-10, **kwargs):
        """
        Project a point x0 or points x0 onto the volume and return
        parametric positions

        Parameters
        ----------
        x0 : array of length 3 or array of size (N, 3)
            Points to embed in the volume. If the points do not
            **actually** lie in the volume, the closest point is
            returned
        nIter : int
            Maximum number of Newton iterations to perform
        eps : float
            Tolerance for the Newton iteration
        u : float or array of len(X0)
            Optional initial guess for u position.
        v : float or array of len(X0)
            Optional initial guess for v position.
        w : float or array of len(X0)
            Optional initial guess for w position.

        Returns
        -------
        u : float or array of length N
            u parametric position of closest point
        v : float or array of length N
            v parametric position of closest point
        w : float or array of length N
            w parametric position of closest point
        D : float or array of length N
            Distance between projected point and closest point
            If the points are 'inside' the volume, D should
            be less than eps.
        """

        x0 = numpy.atleast_2d(x0)

        if "u" in kwargs and "v" in kwargs and "w" in kwargs:
            u = numpy.atleast_1d(kwargs["u"])
            v = numpy.atleast_1d(kwargs["v"])
            w = numpy.atleast_1d(kwargs["w"])
        else:
            u = -1 * numpy.ones(len(x0))
            v = -1 * numpy.ones(len(x0))
            w = -1 * numpy.ones(len(x0))

        if not len(x0) == len(u) == len(v) == len(w):
            raise Error("The length of x0 and u, v, w must be the same")

        # If necessary get brute-force starting point
        if numpy.any(u < 0) or numpy.any(u > 1) or numpy.any(v < 0) or numpy.any(v > 1):
            self.computeData()
            u, v, w = libspline.point_volume_start(x0.real.T, self.udata, self.vdata, self.wdata, self.data.T)
        D = numpy.zeros_like(x0)
        for i in range(len(x0)):
            u[i], v[i], w[i], D[i] = libspline.point_volume(
                x0[i].real,
                self.tu,
                self.tv,
                self.tw,
                self.ku,
                self.kv,
                self.kw,
                self.coef.T,
                nIter,
                eps,
                u[i],
                v[i],
                w[i],
            )

        return u.squeeze(), v.squeeze(), w.squeeze(), D.squeeze()

    def computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """
        # Only recompute if it doesn't exist already
        if self.data is None:
            self.edgeCurves[0].calcInterpolatedGrevillePoints()
            self.udata = self.edgeCurves[0].sdata
            self.edgeCurves[2].calcInterpolatedGrevillePoints()
            self.vdata = self.edgeCurves[2].sdata
            self.edgeCurves[8].calcInterpolatedGrevillePoints()
            self.wdata = self.edgeCurves[8].sdata
            U = numpy.zeros((len(self.udata), len(self.vdata), len(self.wdata)))
            V = numpy.zeros((len(self.udata), len(self.vdata), len(self.wdata)))
            W = numpy.zeros((len(self.udata), len(self.vdata), len(self.wdata)))
            for i in range(len(self.udata)):
                for j in range(len(self.vdata)):
                    for k in range(len(self.wdata)):
                        U[i, j, k] = self.udata[i]
                        V[i, j, k] = self.vdata[j]
                        W[i, j, k] = self.wdata[k]
            self.data = self.getValue(U, V, W)

    def insertKnot(self, direction, s, r):
        """
        Insert a knot into the volume along either u, v, w:

        Parameters
        ----------
        direction : str
            Parameteric direction to insert. Either 'u', 'v', or 'w'
        s : float
            Parametric position along 'direction' to insert
        r : int
        Desired number of times to insert.

        Returns
        -------
        r : int
            The **actual** number of times the knot was inserted.
        """
        if direction not in ["u", "v", "w"]:
            raise Error("Direction must be one of 'u' or 'v' or 'w'")

        s = checkInput(s, "s", float, 0)
        r = checkInput(r, "r", int, 0)
        if s <= 0.0:
            return
        if s >= 1.0:
            return

        # This is relatively inefficient, but we'll do it for
        # simplicity just call insertknot for each slab

        if direction == "u":
            # Insert once to know how many times it was actually inserted
            # so we know how big to make the new coef:
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tu, self.ku, self.coef[:, 0, 0].T)

            newCoef = numpy.zeros((self.nCtlu + actualR, self.nCtlv, self.nCtlw, self.nDim))
            for k in range(self.nCtlvw):
                for j in range(self.nCtlv):
                    actualR, tNew, coefSlice, breakPt = libspline.insertknot(
                        s, r, self.tu, self.ku, self.coef[:, j, k].T
                    )
                    newCoef[:, j, k] = coefSlice[:, 0 : self.nCtlu + actualR].T

            self.tu = tNew[0 : self.nCtlu + self.ku + actualR]
            self.nCtlu = self.nCtlu + actualR

        elif direction == "v":
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tv, self.kv, self.coef[0, :, 0].T)

            newCoef = numpy.zeros((self.nCtlu, self.nCtlv + actualR, self.nCtlw, self.nDim))

            for k in range(self.nCtlw):
                for i in range(self.nCtlu):
                    actualR, tNew, coefSlice, breakPt = libspline.insertknot(
                        s, r, self.tv, self.kv, self.coef[i, :, k].T
                    )
                    newCoef[i, :, k] = coefSlice[:, 0 : self.nCtlv + actualR].T

            self.tv = tNew[0 : self.nCtlv + self.kv + actualR]
            self.nCtlv = self.nCtlv + actualR

        elif direction == "w":
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tw, self.kw, self.coef[0, 0, :].T)

            newCoef = numpy.zeros((self.nCtlu, self.nCtlv, self.nCtlw + actualR, self.nDim))

            for j in range(self.nCtlv):
                for i in range(self.nCtlu):
                    actualR, tNew, coefSlice, breakPt = libspline.insertknot(
                        s, r, self.tw, self.kw, self.coef[i, j, :].T
                    )
                    newCoef[i, j, :] = coefSlice[:, 0 : self.nCtlw + actualR].T

            self.tw = tNew[0 : self.nCtlw + self.kw + actualR]
            self.nCtlw = self.nCtlw + actualR

        self.coef = newCoef

        # break_pt is converted to zero based ordering here!!!
        return actualR, breakPt - 1

    def writeTecplot(self, fileName, vols=True, coef=True, orig=False):
        """Write the volume to a tecplot data file.

        Parameters
        ----------
        fileName : str
            Tecplot filename. Should end in .dat
        vols : bool
            Flag specifiying whether the interpolated volume should
            be used. This is usually True if you want to get an
            approximation of the entire volume.
        coef : bool
            Flag specifiying if the control points are to be plotted
        orig : bool
            Flag specifiying if original data (used for fitting) is
            to be included. If on original data exists, this argument
            is ignored.
        """
        f = openTecplot(fileName, self.nDim)
        if vols:
            self.computeData()
            writeTecplot3D(f, "interpolated", self.data)
        if coef:
            writeTecplot3D(f, "control_pts", self.coef)
        if orig and self.origData:
            writeTecplot3D(f, "orig_data", self.X)
        closeTecplot(f)


# For backwards compatibility, the old curve, surface and volume definitions:
def curve(*args, **kwargs):
    warnings.warn("pySpline.curve has been changed to Curve()")
    return Curve(*args, **kwargs)


def surface(*args, **kwargs):
    warnings.warn("pySpline.surface has been changed to Surface()")
    return Surface(*args, **kwargs)


def volume(*args, **kwargs):
    warnings.warn("pySpline.volume has been changed to Volume()")
    return Volume(*args, **kwargs)

    # ----------------------------------------------------------------------
    #                     Misc Helper Functions
    # ----------------------------------------------------------------------


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
    tu = [0, 0, 1, 1]
    tv = [0, 0, 1, 1]
    tw = [0, 0, 1, 1]
    ku = 2
    kv = 2
    kw = 2

    if len(args) == 1:
        return Volume(coef=args[0], tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    elif len(args) == 2:
        xmin = numpy.array(args[0]).astype("d")
        xmax = numpy.array(args[1]).astype("d")

        xLow = xmin[0]
        xHigh = xmax[0]
        yLow = xmin[1]
        yHigh = xmax[1]
        zLow = xmin[2]
        zHigh = xmax[2]

        coef = numpy.zeros((2, 2, 2, 3))
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
        raise Error(
            "An unknown number of arguments was passed to\
 trilinear  Volume"
        )


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
    if len(args) == 1:
        # One argument passed in ... assume its X
        if len(args[0]) != 4:
            raise Error("A single argument passed to bilinear surface must contain 4 points and be of size (4, 3)")
        coef = numpy.zeros((2, 2, 3))
        coef[0, 0] = args[0][0]
        coef[1, 0] = args[0][1]
        coef[0, 1] = args[0][2]
        coef[1, 1] = args[0][3]
        return Surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)
    else:
        # Assume 4 arguments
        coef = numpy.zeros([2, 2, 3])
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

    if len(args) == 2:
        # Its a two-point type
        return Curve(coef=[args[0], args[1]], k=2, t=[0, 0, 1, 1])
    elif len(args) == 1:
        if len(args[0]) == 2:  # its X
            return Curve(coef=args[0], k=2, t=[0, 0, 1, 1])
        elif "dir" in kwargs:
            # We have point and direction
            if "length" in kwargs:
                x2 = args[0] + kwargs["dir"] / numpy.linalg.norm(kwargs["dir"]) * kwargs["length"]
            else:
                x2 = args[0] + kwargs["dir"]

            return Curve(coef=[args[0], x2], k=2, t=[0, 0, 1, 1])
        else:
            Error("Error: dir must be specified if only 1 argument is given")
