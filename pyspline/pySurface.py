# External modules
import numpy as np
from scipy.sparse import linalg

# Local modules
from . import libspline
from .pyCurve import Curve
from .utils import Error, _assembleMatrix, checkInput, closeTecplot, openTecplot, writeTecplot1D, writeTecplot2D


class Surface(object):

    """
    Create an instance of a b-spline surface. There are two
    ways to initialize the class

    * **Creation**: Create an instance of the Surface class
      directly by supplying the required information: kwargs MUST
      contain the following information: ``ku, kv, tu, tv, coef``.

    * **LMS/Interpolation**: Create an instance of the Surface class by
      using an interpolating spline to given data points or a LMS
      spline. The following keyword argument information is required:

      1. ``ku`` and ``kv`` Spline Orders

      2. ``X`` real array size (Nu, Nv, nDim) of data to fit. **OR**
          1. ``x`` (2D) and ``y`` (2D)  for 2D surface interpolation
          2. ``x`` (3D) and ``y`` (3D) and ``z`` (3) for 3D surface

      3. ``u``, ``v`` real array of size (Nu, Nv). Optional

      4. ``nCtlu``, ``nCtlv`` Specify number of control points. Only for
         LMS fitting.

    Parameters
    ----------
    ku : int
       Spline order in u
    kv : int
       Spline order in v
    nCtlu : int
       Number of control points in u
    nCtlv : int
       Number of control points in v
    coef : array, size (nCtlu, nCtl, nDim)
       b-spline coefficient array.
    tu : array, size(nCtlu + ku)
       knot array in u
    tv : array, size(nCtlv + kv)
       knot array in v
    X : array, size (Nu, Nv, ndim)
       Full data array to fit
    x : array, size (Nu, Nv)
       Just x data to fit/interpolate
    y : array, size (Nu, Nv)
       Just y data to fit/interpolate
    u : array, size (Nu, Nv)
       Explict u parameters to use. Optional.
    v : array, size (Nu, Nv)
       Explict v parameters to use. Optional.
    scaledParams : bool
       default is to use u,v for parameterization. If true use u,v as well.
       If false, use U,V.
    nIter : int
       Number of Hoscheks parater corrections to run

    Notes
    -----
    The orientation of the nodes, edges and faces is the same as the
    **bottom** surface as described in :class:`Volume` documentation.
    """

    def __init__(self, recompute=True, **kwargs):

        self.name = None
        self.edgeCurves = [None, None, None, None]
        self.data = None
        self.udata = None
        self.vdata = None
        if "ku" in kwargs and "kv" in kwargs and "tu" in kwargs and "tv" in kwargs and "coef" in kwargs:
            self.X = None
            self.u = None
            self.v = None
            self.U = None
            self.V = None
            self.ku = checkInput(kwargs["ku"], "ku", int, 0)
            self.kv = checkInput(kwargs["kv"], "kv", int, 0)
            self.coef = checkInput(kwargs["coef"], "coef", float, 3)
            self.nCtlu = self.coef.shape[0]
            self.nCtlv = self.coef.shape[1]
            self.tu = checkInput(kwargs["tu"], "tu", float, 1, (self.nCtlu + self.ku,))
            self.tv = checkInput(kwargs["tv"], "tv", float, 1, (self.nCtlv + self.kv,))
            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.nDim = self.coef.shape[2]
            self.origData = False
            self.setEdgeCurves()
            self.interp = False
            return
        elif "localInterp" in kwargs:
            # Local, non-global cubic interpolation. See The Nurbs
            # Book section 9.3.5 "Local Bicubic Surface Interpolation"
            self.localInterp = True
            self.ku = 4
            self.kv = 4

            # Do some checking on the number of control points
            if "X" in kwargs:
                self.X = np.array(kwargs["X"])
                self.nDim = self.X.shape[2]
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = np.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 3))
                self.X[:, :, 0] = kwargs["x"]
                self.X[:, :, 1] = kwargs["y"]
                self.X[:, :, 2] = kwargs["z"]
                self.nDim = 3
            else:
                raise Error("Error: X (or x, y, z)  MUST be defined for task localInterp!")

            self.origData = True
            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]

            self.u, self.v, self.U, self.V = self.calcParameterization()

            # Contains the T^u_kl, T^v_kl and D^uv_kl values
            Td = np.zeros((self.Nu, self.Nv, 5, 3))

            def getT(Q, s):
                N = len(Q)
                T = np.zeros_like(Q)
                qq = np.zeros_like(Q)
                deltaS = np.zeros(N)
                for i in range(1, N):
                    deltaS[i] = s[i] - s[i - 1]
                    qq[i, :] = Q[i] - Q[i - 1]

                for i in range(1, N - 1):
                    a = deltaS[i] / (deltaS[i] + deltaS[i + 1])
                    T[i] = (1 - a) * qq[i] + a * qq[i + 1]

                # Do the start and end points: (eqn: 9.32, The NURBS book)
                T[0] = 2 * qq[1] / deltaS[1] - T[1]
                T[-1] = 2 * qq[-1] / deltaS[-1] - T[-2]
                for i in range(N):
                    T[i] /= np.linalg.norm(T[i]) + 1e-16
                return T

            def interpKnots(u):
                t = np.zeros(2 * len(u) + 2 + 2)
                t[0:4] = 0.0
                t[-4:] = 1.0

                ii = 4
                for i in range(1, len(u) - 1):
                    t[ii] = u[i]
                    t[ii + 1] = u[i]
                    ii = ii + 2
                return t

            def bezierCoef(Q, T, length, s):
                N = len(Q)
                coef = np.zeros((3 * N - 2, self.nDim))
                for i in range(N):
                    coef[3 * i] = Q[i].copy()

                for i in range(N - 1):
                    a = length * (s[i + 1] - s[i])
                    coef[3 * i + 1] = Q[i] + a / 3.0 * T[i]
                    coef[3 * i + 2] = Q[i + 1] - a / 3.0 * T[i + 1]
                return coef

            def getLength(Q):
                length = 0
                for i in range(len(Q) - 1):
                    length += np.linalg.norm(Q[i + 1] - Q[i])
                return length

            def getD(Q, s):
                N = len(Q)
                D = np.zeros_like(Q)
                dd = np.zeros_like(D)
                deltaS = np.zeros(N)
                for i in range(1, N):
                    deltaS[i] = s[i] - s[i - 1]
                    dd[i] = (Q[i] - Q[i - 1]) / (deltaS[i] + 1e-16)

                for i in range(1, N - 1):
                    a = deltaS[i] / (deltaS[i] + deltaS[i + 1])
                    D[i] = (1 - a) * dd[i] + a * dd[i + 1]

                D[0] = 2 * dd[1] - D[1]
                D[-1] = 2 * dd[-1] - D[-2]
                return D

            # -------- Algorithm A9.5 -------------

            rowLen = np.zeros(self.Nv)
            colLen = np.zeros(self.Nu)

            for j in range(self.Nv):
                # Compute the U-tangent values
                Td[:, j, 0] = getT(self.X[:, j], self.u)
                rowLen[j] = getLength(self.X[:, j])

            for i in range(self.Nu):
                # Compute the U-tangent values
                Td[i, :, 1] = getT(self.X[i, :], self.v)
                colLen[i] = getLength(self.X[i, :])

            self.tu = interpKnots(self.u)
            self.tv = interpKnots(self.v)

            # This contains all the coef...including the ones we will
            # eventually knock out.
            coef = np.zeros((3 * self.Nu - 2, 3 * self.Nv - 2, self.nDim))

            scaledParams = kwargs.pop("scaledParams", True)

            for i in range(self.Nu):
                if scaledParams:
                    coef[3 * i, :] = bezierCoef(self.X[i, :], Td[i, :, 1], colLen[i], self.v)
                else:
                    coef[3 * i, :] = bezierCoef(self.X[i, :], Td[i, :, 1], colLen[i], self.V[i, :])

            for j in range(self.Nv):
                if scaledParams:
                    coef[:, 3 * j] = bezierCoef(self.X[:, j], Td[:, j, 0], rowLen[j], self.u)
                else:
                    coef[:, 3 * j] = bezierCoef(self.X[:, j], Td[:, j, 0], rowLen[j], self.U[:, j])

            # Now compute the cross derivatives, assuming that the uv
            # derivates can be averaged.
            for j in range(self.Nv):
                Td[:, j, 0] *= rowLen[j]
                Td[:, j, 3] = getD(Td[:, j, 0], self.u)
            for i in range(self.Nu):
                Td[i, :, 1] *= colLen[i]
                Td[i, :, 4] = getD(Td[i, :, 1], self.v)
            Td = 0.5 * (Td[:, :, 3] + Td[:, :, 4])
            for i in range(self.Nu):
                for j in range(self.Nv):
                    Td[i, j] /= np.linalg.norm(Td[i, j]) + 1e-16

            for j in range(self.Nv - 1):
                for i in range(self.Nu - 1):
                    gam = (self.u[i + 1] - self.u[i]) * (self.v[j + 1] - self.v[j]) / 9.0
                    ii = 3 * i
                    jj = 3 * j
                    coef[ii + 1, jj + 1] = gam * Td[i, j] + coef[ii, jj + 1] + coef[ii + 1, jj] - coef[ii, jj]

                    coef[ii + 2, jj + 1] = (
                        -gam * Td[i + 1, j] + coef[ii + 3, jj + 1] - coef[ii + 3, jj] + coef[ii + 2, jj]
                    )

                    coef[ii + 1, jj + 2] = (
                        -gam * Td[i, j + 1] + coef[ii + 1, jj + 3] - coef[ii, jj + 3] + coef[ii, jj + 2]
                    )

                    coef[ii + 2, jj + 2] = (
                        gam * Td[i + 1, j + 1] + coef[ii + 2, jj + 3] + coef[ii + 3, jj + 2] - coef[ii + 3, jj + 3]
                    )

            # Coef was just used for convience. We will now extract
            # just the values we need. We can do this with *super*
            # fancy indexing.
            def getIndex(N):
                ind = np.zeros(2 * (N - 1) + 2, "intc")
                ind[0] = 0
                ind[-1] = 3 * N - 3
                ii = 1
                jj = 1
                for _i in range(N - 1):
                    ind[ii] = jj
                    ind[ii + 1] = jj + 1
                    ii += 2
                    jj += 3

                return ind

            coef = coef[:, getIndex(self.Nv), :]
            self.coef = coef[getIndex(self.Nu)]

            self.nCtlu = self.coef.shape[0]
            self.nCtlv = self.coef.shape[1]
            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.origData = True
            self.setEdgeCurves()
            self.interp = True

        else:  # We have LMS/Interpolate
            # Do some checking on the number of control points
            if not (
                "ku" in kwargs
                and "kv" in kwargs
                and (
                    "X" in kwargs
                    or "x" in kwargs
                    or ("x" in kwargs and "y" in kwargs)
                    or ("x" in kwargs and "y" in kwargs and "z" in kwargs)
                )
            ):
                raise ValueError(
                    "ku, kv, and X (or x, or x and y, or x and y and z MUST be defined for task lms or interpolate"
                )

            if "X" in kwargs:
                self.X = np.array(kwargs["X"])
                if len(self.X.shape) == 1:
                    self.nDim = 1
                else:
                    self.nDim = self.X.shape[2]
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = np.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 3))
                self.X[:, :, 0] = kwargs["x"]
                self.X[:, :, 1] = kwargs["y"]
                self.X[:, :, 2] = kwargs["z"]
                self.nDim = 3
            elif "x" in kwargs and "y" in kwargs:
                self.X = np.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 2))
                self.X[:, :, 0] = kwargs["x"]
                self.X[:, :, 1] = kwargs["y"]
                self.nDim = 2
            elif "x" in kwargs:
                self.X = np.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 1))
                self.X[:, :, 0] = kwargs["x"]
                self.nDim = 1
            # enf if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]

            self.ku = checkInput(kwargs["ku"], "ku", int, 0)
            self.kv = checkInput(kwargs["kv"], "kv", int, 0)

            if "nCtlu" in kwargs and "nCtlv" in kwargs:
                self.nCtlu = checkInput(kwargs["nCtlu"], "nCtlu", int, 0)
                self.nCtlv = checkInput(kwargs["nCtlv"], "nCtlv", int, 0)
                self.interp = False
            else:
                self.nCtlu = self.Nu
                self.nCtlv = self.Nv
                self.interp = True

            self.origData = True

            # Sanity Check on Inputs
            if self.nCtlu >= self.Nu:
                self.nCtlu = self.Nu
            if self.nCtlv >= self.Nv:
                self.nCtlv = self.Nv

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku:
                self.ku = self.Nu
            if self.Nv < self.kv:
                self.kv = self.Nv
            if self.nCtlu < self.ku:
                self.ku = self.nCtlu
            if self.nCtlv < self.kv:
                self.kv = self.nCtlv

            if "nIter" in kwargs:
                self.nIter = checkInput(kwargs["nIter"], "nIter", int, 0)
            else:
                self.nIter = 1

            if "u" in kwargs and "v" in kwargs:
                self.u = checkInput(kwargs["u"], "u", float, 1, (self.Nu,))
                self.v = checkInput(kwargs["v"], "v", float, 1, (self.Nv,))
                self.u = self.u / self.u[-1]
                self.v = self.v / self.v[-1]
                [self.V, self.U] = np.meshgrid(self.v, self.u)
            else:
                if self.nDim == 3:
                    self.u, self.v, self.U, self.V = self.calcParameterization()
                else:
                    raise Error(
                        "Automatic parameterization of ONLY available for spatial data in 3 dimensions. "
                        + "Please supply u and v key word arguments otherwise."
                    )

            self.umin = 0
            self.umax = 1
            self.vmin = 0
            self.vmax = 1
            self.calcKnots()
            self.coef = np.zeros((self.nCtlu, self.nCtlv, self.nDim))
            if recompute:
                self.recompute()

    def recompute(self):
        """Recompute the surface if any data has been modified"""

        vals, rowPtr, colInd = libspline.surface_jacobian_wrap(
            self.U.T, self.V.T, self.tu, self.tv, self.ku, self.kv, self.nCtlu, self.nCtlv
        )

        N = _assembleMatrix(vals, colInd, rowPtr, (self.Nu * self.Nv, self.nCtlu * self.nCtlv))

        self.coef = np.zeros((self.nCtlu, self.nCtlv, self.nDim))
        if self.interp:
            # Factorize once for efficiency
            solve = linalg.dsolve.factorized(N)
            for idim in range(self.nDim):
                self.coef[:, :, idim] = solve(self.X[:, :, idim].flatten()).reshape([self.nCtlu, self.nCtlv])
        else:
            solve = linalg.dsolve.factorized(N.transpose() * N)
            for idim in range(self.nDim):
                rhs = N.transpose() * self.X[:, :, idim].flatten()
                self.coef[:, :, idim] = solve(rhs).reshape([self.nCtlu, self.nCtlv])

        self.setEdgeCurves()

    def calcParameterization(self):
        """Compute a spatial parameterization"""

        u = np.zeros(self.Nu, "d")
        U = np.zeros((self.Nu, self.Nv), "d")
        singularSounter = 0
        # loop over each v, and average the 'u' parameter
        for j in range(self.Nv):
            temp = np.zeros(self.Nu, "d")

            for i in range(self.Nu - 1):
                temp[i + 1] = temp[i] + np.linalg.norm(self.X[i, j] - self.X[i + 1, j])

            if temp[-1] == 0:  # We have a singular point
                singularSounter += 1
                temp[:] = 0.0
                U[:, j] = np.linspace(0, 1, self.Nu)
            else:
                temp = temp / temp[-1]
                U[:, j] = temp.copy()

            u += temp  # accumulate the u-parameter calcs for each j

        u = u / (self.Nv - singularSounter)  # divide by the number of 'j's we had

        v = np.zeros(self.Nv, "d")
        V = np.zeros((self.Nu, self.Nv), "d")
        singularSounter = 0
        # loop over each v, and average the 'u' parameter
        for i in range(self.Nu):
            temp = np.zeros(self.Nv, "d")
            for j in range(self.Nv - 1):
                temp[j + 1] = temp[j] + np.linalg.norm(self.X[i, j] - self.X[i, j + 1])

            if temp[-1] == 0:  # We have a singular point
                singularSounter += 1
                temp[:] = 0.0
                V[i, :] = np.linspace(0, 1, self.Nv)
            else:
                temp = temp / temp[-1]
                V[i, :] = temp.copy()

            v += temp  # accumulate the v-parameter calcs for each i

        v = v / (self.Nu - singularSounter)  # divide by the number of 'i's we had

        return u, v, U, V

    def calcKnots(self):
        """Determine the knots depending on if it is inerpolated or
        an LMS fit"""
        if self.interp:
            self.tu = libspline.knots_interp(self.u, np.array([], "d"), self.ku)
            self.tv = libspline.knots_interp(self.v, np.array([], "d"), self.kv)
        else:
            self.tu = libspline.knots_lms(self.u, self.nCtlu, self.ku)
            self.tv = libspline.knots_lms(self.v, self.nCtlv, self.kv)

    def setEdgeCurves(self):
        """Create curve spline objects for each of the edges"""
        self.edgeCurves[0] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, 0])
        self.edgeCurves[1] = Curve(k=self.ku, t=self.tu, coef=self.coef[:, -1])
        self.edgeCurves[2] = Curve(k=self.kv, t=self.tv, coef=self.coef[0, :])
        self.edgeCurves[3] = Curve(k=self.kv, t=self.tv, coef=self.coef[-1, :])

    def getValueCorner(self, corner):
        """Evaluate the spline spline at one of the four corners

        Parameters
        ----------
        corner : int
            Corner index 0<=corner<=3

        Returns
        -------
        value : float
            Spline evaluated at corner
        """
        if corner not in [0, 1, 2, 3]:
            raise Error("Corner must be in range 0..3 inclusive")

        if corner == 0:
            return self.getValue(self.umin, self.vmin)
        elif corner == 1:
            return self.getValue(self.umax, self.vmin)
        elif corner == 2:
            return self.getValue(self.umin, self.vmax)
        elif corner == 3:
            return self.getValue(self.umax, self.vmax)

    def getOrigValueCorner(self, corner):
        """Return the original data for he spline at the corner if it
        exists

        Parameters
        ----------
        corner : int
            Corner index 0<=corner<=3

        Returns
        -------
        value : float
            Original value at corner
        """
        if corner not in range(0, 4):
            raise Error("Corner must be in range 0..3 inclusive")
        if not self.origData:
            raise Error("No original data for this surface")

        if corner == 0:
            return self.X[0, 0]
        elif corner == 1:
            return self.X[-1, 0]
        elif corner == 2:
            return self.X[0, -1]
        elif corner == 3:
            return self.X[-1, -1]

    def getOrigValuesEdge(self, edge):
        """Return the endpoints and the mid-point value for a given edge.

        Parameters
        ----------
        edge : int
            Edge index 0<=edge<=3

        Returns
        -------
        startValue : array size nDim
            Original value at start of edge

        midValue : array size nDim
            Original value at mid point of edge

        endValue : array size nDim
            Original value at end of edge.
        """
        if edge not in range(0, 4):
            raise Error("Edge must be in range 0..3 inclusive")
        if not self.origData:
            raise Error("No original data for this surface")

        if edge == 0:
            if np.mod(self.Nu, 2) == 1:  # Its odd
                mid = (self.Nu - 1) // 2
                return self.X[0, 0], self.X[mid, 0], self.X[-1, 0]
            else:
                Xmid = 0.5 * (self.X[self.Nu // 2, 0] + self.X[self.Nu // 2 - 1, 0])
                return self.X[0, 0], Xmid, self.X[-1, 0]
        elif edge == 1:
            if np.mod(self.Nu, 2) == 1:  # Its odd
                mid = (self.Nu - 1) // 2
                return self.X[0, -1], self.X[mid, -1], self.X[-1, -1]
            else:
                Xmid = 0.5 * (self.X[self.Nu // 2, -1] + self.X[self.Nu // 2 - 1, -1])
                return self.X[0, -1], Xmid, self.X[-1, -1]
        elif edge == 2:
            if np.mod(self.Nv, 2) == 1:  # Its odd
                mid = (self.Nv - 1) // 2
                return self.X[0, 0], self.X[0, mid], self.X[0, -1]
            else:
                Xmid = 0.5 * (self.X[0, self.Nv // 2] + self.X[0, self.Nv // 2 - 1])
                return self.X[0, 0], Xmid, self.X[0, -1]
        elif edge == 3:
            if np.mod(self.Nv, 2) == 1:  # Its odd
                mid = (self.Nv - 1) // 2
                return self.X[-1, 0], self.X[-1, mid], self.X[-1, -1]
            else:
                Xmid = 0.5 * (self.X[-1, self.Nv // 2] + self.X[-1, self.Nv // 2 - 1])
                return self.X[-1, 0], Xmid, self.X[-1, -1]

    def getValueEdge(self, edge, s):
        """
        Evaluate the spline at parametric distance s along edge 'edge'

        Parameters
        ----------
        edge : int
            Edge index 0<=edge<=3

        s : float or array
            Parameter values to evaluate

        Returns
        -------
        values : array size (nDim) or array of size (N,nDim)
            Requested spline evaluated values
        """

        return self.edgeCurves[edge](s)

    def getBasisPt(self, u, v, vals, istart, colInd, lIndex):
        """This function should only be called from pyGeo
        The purpose is to compute the basis function for
        a u, v point and add it to pyGeo's global dPt/dCoef
        matrix. vals, row_ptr, col_ind is the CSR data and
        lIndex in the local -> global mapping for this
        surface"""
        return libspline.getbasisptsurface(u, v, self.tu, self.tv, self.ku, self.kv, vals, colInd, istart, lIndex.T)

    def __call__(self, u, v):
        """
        Equivalant to getValue()
        """
        return self.getValue(u, v)

    def insertKnot(self, direction, s, r):
        """
        Insert a knot into the surface along either u or v.

        Parameters
        ----------
        direction : str
            Parameteric direction to insert. Either 'u' or 'v'.
        s : float
            Parametric position along 'direction' to insert
        r : int
            Desired number of times to insert.

        Returns
        -------
        r : int
            The **actual** number of times the knot was inserted.
        """
        if direction not in ["u", "v"]:
            raise Error("Direction must be one of 'u' or 'v'")

        s = checkInput(s, "s", float, 0)
        r = checkInput(r, "r", int, 0)
        if s <= 0.0:
            return
        if s >= 1.0:
            return

        # This is relatively inefficient, but we'll do it for
        # simplicity just call insertknot for each slice in the
        # v-direction:

        if direction == "u":
            # Insert once to know how many times it was actually inserted
            # so we know how big to make the new coef:
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tu, self.ku, self.coef[:, 0].T)

            newCoef = np.zeros((self.nCtlu + actualR, self.nCtlv, self.nDim))

            for j in range(self.nCtlv):
                actualR, tNew, coefSlice, breakPt = libspline.insertknot(s, r, self.tu, self.ku, self.coef[:, j].T)
                newCoef[:, j] = coefSlice[:, 0 : self.nCtlu + actualR].T

            self.tu = tNew[0 : self.nCtlu + self.ku + actualR]
            self.nCtlu = self.nCtlu + actualR

        elif direction == "v":
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tv, self.kv, self.coef[0, :].T)

            newCoef = np.zeros((self.nCtlu, self.nCtlv + actualR, self.nDim))

            for i in range(self.nCtlu):
                actualR, tNew, coefSlice, breakPt = libspline.insertknot(s, r, self.tv, self.kv, self.coef[i, :].T)
                newCoef[i, :] = coefSlice[:, 0 : self.nCtlv + actualR].T

            self.tv = tNew[0 : self.nCtlv + self.kv + actualR]
            self.nCtlv = self.nCtlv + actualR

        self.coef = newCoef

        # break_pt is converted to zero based ordering here!!!
        return actualR, breakPt - 1

    def splitSurface(self, direction, s):
        """
        Split surface into two surfaces at parametric location s

        Parameters
        ----------
        direction : str
            Parameteric direction along which to split. Either 'u' or 'v'.
        s : float
            Parametric position along 'direction' to split

        Returns
        -------
        surf1 : pySpline.surface
            Lower part of the surface

        surf2 : pySpline.surface
            Upper part of the surface
        """
        if direction not in ["u", "v"]:
            raise Error("Direction must be one of 'u' or 'v'")

        # Special case the bounds: (same for both directions)
        if s <= 0.0:
            return None, Surface(tu=self.tu.copy(), tv=self.tv.copy(), ku=self.ku, kv=self.kv, coef=self.coef.copy())
        if s >= 1.0:
            return Surface(tu=self.tu.copy(), tv=self.tv.copy(), ku=self.ku, kv=self.kv, coef=self.coef.copy()), None

        if direction == "u":
            r, breakPt = self.insertKnot(direction, s, self.ku - 1)
            # Break point is now at the right so we need to adjust the
            # counter to the left
            breakPt = breakPt - self.ku + 2

            tt = self.tu[breakPt]
            # Process knot vectors:
            t1 = np.hstack((self.tu[0 : breakPt + self.ku - 1].copy(), tt)) / tt
            t2 = (np.hstack((tt, self.tu[breakPt:].copy())) - tt) / (1.0 - tt)

            coef1 = self.coef[0:breakPt, :, :].copy()
            coef2 = self.coef[breakPt - 1 :, :, :].copy()

            return (
                Surface(tu=t1, tv=self.tv, ku=self.ku, kv=self.kv, coef=coef1),
                Surface(tu=t2, tv=self.tv, ku=self.ku, kv=self.kv, coef=coef2),
            )
        elif direction == "v":

            r, breakPt = self.insertKnot(direction, s, self.kv - 1)
            # Break point is now at the right so we need to adjust the
            # counter to the left
            breakPt = breakPt - self.kv + 2

            tt = self.tv[breakPt]
            # Process knot vectors:
            t1 = np.hstack((self.tv[0 : breakPt + self.kv - 1].copy(), tt)) / tt
            t2 = (np.hstack((tt, self.tv[breakPt:].copy())) - tt) / (1.0 - tt)

            coef1 = self.coef[:, 0:breakPt, :].copy()
            coef2 = self.coef[:, breakPt - 1 :, :].copy()

            return (
                Surface(tu=self.tu, tv=t1, ku=self.ku, kv=self.kv, coef=coef1),
                Surface(tu=self.tu, tv=t2, ku=self.ku, kv=self.kv, coef=coef2),
            )

    def windowSurface(self, uvLow, uvHigh):
        """Create a surface that is windowed by the rectangular
        parametric range defined by uvLow and uvHigh.

        Parameters
        ----------
        uvLow : list or array of length 2
           (u,v) coordinates at the bottom left corner of the parameteric
           box
        uvHigh : list or array of length 2
           (u,v) coordinates at the top left corner of the parameteric
           box

        Returns
        -------
        surf : pySpline.surface
            A new surface defined only on the interior of uvLow -> uvHigh
        """

        # Do u-low split:
        __, surf = self.splitSurface("u", uvLow[0])

        # Do u-high split (and re-normalize the split coordinate)
        surf, __ = surf.splitSurface("u", (uvHigh[0] - uvLow[0]) / (1.0 - uvLow[0]))

        # Do v-low split:
        __, surf = surf.splitSurface("v", uvLow[1])

        # Do v-high split (and re-normalize the split coordinate)
        surf, __ = surf.splitSurface("v", (uvHigh[1] - uvLow[1]) / (1.0 - uvLow[1]))

        return surf

    def getValue(self, u, v):
        """Evaluate the spline surface at parametric positions u,v. This is the
        main function for spline evaluation.

        Parameters
        ----------
        u : float, array or matrix (rank 0, 1, or 2)
            Parametric u values
        v : float, array or matrix (rank 0, 1, or 2)
            Parametric v values

        Returns
        -------
        values : varies
            Spline evaluation at all points u,v. Shape depend on the
            input. If u,v are scalars, values is array of size nDim. If
            u,v are a 1D list, return is (N,nDim) etc.
        """

        u = np.array(u).T
        v = np.array(v).T
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")

        vals = libspline.eval_surface(
            np.atleast_2d(u), np.atleast_2d(v), self.tu, self.tv, self.ku, self.kv, self.coef.T
        )
        return vals.squeeze().T

    def getDerivative(self, u, v):
        """Evaluate the first derivatvies of the spline surface

        Parameters
        ----------
        u : float
            Parametric u value
        v : float
            Parametric v value

        Returns
        -------
        deriv : array size (2,3)
            Spline derivative evaluation at u,vall points u,v. Shape
            depend on the input.
        """
        u = np.array(u)
        v = np.array(v)
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")
        if not np.ndim(u) == 0:
            raise Error("getDerivative only accepts scalar arguments")

        deriv = libspline.eval_surface_deriv(u, v, self.tu, self.tv, self.ku, self.kv, self.coef.T)
        return deriv.T

    def getSecondDerivative(self, u, v):
        """Evaluate the second derivatvies of the spline surface

        deriv = [ (d^2)/(du^2)    (d^2)/(dudv) ]
                [ (d^2)/(dudv)    (d^2)/(dv^2) ]

        Parameters
        ----------
        u : float
            Parametric u value
        v : float
            Parametric v value

        Returns
        -------
        deriv : array size (2,2,3)
            Spline derivative evaluation at u,vall points u,v. Shape
            depend on the input.
        """
        u = np.array(u)
        v = np.array(v)
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")
        if not np.ndim(u) == 0:
            raise Error("getSecondDerivative only accepts scalar arguments")

        deriv = libspline.eval_surface_deriv2(u, v, self.tu, self.tv, self.ku, self.kv, self.coef.T)
        return deriv.T

    def getBounds(self):
        """Determine the extents of the surface

        Returns
        -------
        xMin : array of length 3
            Lower corner of the bounding box
        xMax : array of length 3
            Upper corner of the bounding box
        """
        if self.nDim != 3:
            raise Error("getBounds is only defined for nDim = 3")

        cx = self.coef[:, :, 0].flatten()
        cy = self.coef[:, :, 1].flatten()
        cz = self.coef[:, :, 2].flatten()

        Xmin = np.array([min(cx), min(cy), min(cz)])
        Xmax = np.array([max(cx), max(cy), max(cz)])

        return Xmin, Xmax

    def projectPoint(self, x0, nIter=25, eps=1e-10, **kwargs):
        """
        Perform a point inversion algorithm. Attempt to find the
        closest parameter values (u,v) to the given points x0.

        Parameters
        ----------
        x0 : array
            A point or list of points in nDim space for which the
            minimum distance to the curve is sought.
        nIter : int
            Maximum number of Newton iterations to perform.
        eps : float
            Desired parameter tolerance.
        u : float or array of length x0
            Optional initial guess for u parameter
        v : float or array of length x0
            Optional initial guess for v parameter

        Returns
        -------
        u : float or array
            Solution to the point inversion. u are the u-parametric
            locations yielding the minimum distance to points x0
        v : float or array
            Solution to the point inversion. v are the v-parametric
            locations yielding the minimum distance to points x0
        D : float or array
            Physical distances between the points and the curve.
            This is simply ||surface(u,v) - X0||_2.
        """

        x0 = np.atleast_2d(x0)
        if "u" in kwargs and "v" in kwargs:
            u = np.atleast_1d(kwargs["u"])
            v = np.atleast_1d(kwargs["v"])
        else:
            u = -1 * np.ones(len(x0))
            v = -1 * np.ones(len(x0))

        if not len(x0) == len(u) == len(v):
            raise Error("The length of x0 and u, v must be the same")

        # If necessary get brute-force starting point
        if np.any(u < 0) or np.any(u > 1) or np.any(v < 0):
            self.computeData()
            u, v = libspline.point_surface_start(x0.T, self.udata, self.vdata, self.data.T)

        D = np.zeros_like(x0)
        for i in range(len(x0)):
            u[i], v[i], D[i] = libspline.point_surface(
                x0[i], self.tu, self.tv, self.ku, self.kv, self.coef.T, nIter, eps, u[i], v[i]
            )

        return u.squeeze(), v.squeeze(), D.squeeze()

    def projectCurve(self, inCurve, nIter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this surface (self) and a curve (inCurve).

        Parameters
        ----------
        inCurve : pySpline.curve objet
           Curve to use
        nIter : int
            Maximum number of Newton iterations to perform.
        eps : float
            Desired parameter tolerance.
        s : float
            Initial solution guess for curve
        u : float
            Initial solution guess for parametric u position
        v : float
            Initial solution guess for parametric v position

        Returns
        -------
        u : float
            Surface parameter u yielding min distance to point x0
        u : float
            Surface parameter v yielding min distance to point x0
        s : float
            Parametric position on curve  yielding min distance to point x0
        D : float
            Minimum distance between this surface and curve.
            is equivalent to ||surface(u,v) - curve(s)||_2.
        """
        u = -1.0
        v = -1.0
        s = -1.0
        if "u" in kwargs:
            u = checkInput(kwargs["u"], "u", float, 0)
        if "v" in kwargs:
            v = checkInput(kwargs["v"], "v", float, 0)
        if "s" in kwargs:
            s = checkInput(kwargs["s"], "s", float, 0)
        nIter = checkInput(nIter, "nIter", int, 0)
        eps = checkInput(eps, "eps", float, 0)

        # If necessary get brute-force starting point
        if np.any(u < 0) or np.any(u > 1) or np.any(v < 0):
            self.computeData()
            inCurve.computeData()
            s, u, v = libspline.curve_surface_start(inCurve.data.T, inCurve.sdata, self.data.T, self.udata, self.vdata)

        return libspline.curve_surface(
            inCurve.t, inCurve.k, inCurve.coef.T, self.tu, self.tv, self.ku, self.kv, self.coef.T, nIter, eps, u, v, s
        )

    def computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """

        # We will base the data on interpolated greville points
        if self.data is None:
            self.edgeCurves[0].calcInterpolatedGrevillePoints()
            self.udata = self.edgeCurves[0].sdata
            self.edgeCurves[2].calcInterpolatedGrevillePoints()
            self.vdata = self.edgeCurves[2].sdata
            [V, U] = np.meshgrid(self.vdata, self.udata)
            self.data = self.getValue(U, V)

    def writeDirections(self, handle, isurf):
        """Write out and indication of the surface direction"""
        if self.nCtlu >= 3 and self.nCtlv >= 3:
            data = np.zeros((4, self.nDim))
            data[0] = self.coef[1, 2]
            data[1] = self.coef[1, 1]
            data[2] = self.coef[2, 1]
            data[3] = self.coef[3, 1]
            writeTecplot1D(handle, "surface%d direction" % (isurf), data)
        else:
            print("Not Enough control points to output direction indicator")

    def writeTecplot(self, fileName, surf=True, coef=True, orig=True, directions=False):
        """
        Write the surface to a tecplot .dat file

        Parameters
        ----------
        fileName : str
            File name for tecplot file. Should have .dat extension
        surf : bool
            Flag to write discrete approximation of the actual surface
        coef: bool
            Flag to write b-spline coefficients
        orig : bool
            Flag to write original data (used for fitting) if it exists
        directions : bool
            Flag to write surface direction visualization
        """
        f = openTecplot(fileName, self.nDim)
        if surf:
            self.computeData()
            writeTecplot2D(f, "interpolated", self.data)
        if coef:
            writeTecplot2D(f, "control_pts", self.coef)
        if orig and self.origData:
            writeTecplot2D(f, "orig_data", self.X)
        if directions:
            self.writeDirections(f, 0)
        closeTecplot(f)

    def writeIGES_directory(self, handle, Dcount, Pcount):
        """
        Write the IGES file header information (Directory Entry Section)
        for this surface
        """
        # A simpler calc based on cmlib definitions The 13 is for the
        # 9 parameters at the start, and 4 at the end. See the IGES
        # 5.3 Manual paraEntries = 13 + Knotsu + Knotsv + Weights +
        # control points
        if self.nDim != 3:
            raise Error("Must have 3 dimensions to write to IGES file")
        paraEntries = 13 + (len(self.tu)) + (len(self.tv)) + self.nCtlu * self.nCtlv + 3 * self.nCtlu * self.nCtlv + 1

        paraLines = (paraEntries - 10) // 3 + 2
        if np.mod(paraEntries - 10, 3) != 0:
            paraLines += 1

        handle.write("     128%8d       0       0       1       0       0       000000001D%7d\n" % (Pcount, Dcount))
        handle.write(
            "     128       0       2%8d       0                               0D%7d\n" % (paraLines, Dcount + 1)
        )
        Dcount += 2
        Pcount += paraLines

        return Pcount, Dcount

    def writeIGES_parameters(self, handle, Pcount, counter):
        """
        Write the IGES parameter information for this surface
        """
        handle.write(
            "%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n"
            % (128, self.nCtlu - 1, self.nCtlv - 1, self.ku - 1, self.kv - 1, Pcount, counter)
        )
        counter += 1
        handle.write("%10d,%10d,%10d,%10d,%10d,          %7dP%7d\n" % (0, 0, 1, 0, 0, Pcount, counter))

        counter += 1
        pos_counter = 0

        for i in range(len(self.tu)):
            pos_counter += 1
            handle.write("%20.12g," % (np.real(self.tu[i])))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in range(len(self.tv)):
            pos_counter += 1
            handle.write("%20.12g," % (np.real(self.tv[i])))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for _i in range(self.nCtlu * self.nCtlv):
            pos_counter += 1
            handle.write("%20.12g," % (1.0))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for j in range(self.nCtlv):
            for i in range(self.nCtlu):
                for idim in range(3):
                    pos_counter += 1
                    handle.write("%20.12g," % (np.real(self.coef[i, j, idim])))
                    if np.mod(pos_counter, 3) == 0:
                        handle.write("  %7dP%7d\n" % (Pcount, counter))
                        counter += 1
                        pos_counter = 0
                    # end if
                # end for
            # end for
        # end for

        # Ouput the ranges
        for i in range(4):
            pos_counter += 1
            if i == 0:
                handle.write("%20.12g," % (np.real(self.umin)))
            if i == 1:
                handle.write("%20.12g," % (np.real(self.umax)))
            if i == 2:
                handle.write("%20.12g," % (np.real(self.vmin)))
            if i == 3:
                # semi-colon for the last entity
                handle.write("%20.12g;" % (np.real(self.vmax)))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            else:  # We have to close it up anyway
                if i == 3:
                    for _j in range(3 - pos_counter):
                        handle.write("%21s" % (" "))
                    # end for
                    pos_counter = 0
                    handle.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                # end if
            # end if
        # end for

        Pcount += 2

        return Pcount, counter

    def writeTin(self, handle):
        """Write the pySpline surface to an open handle in .tin format"""
        handle.write("bspline\n")

        # Sizes and Order
        handle.write("%d, %d, %d, %d, 0\n" % (self.nCtlu, self.nCtlv, self.ku, self.kv))

        # U - Knot Vector
        for i in range(len(self.tu)):
            handle.write("%16.12g, \n" % (self.tu[i]))

        # V - Knot Vector
        for j in range(len(self.tv)):
            handle.write("%16.12g, \n" % (self.tv[j]))

        # Control points:
        for j in range(self.nCtlv):
            for i in range(self.nCtlu):
                handle.write(
                    "%16.12g, %16.12g, %16.12g\n" % (self.coef[i, j, 0], self.coef[i, j, 1], self.coef[i, j, 2])
                )
