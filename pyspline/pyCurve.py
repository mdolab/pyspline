# External modules
import numpy as np
from scipy.sparse import linalg

# Local modules
from . import libspline
from .utils import Error, _assembleMatrix, checkInput, closeTecplot, openTecplot, writeTecplot1D


class Curve(object):
    """
    Create an instance of a b-spline curve. There are two
    ways to initialize the class

    * **Creation**: Create an instance of the spline class
      directly by supplying the required information. kwargs MUST
      contain the following information: ``k, t, coef``.


    * **LMS/Interpolation**: Create an instance of the spline class by
      using an interpolating spline to given data points or a LMS
      spline. The following keyword argument information is required:

      1. ``k`` Spline Order

      2. ``X`` real array size (N, nDim) of data to fit. **OR**
          1. ``x`` (1D) and ``s`` for 1D
          2. ``x`` (1D) and ``y`` (1D) for 2D spatial curve
          3. ``x`` (1D) and ``y``` (1D) and ``z`` 1D for 3D spatial curve

      3. ``s`` real array of size (N). Optional for nDim >= 2

    Parameters
    ----------
    k : {2, 3, 4}
        Order for spline. A spline with order :math:`n` has at most :math:`C^{n-2}` continuity. The :math:`n-1` derivative will be piecewise constant.
    nCtl : int
        Number of control points. If this is specified then LMS will be used instead of interpolation.
    t : array, list
        Knot vector. Used only for creation. Must be size ``nCtl + k``
    coef : array, size (nCtl, nDim)
        Coefficients to use. Only used for creation. The second
        dimension determine the spatial order of the spline
    X : arary, size (N, ndim)
        Full array of data to interpolate/fit
    x : array, size (N)
        Just x coordinates of data to fit
    y : array, size (N)
        Just y coordinates of data to fit
    z : array, size (N)
        Just z coordinates of data to fit
    s : array, size(N)
        Optional parameter to use. Not required for nDim >=2
    nIter : int
        The number of Hoscheks parameter correction iterations to
        run. Only used for LMS fitting.
    weight : array, size(N)
        A array of weighting factors for each fitting point. A value
        of -1 can be used to exactly constrain a point. By default,
        all weights are 1.0
    deriv : array, size (len(derivPtr), nDim)
        Array containing derivative information the
        user wants to use at a particular point. **EXPERIMENTAL**
    derivPtr : int, array, size (nDerivPtr)
        Array of indices pointing to the index of points in X (or
        x,y,z), for which the user has supplied a derivative for in
        ``deriv``. **EXPERIMENTAL**
    derivWeights : array size(nDerivPtr)
        Optional array of the weighting to use for
        derivatives. A value of -1 can be used to exactly constrain a
        derivative. **EXPERIMENTAL**

    Examples
    --------
    >>> x = [0, 0.5, 1.0]
    >>> y = [0, 0.25, 1.0]
    >>> s = [0., 0.5, 1.0]
    >>> # Spatial interpolated seg (k=2 makes this a straight line)
    >>> line_seg = Curve(x=x, y=y, k=2)
    >>> # With explicit parameter values
    >>> line_seg = Curve(x=x, y=y, k=2, s=s)
    >>> #Interpolate parabolic curve
    >>> parabola = Curve(x=x, y=y, k=3)
    >>> #Interpolate parabolic curve with parameter values
    >>> parabola = Curve(x=x, y=y, k=3, s=s)
    """

    def __init__(self, **kwargs):
        self.length = None
        self.gpts = None
        self.data = None
        self.sdata = None
        self.localInterp = False
        # We have provided information to create curve directly
        if "k" in kwargs and "t" in kwargs and "coef" in kwargs:
            self.s = None
            self.X = None
            self.N = None
            self.k = checkInput(kwargs["k"], "k", int, 0)
            self.coef = np.atleast_2d(kwargs["coef"])
            self.nCtl = self.coef.shape[0]
            self.nDim = self.coef.shape[1]
            self.t = checkInput(kwargs["t"], "t", float, 1, (self.nCtl + self.k,))
            self.origData = False
            self.calcGrevillePoints()

        elif "localInterp" in kwargs:
            # Local, non-global interpolation. We could use this for
            # second-order (linear interpolation) but there is not
            # much point, since it would be the same as 'interp'
            # below.
            self.localInterp = True
            self.k = 4
            self.origData = True
            if "X" in kwargs:
                self.X = np.atleast_2d(kwargs["X"])
                if np.ndim(kwargs["X"]) == 1:
                    self.X = np.transpose(self.X)
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = np.vstack([kwargs["x"], kwargs["y"], kwargs["z"]]).T
            elif "x" in kwargs and "y" in kwargs:
                self.X = np.vstack([kwargs["x"], kwargs["y"]]).T
            elif "x" in kwargs:
                self.X = np.transpose(np.atleast_2d(kwargs["x"]))
            # enf if

            self.nDim = self.X.shape[1]
            self.N = self.X.shape[0]

            if "s" in kwargs:
                self.s = checkInput(kwargs["s"], "s", float, 1, (self.N,))
                self.length = 1.0  # Assume it is 0->1
            else:
                if self.nDim <= 1:
                    raise Error("For 1D splines, the basis, 's' must be given.")
                self._getParameterization()

            # Now we have the data we need to generate the local
            # interpolation

            T = np.zeros((self.N, self.nDim))
            # Compute tangents
            qq = np.zeros_like(self.X)
            T = np.zeros_like(self.X)
            deltaS = np.zeros(self.N)
            for i in range(1, self.N):
                deltaS[i] = self.s[i] - self.s[i - 1]
                qq[i, :] = self.X[i] - self.X[i - 1]

            for i in range(1, self.N - 1):
                a = deltaS[i] / (deltaS[i] + deltaS[i + 1])
                T[i] = (1 - a) * qq[i] + a * qq[i + 1]

            # Do the start and end points: (eqn: 9.32, The NURBS book)
            T[0] = 2 * qq[1] / deltaS[1] - T[1]
            T[-1] = 2 * qq[-1] / deltaS[-1] - T[-2]
            # Normalize
            for i in range(self.N):
                T[i] = T[i] / np.linalg.norm(T[i])

            # Final coefficients and t
            self.coef = np.zeros((2 * (self.N - 1) + 2, self.nDim))
            self.t = np.zeros(len(self.coef) + self.k)
            self.nCtl = len(self.coef)
            # End Pts
            self.coef[0] = self.X[0].copy()
            self.coef[-1] = self.X[-1].copy()

            # Interior coefficients
            for i in range(self.N - 1):
                a = self.length * (self.s[i + 1] - self.s[i])
                self.coef[2 * i + 1] = self.X[i] + a / 3.0 * T[i]
                self.coef[2 * i + 2] = self.X[i + 1] - a / 3.0 * T[i + 1]

            # Knots
            self.t[-4:] = 1.0
            u = np.zeros(self.N)
            for i in range(0, self.N - 1):
                u[i + 1] = u[i] + np.linalg.norm(self.coef[2 * i + 2] - self.coef[2 * i + 1])

            for i in range(1, self.N - 1):
                self.t[2 * i + 2] = u[i] / u[self.N - 1]
                self.t[2 * i + 3] = u[i] / u[self.N - 1]

        else:  # lms or interpolate function
            if not ("k" in kwargs and ("X" in kwargs or "x" in kwargs)):
                raise ValueError(
                    "At least spline order, k and X (or x=, y=) MUST be defined for (interpolation) spline creation. "
                    + "nCtl=<number of control points> must be specified for a LMS fit"
                )
            self.origData = True
            if "X" in kwargs:
                self.X = np.atleast_2d(kwargs["X"])
                if np.ndim(kwargs["X"]) == 1:
                    self.X = np.transpose(self.X)
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = np.vstack([kwargs["x"], kwargs["y"], kwargs["z"]]).T
            elif "x" in kwargs and "y" in kwargs:
                self.X = np.vstack([kwargs["x"], kwargs["y"]]).T
            elif "x" in kwargs:
                self.X = np.transpose(np.atleast_2d(kwargs["x"]))
            # enf if

            self.nDim = self.X.shape[1]
            self.N = self.X.shape[0]

            self.k = checkInput(kwargs["k"], "k", int, 0)
            self.t = None
            if "nIter" in kwargs:
                self.nIter = checkInput(kwargs["nIter"], "nIter", int, 0)
            else:
                self.nIter = 1

            if "s" in kwargs:
                self.s = checkInput(kwargs["s"], "s", float, 1, (self.N,))
            else:
                if self.nDim <= 1:
                    raise Error("For 1D splines, the basis, 's' must be given.")
                self._getParameterization()

            if "weights" in kwargs:
                self.weights = checkInput(kwargs["weights"], "weights", float, 1, (self.N,))
            else:
                self.weights = np.ones(self.N)

            if "deriv" in kwargs and "derivPtr" in kwargs:
                self.deriv = checkInput(kwargs["deriv"], "deriv", float, 2)
                self.derivPtr = checkInput(kwargs["derivPtr"], "derivPtr", int, 1, (len(self.deriv),))
            else:
                self.deriv = None
                self.derivPtr = np.array([])

            if "derivWeights" in kwargs and self.deriv is not None:
                self.derivWeights = checkInput(kwargs["derivWeights"], "derivWeights", float, 1, (len(self.derivPtr),))
            else:
                if self.deriv is not None:
                    self.derivWeights = np.ones(len(self.deriv))
                else:
                    self.derivWeights = None

            # We are doing an interpolation...set all weights to 1
            if "nCtl" not in kwargs:
                self.interp = True
                self.weights[:] = 1
                if self.derivWeights is not None:
                    self.derivWeights[:] = 1
            else:
                self.interp = False
                self.nCtl = checkInput(kwargs["nCtl"], "nCtl", int, 0)

            self.recompute(self.nIter, computeKnots=True)

    def recompute(self, nIter, computeKnots=True):
        """
        Run iterations of Hoscheks Parameter Correction on the current curve

        Parameters
        ----------
        nIter : int
            The number of parameter correction iterations to run
        computeKnots : bool
            Flag whether or not the knots should be recomputed
        """

        # Return if we don't have original data to fit
        if not self.origData or self.localInterp:
            return

        # Do the separation between the constrained and unconstrained:
        # u -> unconstrained
        # s -> constrained
        suSelect = np.where(self.weights > 0.0)
        scSelect = np.where(self.weights <= 0.0)
        S = self.X[suSelect]
        su = self.s[suSelect]
        T = self.X[scSelect]
        sc = self.s[scSelect]
        weights = self.weights[np.where(self.weights > 0.0)]

        nu = len(S)
        nc = len(T)

        # And the derivative info
        if self.deriv is not None:
            sduSelect = np.where(self.derivWeights > 0.0)
            sdcSelect = np.where(self.derivWeights <= 0.0)
            S = np.vstack((S, self.deriv[sduSelect]))
            sdu = self.s[self.derivPtr][sduSelect]
            T = np.vstack((T, self.deriv[sdcSelect]))
            sdc = self.s[self.derivPtr][sdcSelect]
            weights = np.append(weights, self.derivWeights[np.where(self.derivWeights > 0.0)])
            ndu = len(sdu)
            ndc = len(sdc)
        else:
            sdu = np.array([], "d")
            sdc = np.array([], "d")
            ndu = 0
            ndc = 0

        W = _assembleMatrix(weights, np.arange(len(weights)), np.arange(len(weights) + 1), (len(weights), len(weights)))

        if self.interp:
            self.nCtl = nu + nc + ndu + ndc
            self.nIter = 1

        # Sanity check to make sure k is ok
        if nu + nc + ndu + ndc < self.k:
            self.k = nu + nc + ndu + ndc

        if computeKnots:
            # Generate the knot vector, if necessary greville points and
            # empty coefficients
            if self.interp:
                self.t = libspline.knots_interp(self.s, self.derivPtr, self.k)
            else:
                self.t = libspline.knots_lms(self.s, self.nCtl, self.k)

        self.calcGrevillePoints()
        self.coef = np.zeros((self.nCtl, self.nDim), "d")

        # Get the 'N' jacobian
        nVals = np.zeros((nu + ndu) * self.k)  # |
        nRowPtr = np.zeros(nu + ndu + 1, "intc")  # | -> CSR formulation
        nColInd = np.zeros((nu + ndu) * self.k, "intc")  # |
        libspline.curve_jacobian_wrap(su, sdu, self.t, self.k, self.nCtl, nVals, nRowPtr, nColInd)
        N = _assembleMatrix(nVals, nColInd, nRowPtr, (nu + ndu, self.nCtl)).tocsc()

        if self.interp:
            # Factorize once for efficiency
            solve = linalg.dsolve.factorized(N)
            for idim in range(self.nDim):
                self.coef[:, idim] = solve(S[:, idim])

            return

        # If we do NOT have an interpolation:
        length = libspline.poly_length(self.X.T)
        for _i in range(nIter):
            su = self.s[suSelect]
            sc = self.s[scSelect]
            if self.deriv is not None:
                sdu = self.s[self.derivPtr][sduSelect]
                sdc = self.s[self.derivPtr][sdcSelect]

            libspline.curve_jacobian_wrap(su, sdu, self.t, self.k, self.nCtl, nVals, nRowPtr, nColInd)
            NTWN = (N.transpose() * W * N).tocsc()  # We need this either way

            if nc + ndc == 0:  # We are doing LMS but no
                # constraints...just a straight weighted
                # LMS

                # Factorize once for efficiency
                solve = linalg.dsolve.factorized(NTWN)
                for idim in range(self.nDim):
                    self.coef[:, idim] = solve(N.transpose() * W * S[:, idim])

            else:
                # Now its more complicated since we have
                # constraints --only works with scipy Sparse
                # matrices

                mVals = np.zeros((nc + ndc) * self.k)  # |
                mRowPtr = np.zeros(nc + ndc + 1, "intc")  # | -> CSR
                mColInd = np.zeros((nc + ndc) * self.k, "intc")  # |

                libspline.curve_jacobian_wrap(sc, sdc, self.t, self.k, self.nCtl, mVals, mRowPtr, mColInd)
                M = _assembleMatrix(mVals, mColInd, mRowPtr, (nc + ndc, self.nCtl))

                # Now we must assemble the constrained jacobian
                # [ N^T*W*T      M^T][P] = [ N^T*W*S]
                # [ M            0  ][R]   [ T      ]

                MT = M.transpose().tocsr()

                jVal, jColInd, jRowPtr = libspline.constr_jac(
                    NTWN.data,
                    NTWN.indptr,
                    NTWN.indices,
                    MT.data,
                    MT.indptr,
                    MT.indices,
                    M.data,
                    M.indptr,
                    M.indices,
                    self.nCtl,
                )

                # Create sparse csr matrix and factorize
                J = _assembleMatrix(jVal, jColInd, jRowPtr, (self.nCtl + nc + ndc, self.nCtl + nc + ndc))

                # Factorize once for efficiency
                solve = linalg.dsolve.factorized(J)
                for idim in range(self.nDim):
                    rhs = np.hstack((N.transpose() * W * S[:, idim], T[:, idim]))
                    self.coef[:, idim] = solve(rhs)[0 : self.nCtl]

            # end if (constr - not constrained

            # Run para correction
            libspline.curve_para_corr(self.t, self.k, self.s, self.coef.T, length, self.X.T)
        # end for (iter loop)

        # Check the RMS
        rms = 0.0
        for idim in range(self.nDim):
            rms += np.linalg.norm(N * self.coef[:, idim] - S[:, idim]) ** 2

        rms = np.sqrt(rms / self.N)

    def _getParameterization(self):
        """Compute a parametrization for the curve based on an
        arc-length formulation
        """
        self.s = np.zeros(self.N, "d")
        for i in range(self.N - 1):
            dist = 0
            for idim in range(self.nDim):
                dist += (self.X[i + 1, idim] - self.X[i, idim]) ** 2
            self.s[i + 1] = self.s[i] + np.sqrt(dist)

        self.length = self.s[-1]
        self.s = self.s / self.s[-1]

    def reverse(self):
        """
        Reverse the direction of this curve
        """
        self.coef = self.coef[::-1, :]
        self.t = 1 - self.t[::-1]

    def insertKnot(self, u, r):
        """
        Insert at knot in the curve at parametric position u

        Parameters
        ----------
        u : float
            Parametric position to split at
        r : int
            Number of time to insert. Should be > 0.

        Returns
        -------
        actualR : int
            The number of times the knot was **actually** inserted
        breakPt : int
            Index in the knot vector of the new knot(s)
        """

        u = checkInput(u, "u", float, 0)
        r = checkInput(r, "r", int, 0)

        if u <= 0:
            return
        if u >= 1.0:
            return

        actualR, tNew, coefNew, breakPt = libspline.insertknot(u, r, self.t, self.k, self.coef.T)
        self.t = tNew[0 : self.nCtl + self.k + actualR]
        self.coef = coefNew[:, 0 : self.nCtl + actualR].T
        self.nCtl = self.nCtl + actualR

        # break_pt is converted to zero based ordering here!!!
        return actualR, breakPt - 1

    def splitCurve(self, u):
        """
        Split the curve at parametric position u. This uses the
        :func:`insertKnot` routine

        Parameters
        ----------
        u : float
            Parametric position to insert knot.

        Returns
        -------
        curve1 : pySpline curve object
            Curve from s=[0, u]

        curve2 : pySpline curve object
            Curve from s=[u, 1]

        Notes
        -----
        curve1 and curve2 may be None if the parameter u is outside
        the range of (0, 1)
        """
        u = checkInput(u, "u", float, 0)

        if u <= 0.0:
            return None, Curve(t=self.t.copy(), k=self.k, coef=self.coef.copy())
        if u >= 1.0:
            return Curve(t=self.t.copy(), k=self.k, coef=self.coef.copy()), None

        r, breakPt = self.insertKnot(u, self.k - 1)
        # Break point is now at the right so we need to adjust the
        # counter to the left
        breakPt = breakPt - self.k + 2

        # Process knot vectors:
        uu = self.t[breakPt]
        t1 = np.hstack((self.t[0 : breakPt + self.k - 1].copy(), uu)) / uu
        t2 = (np.hstack((uu, self.t[breakPt:].copy())) - uu) / (1.0 - uu)

        coef1 = self.coef[0:breakPt, :].copy()
        coef2 = self.coef[breakPt - 1 :, :].copy()

        return Curve(t=t1, k=self.k, coef=coef1), Curve(t=t2, k=self.k, coef=coef2)

    def windowCurve(self, uLow, uHigh):
        """
        Compute the segment of the curve between the two parameter values

        Parameters
        ----------
        uLow : float
            Lower bound for the clip
        uHigh : float
            Upper bound for the clip
        """
        __, c = self.splitCurve(uLow)
        c, __ = c.splitCurve((uHigh - uLow) / (1.0 - uLow))

        return c

    def getLength(self):
        """Compute the length of the curve using the Euclidean Norm

        Returns
        -------
        length : float
           Physical length of curve
        """
        points = self.getValue(self.gpts)
        length = 0
        for i in range(len(points) - 1):
            length += np.linalg.norm(points[i] - points[i + 1])

        return length

    def calcGrevillePoints(self):
        """Calculate the Greville points"""
        self.gpts = np.zeros(self.nCtl)
        for i in range(self.nCtl):
            for n in range(self.k - 1):  # degree loop
                self.gpts[i] += self.t[i + n + 1]
            self.gpts[i] = self.gpts[i] / (self.k - 1)

    def calcInterpolatedGrevillePoints(self):
        """
        Compute greville points, but with additional interpolated knots
        """
        self.calcGrevillePoints()
        s = [self.gpts[0]]
        N = 2
        for i in range(len(self.gpts) - 1):
            for j in range(N):
                s.append((self.gpts[i + 1] - self.gpts[i]) * (j + 1) / (N + 1) + self.gpts[i])
            s.append(self.gpts[i + 1])

        self.sdata = np.array(s)

    def __call__(self, s):
        """
        Equivalent to getValue()
        """
        return self.getValue(s)

    def getValue(self, s):
        """
        Evaluate the spline at parametric position, s

        Parameters
        ----------
        s : float or array
            Parametric position(s) at which to evaluate the curve.

        Returns
        -------
        values : array of size nDim or an array of size (N, 3)
            The evaluated points. If a single scalar s was given, the
            result with be an array of length nDim (or a scalar if
            nDim=1). If a vector of s values were given it will be an
            array of size (N, 3) (or size (N) if ndim=1)
        """

        s = np.array(s).T
        if self.coef.dtype == np.dtype("d"):
            vals = libspline.eval_curve(np.atleast_1d(s), self.t, self.k, self.coef.T)
        else:
            vals = libspline.eval_curve_c(np.atleast_1d(s).astype("D"), self.t, self.k, self.coef.T)

        return vals.squeeze().T

    def getDerivative(self, s):
        """
        Evaluate the derivatie of the spline at parametric position, s

        Parameters
        ----------
        s : float
            Parametric location for derivative

        Returns
        -------
        ds : array
            The first derivative. This is an array of size nDim
        """
        if self.coef.dtype == np.dtype("d"):
            ds = libspline.eval_curve_deriv(s, self.t, self.k, self.coef.T).squeeze()
        else:
            ds = libspline.eval_curve_deriv_c(s, self.t, self.k, self.coef.T).squeeze()

        return ds

    def getSecondDerivative(self, s):
        """
        Evaluate the second derivatie of the spline at parametric
        position, s

        Parameters
        ----------
        s : float
            Parametric location for derivative

        Returns
        -------
        d2s : array
            The second derivative. This is an array of size nDim

        """
        if self.coef.dtype == np.dtype("d"):
            d2s = libspline.eval_curve_deriv2(s, self.t, self.k, self.coef.T).squeeze()
        else:
            d2s = libspline.eval_curve_deriv2_c(s, self.t, self.k, self.coef.T).squeeze()

        return d2s

    def projectPoint(self, x0, nIter=25, eps=1e-10, **kwargs):
        """
        Perform a point inversion algorithm. Attempt to find the
        closest parameter values to the given points x0.

        Parameters
        ----------
        x0 : array
            A point or list of points in nDim space for which the
            minimum distance to the curve is sought.
        nIter : int
            Maximum number of Newton iterations to perform.
        eps : float
            Desired parameter tolerance.

        Returns
        -------
        s : float or array
            Solution to the point inversion. s are the parametric
            locations yielding the minimum distance to points x0
        D : float or array
            Physical distances between the points and the curve.
            This is simply ||curve(s) - X0||_2.
        """
        x0 = np.atleast_2d(x0)
        if "s" in kwargs:
            s = np.atleast_1d(kwargs["s"])
        else:
            s = -1 * np.ones(len(x0))

        if len(x0) != len(s):
            raise Error("projectPoint: The length of x0 and s must be the same")

        # If necessary get brute-force starting point
        if np.any(s < 0) or np.any(s > 1):
            self.computeData()
            s = libspline.point_curve_start(x0.T, self.sdata, self.data.T)

        D = np.zeros_like(x0)
        for i in range(len(x0)):
            s[i], D[i] = libspline.point_curve(x0[i], self.t, self.k, self.coef.T, nIter, eps, s[i])
        return s.squeeze(), D.squeeze()

    def projectCurve(self, inCurve, nIter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this curve (self) and a
        second curve passed in (inCurve)

        Parameters
        ----------
        inCurve : pySpline.curve objet
            Other curve to use
        nIter : int
            Maximum number of Newton iterations to perform.
        eps : float
            Desired parameter tolerance.
        s : float
            Initial guess for curve1 (this curve class)
        t : float
            Initial guess for inCurve (curve passed in )

        Returns
        -------
        s : float
           Parametric position on curve1 (this class)
        t : float
            Parametric position on curve2 (inCurve)
        D : float
            Minimum distance between this curve and inCurve. It
            is equivalent to ||self(s) - inCurve(t)||_2.
        """
        s = -1
        t = -1
        if "s" in kwargs:
            s = checkInput(kwargs["s"], "s", float, 0)
        if "t" in kwargs:
            t = checkInput(kwargs["t"], "t", float, 0)
        eps = checkInput(eps, "eps", float, 0)

        if s < 0 or s > 1 or t < 0 or t > 1:
            self.computeData()
            inCurve.computeData()
            s, t = libspline.curve_curve_start(self.data.T, self.sdata, inCurve.data.T, inCurve.sdata)

        s, t, Diff = libspline.curve_curve(
            self.t, self.k, self.coef.T, inCurve.t, inCurve.k, inCurve.coef.T, nIter, eps, s, t
        )

        return s, t, Diff

    def projectCurveMultiSol(self, inCurve, nIter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this curve (self) and a
        second curve passed in (inCurve). Try to find as many local
        solutions as possible

        Parameters
        ----------
        inCurve : pySpline.curve objet
            Other curve to use
        nIter : int
            Maximum number of Newton iterations to perform.
        eps : float
            Desired parameter tolerance.
        s : float
            Initial guess for curve1 (this curve class)
        t : float
            Initial guess for inCurve (curve passed in )

        Returns
        -------
        s : array
           Parametric position(s) on curve1 (this class)
        t : float
            Parametric position)s_ on curve2 (inCurve)
        D : float
            Minimum distance(s) between this curve and inCurve. It
            is equivalent to ||self(s) - inCurve(t)||_2.
        """
        s = -1
        t = -1
        if "s" in kwargs:
            s = checkInput(kwargs["s"], "s", float, 0)
        if "t" in kwargs:
            t = checkInput(kwargs["t"], "t", float, 0)
        eps = checkInput(eps, "eps", float, 0)

        self.computeData()
        inCurve.computeData()

        uSol = []
        tSol = []
        diff = []
        for i in range(len(self.sdata)):
            for j in range(len(inCurve.sdata)):
                s, t, Diff = libspline.curve_curve(
                    self.t,
                    self.k,
                    self.coef.T,
                    inCurve.t,
                    inCurve.k,
                    inCurve.coef.T,
                    nIter,
                    eps,
                    self.sdata[i],
                    inCurve.sdata[j],
                )

                if np.linalg.norm(Diff) < eps:
                    # Its a solution. Check it it is already in list:
                    if len(uSol) == 0:
                        uSol.append(s)
                        tSol.append(t)
                        diff.append(Diff)

                    for ii in range(len(uSol)):
                        if abs(uSol[ii] - s) < eps and abs(tSol[ii] - t) < eps:
                            pass
                        else:
                            uSol.append(s)
                            tSol.append(t)
                            diff.append(Diff)

        return np.array(uSol), np.array(tSol), np.array(diff)

    def computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """
        # We will base the data on interpolated greville points

        if self.data is None:
            self.calcInterpolatedGrevillePoints()
            self.data = self.getValue(self.sdata)

    def writeTecplot(self, fileName, curve=True, coef=True, orig=True):
        """
        Write the curve to a tecplot .dat file

        Parameters
        ----------
        fileName : str
            File name for tecplot file. Should have .dat extension
        curve : bool
            Flag to write discrete approximation of the actual curve
        coef : bool
            Flag to write b-spline coefficients
        orig : bool
            Flag to write original data (used for fitting) if it exists
        """
        f = openTecplot(fileName, self.nDim)
        if curve:
            self.computeData()
            writeTecplot1D(f, "interpolated", self.data)
        if coef:
            writeTecplot1D(f, "control_pts", self.coef)
        if orig and self.origData:
            writeTecplot1D(f, "orig_data", self.X)

        closeTecplot(f)

    def writeIGES_directory(self, handle, Dcount, Pcount, twoD=False):
        """
        Write the IGES file header information (Directory Entry Section)
        for this curve.

        DO NOT MODIFY ANYTHING HERE UNLESS YOU KNOW **EXACTLY** WHAT
        YOU ARE DOING!

        """

        if self.nDim != 3:
            raise Error("Must have 3 dimensions to write to IGES file")
        paraEntries = 6 + len(self.t) + self.nCtl + 3 * self.nCtl + 5

        paraLines = (paraEntries - 11) // 3 + 2
        if np.mod(paraEntries - 11, 3) != 0:
            paraLines += 1
        if twoD:
            handle.write("     126%8d       0       0       1       0       0       001010501D%7d\n" % (Pcount, Dcount))
            handle.write(
                "     126       0       2%8d       0                               0D%7d\n" % (paraLines, Dcount + 1)
            )
        else:
            handle.write("     126%8d       0       0       1       0       0       000000001D%7d\n" % (Pcount, Dcount))
            handle.write(
                "     126       0       2%8d       0                               0D%7d\n" % (paraLines, Dcount + 1)
            )

        Dcount += 2
        Pcount += paraLines

        return Pcount, Dcount

    def writeIGES_parameters(self, handle, Pcount, counter):
        """Write the IGES parameter information for this curve.

        DO NOT MODIFY ANYTHING HERE UNLESS YOU KNOW **EXACTLY** WHAT
        YOU ARE DOING!
        """
        handle.write(
            "%10d,%10d,%10d,0,0,0,0,                        %7dP%7d\n"
            % (126, self.nCtl - 1, self.k - 1, Pcount, counter)
        )
        counter += 1
        pos_counter = 0

        for i in range(len(self.t)):
            pos_counter += 1
            handle.write("%20.12g," % (np.real(self.t[i])))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0

        for _i in range(self.nCtl):
            pos_counter += 1
            handle.write("%20.12g," % (1.0))
            if np.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0

        for i in range(self.nCtl):
            for idim in range(3):
                pos_counter += 1
                handle.write("%20.12g," % (np.real(self.coef[i, idim])))
                if np.mod(pos_counter, 3) == 0:
                    handle.write("  %7dP%7d\n" % (Pcount, counter))
                    counter += 1
                    pos_counter = 0
        if pos_counter == 1:
            handle.write("%s    %7dP%7d\n" % (" " * 40, Pcount, counter))
            counter += 1
        elif pos_counter == 2:
            handle.write("%s    %7dP%7d\n" % (" " * 20, Pcount, counter))
            counter += 1

        # Ouput the ranges
        handle.write("%20.12g,%20.12g,0.0,0.0,0.0;         " % (np.min(self.t), np.max(self.t)))
        handle.write("  %7dP%7d\n" % (Pcount, counter))
        counter += 1
        Pcount += 2

        return Pcount, counter
