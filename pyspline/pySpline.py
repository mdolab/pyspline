"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# ===========================================================================
# External Python modules
# ===========================================================================
import warnings
import numpy
from scipy import sparse
from scipy.sparse import linalg

# ===========================================================================
# Custom Python modules
# ===========================================================================
from . import libspline

# ===========================================================================
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
    """ Write tecplot zone header depending on spatial dimension"""
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
    """ Close Tecplot file opened with openTecplot()"""
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


# =============================================================================
# pySpline classes
# =============================================================================
class Curve(object):
    """
    Create an instance of a b-spline curve. There are two
    ways to initialize the class

    * **Creation**: Create an instance of the spline class
      directly by supplying the required information. kwargs MUST
      contain the folloiwng information: ``k, t, coef``.


    * **LMS/Interpolation**: Create an instance of the spline class by
      using an interpolating spline to given data points or a LMS
      spline. The following keyword argument information is required:

      1. ``k`` Spline Order

      2. ``X`` real arry size (N, nDim) of data to fit. **OR**
          1. ``x`` (1D) and ``s`` for 1D
          2. ``x`` (1D) and ``y`` (1D) for 2D spatial curve
          3. ``x`` (1D) and ``y``` (1D) and ``z`` 1D for 3D spatial curve

      3. ``s`` real array of size (N). Optional for nDim >= 2

    Parameters
    ----------
    k : int
        Order for spline
    nCtl : int
        Number of control points. Used only by LMS initialization.
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
    >>> # Spatial interpolated seg
    >>> line_seg = pySpline.Curve(x=x, y=y, k=2)
    >>> # With explicit parameter values
    >>> line_seg = pySpline.Curve(x=x, y=y, k=2, s=s)
    >>> #LMS parabolic curve
    >>> parabola = pySpline.Curve(x=x, y=y, k=3)
    >>> #LMS parabolic curve with parameter values
    >>> parabola = pySpline.Curve(x=x, y=y, k=3, s=s)
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
            self.coef = numpy.atleast_2d(kwargs["coef"])
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
                self.X = numpy.atleast_2d(kwargs["X"])
                if numpy.ndim(kwargs["X"]) == 1:
                    self.X = numpy.transpose(self.X)
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = numpy.vstack([kwargs["x"], kwargs["y"], kwargs["z"]]).T
            elif "x" in kwargs and "y" in kwargs:
                self.X = numpy.vstack([kwargs["x"], kwargs["y"]]).T
            elif "x" in kwargs:
                self.X = numpy.transpose(numpy.atleast_2d(kwargs["x"]))
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

            T = numpy.zeros((self.N, self.nDim))
            # Compute tangents
            qq = numpy.zeros_like(self.X)
            T = numpy.zeros_like(self.X)
            deltaS = numpy.zeros(self.N)
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
                T[i] = T[i] / numpy.linalg.norm(T[i])

            # Final coefficients and t
            self.coef = numpy.zeros((2 * (self.N - 1) + 2, self.nDim))
            self.t = numpy.zeros(len(self.coef) + self.k)
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
            u = numpy.zeros(self.N)
            for i in range(0, self.N - 1):
                u[i + 1] = u[i] + numpy.linalg.norm(self.coef[2 * i + 2] - self.coef[2 * i + 1])

            for i in range(1, self.N - 1):
                self.t[2 * i + 2] = u[i] / u[self.N - 1]
                self.t[2 * i + 3] = u[i] / u[self.N - 1]

        else:  # lms or interpolate function
            assert "k" in kwargs and (
                "X" in kwargs or "x" in kwargs
            ), "Error: At least spline order, k and X (or x=, y=) \
MUST be defined for (interpolation) spline creation.\
nCtl=<number of control points> must be specified for a LMS fit"
            self.origData = True
            if "X" in kwargs:
                self.X = numpy.atleast_2d(kwargs["X"])
                if numpy.ndim(kwargs["X"]) == 1:
                    self.X = numpy.transpose(self.X)
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = numpy.vstack([kwargs["x"], kwargs["y"], kwargs["z"]]).T
            elif "x" in kwargs and "y" in kwargs:
                self.X = numpy.vstack([kwargs["x"], kwargs["y"]]).T
            elif "x" in kwargs:
                self.X = numpy.transpose(numpy.atleast_2d(kwargs["x"]))
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
                self.weights = numpy.ones(self.N)

            if "deriv" in kwargs and "derivPtr" in kwargs:
                self.deriv = checkInput(kwargs["deriv"], "deriv", float, 2)
                self.derivPtr = checkInput(kwargs["derivPtr"], "derivPtr", int, 1, (len(self.deriv),))
            else:
                self.deriv = None
                self.derivPtr = numpy.array([])

            if "derivWeights" in kwargs and self.deriv is not None:
                self.derivWeights = checkInput(kwargs["derivWeights"], "derivWeights", float, 1, (len(self.derivPtr),))
            else:
                if self.deriv is not None:
                    self.derivWeights = numpy.ones(len(self.deriv))
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
        suSelect = numpy.where(self.weights > 0.0)
        scSelect = numpy.where(self.weights <= 0.0)
        S = self.X[suSelect]
        su = self.s[suSelect]
        T = self.X[scSelect]
        sc = self.s[scSelect]
        weights = self.weights[numpy.where(self.weights > 0.0)]

        nu = len(S)
        nc = len(T)

        # And the derivative info
        if self.deriv is not None:
            sduSelect = numpy.where(self.derivWeights > 0.0)
            sdcSelect = numpy.where(self.derivWeights <= 0.0)
            S = numpy.vstack((S, self.deriv[sduSelect]))
            sdu = self.s[self.derivPtr][sduSelect]
            T = numpy.vstack((T, self.deriv[sdcSelect]))
            sdc = self.s[self.derivPtr][sdcSelect]
            weights = numpy.append(weights, self.derivWeights[numpy.where(self.derivWeights > 0.0)])
            ndu = len(sdu)
            ndc = len(sdc)
        else:
            sdu = numpy.array([], "d")
            sdc = numpy.array([], "d")
            ndu = 0
            ndc = 0

        W = _assembleMatrix(
            weights, numpy.arange(len(weights)), numpy.arange(len(weights) + 1), (len(weights), len(weights))
        )

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
        self.coef = numpy.zeros((self.nCtl, self.nDim), "d")

        # Get the 'N' jacobian
        nVals = numpy.zeros((nu + ndu) * self.k)  # |
        nRowPtr = numpy.zeros(nu + ndu + 1, "intc")  # | -> CSR formulation
        nColInd = numpy.zeros((nu + ndu) * self.k, "intc")  # |
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
        for i in range(nIter):
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

                mVals = numpy.zeros((nc + ndc) * self.k)  # |
                mRowPtr = numpy.zeros(nc + ndc + 1, "intc")  # | -> CSR
                mColInd = numpy.zeros((nc + ndc) * self.k, "intc")  # |

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
                    rhs = numpy.hstack((N.transpose() * W * S[:, idim], T[:, idim]))
                    self.coef[:, idim] = solve(rhs)[0 : self.nCtl]

            # end if (constr - not constrained

            # Run para correction
            libspline.curve_para_corr(self.t, self.k, self.s, self.coef.T, length, self.X.T)
        # end for (iter loop)

        # Check the RMS
        rms = 0.0
        for idim in range(self.nDim):
            rms += numpy.linalg.norm(N * self.coef[:, idim] - S[:, idim]) ** 2

        rms = numpy.sqrt(rms / self.N)

    def _getParameterization(self):
        """Compute a parametrization for the curve based on an
        arc-length formulation
        """
        self.s = numpy.zeros(self.N, "d")
        for i in range(self.N - 1):
            dist = 0
            for idim in range(self.nDim):
                dist += (self.X[i + 1, idim] - self.X[i, idim]) ** 2
            self.s[i + 1] = self.s[i] + numpy.sqrt(dist)

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
        t1 = numpy.hstack((self.t[0 : breakPt + self.k - 1].copy(), uu)) / uu
        t2 = (numpy.hstack((uu, self.t[breakPt:].copy())) - uu) / (1.0 - uu)

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
            length += numpy.linalg.norm(points[i] - points[i + 1])

        return length

    def calcGrevillePoints(self):
        """Calculate the Greville points"""
        self.gpts = numpy.zeros(self.nCtl)
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

        self.sdata = numpy.array(s)

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

        s = numpy.array(s).T
        if self.coef.dtype == numpy.dtype("d"):
            vals = libspline.eval_curve(numpy.atleast_1d(s), self.t, self.k, self.coef.T)
        else:
            vals = libspline.eval_curve_c(numpy.atleast_1d(s).astype("D"), self.t, self.k, self.coef.T)

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
        if self.coef.dtype == numpy.dtype("d"):
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
        if self.coef.dtype == numpy.dtype("d"):
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
        x0 = numpy.atleast_2d(x0)
        if "s" in kwargs:
            s = numpy.atleast_1d(kwargs["s"])
        else:
            s = -1 * numpy.ones(len(x0))

        if len(x0) != len(s):
            raise Error("projectPoint: The length of x0 and s must be the same")

        # If necessary get brute-force starting point
        if numpy.any(s < 0) or numpy.any(s > 1):
            self.computeData()
            s = libspline.point_curve_start(x0.T, self.sdata, self.data.T)

        D = numpy.zeros_like(x0)
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
            is equilivent to ||self(s) - inCurve(t)||_2.
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
            is equilivent to ||self(s) - inCurve(t)||_2.
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

                if numpy.linalg.norm(Diff) < eps:
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

        return numpy.array(uSol), numpy.array(tSol), numpy.array(diff)

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
        if numpy.mod(paraEntries - 11, 3) != 0:
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
            handle.write("%20.12g," % (numpy.real(self.t[i])))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0

        for i in range(self.nCtl):
            pos_counter += 1
            handle.write("%20.12g," % (1.0))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0

        for i in range(self.nCtl):
            for idim in range(3):
                pos_counter += 1
                handle.write("%20.12g," % (numpy.real(self.coef[i, idim])))
                if numpy.mod(pos_counter, 3) == 0:
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
        handle.write("%20.12g,%20.12g,0.0,0.0,0.0;         " % (numpy.min(self.t), numpy.max(self.t)))
        handle.write("  %7dP%7d\n" % (Pcount, counter))
        counter += 1
        Pcount += 2

        return Pcount, counter


class Surface(object):

    """
    Create an instance of a b-spline surface. There are two
    ways to initialize the class

    * **Creation**: Create an instance of the Surface class
      directly by supplying the required information: kwargs MUST
      contain the folloiwng information: ``ku, kv, tu, tv, coef``.

    * **LMS/Interpolation**: Create an instance of the Surface class by
      using an interpolating spline to given data points or a LMS
      spline. The following keyword argument information is required:

      1. ``ku`` and ``kv`` Spline Orders

      2. ``X`` real arry size (Nu, Nv, nDim) of data to fit. **OR**
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
                self.X = numpy.array(kwargs["X"])
                self.nDim = self.X.shape[2]
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = numpy.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 3))
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
            Td = numpy.zeros((self.Nu, self.Nv, 5, 3))

            def getT(Q, s):
                N = len(Q)
                T = numpy.zeros_like(Q)
                qq = numpy.zeros_like(Q)
                deltaS = numpy.zeros(N)
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
                    T[i] /= numpy.linalg.norm(T[i]) + 1e-16
                return T

            def interpKnots(u):
                t = numpy.zeros(2 * len(u) + 2 + 2)
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
                coef = numpy.zeros((3 * N - 2, self.nDim))
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
                    length += numpy.linalg.norm(Q[i + 1] - Q[i])
                return length

            def getD(Q, s):
                N = len(Q)
                D = numpy.zeros_like(Q)
                dd = numpy.zeros_like(D)
                deltaS = numpy.zeros(N)
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

            rowLen = numpy.zeros(self.Nv)
            colLen = numpy.zeros(self.Nu)

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
            coef = numpy.zeros((3 * self.Nu - 2, 3 * self.Nv - 2, self.nDim))

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
                    Td[i, j] /= numpy.linalg.norm(Td[i, j]) + 1e-16

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
                ind = numpy.zeros(2 * (N - 1) + 2, "intc")
                ind[0] = 0
                ind[-1] = 3 * N - 3
                ii = 1
                jj = 1
                for i in range(N - 1):
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
            assert (
                "ku" in kwargs
                and "kv" in kwargs
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
                    self.nDim = self.X.shape[2]
            elif "x" in kwargs and "y" in kwargs and "z" in kwargs:
                self.X = numpy.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 3))
                self.X[:, :, 0] = kwargs["x"]
                self.X[:, :, 1] = kwargs["y"]
                self.X[:, :, 2] = kwargs["z"]
                self.nDim = 3
            elif "x" in kwargs and "y" in kwargs:
                self.X = numpy.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 2))
                self.X[:, :, 0] = kwargs["x"]
                self.X[:, :, 1] = kwargs["y"]
                self.nDim = 2
            elif "x" in kwargs:
                self.X = numpy.zeros((kwargs["x"].shape[0], kwargs["x"].shape[1], 1))
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
                [self.V, self.U] = numpy.meshgrid(self.v, self.u)
            else:
                if self.nDim == 3:
                    self.u, self.v, self.U, self.V = self.calcParameterization()
                else:
                    raise Error(
                        "Automatic parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise."
                    )

            self.umin = 0
            self.umax = 1
            self.vmin = 0
            self.vmax = 1
            self.calcKnots()
            self.coef = numpy.zeros((self.nCtlu, self.nCtlv, self.nDim))
            if recompute:
                self.recompute()

    def recompute(self):
        """Recompute the surface if any data has been modified"""

        vals, rowPtr, colInd = libspline.surface_jacobian_wrap(
            self.U.T, self.V.T, self.tu, self.tv, self.ku, self.kv, self.nCtlu, self.nCtlv
        )

        N = _assembleMatrix(vals, colInd, rowPtr, (self.Nu * self.Nv, self.nCtlu * self.nCtlv))

        self.coef = numpy.zeros((self.nCtlu, self.nCtlv, self.nDim))
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

        u = numpy.zeros(self.Nu, "d")
        U = numpy.zeros((self.Nu, self.Nv), "d")
        singularSounter = 0
        # loop over each v, and average the 'u' parameter
        for j in range(self.Nv):
            temp = numpy.zeros(self.Nu, "d")

            for i in range(self.Nu - 1):
                temp[i + 1] = temp[i] + numpy.linalg.norm(self.X[i, j] - self.X[i + 1, j])

            if temp[-1] == 0:  # We have a singular point
                singularSounter += 1
                temp[:] = 0.0
                U[:, j] = numpy.linspace(0, 1, self.Nu)
            else:
                temp = temp / temp[-1]
                U[:, j] = temp.copy()

            u += temp  # accumulate the u-parameter calcs for each j

        u = u / (self.Nv - singularSounter)  # divide by the number of 'j's we had

        v = numpy.zeros(self.Nv, "d")
        V = numpy.zeros((self.Nu, self.Nv), "d")
        singularSounter = 0
        # loop over each v, and average the 'u' parameter
        for i in range(self.Nu):
            temp = numpy.zeros(self.Nv, "d")
            for j in range(self.Nv - 1):
                temp[j + 1] = temp[j] + numpy.linalg.norm(self.X[i, j] - self.X[i, j + 1])

            if temp[-1] == 0:  # We have a singular point
                singularSounter += 1
                temp[:] = 0.0
                V[i, :] = numpy.linspace(0, 1, self.Nv)
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
            self.tu = libspline.knots_interp(self.u, numpy.array([], "d"), self.ku)
            self.tv = libspline.knots_interp(self.v, numpy.array([], "d"), self.kv)
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
            if numpy.mod(self.Nu, 2) == 1:  # Its odd
                mid = (self.Nu - 1) // 2
                return self.X[0, 0], self.X[mid, 0], self.X[-1, 0]
            else:
                Xmid = 0.5 * (self.X[self.Nu // 2, 0] + self.X[self.Nu // 2 - 1, 0])
                return self.X[0, 0], Xmid, self.X[-1, 0]
        elif edge == 1:
            if numpy.mod(self.Nu, 2) == 1:  # Its odd
                mid = (self.Nu - 1) // 2
                return self.X[0, -1], self.X[mid, -1], self.X[-1, -1]
            else:
                Xmid = 0.5 * (self.X[self.Nu // 2, -1] + self.X[self.Nu // 2 - 1, -1])
                return self.X[0, -1], Xmid, self.X[-1, -1]
        elif edge == 2:
            if numpy.mod(self.Nv, 2) == 1:  # Its odd
                mid = (self.Nv - 1) // 2
                return self.X[0, 0], self.X[0, mid], self.X[0, -1]
            else:
                Xmid = 0.5 * (self.X[0, self.Nv // 2] + self.X[0, self.Nv // 2 - 1])
                return self.X[0, 0], Xmid, self.X[0, -1]
        elif edge == 3:
            if numpy.mod(self.Nv, 2) == 1:  # Its odd
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

            newCoef = numpy.zeros((self.nCtlu + actualR, self.nCtlv, self.nDim))

            for j in range(self.nCtlv):
                actualR, tNew, coefSlice, breakPt = libspline.insertknot(s, r, self.tu, self.ku, self.coef[:, j].T)
                newCoef[:, j] = coefSlice[:, 0 : self.nCtlu + actualR].T

            self.tu = tNew[0 : self.nCtlu + self.ku + actualR]
            self.nCtlu = self.nCtlu + actualR

        elif direction == "v":
            actualR, tNew, coefNew, breakPt = libspline.insertknot(s, r, self.tv, self.kv, self.coef[0, :].T)

            newCoef = numpy.zeros((self.nCtlu, self.nCtlv + actualR, self.nDim))

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
            t1 = numpy.hstack((self.tu[0 : breakPt + self.ku - 1].copy(), tt)) / tt
            t2 = (numpy.hstack((tt, self.tu[breakPt:].copy())) - tt) / (1.0 - tt)

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
            t1 = numpy.hstack((self.tv[0 : breakPt + self.kv - 1].copy(), tt)) / tt
            t2 = (numpy.hstack((tt, self.tv[breakPt:].copy())) - tt) / (1.0 - tt)

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

        u = numpy.array(u).T
        v = numpy.array(v).T
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")

        vals = libspline.eval_surface(
            numpy.atleast_2d(u), numpy.atleast_2d(v), self.tu, self.tv, self.ku, self.kv, self.coef.T
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
        u = numpy.array(u)
        v = numpy.array(v)
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")
        if not numpy.ndim(u) == 0:
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
        u = numpy.array(u)
        v = numpy.array(v)
        if not u.shape == v.shape:
            raise Error("u and v must have the same shape")
        if not numpy.ndim(u) == 0:
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

        Xmin = numpy.array([min(cx), min(cy), min(cz)])
        Xmax = numpy.array([max(cx), max(cy), max(cz)])

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

        x0 = numpy.atleast_2d(x0)
        if "u" in kwargs and "v" in kwargs:
            u = numpy.atleast_1d(kwargs["u"])
            v = numpy.atleast_1d(kwargs["v"])
        else:
            u = -1 * numpy.ones(len(x0))
            v = -1 * numpy.ones(len(x0))

        if not len(x0) == len(u) == len(v):
            raise Error("The length of x0 and u, v must be the same")

        # If necessary get brute-force starting point
        if numpy.any(u < 0) or numpy.any(u > 1) or numpy.any(v < 0):
            self.computeData()
            u, v = libspline.point_surface_start(x0.T, self.udata, self.vdata, self.data.T)

        D = numpy.zeros_like(x0)
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
        if numpy.any(u < 0) or numpy.any(u > 1) or numpy.any(v < 0):
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
            [V, U] = numpy.meshgrid(self.vdata, self.udata)
            self.data = self.getValue(U, V)

    def writeDirections(self, handle, isurf):
        """Write out and indication of the surface direction"""
        if self.nCtlu >= 3 and self.nCtlv >= 3:
            data = numpy.zeros((4, self.nDim))
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
        if numpy.mod(paraEntries - 10, 3) != 0:
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
            handle.write("%20.12g," % (numpy.real(self.tu[i])))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in range(len(self.tv)):
            pos_counter += 1
            handle.write("%20.12g," % (numpy.real(self.tv[i])))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in range(self.nCtlu * self.nCtlv):
            pos_counter += 1
            handle.write("%20.12g," % (1.0))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for j in range(self.nCtlv):
            for i in range(self.nCtlu):
                for idim in range(3):
                    pos_counter += 1
                    handle.write("%20.12g," % (numpy.real(self.coef[i, j, idim])))
                    if numpy.mod(pos_counter, 3) == 0:
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
                handle.write("%20.12g," % (numpy.real(self.umin)))
            if i == 1:
                handle.write("%20.12g," % (numpy.real(self.umax)))
            if i == 2:
                handle.write("%20.12g," % (numpy.real(self.vmin)))
            if i == 3:
                # semi-colon for the last entity
                handle.write("%20.12g;" % (numpy.real(self.vmax)))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write("  %7dP%7d\n" % (Pcount, counter))
                counter += 1
                pos_counter = 0
            else:  # We have to close it up anyway
                if i == 3:
                    for j in range(3 - pos_counter):
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
                x = checkInput(x, "x", float, 3)
                y = checkInput(x, "y", float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[3], 3))
                self.X[:, :, :, 0] = x
                self.X[:, :, :, 1] = y
                self.nDim = 2
            elif "x" in kwargs:
                x = checkInput(x, "x", float, 3)
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
            raise Error(
                "A single argument passed to bilinear " "surface must contain 4 points and be of " "size (4, 3)"
            )
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


# ==============================================================================
# Class Test
# ==============================================================================
if __name__ == "__main__":
    print("There are two examples in the example directory.")
    print("Look at test_curve.py and test_surf.py for more information")
