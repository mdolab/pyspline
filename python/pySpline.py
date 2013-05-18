'''
pySpline

Contains classes for working with B-spline curves, surfaces and
volumes

Copyright (c) 2009 by G. Kenway
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 24/05/2009$


Developers:
-----------
- Gaetan Kenway (GKK)

History
-------
	v. 1.0 - Initial Class Creation (GKK, 2009)
'''

# ===========================================================================
# Standard Python modules
# ===========================================================================
import sys, warnings

# ===========================================================================
# External Python modules
# ===========================================================================
import numpy

USE_SCIPY = False
try:
    import scipy
    import scipy.sparse.linalg.dsolve
    USE_SCIPY = True
except:
    warnings.warn("pySpline: Scipy could not be importd. Fitting is possible,\
 but may be slow for very large fits")

# ===========================================================================
# Custom Python modules
# ===========================================================================

import pyspline
from mdo_import_helper import mpiPrint, import_modules
exec(import_modules('geo_utils'))

# ===========================================================================
def writeTecplot1D(handle, name, data):
    """A Generic function to write a 1D data to tecplot. 
    Input:
        handle: an open file handle 
        name, str: Name of zone in Tecplot
        data, numpy array: 2D array data to write
    Output:
        None

    """
    nx = data.shape[0]
    ndim = data.shape[1]
    handle.write('Zone T=\"%s\" I=%d\n'%(name, nx))
    handle.write('DATAPACKING=POINT\n')
    for i in xrange(nx):
        for idim in xrange(ndim):
            handle.write('%f '%(data[i, idim]))
        # end for
        handle.write('\n')
    # end for

    return

def writeTecplot2D(handle, name, data):
    """A Generic function to write a 2D data to tecplot. 
    Input:
        handle: an open file handle 
        name, str: Name of zone in Tecplot
        data, numpy array: 2D array data to write
    Output:
        None

    """
    nx = data.shape[0]
    ny = data.shape[1]
    ndim = data.shape[2]
    handle.write('Zone T=\"%s\" I=%d J=%d\n'%(name, nx, ny))
    handle.write('DATAPACKING=POINT\n')
    for j in xrange(ny):
        for i in xrange(nx):
            for idim in xrange(ndim):
                handle.write('%20.16g '%(data[i, j, idim]))
            # end for
            handle.write('\n')
        # end for
    # end for

    return

def writeTecplot3D(handle, name, data):
    """A Generic function to write a 3D data to tecplot. 
    Input:
        handle: an open file handle
        name, str: Name of zone in Tecplot
        data, numpy array: 3D array data to write
    Output:
        None

    """
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    ndim = data.shape[3]
    handle.write('Zone T=\"%s\" I=%d J=%d K=%d\n'%(name, nx, ny, nz))
    handle.write('DATAPACKING=POINT\n')
    for k in xrange(nz):
        for j in xrange(ny):
            for i in xrange(nx):
                for idim in xrange(ndim):
                    handle.write('%f '%(data[i, j, k, idim]))
                # end for
                handle.write('\n')
            # end for 
        # end for
    # end for

    return

def _writeHeader(f, ndim):
    """ Write tecplote zone header depending on spatial dimension"""
    if ndim == 1:
        f.write ('VARIABLES = "X"\n')
    elif ndim == 2:
        f.write ('VARIABLES = "X","Y"\n')
    else:
        f.write ('VARIABLES = "X", "Y","Z"\n')
    # end if

def openTecplot(file_name, ndim):
    """A Generic function to open a Tecplot file to write spatial data.
    Input:
        file_name, str: file name to open
        ndim, int: Number of spatial dimension
    Output:
        f: file handle to file. 
    """
    mpiPrint('Opening ascii Tecplot File: %s'%(file_name))
    f = open(file_name, 'w')
    _writeHeader(f, ndim)

    return f

def closeTecplot(f):
    """ Close Tecplot file opened with openTecplot"""
    f.close()
    
    return 

def _assembleMatrix(data, indices, indptr, shape):
    """
    Generic assemble matrix function to create a CSR matrix if scipy
    is available or a general dense matrix otherwise. 
    Input:
        data, numpy array: Array of values
        indices, numpy array: Array of colume indices
        indptr, numpy array: Array of row pointers
        shape, tuple-like:, Shape of matrix e.g. (5,6)
    Ouput:
        M: Matrix-like array either a sparse or dense matrix
    
        """
    if USE_SCIPY:
        M = scipy.sparse.csr_matrix((data, indices, indptr), shape)
    else:
        M = numpy.zeros(shape,'d')
        # Add the required values
        counter = 0
        for ii in xrange(len(indptr)-1):
            for jj in xrange(indptr[ii], indptr[ii+1]):
                M[ii, indices[counter]] = \
                    data[counter]
                counter += 1
            # end for
        # end for
    # end if

    return M

# =============================================================================
# pySpline classes
# =============================================================================
class curve(object):
    """
    A Generic b-spline curve class. See the __init__ documentation for how to
    create a curve object

    """
    def __init__(self, **kwargs):
        """
        Create an instance of a b-spline curve. There are two
        ways to initialize the class

        Creation: Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            k, integer: Order for spline
            t, real array: Knot vector 
            coef, real array size(Nctl,nDim): Array of control points

        LMS/Interpolation

        Create an instance of the spline class by using an
        interpolating spline to given data points or a LMS spline. The
        following keyword argument information is required:

            k, integer: Order for spline
            X, real, array, size(len(s),nDim): Array of data to fit 
              --> OR  x=<vals>, len(s) ->1D interpolation
              --> OR  x=<vals>,y=<vals>, each of len(s) 
                      -> 2D parametric curve
              --> OR  x=<vals>,y=<vals>,z=<vals> each of len(s) 
                      -> 3D parametric curve
            s, real, array: (OPTIONAL for nDim >=2 ) Array of s values 
            
        For LMS spline:
             Nctl -> (integer) number of control points (Less than len(s)

        Optional Data:
        
        weights:   A array of weighting factors for each fitting point
                   A value of -1 can be used to exactly constrain a point

        deriv and deriv_ptr: deriv is a array which contains
                             derivative information at point indicies
                             defined by deriv_ptr deriv_weights: A
                             array of len deriv_ptr which contains
                             derivative weighting A value of -1 can be
                             used to exactly constrain a derivative

        niter:  The number of Hoscheks parameter corrections to run
        """
        
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
        
        self.length = None
        self.gpts = None
        self.data = None
        self.local_interp = False
        # We have provided information to create curve directly
        if 'k' in kwargs and 't' in kwargs and 'coef' in kwargs: 
            self.s = None
            self.X = None
            self.N = None
            self.k = geo_utils.checkInput(kwargs['k'], 'k', int, 0)
            self.coef = numpy.atleast_2d(kwargs['coef'])
            self.Nctl = self.coef.shape[0]
            self.nDim = self.coef.shape[1]
            self.t = geo_utils.checkInput(
                kwargs['t'], 't', float, 1, self.Nctl+self.k)
            self.orig_data = False
            self._calcGrevillePoints()

        elif 'local_interp' in kwargs:
            # Local, non-global interpolation. We could use this for
            # second-order (linear interpolation) but there is not
            # much point, since it would be the same as 'interp'
            # below. 
            self.local_interp = True
            if 'k' in kwargs:
                mpiPrint('* Warning: k is ignored for local_interp. \
Cubic spline order is always used')
            self.k = 4
            
            self.orig_data = True
            if 'X' in kwargs:
                self.X = numpy.atleast_2d(kwargs['X'])
                if numpy.rank(kwargs['X']) == 1:
                    self.X = numpy.transpose(self.X)
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = numpy.vstack([kwargs['x'], kwargs['y'], kwargs['z']]).T
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = numpy.vstack([kwargs['x'], kwargs['y']]).T
            elif 'x' in kwargs:
                self.X = numpy.transpose(numpy.atleast_2d(kwargs['x']))
            # enf if

            self.nDim = self.X.shape[1]
            self.N = self.X.shape[0]

            if 's' in kwargs:
                self.s = geo_utils.checkInput(
                    kwargs['s'], 's', float, 1, self.N)
            else:
                assert self.nDim > 1, 'Error, pySpline: For 1D splines,\
 the basis, s must be given'
                self._getParameterization()
            # end if

            # Now we have the data we need to generate the local
            # interpolation
     
            T = numpy.zeros((self.N, self.nDim))
            # Compute tangents
                 qq = numpy.zeros_like(self.X)
            T  = numpy.zeros_like(self.X)
            delta_s = numpy.zeros(self.N)
            for i in xrange(1, self.N):
                delta_s[i] = self.s[i]-self.s[i-1]
                qq[i, :] = self.X[i]-self.X[i-1] 
            # end for

            for i in xrange(1, self.N-1):
                a = delta_s[i]/(delta_s[i] + delta_s[i+1])
                T[i] = (1-a)*qq[i] + a*qq[i+1]
            # end for

            # Do the start and end points: (eqn: 9.32, The NURBS book)
            T[0] = 2*qq[1]/delta_s[1] - T[1]
            T[-1] = 2*qq[-1]/delta_s[-1] - T[-2]

            # Normalize
            for i in xrange(self.N):
                T[i] /= numpy.linalg.norm(T[i])
            # end for
                
            # Final coefficients and t
            self.coef = numpy.zeros((2*(self.N-1)+2,self.nDim))
            self.t    = numpy.zeros(len(self.coef) + self.k)
            
            # End Pts
            self.coef[0] = self.X[0].copy()
            self.coef[-1] = self.X[-1].copy()

            # Interior coefficients
            for i in xrange(self.N-1):
                a = self.length*(self.s[i+1]-self.s[i])
                self.coef[2*i+1] = self.X[i] + a/3.0*T[i]
                self.coef[2*i+2] = self.X[i+1] - a/3.0*T[i+1]
            # end for

            # Knots
            self.t[-4:] = 1.0
            u = numpy.zeros(self.N)
            for i in xrange(0,self.N-1):
                u[i+1] = u[i] + numpy.linalg.norm(self.coef[2*i+2]-self.coef[2*i+1])
            # end for

            for i in xrange(1,self.N-1):
                self.t[2*i+2] = u[i]/u[self.N-1]
                self.t[2*i+3] = u[i]/u[self.N-1]
            # end for

        else: #lms or interpolate function
            assert 'k' in kwargs and ('X' in kwargs or 'x' in kwargs), \
                'Error: At least spline order, k and X (or x=, y=) \
MUST be defined for (interpolation) spline creation.\
Nctl=<number of control points> must be specified for a LMS fit'
            self.orig_data = True
            if 'X' in kwargs:
                self.X = numpy.atleast_2d(kwargs['X'])
                if numpy.rank(kwargs['X']) == 1:
                    self.X = numpy.transpose(self.X)
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = numpy.vstack([kwargs['x'], kwargs['y'], kwargs['z']]).T
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = numpy.vstack([kwargs['x'], kwargs['y']]).T
            elif 'x' in kwargs:
                self.X = numpy.transpose(numpy.atleast_2d(kwargs['x']))
            # enf if

            self.nDim = self.X.shape[1]
            self.N = self.X.shape[0]

            self.k = geo_utils.checkInput(kwargs['k'], 'k', int, 0)
            self.t = None
            if 'niter' in kwargs: 
                self.niter = geo_utils.checkInput(
                    kwargs['niter'], 'niter', int, 0)
            else:
                self.niter = 1
            # end if

            if 's' in kwargs:
                self.s = geo_utils.checkInput(
                    kwargs['s'], 's', float, 1, self.N)
            else:
                assert self.nDim > 1, 'Error, pySpline: For 1D splines,\
 the basis, s must be given'
                self._getParameterization()
            # end if
            
            if 'weights' in kwargs:
                self.weights = geo_utils.checkInput(
                    kwargs['weights'], 'weights', float, 1, self.N)
            else:
                self.weights = numpy.ones(self.N)
            # end if

            if 'deriv' in kwargs and 'deriv_ptr' in kwargs:
                self.deriv = geo_utils.checkInput(kwargs['deriv'], float, 2)
                self.deriv_ptr = geo_utils.checkInput(kwargs['deriv_ptr'], 
                                            'deriv_ptr', int, 1, 
                                            len(self.deriv))
            else:
                self.deriv = None
                self.deriv_ptr = numpy.array([])
            # end if
            if 'deriv_weights' in kwargs and self.deriv != None:
                self.deriv_weights = geo_utils.checkInput(
                    kwargs['deriv_weights'], 'deriv_weights', float, 1, 
                    len(self.deriv_ptr))
            else:
                if self.deriv != None:
                    self.deriv_weights = numpy.ones(len(self.deriv))
                else:
                    self.deriv_weights = None
                # end if
            # end if

            # We are doing an interpolation...set all weights to 1
            if not 'Nctl' in kwargs:
                self.interp = True
                self.weights[:] = 1
                if self.deriv_weights != None:
                    self.deriv_weights[:] = 1
                # end if
            else:
                self.interp = False
                self.Nctl = geo_utils.checkInput(kwargs['Nctl'], 'Nctl', int, 0)
            # end if
            self.recompute(self.niter, computeKnots=True)
        return 

    def recompute(self, niter, computeKnots=True):
        """
        Run iterations of Hoscheks Parameter Correction on the current curve
        Input:
            niter, integer: The number of parameter correction iterations to run
        Returns:
            None

            """

        # Return if we don't have original data to fit
        if not self.orig_data or self.local_interp:
            return

        # Do the sparation between the constrained and unconstrained: 
        # u -> unconstrained
        # s -> constrainted
        su_select = numpy.where(self.weights > 0.0)
        sc_select = numpy.where(self.weights <= 0.0)
        S  = self.X[su_select]
        su = self.s[su_select]
        T  = self.X[sc_select]
        sc = self.s[sc_select]
        weights = self.weights[numpy.where(self.weights > 0.0)]

        nu = len(S)
        nc = len(T)
     
        # And the derivative info
        if self.deriv != None:
            sdu_select = numpy.where(self.deriv_weights > 0.0)
            sdc_select = numpy.where(self.deriv_weights <= 0.0)
            S = numpy.vstack((S, self.deriv[sdu_select]))
            sdu = self.s[self.deriv_ptr][sdu_select]
            T = numpy.vstack((T, self.deriv[sdc_select]))
            sdc = self.s[self.deriv_ptr][sdc_select]
            weights = numpy.append(weights, self.deriv_weights[
                    numpy.where(self.deriv_weights > 0.0)])
            ndu = len(sdu)
            ndc = len(sdc)
        else:
            sdu = numpy.array([], 'd')
            sdc = numpy.array([], 'd')
            ndu = 0
            ndc = 0
        # end if
        W = _assembleMatrix(weights, numpy.arange(len(weights)),
                            numpy.arange(len(weights)+1),
                            (len(weights),len(weights)))
            
        if self.interp:
            self.Nctl = nu+nc+ndu+ndc
            self.niter = 1
        # end if

        # Sanity check to make sure k is ok
        if nu+nc+ndu+ndc < self.k:
            self.k = nu+nc+ndu+ndc
        # end if

        if computeKnots:
            # Generate the knot vector, if necessary greville points and
            # empty coefficients
            if self.interp:
                self.t = pyspline.knots_interp(self.s, self.deriv_ptr, self.k)
            else:
                self.t = pyspline.knots_lms(self.s, self.Nctl, self.k)
            # end if
        # end if
        self._calcGrevillePoints()
        self.coef = numpy.zeros((self.Nctl, self.nDim), 'd')

        # Get the 'N' jacobain
        N_vals = numpy.zeros((nu+ndu)*self.k)           # |
        N_row_ptr = numpy.zeros(nu+ndu+1, 'intc')       # | -> CSR formulation
        N_col_ind = numpy.zeros((nu+ndu)*self.k, 'intc')# |
        pyspline.curve_jacobian_wrap(
            su, sdu, self.t, self.k, self.Nctl, N_vals, N_row_ptr, N_col_ind)
        N = _assembleMatrix(N_vals, N_col_ind, N_row_ptr, (nu+ndu, self.Nctl)).tocsc()

        if self.interp:
            if USE_SCIPY:
                # Factorize once for efficiency
                solve = scipy.sparse.linalg.dsolve.factorized( N )
                for idim in xrange(self.nDim):
                    self.coef[:, idim] = solve(S[:, idim])
                # end for
            else:
                for idim in xrange(self.nDim):
                    self.coef[:, idim] = numpy.linalg.solve(N, S[:, idim])
                # end for
            # end if
            return
        # end if

        # If we do NOT have an interpolation:
        length = pyspline.poly_length(self.X.T)
        for i in xrange(niter):
            su = self.s[su_select]
            sc = self.s[sc_select]
            if self.deriv != None:
                sdu = self.s[self.deriv_ptr][sdu_select]
                sdc = self.s[self.deriv_ptr][sdc_select]
            # end if
                
            pyspline.curve_jacobian_wrap(su, sdu, self.t, self.k, self.Nctl, 
                                         N_vals, N_row_ptr, N_col_ind)
            NTWN = (N.transpose()*W*N).tocsc() # We need this either way

            if nc + ndc == 0: # We are doing LMS but no
                              # constraints...just a straight weighted
                              # LMS

                if USE_SCIPY:
                    # Factorize once for efficiency                    

                    solve = scipy.sparse.linalg.dsolve.factorized( NTWN )

                    for idim in xrange(self.nDim):
                        self.coef[:, idim] = solve(N.transpose()*W*S[:, idim])
                    # end for
                else:
                    for idim in xrange(self.nDim):
                        self.coef[:, idim] = numpy.linalg.solve(
                            NTWN, N.transpose()*W*S[:, idim])
                    # end for
                # end if

            else:  # Now its more complicated since we have
                   # constraints --only works with scipy Sparse
                   # matrices

                if USE_SCIPY:
                    M_vals = numpy.zeros((nc+ndc)*self.k)           #|
                    M_row_ptr = numpy.zeros(nc+ndc+1, 'intc')       #| -> CSR
                    M_col_ind = numpy.zeros((nc+ndc)*self.k, 'intc')#|

                    pyspline.curve_jacobian_wrap(sc, sdc, self.t, self.k, 
                                                 self.Nctl, M_vals, M_row_ptr,
                                                 M_col_ind)
                    M = _assembleMatrix(M_vals, M_col_ind, M_row_ptr,
                                        (nc+ndc, self.Nctl))

                    # Now we must assemble the constrained jacobian
                    # [ N^T*W*T      M^T][P] = [ N^T*W*S]
                    # [ M            0  ][R]   [ T      ] 

                    MT   = M.transpose().tocsr()

                    j_val, j_col_ind, j_row_ptr = pyspline.constr_jac(
                        NTWN.data, NTWN.indptr, NTWN.indices, MT.data, 
                        MT.indptr, MT.indices, 
                        M.data, M.indptr, M.indices, self.Nctl)

                    # Create sparse csr matrix and factorize
                    J = _assembleMatrix(j_val, j_col_ind, j_row_ptr, 
                                        (self.Nctl+nc+ndc, self.Nctl+nc+ndc))

                    # Factorize once for efficiency
                    solve = scipy.sparse.linalg.dsolve.factorized( J ) 
                    for idim in xrange(self.nDim):
                        rhs = numpy.hstack((N.transpose()*W*S[:, idim], 
                                            T[:, idim]))
                        self.coef[:, idim] = solve(rhs)[0:self.Nctl]
                    # end for
                else:
                    mpiPrint('Constrained curve fitting is only availble when \
scipy is used.')
                    sys.exit(1)
                # end if
            # end if (constr - not constrained

            # Run para correction
            pyspline.curve_para_corr(self.t, self.k, self.s, 
                                     self.coef.T, length, self.X.T)
        # end for (iter loop)
        # Check the RMS
        rms = 0.0
        for idim in xrange(self.nDim):
            rms += numpy.linalg.norm(N*self.coef[:, idim]-S[:, idim])**2
        # end for
        rms = numpy.sqrt(rms/self.N)
        mpiPrint('Rms is: %f'%(rms), self.NO_PRINT)
      
        return

    def _getParameterization(self):
        """ Compute a paramerization for the curve based on an
        arc-length formulation

        """
        self.s = numpy.zeros(self.N, 'd')
        for i in xrange(self.N-1):
            dist = 0
            for idim in xrange(self.nDim):
                dist += (self.X[i+1, idim] - self.X[i, idim])**2
            # end for
            self.s[i+1] = self.s[i] + numpy.sqrt(dist)
            # end for
        # end for
        self.length = self.s[-1]
        self.s /= self.s[-1]

        return

    def reverse(self):
        """ Reverse the direction of this curve

        """
        self.coef = self.coef[::-1, :]
        self.t    = 1-self.t[::-1]

        return

    def insertKnot(self, u, r):
        """
        Insert at knot in the curve at u
        Required:
            u : Parametric position to split at
            r : The number of times to insert the knot
            """
        u = geo_utils.checkInput(u, 'u', float, 0)
        r = geo_utils.checkInput(r, 'r', int, 0)
        
        if u <= 0:
            return
        if u >= 1.0:
            return
        
        actual_r, t_new, coef_new, break_pt = pyspline.insertknot(
            u, r, self.t, self.k, self.coef.T)
        self.t = t_new[0:self.Nctl+self.k+actual_r]
        self.coef = coef_new[:,0:self.Nctl+actual_r].T
        self.Nctl = self.Nctl + actual_r

        # break_pt is converted to zero based ordering here!!!
        return actual_r, break_pt-1

    def splitCurve(self, u):
        """ 
        Split the curve at parametric position
        Input:
            u : Parametric position for split
        Returns:
            curve1, curve2: The two split spline curve
        """
        u = geo_utils.checkInput(u, 'u', float, 0)

        if u <= 0.0:
            return None, curve(t=self.t.copy(), k=self.k, 
                                  coef=self.coef.copy())
        if u >= 1.0:
            return curve(t=self.t.copy(), k=self.k, 
                                  coef=self.coef.copy()), None

        
        r, break_pt = self.insertKnot(u, self.k-1)

        if r == self.k -1:
            # A new knot was inserted, so index of relevant place
            # in knot vector is
            break_pt = break_pt + 1
        else:
            break_pt = break_pt
        # end if

        uu = self.t[break_pt]

        # Process knot vectors:
        t1 = numpy.hstack((self.t[0:break_pt+self.k-1].copy(),uu))/uu
        t2 = (numpy.hstack((uu,self.t[break_pt:].copy()))-uu)/(1.0-uu)
            
        coef1 = self.coef[0:break_pt,:].copy()
        coef2 = self.coef[break_pt-1:,:].copy()

        return \
            curve(t=t1, k=self.k, coef=coef1),\
            curve(t=t2, k=self.k, coef=coef2)

    def windowCurve(self, uLow, uHigh):
        dummy, curve = self.splitCurve(uLow)
        curve, dummy = curve.splitCurve((uHigh-uLow)/(1.0-uLow))

        return curve

    def getLength(self):
        """ Compute the length of the curve using the Eculdian Norm
        Input:
            None
        Returns:
            Length, real scalar: Lenght of curve
            """
        points = self.getValue(self.gpts)
        length = 0
        for i in xrange(len(points)-1):
            length += geo_utils.e_dist(points[i], points[i+1])
        # end for

        return length

    def _calcGrevillePoints(self):
        """Calculate the Greville points"""

        self.gpts = numpy.zeros(self.Nctl)
        for i in xrange(self.Nctl):
            for n in xrange(self.k-1): #degree loop
                self.gpts[i] += self.t[i+n+1]
            # end for
            self.gpts[i] /= (self.k-1)
        # end for

        return
    
    def _calcInterpolatedGrevillePoints(self):
        self._calcGrevillePoints()
        s = [self.gpts[0]]
        N = 2
        for i in xrange(len(self.gpts)-1):
            for j in xrange(N):
                s.append( (self.gpts[i+1]-self.gpts[i])*(j+1)/(N+1) + self.gpts[i])
            s.append(self.gpts[i+1])
        # end for

        self.sdata = numpy.array(s)

        return

    def __call__(self, s):
        """
        Equivilant to getValue()
        """
        return self.getValue(s)
    
    def getValue(self, s):
        """
        Evalue the spline at parametric position, s
        Input:
           s: s can be a scalar or vector of values
        Returns:
           values: An array of shape(len(s), nDim)) if s is an array
                   or array of len(nDim) if nDim == 1
                   """
     
        s = numpy.array(s).T
        if self.coef.dtype == numpy.dtype('d'):
            vals = pyspline.eval_curve(numpy.atleast_1d(s), 
                                       self.t, self.k, self.coef.T)
        else:
            vals = pyspline.eval_curve_c(numpy.atleast_1d(s), 
                                         self.t, self.k, self.coef.T)
        # end if

        return vals.squeeze().T
        
    def getDerivative(self, s):
        """
        Evalue the derivatie of the spline at parametric position, s
        Input:
            s, scalar: A real scalar
        Returns:
            derivative: The first derivative [ndim] vector

            """
        if self.coef.dtype == numpy.dtype('d'):
            derivative = pyspline.eval_curve_deriv(
                s, self.t, self.k, self.coef.T).squeeze()
        else:
            derivative = pyspline.eval_curve_deriv_c(
                s, self.t, self.k, self.coef.T).squeeze()
        # end if

        return derivative
        
    def getSecondDerivative(self, s):
        """
        Evalue the 2nd derivatie of the spline at parametric position, s
        Input:
            s, scalar: A real scalar
        Returns:
            derivative: The first derivative [ndim] vector
            
            """
        if self.coef.dtype == numpy.dtype('d'):
            derivative = pyspline.eval_curve_deriv2(
                s, self.t, self.k, self.coef.T).squeeze()
        else:
            derivative = pyspline.eval_curve_deriv2_c(
                s, self.t, self.k, self.coef.T).squeeze()
        # end if

        return derivative

    def projectPoint(self, x0, Niter=25, eps=1e-10, **kwargs):
        """Project a point x0 or (points x0) onto the curve and return
        parametric position(s)

        curve: 
        Required Arguments:
            x0   : A point(s) in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            s    : Initial guess for position
        Returns:
            s: Parametric position(s) of closest point on curve
            D: Distance(s) from point to curve

        """
        x0 = numpy.atleast_2d(x0)
        if 's' in kwargs: 
            s = numpy.atleast_1d(kwargs['s'])
        else:
            s = -1*numpy.ones(len(x0))
        # end if

        assert len(x0) == len(s), 'projectPoint: The length of x0\
 and s must be the same'

        # If necessary get brute-force starting point
        if numpy.any(s<0) or numpy.any(s>1):
            self._computeData()
            s = pyspline.point_curve_start(x0.T, self.sdata, self.data.T)

        D = numpy.zeros_like(x0)
        for i in xrange(len(x0)):
            s[i], D[i] = pyspline.point_curve(x0[i], self.t, self.k, self.coef.T, 
                                              Niter, eps, s[i])
        # end for

        return s.squeeze(), D.squeeze()

    def projectCurve(self, in_curve, Niter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this curve (self) and a second
        curve passed in (curve)
        Required Arguments:
            curve: A pyspline curve class to do the projection with
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            s    : Initial guess for curve1 (this curve class)
            t    : Initial guess for curve2 (cuve passed in)
        Returns:
            s    : Parametric position on curve1 (this class)
            t    : Parametric position on curve2 (curve passed in)
            D    : Distance between curve1(s) and curve2(t)
            """
        s = -1
        t = -1
        if 's' in kwargs:
            s = geo_utils.checkInput(kwargs['s'], 's', float, 0)
        if 't' in kwargs:
            t = geo_utils.checkInput(kwargs['t'], 't', float, 0)
        eps   = geo_utils.checkInput(eps, 'eps', float, 0)

        if s < 0 or s > 1 or t < 0 or t >1:
            self._computeData()
            in_curve._computeData()
            s, t = pyspline.curve_curve_start(
                self.data.T, self.sdata, in_curve.data.T, in_curve.sdata)

        # end if
        s, t, Diff = pyspline.curve_curve(
            self.t, self.k, self.coef.T, in_curve.t, in_curve.k, 
            in_curve.coef.T, Niter, eps, s, t)
        
        return s, t, Diff
        
    def projectCurveMultiSol(self, in_curve, Niter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this curve (self) and a
        second curve passed in (curve). This tries to find more than
        one solution if there are more than one.

        Required Arguments:
            curve: A pyspline curve class to do the projection with
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            s    : Initial guess for curve1 (this curve class)
            t    : Initial guess for curve2 (cuve passed in)
        Returns:
            s    : Parametric position on curve1 (this class)
            t    : Parametric position on curve2 (curve passed in)
            D    : Distance between curve1(s) and curve2(t)
            """
        s = -1
        t = -1
        if 's' in kwargs:
            s = geo_utils.checkInput(kwargs['s'], 's', float, 0)
        if 't' in kwargs:
            t = geo_utils.checkInput(kwargs['t'], 't', float, 0)
        eps   = geo_utils.checkInput(eps, 'eps', float, 0)

        self._computeData()
        in_curve._computeData()

        u_sol = []
        t_sol = []
        diff  = []
        for i in xrange(len(self.sdata)):
            for j in xrange(len(in_curve.sdata)):
                s, t, Diff = pyspline.curve_curve(
                    self.t, self.k, self.coef.T, in_curve.t, in_curve.k, 
                    in_curve.coef.T, Niter, eps, self.sdata[i], in_curve.sdata[j])

                if numpy.linalg.norm(Diff) < eps:
                    # Its a solution. Check it it is already in list:
                    if len(u_sol) == 0:
                        u_sol.append(s)
                        t_sol.append(t)
                        diff.append(Diff)
                    # end if

                    for ii in xrange(len(u_sol)):
                        if abs(u_sol[ii] - s) < eps and abs(t_sol[ii]-t) < eps:
                            pass
                        else:
                            u_sol.append(s), t_sol.append(t), diff.append(Diff)
                        # end if
                    # end for
            # end for
        # end for

        return numpy.array(u_sol), numpy.array(t_sol), numpy.array(diff)
        
    def _computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """
        # We will base the data on interpolated greville points

        if self.data is None:
            self._calcInterpolatedGrevillePoints()
            self.data = self.getValue(self.sdata)
        # end if
    
    def writeTecplot(self, file_name, curve=True, coef=True, orig=True):
        """
        Write the cuve to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            curve: boolean flag to write the interpolated curve (default=True)
            coef : boolean flag to write the control points (defaut=True)
            orig : boolean flag to write the original data (default=True)
            size : A distance to determine the output resolution of 
                   interpolated curve
        Returns:
            None
            """

        f = openTecplot(file_name, self.nDim)
        if curve:
            self._writeTecplotCurve(f)
        if coef:
            writeTecplot1D(f, 'control_pts', self.coef)        
        if orig and self.orig_data:
            writeTecplot1D(f, 'orig_data', self.X)
        
        closeTecplot(f)

        return 

    def _writeTecplotCurve(self, handle, **kwargs):
        """
        Write the interpolated curve to a handle
        """
        self._computeData()
        writeTecplot1D(handle, 'interpolated', self.data)

        return         
  
class surface(object):

    def __init__(self, recompute=True, **kwargs):
        """
        Create an instance of a b-spline surface. There are two
        ways to initialize the class

        Creation: Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu, Nctlv, nDim): Array of control points

        LMS/Interpolation
        Create an instance of the spline class by using an interpolating spline to given data points. **kwarg
        or a LMS spline. The following information is required:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            X, real, array, size(len(u), len(v), nDim): Array of data to fit 
              --> OR  x=<vals>, (len(u), len(v) -2D interpolation
              --> OR  x=<vals>, y=<vals>, each of size (len(u), len(v) 
              --> OR  x=<vals>, y=<vals>, z=<vals> each of len(s) ->3D parametric surface
            u, v, real, arrays: (OPTIONAL for nDim == 3 ) Arrays of u and v values
            
        For LMS spline:
            Nctlu -> (integer) number of control points in u
            Nctlv -> (integer) number of control points in v

        Optional Data:
           u    : array of u parameter locations (optional for ndim == 3)
           v    : array of v parameter locations (optional for ndim == 3)
           niter:  The number of Hoschek\'s parameter corrections to run
           """
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      

        self.edge_curves = [None, None, None, None]
        self.data = None
        if 'ku' in kwargs and 'kv' in kwargs and 'tu' in kwargs and 'tv' in \
                kwargs and 'coef' in kwargs:
            self.X = None
            self.u = None
            self.v = None
            self.U = None
            self.V = None
            self.ku = geo_utils.checkInput(kwargs['ku'], 'ku', int, 0)
            self.kv = geo_utils.checkInput(kwargs['kv'], 'kv', int, 0)
            self.coef = geo_utils.checkInput(kwargs['coef'], 'coef', float, 3)
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.tu = geo_utils.checkInput(kwargs['tu'], 'tu', float, 1, 
                                 self.Nctlu+self.ku)
            self.tv = geo_utils.checkInput(kwargs['tv'], 'tv', float, 1,
                                 self.Nctlv+self.kv)
            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.nDim = self.coef.shape[2]
            self.orig_data = False
            self._setEdgeCurves()
            self.interp =  False
            return
        else: # We have LMS/Interpolate
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and \
                ( 'X' in kwargs or 'x' in kwargs or \
                      ( 'x' in kwargs and 'y' in kwargs) or \
                      ( 'x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
                      'Error: ku, kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate'

            if 'X' in kwargs:
                self.X  = numpy.array(kwargs['X'])
                if len(self.X.shape) == 1:
                    self.nDim = 1
                else:
                    self.nDim = self.X.shape[2]
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = numpy.zeros((
                        kwargs['x'].shape[0], kwargs['x'].shape[1], 3))
                self.X[:, :, 0] = kwargs['x']
                self.X[:, :, 1] = kwargs['y']
                self.X[:, :, 2] = kwargs['z']
                self.nDim = 3
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = numpy.zeros((
                        kwargs['x'].shape[0], kwargs['x'].shape[1], 2))
                self.X[:, :, 0] = kwargs['x']
                self.X[:, :, 1] = kwargs['y']
                self.nDim = 2
            elif 'x' in kwargs:
                self.X = numpy.zeros((
                        kwargs['x'].shape[0], kwargs['x'].shape[1], 1))
                self.X[:, :, 0] = kwargs['x']
                self.nDim = 1
            # enf if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]

            self.ku = geo_utils.checkInput(kwargs['ku'], 'ku', int, 0)
            self.kv = geo_utils.checkInput(kwargs['kv'], 'kv', int, 0)

            if 'Nctlu' in kwargs and 'Nctlv' in kwargs:
                self.Nctlu = geo_utils.checkInput(
                    kwargs['Nctlu'], 'Nctlu', int, 0)
                self.Nctlv = geo_utils.checkInput(
                    kwargs['Nctlv'], 'Nctlv', int, 0)
                self.interp = False
            else:
                self.Nctlu = self.Nu
                self.Nctlv = self.Nv
                self.interp = True
                
            self.orig_data = True

            # Sanity Check on Inputs
            if self.Nctlu >= self.Nu:
                self.Nctlu = self.Nu
            if self.Nctlv >= self.Nv:
                self.Nctlv = self.Nv

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku:
                self.ku = self.Nu
            if self.Nv < self.kv:
                self.kv = self.Nv
            if self.Nctlu < self.ku:
                self.ku = self.Nctlu
            if self.Nctlv < self.kv:
                self.kv = self.Nctlv

            if 'niter' in kwargs:
                self.niter = geo_utils.checkInput(
                    kwargs['niter'], 'niter', int, 0)
            else:
                self.niter = 1
            # end if

            if 'u' in kwargs and 'v' in kwargs:
                self.u = geo_utils.checkInput(
                    kwargs['u'], 'u', float, 1, self.Nu)
                self.v = geo_utils.checkInput(
                    kwargs['v'], 'v', float, 1, self.Nv)
                self.u /= self.u[-1]
                self.v /= self.v[-1]
                [self.V, self.U] = numpy.meshgrid(self.v, self.u)
            else:
                if self.nDim == 3:
                    self.u, self.v, self.U, self.V = \
                        self._calcParameterization()
                else:
                    mpiPrint('Automatic parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise.')
                    sys.exit(1)
                # end if
            # end if
            self.umin = 0
            self.umax = 1
            self.vmin = 0
            self.vmax = 1
            self._calcKnots()
            self.coef = numpy.zeros((self.Nctlu, self.Nctlv, self.nDim))
            if recompute:
                self.recompute()

        return

    def recompute(self):
        """Recompute the surface if any data has been modified
         Required:
             None
         Returns:
             None
             """
        vals, row_ptr, col_ind = pyspline.surface_jacobian_wrap(\
            self.U.T, self.V.T, self.tu, self.tv, self.ku, self.kv, 
            self.Nctlu, self.Nctlv)
        
        N = _assembleMatrix(vals, col_ind, row_ptr, (self.Nu*self.Nv, 
                                                     self.Nctlu*self.Nctlv))

        self.coef = numpy.zeros((self.Nctlu, self.Nctlv, self.nDim))
        if self.interp:
            if USE_SCIPY:
                # Factorize once for efficiency
                solve = scipy.sparse.linalg.dsolve.factorized( N )
                for idim in xrange(self.nDim):
                    self.coef[:, :, idim] = \
                        solve(self.X[:, :, idim].flatten()).reshape(
                        [self.Nctlu, self.Nctlv])
                # end for
            else:
                for idim in xrange(self.nDim):
                    self.coef[:, :, idim] = numpy.linalg.solve(
                        N, self.X[:, :, idim].flatten()).reshape(
                        [self.Nctlu, self.Nctlv])
                # end for
            # end if
        else:
            if USE_SCIPY:
                solve = scipy.sparse.linalg.dsolve.factorized(N.transpose()*N)
                for idim in xrange(self.nDim):
                    rhs  = N.transpose()*self.X[:, :, idim].flatten()
                    self.coef[:, :, idim] = \
                        solve(rhs).reshape([self.Nctlu, self.Nctlv])
                # end for
            else:
                NTN = (N.transpose()*N)
                for idim in xrange(self.nDim):
                    rhs = N.transpose()*self.X[:, :, idim].flatten()
                    self.coef[:, :, idim] = numpy.linalg.solve(
                        NTN, rhs).reshape([self.Nctlu, self.Nctlv])
                # end for
            # end if
        # end if
        self._setEdgeCurves()

        return

    def _calcParameterization(self):
        u = numpy.zeros(self.Nu, 'd')
        U = numpy.zeros((self.Nu, self.Nv), 'd')
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for j in xrange(self.Nv): 
            temp = numpy.zeros(self.Nu, 'd')

            for i in xrange(self.Nu-1):
                temp[i+1] = temp[i] + geo_utils.e_dist(self.X[i, j], self.X[i+1, j])
            # end for

            if temp[-1] == 0: # We have a singular point
                singular_counter += 1
                temp[:] = 0.0
                U[:, j] = numpy.linspace(0, 1, self.Nu)
            else:
                temp /= temp[-1]
                U[:, j] = temp.copy()
            # end if

            u += temp #accumulate the u-parameter calcs for each j
        # end for 
        u = u/(self.Nv-singular_counter) #divide by the number of 'j's we had
        
        v = numpy.zeros(self.Nv, 'd')
        V = numpy.zeros((self.Nu, self.Nv), 'd')
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for i in xrange(self.Nu): 
            temp = numpy.zeros(self.Nv, 'd')
            for j in xrange(self.Nv-1):
                temp[j+1] = temp[j] + geo_utils.e_dist(self.X[i, j], self.X[i, j+1])
            # end for
            if temp[-1] == 0: #We have a singular point
                singular_counter += 1
                temp[:] = 0.0
                V[i, :] = numpy.linspace(0, 1, self.Nv)
            else:
                temp /= temp[-1]
                V[i, :] = temp.copy()
            #end if 

            v += temp #accumulate the v-parameter calcs for each i
        # end for 
        v = v/(self.Nu-singular_counter) #divide by the number of 'i's we had

        return u, v, U, V

    def _calcKnots(self):
        if self.interp:
            self.tu = pyspline.knots_interp(self.u, numpy.array([], 'd'),
                                            self.ku)
            self.tv = pyspline.knots_interp(self.v, numpy.array([], 'd'), 
                                            self.kv)
        else:
            self.tu = pyspline.knots_lms(self.u, self.Nctlu, self.ku)
            self.tv = pyspline.knots_lms(self.v, self.Nctlv, self.kv)
        # end if
            
        return

    def _setEdgeCurves(self):
        """Create curve spline objects for each of the edges"""
        self.edge_curves[0] = curve(k=self.ku, t=self.tu, 
                                    coef=self.coef[:, 0])
        self.edge_curves[1] = curve(k=self.ku, t=self.tu,
                                    coef=self.coef[:, -1])
        self.edge_curves[2] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[0, :])
        self.edge_curves[3] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[-1, :])

        return

    def getValueCorner(self, corner):
        """Get the value of the spline on corner 
        Requred:
            corner: corner index=0, 1, 2 or 3 
        Returns:
            value: Surface value on corner
            """
        assert corner in [0, 1, 2, 3], 'Error, getValueCorner:\
 Corner must be in range 0->3'
        if corner == 0:
            return self.getValue(self.umin, self.vmin)
        elif corner == 1:
            return self.getValue(self.umax, self.vmin)
        elif corner == 2:
            return self.getValue(self.umin, self.vmax)
        elif corner == 3:
            return self.getValue(self.umax, self.vmax)
        #end if

    def getOrigValueCorner(self, node):
        """ 
        Get the values of the original data on corner.
        Required Arguments:
           node: index of conrner = 0, 1, 2, 3
        Returns:
           value: Value at corner
           """
        assert node in [0, 1, 2, 3] and self.orig_data == True, \
 'Error, getOrigValueCorner: No original data for this surface or node\
 is not in range 0->3'

        if node == 0:
            return self.X[0, 0]
        elif node == 1:
            return self.X[-1, 0]
        elif node == 2:
            return self.X[0, -1]
        elif node == 3:
            return self.X[-1, -1]

    def getOrigValuesEdge(self, edge):
        """Get the two end points and midpoint values of the original data
        on edge.
        Required:
           edge: edge index = 0, 1, 2, 3
        Returns:
            start value, mid point value, end_point value
            """
        assert edge in [0, 1, 2, 3] and self.orig_data == True, \
'Error, getOrigValuesEdge: No original data for this surface or edge is\
 not in range 0->3'

        if edge == 0:
            if numpy.mod(self.Nu, 2) == 1: # Its odd
                mid = (self.Nu-1)/2
                return self.X[0, 0], self.X[mid, 0], self.X[-1, 0]
            else:
                Xmid = 0.5 *(self.X[self.Nu/2, 0] + self.X[self.Nu/2 - 1, 0])
                return self.X[0, 0], Xmid, self.X[-1, 0]
        elif edge == 1:
            if numpy.mod(self.Nu, 2) == 1: # Its odd
                mid = (self.Nu-1)/2
                return self.X[0, -1], self.X[mid, -1], self.X[-1, -1]
            else:
                Xmid = 0.5 *(self.X[self.Nu/2, -1] + self.X[self.Nu/2 - 1, -1])
                return self.X[0, -1], Xmid, self.X[-1, -1]
        elif edge == 2:
            if numpy.mod(self.Nv, 2) == 1: # Its odd
                mid = (self.Nv-1)/2
                return self.X[0, 0], self.X[0, mid], self.X[0, -1]
            else:
                Xmid = 0.5 *(self.X[0, self.Nv/2] + self.X[0, self.Nv/2 - 1])
                return self.X[0, 0], Xmid, self.X[0, -1]
        elif edge == 3:
            if numpy.mod(self.Nv, 2) == 1: # Its odd
                mid = (self.Nv-1)/2
                return self.X[-1, 0], self.X[-1, mid], self.X[-1, -1]
            else:
                Xmid = 0.5 *(self.X[-1, self.Nv/2] + self.X[-1, self.Nv/2 - 1])
                return self.X[-1, 0], Xmid, self.X[-1, -1]
        # end if

    def getValueEdge(self, edge, s):
        return self.edge_curves[edge](s)

    def _getBasisPt(self, u, v, vals, istart, col_ind, l_index):
        # This function should only be called from pyGeo
        # The purpose is to compute the basis function for 
        # a u, v point and add it to pyGeo's global dPt/dCoef
        # matrix. vals, row_ptr, col_ind is the CSR data and 
        # l_index in the local -> global mapping for this 
        # surface
        return pyspline.getbasisptsurface(u, v, self.tu, self.tv, self.ku, 
                                          self.kv, vals, col_ind, istart,
                                          l_index.T)
                                
    def __call__(self, u, v):
        """
        Equivalant to getValue
        """
        return self.getValue(u, v)

    def insertKnot(self, direction, s, r):
        """
        Insert at knot in tu at u
        Required:
            direction : 'u' or 'v' The parametric direction for the split
            u : Parametric position to split at
            r : The number of times to insert the knot
        Returns: r, the number of times the knot was actually added
        """
        assert direction in ['u','v'], 'pySpline.surface.splitSurace: direction must be one of \
\'u\' or \'v\''

        s = geo_utils.checkInput(s, 's', float, 0)
        r = geo_utils.checkInput(r, 'r', int, 0)
        if s <= 0.0:
            return 
        if s >= 1.0:
            return
        
        # This is relatively inefficient, but we'll do it for
        # simplicity just call insertknot for each slice in the
        # v-direction:
        
        if direction == 'u':
            # Insert once to know how many times it was actually inserted
            # so we know how big to make the new coef:
            actual_r, t_new, coef_new, break_pt = pyspline.insertknot(
                s, r, self.tu, self.ku, self.coef[:,0].T)

            new_coef = numpy.zeros((self.Nctlu + actual_r, self.Nctlv, self.nDim))
            
            for j in xrange(self.Nctlv):
                actual_r, t_new, coef_slice, break_pt = pyspline.insertknot(
                    s, r, self.tu, self.ku, self.coef[:,j].T)
                new_coef[:,j] = coef_slice[:,0:self.Nctlu+actual_r].T

            self.tu = t_new[0:self.Nctlu+self.ku+actual_r]
            self.Nctlu = self.Nctlu + actual_r

        elif direction == 'v':
            actual_r, t_new, coef_new, break_pt = pyspline.insertknot(
                s, r, self.tv, self.kv, self.coef[0,:].T)

            new_coef = numpy.zeros((self.Nctlu, self.Nctlv+actual_r,self.nDim))

            for i in xrange(self.Nctlu):
                actual_r, t_new, coef_slice, break_pt = pyspline.insertknot(
                    s, r, self.tv, self.kv, self.coef[i,:].T)
                new_coef[i,:] = coef_slice[:,0:self.Nctlv+actual_r].T

            self.tv = t_new[0:self.Nctlv+self.kv+actual_r]
            self.Nctlv = self.Nctlv + actual_r
        # end if

        self.coef = new_coef

        # break_pt is converted to zero based ordering here!!!
        return actual_r, break_pt-1

    def splitSurface(self, direction, t):
        """
        Split surface into two surfaces at parametric location u
        Required: 
            direction : 'u' or 'v' The parametric direction for the split
            t : Parametric position to split at

        Returns: surf1 and surf2, two pyspline surfaces. surf1
        is the lower part and surf2 is the upper part. 
        """

        assert direction in ['u','v'], 'pySpline.surface.splitSurace: direction must be one of \
\'u\' or \'v\''

        # Special case the bounds: (same for both directions)
        if t<= 0.0:
            return None, surface(tu=self.tu.copy(), tv=self.tv.copy(),
                                 ku=self.ku, kv=self.kv, coef=self.coef.copy())
        if t >= 1.0:
            return surface(tu=self.tu.copy(), tv=self.tv.copy(),
                           ku=self.ku, kv=self.kv, coef=self.coef.copy()), None
        
        if direction == 'u':

            r, break_pt = self.insertKnot(direction, t, self.ku-1)

            if r == self.ku -1:
                # A new knot was inserted, so index of relevant place
                # in knot vector is
                break_pt = break_pt + 1
            else:
                break_pt = break_pt
            # end if

            tt = self.tu[break_pt]
            # Process knot vectors:
            t1 = numpy.hstack((self.tu[0:break_pt+self.ku-1].copy(),tt))/tt
            t2 = (numpy.hstack((tt,self.tu[break_pt:].copy()))-tt)/(1.0-tt)
            
            coef1 = self.coef[0:break_pt,:,:].copy()
            coef2 = self.coef[break_pt-1:,:,:].copy()

            return \
                surface(tu=t1, tv=self.tv, ku=self.ku, kv=self.kv, coef=coef1),\
                surface(tu=t2, tv=self.tv, ku=self.ku, kv=self.kv, coef=coef2)
        elif direction == 'v':

            r, break_pt = self.insertKnot(direction, t, self.kv-1)

            if r == self.kv-1:
                # A new knot was inserted, so index of relevant place
                # in knot vector is
                break_pt = break_pt + 1
            else:
                break_pt = break_pt
            # end if

            tt = self.tv[break_pt]
            # Process knot vectors:
            t1 = numpy.hstack((self.tv[0:break_pt+self.kv-1].copy(),tt))/tt
            t2 = (numpy.hstack((tt,self.tv[break_pt:].copy()))-tt)/(1.0-tt)
            
            coef1 = self.coef[:, 0:break_pt,:].copy()
            coef2 = self.coef[:, break_pt-1:,:].copy()

            return \
                surface(tu=self.tu, tv=t1, ku=self.ku, kv=self.kv, coef=coef1),\
                surface(tu=self.tu, tv=t2, ku=self.ku, kv=self.kv, coef=coef2)
        # end if

    def windowSurface(self, uvLow, uvHigh):
        """Create a surface that is windowed by the rectangular
        parametric range defined by uvLow and uvHigh. 
        Required: 

            uvLow: array or list of length two defining u-v coordintes
            of bottom left corner of parametric box

            uvHihg: array of list of length two definign u-v coordines
            of top right corner of parametric box

        Returns: pySpline surface defined on box interior

        """

        # Do u-low split:
        dummy, surf = self.splitSurface('u',uvLow[0])

        # Do u-high split (and re-normalize the split coordinate)
        surf, dummy = surf.splitSurface('u',(uvHigh[0]-uvLow[0])/(1.0-uvLow[0]))

        # Do v-low split:
        dummy, surf = surf.splitSurface('v',uvLow[1])

        # Do v-high split (and re-normalize the split coordinate)
        surf, dummy = surf.splitSurface('v',(uvHigh[1]-uvLow[1])/(1.0-uvLow[1]))

        return surf
        
    def getValue(self, u, v):
        """Get the value at the surface point(s) u, v
        Required Arguments:
           u, v: u and v can be a scalar, vector or matrix of values
        Returns:
           values: An array of size (shape(u), nDim)
           """
        
        u = numpy.array(u).T
        v = numpy.array(v).T
        assert u.shape == v.shape, 'Error, getValue: u and v must have\
 the same shape'
        
        vals = pyspline.eval_surface(numpy.atleast_2d(u), numpy.atleast_2d(v),
                                     self.tu, self.tv, self.ku, self.kv, 
                                     self.coef.T)
        return vals.squeeze().T
                                
    def getDerivative(self, u, v):
        """Get the derivative at the surface point(s) u, v
        Required Arguments:
           u, v: u, v can be a scalar, vector or matrix of values
        Returns:
           values: An array of size (shape(u), 2, ndim)
           """
        u = numpy.array(u)
        v = numpy.array(v)
        assert u.shape == v.shape, 'Error, getDerivative: u and v\
 must have the same shape'
        assert numpy.rank(u) == 0, 'Error: getDerivative only accepts\
 scalar arguments'
        deriv = pyspline.eval_surface_deriv(u, v, self.tu, self.tv, self.ku,
                                            self.kv, self.coef.T)

        return deriv.T
        
    def getSecondDerivative(self, u, v):
        """Get the second derivative matrix at point(s) u, v
        [ (d^2)/(du^2)    (d^2)/(dudv) ]
        [ (d^2)/(dudv)    (d^2)/(dv^2) ]
        Required Arguments:
            u, v: u, v can be a scalar, vector or matrix of values
        Returns:
           values: An array of size (shape(u), 2, 2, ndim)
           """
        u = numpy.array(u)
        v = numpy.array(v)
        assert u.shape == v.shape, 'Error, getSecondDerivative: u and v\
 must have the same shape'
        assert numpy.rank(u) == 0, 'Erorr getSecondDerivative only accepts\
 scalar arguments'

        deriv = pyspline.eval_surface_deriv2(u, v, self.tu, self.tv, self.ku, 
                                             self.kv, self.coef.T)
        
        return deriv.T

    def getBounds(self):
        """Determine the extents of the surface
        Required: 
            None:
        Returns:
            xmin, xmax: xmin is the lowest x, y, z point and xmax the highest
            """
        assert self.nDim == 3, 'getBounds is only defined for nDim = 3'
        cx = self.coef[:, :, 0].flatten()
        cy = self.coef[:, :, 1].flatten()
        cz = self.coef[:, :, 2].flatten()

        Xmin = numpy.array([min(cx), min(cy), min(cz)])
        Xmax = numpy.array([max(cx), max(cy), max(cz)])

        return Xmin, Xmax
    
    def projectPoint(self, x0, Niter=25, eps=1e-10, **kwargs):
        """
        Project a point x0 onto the surface and return parametric position
        curve: 
        Required Arguments:
            x0   : A point or points in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            u    : Initial guess for u position
            v    : Initial guess for v position
        Returns:
            u    : Parametric position(s) u on surface
            v    : Parametric position(s) v on surface
            D    : Distance(s) between point(s) and surface(u, v)
        """

        x0 = numpy.atleast_2d(x0)
        if 'u' in kwargs and 'v' in kwargs:
            u = numpy.atleast_1d(kwargs['u'])
            v = numpy.atleast_1d(kwargs['v'])
        else:
            u = -1*numpy.ones(len(x0))
            v = -1*numpy.ones(len(x0))
        # end if

        assert len(x0) == len(u) == len(v), 'surface projectPoint: The length of x0\
 and u,v,w must be the same'

        # If necessary get brute-force starting point
        if numpy.any(u<0) or numpy.any(u>1) or numpy.any(v<0):
            self._computeData()
            u,v = pyspline.point_surface_start(
                x0.T, self.udata, self.vdata, self.data.T)
        # end if
       
        D = numpy.zeros_like(x0)
        for i in xrange(len(x0)):
            u[i], v[i], D[i] = pyspline.point_surface(
                x0[i], self.tu, self.tv, self.ku, self.kv, self.coef.T, 
                Niter, eps, u[i], v[i])

        return u.squeeze(), v.squeeze(), D.squeeze()

    def projectCurve(self, in_curve, Niter=25, eps=1e-10, **kwargs):
        """
        Find the minimum distance between this surface and a curve
        Required Arguments:
            curve: A pyspline curve class to do the projection with
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            u    : Initial guess for u parameter on surface
            v    : Initial guess for v parameter on surface
            s    : Initial guess for s parameter on curve
        Returns:
            u    : Parametric position u on surface
            v    : Parametric position v on surface
            s    : Parametric position s on curve
            D    : Distance between curve(s) and surface(u, v)
            """      
        u = -1.0
        v = -1.0
        s = -1.0
        if 'u' in kwargs: 
            u = geo_utils.checkInput(kwargs['u'], 'u', float, 0)
        if 'v' in kwargs: 
            v = geo_utils.checkInput(kwargs['v'], 'v', float, 0)
        if 's' in kwargs: 
            s = geo_utils.checkInput(kwargs['s'], 's', float, 0)
        Niter = geo_utils.checkInput(Niter, 'Niter', int, 0)
        eps   = geo_utils.checkInput(eps, 'eps', float, 0)

        # If necessary get brute-force starting point
        if numpy.any(u<0) or numpy.any(u>1) or numpy.any(v<0):
            self._computeData()
            in_curve._computeData()
            s, u, v = pyspline.curve_surface_start(
                in_curve.data.T, in_curve.sdata, self.data.T, 
                self.udata, self.vdata)
        # end if

        return pyspline.curve_surface(\
            in_curve.t, in_curve.k, in_curve.coef.T, self.tu, self.tv, \
                self.ku, self.kv, self.coef.T, Niter, eps, u, v, s)
   
    def _computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """
        # We will base the data on interpolated greville points

        if self.data is None:
            self.edge_curves[0]._calcInterpolatedGrevillePoints()
            self.udata = self.edge_curves[0].sdata
            self.edge_curves[2]._calcInterpolatedGrevillePoints()
            self.vdata = self.edge_curves[2].sdata
            [V, U] = numpy.meshgrid(self.vdata, self.udata)
            self.data = self.getValue(U,V)

        # end if
    
    def _writeTecplotSurface(self, handle):
        """Output this surface\'s data to a open file handle \'handle\' """
        
        self._computeData()
        writeTecplot2D(handle, 'interpolated', self.data)
        
        return
                    
    def _writeDirections(self, handle, isurf):
        """Write out and indication of the surface direction"""
        if self.Nctlu >= 3 and self.Nctlv >= 3:
            data = numpy.zeros((4, self.nDim))
            data[0] = self.coef[1, 2]
            data[1] = self.coef[1, 1]
            data[2] = self.coef[2, 1]
            data[3] = self.coef[3, 1]
            writeTecplot1D(handle, 'surface%d direction'%(isurf), data)
        else:
            mpiPrint('Not Enough control points to output direction indicator')
        #end if

        return 

    def writeTecplot(self, file_name, surfs=True, coef=True, orig=True, dir=False):
        """Write the surface to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            surfs: Boolean to write interpolated surfaces (default=True)
            coef : Boolean to write coefficients (default=True)
            orig : Boolean to write original data (default=True)
            dir  : Boolean to write out surface direction indicators 
                   (default=False)
            """
        f = openTecplot(file_name, self.nDim)

        if surfs:
            self._writeTecplotSurface(f)
        if coef:
            writeTecplot2D(f, 'control_pts', self.coef)
        if orig and self.orig_data:
            writeTecplot2D(f, 'orig_data', self.X)
        if dir:
            self._writeDirections(f, 0)
        closeTecplot(f)

    def _writeIGES_directory(self, handle, Dcount, Pcount):
        """
        Write the IGES file header information (Directory Entry Section)
        for this surface
        """
        # A simplier Calc based on cmlib definations The 13 is for the
        # 9 parameters at the start, and 4 at the end. See the IGES
        # 5.3 Manual paraEntries = 13 + Knotsu + Knotsv + Weights +
        # control points
        assert self.nDim == 3, 'Must have 3 dimensions to write to IGES file'
        paraEntries = 13 + (len(self.tu)) + (len(self.tv)) + \
            self.Nctlu*self.Nctlv + 3*self.Nctlu*self.Nctlv+1

        paraLines = (paraEntries-10) / 3 + 2
        if numpy.mod(paraEntries-10, 3) != 0:
            paraLines += 1

        handle.write('     128%8d       0       0       1       0       0       000000001D%7d\n'%(Pcount, Dcount))
        handle.write('     128       0       2%8d       0                               0D%7d\n'%(paraLines, Dcount+1))
        Dcount += 2
        Pcount += paraLines

        return Pcount , Dcount

    def _writeIGES_parameters(self, handle, Pcount, counter):
        """
        Write the IGES parameter information for this surface
        """

        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'\
                         %(128, self.Nctlu-1, self.Nctlv-1, \
                               self.ku-1, self.kv-1, Pcount, counter))
        counter += 1
        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'\
                         %(0, 0, 1, 0, 0, Pcount, counter))
        counter += 1
        pos_counter = 0

        for i in xrange(len(self.tu)):
            pos_counter += 1
            handle.write('%20.12g,'%(numpy.real(self.tu[i])))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(len(self.tv)):
            pos_counter += 1
            handle.write('%20.12g,'%(numpy.real(self.tv[i])))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(self.Nctlu*self.Nctlv):
            pos_counter += 1
            handle.write('%20.12g,'%(1.0))
            if numpy.mod(pos_counter, 3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount, counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for 

        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                for idim in xrange(3):
                    pos_counter += 1
                    handle.write('%20.12g,'%(numpy.real(
                                self.coef[i, j, idim])))
                    if numpy.mod(pos_counter, 3) == 0:
                        handle.write('  %7dP%7d\n'%(Pcount, counter))
                        counter += 1
                        pos_counter = 0
                    # end if
                # end for
            # end for
        # end for
                        
        # Ouput the ranges
        for  i in xrange(4):
            pos_counter += 1
            if i == 0:
                handle.write('%20.12g,'%(numpy.real(self.umin)))
            if i == 1:
                handle.write('%20.12g,'%(numpy.real(self.umax)))
            if i == 2:
                handle.write('%20.12g,'%(numpy.real(self.vmin)))
            if i == 3:
                # semi-colon for the last entity
                handle.write('%20.12g;'%(numpy.real(self.vmax)))
            if numpy.mod(pos_counter, 3)==0:
                handle.write('  %7dP%7d\n'%(Pcount, counter))
                counter += 1
                pos_counter = 0
            else: # We have to close it up anyway
                if i == 3:
                    for j  in xrange(3-pos_counter):
                        handle.write('%21s'%(' '))
                    # end for
                    pos_counter = 0
                    handle.write('  %7dP%7d\n'%(Pcount, counter))
                    counter += 1
                # end if
            # end if
        # end for

        Pcount += 2

        return Pcount, counter

    def writeTin(self, handle):
        '''Write the pySpline surface to an open handle in .tin format'''
        handle.write('bspline\n')

        # Sizes and Order
        handle.write('%d,%d,%d,%d,0\n'%(self.Nctlu, self.Nctlv, self.ku, self.kv))

        # U - Knot Vector
        for i in xrange(len(self.tu)):
            handle.write('%16.12g,\n'%(self.tu[i]))
        # end for

        # V - Knot Vector
        for j in xrange(len(self.tv)):
            handle.write('%16.12g,\n'%(self.tv[j]))
        # end for

        # Control points:
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                handle.write('%16.12g,%16.12g,%16.12g\n'%(self.coef[i, j, 0], 
                                                          self.coef[i, j, 1], 
                                                          self.coef[i, j, 2]))
            # end for
        # end for

        return 

class volume(object):

    def __init__(self, recompute=True, **kwargs):
        """
        Create an instance of a b-spline surface. There are two
        ways to initialize the class

        Creation: Create an instance of the volume spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            kw, integer: Order for spline in w
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            tw, real array: Knot vector for w
            coef, real array size(Nctlu, Nctlv, Nctlw, nDim): Array of
            control points

        LMS/Interpolation

        Create an instance of the Volume spline class by using an
        interpolating (or lms) spline to given data points. The
        following information is required:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            kw, integer: Order for spline in w
            X, real, array, size(len(u), len(v), len(w), nDim): Array
            of data to fit
            u, v, w real, arrays: (OPTIONAL for nDim == 3 ) Arrays of
            u, v, w values


        Optional: recompute: Specifies whether the actual fitting is complted.
                             Default is True


           NODES      |           EDGES          |           FACES
       6             7|             5          |                   
       #-------------#|       #-------------#  |          #-------------#
      /             / |      /|            /|  |         /|            /|
     /             /  |     / |           / |  |        / |           / |
    /             /   |   6/  |         7/  |  |       /  |   1      /  |
   /             /    |   /   |10       /   |11|      /   |      ---------- 5
  /             /     |  /    |    4   /    |  |     /    |        /    |(back)
 #-------------#      | #-------------#     |  |    #-------------#     |
 4             5      | |     |       |     |  |    |     |       |     |
                      | |     |       |     |  |    |     |       |     | <-3 
       2             3| |     |   1   |     |  |2-> |     |       |     |  
       #-------------#| |     #-------|-----#  |    |     #-------|-----#
      /             / | |8   /        |9   /   |4 ----------      |    /
     /             /  | |   /         |   /    |    |   /         |   /
    /             /   | |  /2         |  /3    |    |  /      0   |  /
   /             /    | | /           | /      |    | /           | /
  /             /     | |/            |/       |    |/            |/
 #-------------#      | #-------------#        |    #-------------#
 0             1      |         0              |

  """
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
            
        self.face_surfaces = [None, None, None, None, None, None]
        self.edge_curves = [None, None, None, None, None, None,
                            None, None, None, None, None, None]
        self.data = None
        if 'ku' in kwargs and 'kv' in kwargs and 'kw' in kwargs and \
                'tu' in kwargs and 'tv' in kwargs and 'tw' in kwargs and \
                'coef' in kwargs:
            self.X = None
            self.u = None
            self.v = None
            self.w = None
            self.U = None
            self.V = None
            self.W = None
            self.interp = False

            self.ku = geo_utils.checkInput(kwargs['ku'], 'ku', int, 0)
            self.kv = geo_utils.checkInput(kwargs['kv'], 'kv', int, 0)
            self.kw = geo_utils.checkInput(kwargs['kw'], 'kw', int, 0)
            self.coef = geo_utils.checkInput(kwargs['coef'], 'coef', float, 4)
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.Nctlw = self.coef.shape[2]
            self.nDim  = self.coef.shape[3]
            self.tu = geo_utils.checkInput(kwargs['tu'], 'tu', float, 1, 
                                 self.Nctlu+self.ku)
            self.tv = geo_utils.checkInput(kwargs['tv'], 'tv', float, 1, 
                                 self.Nctlv+self.kv)
            self.tw = geo_utils.checkInput(kwargs['tw'], 'tw', float, 1, 
                                 self.Nctlw+self.kw)

            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.wmin = self.tw[0]
            self.wmax = self.tw[-1]
            self.orig_data = False
            self._setFaceSurfaces()
            self._setEdgeCurves()
            self.faceBCs = [None, None, None, None, None, None]
        else: # We have LMS/Interpolate
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and 'kw' in kwargs and \
                ( 'X' in kwargs or 'x' in kwargs or \
                      ( 'x' in kwargs and 'y' in kwargs) or \
                      ( 'x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
                      'Error: ku, kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate'

            if 'X' in kwargs:
                self.X  = numpy.array(kwargs['X'])
                if len(self.X.shape) == 1:
                    self.nDim = 1
                else:
                    self.nDim = self.X.shape[3]
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                x = geo_utils.checkInput(kwargs['x'], 'x', float, 3)
                y = geo_utils.checkInput(kwargs['y'], 'y', float, 3)
                z = geo_utils.checkInput(kwargs['z'], 'z', float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[2], 3))
                self.X[:, :, :, 0] = x
                self.X[:, :, :, 1] = y
                self.X[:, :, :, 2] = z
                self.nDim = 3
            elif 'x' in kwargs and 'y' in kwargs:
                x = geo_utils.checkInput(x, 'x', float, 3)
                y = geo_utils.checkInput(x, 'y', float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[3], 3))
                self.X[:, :, :, 0] = x
                self.X[:, :, :, 1] = y
                self.nDim = 2
            elif 'x' in kwargs:
                x = geo_utils.checkInput(x, 'x', float, 3)
                self.X = numpy.zeros((x.shape[0], x.shape[1], x.shape[3], 3))
                self.X[:, :, :, 0] = kwargs['x']
                self.nDim = 1
            # enf if

            if 'faceBCs' in kwargs:
                self.faceBCs = kwargs['faceBCs']
            else:
                self.faceBCs = [None, None, None, None, None, None]
            # end if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.Nw = self.X.shape[2]
            self.ku = geo_utils.checkInput(kwargs['ku'], 'ku', int, 0)
            self.kv = geo_utils.checkInput(kwargs['kv'], 'kv', int, 0)
            self.kw = geo_utils.checkInput(kwargs['kw'], 'kw', int, 0)

            if 'Nctlu' in kwargs and 'Nctlv' in kwargs and 'Nctlw' in kwargs:
                self.Nctlu = geo_utils.checkInput(
                    kwargs['Nctlu'], 'Nctlu', int, 0)
                self.Nctlv = geo_utils.checkInput(
                    kwargs['Nctlv'], 'Nctlv', int, 0)
                self.Nctlw = geo_utils.checkInput(
                    kwargs['Nctlw'], 'Nctlw', int, 0)

                self.interp = False
            else:
                self.Nctlu = self.Nu
                self.Nctlv = self.Nv
                self.Nctlw = self.Nw
                self.interp = True
                
            self.orig_data = True


            # Sanity Check on Inputs
            if self.Nctlu >= self.Nu:
                self.Nctlu = self.Nu
            if self.Nctlv >= self.Nv:
                self.Nctlv = self.Nv
            if self.Nctlw >= self.Nw:
                self.Nctlw = self.Nw

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku: 
                self.ku = self.Nu
            if self.Nv < self.kv: 
                self.kv = self.Nv
            if self.Nw < self.kw: 
                self.kw = self.Nw
            if self.Nctlu < self.ku: 
                self.ku = self.Nctlu
            if self.Nctlv < self.kv: 
                self.kv = self.Nctlv
            if self.Nctlw < self.kw: 
                self.kw = self.Nctlw

          

            if 'niter' in kwargs:
                self.niter = kwargs['niter']
            else:
                self.niter = 1
            # end if

            if 'u' in kwargs and 'v' in kwargs and 'w' in kwargs:
                self.u = geo_utils.checkInput(kwargs['u'], 'u', float, 1)
                self.v = geo_utils.checkInput(kwargs['v'], 'v', float, 1)
                self.w = geo_utils.checkInput(kwargs['w'], 'w', float, 1)
            else:
                if self.nDim == 3:
                    self._calcParameterization()
                else:
                    mpiPrint('Automatic parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise.')
                    sys.exit(1)
                # end if
            # end if
            self.umin = 0
            self.umax = 1
            self.vmin = 0
            self.vmax = 1
            self.wmin = 0
            self.wmax = 1
            self._calcKnots()
            self._setCoefSize()

            if recompute:
                self.recompute()
        # end if

    def recompute(self):
        """Recompute the volume if any driving data has been modified
         Required:
             None
         Returns:
             None
             """
        self._setCoefSize()

        vals, row_ptr, col_ind = pyspline.volume_jacobian_wrap(\
            self.U, self.V, self.W, self.tu, self.tv, self.tw, self.ku,\
                self.kv, self.kw, self.Nctlu, self.Nctlv, self.Nctlw)
        
        N = _assembleMatrix(vals, col_ind, row_ptr, 
                             (self.Nu*self.Nv*self.Nw,
                              self.Nctlu*self.Nctlv*self.Nctlw))
        if self.interp:
            if USE_SCIPY:
                # Factorize once for efficiency
                solve = scipy.sparse.linalg.dsolve.factorized( N )
                for idim in xrange(self.nDim):
                    self.coef[:, :, :, idim] =  solve(
                        self.X[:, :, :, idim].flatten()).reshape(
                        [self.Nctlu, self.Nctlv, self.Nctlw])
                # end for
            else:
                for idim in xrange(self.nDim):
                    self.coef[:, :, :, idim] = numpy.linalg.solve(
                        N, self.X[:, :, :, idim].flatten()).reshape(
                        [self.Nctlu, self.Nctlv, self.Nctlw])
                # end for
            # end if
        else:
            if USE_SCIPY:
                solve = scipy.sparse.linalg.dsolve.factorized(N.transpose()*N)
                for idim in xrange(self.nDim):
                    rhs  = N.transpose()*self.X[:, :, :, idim].flatten()
                    self.coef[:, :, :, idim] = solve(rhs).reshape(
                        [self.Nctlu, self.Nctlv, self.Nctlw])
                # end for
            else:
                NTN = (N.transpose()*N)
                for idim in xrange(self.nDim):
                    rhs = N.transpose()*self.X[:, :, idim].flatten()
                    self.coef[:, :, :, idim] = numpy.linalg.solve(
                        NTN, rhs).reshape(
                        [self.Nctlu, self.Nctlv.self.Nctlw])
                # end for
            # end if
        # end if
     
        self._setFaceSurfaces()
        self._setEdgeCurves()
    
        return
    
    def _setCoefSize(self):
        self.coef = numpy.zeros((self.Nctlu, self.Nctlv, self.Nctlw, self.nDim))
        return 

    def _calcParameterization(self):

        S, u, v, w = pyspline.para3d(self.X.T)
        S = S.T
        self.u = u
        self.v = v
        self.w = w

        self.U = numpy.asarray(S[:, :, :, 0], order='c')
        self.V = numpy.asarray(S[:, :, :, 1], order='c')
        self.W = numpy.asarray(S[:, :, :, 2], order='c')

        return

    def _calcKnots(self):
        if self.interp:
            self.tu = pyspline.knots_interp(self.u, numpy.array([], 'd'),
                                            self.ku)
            self.tv = pyspline.knots_interp(self.v, numpy.array([], 'd'),
                                            self.kv)
            self.tw = pyspline.knots_interp(self.w, numpy.array([], 'd'),
                                            self.kw)
        else:
            self.tu = pyspline.knots_lms(self.u, self.Nctlu, self.ku)
            self.tv = pyspline.knots_lms(self.v, self.Nctlv, self.kv)
            self.tw = pyspline.knots_lms(self.w, self.Nctlw, self.kw)
        # end if
            
        return

    def getValueCorner(self, corner):
        """Get the value of the spline on corner 
        Requred:
            corner: corner index=0, 1, 2 or 3 
        Returns:
            value: Surface value on corner
            """
        assert corner in range(0, 8), 'Error, getValueCorner: Corner must\
 be in range 0->8'
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
        # end if

        return val


    def getOrigValueCorner(self, corner):
        """Return the original value on corner corner"""
        assert corner in [0, 1, 2, 3, 4, 5, 6, 7], 'Error: corner must\
 be 0 to 7'
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
        # end if

        return val

    def getOrigValuesFace(self, face):
        """Return an array of length 8*ndim which cooresponds to the
        the four corners (node 0-3) and the four midpoints (on edges
        0-3) for face face
        Required: 
            face: integer (0-5)
        Returns:
            coordinates: size(8, ndim)
            """

        assert face in [0, 1, 2, 3, 4, 5], 'Error: face must be 0 to 5'
        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu-1)/2, (self.Nu-1)/2]
        else:
            midu = [self.Nu/2, self.Nu/2-1]
        # end if

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv-1)/2, (self.Nv-1)/2]
        else:
            midv = [self.Nv/2, self.Nv/2-1]
        # end if

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw-1)/2, (self.Nw-1)/2]
        else:
            midw = [self.Nw/2, self.Nw/2-1]
        # end if

        if   face == 0:
            values = [self.X[0, 0, 0], self.X[-1, 0, 0], self.X[0, -1, 0],
                      self.X[-1, -1, 0], 
                      0.5*(self.X[midu[0], 0, 0 ] + self.X[midu[1], 0, 0]), 
                      0.5*(self.X[midu[0], -1, 0] + self.X[midu[1], -1, 0]), 
                      0.5*(self.X[0, midv[0], 0 ] + self.X[0, midv[1], 0]), 
                      0.5*(self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0])]
        elif face == 1:
            values = [self.X[0, 0, -1], self.X[-1, 0, -1], self.X[0, -1, -1],
                      self.X[-1, -1, -1], 
                      0.5*(self.X[midu[0], 0, -1 ] + self.X[midu[1], 0, -1]), 
                      0.5*(self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1]), 
                      0.5*(self.X[0, midv[0], -1 ] + self.X[0, midv[1], -1]), 
                      0.5*(self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1])]
        elif face == 2:
            values = [self.X[0, 0, 0], self.X[0, -1, 0], self.X[0, 0, -1],
                      self.X[0, -1, -1], 
                      0.5*(self.X[0, midv[0], 0] + self.X[0, midv[1], 0 ]), 
                      0.5*(self.X[0, midv[0], -1] + self.X[0, midv[1], -1]), 
                      0.5*(self.X[0, 0, midw[0] ] + self.X[0, 0, midw[1] ]), 
                      0.5*(self.X[0, -1, midw[0]] + self.X[0, -1, midw[1]])]
        elif face == 3:
            values = [self.X[-1, 0, 0], self.X[-1, -1, 0], self.X[-1, 0, -1],
                      self.X[-1, -1, -1], 
                      0.5*(self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0 ]), 
                      0.5*(self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1]), 
                      0.5*(self.X[-1, 0, midw[0] ] + self.X[-1, 0, midw[1] ]), 
                      0.5*(self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]])]
        elif face == 4:
            values = [self.X[0, 0, 0], self.X[-1, 0, 0], self.X[0, 0, -1],
                      self.X[-1, 0, -1], 
                      0.5*(self.X[midu[0], 0, 0 ] + self.X[midu[1], 0, 0 ]), 
                      0.5*(self.X[midu[0], 0, -1] + self.X[midu[1], 0, -1]), 
                      0.5*(self.X[0, 0, midw[0] ] + self.X[0, 0, midw[1] ]), 
                      0.5*(self.X[-1, 0, midw[0]] + self.X[-1, 0, midw[1]])]
        elif face == 5:
            values = [self.X[0, -1, 0], self.X[-1, -1, 0], self.X[0, -1, -1],
                      self.X[-1, -1, -1], 
                      0.5*(self.X[midu[0], -1, 0 ] + self.X[midu[1], -1, 0 ]), 
                      0.5*(self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1]), 
                      0.5*(self.X[0, -1, midw[0] ] + self.X[0, -1, midw[1] ]), 
                      0.5*(self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]])]
        # end if
        return numpy.array(values)


    def getMidPointEdge(self, edge):

        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu-1)/2, (self.Nu-1)/2]
        else:
            midu = [self.Nu/2, self.Nu/2-1]
        # end if

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv-1)/2, (self.Nv-1)/2]
        else:
            midv = [self.Nv/2, self.Nv/2-1]
        # end if

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw-1)/2, (self.Nw-1)/2]
        else:
            midw = [self.Nw/2, self.Nw/2-1]
        # end if

        if   edge == 0:
            val = (self.X[midu[0], 0, 0] + self.X[midu[1], 0, 0])
        elif edge == 1:
            val = (self.X[midu[0], -1, 0] + self.X[midu[1], -1, 0])
        elif edge == 2:
            val = (self.X[0, midv[0], 0] + self.X[0, midv[1], 0])
        elif edge == 3:
            val = (self.X[-1, midv[0], 0] + self.X[-1, midv[1], 0])
        elif edge == 4:
            val = (self.X[midu[0], 0, -1] + self.X[midu[1], 0, -1])
        elif edge == 5:
            val = (self.X[midu[0], -1, -1] + self.X[midu[1], -1, -1])
        elif edge == 6:
            val = (self.X[0, midv[0], -1] + self.X[0, midv[1], -1])
        elif edge == 7:
            val = (self.X[-1, midv[0], -1] + self.X[-1, midv[1], -1])
        elif edge == 8:
            val = (self.X[0, 0, midw[0]] + self.X[0, 0, midw[1]])
        elif edge == 9:
            val = (self.X[-1, 0, midw[0]] + self.X[-1, 0, midw[1]])
        elif edge == 10:
            val = (self.X[0, -1, midw[0]] + self.X[0, -1, midw[1]])
        elif edge == 11:
            val = (self.X[-1, -1, midw[0]] + self.X[-1, -1, midw[1]])

        return val


    def getMidPointFace(self, face):
        """Get the midpoint of the face
        on edge.
        Required:
           edge: face index = 0, 1, 2, 3, 4, 5
        Returns:
            midpoint 
            """
        assert face in [0, 1, 2, 3, 4, 5] and self.orig_data == True, 'Error,\
 getMidPointFace: No original data for this surface or face is not in\
 range 0->5'

        if numpy.mod(self.Nu, 2) == 1:
            midu = [(self.Nu-1)/2, (self.Nu-1)/2]
        else:
            midu = [self.Nu/2, self.Nu/2-1]
        # end if

        if numpy.mod(self.Nv, 2) == 1:
            midv = [(self.Nv-1)/2, (self.Nv-1)/2]
        else:
            midv = [self.Nv/2, self.Nv/2-1]
        # end if

        if numpy.mod(self.Nw, 2) == 1:
            midw = [(self.Nw-1)/2, (self.Nw-1)/2]
        else:
            midw = [self.Nw/2, self.Nw/2-1]
        # end if

        if   face == 0:
            val = 0.25*(
                self.X[midu[0], midv[0], 0] + self.X[midu[1], midv[0], 0] +
                self.X[midu[0], midv[1], 0] + self.X[midu[1], midv[1], 0])
        elif face == 1:
            val = 0.25*(
                self.X[midu[0], midv[0], -1] + self.X[midu[1], midv[0], -1] +
                self.X[midu[0], midv[1], -1] + self.X[midu[1], midv[1], -1])
        elif face == 2:
            val = 0.25*(
                self.X[0, midv[0], midw[0]] + self.X[0, midv[1], midw[0]] +
                self.X[0, midv[0], midw[1]] + self.X[0, midv[1], midw[1]])
        elif face == 3:
            val = 0.25*(
                self.X[-1, midv[0], midw[0]] + self.X[-1, midv[1], midw[0]] +
                self.X[-1, midv[0], midw[1]] + self.X[-1, midv[1], midw[1]])
        elif face == 4:
            val = 0.25*(
                self.X[midu[0], 0, midw[0]] + self.X[midu[1], 0, midw[0]] +
                self.X[midu[0], 0, midw[1]] + self.X[midu[1], 0, midw[1]])
        elif face == 5:
            val = 0.25*(
                self.X[midu[0], -1, midw[0]] + self.X[midu[1], -1, midw[0]] +
                self.X[midu[0], -1, midw[1]] + self.X[midu[1], -1, midw[1]])
        # end if
        return val
 
    def _setFaceSurfaces(self):
        """Create face spline objects for each of the faces"""
      
        self.face_surfaces[0] = surface(
            ku=self.ku, kv=self.kv, tu=self.tu, tv=self.tv, 
            coef=self.coef[:, :, 0, :])
        self.face_surfaces[1] = surface(
            ku=self.ku, kv=self.kv, tu=self.tu, tv=self.tv, 
            coef=self.coef[:, :, -1, :])
        self.face_surfaces[2] = surface(
            ku=self.ku, kv=self.kw, tu=self.tu, tv=self.tw, 
            coef=self.coef[:, 0, :, :])
        self.face_surfaces[3] = surface(
            ku=self.ku, kv=self.kw, tu=self.tu, tv=self.tw, 
            coef=self.coef[:, -1, :, :])
        self.face_surfaces[4] = surface(
            ku=self.kv, kv=self.kw, tu=self.tv, tv=self.tw, 
            coef=self.coef[0, :, :, :])
        self.face_surfaces[5] = surface(
            ku=self.kv, kv=self.kw, tu=self.tv, tv=self.tw, 
            coef=self.coef[-1, :, :, :])
        
        return

    def _setEdgeCurves(self):
        """Create edge spline objects for each edge"""
      
        self.edge_curves[0] = curve(k=self.ku, t=self.tu, 
                                    coef=self.coef[:, 0, 0])
        self.edge_curves[1] = curve(k=self.ku, t=self.tu,
                                    coef=self.coef[:, -1, 0])
        self.edge_curves[2] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[0, :, 0])
        self.edge_curves[3] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[-1, :, 0])
        self.edge_curves[4] = curve(k=self.ku, t=self.tu, 
                                    coef=self.coef[:, 0, -1])
        self.edge_curves[5] = curve(k=self.ku, t=self.tu, 
                                    coef=self.coef[:, -1, -1])
        self.edge_curves[6] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[0, :, -1])
        self.edge_curves[7] = curve(k=self.kv, t=self.tv, 
                                    coef=self.coef[-1, :, -1])
        self.edge_curves[8] = curve(k=self.kw, t=self.tw, 
                                    coef=self.coef[0, 0, :])
        self.edge_curves[9] = curve(k=self.kw, t=self.tw, 
                                    coef=self.coef[-1, 0, :])
        self.edge_curves[10] = curve(k=self.kw, t=self.tw, 
                                     coef=self.coef[0, -1, :])
        self.edge_curves[11] = curve(k=self.kw, t=self.tw, 
                                     coef=self.coef[-1, -1, :])

        return 

    def _getBasisPt(self, u, v, w, vals, istart, col_ind, l_index):
        # This function should only be called from pyBlock The purpose
        # is to compute the basis function for a u, v, w point and add
        # it to pyBlcoks's global dPt/dCoef
        # matrix. vals, row_ptr, col_ind is the CSR data and l_index in
        # the local -> global mapping for this volume

        return pyspline.getbasisptvolume(u, v, w, self.tu, self.tv, self.tw, 
                                         self.ku, self.kv, self.kw, vals, 
                                         col_ind, istart, l_index.T)
    
    def __call__(self, u, v, w):
        """
        Equivalant to getValue
        """
        return self.getValue(u, v, w)

    def getValue(self, u, v, w):
        """Get the value at the volume points(s) u, v, w
        Required Arguments:
            u, v, w: u, w and w can be a scalar, vector or matrix of values
        Returns:
           values: An array of size (shape(u), nDim)
           """
        u = numpy.atleast_3d(u).T
        v = numpy.atleast_3d(v).T
        w = numpy.atleast_3d(w).T

        assert u.shape == v.shape == w.shape, 'Error, getValue: u and v must\
 have the same shape'

        vals = pyspline.eval_volume(u, v, w, self.tu, self.tv, self.tw, 
                                    self.ku, self.kv, self.kw, self.coef.T)
        return vals.squeeze().T

    def getValueEdge(self, edge, s):
        """Get the value at the volume points(s) u, v, w
        Required Arguments:
            u, v, w: u, w and w can be a scalar, vector or matrix of values
        Returns:
           values: An array of size (shape(u), nDim)
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

        assert u.shape == v.shape == w.shape, 'Error, getValue: u and v\
 must have the same shape'

        vals = pyspline.eval_volume(u, v, w, self.tu, self.tv, self.tw, 
                                    self.ku, self.kv, self.kw, self.coef.T)
        return vals.squeeze().T

    def getBounds(self):
        """Determine the extents of the volume
        Required: 
            None:
        Returns:
            xmin, xmax: xmin is the lowest x, y, z point and xmax the highest
            """
        assert self.nDim == 3, 'getBounds is only defined for nDim = 3'
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

    def projectPoint(self, x0, Niter=25, eps=1e-10, **kwargs):
        """
        Project a point x0 onto the volume and return parametric position

        Required Arguments:
            x0   : A point or points in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps  : Tolerance
            u    : Initial guess for u position
            v    : Initial guess for v position
            w    : Initial guess for w position 
        Returns:
            u    : Parametric position(s) u on surface
            v    : Parametric position(s) v on surface
            w    : Parametric position(s) w on surface
            D    : Distance(s) between point(s) and surface(u, v, w)
        """
      
        x0 = numpy.atleast_2d(x0)
        if 'u' in kwargs and 'v' in kwargs and 'w' in kwargs:
            u = numpy.atleast_1d(kwargs['u'])
            v = numpy.atleast_1d(kwargs['v'])
            w = numpy.atleast_1d(kwargs['w'])
        else:
            u = -1*numpy.ones(len(x0))
            v = -1*numpy.ones(len(x0))
            w = -1*numpy.ones(len(x0))
        # end if

        assert len(x0) == len(u) == len(v) == len(w), 'volume projectPoint: The length of x0\
 and u,v,w must be the same'

        # If necessary get brute-force starting point
        if numpy.any(u<0) or numpy.any(u>1) or numpy.any(v<0) or numpy.any(v>1):
            self._computeData()
            u,v,w = pyspline.point_volume_start(x0.T, self.udata, self.vdata, self.wdata, 
                                                self.data.T)
        # end if
        D = numpy.zeros_like(x0)
        for i in xrange(len(x0)):
            u[i], v[i], w[i], D[i] = pyspline.point_volume(
                x0[i], self.tu, self.tv, self.tw, self.ku, self.kv, self.kw, self.coef.T, 
                Niter, eps, u[i], v[i], w[i])

        return u.squeeze(), v.squeeze(), w.squeeze(), D.squeeze()

    def _computeData(self):
        """
        Compute discrete data that is used for the Tecplot
        Visualization as well as the data for doing the brute-force
        checks
        """
        # We will base the data on interpolated greville points

        if self.data is None:
            self.edge_curves[0]._calcInterpolatedGrevillePoints()
            self.udata = self.edge_curves[0].sdata
            self.edge_curves[2]._calcInterpolatedGrevillePoints()
            self.vdata = self.edge_curves[2].sdata
            self.edge_curves[8]._calcInterpolatedGrevillePoints()
            self.wdata = self.edge_curves[8].sdata
            U = numpy.zeros((len(self.udata),len(self.vdata),len(self.wdata)))
            V = numpy.zeros((len(self.udata),len(self.vdata),len(self.wdata)))
            W = numpy.zeros((len(self.udata),len(self.vdata),len(self.wdata)))
            for i in xrange(len(self.udata)):
                for j in xrange(len(self.vdata)):
                    for k in xrange(len(self.wdata)):
                        U[i,j,k] = self.udata[i]
                        V[i,j,k] = self.vdata[j]
                        W[i,j,k] = self.wdata[k]
                    # end for
                # end for
            # end for
            self.data = self.getValue(U,V,W)
        # end if

        return
        
    def _writeTecplotVolume(self, handle):
        """Output this volume\'s data to a open file handle \'handle\' """

        self._computeData()
        writeTecplot3D(handle, 'interpolated', self.data)
        
        return

    def writeTecplot(self, file_name, vols=True, coef=True, orig=False):

        """Write the surface to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            surfs: Boolean to write interpolated surfaces (default=True)
            coef : Boolean to write coefficients (default=True)
            orig : Boolean to write original data (default=True)
            dir  : Boolean to write out surface direction indicators 
                   (default=False)

            """
        f = openTecplot(file_name, self.nDim)
        if vols:
            self._writeTecplotVolume(f)
        if coef:
            writeTecplot3D(f, 'control_pts', self.coef)
        if orig and self.orig_data:
            writeTecplot3D(f, 'orig_data', self.X)
        closeTecplot(f)

        return

    def _writeBvol(self, handle, binary):
        """Write the data to an open file handle"""
        # Initial Integers
        init_values = numpy.array([self.Nctlu, self.Nctlv, self.Nctlw, 
                             self.ku, self.kv, self.kw])
        if binary:
            init_values.tofile(handle, sep="")

            # Knot Vectors
            self.tu.tofile(handle, sep="")
            self.tv.tofile(handle, sep="")
            self.tw.tofile(handle, sep="")

            # Control Points
            self.coef.flatten().tofile(handle)
        else:
            init_values.tofile(handle, sep="\n", format="%d")
            handle.write('\n')

            # Knot Vectors
            self.tu.tofile(handle, sep="\n", format="%f")
            handle.write('\n')
            self.tv.tofile(handle, sep="\n", format="%f")
            handle.write('\n')
            self.tw.tofile(handle, sep="\n", format="%f")
            handle.write('\n')

            # Control Points
            self.coef.flatten().tofile(handle, sep="\n", format="%f")
            handle.write('\n')
        # end if
        
        return 

# ----------------------------------------------------------------------
#                     Misc Helper Functions
# ----------------------------------------------------------------------

def trilinear_volume(*args):
    """This is a short-cut function to create a trilinear b-spline volume
    Args can contain:
        X: array of size(8, 3) which contains the corners of the box is
        coordinate order (i, j, k)
      
        xmin, xmax: The lower extreme and upper extreme corners of the box
        
        """
    tu = [0, 0, 1, 1]
    tv = [0, 0, 1, 1]
    tw = [0, 0, 1, 1]
    ku = 2
    kv = 2
    kw = 2

    if len(args) == 1:
        return volume(coef=args[0], tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    elif len(args) == 2:
        xmin = numpy.array(args[0]).astype('d')
        xmax = numpy.array(args[1]).astype('d')

        x_low  = xmin[0]
        x_high = xmax[0]
        y_low  = xmin[1]
        y_high = xmax[1]
        z_low  = xmin[2]
        z_high = xmax[2]

        coef = numpy.zeros((2, 2, 2, 3))
        coef[0, 0, 0, :] = [x_low, y_low, z_low]
        coef[1, 0, 0, :] = [x_high, y_low, z_low]
        coef[0, 1, 0, :] = [x_low, y_high, z_low]
        coef[1, 1, 0, :] = [x_high, y_high, z_low]
        coef[0, 0, 1, :] = [x_low, y_low, z_high]
        coef[1, 0, 1, :] = [x_high, y_low, z_high]
        coef[0, 1, 1, :] = [x_low, y_high, z_high]
        coef[1, 1, 1, :] = [x_high, y_high, z_high]

        return volume(coef=coef, tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    else:
        mpiPrint('Error: An unknown number of arguments was passed to\
 trilinear  volume')
        sys.exit(1)
    # end if
    
    return

def bilinear_surface(*args):
    """This is short-cut function to create a bilinear surface
    Args can contain:
        x: array of size(4, 3) The four corners of the array arranged in
        the coordinate direction orientation:

        2          3
        /----------\
        |          |
        |          |
        |          |
        \----------/
        0          1
    
   OR

   Args can contain pt1, pt2, pt3, pt4 is CCW Ordering

        3          2
        /----------\
        |          |
        |          |
        |          |
        \----------/
        0          1
        """
    if len(args) == 1:
        # One argument passed in ... assume its X
        assert len(args[0]) == 4, 'Error: a single argument passed to bilinear\
 surface must contain 4 points and be of size (4, 3)'
        coef = numpy.zeros((2, 2, 3))
        coef[0, 0] = args[0][0]
        coef[1, 0] = args[0][1]
        coef[0, 1] = args[0][2]
        coef[1, 1] = args[0][3]
        return surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)
    else:
        # Assume 4 arguments
        coef = numpy.zeros([2, 2, 3])
        coef[0, 0] = args[0]
        coef[1, 0] = args[1]
        coef[0, 1] = args[3]
        coef[1, 1] = args[2]
        return surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)
    # end if

def line(*args, **kwargs):
    """This is a short cut function to create a line curve

    Args can contain:
       X: array of size(2, ndim) The two end points
       
       OR:
       x1, x2: The two end points (each of size ndim)

       OR:
       x1, dir=direction
       x1 and the keyword argument direction

       OR: 
       x1, dir=direction, length=length
       x1, direction and specific length
       """
    if len(args) == 2:
        # Its a two-point type
        return curve(coef=[args[0], args[1]], k=2, t=[0, 0, 1, 1])
    elif len(args) == 1:
        if len(args[0]) == 2: # its X
            return curve(coef=args[0], k=2, t=[0, 0, 1, 1])
        elif 'dir' in kwargs:
            # We have point and direction
            if 'length' in kwargs:
                x2 = args[0] + kwargs['dir']/numpy.linalg.norm(
                    kwargs['dir'])*kwargs['length']
            else:
                x2 = args[0] + kwargs['dir']
            # end if
            return curve(coef=[args[0], x2], k=2, t=[0, 0, 1, 1])
        else:
            mpiPrint('Error: dir must be specified if only 1 argument is given')
        # end if
    # end if

#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
    print 'There are two examples in the example directory.'
    print 'Look at test_curve.py and test_surf.py for more information'

