'''
pySpline

Contains an class functions for working with B-spline curves and surfaces

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


# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, time, copy, pdb

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,zeros,sqrt,array,reshape,meshgrid,mod,floor,\
    ones,vstack,real,where,arange,append,hstack,mgrid
from numpy.linalg import norm

import scipy
from scipy import sparse,io
try:
    from scipy.sparse.linalg.dsolve import factorized # Version 0.8 of scipy
except:
    from scipy.linsolve import factorized

# =============================================================================
# Custom Python modules
# =============================================================================

import pyspline
from mdo_import_helper import *
USE_TECIO = pyspline.tecplot_test()

# =============================================================================
def e_dist(x1,x2):
    '''Get the eculidean distance between two points in nDim space'''
    total = 0.0
    for idim in xrange(len(x1)):
        total += (x1[idim]-x2[idim])**2
    # end for
    return sqrt(total)


def _writeTecplot1D(handle,name,data):
    '''A Generic write tecplot zone to file'''
    if handle != None:
        nx = data.shape[0]
        ndim = data.shape[1]
        handle.write('Zone T=\"%s\" I=%d\n'%(name,nx))
        handle.write('DATAPACKING=POINT\n')
        for i in xrange(nx):
            for idim in xrange(ndim):
                handle.write('%f '%(data[i,idim]))
            # end for
            handle.write('\n')
        # end for
    else:
        pyspline.i_ordered(name,data)
    # end if
    return

def _writeTecplot2D(handle,name,data):
    '''A Generic write tecplot zone to file'''
    if handle != None:
        nx = data.shape[0]
        ny = data.shape[1]
        ndim = data.shape[2]
        handle.write('Zone T=\"%s\" I=%d J=%d\n'%(name,nx,ny))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(ny):
            for i in xrange(nx):
                for idim in xrange(ndim):
                    handle.write('%f '%(data[i,j,idim]))
                # end for
                handle.write('\n')
            # end for
        # end for
    else:
        pyspline.ij_ordered(name,data)
    # end if
    return

def _writeTecplot3D(handle,name,data):
    '''A Generic write tecplot zone to file'''
    if handle != None:
        nx = data.shape[0]
        ny = data.shape[1]
        nz = data.shape[2]
        ndim = data.shape[3]
        handle.write('Zone T=\"%s\" I=%d J=%d K=%d\n'%(name,nx,ny,nz))
        handle.write('DATAPACKING=POINT\n')
        for k in xrange(nz):
            for j in xrange(ny):
                for i in xrange(nx):
                    for idim in xrange(ndim):
                        handle.write('%f '%(data[i,j,k,idim]))
                    # end for
                    handle.write('\n')
                # end for 
            # end for
        # end for
    else:
        pyspline.ijk_ordered(name,data)
    # end if
    return

def _writeHeader(f,ndim):
    if ndim == 1:
        f.write ('VARIABLES = "X"\n')
    elif ndim == 2:
        f.write ('VARIABLES = "X","Y"\n')
    else:
        f.write ('VARIABLES = "X", "Y","Z"\n')
    # end if

def openTecplot(file_name,ndim,tecio=USE_TECIO):
    mpiPrint(' ')
    if tecio:
        mpiPrint('Opening binary Tecplot File: %s'%(file_name))
        pyspline.open_tecplot(file_name,ndim)
        f = None
    else:
        mpiPrint('Opening ascii Tecplot File: %s'%(file_name))
        f = open(file_name,'w')
        _writeHeader(f,ndim)
    # end if
    return f

def closeTecplot(f):
    if f == None:
        pyspline.close_tecplot()
    else:
        f.close()
    #end if
    
    return 

# =============================================================================
# pySpline classes
# =============================================================================
  
class surface(object):

    def __init__(self,recompute=True,*args,**kwargs):
        '''
        Create an instance of a b-spline surface. There are two
        ways to initialize the class

        Creation: Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control points

        LMS/Interpolation
        Create an instance of the spline class by using an interpolating spline to given data points. **kwarg
        or a LMS spline. The following information is required:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            X, real, array, size(len(u),len(v),nDim): Array of data to fit 
              --> OR  x=<vals>, (len(u),len(v) -2D interpolation
              --> OR  x=<vals>,y=<vals>, each of size (len(u),len(v) 
              --> OR  x=<vals>,y=<vals>,z=<vals> each of len(s) ->3D parametric surface
            u,v, real, arrays: (OPTIONAL for nDim == 3 ) Arrays of u and v values
            
        For LMS spline:
            Nctlu -> (integer) number of control points in u
            Nctlv -> (integer) number of control points in v

        Optional Data:
           u    : array of u parameter locations (optional for ndim == 3)
           v    : array of v parameter locations (optional for ndim == 3)
           niter:  The number of Hoschek\'s parameter corrections to run
           '''
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
        
        if 'ku' in kwargs and 'kv' in kwargs and 'tu' in kwargs and 'tv' in \
                kwargs and 'coef' in kwargs:
            self.X = None
            self.u = None
            self.v = None
            self.U = None
            self.V = None
            self.ku = int(kwargs['ku'])
            self.kv = int(kwargs['kv'])
            self.tu = array(kwargs['tu'],'d')
            self.tv = array(kwargs['tv'],'d')
            self.coef = array(kwargs['coef'],'d')
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.nDim = self.coef.shape[2]
            self.orig_data = False
            self._setEdgeCurves()
            return
        else: # We have LMS/Interpolate
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and \
                ( 'X' in kwargs or 'x' in kwargs or \
                      ( 'x' in kwargs and 'y' in kwargs) or \
                      ( 'x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
                      'Error: ku,kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate'

            if 'X' in kwargs:
                self.X  = array(kwargs['X'])
                if len(self.X.shape) == 1:
                    self.nDim =1
                else:
                    self.nDim = self.X.shape[2]
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = zeros((kwargs['x'].shape[0],kwargs['x'].shape[1],3))
                self.X[:,:,0] = kwargs['x']
                self.X[:,:,1] = kwargs['y']
                self.X[:,:,2] = kwargs['z']
                self.nDim = 3
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = zeros((kwargs['x'].shape[0],kwargs['x'].shape[1],2))
                self.X[:,:,0] = kwargs['x']
                self.X[:,:,1] = kwargs['y']
                self.nDim = 2
            elif 'x' in kwargs:
                self.X = zeros((kwargs['x'].shape[0],kwargs['x'].shape[1],1))
                self.X[:,:,0] = kwargs['x']
                self.nDim = 1
            # enf if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]

            self.ku = int(kwargs['ku'])
            self.kv = int(kwargs['kv'])


            if 'Nctlu' in kwargs and 'Nctlv' in kwargs:
                self.Nctlu = kwargs['Nctlu']
                self.Nctlv = kwargs['Nctlv']
                self.interp = False
            else:
                self.Nctlu = self.Nu
                self.Nctlv = self.Nv
                self.interp = True
                
            self.orig_data = True

            # Sanity Check on Inputs
            if self.Nctlu >= self.Nu:  self.Nctlu = self.Nu
            if self.Nctlv >= self.Nv:  self.Nctlv = self.Nv

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku: self.ku = self.Nu
            if self.Nv < self.kv: self.kv = self.Nv
            if self.Nctlu < self.ku: self.ku = self.Nctlu
            if self.Nctlv < self.kv: self.kv = self.Nctlv

            if 'niter' in kwargs:
                self.niter = kwargs['niter']
            else:
                self.niter = 1
            # end if

            if 'u' in kwargs and 'v' in kwargs:
                self.u = array(kwargs['u'],'f')
                self.v = array(kwargs['v'],'f')
                self.u/=self.u[-1]
                self.v/=self.v[-1]
                [self.V,self.U] = meshgrid(self.v,self.u)
            else:
                if self.nDim == 3:
                    self.u,self.v,self.U,self.V = self._calcParameterization()
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
            self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim))
            if recompute:
                self.recompute()

        return

    def recompute(self):
        '''Recompute the surface if any data has been modified
         Required:
             None
         Returns:
             None
             '''
        vals,row_ptr,col_ind = pyspline.surface_jacobian_wrap(\
            self.U,self.V,self.tu,self.tv,self.ku,self.kv,self.Nctlu,self.Nctlv)
        N = sparse.csr_matrix((vals,col_ind,row_ptr),
                             [self.Nu*self.Nv,self.Nctlu*self.Nctlv])
        if self.interp:
            solve = factorized( N )
            for idim in xrange(self.nDim):
                self.coef[:,:,idim] =  solve(self.X[:,:,idim].flatten()).reshape([self.Nctlu,self.Nctlv])
            # end for
        else:
            solve = factorized (N.transpose()*N)
            for idim in xrange(self.nDim):
                rhs  = N.transpose()*self.X[:,:,idim].flatten()
                self.coef[:,:,idim] = solve(rhs).reshape([self.Nctlu,self.Nctlv])
            # end for
        # end if
        self._setEdgeCurves()

        return

    def _calcParameterization(self):
        u = zeros(self.Nu,'d')
        U = zeros((self.Nu,self.Nv),'d')
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for j in xrange(self.Nv): 
            temp = zeros(self.Nu,'d')

            for i in xrange(self.Nu-1):
                temp[i+1] = temp[i] + e_dist(self.X[i,j],self.X[i+1,j])
            # end for

            if temp[-1] == 0: # We have a singular point
                singular_counter += 1
                temp[:] = 0.0
                U[:,j] = linspace(0,1,self.Nu)
            else:
                temp /= temp[-1]
                U[:,j] = temp.copy()
            # end if

            u += temp #accumulate the u-parameter calcs for each j
        # end for 
        u =u/(self.Nv-singular_counter) #divide by the number of 'j's we had
        
        v = zeros(self.Nv,'d')
        V = zeros((self.Nu,self.Nv),'d')
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for i in xrange(self.Nu): 
            temp = zeros(self.Nv,'d')
            for j in xrange(self.Nv-1):
                temp[j+1] = temp[j] + e_dist(self.X[i,j],self.X[i,j+1])
            # end for
            if temp[-1] == 0: #We have a singular point
                singular_counter += 1
                temp[:] = 0.0
                V[i,:] = linspace(0,1,self.Nv)
            else:
                temp /= temp[-1]
                V[i,:] = temp.copy()
            #end if 

            v += temp #accumulate the v-parameter calcs for each i
        # end for 
        v = v/(self.Nu-singular_counter) #divide by the number of 'i's we had

        return u,v,U,V

    def _calcKnots(self):
        if self.interp:
            self.tu = pyspline.knots_interp(self.u,array([],'d'),self.ku)
            self.tv = pyspline.knots_interp(self.v,array([],'d'),self.kv)
        else:
            self.tu = pyspline.knots_lms(self.u,self.Nctlu,self.ku)
            self.tv = pyspline.knots_lms(self.v,self.Nctlv,self.kv)
        # end if
            
        return

    def _setEdgeCurves(self):
        '''Create curve spline objects for each of the edges'''
        self.edge_curves = [None,None,None,None]
        self.edge_curves[0] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,0])
        self.edge_curves[1] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,-1])
        self.edge_curves[2] = curve('create',k=self.kv,t=self.tv,coef=self.coef[0,:])
        self.edge_curves[3] = curve('create',k=self.kv,t=self.tv,coef=self.coef[-1,:])

        return

    def getValueCorner(self,corner):
        '''Get the value of the spline on corner 
        Requred:
            corner: corner index=0,1,2 or 3 
        Returns:
            value: Surface value on corner
            '''
        assert corner in [0,1,2,3],'Error, getValueCorner: Corner must be in range 0->3'
        if corner == 0:
            return self.getValue(self.umin,self.vmin)
        elif corner == 1:
            return self.getValue(self.umax,self.vmin)
        elif corner == 2:
            return self.getValue(self.umin,self.vmax)
        elif corner ==3:
            return self.getValue(self.umax,self.vmax)
        #end if

    def getOrigValuesEdge(self,edge):
        '''Get the two end points and midpoint values of the original data
        on edge.
        Required:
           edge: edge index = 0,1,2,3
        Returns:
            start value, mid point value, end_point value
            '''
        assert edge in [0,1,2,3] and self.orig_data == True,'Error, getOrigValuesEdge: No \
original data for this surface or edge is not in range 0->3'

        if edge == 0:
            if mod(self.Nu,2) == 1: # Its odd
                mid = (self.Nu-1)/2
                return self.X[0,0],self.X[mid,0],self.X[-1,0]
            else:
                Xmid = 0.5 *(self.X[self.Nu/2,0] + self.X[self.Nu/2 - 1,0])
                return self.X[0,0],Xmid,self.X[-1,0]
        elif edge == 1:
            if mod(self.Nu,2) == 1: # Its odd
                mid = (self.Nu-1)/2
                return self.X[0,-1],self.X[mid,-1],self.X[-1,-1]
            else:
                Xmid = 0.5 *(self.X[self.Nu/2,-1] + self.X[self.Nu/2 - 1,-1])
                return self.X[0,-1],Xmid,self.X[-1,-1]
        elif edge == 2:
            if mod(self.Nv,2) == 1: # Its odd
                mid = (self.Nv-1)/2
                return self.X[0,0],self.X[0,mid],self.X[0,-1]
            else:
                Xmid = 0.5 *(self.X[0,self.Nv/2] + self.X[0,self.Nv/2 - 1])
                return self.X[0,0],Xmid,self.X[0,-1]
        elif edge == 3:
            if mod(self.Nv,2) == 1: # Its odd
                mid = (self.Nv-1)/2
                return self.X[-1,0],self.X[-1,mid],self.X[-1,-1]
            else:
                Xmid = 0.5 *(self.X[-1,self.Nv/2] + self.X[-1,self.Nv/2 - 1])
                return self.X[-1,0],Xmid,self.X[-1,-1]
        # end if

    def getOrigValueCorner(self,node):
        ''' 
        Get the values of the original data on corner.
        Required Arguments:
           node: index of conrner = 0,1,2,3
        Returns:
           value: Value at corner
           '''
        assert node in [0,1,2,3] and self.orig_data == True,'Error, getOrigValueCorner: No \
original data for this surface or node is not in range 0->3'

        if node == 0:
            return self.X[0,0]
        elif node == 1:
            return self.X[-1,0]
        elif node == 2:
            return self.X[0,-1]
        elif node == 3:
            return self.X[-1,-1]

    def _getBasisPt(self,u,v,vals,istart,col_ind,l_index):
        # This function should only be called from pyGeo
        # The purpose is to compute the basis function for 
        # a u,v point and add it to pyGeo's global dPt/dCoef
        # matrix. vals,row_ptr,col_ind is the CSR data and 
        # l_index in the local -> global mapping for this 
        # surface
        return pyspline.getbasisptsurface(u,v,self.tu,self.tv,self.ku,self.kv,vals,
                                          col_ind,istart,l_index)
                                
    def __call__(self,u,v):
        '''
        Equivalant to getValue
        '''
        return self.getValue(u,v)

    def getValue(self,u,v):
        '''Get the value at the surface point(s) u,v
        Required Arguments:
           u,v: u and v can be a scalar,vector or matrix of values
        Returns:
           values: An array of size (shape(u), nDim)
           '''
        u = array(u)
        v = array(v)
        assert u.shape == v.shape,'Error, getValue: u and v must have the same shape'
        if len(u.shape) == 0:
            return pyspline.eval_surface(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 1:
            return pyspline.eval_surface_v(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 2:
            return pyspline.eval_surface_m(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
                                
    def getDerivative(self,u,v):
        '''Get the derivative at the surface point(s) u,v
        Required Arguments:
           u,v: u,v can be a scalar,vector or matrix of values
        Returns:
           values: An array of size (shape(u), 2, ndim)
           '''
        u = array(u)
        v = array(v)
        assert u.shape == v.shape,'Error, getDerivative: u and v must have the same shape'
        if len(u.shape) == 0:
            return pyspline.eval_surface_deriv(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 1:
            return pyspline.eval_surface_deriv_v(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 2:
            return pyspline.eval_surface_deriv_m(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        
    def getNormal(self,u,v):
        '''Get the normalized normal at the surface point(s) u,v
        Required Arguments:
           u: u can be a scalar,vector or matrix of values
           v: v can be a scalar,vector or matrix of values
        Returns:
           values: An array of size (shape(u), 3)
           '''
        u = array(u)
        v = array(v)
        assert self.nDim == 3,'Error, getNormal is only defined for three spatial dimensions'
        assert u.shape == v.shape,'Error, getNormal: u and v must have the same shape'
        if len(u.shape) == 0:
            return pyspline.eval_surf_normal(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 1:
            return pyspline.eval_surf_normal_v(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 2:
            return pyspline.eval_surf_normal_m(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        # end if

    def getSecondDerivative(self,u,v):
        '''Get the second derivative matrix at point(s) u,v
        [ (d^2)/(du^2)    (d^2)/(dudv) ]
        [ (d^2)/(dudv)    (d^2)/(dv^2) ]
        Required Arguments:
            u,v: u,v can be a scalar,vector or matrix of values
        Returns:
           values: An array of size (shape(u), 2,2,ndim)
           '''
        u = array(u)
        v = array(v)
        assert u.shape == v.shape,'Error, getSecondDerivative: u and v must have the same shape'
        if len(u.shape) == 0:
            return pyspline.eval_surface_deriv2(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 1:
            return pyspline.eval_surface_deriv2_v(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        elif len(u.shape) == 2:
            return pyspline.eval_surface_deriv2_m(u,v,self.tu,self.tv,self.ku,self.kv,self.coef)
        # end if
 
    def getBounds(self):
        '''Determine the extents of the surface
        Required: 
            None:
        Returns:
            xmin,xmax: xmin is the lowest x,y,z point and xmax the highest
            '''
        assert self.nDim == 3,'getBounds is only defined for nDim = 3'
        cx = self.coef[:,:,0].flatten()
        cy = self.coef[:,:,1].flatten()
        cz = self.coef[:,:,2].flatten()

        Xmin = zeros(self.nDim)
        Xmin[0] = min(cx)
        Xmin[1] = min(cy)
        Xmin[2] = min(cz)

        Xmax = zeros(self.nDim)
        Xmax[0] = max(cx)
        Xmax[1] = max(cy)
        Xmax[2] = max(cz)

        return Xmin,Xmax
    
    def projectPoint(self,x0,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''
        Project a point x0 onto the surface and return parametric position
        curve: 
        Required Arguments:
            x0   : A point or points in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps1 : Intersection/Relative change tolerance
            eps2 : Cosine convergence measure tolerance
            u    : Initial guess for u position
            v    : Initial guess for v position
        Returns:
            u    : Parametric position(s) u on surface
            v    : Parametric position(s) v on surface
            D    : Distance(s) between point(s) and surface(u,v)
        '''
        # We will use a starting point u0,v0 if given
        x0 = array(x0)

        if len(x0.shape) == 1: # Just 1 value passed in 
            if 'u' in kwargs:
                u = array([kwargs['u']])
            else:
                u = array([-1.0])
            # end if
            if 'v' in kwargs:
                v = array([kwargs['v']])
            else:
                v = array([-1.0])
            # end if
            result = pyspline.point_surface(\
                array([x0]),self.tu,self.tv,self.ku,self.kv,self.coef,Niter,eps1,eps2,u,v)
            return result[0][0],result[1][0],result[2][0] #u,v,D
        else:
            if 'u' in kwargs:
                u = kwargs['u']
            else:
                u = -1.0*ones(len(x0))
            # end if
            if 'v' in kwargs:
                v = kwargs['v']
            else:
                v = -1.0*ones(len(x0))
            # end if
            result = pyspline.point_surface(\
                array([x0]),self.tu,self.tv,self.ku,self.kv,self.coef,Niter,eps1,eps2,u,v)
            return result[0],result[1],result[2] #u,v,D
        # end if

    def projectCurve(self,curve,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''
        Find the minimum distance between this surface and a curve
        Required Arguments:
            curve: A pyspline curve class to do the projection with
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps1 : Intersection/Relative change tolerance
            eps2 : Cosine convergence measure tolerance
            u    : Initial guess for u parameter on surface
            v    : Initial guess for v parameter on surface
            s    : Initial guess for s parameter on curve
        Returns:
            u    : Parametric position u on surface
            v    : Parametric position v on surface
            s    : Parametric position s on curve
            D    : Distance between curve(s) and surface(u,v)
            '''      
        u = -1.0
        v = -1.0
        s = -1.0
        if 'u' in kwargs: u = kwargs['u']
        if 'v' in kwargs: v = kwargs['v']
        if 's' in kwargs: s = kwargs['s']
        return pyspline.curve_surface(\
            curve.t,curve.k,curve.coef,self.tu,self.tv,\
                self.ku,self.kv,self.coef,Niter,eps1,eps2,u,v,s)
   
    def _writeTecplotOrigData(self,handle):
        if self.orig_data:
            _writeTecplot2D(handle,'orig_data',self.X)
        # end if

        return

    def _writeTecplotSurface(self,handle):
        '''Output this surface\'s data to a open file handle \'handle\' '''
        
        Nx = self.Nctlu*(self.ku)*2
        Ny = self.Nctlv*(self.kv)*2
        u_plot = 0.5*(1-cos(linspace(0,pi,Nx)))
        v_plot = 0.5*(1-cos(linspace(0,pi,Ny)))
        [V_plot,U_plot] = meshgrid(v_plot,u_plot)
        X = self.getValue(U_plot,V_plot)
        _writeTecplot2D(handle,'interpolated',X)
        
        return
                    
    def _writeTecplotCoef(self,handle):
        '''Write the Spline coefficients to handle'''
        _writeTecplot2D(handle,'control_pts',self.coef)

        return

    def _writeTecplotEdge(self,handle,edge,*args,**kwargs):
        '''Dump out a linear zone along edge used for visualizing edge connections'''
        self.edge_curves[edge].writeTecplotEdge(handle,*args,**kwargs)

        return

    def _writeDirections(self,handle,isurf):
        '''Write out and indication of the surface direction'''
        if self.Nctlu >= 3 and self.Nctlv >=3:
            data = zeros((4,self.nDim))
            data[0] = self.coef[1,2]
            data[1] = self.coef[1,1]
            data[2] = self.coef[2,1]
            data[3] = self.coef[3,1]
            _writeTecplot1D(handle,'surface%d direction'%(isurf),data)
        else:
            mpiPrint('Not Enough control points to output direction indicator')
        #end if

        return 

    def writeTecplot(self,file_name,surfs=True,coef=True,orig=True,dir=False,tecio=USE_TECIO):
        '''Write the surface to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            surfs: Boolean to write interpolated surfaces (default=True)
            coef : Boolean to write coefficients (default=True)
            orig : Boolean to write original data (default=True)
            dir  : Boolean to write out surface direction indicators (default=False)
            '''
        f = openTecplot(file_name,self.nDim,tecio)

        if surfs:
            self._writeTecplotSurface(f)
        if coef:
            self._writeTecplotCoef(f)
        if orig:
            self._writeTecplotOrigData(f)
        if dir:
            self._writeDirections(f)
        closeTecplot(f)

    def _writeIGES_directory(self,handle,Dcount,Pcount):
        '''
        Write the IGES file header information (Directory Entry Section)
        for this surface
        '''
        # A simplier Calc based on cmlib definations The 13 is for the
        # 9 parameters at the start, and 4 at the end. See the IGES
        # 5.3 Manual paraEntries = 13 + Knotsu + Knotsv + Weights +
        # control points
        assert self.nDim == 3, 'Must have 3 dimensions to write to IGES file'
        paraEntries = 13 + (len(self.tu)) + (len(self.tv)) + \
            self.Nctlu*self.Nctlv + 3*self.Nctlu*self.Nctlv+1

        paraLines = (paraEntries-10) / 3 + 2
        if mod(paraEntries-10,3) != 0: paraLines += 1

        handle.write('     128%8d       0       0       1       0       0       000000001D%7d\n'%(Pcount,Dcount))
        handle.write('     128       0       2%8d       0                               0D%7d\n'%(paraLines,Dcount+1))
        Dcount += 2
        Pcount += paraLines

        return Pcount , Dcount

    def _writeIGES_parameters(self,handle,Pcount,counter):
        '''
        Write the IGES parameter information for this surface
        '''

        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'\
                         %(128,self.Nctlu-1,self.Nctlv-1,\
                               self.ku-1,self.kv-1,Pcount,counter))
        counter += 1
        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'\
                         %(0,0,1,0,0,Pcount,counter))
        counter += 1
        pos_counter = 0

        for i in xrange(len(self.tu)):
            pos_counter += 1
            handle.write('%20.12g,'%(real(self.tu[i])))
            if mod(pos_counter,3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(len(self.tv)):
            pos_counter += 1
            handle.write('%20.12g,'%(real(self.tv[i])))
            if mod(pos_counter,3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(self.Nctlu*self.Nctlv):
            pos_counter += 1
            handle.write('%20.12g,'%(1.0))
            if mod(pos_counter,3) == 0:
                handle.write('  %7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for 

        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                for idim in xrange(3):
                    pos_counter += 1
                    handle.write('%20.12g,'%(real(self.coef[i,j,idim])))
                    if mod(pos_counter,3) == 0:
                        handle.write('  %7dP%7d\n'%(Pcount,counter))
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
                handle.write('%20.12g,'%(real(self.umin)))
            if i == 1:
                handle.write('%20.12g,'%(real(self.umax)))
            if i == 2:
                handle.write('%20.12g,'%(real(self.vmin)))
            if i == 3:
                # semi-colon for the last entity
                handle.write('%20.12g;'%(real(self.vmax)))
            if mod(pos_counter,3)==0:
                handle.write('  %7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            else: # We have to close it up anyway
                if i ==3:
                    for j  in xrange(3-pos_counter):
                        handle.write('%21s'%(' '))
                    # end for
                    pos_counter = 0
                    handle.write('  %7dP%7d\n'%(Pcount,counter))
                    counter += 1
                # end if
            # end if
        # end for

        Pcount += 2

        return Pcount,counter

class curve(object):

    def __init__(self,*args,**kwargs):
        '''
        Create an instance of a b-spline curve. There are two
        ways to initialize the class

        Creation: Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            k, integer: Order for spline
            t, real array: Knot vector 
            coef, real array size(Nctl,nDim): Array of control points

        LMS/Interpolation
        Create an instance of the spline class by using an interpolating spline to given data points. **kwarg
        or a LMS spline. The following information is required:

            k, integer: Order for spline
            X, real, array, size(len(s),nDim): Array of data to fit 
              --> OR  x=<vals>, len(s) ->1D interpolation
              --> OR  x=<vals>,y=<vals>, each of len(s) -> 2D parametric curve
              --> OR  x=<vals>,y=<vals>,z=<vals> each of len(s) ->3D parametric curve
            s, real, array: (OPTIONAL for nDim >=2 ) Array of s values 
            
        For LMS spline:
             Nctl -> (integer) number of control points

        Optional Data:
        
        weights:   A array of weighting factors for each fitting point
                   A value of -1 can be used to exactly constrain a point
        deriv and deriv_ptr: deriv is a array which contains derivative
                             information at point indicies defined by deriv_ptr
        deriv_weights: A array of len deriv_ptr which contains derivative weighting
                       A value of -1 can be used to exactly constrain a derivative
        niter:  The number of Hoschek's parameter corrections to run
        '''
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
        
        if 'k' in kwargs and 't' in kwargs and 'coef' in kwargs: # We have a create class
            self.s = None
            self.X = None
            self.N = None
            self.k = int(kwargs['k'] )
            self.t = array(kwargs['t'])
            self.coef = array(kwargs['coef'])
            self.Nctl = self.coef.shape[0]
            self.orig_data = False
            self.smin = self.t[0]
            self.smax = self.t[-1]
            if len(self.coef.shape) == 1:
                self.nDim = 1
            else:
                self.nDim = self.coef.shape[1]
            self._calcGrevillePoints()
        else: #lms or interpolate function
            assert 'k' in kwargs and ('X' in kwargs or 'x' in kwargs),\
                'Error: At least spline orderk, k and X (or x=,y=) MUST be defined for (interpolation) spline creation.\
Nctl=<number of control points> must be specified for a LMS fit'
            if 'X' in kwargs:
                self.X  = array(kwargs['X'])
                if len(self.X.shape) == 1:
                    self.nDim =1
                else:
                    self.nDim = self.X.shape[1]
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = vstack([kwargs['x'],kwargs['y'],kwargs['z']]).T          
                self.nDim = 3
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = vstack([kwargs['x'],kwargs['y']]).T
                self.nDim = 2
            elif 'x' in kwargs:
                self.X = array(kwargs['x']).reshape((len(kwargs['x']),1))
                self.nDim = 1
            # enf if
            self.X = self.X.astype('d') # Make sure its real
            self.k = int(kwargs['k'])
            self.N = len(self.X)
            self.t = None
            if 'niter' in kwargs: 
                self.niter = kwargs['niter']
            else:
                self.niter = 1
            # end if

            if 's' in kwargs:
                self.s = array(kwargs['s'])
            else:
                assert self.nDim > 1,'Error, pySpline: For 1D splines, the basis, s must be given'
                self._getParameterization()
            # end if
            self.smin = self.s[0]
            self.smax = self.s[-1]
            self.orig_data = True
            
            if 'weights' in kwargs:
                self.weights = array(kwargs['weights'])
            else:
                self.weights = ones(self.N)
            # end if

            if 'deriv' in kwargs and 'deriv_ptr' in kwargs:
                self.deriv = array(kwargs['deriv'])
                self.deriv_ptr = array(kwargs['deriv_ptr'])
            else:
                self.deriv = None
                self.deriv_ptr = array([])
            # end if
            if 'deriv_weights' in kwargs and self.deriv != None:
                self.deriv_weights = array(kwargs['deriv_weights'])
            else:
                if self.deriv != None:
                    self.deriv_weights = ones(len(self.deriv))
                else:
                    self.deriv_weights = None
                # end if
            # end if
            if not 'Nctl' in kwargs: # We are doing an interpolation...set all weights to 1
                self.interp = True
                self.weights[:] = 1
                if self.deriv_weights != None:
                    self.deriv_weights[:] = 1
                # end if
            else:
                self.interp = False
                self.Nctl = int(kwargs['Nctl'])
            # end if
            self.recompute(self.niter)
        return

    def recompute(self,niter):
        '''
        Run iterations of Hoschek's Parameter Correction
        Required:
            niter: The number of parameter correction iterations to run
        Returns:
            None
            '''
        # Now do the sparation between the constrained and unconstrained
        su_select = where(self.weights > 0.0)
        sc_select = where(self.weights <= 0.0)
        S  = self.X[su_select]
        su = self.s[su_select]
        T  = self.X[sc_select]
        sc = self.s[sc_select]
        weights = self.weights[where(self.weights > 0.0)]

        nu = len(S)
        nc = len(T)
     
        # And the derivative info
        if self.deriv != None:
            sdu_select = where(self.deriv_weights > 0.0)
            sdc_select = where(self.deriv_weights <= 0.0)
            S = vstack((S,self.deriv[sdu_select]))
            sdu = self.s[self.deriv_ptr][sdu_select]
            T = vstack((T,self.deriv[sdc_select]))
            sdc = self.s[self.deriv_ptr][sdc_select]
            weights = append(weights,self.deriv_weights[where(self.deriv_weights > 0.0)])
            ndu = len(sdu)
            ndc = len(sdc)
        else:
            sdu = array([],'d')
            sdc = array([],'d')
            ndu = 0
            ndc = 0
        # end if
        W = sparse.csr_matrix((weights,arange(len(weights)),arange(len(weights)+1))) # Diagonal
            
        if self.interp:
            self.Nctl = nu+nc+ndu+ndc
            self.niter = 1
        # end if

        # Sanity check to make sure k is ok
        if nu+nc+ndu+ndc < self.k:
            self.k = nu+nc+ndu+ndc
        # end if

        # Generate the knot vector,greville points and empty coefficients
        if self.interp:
            if self.t == None:
                self.t = pyspline.knots_interp(self.s,self.deriv_ptr,self.k)
            # end if
        else:
            if self.t == None:
                self.t = pyspline.knots_lms(self.s,self.Nctl,self.k)
            # end if
        # end if
        self._calcGrevillePoints()
        self.coef = zeros((self.Nctl,self.nDim),'d')

        # Get the 'N' jacobain
        N_vals = zeros((nu+ndu)*self.k)          # |
        N_row_ptr = zeros(nu+ndu+1,'intc')       # | ->  Standard CSR formulation
        N_col_ind = zeros((nu+ndu)*self.k,'intc')# |
        pyspline.curve_jacobian_wrap(su,sdu,self.t,self.k,self.Nctl,N_vals,N_row_ptr,N_col_ind)
        N = sparse.csr_matrix((N_vals,N_col_ind,N_row_ptr),[nu+ndu,self.Nctl])
        if self.interp:
            solve = factorized( N ) # Factorize once for efficiency
            for idim in xrange(self.nDim):
                self.coef[:,idim] = solve(S[:,idim])
            # end for
            return
        # end if
                
        length = pyspline.poly_length(self.X)
        for iter in xrange(niter):
            su = self.s[su_select]
            sc = self.s[sc_select]
            if self.deriv != None:
                sdu = self.s[self.deriv_ptr][sdu_select]
                sdc = self.s[self.deriv_ptr][sdc_select]
            # end if

            pyspline.curve_jacobian_wrap(su,sdu,self.t,self.k,self.Nctl,N_vals,N_row_ptr,N_col_ind)
            NTWN = N.transpose()*W*N # We need this either way

            if nc + ndc == 0: # We are doing LMS but no constraints...just a straight weighted LMS
                solve = factorized( NTWN ) # Factorize once for efficiency
                for idim in xrange(self.nDim):
                    self.coef[:,idim] = solve(N.transpose()*W*S[:,idim])
                # end for
            else:  # Now its more complicated since we have constraints
                M_vals = zeros((nc+ndc)*self.k)          #|
                M_row_ptr = zeros(nc+ndc+1,'intc')       #| -> Standard CSR formulation
                M_col_ind = zeros((nc+ndc)*self.k,'intc')#|
                pyspline.curve_jacobian_wrap(sc,sdc,self.t,self.k,self.Nctl,M_vals,M_row_ptr,M_col_ind)
                M = sparse.csr_matrix((M_vals,M_col_ind,M_row_ptr),[nc+ndc,self.Nctl])
                # Now we must assemble the constrained jacobian
                # [ N^T*W*T      M^T][P] = [ N^T*W*S]
                # [ M            0  ][R]   [ T      ] 

                NTWN = NTWN.tocsr() # Convert to CSR for the constr_jac function
                MT   = M.transpose().tocsr() # Convert to CSR for the constr_jac function
                j_val,j_col_ind,j_row_ptr = pyspline.constr_jac(
                    NTWN.data,NTWN.indptr,NTWN.indices,MT.data,MT.indptr,MT.indices,
                    M.data,M.indptr,M.indices,self.Nctl)
                # Create sparse csr matrix and factorize
                J = sparse.csr_matrix((j_val,j_col_ind,j_row_ptr),[self.Nctl+nc+ndc,self.Nctl+nc+ndc])
                solve = factorized( J )
                for idim in xrange(self.nDim):
                    rhs = hstack((N.transpose()*W*S[:,idim],T[:,idim]))
                    self.coef[:,idim] = solve(rhs)[0:self.Nctl]
                # end for
            # end if (constr - not constrained

            # Run para correction
            self.s = pyspline.curve_para_corr(self.t,self.k,self.s.copy(),self.coef,length,self.X)
        # end for (iter loop)
        # Check the RMS
        rms = 0.0
        for idim in xrange(self.nDim):
            rms += norm(N*self.coef[:,idim]-S[:,idim])**2
        # end for
        rms = sqrt(rms/self.N)
        mpiPrint('Rms is: %f'%(rms),self.NO_PRINT)
      
        return

    def _getParameterization(self):
        # We need to parameterize the curve
        self.s = zeros(self.N,'d');
        for i in xrange(self.N-1):
            dist = 0
            for idim in xrange(self.nDim):
                dist += (self.X[i+1,idim] - self.X[i,idim])**2
                # end for
            self.s[i+1] = self.s[i] + sqrt(dist)
            # end for
        # end for
        self.length = self.s[-1]
        self.s /= self.s[-1]

        return
     
    def insertKnot(self,u,r):
        '''
        Insert at knot at u
        Required:
            u : Parametric position to split at
            r : The number of times to insert the knot
            '''
        t_new,coef_new,break_pt = pyspline.insertknot(u,r,self.t,self.k,self.coef)
        self.t = t_new
        self.coef = coef_new
        self.Nctl = self.Nctl + r
        
        return

    def splitCurve(self,u):
        ''' 
        Split the curve at parametric position
        Required:
            u : Parametric position for split
        Optional: 
            None
        Returns:
            curve1,curve2: The two split spline curve
        '''
        # First determine if u is very close to an existing knot
        ileft,mflag = pyspline.intrv(self.t,u,1)
        if abs(u-self.t[ileft-1]) < 0.005:
            u = self.t[ileft-1]
        # end if

        r,t_new,coef_new,break_pt = pyspline.insertknot(u,self.k-1,self.t,self.k,self.coef)
        # r is the number of time the knot was actually added
        s = self.k-1-r # Back out the multiplicity of the point
        break_pt = break_pt - s
        # -------- Process the Knot Vectors--------
        t1 = zeros(break_pt + self.k)
        t2 = zeros(self.Nctl+2*self.k-break_pt-s)
        t1[0:break_pt] = t_new[0:break_pt]
        t1[break_pt:] = self.k*[u]
        if s== 0:
            t2[1:] = t_new[break_pt:]
        else:
            t2[1:] = t_new[break_pt:-s]
        t2[0]  = u

        # Now Normalize
        t1 = t1/t1[-1]
        t2 = (t2-u)/(1-u)

        # ------- Proces the Coefficient Arrays
        coef1 = zeros((break_pt,self.nDim))
        coef2 = zeros((len(t2)-self.k,self.nDim))
        coef1[:,:] = coef_new[0:break_pt,:]
        if s==0:
            coef2[:,:] = coef_new[break_pt-1:,:]        
        else:
            coef2[:,:] = coef_new[break_pt-1:-s,:]        
        # end if
        curve1 = curve('create',k=self.k,t=t1,coef=coef1)
        curve2 = curve('create',k=self.k,t=t2,coef=coef2)

        return curve1,curve2

    def getLength(self):
        '''
        Compute the length of the curve using the Eculdian Norm
        Required Argument: 
            None
        Returns:
            Length: Lenght of curve
            '''
        points = self.getValue(self.gpts)# These are physical points
        length = 0
        for i in xrange(len(points)-1):
            length += e_dist(points[i],points[i+1])
        # end for

        return length

    def _calcGrevillePoints(self):
        '''Calculate the Greville points'''
        # Interpolate
        self.gpts = zeros(self.Nctl)
        for i in xrange(self.Nctl):
            for n in xrange(self.k-1): #degree loop
                self.gpts[i] += self.t[i+n+1]
            # end for
            self.gpts[i] /= (self.k-1)
        # end for

        return

    def __call__(self,s):
        '''
        Equivilant to getValue
        '''
        return self.getValue(s)
    
    def getValue(self,s):
        '''
        Get the value of the spline
        Required Arguments:
           s: s can be a scalar or vector of values
        Returns:
           values: An array of shape(len(s),nDim)) if s is an array
                   or array of len(nDim) if nDim == 1
                   '''
        s = array(s)
        if len(s.shape) == 0:
            return pyspline.eval_curve(s,self.t,self.k,self.coef)
        elif len(s.shape) == 1:
            return  pyspline.eval_curve_v(s,self.t,self.k,self.coef)
        # end if

    def getDerivative(self,s):
        '''
        Get the value of the derivatve of the spline
        Required Arguments:
           s: s can be a scalar or vector of values
        Returns:
           values: An array of shape(len(s),nDim)) if s is an array
                   or array of len(nDim) if nDim == 1
                   '''
        s = array(s)
        if len(s.shape) == 0:
            return pyspline.eval_curve_deriv(s,self.t,self.k,self.coef)
        elif len(s.shape) == 1:
            return pyspline.eval_curve_deriv_v(s,self.t,self.k,self.coef)

    def projectPoint(self,x0,Niter=20,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''Project a point x0 or (points x0) onto the curve and return parametric
        position(s)
        curve: 
        Required Arguments:
            x0   : A point(s) in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps1 : Intersection/Relative change tolerance
            eps2 : Cosine convergence measure tolerance
            s    : Initial guess for position
        Returns:
            s: Parametric position(s) of closest point on curve
            D: Distance(s) from point to curve
        '''
        x0 = array(x0) # Ensure x0 is array

        if len(x0.shape) == 1: # Just 1 value pased in
            if 's' in kwargs: 
                s=kwargs['s']
            else:
                s=[-1.0]
            # end if
            result = pyspline.point_curve(array([x0]),self.t,self.k,self.coef,
                                          Niter,eps1,eps2,s)
            return result[0][0],result[1][0] # s,D
        else:
            if 's' in kwargs: 
                s=kwargs['s']
            else:
                s=-1.0*ones(len(x0))
            # end if
            result = pyspline.point_curve(x0,self.t,self.k,self.coef,
                                          Niter,eps1,eps2,s)
            return result[0],result[1] # s,D
        # end if

    def projectCurve(self,curve,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''
        Find the minimum distance between this curve (self) and a second
        curve passed in (curve)
        Required Arguments:
            curve: A pyspline curve class to do the projection with
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps1 : Intersection/Relative change tolerance
            eps2 : Cosine convergence measure tolerance
            s    : Initial guess for curve1 (this curve class)
            t    : Initial guess for curve2 (cuve passed in)
        Returns:
            s    : Parametric position on curve1 (this class)
            t    : Parametric position on curve2 (curve passed in)
            D    : Distance between curve1(s) and curve2(t)
            '''
        s = -1.0
        t = -1.0
        if 's' in kwargs:  s = kwargs['s']
        if 't' in kwargs:  t = kwargs['t']
        return pyspline.curve_curve(self.t,self.k,self.coef,
                                    curve.t,curve.k, curve.coef,
                                    Niter,eps1,eps2,s,t)

    def writeTecplot(self,file_name,curve=True,coef=True,orig=True,tecio=USE_TECIO,*args,**kwargs):
        '''
        Write the cuve to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            curve: boolean flag to write the interpolated curve (default=True)
            coef : boolean flag to write the control points (defaut=True)
            orig : boolean flag to write the original data (default=True)
            size : A distance to determine the output resolution of interpolated curve
        Returns:
            None
            '''

        f = openTecplot(file_name,self.nDim,tecio)
        if curve:
            self._writeTecplotCurve(f,size=0.1)
        if coef:
            self._writeTecplotCoef(f)
        if orig:
            self._writeTecplotOrigData(f)
        closeTecplot(f)

        return 

    def _writeTecplotCoef(self,handle,*args,**kwargs):
        '''
        Write the Spline coefficients to handle
        '''
        _writeTecplot1D(handle,'control_pts',self.coef)

        return

    def _writeTecplotOrigData(self,handle,*args,**kwargs):
        '''
        Write the original data to a handle
        '''
        if self.orig_data:
            _writeTecplot1D(handle,'orig_data',self.X)
        # end if

        return

    def _writeTecplotCurve(self,handle,*args,**kwargs):
        '''
        Write the interpolated curve to a handle
        '''
        if 'size' in kwargs:
            length = self.getLength()
            n=int(floor(real(length/kwargs['size'])))
            X = self.getValue(linspace(0,1,n))
        else:
            if not self.s == None:
                X = self.getValue(self.s)
            else:
                s = 0.5*(1-cos(linspace(0,pi,self.k*self.Nctl+2)))
                X = self.getValue(s)
            # end if
        # end if

        _writeTecplot1D(handle,'interpolated',X)

        return 
  
class volume(object):

    def __init__(self,recompute=True,*args,**kwargs):
        '''
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
            coef, real array size(Nctlu,Nctlv,Nctlw,nDim): Array of control points

        LMS/Interpolation
        Create an instance of the Volume spline class by using an interpolating (or lms)
        spline to given data points. The following information is required:

            ku, integer: Order for spline in u
            kv, integer: Order for spline in v
            kw, integer: Order for spline in w
            X, real, array, size(len(u),len(v),len(w),nDim): Array of data to fit 
            u,v,w real, arrays: (OPTIONAL for nDim == 3 ) Arrays of u,v,w values


        Optional: recompute: Specifies whether the actual fitting is complted.
                             Default is True


              NODES       |           EDGES          |           FACES
        6             7 |             5            |                   
        #-------------# |       #-------------#    |          #-------------#
       /             /  |      /|            /|    |         /|            /|
      /             /   |     / |           / |    |        / |           / |
     /             /    |   6/  |         7/  |    |       /  |   1      /  |
    /             /     |   /   |10       /   |11  |      /   |      ---------- 5
   /             /      |  /    |    4   /    |    |     /    |        /    |      
  #-------------#       | #-------------#     |    |    #-------------#     |
  4             5       | |     |       |     |    |    |     |       |     |
                        | |     |       |     |    |    |     |       |     | <-3
        2             3 | |     |   1   |     |    |2-> |     |       |     |  
        #-------------# | |     #-------|-----#    |    |     #-------|-----#
       /             /  | |8   /        |9   /     |4 ----------      |    /
      /             /   | |   /         |   /      |    |   /         |   /
     /             /    | |  /2         |  /3      |    |  /      0   |  /
    /             /     | | /           | /        |    | /           | /
   /             /      | |/            |/         |    |/            |/
  #-------------#       | #-------------#          |    #-------------#
  0             1       |         0                |

  '''
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
            
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
            self.ku = int(kwargs['ku'])
            self.kv = int(kwargs['kv'])
            self.kw = int(kwargs['kw'])
            self.tu = array(kwargs['tu'],'d')
            self.tv = array(kwargs['tv'],'d')
            self.tw = array(kwargs['tv'],'d')
            self.coef = array(kwargs['coef'],'d')
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.Nctlw = self.coef.shape[2]
            self.nDim  = self.coef.shape[3]
            self.umin = self.tu[0]
            self.umax = self.tu[-1]
            self.vmin = self.tv[0]
            self.vmax = self.tv[-1]
            self.wmin = self.tw[0]
            self.wmax = self.tw[-1]
            self.orig_data = False
            self._setFaceSurfaces()
            self._setEdgeCurves()
            self.faceBCs = [None,None,None,None,None,None]
        else: # We have LMS/Interpolate
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and 'kw' in kwargs and \
                ( 'X' in kwargs or 'x' in kwargs or \
                      ( 'x' in kwargs and 'y' in kwargs) or \
                      ( 'x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
                      'Error: ku,kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate'

            if 'X' in kwargs:
                self.X  = array(kwargs['X'])
                if len(self.X.shape) == 1:
                    self.nDim =1
                else:
                    self.nDim = self.X.shape[3]
            elif 'x' in kwargs and 'y' in kwargs and 'z'in kwargs:
                self.X = zeros((kwargs['x'].shape[0],kwargs['x'].shape[1],kwargs['x'].shape[2],3))
                self.X[:,:,:,0] = kwargs['x']
                self.X[:,:,:,1] = kwargs['y']
                self.X[:,:,:,2] = kwargs['z']
                self.nDim = 3
            elif 'x' in kwargs and 'y' in kwargs:
                self.X = zeros((kwargs['x'].shape[0],kwargs['x'].shape[1],2))
                self.X[:,:,:,0] = kwargs['x']
                self.X[:,:,:,1] = kwargs['y']
                self.nDim = 2
            elif 'x' in kwargs:
                self.X = zeros((kwargs['x'].shape[0],1))
                self.X[:,:,:,0] = kwargs['x']
                self.nDim = 1
            # enf if

            if 'faceBCs' in kwargs:
                self.faceBCs = kwargs['faceBCs']
            else:
                self.faceBCs = [None,None,None,None,None,None]
            # end if

            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.Nw = self.X.shape[2]
            self.ku = int(kwargs['ku'])
            self.kv = int(kwargs['kv'])
            self.kw = int(kwargs['kv'])

            if 'Nctlu' in kwargs and 'Nctlv' in kwargs:
                self.Nctlu = kwargs['Nctlu']
                self.Nctlv = kwargs['Nctlv']
                self.Nctlw = kwargs['Nctlw']
                self.interp = False
            else:
                self.Nctlu = self.Nu
                self.Nctlv = self.Nv
                self.Nctlw = self.Nw
                self.interp = True
                
            self.orig_data = True

            # Sanity Check on Inputs
            if self.Nctlu >= self.Nu:  self.Nctlu = self.Nu
            if self.Nctlv >= self.Nv:  self.Nctlv = self.Nv
            if self.Nctlw >= self.Nw:  self.Nctlw = self.Nw

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku: self.ku = self.Nu
            if self.Nv < self.kv: self.kv = self.Nv
            if self.Nw < self.kw: self.kw = self.Nw
            if self.Nctlu < self.ku: self.ku = self.Nctlu
            if self.Nctlv < self.kv: self.kv = self.Nctlv
            if self.Nctlw < self.kv: self.kv = self.Nctlw

            if 'niter' in kwargs:
                self.niter = kwargs['niter']
            else:
                self.niter = 1
            # end if

            if 'u' in kwargs and 'v' in kwargs and 'w' in kwargs:
                self.u = kwargs['u']
                self.v = kwargs['v']
                self.w = kwargs['w']
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
        '''Recompute the volume if any driving data has been modified
         Required:
             None
         Returns:
             None
             '''
        self._setCoefSize()
        
        vals,row_ptr,col_ind = pyspline.volume_jacobian_wrap(\
            self.U,self.V,self.W,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,\
                self.Nctlu,self.Nctlv,self.Nctlw)
        
        N = sparse.csr_matrix((vals,col_ind,row_ptr),
                              [self.Nu*self.Nv*self.Nw,self.Nctlu*self.Nctlv*self.Nctlw])
     
        if self.interp:
            solve = factorized( N )
            for idim in xrange(self.nDim):
                self.coef[:,:,:,idim] =  solve(self.X[:,:,:,idim].flatten()).reshape([self.Nctlu,self.Nctlv,self.Nctlw])
            # end for
            
        else:
            solve = factorized (N.transpose()*N)
            for idim in xrange(self.nDim):
                rhs  = N.transpose()*self.X[:,:,:,idim].flatten()
                self.coef[:,:,:,idim] = solve(rhs).reshape([self.Nctlu,self.Nctlv,self.Nctlw])
            # end for
        # end if
        self._setFaceSurfaces()
        self._setEdgeCurves()
        return
    
    def _setCoefSize(self):
        self.coef = zeros((self.Nctlu,self.Nctlv,self.Nctlw,self.nDim))
        return 

    def _calcParameterization(self):
        S,u,v,w = pyspline.para3d(self.X)
        self.u = u
        self.v = v
        self.w = w
        self.U = S[:,:,:,0]
        self.V = S[:,:,:,1]
        self.W = S[:,:,:,2]

        return

    def _calcKnots(self):
        if self.interp:
            self.tu = pyspline.knots_interp(self.u,array([],'d'),self.ku)
            self.tv = pyspline.knots_interp(self.v,array([],'d'),self.kv)
            self.tw = pyspline.knots_interp(self.w,array([],'d'),self.kw)
        else:
            self.tu = pyspline.knots_lms(self.u,self.Nctlu,self.ku)
            self.tv = pyspline.knots_lms(self.v,self.Nctlv,self.kv)
            self.tw = pyspline.knots_lms(self.w,self.Nctlw,self.kw)
        # end if
            
        return

    def getOrigValueCorner(self,corner):
        '''Return the original value on corner corner'''
        if corner == 0:
            return self.X[0,0,0]
        elif corner == 1:
            return self.X[-1,0,0]
        elif corner == 2:
            return self.X[0,-1,0]
        elif corner == 3:
            return self.X[-1,-1,0]
        elif corner == 4:
            return self.X[0,0,-1]
        elif corner == 5:
            return self.X[-1,0,-1]
        elif corner == 6:
            return self.X[0,-1,-1]
        elif corner == 7:
            return self.X[-1,-1,-1]

    def getOrigValuesFace(self,face):
        '''Return an array of length 8*ndim which cooresponds to the
        the four corners (node 0-3) and the four midpoints (on edges
        0-3) for face face
        Required: 
            face: integer (0-5)
        Returns:
            coordinates: size(8,ndim)
            '''
        if mod(self.Nu,2) == 1:
            midu = [(self.Nu-1)/2,(self.Nu-1)/2]
        else:
            midu = [self.Nu/2,self.Nu/2-1]
        # end if

        if mod(self.Nv,2) == 1:
            midv = [(self.Nv-1)/2,(self.Nv-1)/2]
        else:
            midv = [self.Nv/2,self.Nv/2-1]
        # end if

        if mod(self.Nw,2) == 1:
            midw = [(self.Nw-1)/2,(self.Nw-1)/2]
        else:
            midw = [self.Nw/2,self.Nw/2-1]
        # end if

        if   face == 0:
            values = [self.X[0,0,0],self.X[-1,0,0],self.X[0,-1,0],self.X[-1,-1,0],
                      0.5*(self.X[midu[0],0,0 ] + self.X[midu[1],0,0]),
                      0.5*(self.X[midu[0],-1,0] + self.X[midu[1],-1,0]),
                      0.5*(self.X[0,midv[0],0 ] + self.X[0,midv[1],0]),
                      0.5*(self.X[-1,midv[0],0] + self.X[-1,midv[1],0])]
        elif face == 1:
            values = [self.X[0,0,-1],self.X[-1,0,-1],self.X[0,-1,-1],self.X[-1,-1,-1],
                      0.5*(self.X[midu[0],0,-1 ] + self.X[midu[1],0,-1]),
                      0.5*(self.X[midu[0],-1,-1] + self.X[midu[1],-1,-1]),
                      0.5*(self.X[0,midv[0],-1 ] + self.X[0,midv[1],-1]),
                      0.5*(self.X[-1,midv[0],-1] + self.X[-1,midv[1],-1])]
        elif face == 2:
            values = [self.X[0,0,0],self.X[0,-1,0],self.X[0,0,-1],self.X[0,-1,-1],
                      0.5*(self.X[0,midv[0], 0] + self.X[0,midv[1],0 ]),
                      0.5*(self.X[0,midv[0],-1] + self.X[0,midv[1],-1]),
                      0.5*(self.X[0,0,midw[0] ] + self.X[0,0,midw[1] ]),
                      0.5*(self.X[0,-1,midw[0]] + self.X[0,-1,midw[1]])]
        elif face == 3:
            values = [self.X[-1,0,0],self.X[-1,-1,0],self.X[-1,0,-1],self.X[-1,-1,-1],
                      0.5*(self.X[-1,midv[0], 0] + self.X[-1,midv[1],0 ]),
                      0.5*(self.X[-1,midv[0],-1] + self.X[-1,midv[1],-1]),
                      0.5*(self.X[-1,0,midw[0] ] + self.X[-1,0,midw[1] ]),
                      0.5*(self.X[-1,-1,midw[0]] + self.X[-1,-1,midw[1]])]
        elif face == 4:
            values = [self.X[0,0,0],self.X[-1,0,0],self.X[0,0,-1],self.X[-1,0,-1],
                      0.5*(self.X[midu[0],0,0 ] + self.X[midu[1],0,0 ]),
                      0.5*(self.X[midu[0],0,-1] + self.X[midu[1],0,-1]),
                      0.5*(self.X[0,0,midw[0] ] + self.X[0,0,midw[1] ]),
                      0.5*(self.X[-1,0,midw[0]] + self.X[-1,0,midw[1]])]
        elif face == 5:
            values = [self.X[0,-1,0],self.X[-1,-1,0],self.X[0,-1,-1],self.X[-1,-1,-1],
                      0.5*(self.X[midu[0],-1,0 ] + self.X[midu[1],-1,0 ]),
                      0.5*(self.X[midu[0],-1,-1] + self.X[midu[1],-1,-1]),
                      0.5*(self.X[0,-1,midw[0] ] + self.X[0,-1,midw[1] ]),
                      0.5*(self.X[-1,-1,midw[0]] + self.X[-1,-1,midw[1]])]
        # end if
        return array(values)


    def getMidPointFace(self,face):
        '''Get the midpoint of the face
        on edge.
        Required:
           edge: face index = 0,1,2,3,4,5
        Returns:
            midpoint 
            '''
        assert face in [0,1,2,3,4,5] and self.orig_data == True,'Error, getMidPointFace: No \
original data for this surface or face is not in range 0->5'

        if mod(self.Nu,2) == 1:
            midu = [(self.Nu-1)/2,(self.Nu-1)/2]
        else:
            midu = [self.Nu/2,self.Nu/2-1]
        # end if

        if mod(self.Nv,2) == 1:
            midv = [(self.Nv-1)/2,(self.Nv-1)/2]
        else:
            midv = [self.Nv/2,self.Nv/2-1]
        # end if

        if mod(self.Nw,2) == 1:
            midw = [(self.Nw-1)/2,(self.Nw-1)/2]
        else:
            midw = [self.Nw/2,self.Nw/2-1]
        # end if

        if   face == 0:
            val = (self.X[midu[0],midv[0],0] + self.X[midu[1],midv[0],0] +
                   self.X[midu[0],midv[1],0] + self.X[midu[1],midv[1],0])/4.0
        elif face == 1:
            val = (self.X[midu[0],midv[0],-1] + self.X[midu[1],midv[0],-1] +
                   self.X[midu[0],midv[1],-1] + self.X[midu[1],midv[1],-1])/4.0
        elif face == 2:
            val = (self.X[0,midv[0],midw[0]] + self.X[0,midv[1],midw[0]] +
                   self.X[0,midv[0],midw[1]] + self.X[0,midv[1],midw[1]])/4.0
        elif face == 3:
            val = (self.X[-1,midv[0],midw[0]] + self.X[-1,midv[1],midw[0]] +
                   self.X[-1,midv[0],midw[1]] + self.X[-1,midv[1],midw[1]])/4.0
        elif face == 4:
            val = (self.X[midu[0],0,midw[0]] + self.X[midu[1],0,midw[0]] +
                   self.X[midu[0],0,midw[1]] + self.X[midu[1],0,midw[1]])/4.0
        elif face == 5:
            val = (self.X[midu[0],-1,midw[0]] + self.X[midu[1],-1,midw[0]] +
                   self.X[midu[0],-1,midw[1]] + self.X[midu[1],-1,midw[1]])/4.0
        # end if
        return val
 
    def _setFaceSurfaces(self):
        '''Create face spline objects for each of the faces'''
        self.face_surfaces = [None,None,None,None,None,None]
        self.face_surfaces[0] = surface('create',ku=self.ku,kv=self.kv,tu=self.tu,
                                      tv=self.tv,coef=self.coef[:,:,0,:])
        self.face_surfaces[1] = surface('create',ku=self.ku,kv=self.kv,tu=self.tu,
                                      tv=self.tv,coef=self.coef[:,:,-1,:])
        self.face_surfaces[2] = surface('create',ku=self.ku,kv=self.kw,tu=self.tu,
                                      tv=self.tw,coef=self.coef[:,0,:,:])
        self.face_surfaces[3] = surface('create',ku=self.ku,kv=self.kw,tu=self.tu,
                                      tv=self.tw,coef=self.coef[:,-1,:,:])
        self.face_surfaces[4] = surface('create',ku=self.kv,kv=self.kw,tu=self.tv,
                                      tv=self.tw,coef=self.coef[0,:,:,:])
        self.face_surfaces[5] = surface('create',ku=self.kv,kv=self.kw,tu=self.tv,
                                      tv=self.tw,coef=self.coef[-1,:,:,:])
        
        return

    def _setEdgeCurves(self):
        '''Create edge spline objects for each edge'''
        self.edge_curves = [None,None,None,None,None,None,None,None,None,None,None,None]
        self.edge_curves[0] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,0,0,:])
        self.edge_curves[1] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,-1,0,:])
        self.edge_curves[2] = curve('create',k=self.kv,t=self.tv,coef=self.coef[0,:,0,:])
        self.edge_curves[3] = curve('create',k=self.kv,t=self.tv,coef=self.coef[-1,:,0,:])
        self.edge_curves[4] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,0,-1,:])
        self.edge_curves[5] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,-1,-1,:])
        self.edge_curves[6] = curve('create',k=self.kv,t=self.tv,coef=self.coef[0,:,-1,:])
        self.edge_curves[7] = curve('create',k=self.kv,t=self.tv,coef=self.coef[-1,:,-1,:])
        self.edge_curves[8] = curve('create',k=self.kw,t=self.tw,coef=self.coef[0,0,:,:])
        self.edge_curves[9] = curve('create',k=self.kw,t=self.tw,coef=self.coef[-1,0,:,:])
        self.edge_curves[10] = curve('create',k=self.kw,t=self.tw,coef=self.coef[0,-1,:,:])
        self.edge_curves[11] = curve('create',k=self.kw,t=self.tw,coef=self.coef[-1,-1,:,:])
        return 

    def _getBasisPt(self,u,v,w,vals,istart,col_ind,l_index):
        # This function should only be called from pyBlock The purpose
        # is to compute the basis function for a u,v,w point and add
        # it to pyBlcoks's global dPt/dCoef
        # matrix. vals,row_ptr,col_ind is the CSR data and l_index in
        # the local -> global mapping for this volume

        return pyspline.getbasisptvolume(u,v,w,self.tu,self.tv,self.tw,
                                         self.ku,self.kv,self.kw,vals,
                                         col_ind,istart,l_index)
    
    def __call__(self,u,v,w):
        '''
        Equivalant to getValue
        '''
        return self.getValue(u,v,w)

    def getValue(self,u,v,w):
        '''Get the value at the volume points(s) u,v,w
        Required Arguments:
            u,v,w: u,w and w can be a scalar,vector or matrix of values
        Returns:
           values: An array of size (shape(u), nDim)
           '''
        u = array(u)
        v = array(v)
        w = array(w)

        assert u.shape == v.shape == w.shape,'Error, getValue: u and v must have the same shape'
        if len(u.shape) == 0:
            return pyspline.eval_volume(u,v,w,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,self.coef)
        elif len(u.shape) == 1:
            return pyspline.eval_volume_v(u,v,w,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,self.coef)
        elif len(u.shape) == 2:
            return pyspline.eval_volume_m(u,v,w,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,self.coef)
        elif len(u.shape) == 3:
            return pyspline.eval_volume_t(u,v,w,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,self.coef)
        # end if

    def getBounds(self):
        '''Determine the extents of the volume
        Required: 
            None:
        Returns:
            xmin,xmax: xmin is the lowest x,y,z point and xmax the highest
            '''
        assert self.nDim == 3,'getBounds is only defined for nDim = 3'
        cx = self.coef[:,:,:,0].flatten()
        cy = self.coef[:,:,:,1].flatten()
        cz = self.coef[:,:,:,2].flatten()

        Xmin = zeros(self.nDim)
        Xmin[0] = min(cx)
        Xmin[1] = min(cy)
        Xmin[2] = min(cz)

        Xmax = zeros(self.nDim)
        Xmax[0] = max(cx)
        Xmax[1] = max(cy)
        Xmax[2] = max(cz)

        return Xmin,Xmax

    def projectPoint(self,x0,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''
        Project a point x0 onto the volume and return parametric position

        Required Arguments:
            x0   : A point or points in nDim space for projection
        Optional Arguments:
            Niter: Maximum number of newton iterations
            eps1 : Intersection/Relative change tolerance
            eps2 : Cosine convergence measure tolerance
            u    : Initial guess for u position
            v    : Initial guess for v position
            w    : Initial guess for w position 
        Returns:
            u    : Parametric position(s) u on surface
            v    : Parametric position(s) v on surface
            w    : Parametric position(s) w on surface
            D    : Distance(s) between point(s) and surface(u,v,w)
        '''
        # We will use a starting point u0,v0 if given
        x0 = array(x0)

        if len(x0.shape) == 1: # Just 1 value passed in 
            if 'u' in kwargs:
                u = array([kwargs['u']])
            else:
                u = array([.5])
            # end if
            if 'v' in kwargs:
                v = array([kwargs['v']])
            else:
                v = array([.5])
            # end if
            if 'w' in kwargs:
                w = array([kwargs['w']])
            else:
                w = array([.5])
            # end if
            result = pyspline.point_volume(array([x0]),self.tu,self.tv,self.tw,
                                           self.ku,self.kv,self.kw,
                                           self.coef,Niter,eps1,eps2,u,v,w)
            return result[0][0],result[1][0],result[2][0],result[3][0] #u,v,D
        else:
            if 'u' in kwargs:
                u = array(kwargs['u'])
            else:
                u = 0.5*ones(len(x0))
            # end if
            if 'v' in kwargs:
                v = array(kwargs['v'])
            else:
                v = 0.5*ones(len(x0))
            # end if
            if 'w' in kwargs:
                w = array(kwargs['w'])
            else:
                w = 0.5*ones(len(x0))
            # end if
            result = pyspline.point_volume(\
                x0,self.tu,self.tv,self.tw,self.ku,self.kv,self.kw,\
                    self.coef,Niter,eps1,eps2,u,v,w)
            return result[0],result[1],result[2],result[3]
        # end if


    def _writeTecplotOrigData(self,handle):
        if self.orig_data:
            _writeTecplot3D(handle,'orig_data',self.X)
        # end if

        return
    def _writeTecplotVolume(self,handle):
        '''Output this volume\'s data to a open file handle \'handle\' '''
        
        if self.U ==None:
            # This works really well actually
            Nx = self.Nctlu*self.ku+1
            Ny = self.Nctlv*self.kv+1
            Nz = self.Nctlw*self.kw+1
            u_plot = 0.5*(1-cos(linspace(0,pi,Nx)))
            v_plot = 0.5*(1-cos(linspace(0,pi,Ny)))
            w_plot = 0.5*(1-cos(linspace(0,pi,Nz)))
        
            # Dump re-interpolated surface
            W_plot = zeros((Nx,Ny,Nz))
            V_plot = zeros((Nx,Ny,Nz))
            U_plot = zeros((Nx,Ny,Nz))
            for i in xrange(Nx):
                [V_plot[i,:,:],U_plot[i,:,:]] = meshgrid(v_plot,u_plot)
                W_plot[i,:,:] = w_plot[i]
            # end for
            values = self.getValue(U_plot,V_plot,W_plot)
        else:
            values = self.getValue(self.U,self.V,self.W)

        _writeTecplot3D(handle,'interpolated',values)
        
        return

    def _writeTecplotCoef(self,handle):
        '''Write the Spline coefficients to handle'''
        _writeTecplot3D(handle,'control_pts',self.coef)

        return

    def writeTecplot(self,file_name,vols=True,coef=True,orig=False,tecio=USE_TECIO):
        '''Write the surface to a tecplot dat file
        Required Arguments:
            file_name: The output file name
        Optional Arguments:
            surfs: Boolean to write interpolated surfaces (default=True)
            coef : Boolean to write coefficients (default=True)
            orig : Boolean to write original data (default=True)
            dir  : Boolean to write out surface direction indicators (default=False)

            '''
        f = openTecplot(file_name,self.nDim,tecio)
        if vols:
            self._writeTecplotVolume(f)
        if coef:
            self._writeTecplotCoef(f)
        if orig:
            self._writeTecplotOrigData(f)
        closeTecplot(f)

    def _writeBvol(self,handle,binary):
        '''Write the data to an open file handle'''
        # Initial Integers
        init_values = array([self.Nctlu,self.Nctlv,self.Nctlw,
                             self.ku,self.kv,self.kw])
        if binary:
            init_values.tofile(handle,sep="")
            # Knot Vectors
            self.tu.tofile(handle,sep="")
            self.tv.tofile(handle,sep="")
            self.tw.tofile(handle,sep="")
            # Control Points
            self.coef.flatten().tofile(handle)
        else:
            init_values.tofile(handle,sep="\n",format="%d")
            handle.write('\n')
            # Knot Vectors
            self.tu.tofile(handle,sep="\n",format="%f")
            handle.write('\n')
            self.tv.tofile(handle,sep="\n",format="%f")
            handle.write('\n')
            self.tw.tofile(handle,sep="\n",format="%f")
            handle.write('\n')
            # Control Points
            self.coef.flatten().tofile(handle,sep="\n",format="%f")
            handle.write('\n')
        # end if
        
        return 
        
# ----------------------------------------------------------------------
#                     Misc Helper Functions
# ----------------------------------------------------------------------

def trilinear_volume(*args,**kwargs):
    '''This is a short-cut function to create a trilinear b-spline volume
    Args can contain:
        X: array of size(8,3) which contains the corners of the box is
        coordinate order (i,j,k)
      
        xmin,xmax: The lower extreme and upper extreme corners of the box
        
        '''
    tu = [0,0,1,1]
    tv = [0,0,1,1]
    tw = [0,0,1,1]
    ku = 2
    kv = 2
    kw = 2
    if len(args) == 0:
        return volume(coef=args[0],tu=tu,tv=tv,tw=tw,ku=ku,kv=kv,kw=kw)
    elif len(args) == 2:
        xmin = args[0]
        xmax = args[1]

        x_low  = xmin[0]
        x_high = xmax[0]
        y_low  = xmin[1]
        y_high = xmax[1]
        z_low  = xmin[2]
        z_high = xmax[2]

        coef = zeros((2,2,2,3))
        coef[0,0,0,:] = [x_low,y_low,z_low]
        coef[1,0,0,:] = [x_high,y_low,z_low]
        coef[0,1,0,:] = [x_low,y_high,z_low]
        coef[1,1,0,:] = [x_high,y_high,z_low]
        coef[0,0,1,:] = [x_low,y_low,z_high]
        coef[1,0,1,:] = [x_high,y_low,z_high]
        coef[0,1,1,:] = [x_low,y_high,z_high]
        coef[1,1,1,:] = [x_high,y_high,z_high]
        return volume(coef=coef,tu=tu,tv=tv,tw=tw,ku=ku,kv=kv,kw=kw)
    else:
        mpiPrint('Error: An unknown number of arguments was passed to trilinear\
 volume')
        sys.exit(1)
    # end if


def bilinear_surface(*args,**kwargs):
    '''This is short-cut function to create a bilinear surface
    Args can contain:
        x: array of size(4,3) The four corners of the array arranged in
        the coordinate direction orientation:

        2          3
        /----------\
        |          |
        |          |
        |          |
        \----------/
        0          1
    
   OR

   Args can contain pt1,pt2,pt3,pt4 is CCW Ordering

        3          2
        /----------\
        |          |
        |          |
        |          |
        \----------/
        0          1
        '''
    if len(args) == 1:
        # One argument passed in ... assume its X
        assert len(args[0]) == 4,'Error: a single argument passed to bilinear\
 surface must contain 4 points and be of size (4,3)'
        coef = zeros((2,2,3))
        coef[0,0] = args[0][0]
        coef[1,0] = args[0][1]
        coef[0,1] = args[0][2]
        coef[1,1] = args[0][3]
        return surface(coef=coef,tu=[0,0,1,1],tv=[0,0,1,1],ku=2,kv=2)
    else:
        # Assume 4 arguments
        coef = zeros([2,2,3])
        coef[0,0] = args[0]
        coef[0,1] = args[1]
        coef[1,0] = args[3]
        coef[1,1] = args[2]
        return surface(coef=coef,tu=[0,0,1,1],tv=[0,0,1,1],ku=2,kv=2)
    # end if

def line(*args,**kwargs):
    '''This is a short cut function to create a line curve

    Args can contain:
       X: array of size(2,ndim) The two end points
       
       OR:
       x1,x2: The two end points (each of size ndim)

       OR:
       x1,dir=direction
       x1 and the keyword argument direction

       OR: 
       x1,dir=direction,length=length
       x1, direction and specific length
       '''
    if len(args) == 2:
        # Its a two-point type
        return curve(coef=[args[0],args[1]],k=2,t=[0,0,1,1])
    elif len(args) == 1:
        if len(args[0]) == 2: # its X
            return curve(coef=args[0],k=2,t=[0,0,1,1])
        elif 'dir' in kwargs:
            # We have point and direction
            if 'length' in kwargs:
                x2 = args[0] + kwargs['dir']/norm(kwargs['dir'])*kwargs['length']
            else:
                x2 = args[0] + kwargs['dir']
            # end if
            return curve(coef=[args[0],x2],k=2,t=[0,0,1,1])
        else:
            mpiPrint('Error: dir must be specified if only 1 argument is given')
        # end if
    # end if


#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
    print 'There are two examples in the example directory.'
    print 'Look at test_curve.py and test_surf.py for more informatin'
