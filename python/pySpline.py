'''
pySpline

Contains an relatively thin interface to the cmlib spline functions

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

__version__ = '$Revision: $'

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, time, copy, pdb

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,zeros,sqrt,array,reshape,meshgrid,mod,floor,\
    ones,vstack,real,where,arange,append,hstack

import scipy
from scipy import sparse
from scipy.sparse.linalg.dsolve import factorized
scipy.sparse.linalg.use_solver(useUmfpack=False)

from numpy.linalg import lstsq
import pyspline

from mdo_import_helper import *

# =============================================================================
def e_dist(x1,x2):
    '''Get the eculidean distance between two points in nDim space'''
    total = 0.0
    for idim in xrange(len(x1)):
        total += (x1[idim]-x2[idim])**2
    # end for
    return sqrt(total)

# =============================================================================
# pySpline class
# =============================================================================

class surface():

    def __init__(self,*args,**kwargs):

        '''Create an instance of a b-spline surface. There are three ways to
        initialize the class as determined by the first argument:

        args[0] = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for u
            kv, integer: Order for v
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control point 
                                                     values

        args[0] = \'interpolate\': Create an instance of the spline class
        by using an interpolating spline to given data points. **kwarg
        MUST contain the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: list of u values (Optional for ndim=3,required otherwise)
            v, real, array: list of v values (Optional for ndim=3,required otherwise)
            X, real, array, size(len(u),len(v),nDim):Array of data points to fit
            --> OR  x=array, size(len(u),len(v)) -> two dimensional interpolation
            --> OR  x=array,y=array, each of size(len(u),len(v))
            --> OR  x=array,y=array,z=array, each of size (len(u),len(v))

        args[0] \'lms\': Create an instance of the spline class using a
        Least-Mean-Squares fit to the spline. . **kwargs MUST contain
        the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            Nctlu, integer: Number of control points for u
            Nctlv, integer: Number of control points for v
            u, real, array: list of u values (Optional for ndim=3,required otherwise)
            v, real, array: list of v values (Optional for ndim=3,required otherwise)
            X, real, array, size(len(u),len(v),nDim):Array of data points to fit
            --> OR  x=array, size(len(u),len(v)) -> two dimensional interpolation
            --> OR  x=array,y=array, each of size(len(u),len(v))
            --> OR  x=array,y=array,z=array, each of size (len(u),len(v))
'''
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      
        task = args[0]
        mpiPrint('pySpline Type: %s. '%(task),self.NO_PRINT)
        
        # Defaults
        self.X    = None
        self.u    = None
        self.v    = None
        self.Nu   = None
        self.Nv   = None
        self.orig_data = False

        if task == 'create':
            assert 'ku' in kwargs and 'kv' in kwargs  and 'tu' in kwargs \
                and 'tv' in kwargs and 'coef' in kwargs,'Error: ku,kv,tu,tv,\
 and coef MUST be defined for task=create'
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
            self._setEdgeCurves()
            return
     
        elif task == 'interpolate' or task == 'lms':
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and \
                ( 'X' in kwargs or 'x' in kwargs or \
                      ( 'x' in kwargs and 'y' in kwargs) or \
                      ( 'x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
                      'Error: ku,kv, and X (or x or x and y or x and y and z \
MUST be defined for task lms or interpolate'
            if task=='lms':
                assert 'Nctlu' in kwargs and 'Nctlv' in kwargs,'Nctlu and Nctlv\
 MUST be defined for task lms'
            
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
            if task == 'interpolate':
                self.Nctlu = self.X.shape[0]
                self.Nctlv = self.X.shape[1]
            else:
                self.Nctlu = kwargs['Nctlu']
                self.Nctlv = kwargs['Nctlv']
            # end if
                
            self.orig_data = True
            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]

            self.ku = int(kwargs['ku'])
            self.kv = int(kwargs['kv'])

            # Sanity Check on Inputs
            if self.Nctlu >= self.Nu:  self.Nctlu  = self.Nu
            if self.Nctlv >= self.Nv:  self.Nctlv = self.Nv

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku: self.ku = self.Nu
            if self.Nv < self.kv: self.kv = self.Nv
            if self.Nctlu < self.ku: self.ku = self.Nctlu
            if self.Nctlv < self.kv: self.kv = self.Nctlv

            if 'rel_tol' in kwargs:
                self.rel_tol = kwargs['rel_tol']
            else:
                self.rel_tol = 5e-4
            # end if

            if 'niter' in kwargs:
                self.niter = kwargs['niter']
            else:
                self.niter = 1000
            # end if

            if 'u' in kwargs and 'v' in kwargs:
                u = kwargs['u']
                v = kwargs['v']
            else:
                if self.nDim == 3:
                    u,v = self._calcParameterization()
                else:
                    mpiPrint('Automatic parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise.')
                    sys.exit(1)
                # end if
            # end if
            self.umin = u[0]
            self.umax = u[-1]
            self.vmin = v[0]
            self.vmax = v[-1]
            
            self.u = u
            self.v = v
            
            self._calcKnots()
            [self.V,self.U] = meshgrid(v,u)
#             self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim),'d')
#             self.niter = 1      
#             self.U,self.V,self.coef, rms = pyspline.compute_surface(\
#                 self.X,self.U,self.V,self.tu,self.tv,self.ku,self.kv,self.coef,\
#                     self.niter,self.rel_tol)
  
            self.recompute()
            #mpiPrint('Done Fitting Surface')

        else:
            mpiPrint('Error: The first argument must be \'create\', \'lms\' or \'interpolate\'')
            sys.exit(0)
        # end if (init type)
        return

    def _calcKnots(self):
        self.tu = pyspline.knots(self.u,self.Nctlu,self.ku)
        self.tv = pyspline.knots(self.v,self.Nctlv,self.kv)


    def recompute(self):
        '''Recompute the curve if the knot vector has been modified'''
        # --- Direct Implementation
        Jac = pyspline.surface_jacobian_wrap2(self.U,self.V,self.tu,self.tv,
                                             self.ku,self.kv,self.Nctlu,
                                             self.Nctlv)
        self.coef = lstsq(Jac,self.X.reshape([self.Nu*self.Nv,self.nDim]))[0].reshape([self.Nctlu,self.Nctlv,self.nDim])

      #   # ---- Sparse Implemention

#         vals,row_ptr,col_ind = pyspline.surface_jacobian_wrap(\
#             self.U,self.V,self.tu,self.tv,self.ku,self.kv,self.Nctlu,self.Nctlv)
#         Jac2 = sparse.csr_matrix((vals,col_ind,row_ptr),
#                                  shape=[self.Nu*self.Nv,self.Nctlu*self.Nctlv])
#         self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim))
#         for idim in xrange(self.nDim):
#             J = Jac2.transpose()*Jac2
#             rhs  = Jac2.transpose()*self.X[:,:,idim].flatten()
#             self.coef[:,:,idim] =  linsolve.spsolve(J,rhs).reshape([self.Nctlu,self.Nctlv])

        self._setEdgeCurves()
        return

    def _calcParameterization(self):
        u = zeros(self.Nu,'d')
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
            else:
                temp /= temp[-1]
            # end if

            u += temp #accumulate the u-parameter calcs for each j
        # end for 
        u =u/(self.Nv-singular_counter) #divide by the number of 'j's we had
        
        v = zeros(self.Nv,'d')
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
            else:
                temp /= temp[-1]
            #end if 

            v += temp #accumulate the v-parameter calcs for each i
        # end for 
        v = v/(self.Nu-singular_counter) #divide by the number of 'i's we had

        return u,v

    def _setEdgeCurves(self):
        '''Create curve spline objects for each of the edges'''
        self.edge_curves = [None,None,None,None]
        self.edge_curves[0] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,0])
        self.edge_curves[1] = curve('create',k=self.ku,t=self.tu,coef=self.coef[:,-1])
        self.edge_curves[2] = curve('create',k=self.kv,t=self.tv,coef=self.coef[0,:])
        self.edge_curves[3] = curve('create',k=self.kv,t=self.tv,coef=self.coef[-1,:])

        return

    def getValueCorner(self,corner):
        '''Get the value of the spline on corner i,i=0,1,2 or 3 '''
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
        ''' Get the two end points and midpoint values of the original data on edge. edge = 0,1,2,3'''
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
        ''' Get the values of the original data on cornere. conrner = 0,1,2,3'''
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
        return pyspline.getbasispt(u,v,self.tu,self.tv,self.ku,self.kv,vals,
                                   col_ind,istart,l_index)
                                
    def __call__(self,u,v):
        '''Equivalant to getValue'''
        return self.getValue(u,v)

    def getValue(self,u,v):
        '''Get the value of the spline at point(s) u,v
        u and v can be scalars,vectors or arrays. They must be the same size'''
        #Assure u,v are numpy arrays
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
        '''Get the value of the derivative of spline at point(s) u,v
        u and v can be scalars,vectors or arrays. They must be the same size'''
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
        u and v can be scalars,vectors or arrays. They must be the same size'''
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
        '''Get the value of the second derivative (d2u2,d2v2,d2uv) of spline at point(s) u,v
        u and v can be scalars,vectors or arrays. They must be the same size'''
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
 
    def projectPoint(self,x0,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''Project a point x0 onto the surface. i.e. Find the point on the
        surface that minimizes the distance from x0 to
        surface(u,v).
        Starting points can be given by
        u=<val>
        v=<val>
        '''
        # We will use a starting point u0,v0 if given
        assert len(x0) == self.nDim,'Error: x0 must have same spatial dimension as surface'
        u=-1
        v=-1
        if 'u' in kwargs: u = kwargs['u']
        if 'v' in kwargs: v = kwargs['v']
        return pyspline.point_surface(x0,self.tu,self.tv,self.ku,self.kv,self.coef,
                                              Niter,eps1,eps2,u,v)
      
    def projectCurve(self,curve,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''Project a curve object curve onto this curve. This function
        finds the parametric values of the shortest distances between
        curves. If the curves intersection it will return ONE of the
        intersections
        Starting points can be given with 
        u=<val>(surface)
        v=<val>(surface)
        s=<val> (curve)'''
        # We will use a starting point u0,v0 if given
        u = -1.0
        v = -1.0
        s = -1.0
        if 'u' in kwargs: u = kwargs['u']
        if 'v' in kwargs: v = kwargs['v']
        if 's' in kwargs: s = kwargs['s']
        return pyspline.curve_surface(\
            curve.t,curve.k,curve.coef,self.tu,self.tv,\
                self.ku,self.kv,self.coef,Niter,eps1,eps2,u,v,s)
   
    def _writeTecplot2D(self,handle,name,data):
         '''A Generic write tecplot zone to file'''
         nx = data.shape[0]
         ny = data.shape[1]
         handle.write('Zone T=%s I=%d J=%d\n'%(name,nx,ny))
         handle.write('DATAPACKING=POINT\n')
         for j in xrange(ny):
             for i in xrange(nx):
                 handle.write('%f %f %f \n'%(data[i,j,0],data[i,j,1],data[i,j,2]))
             # end for
         # end for 
         
         return

    def _writeTecplotOrigData(self,handle):
        if self.orig_data:
            self._writeTecplot2D(handle,'orig_data',self.X)
        # end if

        return

    def _writeTecplotSurface(self,handle,size=None):
        '''Output this surface\'s data to a open file handle \'handle\' '''
        
        # The size should be based on the knot vectors
               
        MAX_SIZE = 100
        MIN_SIZE = 5
                
        if size == None:
            # This works really well actually
            u_plot = 0.5*(1-cos(linspace(0,pi,25)))
            v_plot = 0.5*(1-cos(linspace(0,pi,25)))
        else:
            # Cheaply calculate the length of each size of the surface
            # to determine its lenght and then use the size parameter 
            # to determine the number of points to use

            u = array([0,0.5,1,1,1,0.5,0,0])
            v = array([0,0,0,0.5,1,1,1,0.5])
            val = self.getValue(u,v)

            u_len = (e_dist(val[0],val[1])+e_dist(val[1],val[2])+
                     e_dist(val[4],val[5])+e_dist(val[5],val[6]))/2
            
            v_len = (e_dist(val[2],val[3])+e_dist(val[3],val[4])+
                     e_dist(val[6],val[7])+e_dist(val[7],val[0]))/2
            
            nu=int(floor(real(u_len/size)))
            nv=int(floor(real(v_len/size)))
            
            if nu > MAX_SIZE: nu = MAX_SIZE
            if nu < MIN_SIZE: nu = MIN_SIZE
            if nv > MAX_SIZE: nv = MAX_SIZE
            if nv < MIN_SIZE: nv = MIN_SIZE
            
            u_plot = linspace(self.umin,self.umax,nu)
            v_plot = linspace(self.vmin,self.vmax,nv)
        # end if

        # Dump re-interpolated surface
        [V_plot,U_plot] = meshgrid(v_plot,u_plot)
        X = self.getValue(U_plot,V_plot)
        self._writeTecplot2D(handle,'interpolated',X)
    
        return
                    
    def _writeTecplotCoef(self,handle):
        '''Write the Spline coefficients to handle'''
        self._writeTecplot2D(handle,'control_pts',self.coef)

        return

    def _writeTecplotEdge(self,handle,edge,*args,**kwargs):
        '''Dump out a linear zone along edge used for visualizing edge connections'''
        self.edge_curves[edge].writeTecplotEdge(handle,*args,**kwargs)

        return

    def _writeDirections(self,handle,isurf):
        '''Write out and indication of the surface direction'''
        handle.write('Zone T=\"surface%d direction\" I=4\n'%(isurf))
        if self.Nctlu >= 3 and self.Nctlv >=3:
            handle.write('%f,%f,%f \n'%(self.coef[1,2,0],self.coef[1,2,1],self.coef[1,2,2]))
            handle.write('%f,%f,%f \n'%(self.coef[1,1,0],self.coef[1,1,1],self.coef[1,1,2]))
            handle.write('%f,%f,%f \n'%(self.coef[2,1,0],self.coef[2,1,1],self.coef[2,1,2]))
            handle.write('%f,%f,%f \n'%(self.coef[3,1,0],self.coef[3,1,1],self.coef[3,1,2]))
        else:
            mpiPrint('Not Enough control points to output direction indicator')
        #end if

        return 

    def writeTecplot(self,file_name,surfs=True,coef=True,orig=True,dir=False):
        '''Write the surface to tecplot'''
        f = open(file_name,'w')
        if surfs:
            self._writeTecplotSurface(f)
        if coef:
            self._writeTecplotCoef(f)
        if orig:
            self._writeTecplotOrigData(f)
        if dir:
            self._writeDirections(f,0)
        f.close()

    def _writeIGES_directory(self,handle,Dcount,Pcount):
        '''Write the IGES file header information (Directory Entry Section)
        for this surface'''
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
        '''Write the IGES parameter information for this surface'''

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

class curve():

    def __init__(self,*args,**kwargs):
        '''
        Create an instance of a b-spline curve. There are three
        ways to initialize the class as determined by the first argument

        args[0] = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            k, integer: Order for spline
            Nctu, integer: Number of u control points
            t, real array: Knot vector 
            coef, real array size(Nctl,nDim+1): Array of control points and weights

        arg[0] = \'interpolate\': Create an instance of the spline class
        by using an interpolating spline to given data points. **kwarg
        MUST contain the following information:

            k, integer: Order for spline
            X, real, array, size(len(s),nDim): Array of data to fit 
              --> OR  x=<vals>, len(s) ->1D interpolation
              --> OR  x=<vals>,y=<vals>, each of len(s) -> 2D parametric curve
              --> OR  x=<vals>,y=<vals>,z=<vals> each of len(s) ->3D parametric curve
            s, real, array: (OPTIONAL for nDim >=2 ) Array of s values 

        arg[0] = \'lms\': Create an instance of the spline class using a lms spline.
       **kwargs MUST contain the following information:
            k, integer: Order for spline
            s, real, array: (OPTIONAL for nDIM >=2) Array of s values 
            X, real, array, size(len(s),nDim): Array of data to fit OR
                 x=<vals>,y=<vals>,z=<vals> where each is of length s
            Nctl, integer: The number of control points
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
                self.X = array(kwargs['x'])
                self.nDim = 1
            # enf if
            self.k = int(kwargs['k'])
            self.N = len(self.X)
            if 'niter' in kwargs: 
                self.niter = kwargs['niter']
            else:
                self.niter = 10
            # end if

            if 'rel_tol' in kwargs:
                self.rel_tol = kwargs['rel_tol']
            else:
                self.rel_tol = 5e-4
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
                weights = array(kwargs['weights'])
            else:
                weights = ones(self.N)
            # end if

            if 'deriv' in kwargs and 'deriv_ptr' in kwargs:
                deriv = array(kwargs['deriv'])
                deriv_ptr = array(kwargs['deriv_ptr'])
            else:
                deriv = None
                deriv_ptr = array([])
            # end if
            if 'deriv_weights' in kwargs and deriv != None:
                deriv_weights = array(kwargs['deriv_weights'])
            else:
                if deriv != None:
                    deriv_weights = ones(len(deriv))
                else:
                    deriv_weights = None
                # end if
            # end if
            if not 'Nctl' in kwargs: # We are doing an interpolation...set all weights to -1
                weights[:] = 1
                if deriv_weights != None:
                    deriv_weights[:] = 1
                # end if
            # end if
            # Now do the sparation between the constrained and unconstrained
            S  = self.X[where(weights > 0.0)]
            su = self.s[where(weights > 0.0)]
            T  = self.X[where(weights <=  0.0)]
            sc = self.s[where(weights <= 0.0)]
            weights = weights[where(weights > 0.0)]
            nu = len(S)
            nc = len(T)
     
            # And the derivative info
            if deriv != None:
                S = vstack((S,deriv[where(deriv_weights > 0.0)]))
                sdu = self.s[deriv_ptr][where(deriv_weights > 0.0)]
                T = vstack((T,deriv[where(deriv_weights <= 0.0)]))
                sdc = self.s[deriv_ptr][where(deriv_weights <= 0.0)] 
                weights = append(weights,deriv_weights[where(deriv_weights > 0.0)])
                ndu = len(sdu)
                ndc = len(sdc)
            else:
                sdu = array([],'d')
                sdc = array([],'d')
                ndu = 0
                ndc = 0
            # end if
            print 'weights:',weights
            W = sparse.csr_matrix((weights,arange(len(weights)),arange(len(weights)+1))) # Diagonal
            
            if 'Nctl' in kwargs:
                self.Nctl = int(kwargs['Nctl'])
                if self.Nctl > nu+nc+ndu+ndc:
                    self.Nctl = nu+nc+ndu+ndc
                    interp = True
            else:
                self.Nctl = nu+nc+ndu+ndc
                self.niter = 1
                interp = True
            # end if
            # Sanity check to make sure k is ok
            if nu+nc+ndu+ndc < self.k:
                self.k = nu+nc+ndu+ndc
            # end if

            # Generate the knot vector,greville points and empty coefficients
            if interp:
                print 'Here'
                self.t = pyspline.knots_interp(self.s,deriv_ptr,self.k)
            else:
                self.t = pyspline.knots_lms(self.s,self.Nctl,self.k)
            # end if
            self._calcGrevillePoints()
            self.coef = zeros((self.Nctl,self.nDim),'d')

            # Get the 'N' jacobain
            N_vals = zeros((nu+ndu)*self.k)          # |
            N_row_ptr = zeros(nu+ndu+1,'intc')       # | ->  Standard CSR formulation
            N_col_ind = zeros((nu+ndu)*self.k,'intc')# |
            print 'su,sdu',su,sdu
            print 'nu,ndu:',nu,ndu
            print 'knots:',self.t
            pyspline.curve_jacobian_wrap(su,sdu,self.t,self.k,self.Nctl,N_vals,N_row_ptr,N_col_ind)
            N = sparse.csr_matrix((N_vals,N_col_ind,N_row_ptr),shape=[nu+ndu,self.Nctl])

            # If there is no constrained data...do a regular fit
            if N.shape[0] == N.shape[1]: # Square N...will ONLY happen for pure interpolation 
                print N.todense()
                solve = factorized( N ) # Factorize once for efficiency
                for idim in xrange(self.nDim):
                    self.coef[:,idim] = solve(S[:,idim])
                # end for
                return

            NTWN = N.transpose()*W*N # We need this either way
            if nc + ndc == 0:
                solve = factorized( NTWN ) # Factorize once for efficiency
                for idim in xrange(self.nDim):
                    self.coef[:,idim] = solve(N.transpose()*W*S[:,idim])
                # end for
            else:  # Now its more complicated since we have constraints
                M_vals = zeros((nc+ndc)*self.k)          #|
                M_row_ptr = zeros(nc+ndc+1,'intc')       #| -> Standard CSR formulation
                M_col_ind = zeros((nc+ndc)*self.k,'intc')#|
                pyspline.curve_jacobian_wrap(sc,sdc,self.t,self.k,self.Nctl,M_vals,M_row_ptr,M_col_ind)
                M = sparse.csr_matrix((M_vals,M_col_ind,M_row_ptr),shape=[nc+ndc,self.Nctl])
                # Now we must assemble the constrained jacobian
                # [ N^T*W*T      M^T][P] = [ N^T*W*S]
                # [ M            0  ][R]   [ T      ] 

                NTWN = NTWN.tocsr() # Convert to CSR for the constr_jac function
                MT   = M.transpose().tocsr() # Convert to CSR for the constr_jac function
                j_val,j_col_ind,j_row_ptr = pyspline.constr_jac(
                    NTWN.data,NTWN.indptr,NTWN.indices,MT.data,MT.indptr,MT.indices,
                    M.data,M.indptr,M.indices,self.Nctl)
                # Create sparse csr matrix and factorize
                J = sparse.csr_matrix((j_val,j_col_ind,j_row_ptr),shape=[self.Nctl+nc+ndc,self.Nctl+nc+ndc])
                solve = factorized( J )
                for idim in xrange(self.nDim):
                    rhs = hstack((N.transpose()*W*S[:,idim],T[:,idim]))
                    self.coef[:,idim] = solve(rhs)[0:self.Nctl]
                # end for
            # end if (constr - not constrained
        # end if
        return

    def runParameterCorrection(self,niter,rel_tol=1e-5):
        '''Run niter of the Hoschek's parameter correction algorithim'''
        # Run more parameter correction
        self.s,self.coef,rms = pyspline.compute_curve(\
            self.s,self.X,self.t,self.k,self.coef,niter,self.rel_tol)
        mpiPrint('ReFitted with rms: %f'%rms,self.NO_PRINT)

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
        '''Insert at knot at u'''
        t_new,coef_new,break_pt = pyspline.insertknot(u,r,self.t,self.k,self.coef)
        self.t = t_new
        self.coef = coef_new
        self.Nctl = self.Nctl + r
        
        return

    def splitCurve(self,u):
        ''' Split the curve at knot position u and return the two new
        curves'''
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
        '''Compute the length of the curve using the Eculdian Norm'''
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
        '''Equivilant to getValue'''
        return self.getValue(s)
    
    def getValue(self,s):
        '''Get The value of the spline
        s can be a scalar or vector of values'''
        s = array(s)
        if len(s.shape) == 0:
            return pyspline.eval_curve(s,self.t,self.k,self.coef)
        elif len(s.shape) == 1:
            return pyspline.eval_curve_v(s,self.t,self.k,self.coef)

    def getDerivative(self,s):
        '''Get the value of the derivatve of the spline
        s can be a scalar or vector of values'''
        s = array(s)
        if len(s.shape) == 0:
            return pyspline.eval_curve_deriv(s,self.t,self.k,self.coef)
        elif len(s.shape) == 1:
            return pyspline.eval_curve_deriv_v(s,self.t,self.k,self.coef)

    def projectPoint(self,x0,Niter=20,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''Project a point x0 onto the curve and return parametric position'''

        # eps1 is tolerance for point on curve and for delta change in parameter
        # eps2 is tolerance for cosine convergence test (point off curve)

        assert len(x0) == self.nDim,'Dimension of x0 and the dimension of\
 spline must be the same'
        s = -1.0
        if 's' in kwargs: s=kwargs['s']
        return pyspline.point_curve(x0,self.t,self.k,self.coef,Niter,eps1,eps2,s)

    def projectCurve(self,curve,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):
        '''Find the minimum distance between this curve (self) and a second
        curve passed in (curve)'''

        # eps1 is tolerance for intersection
        # eps2 is tolerance for cosine convergence test (non-intersecting)

        s = -1.0
        t = -1.0
        if 's' in kwargs:  s = kwargs['s']
        if 't' in kwargs:  t = kwargs['t']
        return pyspline.curve_curve(self.t,self.k,self.coef,
                                    curve.t,curve.k, curve.coef,
                                    Niter,eps1,eps2,s,t)

    def writeTecplot(self,file_name,curve=True,coef=True,orig=True,*args,**kwargs):
        '''Write the cuve to tecplot'''
        f = open(file_name,'w')
        if curve:
            self._writeTecplotCurve(f,*args,**kwargs)
        if coef:
            self._writeTecplotCoef(f,*args,**kwargs)
        if orig:
            self._writeTecplotOrigData(f,*args,**kwargs)
        f.close()

    def _writeTecplot1D(self,handle,name,data):
        '''A Generic write tecplot zone to file'''
        n = data.shape[0]
        handle.write('Zone T=%s I=%d \n'%(name,n))
        handle.write('DATAPACKING=POINT\n')
        for i in xrange(n):
            for idim in xrange(self.nDim):
                handle.write('%f '%data[i,idim])
            # end for
            handle.write('\n')
        # end for
        return

    def _writeTecplotCoef(self,handle,*args,**kwargs):
        '''Write the Spline coefficients to handle'''
        self._writeTecplot1D(handle,'control_pts',self.coef)
        return

    def _writeTecplotOrigData(self,handle,*args,**kwargs):
        '''Write the original data to a handle'''
        if self.orig_data:
            self._writeTecplot1D(handle,'orig_data',self.X)
        # end if

        return
    def _writeTecplotCurve(self,handle,*args,**kwargs):
        '''Write the interpolated curve to a handle'''
        if 'size' in kwargs:
            length = self.getLength()
            n=int(floor(real(length/kwargs['size'])))
            X = self.getValue(linspace(0,1,n))
        else:
            if not self.s == None:
                X = self.getValue(self.s)
            else:
                s = 0.5*(1-cos(linspace(0,pi,25)))
                X = self.getValue(s)
            # end if
        # end if
        self._writeTecplot1D(handle,'interpolated',X)

        return 

#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
    print 'There are two examples in the example directory.'
    print 'Look at test_curve.py and test_surf.py for more informatin'
