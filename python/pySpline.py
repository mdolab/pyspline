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
import os, sys, string, time, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,zeros,sqrt,array,reshape,meshgrid,mod,floor,\
vstack,real

import numpy.linalg
import pyspline

from mdo_import_helper import *

# =============================================================================
def e_dist(x1,x2):
    '''Get the eculidean distance between two points in e3'''
    return sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
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
            return
     
        elif task == 'interpolate' or task == 'lms':
            # Do some checking on the number of control points
            assert 'ku' in kwargs and 'kv' in kwargs and \
                ('X' in kwargs or 'x' in kwargs or \
                     ('x' in kwargs and 'y' in kwargs) or \
                     ('x' in kwargs and 'y' in kwargs and 'z' in kwargs)), \
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
            if self.Nctlu >= self.Nu:
                self.Nctlu  = self.Nu 
            if self.Nctlv >= self.Nv:
                self.Nctlv = self.Nv
            # end if

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku:
                self.ku = self.Nu
            # end if

            if self.Nv < self.kv:
                self.kv = self.Nv
            # end if

            if self.Nctlu < self.ku:
                self.ku = self.Nctlu
            # end if

            if self.Nctlv < self.kv:
                self.kv = self.Nctlv
            # end if

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
            self.tu = pyspline.knots(u,self.Nctlu,self.ku)
            self.tv = pyspline.knots(v,self.Nctlv,self.kv)

            #self._calcKnots()
            [self.V,self.U] = meshgrid(v,u)
            self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim),'d')
            self.niter = 1
            print self.Nctlu,self.Nctlv
            self.U,self.V,self.coef, rms = pyspline.compute_surface(\
                self.X,self.U,self.V,self.tu,self.tv,self.ku,self.kv,self.coef,\
                    self.niter,self.rel_tol)
            print 'RMS is:',rms

            # Now do it dense
            Jac = pyspline.surface_jacobian_linear2(self.U,self.V,self.tu,self.tv,
                                                    self.ku,self.kv,self.Nctlu,
                                                    self.Nctlv)
            coef2 = zeros((self.Nctlu,self.Nctlv,self.nDim))
            for idim in xrange(self.nDim):
                coef2[:,:,idim] = numpy.linalg.lstsq(Jac,self.X[:,:,idim].flatten())[0].reshape((self.Nctlu,self.Nctlv))
            # end for

            print 'x - Comparison:'
            print self.coef[:,:,0] - coef2[:,:,0]
            print 'y - Comparison:'
            print self.coef[:,:,1] - coef2[:,:,2]
            print 'z - Comparison:'
            print self.coef[:,:,1] - coef2[:,:,2]
            self.coef = coef2



        else:
            mpiPrint('Error: The first argument must be \'create\', \'lms\' or \'interpolate\'')
            sys.exit(0)
        # end if (init type)
        return

    def _calcKnots(self):
        self.tu = pyspline.knots(u,self.Nctlu,self.ku)
        self.tv = pyspline.knots(v,self.Nctlv,self.kv)

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
        self.edge_curves = []
        self.edge_curves.append(curve('create',k=self.ku,t=self.tu,coef=self.coef[:,0]))
        self.edge_curves.append(curve('create',k=self.ku,t=self.tu,coef=self.coef[:,-1]))
        self.edge_curves.append(curve('create',k=self.kv,t=self.tv,coef=self.coef[0,:]))
        self.edge_curves.append(curve('create',k=self.kv,t=self.tv,coef=self.coef[-1,:]))

        return

    def getValueCorner(self,corner):
        '''Get the value of the spline on corner i '''
        assert corner in [0,1,2,3],'Error, getValueCorner: Corner must be in range 0->3'
        if corner == 0:
            return self.getValue(self.umin,self.vmin)
        elif corner == 1:
            return self.getValue(self.umax,self.vmin)
        elif corner == 2:
            return self.getValue(self.umin,self.vamx)
        elif corner ==3:
            return self.getValue(self.umax,self.vmax)
        #end if

    def getOrigValuesEdge(self,edge):
        ''' Get the values of the original data on edge. edge = 0,1,2,3'''
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
        ''' Get the values of the original data on edge. edge = 0,1,2,3'''
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
                                
    def __call__(self,u,v):
        return self.getValue(u,v)

    def getValue(self,u,v):
        '''Get the value of the spline at point(s) u,v'''
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
        '''Get the value of the derivative of spline at point(s) u,v'''
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
        '''Get the normalized normal at the surface point(s) u,v'''
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
        '''Get the value of the second derivative (d2u2,d2v2,d2uv) of spline at point(s) u,v'''
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
        surface(u,v). See the fortran function project_point_surface.90'''

        # We will use a starting point u0,v0 if given
        assert len(x0) == self.nDim,'Error: x0 must have same spatial dimension as surface'
        u=-1
        v=-1
        if 'u' in kwargs: u = kwargs['u']
        if 'v' in kwargs: v = kwargs['v']
        return pyspline.point_surface(x0,self.tu,self.tv,self.ku,self.kv,self.coef,
                                              Niter,eps1,eps2,u,v)
      

    def projectCurve(self,curve,Niter=25,eps1=1e-6,eps2=1e-6,*args,**kwargs):

        '''Project a point x0 onto the surface. i.e. Find the point on the
        surface that minimizes the distance from x0 to
        surface(u,v). See the fortran function project_point_surface.90'''

        # We will use a starting point u0,v0 if given
        u = -1.0
        v = -1.0
        s = -1.0
        if 'u' in kwargs: u = kwargs['u']
        if 'v' in kwargs: v = kwargs['v']
        if 's' in kwargs: s = kwargs['s']
        return pyspline.curve_surface(curve.t,curve.k,curve.coef,
                                      self.tu,self.tv,self.ku,self.kv,self.coef,
                                      Niter,eps1,eps2,u,v,s)
   
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
            val = self.getValueV(u,v)

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
            
            u_plot = linspace(self.range[0],self.range[1],nu)
            v_plot = linspace(self.range[2],self.range[3],nv)
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

    def writeIGES_directory(self,handle,Dcount,Pcount):

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

    def writeIGES_parameters(self,handle,Pcount,counter):
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
            X, real, array, size(len(s),nDim): Array of data to fit OR
                 x=<vals>,y=<vals>,z=<vals> where each is of length s
            s, real, array: (OPTIONAL for nDim >=2 ) Array of s values 

        arg[0] = \'lms\': Create an instance of the spline class using a lms spline.
       **kwargs MUST contain the following information:
            k, integer: Order for spline
            s, real, array: (OPTIONAL for nDIM >=2) Array of s values 
            X, real, array, size(len(s),nDim): Array of data to fit OR
                 x=<vals>,y=<vals>,z=<vals> where each is of length s
            Nctl, integer: The number of control points
    '''
        
        task = args[0]
        if task == 'create':
            assert 'k' in kwargs and 't' in kwargs and 'coef' in kwargs and 'range' in kwargs, \
                'Error: k,t,coef, and range MUST be defined for task=\'create\''
            
            self.s = None
            self.X = None
            self.N = None
            self.k = int(kwargs['k'] )
            self.t = array(kwargs['t'])
            self.coef = array(kwargs['coef'])
            self.Nctl = self.coef.shape[0]
            self.orig_data = False
            self.smin = kwargs['range'][0]
            self.smax = kwargs['range'][1]
            self.nDim = self.coef.shape[1]-1
            self._calcGrevillePts()
            return
             
        elif task=='interpolate' or task == 'lms':
            if task == 'interpolate':
                assert 'k' in kwargs and ('X' in kwargs or 'x' in kwargs),\
                'Error: k, and X (or x=,y=) MUST be defined for task \'interpolate\''
            else:
                assert 'k' in kwargs and ('X' in kwargs or 'x' in kwargs) \
                and 'Nctl' in kwargs ,'Error: k, X and Nctl MUST be defined for task lms'
            # end if

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

            self.N = self.X.shape[0]
            self.k = kwargs['k']
            assert self.k>=2 and self.k<9,'The value of k must be between 2 and 9'
            
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
            
            if task == 'interpolate':
                self.Nctl = self.N      
                self.niter = 1
            else:
                self.Nctl = kwargs['Nctl']
                if self.Nctl > self.N:
                    self.Nctl  = self.N
                # end if
            # end if

            # Sanity check to make sure k is less than N
            if self.N <= self.k:
                self.k = self.N
            # end if

            # Generate the knot vector
            self.t = pyspline.knots(self.s,self.Nctl,self.k)
            self.coef = zeros((self.Nctl,self.nDim),'d')
            self.s,self.coef,rms = pyspline.compute_curve(self.s,self.X,self.t,self.k,self.coef,
                                                   self.niter,self.rel_tol)
            mpiPrint('Fitted with rms:%f'%rms)
            self._calcGrevillePoints()

    def runParameterCorrection(self,niter,rel_tol=1e-5):
        # Run more parameter correction
        self.s,self.coef,rms = pyspline.compute_curve(self.s,self.X,self.t,self.k,self.coef,
                                               niter,self.rel_tol)
        mpiPrint('ReFitted with rms: %f'%rms)

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
      
    def getLength(self):
        # Get the length of the actual Curve
        # Use the greville points for the positions
        
        # We should do this with exact gaussian integration 

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
        return self.getValue(s)
    
    def getValue(self,s):
        '''Get The value of the spline'''
        s = array(s)
        if len(s.shape) == 0:
            return pyspline.eval_curve(s,self.t,self.k,self.coef)
        elif len(s.shape) == 1:
            return pyspline.eval_curve_v(s,self.t,self.k,self.coef)

    def getDerivative(self,s):
        '''Get the value of the derivatve of the spline'''
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
        if self.orig_data:
            self._writeTecplot1D(handle,'orig_data',self.X)
        # end if

        return
    def _writeTecplotCurve(self,handle,*args,**kwargs):
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
	
    # Run a Simple Test Case
    print 'Testing pyspline...\n'
    print 'There is an example in the ./example directory'



# Functions to Eliminate Possibly
#   def checkCoef(self):
#         '''This function determines if there are crossed control points which
#         could result in an invalid mesh'''
#         # This function can easily migrate to fortran
#         counter = 0
#         # Loop over number of "elements or patches"
#         ustart = 0
#         uend = self.Nctlu - 1
#         vstart = 0
#         vend = self.Nctlv - 1


#         for i in xrange(ustart,uend):
#             for j in xrange(vstart,vend):
#                 # First find normal of "patch" by taking cross product
#                 # of the diagonals
#                 d1 = self.coef[i+1,j+1]-self.coef[i,j]
#                 d2 = self.coef[i,j+1]-self.coef[i+1,j]
#                 normal = cross(d1,d2)
#                 normal /= sqrt(dot(normal,normal))
#                 if sqrt(dot(normal,normal)) == 0:
#                     print 'ha ha '
#                     sys.exit(0)
#                 v1 = self.coef[i+1,j  ]-self.coef[i  ,j  ]
#                 v2 = self.coef[i+1,j+1]-self.coef[i+1,j  ]
#                 v3 = self.coef[i  ,j+1]-self.coef[i+1,j+1]
#                 v4 = self.coef[i  ,j  ]-self.coef[i  ,j+1]
#                 if sqrt(dot(v1,v1)) == 0:
#                     v1 = zeros(3)
#                 else:
#                     v1 /=sqrt(dot(v1,v1))
#                 if sqrt(dot(v2,v2)) == 0:
#                     v2 = zeros(3)
#                 else:
#                     v2 /=sqrt(dot(v2,v2))
#                 if sqrt(dot(v3,v3)) == 0:
#                     v3 = zeros(3)
#                 else:
#                     v3 /=sqrt(dot(v3,v3))
#                 if sqrt(dot(v4,v4)) == 0:
#                     v4 = zeros(3)
#                 else:
#                     v4 /=sqrt(dot(v4,v4))
          

#                 # Now cross normal with all the vectors
#                 values = zeros(4)
#                 values[0] = dot(cross(v1,v2),normal)
#                 values[1] = dot(cross(v2,v3),normal)
#                 values[2] = dot(cross(v3,v4),normal)
#                 values[3] = dot(cross(v4,v1),normal)
# #                print values
#                 for ii in xrange(4):
#                     if abs(values[ii]) < 1e-3:
#                         values[ii] = 0
#                     # end if
#                 # end for

#                 if values[0] >=0 and values[1] >= 0 and values[2]>=0 and values[3]>=0:
#                     pass
#                 else:
#                     if values[0]<=0 and values[1]<=0 and values[2]<=0 and values[3]<=0:
#                         pass
#                     else:
#                         print 'values:',values,i,j
# #                         print 'normal:',normal
# #                         print 'v:',v1,v2,v3,v4
#                         counter += 1
#                     # end if
#                 # end if
#             # end for
#         # end for
#         return counter 


#    if 'dx1' and 'dx2' in kwargs:
#                 # We have defined a tangent vector at each end

#                 #print 'dx1 and dx2 defined'
                
#                 self.Nctl = self.N + 2
#                 dx1 = kwargs['dx1']
#                 dx2 = kwargs['dx2']
#                 ibcl = 1
#                 ibcr = 1
#                 kntopt = 1
#                 assert len(dx1)==len(dx2)==self.nDim,'The length of the \
# derivative vectors must match the spatial dimension of the curve'

#                 if self.nDim > 1:
#                     self.coef = zeros((self.Nctl,self.nDim),'d')

#                     for idim in xrange(self.nDim):
#                         fbcl = dx1[idim]
#                         fbcr = dx2[idim]

#                         self.t,self.coef[:,idim],k = \
#                             pyspline.bint4(self.s,self.X[:,idim],\
#                                                     ibcl,ibcr,fbcl,fbcr,kntopt)

#                 else:
#                     print 'do something'
#                 # end if




#             # Calculate the knot vector and Jacobian

#             self._calcJacobian()

#             # Lets do a lms 
#             timeA = time.time()
#             self.coef = zeros([self.Nctlu,self.Nctlv,self.nDim],'d')

#             res = numpy.linalg.lstsq(self.J,self.X.reshape[self.Nu*self.Nv,self.nDim])
#             for idim in xrange(self.nDim):
#                 self.coef[:,:,idim] = reshape(res[0][idim],self.Nctlu,self.Nctlv)
#             # end for
                
#             mpiPrint(' LMS Fit Time: %6.5f s\n'%(time.time()-timeA),self.NO_PRINT)
            
#             # Finally, set edges as linear splines
#             self._setEdgeCurves()



#         self.coef = zeros((self.Nctlu*self.Nctlv,self.nDim))
        
#         if USE_PETSC:

#             J = PETSc.Mat()
#             if PETSC_MAJOR_VERSION == 1:
#                 J.createAIJ([self.Nu*self.Nv,self.Nctlu*self.Nctlv],nnz=self.ku*self.kv,comm=PETSc.COMM_SELF)
#             elif PETSC_MAJOR_VERSION == 0:
#                 J.createSeqAIJ([self.Nu*self.Nv,self.Nctlu*self.Nctlv],nz=self.ku*self.kv)

#             rhs = PETSc.Vec().createSeq(self.Nu*self.Nv)
#             sol = PETSc.Vec().createSeq(self.Nctlu*self.Nctlv)

#             ksp = PETSc.KSP()
#             ksp.create(PETSc.COMM_WORLD)
#             ksp.getPC().setType('none')
#             ksp.setType('lsqr')
#             ksp.setConvergenceTest(self._converge_test)
#             nrow = self.Nu*self.Nv
#             ncol = self.Nctlu*self.Nctlv
#             for iter in xrange(niter):
#                 J.zeroEntries()
#                 rows,cols,vals = pyspline.surface_jacobian_linear(
#                     self.U,self.V,self.tu,self.tv,self.ku,self.kv,
#                     self.Nctlu,self.Nctlv)

#                 for i in xrange(len(rows)):
#                     J.setValue(int(rows[i]),int(cols[i]),vals[i])#,PETSc.InsertMode.INSERT_VALUES)
                    
#                 # end for
#                 J.assemble()
#                 ksp.setOperators(J)

#                 for idim in xrange(self.nDim): # Solve for each RHS
#                     rhs[:] = self.X[:,:,idim].flatten()
#                     ksp.solve(rhs, sol)
#                     self.coef[:,idim] = sol[:]

#                 # end for
#                 self.U,self.V,rms = pyspline.surface_para_corr(
#                     self.tu,self.tv,self.ku,self.kv,self.U,self.V,\
#                     self.coef.reshape((self.Nctlu,self.Nctlv,self.nDim)),self.X)

#                 # Convergence Checks
#                 if iter == 0:
#                     self.old_rms = rms
#                     rms0 = rms
#                     if rms < 1e-12:
#                         break
#                 else:
#                     rel_change = abs((rms-self.old_rms)/self.old_rms)
#                     if rel_change > 1:
#                         print 'Fucked up'
#                         break
#                     self.old_rms = rms
#                     if rel_change < self.rel_tol or rms <1e-12:
#                         break
#                     # end if
                    
#                 # end if
#             # end for
#         else: # Use Numpy
#                pass

#         self.coef = self.coef.reshape((self.Nctlu,self.Nctlv,self.nDim))
#         print 'Succefully fit: rms0 = %g, rms_final = %g,iterations %d:'%(rms0,rms,iter+1)
        
#        return

#     def _converge_test(self,ksp,iter,rnorm):
#         # Do a relative norm and a happy breakdown absolute norm
#         if iter== 0:
#             self.old_rnorm = rnorm
#             if rnorm < 1e-12:
#                 return 1
#             # end if 
#         else:
#             val = abs((rnorm-self.old_rnorm)/self.old_rnorm)
#             self.old_rnorm = rnorm
#             if val < 1e-12 or rnorm < 1e-12:
#                 return 1
#             # end if
#         # end if

    # If Using PETSc
#             self.coef = zeros((self.Nctl,self.nDim))
                 
#             if USE_PETSC:
#                 J = PETSc.Mat()
#                 if PETSC_MAJOR_VERSION == 1:
#                     J.createAIJ([self.N,self.Nctl],nnz=self.k,comm=PETSc.COMM_SELF)
#                 elif PETSC_MAJOR_VERSION == 0:
#                     J.createSeqAIJ([self.N,self.Nctl],nz=self.k)
#                 # end if
                
#                 rhs = PETSc.Vec().createSeq(self.N)
#                 temp = PETSc.Vec().createSeq(self.Nctl)
                
#                 ksp = PETSc.KSP()
#                 ksp.create(PETSc.COMM_WORLD)
#                 ksp.getPC().setType('none')
#                 ksp.setType('lsqr')
#                 ksp.setConvergenceTest(self._converge_test)
#                 for j in xrange(self.niter):
#                     J.zeroEntries()
#                     rows,cols,vals = pyspline.curve_jacobian_linear(
#                         self.t,self.k,self.s,self.Nctl)
#                     for i in xrange(len(rows)):
#                         J.setValue(int(rows[i]),int(cols[i]), vals[i],PETSc.InsertMode.INSERT_VALUES)
#                     # end for
#                     J.assemble()
#                     ksp.setOperators(J)

#                     for idim in xrange(self.nDim): # Solve for each RHS
#                         rhs[:] = self.X[:,idim]
#                         ksp.solve(rhs, temp)
#                         self.coef[:,idim] = temp
#                     # end for

#                     # Run Hoschek Parameter Correction
#                     self.s,rms = pyspline.curve_para_corr(
#                         self.t,self.k,self.s,self.coef,length,self.X)

#                     # Convergence Checks
#                     if j == 0:
#                         self.old_rms = rms
#                         rms0 = rms
#                     else:
#                         val = abs((rms-self.old_rms)/self.old_rms)
#                         self.old_rms = rms
#                         if val < self.rel_tol or rms <1e-12:
#                             break
#                         # end if
#                     # end if
#                 # end for
#             else: # Use Numpy
#                 print 'Not implemented yet'
#                 pass

#             # end if
#         # end if
       
#         print 'Succefully fit: rms0 = %g, rms_final = %g,iterations %d:'%(rms0,rms,j+1)
                 
#         return


#     def _converge_test(self,ksp,iter,rnorm):
#         # Do a relative norm and a happy breakdown absolute norm
#         if iter== 0:
#             self.old_rnorm = rnorm
#             if rnorm < 1e-12:
#                 return 1
#             # end if 
#         else:
#             val = abs((rnorm-self.old_rnorm)/self.old_rnorm)
#             self.old_rnorm = rnorm
#             if val < 1e-12 or rnorm < 1e-12:
#                 return 1
#             # end if
#         # end if
    # --- From curve class ----
   #  def writeTecplot(self,handle):
#         '''Output this line\'s data to a open file handle \'handle\' '''

#         if self.orig_data:
#             handle.write('Zone T=%s I=%d \n'%('orig_data',self.N))
#             handle.write('DATAPACKING=POINT\n')
#             for i in xrange(self.N):
#                 for idim in xrange(self.nDim):
#                     handle.write('%f '%(real(self.X[i,idim])))
#                 # end for
#                 handle.write('\n')
#         # end if

#         if self.orig_data:
#             s_plot = self.s
#         else:
#             s_plot = linspace(self.range[0],self.range[1],25)
#         # end if 
        
#         # Dump re-interpolated spline
#         handle.write('Zone T=%s I=%d \n'%('interpolated',len(s_plot)))
#         handle.write('DATAPACKING=POINT\n')
#         for i in xrange(len(s_plot)):
#             for idim in xrange(self.nDim):
#                 handle.write('%f '%(\
#                         pyspline_real.bvalu(self.t,\
#                                            self.coef[:,idim],\
#                                            self.k,0,s_plot[i])))
#             # end for 
#             handle.write('\n')
#         # end for 

#         # Dump Control Points (Always have these :-) ) 
#         handle.write('Zone T=%s I = %d\n'%('control_pts',self.Nctl))
#         handle.write('DATAPACKING=POINT\n')

#         for i in xrange(self.Nctl):
#             for idim in xrange(self.nDim):
#                 handle.write('%f '%(real(self.coef[i,idim])))
#             # end for
#             handle.write('\n')
#         # end for

#         return
#     def findUV(self,x0,r,u0=0.5,v0=0.5,NO_PRINT=False):

  # Generic algorithim -- Min distance between curve and surface
                                
