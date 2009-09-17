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
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot, cross, mod, floor, mat

import numpy.linalg
from numpy.linalg import lstsq

import pyspline_cs
import pyspline as pyspline_real
# =============================================================================
# pySpline class
# =============================================================================

class surf_spline():

    def __init__(self,task='create',*args,**kwargs):

        '''Create an instance of a b-spline surface. There are three ways to
        initialize the class as determined by the task flag:

        task = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for u
            kv, integer: Order for v
            Nctlu, integer: Number of u control points
            Nctlv, integer: Number of v control points
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control point 
                                                     values

        task = \'interpolate\': Create an instance of the spline class
        by using an interpolating spline to given data points. **kwarg
        MUST contain the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: list of u values 
            v, real, array: list of v values
            X, real, array, size(len(u),len(v),nDim):Array of data points to fit

        task = \'lms\': Create an instance of the spline class using a
        Least-Mean-Squares fit to the spline. . **kwargs MUST contain
        the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: list of u values 
            v, real, array: list of v values
            X, real, array, size(len(u),len(v),nDim):Array of data points to fit
'''
        if 'no_print' in kwargs:
            self.NO_PRINT = kwargs['no_print']
        else:
            self.NO_PRINT = False
        # end if      

        if not self.NO_PRINT:
            sys.stdout.write('pySpline Type: %s. '%(task))

        self.pyspline_real = pyspline_real # Real version of spline library
        self.pyspline_cs = pyspline_cs # Complex version of spline library
        # sef.
        if 'complex' in kwargs:
            self.pyspline = self.pyspline_cs # Use the complex version for everything
            self.dtype = 'D'
        else:
            self.pyspline = self.pyspline_real
            self.dtype = 'd'
        # end if

        # Defaults
        self.task = task
        self.X    = None
        self.u    = None
        self.v    = None
        self.Nu   = None
        self.Nv   = None
        self.orig_data = False

        if task == 'create':
            assert 'ku' in kwargs and 'kv' in kwargs  and 'tu' in kwargs \
                and 'tv' in kwargs and 'coef' in kwargs and 'range' in kwargs,\
                'Error: ku,kv,tu,tv,coef and range MUST be defined for task=\
\'create\''
            if not self.NO_PRINT:
                sys.stdout.write('\n')
            self.ku = kwargs['ku'] 
            self.kv = kwargs['kv']
            self.tu = array(kwargs['tu'],self.dtype)
            self.tv = array(kwargs['tv'],self.dtype)
            self.coef = array(kwargs['coef'],self.dtype)
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.Nctl  = self.Nctlu*self.Nctlv
            self.range = array(kwargs['range'],self.dtype)
            self.nDim = self.coef.shape[2]
            return
     
        elif task == 'interpolate' or task == 'lms':
            # Do some checking on the number of control points
            if task == 'lms':
                assert 'ku' in kwargs and 'kv' in kwargs and \
                    'Nctlu' in kwargs and 'Nctlv' in kwargs and 'X' in kwargs, \
                    'Error: ku,kv,Nctlu,Nctlv and X MUST be defined for task lms'
                self.Nctlu = kwargs['Nctlu']
                self.Nctlv = kwargs['Nctlv']
            else:
                assert 'ku' in kwargs and 'kv' in kwargs and \
                   'X' in kwargs,'Error: ku,kv,and X MUST be defined for task \
interpolate'
                self.Nctlu = kwargs['X'].shape[0]
                self.Nctlv = kwargs['X'].shape[1]
            # end if
                
            self.task = task
            self.X  = kwargs['X']
            self.orig_data = True
            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.nDim  = self.X.shape[2]
           
            self.ku = kwargs['ku']
            self.kv = kwargs['kv']

            # Sanity Check on Inputs
            if self.Nctlu > self.Nu:
                self.Nctlu  = self.Nu
                print 'Warning: Number of control points in u has been capped \
to Nu: Nctlu = %d'%self.Nctlu
            if self.Nctlv > self.Nv:
                self.Nctlv = self.Nv
                print 'Warning: Number of control points in v has been capped \
to Nv: Nctlv = %d'%self.Nctlv
            # end if

            self.Nctl  = self.Nctlu*self.Nctlv

            # Sanity check to make sure k is less than N
            if self.Nu < self.ku:
                self.ku = self.Nu
                print 'Warning: ku has been set to Nu. ku is now %d'%(self.ku)
            if self.Nv < self.kv:
                print 'Warning: kv has been set to Nv. kv is now %d'%(self.kv)
                self.kv = self.Nv

            if 'u' in kwargs and 'v' in kwargs:
                self.u = kwargs['u']
                self.v = kwargs['v']
            else:
                if self.nDim == 3:
                    self._calcParameterization()
                else:
                    print 'Automatric parameterization of ONLY available\
 for spatial data in 3 dimensions. Please supply u and v key word arguments\
 otherwise.'
                    sys.exit(1)
                # end if
            # end if
            
            self.range = array([self.u[0],self.u[-1],self.v[0],self.v[-1]])

            # Set a linear version of U and V
            [self.V, self.U] = meshgrid(self.v,self.u)

           #Calculate the knot vector and Jacobian
            if not self.NO_PRINT:
                sys.stdout.write(' Calculating: knots, ')
            self._calcKnots()
            if not self.NO_PRINT:
                sys.stdout.write(' jacobian, ')
            self._calcJacobian()

#             # Lets do a lms 
            timeA = time.time()
            self.coef = zeros([self.Nctlu,self.Nctlv,self.nDim],self.dtype)
            for idim in xrange(self.nDim):
                self.coef[:,:,idim] =\
                    reshape(lstsq(self.J,self.X[:,:,idim].flatten())[0]\
                                ,[self.Nctlu,self.Nctlv])
            if not self.NO_PRINT:
                sys.stdout.write(' LMS Fit Time: %6.5f s\n'%(time.time()-timeA))
            
            return 

        print 'Error: task is not understood. Type must be \'lms\',\
\'interpolate\' or \'create\''
        sys.exit(1)
        
        return

    def recompute(self):
        '''Recompute the surface if the knot vector has changed:'''

        assert self.orig_data,'Only surface initialization with original \
data can be recomputed'

        if self.task == 'lms' or self.task == 'interpolate':
            self._calcJacobian()
            self.coef = zeros([self.Nctlu,self.Nctlv,self.nDim],self.dtype)
            for idim in xrange(self.nDim):
                self.coef[:,:,idim] =\
                    reshape(lstsq(self.J,self.X[:,:,idim].flatten())[0]\
                                ,[self.Nctlu,self.Nctlv])
            # end for
        # end if

    def _calcParameterization(self):

        u = zeros(self.Nu,self.dtype)
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for j in xrange(self.Nv): 
            temp = zeros(self.Nu,self.dtype)

            for i in xrange(self.Nu-1):
                temp[i+1] = temp[i] + sqrt((self.X[i+1,j,0]-self.X[i,j,0])**2 +\
                                           (self.X[i+1,j,1]-self.X[i,j,1])**2 +\
                                           (self.X[i+1,j,2]-self.X[i,j,2])**2)
            # end for

            if temp[-1] == 0: # We have a singular point
                singular_counter += 1
                temp[:] = 0.0
            else:
                temp /= temp[-1]
            # end if

            u += temp #accumulate the u-parameter calcs for each j
            
        # end for 
        u/=(self.Nv-singular_counter) #divide by the number of 'j's we had
        self.u = u
        
        v = zeros(self.Nv,self.dtype)
        singular_counter = 0
        # loop over each v, and average the 'u' parameter 
        for i in xrange(self.Nu): 
            temp = zeros(self.Nv,self.dtype)
            for j in xrange(self.Nv-1):
                temp[j+1] = temp[j] + sqrt((self.X[i,j+1,0]-self.X[i,j,0])**2 +\
                                           (self.X[i,j+1,1]-self.X[i,j,1])**2 +\
                                           (self.X[i,j+1,2]-self.X[i,j,2])**2)
            # end for
            if temp[-1] == 0: #We have a singular point
                singular_counter += 1
                temp[:] = 0.0
            else:
                temp /= temp[-1]
            #end if 

            v += temp #accumulate the v-parameter calcs for each i
        # end for 
        v/=(self.Nu-singular_counter) #divide by the number of 'i's we had
        self.v = v

        return

    def _calcKnots(self):

        '''Find the initial knot vector for this problem'''

        # ADD COMPLEX VERSION OF KNOTS
        #self.tu = self.pyspline.knots(self.u,self.Nctlu,self.ku)
        #self.tv = self.pyspline.knots(self.v,self.Nctlv,self.kv)

        self.tu = pyspline_real.knots(self.u,self.Nctlu,self.ku)
        self.tv = pyspline_real.knots(self.v,self.Nctlv,self.kv)

        return

    def _calcJacobian(self):

        # Calculate the jacobian J, for fixed t and s

        self.J = zeros([self.Nu*self.Nv,self.Nctlu*self.Nctlv],self.dtype)
        ctl = zeros([self.Nctlu,self.Nctlv],self.dtype)
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                ctl[i,j] += 1
                self.J[:,i*self.Nctlv + j] = \
                    self.pyspline.b2valv(\
                    self.U.flatten(),self.V.flatten(),0,0,\
                        self.tu,self.tv,self.ku,self.kv,ctl)
                ctl[i,j] -= 1
                # end for
            # end for 
        # end for
        return

    def checkCoef(self):
        '''This function determines if there are crossed control points which
        could result in an invalid mesh'''
        # This function can easily migrate to fortran
        counter = 0
        # Loop over number of "elements or patches"
        ustart = 0
        uend = self.Nctlu - 1
        vstart = 0
        vend = self.Nctlv - 1


        for i in xrange(ustart,uend):
            for j in xrange(vstart,vend):
                # First find normal of "patch" by taking cross product
                # of the diagonals
                d1 = self.coef[i+1,j+1]-self.coef[i,j]
                d2 = self.coef[i,j+1]-self.coef[i+1,j]
                normal = cross(d1,d2)
                normal /= sqrt(dot(normal,normal))
                if sqrt(dot(normal,normal)) == 0:
                    print 'ha ha '
                    sys.exit(0)
                v1 = self.coef[i+1,j  ]-self.coef[i  ,j  ]
                v2 = self.coef[i+1,j+1]-self.coef[i+1,j  ]
                v3 = self.coef[i  ,j+1]-self.coef[i+1,j+1]
                v4 = self.coef[i  ,j  ]-self.coef[i  ,j+1]
                if sqrt(dot(v1,v1)) == 0:
                    v1 = zeros(3)
                else:
                    v1 /=sqrt(dot(v1,v1))
                if sqrt(dot(v2,v2)) == 0:
                    v2 = zeros(3)
                else:
                    v2 /=sqrt(dot(v2,v2))
                if sqrt(dot(v3,v3)) == 0:
                    v3 = zeros(3)
                else:
                    v3 /=sqrt(dot(v3,v3))
                if sqrt(dot(v4,v4)) == 0:
                    v4 = zeros(3)
                else:
                    v4 /=sqrt(dot(v4,v4))
          

                # Now cross normal with all the vectors
                values = zeros(4)
                values[0] = dot(cross(v1,v2),normal)
                values[1] = dot(cross(v2,v3),normal)
                values[2] = dot(cross(v3,v4),normal)
                values[3] = dot(cross(v4,v1),normal)
#                print values
                for ii in xrange(4):
                    if abs(values[ii]) < 1e-3:
                        values[ii] = 0
                    # end if
                # end for

                if values[0] >=0 and values[1] >= 0 and values[2]>=0 and values[3]>=0:
                    pass
                else:
                    if values[0]<=0 and values[1]<=0 and values[2]<=0 and values[3]<=0:
                        pass
                    else:
                        print 'values:',values,i,j
#                         print 'normal:',normal
#                         print 'v:',v1,v2,v3,v4
                        counter += 1
                    # end if
                # end if
            # end for
        # end for
        return counter 


    def calcPtDeriv(self,u,v,i,j):
        '''Calc the derivative of point u,v at control point i,j'''
        coef = zeros((self.Nctlu,self.Nctlv))
        coef[i,j] = 1.0
        x = self.pyspline.b2val(\
            u,v,0,0,self.tu,self.tv,self.ku,self.kv,coef)
        
        return x

    def getValueEdge(self,edge,s):
        '''Get the value of the spline on edge, edge=0,1,2,3'''

        if edge == 0:
            return self.getValue(s,self.range[2])
        elif edge == 1:
            return self.getValue(s,self.range[3])
        elif edge == 2:
            return self.getValue(self.range[0],s)
        elif edge ==3:
            return self.getValue(self.range[1],s)
        else:
            print 'Edge must be between 0 and 3'
            sys.exit(1)
            return
        #end if 

    def getValueCorner(self,corner):
        '''Get the value of the spline on corner i '''

        if corner == 0:
            return self.getValue(self.range[0],self.range[2])
        elif corner == 1:
            return self.getValue(self.range[1],self.range[2])
        elif corner == 2:
            return self.getValue(self.range[1],self.range[3])
        elif corner ==3:
            return self.getValue(self.range[0],self.range[3])
        else:
            print 'Corner must be between 0 and 3'
            sys.exit(1)
            return
        #end if

    def getOrigValuesEdge(self,edge):
        ''' Get the values of the original data on edge. edge = 0,1,2,3'''
        if self.orig_data == False:
            print 'Error: No original data exists for this surface. The \
initialization type for this spline class was \'create\''
            sys.exit(1)
        # end if

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
        else:
            print 'Error: edge must be between 0 and 3'
            sys.exit(1)
        # end if

    def getOrigValueCorner(self,node):
        ''' Get the values of the original data on edge. edge = 0,1,2,3'''
        if self.orig_data == False:
            print 'Error: No original data exists for this surface. The \
initialization type for this spline class was \'create\''
            sys.exit(1)
        # end if

        if node == 0:
            return self.X[0,0]
        elif node == 1:
            return self.X[-1,0]
        elif node == 2:
            return self.X[0,-1]
        elif node == 3:
            return self.X[-1,-1]

    def checkDegenerateEdge(self,edge):
        
        '''Check to see if the values on edge \'edge\' are degenerate'''
        
        degen = False

        if edge == 0:
            values = self.X[:,0]
        elif edge == 1:
            values = self.X[:,-1]
        elif edge == 2:
            values = self.X[0,:]
        elif edge == 3:
            values = self.X[-1,:]
        else:
            print 'Error: Edge must be between 0 and 3'
            sys.exit(1)
        # end if

        # Now loop over edge to see if its degenerate
        N = len(values)
        length = zeros(N)
        for i in xrange(N-1):
            length[i+1] = length[i] + self._e_dist(values[i+1],values[i])
        # end for
        
        if abs(length[-1]) < 1e-12:
            degen = True
        # end if

        return degen,values[0]


    def getNormal(self,u,v):
        '''Get the normalized normal at the surface point u,v'''
        if self.nDim == 3:
            du,dv = self.getDerivative(u,v)

            # Now cross and normalize
            n = zeros(3,self.dtype)
            n[0] = du[1]*dv[2]-du[2]*dv[1]
            n[1] = du[2]*dv[0]-du[0]*dv[2]
            n[2] = du[0]*dv[1]-du[1]*dv[0]
            
            n/= sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
        else:
            print 'Warning: getNormal is only defined for 3 spatial dimension'
            sys.exit(1)
        # end if
        return n

    def getValue(self,u,v,x=None):
        
        '''Get the value of the spline at point u,v'''
        if x == None:
            x = zeros((self.nDim),self.dtype)
        # end if
        for idim in xrange(self.nDim):
            x[idim] = self.pyspline.b2val(\
                u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        
        return x

    def getValueV(self,u,v,x=None):
        '''Get the value of a spline at vector of points u,v'''
        # Note: If the user is passing in x, it must already be the
        # right shape/size
        assert u.shape == v.shape, 'u and v must be the same length'
        if x == None:
            x = zeros((len(u),self.nDim),self.dtype)
        # end if

        for idim in xrange(self.nDim):
            x[:,idim] = self.pyspline.b2valv(\
                u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return x

    def getValueM(self,u,v,x=None):
        '''Get the value of a spline at matrix of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same shape'
        if x == None:
            x = zeros((u.shape[0],u.shape[1],self.nDim),self.dtype)
        # end if
        for idim in xrange(self.nDim):
            x[:,:,idim] = self.pyspline.b2valm(\
                u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return x

    def getDerivative(self,u,v):
        
        '''Get the value of the derivative of spline at point u,v'''
        du = zeros(self.nDim,self.dtype)
        dv = zeros(self.nDim,self.dtype)
        for idim in xrange(self.nDim):
            du[idim] = self.pyspline.b2val(\
                u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            dv[idim] = self.pyspline.b2val(\
                u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        return du,dv


    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((self.nDim,2),self.dtype)
        
        for idim in xrange(self.nDim):
            J[idim,0] = self.pyspline.b2val(\
                u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            J[idim,1] = self.pyspline.b2val(\
                u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return J

    def projectPoint(self,x0,u0=0.5,v0=0.5,Niter=25,tol=1e-6):

        '''Project a point x0 onto the surface. i.e. Find the point on the
        surface that minimizes the distance from x0 to
        surface(u,v). See the fortran function project point.'''

        # We will use a starting point u0,v0 if given

        u0,v0,D,converged = self.pyspline.projectpoint(\
            self.coef,self.ku,self.kv,self.tu,self.tv,x0,u0,v0,Niter,tol)

        return u0,v0,D,converged

    def findUV(self,x0,r,u0=0.5,v0=0.5):
        ''' Try to find the parametric u-v coordinate of the spline
        which coorsponds to the intersection of the directed vector 
        v = x0 + r*s where x0 is a basepoint, r  is a direction vector
        and s is the distance along the vector.
        
        If possible, both intersections are attempted with the first 
        coorsponding to the first intersection in the direction of r
        
        Input: 
        
        x0: array, length 3: The base of the vector
        r : array, length 3: The direction of the vector

        u0: scalar: The guess for the u coordinate of the intersection
        v0: scalar: The guess for the v coordinate of the intersection

        '''
        maxIter = 25

        u = u0
        v = v0
        s = 0 #Start at the basepoint

        for iter in xrange(maxIter):

            #just in case, force crop u,v to [-1,1]
            if u<-1: u = -1
            if u>1 : u =  1
            if v<-1: v = -1
            if v>1 : v =  1

            x = self.getValue(u,v) #x contains the x,y,z coordinates 

            f = mat(zeros((3,1),self.dtype))
            f[0] = x[0]-(x0[0]+r[0]*s)
            f[1] = x[1]-(x0[1]+r[1]*s)
            f[2] = x[2]-(x0[2]+r[2]*s)

            J = self.getJacobian(u,v)
            A = mat(zeros((3,3),self.dtype))
            A[:,0:2] = J
            A[0,2]   = -r[0]
            A[1,2]   = -r[1]
            A[2,2]   = -r[2]

            x_up = numpy.linalg.solve(A,-f)

            # Add a little error checking here:
            
            if u + x_up[0] < -1 or u + x_up[0] > 1 or \
               v + x_up[1] < -1 or v + x_up[1] > 1:
                #Cut the size of the step in 1/2
                x_up /= 2

            u = u + x_up[0]
            v = v + x_up[1]
            s = s + x_up[2]

            if numpy.linalg.norm(x_up) < 1e-12:
                return u,v,x

        # end for

        print 'Warning: Newton Iteration for u,v,s did not converge:'
        print 'u = %f, v = %f, s = %f\n'%(u,v,s)
        print 'Norm of update:',numpy.linalg.norm(x_up)

        return u,v,x


    def writeTecplotSurface(self,handle,size=None):
        '''Output this surface\'s data to a open file handle \'handle\''''

        MAX_SIZE = 100
        MIN_SIZE = 5

        if self.orig_data:
            handle.write('Zone T=%s I=%d J=%d\n'%('orig_data',self.Nu,self.Nv))
            handle.write('DATAPACKING=POINT\n')
            for j in xrange(self.Nv):
                for i in xrange(self.Nu):
                    handle.write('%f %f %f \n'%(\
                            self.X[i,j,0].astype('d'),self.X[i,j,1].astype('d')\
                                ,self.X[i,j,2].astype('d')))
                # end for
            # end for 
        # end if
                
        if size == None:
            u_plot = linspace(self.range[0],self.range[1],50).astype('d')
            v_plot = linspace(self.range[2],self.range[3],50).astype('d')
        else:

            # Cheaply calculate the length of each size of the surface
            # to determine its lenght and then use the size parameter 
            # to determine the number of points to use

            u = array([0,0.5,1,1,1,0.5,0,0])
            v = array([0,0,0,0.5,1,1,1,0.5])
            val = self.getValueV(u,v)

            u_len = (self._e_dist(val[0],val[1])+self._e_dist(val[1],val[2])+\
                         self._e_dist(val[4],val[5])+\
                         self._e_dist(val[5],val[6]))/2
            
            v_len = (self._e_dist(val[2],val[3])+self._e_dist(val[3],val[4])+\
                         self._e_dist(val[6],val[7])+\
                         self._e_dist(val[7],val[0]))/2
            
            nu=int(floor(real(u_len/size)))
            nv=int(floor(real(v_len/size)))
            
            if nu > MAX_SIZE: nu = MAX_SIZE
            if nu < MIN_SIZE: nu = MIN_SIZE
            if nv > MAX_SIZE: nv = MAX_SIZE
            if nv < MIN_SIZE: nv = MIN_SIZE
            
            u_plot = linspace(self.range[0],self.range[1],nu).astype('d')
            v_plot = linspace(self.range[2],self.range[3],nv).astype('d')

        # Dump re-interpolated surface
        handle.write('Zone T=%s I=%d J = %d\n'\
                         %('interpolated',len(u_plot),len(v_plot)))
        handle.write('DATAPACKING=POINT\n')

        [V_plot,U_plot] = meshgrid(v_plot,u_plot)
        X = self.getValueM( U_plot,V_plot)
        for j in xrange(len(v_plot)):
            for i in xrange(len(u_plot)):
                    handle.write('%f %f %f \n'%(X[i,j,0],X[i,j,1],X[i,j,2]))
                # end if
            # end if
        # end if
        #sys.exit(0)
        # ---------------------------------------------------------------------
        
        # Dump Control Points (Always have these :-) ) 
        handle.write('Zone T=%s I=%d J = %d\n'\
                         %('control_pts',self.Nctlu,self.Nctlv))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                handle.write('%f %f %f \n'%(self.coef[i,j,0].astype('d'),\
                                                self.coef[i,j,1].astype('d'),\
                                                self.coef[i,j,2].astype('d')))
            # end for
        # end for 

        return


    def writeTecplotEdge(self,handle,edge,*args,**kwargs):
        '''Dump out a linear zone along edge used for visualizing edge
        connections'''

        N = 25
        if 'name' in kwargs:
            handle.write('Zone T=%s I=%d\n'%(kwargs['name'],N))
        else:
            handle.write('Zone T=%s I=%d\n'%('edge',N))
        handle.write('DATAPACKING=POINT\n')
        s = linspace(0,1,N)
        for i in xrange(N):
            value = self.getValueEdge(edge,s[i]).astype('d')
            handle.write('%f %f %f \n'%(value[0],value[1],value[2]))
        # end for
        return

    def writeDirections(self,handle,isurf):
        '''Write out and indication of the surface direction'''
        handle.write('Zone T=\"surface%d direction\" I=4\n'%(isurf))
        if self.Nctlu >= 3 and self.Nctlv >=3:
            handle.write('%f,%f,%f \n'%(self.coef[1,2,0],self.coef[1,2,1],self.coef[1,2,2]))
            handle.write('%f,%f,%f \n'%(self.coef[1,1,0],self.coef[1,1,1],self.coef[1,1,2]))
            handle.write('%f,%f,%f \n'%(self.coef[2,1,0],self.coef[2,1,1],self.coef[2,1,2]))
            handle.write('%f,%f,%f \n'%(self.coef[3,1,0],self.coef[3,1,1],self.coef[3,1,2]))
        else:
            print 'Not Enough control points to output direction indicator'
        #end if
        return 
                     



    def writeIGES_directory(self,handle,Dcount,Pcount):

        '''Write the IGES file header information (Directory Entry Section)
        for this surface'''
        # A simplier Calc based on cmlib definations The 13 is for the
        # 9 parameters at the start, and 4 at the end. See the IGES
        # 5.3 Manual paraEntries = 13 + Knotsu + Knotsv + Weights +
        # control points
        assert self.nDim == 3, 'Must have 3 dimensions to write to IGES file'
        paraEntries = 13 + (len(self.tu)) + (len(self.tv)) +\
            self.Nctlu*self.Nctlv + 3*self.Nctlu*self.Nctlv+1

        paraLines = paraEntries / 5
        if mod(paraEntries,5) != 0: paraLines += 1

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
            handle.write('%12.6g,'%(real(self.tu[i])))
            if mod(pos_counter,5) == 0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(len(self.tv)):
            pos_counter += 1
            handle.write('%12.6g,'%(real(self.tv[i])))
            if mod(pos_counter,5) == 0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(self.Nctlu*self.Nctlv):
            pos_counter += 1
            handle.write('%12.6g,'%(1.0))
            if mod(pos_counter,5) == 0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for 

        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                for idim in xrange(3):
                    pos_counter += 1
                    handle.write('%12.6g,'%(real(self.coef[i,j,idim])))
                    if mod(pos_counter,5) == 0:
                        handle.write('%7dP%7d\n'%(Pcount,counter))
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
                handle.write('%12.6g,'%(real(self.tu[0])))
            if i == 1:
                handle.write('%12.6g,'%(real(self.tu[-1])))
            if i == 2:
                handle.write('%12.6g,'%(real(self.tv[0])))
            if i == 3:
                # semi-colon for the last entity
                handle.write('%12.6g;'%(real(self.tv[-1]))) 
            if mod(pos_counter,5)==0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            else: # We have to close it up anyway
                if i ==3:
                    for j  in xrange(5-pos_counter):
                        handle.write('%13s'%(' '))
                    # end for
                    pos_counter = 0
                    handle.write('%7dP%7d\n'%(Pcount,counter))
                    counter += 1
                # end if
            # end if
        # end for

        Pcount += 2
        return Pcount,counter

    def _e_dist(self,x1,x2):
        '''Get the eculidean distance between two points'''
        return sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)


class linear_spline():

    def __init__(self,task='create',*args,**kwargs):

        '''Create an instance of a b-spline surface. There are three ways to
        initialize the class as determined by the task flag:

        task = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for u
            kv, integer: Order for v
            Nctlu, integer: Number of u control points
            Nctlv, integer: Number of v control points
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control points

        task = \'interpolate\': Create an instance of the spline class
        by using an nterpolating spline to given data points. **kwarg
        MUST contain the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: Array of u values 
            v, real, array: Array of v values
            X, real, array, size(len(u),len(v),nDim): Array of data to fit
'''
        #print 'pySpline Class Initialization Type: %s'%(task)

        self.pyspline_real = pyspline_real
        self.pyspline_cs   = pyspline_cs
        self.dtype = 'd'
        if 'complex' in kwargs:
            self.pyspline = pyspline_cs
            self.dtype = 'D'
        else:
            self.pyspline = pyspline_real
            self.dtype = 'd'
        # end if

        if task == 'create':
            assert 'k' in kwargs and 't' in kwargs and \
                'coef' in kwargs and 'range' in kwargs, \
                'Error: k,t,coef, and range MUST be defined for task=\'create\''
            
            self.s = None
            self.X = None
            self.N = None
            self.k = kwargs['k'] 
            self.t = kwargs['t']
            self.coef = kwargs['coef']
            self.Nctl = self.coef.shape[0]
            self.orig_data = False
            self.range = kwargs['range']
            self.nDim = self.coef.shape[1]

            return
             
        if task == 'interpolate':
            assert 'k' in kwargs and 'X' in kwargs, \
                'Error: k, and X MUST be defined for task \'interpolate\''

            self.X  = kwargs['X']
            if len(self.X.shape) == 1:
                self.nDim = 1
            else:
                self.nDim  = self.X.shape[1]
            # end if

            self.N = self.X.shape[0]
            self.k = kwargs['k']
            
            if 's' in kwargs:
                self.s = kwargs['s']
                if self.nDim > 1:
                    self._getLength()
                # end if
            else:
                if self.nDim > 1:
                    self._getParameterization()
                else:
                    print 'Erorr: For 1D mapping splines, the basis, s,\
 must be given as input'
                    sys.exit(1)
                # end if
            # end if
            self.range = array([self.s[0],self.s[-1]])
            self.orig_data = True
            
            if 'dx1' and 'dx2' in kwargs:
                # We have defined a tangent vector at each end

                #print 'dx1 and dx2 defined'
                
                self.Nctl = self.N + 2
                dx1 = kwargs['dx1']
                dx2 = kwargs['dx2']
                ibcl = 1
                ibcr = 1
                kntopt = 1
                assert len(dx1)==len(dx2)==self.nDim,'The length of the \
derivative vectors must match the spatial dimension of the curve'

                if self.nDim > 1:
                    self.coef = zeros((self.Nctl,self.nDim),self.dtype)

                    for idim in xrange(self.nDim):
                        fbcl = dx1[idim]
                        fbcr = dx2[idim]
#                         print 'input: s:',self.s
#                         print 'input: X:',self.X[:,idim]
#                         print 'fbcl:',fbcl
#                         print 'fbcr:',fbcr
                       
                        #self.t,self.coef[:,idim],n,k = \
                        #    self.pyspline.bint4(self.s,self.X[:,idim],\
                                   #                 ibcl,ibcr,fbcl,fbcr,kntopt)

                        self.t,self.coef[:,idim],k = \
                            self.pyspline.bint4(self.s,self.X[:,idim],\
                                                    ibcl,ibcr,fbcl,fbcr,kntopt)

                else:
                    print 'do something'
                # end if


                
        # t,bcoef,k = bint4(x,y,ibcl,ibcr,fbcl,fbcr,kntopt,[ndata,n,w])
        
            else: # Do a regular fit
                
                self.Nctl = self.N


                # Sanity check to make sure k is less than N
                if self.N <= self.k:
                    self.k = self.N
                # end if

                # Generate the knot vector
                self.t = self.pyspline.bknot(self.s,self.k)

                if self.nDim > 1:
                    self.coef = zeros((self.Nctl,self.nDim),self.dtype)
                    for idim in xrange(self.nDim):
                        self.coef[:,idim]= \
                            self.pyspline.bintk(self.s,self.X[:,idim],self.t,self.k)
                else:
                    self.coef = self.pyspline.bintk(self.s,self.X,self.t,self.k)
                # end if
                
            #end for
                
        if task == 'lms':
            assert 'k' in kwargs and 'X' in kwargs and 'Nctl' in kwargs , \
                'Error: k, X and Nctl MUST be defined for task \'interpolate\''
           
            self.X  = kwargs['X']
            if len(self.X.shape) == 1:
                self.nDim = 1
            else:
                self.nDim  = self.X.shape[1]
            # end if

            self.N = self.X.shape[0]
            self.k = kwargs['k']
            self.Nctl = kwargs['Nctl']

            if 's' in kwargs:
                self.s = kwargs['s']
                if self.nDim > 1:
                    self._getLength()
                # end if
            else:
                if self.nDim > 1:
                    self._getParameterization()
                else:
                    print 'Erorr: For 1D mapping splines, the basis, s,\
 must be given as input'
                    sys.exit(1)
                # end if
            # end if
                
            self.range = array([self.s[0],self.s[-1]])
            self.orig_data = True

            # Sanity check to make sure k is less than N
            if self.N <= self.k:
                self.k = self.N
            # end if

            if self.Nctl > self.N:
                self.Nctl  = self.N

            # Generate the knot vector
            self.t = self.pyspline_real.knots(self.s,self.Nctl,self.k)
            if self.nDim > 1:
                # Calculate the Jacobian
                self._calcJacobian()
                self.coef = zeros((self.Nctl,self.nDim),self.dtype)
                for idim in xrange(self.nDim):
                    self.coef[:,idim] = lstsq(self.J,self.X[:,idim])[0]
                # end for
            else:
                self.coef = lstsq(self.J,self.X)[0]
            # end if

        return


    def projectPoint(self,x0,s0=0,Niter=20,tol=1e-8):
        '''Project a point x0 onto the curve and return parametric position
        giving minimum distance. This should also work if point is
        already on the curve as well'''

        # We will use a starting point s0 if given

        # Check we have the same dimension:
        assert len(x0) == self.nDim,'Dimension of x0 and the dimension of\
 spline must be the same'
        converged = False

        for i in xrange(Niter):
            #print 's0:',s0
            D = x0 - self.getValue(s0)
            Ydot = self.getDerivative(s0)
            update = dot(D,Ydot)/(sqrt(dot(Ydot,Ydot)))/self.length

            # Check to see if update went too far

           #  if s0+update > self.range[1]:
#                 update = self.range[1]-s0
#             elif s0+update< self.range[0]:
#                 # Only put the update up to the end
#                 update = self.range[0]-s0
                                
            D2 = x0-self.getValue(s0+update)            
            if abs(dot(D2,D2)) > abs(dot(D,D)):
                update /= 2

            if abs(update)<tol:
                s0 += update
                converged=True
                D = x0-self.getValue(s0) # Get the final Difference
                break
            else:
                s0 += update
                # end if
            # end if
        # end for
        
        return s0,D,converged,update


    def _getParameterization(self):
        # We need to parameterize the curve
        self.s = zeros(self.N,self.dtype);
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

    def _getLength(self):
        # We need to the length of the curve
        s = zeros(self.N,self.dtype);
        for i in xrange(self.N-1):
            dist = 0
            for idim in xrange(self.nDim):
                dist += (self.X[i+1,idim] - self.X[i,idim])**2
                # end for
                s[i+1] = s[i] + sqrt(dist)
            # end for
        # end for
        self.length = s[-1]
        return

    def __call__(self,s):

        if self.nDim == 1:
            x = self.pyspline.bvalu(self.t,self.coef,self.k,0,s)
        else:
            x = zeros((self.nDim),self.dtype)
            for idim in xrange(self.nDim):
                x[idim] = self.pyspline.bvalu(\
                    self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

    def minDistance(self,curve,s=0,t=0,Niter=25,tol=1e-6):

        '''Find the minimum distance between this curve (self) and a second
        curve passed in (curve)'''
        return self.pyspline_real.mincurvedistance(self.t,self.k,self.coef,
                                         curve.t,curve.k,curve.coef,
                                         s,t,Niter,tol)

   
    def getValue(self,s,x=None):
        
        '''Get the value of the spline at point u,v'''
        if self.nDim == 1:
            x = self.pyspline.bvalu(self.t,self.coef,self.k,0,s)
        else:
            if ( x == None ):
                x = zeros((self.nDim),self.dtype)
            # end
            
            for idim in xrange(self.nDim):
                x[idim] = self.pyspline.bvalu(\
                    self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

    def getValueV(self,s):
        '''Get the value of a spline at vector of points s'''
        if self.nDim == 1:
            x = self.pyspline.bvaluv(self.t,self.coef,self.k,0,s)
        else:
            x = zeros((len(s),self.nDim),self.dtype)
            for idim in xrange(self.nDim):
                x[:,idim] = self.pyspline.bvaluv(\
                    self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

 
    def getDerivative(self,s):
        
        '''Get the value of the derivative of spline at point u,v'''
        x = zeros(self.nDim,self.dtype)
        for idim in xrange(self.nDim):
            x[idim] = self.pyspline.bvalu(self.t,self.coef[:,idim],self.k,1,s)
            
        return x        


    def _calcJacobian(self):
        
        # Calculate the jacobian J, for fixed t and s
        self.J = zeros((self.N,self.Nctl),self.dtype)
        ctl = zeros((self.Nctl),self.dtype)
        
        for i in xrange(self.Nctl):
            ctl[i] += 1
            self.J[:,i] = self.pyspline.bvaluv(self.t,ctl,self.k,0,self.s)
            ctl[i] -= 1
        # end for
            
        return


    def writeTecplot(self,handle):
        '''Output this line\'s data to a open file handle \'handle\''''

        if self.orig_data:
            handle.write('Zone T=%s I=%d \n'%('orig_data',self.N))
            handle.write('DATAPACKING=POINT\n')
            for i in xrange(self.N):
                for idim in xrange(self.nDim):
                    handle.write('%f '%(real(self.X[i,idim])))
                # end for
                handle.write('\n')
        # end if

        if self.orig_data:
            s_plot = self.s
        else:
            s_plot = linspace(self.range[0],self.range[1],25).astype('d')
        # end if 
        
        # Dump re-interpolated spline
        handle.write('Zone T=%s I=%d \n'%('interpolated',len(s_plot)))
        handle.write('DATAPACKING=POINT\n')
        for i in xrange(len(s_plot)):
            for idim in xrange(self.nDim):
                handle.write('%f '%(\
                        pyspline_real.bvalu(self.t.astype('d'),\
                                           self.coef[:,idim].astype('d'),\
                                           self.k,0,s_plot[i])))
            # end for 
            handle.write('\n')
        # end for 

        # Dump Control Points (Always have these :-) ) 
        handle.write('Zone T=%s I = %d\n'%('control_pts',self.Nctl))
        handle.write('DATAPACKING=POINT\n')

        for i in xrange(self.Nctl):
            for idim in xrange(self.nDim):
                handle.write('%f '%(real(self.coef[i,idim])))
            # end for
            handle.write('\n')
        # end for

        return


#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pyspline...\n'
    print 'There is an example in the ./example directory'


       

#     def __objcon(self,x):
#         '''Get the rms error for the given set of design variables'''
#         # Unpack the x-values
#         Bcon  = self.Bcon
#         ctl = self.__unpack_x(x)

#         total = 0.0

#         for idim in xrange(self.nDim):
#             total += sum((dot(self.J,ctl[:,:,idim].flatten()) - self.X[:,:,idim].flatten())**2)
#         # end for 
#         fcon = dot(Bcon,x)
#         index = 4*self.nDim + 2*self.Nctlv*self.nDim

#        #  # Calculate the LE constraint
#         for j in xrange(self.Nctlv):

#             A = ctl[0,0,j,:] # Root LE (upper)
#             B = ctl[0,1,j,:] # Root LE upper +1
#             C = ctl[1,-2,j,:]# Root LE lower +1

#             # Area = 0.5*abs( xA*yC - xAyB + xByA - xByC + xCyB - xCyA )

#             A1 = A[0]*C[1] - A[0]*B[1] + B[0]*A[1] -B[0]*C[1] + C[0]*B[1] - C[0]*A[1]
#             A2 = A[1]*C[2] - A[1]*B[2] + B[1]*A[2] -B[1]*C[2] + C[1]*B[2] - C[1]*A[2]
#             A3 = A[0]*C[2] - A[0]*B[2] + B[0]*A[2] -B[0]*C[2] + C[0]*B[2] - C[0]*A[2]

#             fcon[index:index+3] = array([A1,A2,A3])
#             index += 3
#         # end for

#         return total,fcon,False


#     def __sens(self,x,f_obj,f_con):

#         ndv = len(x)
#         g_obj = zeros(ndv)
#          # Unpack the x-values
#         ctl = self.__unpack_x(x)
#         N = self.Nctlu*self.Nctlv

#         for idim in xrange(self.nDim):
#             g_obj[self.nDim*N + idim*N :  self.nDim*N + idim*N + N] = \
#                 2*dot(dot(self.J,ctl[:,:,idim].flatten())-self.X[:,:,idim].flatten(),self.J)
#         # end for 

#         g_con = self.Bcon
#         h = 1.0e-40j
#         x = array(x,'D')

#         for i in xrange(ndv):
#             index = 4*self.nDim + 2*self.Nctlv*self.nDim
#             x[i] += h
#             ctl = self.__unpack_x(x,'D')
#             for j in xrange(self.Nctlv):
#                 A = ctl[0,0,j,:] # Root LE (upper)
#                 B = ctl[0,1,j,:] # Root LE upper +1
#                 C = ctl[1,-2,j,:]# Root LE lower +1

#                 # Area = 0.5*abs( xA*yC - xAyB + xByA - xByC + xCyB - xCyA )

#                 A1 = A[0]*C[1] - A[0]*B[1] + B[0]*A[1] -B[0]*C[1] + C[0]*B[1] - C[0]*A[1]
#                 A2 = A[1]*C[2] - A[1]*B[2] + B[1]*A[2] -B[1]*C[2] + C[1]*B[2] - C[1]*A[2]
#                 A3 = A[0]*C[2] - A[0]*B[2] + B[0]*A[2] -B[0]*C[2] + C[0]*B[2] - C[0]*A[2]

#                 g_con[index:index+3,i] = imag(array([A1,A2,A3]))/imag(h)
#                 index += 3
#             # end for
#             x[i] -= h
#         # end for
#         return g_obj,g_con,False

#     def __unpack_x(self,x,dtype='d'):
#         ctl = zeros((self.Nctlu,self.Nctlv,self.nDim),dtype)
#         N = self.Nctlu*self.Nctlv
#         for idim in xrange(self.nDim):
#             ctl[:,:,idim] = reshape(x[self.nDim*N + idim*N : self.nDim*N + idim*N + N],[self.Nctlu,self.Nctlv])
#         # end for
#         return ctl


