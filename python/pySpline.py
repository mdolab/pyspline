#!/usr/local/bin/python
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
import os, sys, string

# =============================================================================
# External Python modules
# =============================================================================
from numpy import array, mat, vstack, zeros, ones, linspace, meshgrid, floor, \
    sqrt, cos, sin, pi

import scipy.linalg


# =============================================================================
# Extension modules
# =============================================================================
import pyspline 
import pyspline_cs

# =============================================================================
# pySpline class
# =============================================================================
class spline_interp():
	
    '''
    Spline Object Class (Interpolating Spline)
    '''
    
    def __init__(self, u,v,x,y,z,ku=4,kv=4,*args, **kwargs):
        
        '''Initialize the class with the data points x,y,z defined at 
        parametric positions u and v. We will find an interpolating spline with 
        knot vector values placed according to the spacing of u and v.
        Spline order can be set in both u and v are are cubic (k=4) by 
        default.'''

        self.tu_x,self.tv_x,self.bcoef_x = pyspline.b2ink(u,v,x,ku,kv)
        self.tu_y,self.tv_y,self.bcoef_y = pyspline.b2ink(u,v,y,ku,kv)
        self.tu_z,self.tv_z,self.bcoef_z = pyspline.b2ink(u,v,z,ku,kv)
        self.ku = ku
        self.kv = kv
        self.x0 = x
        self.y0 = y
        self.z0 = z
        self.u0 = u
        self.v0 = v
        
        return
	

    def getValue(self,u,v):
        
        '''Get the value of the spline at point u,v'''
        
        x = pyspline.b2val(u,v,0,0,self.tu_x,self.tv_x,self.ku,self.kv,self.bcoef_x)
        y = pyspline.b2val(u,v,0,0,self.tu_y,self.tv_y,self.ku,self.kv,self.bcoef_y)
        z = pyspline.b2val(u,v,0,0,self.tu_z,self.tv_z,self.ku,self.kv,self.bcoef_z)
        
        return array([x,y,z],float)


    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((3,2))
        J[0,0] = pyspline.b2val(u,v,1,0,self.tu_x,self.tv_x,self.ku,self.kv,self.bcoef_x)
        J[0,1] = pyspline.b2val(u,v,0,1,self.tu_x,self.tv_x,self.ku,self.kv,self.bcoef_x)

        J[1,0] = pyspline.b2val(u,v,1,0,self.tu_y,self.tv_y,self.ku,self.kv,self.bcoef_y)
        J[1,1] = pyspline.b2val(u,v,0,1,self.tu_y,self.tv_y,self.ku,self.kv,self.bcoef_y)

        J[2,0] = pyspline.b2val(u,v,1,0,self.tu_z,self.tv_z,self.ku,self.kv,self.bcoef_z)
        J[2,1] = pyspline.b2val(u,v,0,1,self.tu_z,self.tv_z,self.ku,self.kv,self.bcoef_z)

        return J

    def findUV(self,x0,r,u0,v0):
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

        maxIter = 50

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

            f = mat(zeros((3,1)))
            f[0] = x[0]-(x0[0]+r[0]*s)
            f[1] = x[1]-(x0[1]+r[1]*s)
            f[2] = x[2]-(x0[2]+r[2]*s)

            J = self.getJacobian(u,v)
            A = mat(zeros((3,3)))
            A[:,0:2] = J
            A[0,2]   = -r[0]
            A[1,2]   = -r[1]
            A[2,2]   = -r[2]

            x_up = scipy.linalg.solve(A,-f)

            # Add a little error checking here:
            
            if u + x_up[0] < -1 or u + x_up[0] > 1 or \
               v + x_up[1] < -1 or v + x_up[1] > 1:
                #Cut the size of the step in 1/2
                x_up /= 2

            u = u + x_up[0]
            v = v + x_up[1]
            s = s + x_up[2]

            if scipy.linalg.norm(x_up) < 1e-12:
                return u,v,x

        # end for

        print 'Warning: Newton Iteration for u,v,s did not converge:'
        print 'u = %f, v = %f, s = %f\n'%(u,v,s)
        print 'Norm of update:',scipy.linalg.norm(x_up)

        return u,v,x


class spline_lms():
	
    '''
    Spline Object Class (Least Mean Squares)
    '''
    
    def __init__(self, u,v,x,y,z,nctlu=None,nctlv=10,w=None,ku=4,kv=4,tv_type='cosine',*args, **kwargs):        
        '''Initialize the class with the data points x,y,z defined at 
        parametric positions u and v. We will find an interpolating spline with 
        knot vector values placed according to the spacing of u and v.
        Spline order can be set in both u and v are are cubic (k=4) by 
        default.'''

        # NOTE: fitpack uses the  polynomial order, NOT the spline order.
        # So a cubic spline is 3 and NOT 4. However, This will be adjust in
        # this class such that the user sees a consistent interface

        if not(w):
            w  = ones(x.shape)
        # end if

        # Lets do something with the weighting...we want to match the le more than elsewhere
       # w[:,0] = 100
       # w[:,-1] = 100
#        w[:,1] = 1000
#         w[:,2] = 100
#         w[:,3] = 50
#         w[:,4] = 10
#         w[:,5] = 5

        #print 'w',w


        [V,U] = meshgrid(v,u)
        
        # nctl is the number of control points in the chord-wise direction
        if not(nctlu):
            nctlu = len(u) # Use the number of u-points as number of ctrl pts
        # end if
        
        # Options
        iopt = -1 #least squares surface
        ku -= 1
        kv -= 1
        
        nu = nctlu + ku + 1 # number of knots in u
        nv = nctlv + kv + 1 # number of knots in v
        smooth = 5e-5 # Only required for iopt = 0

        # Knot Stuff

        m = len(x.flatten()) # number of data points

        # Knot estimates
        nuest=nu # estimate of max number of knots (u)
        nvest=nv # estimate of max number of knots (v)

        # U-Knots
        t = pyspline.bknot(u,ku+1) #cmlib k sense
        tu = zeros(max(nuest,nvest))
        tu[0:nu] = t
    
        # V-Knots
        tv = zeros(max(nuest,nvest))
        tv[0:kv+1] = 0
        tv[nv-kv-1:nv] = 1

        if tv_type =='cosine': # do a cosine knot vector
            tv[kv:nv-kv] = 0.5*(1-cos(linspace(0,pi,nv-2*kv)))
            #tv[kv:nv-kv] = linspace(0,1,nv-2*kv)**3
        elif tv_type == 'linear':
            tv[kv:nv-kv] = linspace(0,1,nv-2*kv)
        else:
            print 'Unknown type of knot vector...using linear distribution'
            tv[kv:nv-kv] = linspace(0,1,nv-2*kv)
        # end if

        # Not do the LMS Fits
            
 #        print 'tu:',tu
#         print 'tv:',tv
#         print 'nu:',nu
#         print 'nv:',nv

#         print 'nctlu:',nctlu
#         print 'nctlv:',nctlv
#         print 'nuest:',nuest
#         print 'nvest:',nvest
#         print 'u:',u
#         print 'v:',v
#         print 'U:',U
#         print 'x:',z

#        sys.exit(0)

        nu,tu,nv,tv,cx,fpx,ierx = \
            pyspline.surfit(iopt,U.flatten(),V.flatten(),x.flatten(),w.flatten()\
                                ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
        nu,tu,nv,tv,cy,fpy,iery = \
            pyspline.surfit(iopt,U.flatten(),V.flatten(),y.flatten(),w.flatten()\
                                ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
        nu,tu,nv,tv,cz,fpz,ierz = \
            pyspline.surfit(iopt,U.flatten(),V.flatten(),z.flatten(),w.flatten()\
                                ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
        print 'Fitting Error: x,y,z:',fpx,fpy,fpz

        if ierx!=0 or iery !=0 or ierz != 0:
            print 'There was an error in fitting one of the surface dimensions.'
            print 'The errors codes are:'
            print 'X-cordinate error:',ierx
            print 'Y-cordinate error:',iery
            print 'Z-cordinate error:',ierz
            print 'See surfit.f in src/fitpack for an explination of error messages'
        #end if

        self.nctlu = nctlu
        self.nctlv = nctlv
        self.tu = tu[0:nu]
        self.tv = tv[0:nv]
        self.cx = cx[0:(nu-kv-1)*(nv-kv-1)]
        self.cy = cy[0:(nu-kv-1)*(nv-kv-1)]
        self.cz = cz[0:(nu-kv-1)*(nv-kv-1)]

        self.ku = ku
        self.kv = kv
        self.x0 = x
        self.y0 = y
        self.z0 = z
        self.u0 = u
        self.v0 = v

        return
	

    def getValue(self,u,v):
        
        '''Get the value of the spline at point u,v'''
        x,ierx = pyspline.bispev(self.tu,self.tv, self.cx, self.ku,self.kv,u,v)
        y,iery = pyspline.bispev(self.tu,self.tv, self.cy, self.ku,self.kv,u,v)
        z,ierz = pyspline.bispev(self.tu,self.tv, self.cz, self.ku,self.kv,u,v)

        if ierx!=0 or iery!=0 or ierz != 0:
            print 'There was an error in spline_lms.getValue' 
            print 'The error codes are:'
            print 'X-cordinate error:',ierx
            print 'Y-cordinate error:',iery
            print 'Z-cordinate error:',ierz
            print 'See bispev.f in src/fitpack for an explination of error messages'
        #end if

        return array([x[0],y[0],z[0]],float)


    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((3,2))
        J[0,0] = pyspline.b2val(u,v,1,0,self.tu_x,self.tv_x,self.ku,self.kv,self.bcoef_x)
        J[0,1] = pyspline.b2val(u,v,0,1,self.tu_x,self.tv_x,self.ku,self.kv,self.bcoef_x)

        J[1,0] = pyspline.b2val(u,v,1,0,self.tu_y,self.tv_y,self.ku,self.kv,self.bcoef_y)
        J[1,1] = pyspline.b2val(u,v,0,1,self.tu_y,self.tv_y,self.ku,self.kv,self.bcoef_y)

        J[2,0] = pyspline.b2val(u,v,1,0,self.tu_z,self.tv_z,self.ku,self.kv,self.bcoef_z)
        J[2,1] = pyspline.b2val(u,v,0,1,self.tu_z,self.tv_z,self.ku,self.kv,self.bcoef_z)

        return J

    def findUV(self,x0,r,u0,v0):
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

        maxIter = 50

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

            f = mat(zeros((3,1)))
            f[0] = x[0]-(x0[0]+r[0]*s)
            f[1] = x[1]-(x0[1]+r[1]*s)
            f[2] = x[2]-(x0[2]+r[2]*s)

            J = self.getJacobian(u,v)
            A = mat(zeros((3,3)))
            A[:,0:2] = J
            A[0,2]   = -r[0]
            A[1,2]   = -r[1]
            A[2,2]   = -r[2]

            x_up = scipy.linalg.solve(A,-f)

            # Add a little error checking here:
            
            if u + x_up[0] < -1 or u + x_up[0] > 1 or \
               v + x_up[1] < -1 or v + x_up[1] > 1:
                #Cut the size of the step in 1/2
                x_up /= 2

            u = u + x_up[0]
            v = v + x_up[1]
            s = s + x_up[2]

            if scipy.linalg.norm(x_up) < 1e-12:
                return u,v,x

        # end for

        print 'Warning: Newton Iteration for u,v,s did not converge:'
        print 'u = %f, v = %f, s = %f\n'%(u,v,s)
        print 'Norm of update:',scipy.linalg.norm(x_up)

        return u,v,x




#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pySpline...\n'
    print 'There is an example in the ./example directory'




