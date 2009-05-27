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
from numpy import array,mat,vstack,zeros

import scipy.linalg


# =============================================================================
# Extension modules
# =============================================================================
import pyspline 
import pyspline_cs

# =============================================================================
# pySpline class
# =============================================================================
class spline():
	
    '''
    Spline Object Class 
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



#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pySpline...\n'




