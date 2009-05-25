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
        positions u and v. We will find an interpolating spline with 
        knot vector values placed according to the spacing of u and v.
        Spline order can be set in both u and v are are cubic by default.'''


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

    def findUV(self,x0,z0,u0,v0):
        ''' Try to find the u-v coordinate of the curved spline coorsponding
        to the point x0,z0 as given'''

        #TODO: Add in error checking. u,v must be bound in [-1,1].
        # if step is too large, make take smaller steps
        
        maxIter = 15

        u = u0
        v = v0
        
        for iter in xrange(maxIter):
            x = self.getValue(u,v) #x contains the x,y,z coordinates 
            f = mat([x[0]-x0,x[2]-z0])
            J = self.getJacobian(u,v)
            A = mat(vstack((J[0,:],J[2,:])))
            x_up = scipy.linalg.solve(A,-f.T)

            u = u + x_up[0]
            v = v + x_up[1]
        # end for

        return u,v

#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pySpline...\n'




