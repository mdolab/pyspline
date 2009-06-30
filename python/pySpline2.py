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
import os, sys, string, time, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot, cross, mod

import numpy.linalg

#from pyOpt_optimization import Optimization
#from pySNOPT import SNOPT
import pyspline
import pyspline_cs

# =============================================================================
# pySpline class
# =============================================================================

class surf_spline():

    def __init__(self,task='create',*args,**kwargs):

        '''Create an instance of a b-spline surface. There are three ways to initialize 
        the class as determined by the task flag:

        task = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for u
            kv, integer: Order for v
            Nctlu, integer: Number of u control points
            Nctlv, integer: Number of v control points
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control point values

        task = \'interpolate\': Create an instance of the spline class
        by using an interpolating spline to given data points. **kwarg
        MUST contain the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: list of u values 
            v, real, array: list of v values
            X, real, array, size(len(u),len(v),nDim): Array of data points to fit

        task = \'lms\': Create an instance of the spline class using a
        Least-Mean-Squares fit to the spline. . **kwargs MUST contain
        the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: list of u values 
            v, real, array: list of v values
            X, real, array, size(len(u),len(v),nDim): Array of data points to fit
      
            
'''
        print 'pySpline Class Initialization Type: %s'%(task)

        if task == 'create':
            assert 'ku' in kwargs and 'kv' in kwargs  and 'tu' in kwargs \
                and 'tv' in kwargs and 'coef' in kwargs and 'range' in kwargs, \
                'Error: ku,kv,tu,tv,coef and range MUST be defined for task=\'create\''
            
            self.u = None
            self.v = None
            self.X = None
            self.Nu = None
            self.Nv = None
            self.ku = kwargs['ku'] 
            self.kv = kwargs['kv']
            self.tu = kwargs['tu']
            self.tv = kwargs['tv']
            self.coef = kwargs['coef']
            self.Nctlu = self.coef.shape[0]
            self.Nctlv = self.coef.shape[1]
            self.orig_data = False
            self.range = kwargs['range']
            self.nDim = self.coef.shape[2]
            return
     
        if task == 'interpolate':
            
            assert 'ku' in kwargs and 'kv' and 'X' in kwargs,\
                'Error: ku,kv,u,v and X MUST be defined for task \'interpolate\''

            self.X  = kwargs['X']
            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.nDim  = self.X.shape[2]
            self.Nctlu = self.Nu
            self.Nctlv = self.Nv
            self.ku = kwargs['ku']
            self.kv = kwargs['kv']


            if 'u' in kwargs and 'v' in kwargs:
                self.u = kwargs['u']
                self.v = kwargs['v']
            else:
                self._calcParameterization()
            # end if
                
            self.orig_data = True
            
            self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim))
            self.range = array([self.u[0],self.u[-1],self.v[0],self.v[-1]])


            #Sanity check to make sure k is less than N
            if self.Nu <= self.ku:
                self.ku = self.Nu-1
            if self.Nv <= self.kv:
                self.kv = self.Nv-1

            for idim in xrange(self.nDim):
                self.tu,self.tv,self.coef[:,:,idim]= pyspline.b2ink(self.u,self.v,self.X[:,:,idim],self.ku,self.kv)
            #end for
                
            return
           
        elif task == 'lms':
            # Do some checking on the number of control points

            assert 'ku' in kwargs and 'kv' in kwargs and \
                   'Nctlu' in kwargs and 'Nctlv' in kwargs and 'X' in kwargs, \
                   'Error: ku,kv,Nctlu, Nctlv and X MUST be defined for task \'lms\''

            self.X  = kwargs['X']
            self.Nu = self.X.shape[0]
            self.Nv = self.X.shape[1]
            self.nDim  = self.X.shape[2]
            self.Nctlu = kwargs['Nctlu']
            self.Nctlv = kwargs['Nctlv']
            self.ku = kwargs['ku']
            self.kv = kwargs['kv']

            if 'u' in kwargs and 'v' in kwargs:
                self.u = kwargs['u']
                self.v = kwargs['v']
            else:
                self._calcParameterization()
            # end if
                
            self.orig_data = True
            self.range = array([self.u[0],self.u[-1],self.v[0],self.v[-1]])

            # Sanity Check on Inputs
            if self.Nctlu > self.Nu:
                self.Nctlu  = self.Nu
            if self.Nctlv > self.Nv:
                self.Nctlv = self.Nv

            if self.Nu <= self.ku:
                self.ku = self.Nu-1
            if self.Nv <= self.kv:
                self.kv = self.Nv-1

           #Calculate the knot vector and Jacobian
            self._calcKnots()
            self._calcJacobian()

#             # Lets do a lms 
            timeA = time.time()
            self.coef = pyspline.fit_surf(self.Nu,self.Nv,self.Nctlu,self.Nctlv,self.J,self.X)
            print 'LMS Fit Time:',time.time()-timeA
            
            return 

        print 'Error: task is not understood. Type must be \'lms\',\'interpolate\' or \'create\''
        sys.exit(0)
        return

    def _calcParameterization(self):

        u = (zeros(self.Nu))
        singular_counter = 0

        for j in xrange(self.Nv): #loop over each v, and average the 'u' parameter 
            temp = zeros(self.Nu,'d')

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
        
        v = zeros(self.Nv)
        singular_counter = 0
        for i in xrange(self.Nu): #loop over each v, and average the 'u' parameter 
            temp = zeros(self.Nv)
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

        self.tu = pyspline.knots(self.u,self.Nctlu,self.ku)
        self.tv = pyspline.knots(self.v,self.Nctlv,self.kv)
        
        return

    def _calcJacobian(self):
        
        # Calculate the jacobian J, for fixed t and s
        h = 1.0e-40j
        self.J = zeros([self.Nu*self.Nv,self.Nctlu*self.Nctlv])
        ctl = zeros([self.Nctlu,self.Nctlv],'D')

        [V, U] = meshgrid(self.v,self.u)
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                ctl[i,j] += h
                val = pyspline_cs.b2valv(U.flatten(),V.flatten(),0,0,self.tu,self.tv,self.ku,self.kv,ctl)
                ctl[i,j] -= h    
                self.J[:,i*self.Nctlv + j] = imag(val)/imag(h)
                # end for
            # end for 
        # end for
        return

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

    def getValueEdge(self,edge,s):
        '''Get the value of the spline on edge, edge=0,1,2,3 where
        edges are oriented in the standard counter-clockwise fashion.'''

        if edge == 0:
            return self.getValue(s,0)
        elif edge == 1:
            return self.getValue(1,s)
        elif edge == 2:
            return self.getValue(1-s,1)
        elif edge ==3:
            return self.getValue(0,1-s)
        else:
            print 'Edge must be between 0 and 3'
            sys.exit(1)
            return
        #end if 

    def getValue(self,u,v):
        
        '''Get the value of the spline at point u,v'''
        x = zeros([self.nDim])
        for idim in xrange(self.nDim):
            x[idim] = pyspline.b2val(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        
        return x

    def getValueV(self,u,v):
        '''Get the value of a spline at vector of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same length'
        x = zeros(len(u),self.nDim)
        for idim in xrange(self.nDim):
            x[:,idim] = pyspline.b2valv(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        return x

    def getValueM(self,u,v):
        '''Get the value of a spline at matrix of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same shape'
        x = zeros((u.shape[0],u.shape[1],self.nDim))
        for idim in xrange(self.nDim):
            x[:,:,idim] = pyspline.b2valm(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return x

    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((self.nDim,2))
        
        for idim in xrange(self.nDim):
            J[idim,0] = pyspline.b2val(u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            J[idim,1] = pyspline.b2val(u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

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

    def writeTecplot(self,handle):
        '''Output this surface\'s data to a open file handle \'handle\''''

        if self.orig_data:
            handle.write('Zone T=%s I=%d J = %d\n'%('orig_data',self.Nu,self.Nv))
            handle.write('DATAPACKING=POINT\n')
            for j in xrange(self.Nv):
                for i in xrange(self.Nu):
                    handle.write('%f %f %f \n'%(self.X[i,j,0],self.X[i,j,1],self.X[i,j,2]))
                # end for
            # end for 
        # end if
                
        if self.orig_data:
            u_plot = self.u
            v_plot = self.v
        else:
            u_plot = linspace(self.range[0],self.range[1],25)
            v_plot = linspace(self.range[2],self.range[3],25)
        # end if 
      

        # Dump re-interpolated surface
        handle.write('Zone T=%s I=%d J = %d\n'%('interpolated',len(u_plot),len(v_plot)))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(len(v_plot)):
            for i in xrange(len(u_plot)):
                for idim in xrange(self.nDim):
                    handle.write('%f '%(pyspline.b2val(u_plot[i],v_plot[j],0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])))
                # end for 
                handle.write('\n')
            # end for
        # end for 

        # Dump Control Points (Always have these :-) ) 
        handle.write('Zone T=%s I=%d J = %d\n'%('control_pts',self.Nctlu,self.Nctlv))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                handle.write('%f %f %f \n'%(self.coef[i,j,0],self.coef[i,j,1],self.coef[i,j,2]))
            # end for
        # end for 


    def writeIGES_directory(self,handle,Dcount,Pcount):

        '''Write the IGES file header information (Directory Entry Section) for this surface'''
        # A simplier Calc based on cmlib definations
        # The 13 is for the 9 parameters at the start, and 4 at the end. See the IGES 5.3 Manual
        #paraEntries = 13 + Knotsu              + Knotsv               + Weights               + control points
        assert self.nDim == 3, 'Must have 3 dimensions to write to IGES file'
        paraEntries = 13 + (len(self.tu)) + (len(self.tv)) +   self.Nctlu*self.Nctlv + 3*self.Nctlu*self.Nctlv+1

        paraLines = paraEntries / 5
        if mod(paraEntries,5) != 0: paraLines += 1

        handle.write('     128%8d       0       0       1       0       0       000000001D%7d\n'%(Pcount,Dcount))
        handle.write('     128       0       2%8d       0                               0D%7d\n'%(paraLines,Dcount+1))
        Dcount += 2
        Pcount += paraLines
        return Pcount , Dcount


    def writeIGES_parameters(self,handle,Pcount,counter):
        '''Write the IGES parameter information for this surface'''

        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'%(128,self.Nctlu-1,self.Nctlv-1,self.ku-1,self.kv-1,Pcount,counter))
        counter += 1
        handle.write('%12d,%12d,%12d,%12d,%12d,%7dP%7d\n'%(0,0,1,0,0,Pcount,counter))
        counter += 1
        
        pos_counter = 0


        for i in xrange(len(self.tu)):
            pos_counter += 1
            handle.write('%12.6g,'%(self.tu[i]))
            if mod(pos_counter,5) == 0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0


        for i in xrange(len(self.tv)):
            pos_counter += 1
            handle.write('%12.6g,'%(self.tv[i]))
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
                    handle.write('%12.6g,'%(self.coef[i,j,idim]))
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
                handle.write('%12.6g,'%(self.tu[0]))
            if i == 1:
                handle.write('%12.6g,'%(self.tu[-1]))
            if i == 2:
                handle.write('%12.6g,'%(self.tv[0]))
            if i == 3:
                handle.write('%12.6g;'%(self.tv[-1])) # semi-colon for the last entity
            if mod(pos_counter,5)==0:
                handle.write('%7dP%7d\n'%(Pcount,counter))
                counter += 1
                pos_counter = 0
            else: # We have to close it up anyway
                if i ==3:
                    for j  in xrange(5-pos_counter):
                        handle.write('%13s'%(' '))
                        
                    pos_counter = 0
                    handle.write('%7dP%7d\n'%(Pcount,counter))
                    counter += 1
        Pcount += 2
        return Pcount,counter


class linear_spline():

    def __init__(self,task='create',*args,**kwargs):

        '''Create an instance of a b-spline surface. There are three ways to initialize 
        the class as determined by the task flag:

        task = \'create\': Create an instance of the spline class
        directly by supplying the required information. **kwargs MUST
        contain the folloiwng information:

            ku, integer: Order for u
            kv, integer: Order for v
            Nctlu, integer: Number of u control points
            Nctlv, integer: Number of v control points
            tu, real array: Knot vector for u
            tv, real array: Knot vector for v
            coef, real array size(Nctlu,Nctlv,nDim): Array of control point values

        task = \'interpolate\': Create an instance of the spline class
        by using an nterpolating spline to given data points. **kwarg
        MUST contain the following information:

            ku, integer: Order for u
            kv, integer: Order for v
            u, real, array: Array of u values 
            v, real, array: Array of v values
            X, real, array, size(len(u),len(v),nDim): Array of data points to fit
'''
        print 'pySpline Class Initialization Type: %s'%(task)

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
            self.nDim  = self.X.shape[1]
            self.N = self.X.shape[0]
            self.k = kwargs['k']
            
            if 's' in kwargs:
                self.s = kwargs['s']
            else:
                self._getParameterization()
            # end if
                
            self.Nctl = self.N
            self.coef = zeros((self.Nctl,self.nDim))
            self.range = array([self.s[0],self.s[-1]])
            self.orig_data = True

            # Sanity check to make sure k is less than N
            if self.N <= self.k:
                self.k = self.N-1
            # end if
    
            # Generate the knot vector
            self.t = pyspline.bknot(self.s,self.k)

            for idim in xrange(self.nDim):
                self.coef[:,idim]= pyspline.bintk(self.s,self.X[:,idim],self.t,self.k)
            #end for
                
            return


    def _getParameterization(self):
        # We need to parameterize the curve
        self.s = zeros(self.N);
        for i in xrange(self.N-1):
            dist = 0
            for idim in xrange(self.nDim):
                dist += (self.X[i+1,idim] - self.X[i,idim])**2
                # end for
                self.s[i+1] = self.s[i] + sqrt(dist)
            # end for
        # end for
        return
   
    def getValue(self,s):
        
        '''Get the value of the spline at point u,v'''
        x = zeros([self.nDim])
        for idim in xrange(self.nDim):
            x[idim] = pyspline.bvalu(self.t,self.coef[:,idim],self.k,0,s)

        return x
        

    def writeTecplot(self,handle):
        '''Output this line\'s data to a open file handle \'handle\''''

        if self.orig_data:
            handle.write('Zone T=%s I=%d \n'%('orig_data',self.N))
            handle.write('DATAPACKING=POINT\n')
            for i in xrange(self.N):
                for idim in xrange(self.nDim):
                    handle.write('%f '%(self.X[i,idim]))
                # end for
                handle.write('\n')
        # end if

        if self.orig_data:
            s_plot = self.s
        else:
            s_plot = linspace(self.range[0],self.range[1],25)
        # end if 

        # Dump re-interpolated spline
        handle.write('Zone T=%s I=%d \n'%('interpolated',len(s_plot)))
        handle.write('DATAPACKING=POINT\n')
        for i in xrange(len(s_plot)):
            for idim in xrange(self.nDim):
                handle.write('%f '%(pyspline.bvalu(self.t,self.coef[:,idim],self.k,0,s_plot[i])))
            # end for 
            handle.write('\n')
        # end for 

        # Dump Control Points (Always have these :-) ) 
        handle.write('Zone T=%s I = %d\n'%('control_pts',self.N))
        handle.write('DATAPACKING=POINT\n')
        for i in xrange(self.N):
            for idim in xrange(self.nDim):
                handle.write('%f '%(self.coef[i,idim]))
            # end for
            handle.write('\n')
        # end for

        return



#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pySpline...\n'
    print 'There is an example in the ./example directory'




