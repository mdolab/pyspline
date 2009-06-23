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
import os, sys, string, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot, cross, mod

import scipy.linalg
from pyOpt_optimization import Optimization
from pySNOPT import SNOPT
import pyspline
import pyspline_cs

# =============================================================================
# pySpline class
# =============================================================================

class spline():

    def __init__(self,u,v,X,task='create',Nctlu = 30, Nctlv = 30, ku=4,kv=4,*args, **kwargs):
        print 'Task: %s'%(task)

        self.u0 = u
        self.v0 = v        
        self.x0 = X
        self.Nu = len(u)
        self.Nv = len(v)     
        if task == 'interpolate':
            self.Nctlu = len(u)
            self.Nctlv = len(v)
            self.ku = ku
            self.kv = kv
          
            self.coef = zeros((self.Nctlu,self.Nctlv,3))

            if self.Nu <= ku:
                self.ku = self.Nu-1
            if self.Nv <= kv:
                self.kv = self.Nv-1

            for idim in xrange(3):
                self.tu,self.tv,self.coef[:,:,idim]= pyspline.b2ink(u,v,X[:,:,idim],self.ku,self.kv)
            #end for
           
        elif task == 'lms':
            # Do some checking on the number of control points
            print 'Nu,Nv:',self.Nu,self.Nv
            self.Nctlu = Nctlu
            self.Nctlv = Nctlv
            self.ku = ku
            self.kv = kv
            if self.Nctlu > self.Nu:
                self.Nctlu  = self.Nu
            if self.Nctlv > self.Nv:
                self.Nctlv = self.Nv

            if self.Nu <= ku:
                ku = self.Nu-1
                self.ku = ku
            if self.Nv <= kv:
                kv = self.Nv-1
                self.kv = kv

#            print 'Nu,Nv,Nctlu,Nctlv,ku,kv:',self.Nu,self.Nv,self.Nctlu,self.Nctlv,self.ku,self.kv

            #Calculate the knot vector
            self.__calcKnots()
#            print 'knots: tu,tv:',self.tu,self.tv

            self.__calcJacobian()
#            print 'jacobian shape:',self.J.shape

#             # Lets do a lms 
            print 'before lms:'
            timeA = time.time()
            print self.J.shape
            ctl = pyspline.fit_surf(self.Nu,self.Nv,self.Nctlu,self.Nctlv,self.J,X)
            print 'LMS Fit Time:',time.time()-timeA
            self.coef = ctl
        else:
            print 'Error: fit_type is not understood. Type must be \'lms\' or \'interpolate\''
            sys.exit(0)
        return

    def __calcKnots(self):

        '''Find the initial knot vector for this problem'''

        # ==============================
        # Find the Initial Knot Vectors
        # ==============================
        
        # U knots
        tu= zeros(self.Nctlu + self.ku)
        tu[self.ku-1:self.Nctlu+1]  = 0.5*(1-cos(linspace(0,pi,self.Nctlu-self.ku+2)))
        tu[0:self.ku] = 0
        tu[self.Nctlu:self.Nctlu+self.ku] = 1.0#+0.1*(tu0[Nctlu]-tu0[Nctlu-1])
            
        # V Knots

        tv= zeros(self.Nctlv + self.kv)
        tv[self.kv-1:self.Nctlv+1] =  linspace(0,1,self.Nctlv-self.kv+2)
            #tv0[kv-1:Nctlv+1] = v[0,1:-1]
            #          tv0[kv-1:Nctlv+1]  = 0.5*(1-cos(linspace(0,pi,Nctlv-kv+2)))
        tv[0:self.kv] = 0
        tv[self.Nctlv:self.Nctlv+self.kv] = 1.0#+0.1*(tv0[Nctlv]-tv0[Nctlv-1])
        
        self.tu = tu
        self.tv = tv

        return

    def __calcJacobian(self):
        
        # Calculate the jacobian J, for fixed t and s
        h = 1.0e-40j
        self.J = zeros([self.Nu*self.Nv,self.Nctlu*self.Nctlv])
        ctl = zeros([self.Nctlu,self.Nctlv],'D')

        V = zeros((self.Nu,self.Nv))
        U = zeros((self.Nu,self.Nv))
        #print 'Nu,Nv,Nctlu,Nctlv,ku,kv:',self.Nu,self.Nv,self.Nctlu,self.Nctlv,self.ku,self.kv
        #print 'u0,v0:',self.u0,self.v0
        [V, U] = meshgrid(self.v0,self.u0)
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                ctl[i,j] += h
                val = pyspline_cs.b2valv(U.flatten(),V.flatten(),0,0,self.tu,self.tv,self.ku,self.kv,ctl)
                ctl[i,j] -= h    
                self.J[:,i*self.Nctlv + j] = imag(val)/imag(h)
                # end for
            # end for 
        # end for


    def __objcon(self,x):
        '''Get the rms error for the given set of design variables'''
        # Unpack the x-values
        Bcon  = self.Bcon
        ctl = self.__unpack_x(x)

        total = 0.0

        for idim in xrange(3):
            total += sum((dot(self.J,ctl[:,:,idim].flatten()) - self.x0[:,:,idim].flatten())**2)
        # end for 
        fcon = dot(Bcon,x)
        index = 4*3 + 2*self.Nctlv*3

       #  # Calculate the LE constraint
        for j in xrange(self.Nctlv):

            A = ctl[0,0,j,:] # Root LE (upper)
            B = ctl[0,1,j,:] # Root LE upper +1
            C = ctl[1,-2,j,:]# Root LE lower +1

            # Area = 0.5*abs( xA*yC - xAyB + xByA - xByC + xCyB - xCyA )

            A1 = A[0]*C[1] - A[0]*B[1] + B[0]*A[1] -B[0]*C[1] + C[0]*B[1] - C[0]*A[1]
            A2 = A[1]*C[2] - A[1]*B[2] + B[1]*A[2] -B[1]*C[2] + C[1]*B[2] - C[1]*A[2]
            A3 = A[0]*C[2] - A[0]*B[2] + B[0]*A[2] -B[0]*C[2] + C[0]*B[2] - C[0]*A[2]

            fcon[index:index+3] = array([A1,A2,A3])
            index += 3
        # end for

        return total,fcon,False


    def __sens(self,x,f_obj,f_con):

        ndv = len(x)
        g_obj = zeros(ndv)
         # Unpack the x-values
        ctl = self.__unpack_x(x)
        N = self.Nctlu*self.Nctlv

        for idim in xrange(3):
            g_obj[3*N + idim*N :  3*N + idim*N + N] = \
                2*dot(dot(self.J,ctl[:,:,idim].flatten())-self.x0[:,:,idim].flatten(),self.J)
        # end for 

        g_con = self.Bcon
        h = 1.0e-40j
        x = array(x,'D')

        for i in xrange(ndv):
            index = 4*3 + 2*self.Nctlv*3
            x[i] += h
            ctl = self.__unpack_x(x,'D')
            for j in xrange(self.Nctlv):
                A = ctl[0,0,j,:] # Root LE (upper)
                B = ctl[0,1,j,:] # Root LE upper +1
                C = ctl[1,-2,j,:]# Root LE lower +1

                # Area = 0.5*abs( xA*yC - xAyB + xByA - xByC + xCyB - xCyA )

                A1 = A[0]*C[1] - A[0]*B[1] + B[0]*A[1] -B[0]*C[1] + C[0]*B[1] - C[0]*A[1]
                A2 = A[1]*C[2] - A[1]*B[2] + B[1]*A[2] -B[1]*C[2] + C[1]*B[2] - C[1]*A[2]
                A3 = A[0]*C[2] - A[0]*B[2] + B[0]*A[2] -B[0]*C[2] + C[0]*B[2] - C[0]*A[2]

                g_con[index:index+3,i] = imag(array([A1,A2,A3]))/imag(h)
                index += 3
            # end for
            x[i] -= h
        # end for
        return g_obj,g_con,False

    def __unpack_x(self,x,dtype='d'):
        ctl = zeros((self.Nctlu,self.Nctlv,3),dtype)
        N = self.Nctlu*self.Nctlv
        for idim in xrange(3):
            ctl[:,:,idim] = reshape(x[3*N + idim*N : 3*N + idim*N + N],[self.Nctlu,self.Nctlv])
        # end for
        return ctl

    def getValue(self,u,v):
        
        '''Get the value of the spline at point u,v'''
        x = zeros([3])
        for idim in xrange(3):
            x[idim] = pyspline.b2val(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        
        return x


    def getValueV(self,u,v):
        '''Get the value of a spline at vector of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same length'
        x = zeros(len(u),3)
        for idim in xrange(3):
            x[:,idim] = pyspline.b2valv(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        return x

    def getValueM(self,u,v):
        '''Get the value of a spline at matrix of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same shape'
        x = zeros((u.shape[0],u.shape[1],3))
        for idim in xrange(3):
            x[:,:,idim] = pyspline.b2valm(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return x

    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((3,2))
        
        for idim in xrange(3):
            J[idim,0] = pyspline.b2val(u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            J[idim,1] = pyspline.b2val(u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.bcoef_x)

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

    def writeTecplot(self,handle=None,file_name=None):
        
        if not(handle) and not(file_name):
            print 'A file handle OR file_name must be given'
            sys.exit(0)

        if handle:
            pass
        else:
            f = open(file_name,'w')
            handle.write ('VARIABLES = "X", "Y","Z"\n')

        handle.write('Zone T=%s I=%d J = %d\n'%('orig_data',self.Nu,self.Nv))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(self.Nv):
            for i in xrange(self.Nu):
                handle.write('%f %f %f \n'%(self.x0[i,j,0],self.x0[i,j,1],self.x0[i,j,2]))
            # end for
        # end for 
        u_plot = linspace(0,1,25)
        u_plot = self.u0
        #u_plot = 0.5*(1-cos(linspace(0,pi,50)))
        v_plot = linspace(0,1,25)
        v_plot = self.v0
        # Dump re-interpolated surface
        handle.write('Zone T=%s I=%d J = %d\n'%('interpolated',len(u_plot),len(v_plot)))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(len(v_plot)):
            for i in xrange(len(u_plot)):
                for idim in xrange(3):
                    handle.write('%f '%(pyspline.b2val(u_plot[i],v_plot[j],0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])))
                # end for 
                handle.write('\n')
            # end for
        # end for 

        # Dump Control Points
        handle.write('Zone T=%s I=%d J = %d\n'%('control_pts',self.Nctlu,self.Nctlv))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                handle.write('%f %f %f \n'%(self.coef[i,j,0],self.coef[i,j,1],self.coef[i,j,2]))
            # end for
        # end for 


    def writeIGES(self,file_name):

        f = open(file_name,'w')
        f.write('                                                                        S      1\n')
        f.write('1H,,1H;,7H128-000,11H128-000.IGS,9H{unknown},9H{unknown},16,6,15,13,15, G      1\n')
        f.write('7H128-000,1.,1,4HINCH,8,0.016,15H19970830.165254,0.0001,0.,             G      2\n')
        f.write('21Hdennette@wiz-worx.com,23HLegacy PDD AP Committee,11,3,               G      3\n')
        f.write('13H920717.080000,23HMIL-PRF-28000B0,CLASS 1;                            G      4\n')
        
        f.write('     128       1       0       1       0       0       0       000000001D      1\n')
        f.write('     128       0       2     128       0                                D      2\n')
        counter = 1
        f.write('%3d,%4d,%4d,%4d,%4d,%1d,%1d,%1d,%1d,%1d,%37s1P%7d\n'%(128,self.Nu-1,self.Nv-1,self.ku-1,self.kv-1,0,0,1,0,0,' ',counter))
        counter += 1
        # Now do the knots

        pos_counter = 0

        for i in xrange(len(self.tu)):
            pos_counter += 1
            f.write('%12.6g,'%(self.tu[i]))
            if mod(pos_counter,5) == 0:
                f.write('      1P%7d\n'%(counter))
                counter += 1
                pos_counter = 0


        for i in xrange(len(self.tv)):
            pos_counter += 1
            f.write('%12.6g,'%(self.tv[i]))
            if mod(pos_counter,5) == 0:
                f.write('      1P%7d\n'%(counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for

        for i in xrange(self.Nu*self.Nv):
            pos_counter += 1
            f.write('%12.6g,'%(1.0))
            if mod(pos_counter,5) == 0:
                f.write('      1P%7d\n'%(counter))
                counter += 1
                pos_counter = 0
            # end if
        # end for 
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                for idim in xrange(3):
                    pos_counter += 1
                    f.write('%12.6g,'%(self.coef[0,i,j,idim]))
                    if mod(pos_counter,5) == 0:
                        f.write('      1P%7d\n'%(counter))
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
                f.write('%12.6g,'%(self.u0[0,0]))
            if i == 1:
                f.write('%12.6g,'%(self.u0[0,-1]))
            if i == 2:
                f.write('%12.6g,'%(self.v0[0,0]))
            if i == 3:
                f.write('%12.6g;'%(self.v0[0,-1])) # semi-colon for the last entity
            if mod(pos_counter,5)==0:
                f.write('      1P%7d\n'%(counter))
                counter += 1
                pos_counter = 0
            else: # We have to close it up anyway
                if i ==3:
                    for j  in xrange(5-pos_counter):
                        f.write('%13s'%(' '))
                    f.write('      1P%7d\n'%(counter))
                    pos_counter = 0

                    
        # Output the Terminate Statement
        f.write('S%7dG%7dD%7dP%7d%40sT%7s\n'%(1,4,2,counter,' ',' '))

        f.close()




#==============================================================================
# Class Test
#==============================================================================
if __name__ == '__main__':
	
    # Run a Simple Test Case
    print 'Testing pySpline...\n'
    print 'There is an example in the ./example directory'




