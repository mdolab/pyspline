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
from numpy.linalg import lstsq

import pyspline_cs
import pyspline

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
        sys.stdout.write('pySpline Type: %s. '%(task))

        # This is the optional additional info required for use with pyGeo
        self.edge_con = [ [], [], [], []] #Is edge connected to another edge 
        self.master_edge = [True,True,True,True] # Is edge a Master Edge?...it is unless it its driven
        self.dir = [None,None,None,None]
        self.edge_type = [None,None,None,None]
        self.node_con = [[],[],[],[]]
        self.master_node = [None,None,None,None] # Is each node a master node?

        self.Nu_free = None
        self.Nv_free = None
        self.N_free  = None
        self.Nctlu_free = None
        self.Nctlv_free = None
        self.Nctl_free  = None
        self.links = None
        self.ref_axis_dir = None
        
        if task == 'create':
            assert 'ku' in kwargs and 'kv' in kwargs  and 'tu' in kwargs \
                and 'tv' in kwargs and 'coef' in kwargs and 'range' in kwargs, \
                'Error: ku,kv,tu,tv,coef and range MUST be defined for task=\'create\''
            sys.stdout.write('\n')
            self.u = None
            self.v = None
            self.X = None
            self.Nu = None
            self.Nv = None
            self.ku = kwargs['ku'] 
            self.kv = kwargs['kv']
            self.tu = array(kwargs['tu'])
            self.tv = array(kwargs['tv'])
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

            # Set a linear version of U and V
            [self.V, self.U] = meshgrid(self.v,self.u)
            self.orig_data = True
            
            self.coef = zeros((self.Nctlu,self.Nctlv,self.nDim),'D')
            self.range = array([self.u[0],self.u[-1],self.v[0],self.v[-1]])

            #Sanity check to make sure k is less than N
            if self.Nu <= self.ku:
                self.ku = self.Nu
            if self.Nv <= self.kv:
                self.kv = self.Nv
            
            timeA = time.time()
            for idim in xrange(self.nDim):
                self.tu,self.tv,self.coef[:,:,idim]= pyspline_cs.b2ink(self.u,self.v,self.X[:,:,idim],self.ku,self.kv)
            #end for
            
            sys.stdout.write(' Interpolate Time: %6.5f s\n'%(time.time()-timeA))
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
            self.orig_data = True

            # Sanity Check on Inputs
            if self.Nctlu > self.Nu:
                self.Nctlu  = self.Nu
            if self.Nctlv > self.Nv:
                self.Nctlv = self.Nv
            # end if
            
            if self.Nu <= self.ku:
                self.ku = self.Nu
            if self.Nv <= self.kv:
                self.kv = self.Nv

            if 'u' in kwargs and 'v' in kwargs:
                self.u = kwargs['u']
                self.v = kwargs['v']
            else:
                self._calcParameterization()
            # end if
            
            self.range = array([self.u[0],self.u[-1],self.v[0],self.v[-1]])

            # Set a linear version of U and V
            [self.V, self.U] = meshgrid(self.v,self.u)

           #Calculate the knot vector and Jacobian
            sys.stdout.write(' Calculating: knots, ')


            self._calcKnots()

            sys.stdout.write(' jacobian, ')
            self._calcJacobian()

#             # Lets do a lms 
            timeA = time.time()
            self.coef = zeros([self.Nctlu,self.Nctlv,self.nDim],'D')
            for idim in xrange(self.nDim):
                self.coef[:,:,idim] = reshape(lstsq(self.J,self.X[:,:,idim].flatten())[0],[self.Nctlu,self.Nctlv])

            sys.stdout.write(' LMS Fit Time: %6.5f s\n'%(time.time()-timeA))
            
            return 

        print 'Error: task is not understood. Type must be \'lms\',\'interpolate\' or \'create\''
        sys.exit(0)
        return

    def _calcParameterization(self):

        u = zeros(self.Nu,'D')
        singular_counter = 0

        for j in xrange(self.Nv): #loop over each v, and average the 'u' parameter 
            temp = zeros(self.Nu,'D')

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
        
        v = zeros(self.Nv,'D')
        singular_counter = 0
        for i in xrange(self.Nu): #loop over each v, and average the 'u' parameter 
            temp = zeros(self.Nv,'D')
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

    def _calcNFree(self):
        
        '''Internal function to calculate how many "free" control points we
        have. i.e. the master edges'''
        Nctlu = self.Nctlu
        Nctlv = self.Nctlv
        Nu    = self.Nu
        Nv    = self.Nv
        
        if self.master_edge[0] == False:
            Nv -= 1
            Nctlv -= 1
            
        if self.master_edge[1] == False:
            Nv -= 1
            Nctlv -= 1

        if self.master_edge[2] == False:
            Nu -= 1
            Nctlu -= 1

        if self.master_edge[3] == False:
            Nu -= 1
            Nctlu -= 1
            
        self.Nu_free = Nu
        self.Nv_free = Nv
        self.N_free = Nu*Nv
        
        self.Nctlu_free = Nctlu
        self.Nctlv_free = Nctlv
        self.Nctl_free = Nctlu*Nctlv
        return

    def _getCropData(self,X):
        '''This internal function gets the original data coorsponding but
        omits the data along edges which are not masters'''

        return eval('X'+self.slice_string)

    def getFreeCtl(self):

        '''Return Free control Points'''
        # This requires some explaining: the getFreeIndex function
        # generates the pythonic index string which coorsponds to the
        # section of the coefficient martrix which is free or
        # driving. Then to return these values we make a string that
        # contains 'self.coef' and then the slice_string which is
        # something like [:,1:-1]. This is only possible in an
        # interpreted language
        return eval('self.coef'+self.slice_string)


    def setFreeCtl(self,coef):

        '''Sets Free control Points'''
        # See getFreeCtl for how this works
        
        return eval('self.coef'+self.slice_string,'=coef')

    
    def _getFreeIndex(self):
        '''This internal function gets string slice which tells
        various functions how to slice out control points or data that
        coorspond to master edges'''
        # Use a Binary Tree-Type Search -> Faster -> log n scaling 

        slice_string = ''

        if self.master_edge[0] == True:
            if self.master_edge[1] == True:
                if self.master_edge[2] == True:
                    if self.master_edge[3] == True:
                        #[True,True,True,True]:
                        slice_string = '[:,:]'
                    else:
                        #[True,True,True,False]:
                         slice_string = '[0:-1,:]'
                    # end if
                else:
                    if self.master_edge[3] == True:
                        #[True,True,False,True]:
                        slice_string = '[1:,:]'
                    else:
                        #[True,True,False,False]:
                        slice_string = '[1:-1,:]'
                    # end if
                # end if
            else:
                if self.master_edge[2] == True:

                    if self.master_edge[3] == True:
                        #[True,False,True,True]:
                        slice_string = '[:,0:-1]'
                    else:
                        #[True,False,True,False]:
                        slice_string = '[0:-1,0:-1]'
                    # end if
                else:
                    if self.master_edge[3] == True:
                        #[True,False,False,True]:
                        slice_string = '[1:,0:-1]'
                    else:
                        #[True,False,False,False]:
                        slice_string = '[1:-1,0:-1]'
                    # end if
                # end if
            # end if
        else:
            if self.master_edge[1] == True:
                if self.master_edge[2] == True:
                    if self.master_edge[3] == True:
                        #[False,True,True,True]:
                        slice_string = '[:,1:]'
                    else:
                        #[False,True,True,False]:
                        slice_string = '[0:-1,1:]'
                    # end if
                else:
                    if self.master_edge[3] == True:
                        # [False,True,False,True]:
                        slice_string = '[1:,1]'
                    else:
                        # [False,True,False,False]:
                        slice_string = '[1:-1,1:]'
                    # end if
                # end if
            else:
                if self.master_edge[2] == True:
                    if self.master_edge[3] == True:
                        #[False,False,True,True]:
                        slice_string = '[:,1:-1]'
                    else:
                        # [False,False,True,False]:
                        slice_string = '[0:-1:,1:-1]'
                    # end if
                else:
                    if self.master_edge[3] == True:
                        #[False,False,False,True]:
                        slice_string = '[1:,1:-1]'
                    else:
                        # [False,False,False,False]:
                        slice_string = '[1:-1,1:-1]'
                    # end if
                # end if
            # end if
        #end if

        self.slice_string = slice_string
        return


    def _calcKnots(self):

        '''Find the initial knot vector for this problem'''

        self.tu = pyspline.knots(self.u,self.Nctlu,self.ku)
        self.tv = pyspline.knots(self.v,self.Nctlv,self.kv)
        
        return

    def _calcJacobian(self):
        
        # Calculate the jacobian J, for fixed t and s
        self.J = zeros([self.Nu*self.Nv,self.Nctlu*self.Nctlv],'D')
        ctl = zeros([self.Nctlu,self.Nctlv],'D')

        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                ctl[i,j] += 1
                self.J[:,i*self.Nctlv + j] = pyspline_cs.b2valv(self.U.flatten(),self.V.flatten(),0,0,self.tu,self.tv,self.ku,self.kv,ctl)
                ctl[i,j] -= 1
                # end for
            # end for 
        # end for
        return

    def checkCtl(self,i,j):
        ''' This function checks to see if a control point is a master. i.e. Internal or on a master edge'''

        if i > 0 and i < self.Nctlu-1 and j > 0 and j < self.Nctlv-1:
            #This is the internal case...
            return 0,None,None

        # Now we need to check if its on a master edge (but not a corner)
        
        if j == 0 and i > 0 and i < self.Nctlu-1:
            if self.master_edge[0]:
                return 0,[self.edge_con[0],i,self.dir[0],self.edge_type[0]],None  #Face,Edge,Node
            else:
                return 1,None,None
        
        if j == self.Nctlv-1 and i > 0 and i < self.Nctlu-1:
            if self.master_edge[1]:
                return 0,[self.edge_con[1],i,self.dir[1],self.edge_type[1]],None  #Face,Edge,Node
            else:
                return 1,None,None

        if i == 0 and j > 0 and j < self.Nctlv-1:
            if self.master_edge[2]:
                return 0,[self.edge_con[2],j,self.dir[2],self.edge_type[2]],None  #Face,Edge,Node
            else:
                return 1,None,None
        if i == self.Nctlu-1 and j > 0 and j < self.Nctlv-1 :
            if self.master_edge[3]:
                return 0,[self.edge_con[3],j,self.dir[3],self.edge_type[3]],None  #Face,Edge,Node
            else:
                return 1,None,None
        
        # Now check Corners...Must be on BOTH master edges

        if j == 0 and i == 0:
            if self.master_node[0]:
                return 0,None,self.node_con[0]
            else:
                return 1,None,None
            
        if j == 0 and i == self.Nctlu-1:
            if self.master_node[1]:
                return 0,None,self.node_con[1]
            else:
                return 1,None,None
            
        if j == self.Nctlv-1 and i == 0:
            if self.master_node[2]:
                return 0,None,self.node_con[2]
            else:
                return 1, None, None

        if j == self.Nctlv-1 and i == self.Nctlu-1:
            if self.master_node[3]:
                return 0,None,self.node_con[3]
            else:
                return 1, None, None
        # else:


    def _calcCtlDeriv(self,i,j):
        '''Calculate the derivative wrt control point i,j'''
        ctl = zeros([self.Nctlu,self.Nctlv],'D')
        ctl[i,j] += 1
        
        return self._getCropData(pyspline_cs.b2valm(self.U,self.V,0,0,self.tu,self.tv,self.ku,self.kv,ctl)).flatten()

    def _calcCtlDerivNode(self,node):
        '''Calculate the derivative wrt control point i,j'''
        ctl = zeros([self.Nctlu,self.Nctlv],'D')
        if node == 0:
            ctl[0,0] += 1
        elif node == 1:
            ctl[-1,0] += 1
        elif node == 2:
            ctl[0,-1] += 1
        else:
            ctl[-1,-1] += 1
        
        return self._getCropData(pyspline_cs.b2valm(self.U,self.V,0,0,self.tu,self.tv,self.ku,self.kv,ctl)).flatten()

    def _calcCtlDerivEdge(self,edge,index,dir):
        '''Calculate the derivative wrt control point i,j'''
        ctl = zeros([self.Nctlu,self.Nctlv],'D')
        if edge == 0:
            if dir==1:
                ctl[index,0] += 1
            else:
                ctl[self.Nctlu-index-1,0] += 1

        elif edge == 1:
            if dir ==1:
                ctl[index,-1] += 1
            else:
                ctl[self.Nctlu-index-1,-1] += 1

        elif edge == 2:
            if dir == 1:
                ctl[0,index] += 1
            else:
                ctl[0,self.Nctlv-index-1] += 1
        else:
            if dir == 1:
                ctl[-1,index] += 1
            else:
                clt[-1,self.Nctlv-index-1] += 1
        
        
        return self._getCropData(pyspline_cs.b2valm(self.U,self.V,0,0,self.tu,self.tv,self.ku,self.kv,ctl)).flatten()


    def associateRefAxis(self,ref_axis):#,section=None):
        '''Create links between ref_axis and this patch'''

        #print 'associateRefAxis'
        # First we Need to deduce along which direction (u or v) the
        # reference axis is directed.  First estimate Over what
        # portion the surface and ref axis co-inside

        # Find the s range for the 4 corners

        p = zeros(4,'D')
        
        p[0],conv,D,esp = ref_axis.xs.projectPoint(self.coef[0 , 0])
        p[1],conv,D,esp = ref_axis.xs.projectPoint(self.coef[0 ,-1])
        p[2],conv,D,esp = ref_axis.xs.projectPoint(self.coef[-1, 0])
        p[3],conv,D,esp = ref_axis.xs.projectPoint(self.coef[-1,-1])

        min_s = min(p)
        max_s = max(p)
       
        #print 'p is:',p
       

        # Now we know over what portion of the ref axis we are dealing
        # with.  Take 3 Normal Vectors
        N = 3
        sn = linspace(min_s,max_s,N)
        dn = zeros((N,3),'D')
        for i in xrange(N):
            dn[i,:] = ref_axis.xs.getDerivative(sn[i])


        # Now Do two tests: Take three points in u and test 3 groups
        # against dn and take three points in v and test the  3 groups
        # again
    
        u_dot_tot = 0
        s = linspace(0,1,N)
    
        for i in xrange(N):
            for n in xrange(N):
                du,dv = self.getDerivative(s[i],s[n])
                u_dot_tot += abs(dot(du,dn[n,:]))
            # end for
        # end for

        v_dot_tot = 0
        for j in xrange(N):
            for n in xrange(N):
                du,dv = self.getDerivative(s[n],s[j])
                v_dot_tot += abs(dot(dv,dn[n,:]))
            # end for
        # end for
        
        self.links = []

        
        if v_dot_tot > u_dot_tot:
            #print 'along v:'
            for j in xrange(self.Nctlv):
                # Create a line (k=2 spline) for the control points along U
                
                ctl_line = linear_spline('lms',X=self.coef[:,j],Nctl=2,k=2)
                # Now find the minimum distance between midpoint and the ref_axis
                s,D,conv,eps = ref_axis.xs.projectPoint(ctl_line.getValue(1.0))
                if s > max_s:
                    s = max_s
                if s < min_s:
                    s = min_s
            #    print 's is:',s
                # Now loop over the control points to set links
                base_point = ref_axis.xs.getValue(s)
                M = ref_axis.getRotMatrixGlobalToLocal(s)
                
                for i in xrange(self.Nctlu):
                    D = self.coef[i,j] - base_point
                    D = dot(M,D) #Rotate to local frame
                    self.links.append([s,D])
                # end if
            # end for
            self.ref_axis_dir = 1 # Ref axis is aligned along v
        else:
            #print 'along u:'
            for i in xrange(self.Nctlu):
                # Create a line (k=2 spline) for the control points along V

                ctl_line = linear_spline('lms',X=self.coef[i,:],Nctl=2,k=2)
                # Now find the minimum distance between midpoint and the ref_axis
                s,D,conv,eps = ref_axis.xs.projectPoint(ctl_line.getValue(0.5))
                if s > max_s:
                    s = max_s
                if s < min_s:
                    s = min_s
                # Now loop over the control points to set links
                base_point = ref_axis.xs.getValue(s)
                M = ref_axis.getRotMatrixGlobalToLocal(s)
                
                for j in xrange(self.Nctlv):
                    D = self.coef[i,j] - base_point
                    D = dot(M,D) #Rotate to local frame
                    self.links.append([s,D])
                # end if
            # end for
            self.ref_axis_dir = 0 # Ref axis is aligned along u


    def associateRefAxis2(self,ref_axis,s):
        '''Create links between ref_axis and this patch.'''

        self.ref_axis_dir = 1 #since it is ONLY called from Lifting Surface
        self.links = []
        assert len(s) == self.Nctlv, 's must be of length Nctlv'
        for j in xrange(self.Nctlv):
            # Now loop over the control points to set links
            base_point = ref_axis.xs.getValue(s[j])
            M = ref_axis.getRotMatrixGlobalToLocal(s[j])

            for i in xrange(self.Nctlu):
                D = self.coef[i,j] - base_point
                D = dot(M,D) #Rotate to local frame
                self.links.append([s[j],D])
            # end if
        # end for

    def update(self,ref_axis):
        '''Update the control points with ref_axis:'''
        counter = 0
        print 'update from cs'
        if self.ref_axis_dir == 1:
            for j in xrange(self.Nctlv):
                s = self.links[counter][0]
                M = ref_axis.getRotMatrixLocalToGloabl(s)
                X_base = ref_axis.xs.getValue(s)
                for i in xrange(self.Nctlu):
                   self.coef[i,j,:] = X_base + dot(M,self.links[counter][1])*ref_axis.scales(s)
                   counter += 1
                # end for
            # end for
        else:
            for i in xrange(self.Nctlu):
                s = self.links[counter][0]
                M = ref_axis.getRotMatrixLocalToGloabl(s)
                X_base = ref_axis.xs.getValue(s)
                
                for j in xrange(self.Nctlv):
                    self.coef[i,j,:] = X_base + dot(M,self.links[counter][1])*ref_axis.scales(s)
                    counter += 1
                # end for
            # end for

    def updateSurfacePoints(self,delta):
        '''Update the control points on surface deltas normal to the surface'''
        update = zeros((self.Nctlu,self.Nctlv),'D') 
        exec('update'+self.slice_string+' = delta')

# ------------- Python Code ------------------------------------------------

#  for i in xrange(self.Nctlu):
#             for j in xrange(self.Nctlv):
#                 # The coordinates where we get the normal are NOT
#                 # stictly correct...the last ones will run out into
#                 # the 1's at the end of the knot vector
                
#                 n = self.getNormal(self.tu[self.ku-2 + i],self.tv[self.kv-2+j])

#                 self.coef[i,j] -= update[i,j] * n
#                 # The MINUS is because the surface normals MUST point inward
                
#             # end for
#         # end for

# -------------------------------------------------------------------------

        # Call the fortran function instead
        self.coef = pyspline_cs.updatesurfacepoints(self.coef,update,self.tu,self.tv,self.ku,self.kv)

        return

    def getValueEdge(self,edge,s):
        '''Get the value of the spline on edge, edge=0,1,2,3 where
        edges are oriented in the standard counter-clockwise fashion.'''

        if edge == 0:
            return self.getValue(s,0)
        elif edge == 1:
            return self.getValue(s,1)
        elif edge == 2:
            return self.getValue(0,s)
        elif edge ==3:
            return self.getValue(1,s)
        else:
            print 'Edge must be between 0 and 3'
            sys.exit(1)
            return
        #end if 

    def getValueCorner(self,corner):
        '''Get the value of the spline on corner i where nodes are oriented in
        the standard counter-clockwise fashion.'''

        if corner == 0:
            return self.getValue(0,0)
        elif corner == 1:
            return self.getValue(1,0)
        elif corner == 2:
            return self.getValue(1,1)
        elif corner ==3:
            return self.getValue(0,1)
        else:
            print 'Corner must be between 0 and 3'
            sys.exit(1)
            return
        #end if 

    def getCoefEdge(self,edge):
        '''Get Coef along edge edge'''

        if   edge == 0:
            return self.coef[:,0,:]
        elif edge == 1:
            return self.coef[:,-1,:]
        elif edge == 2:
            return self.coef[0,:,:]
        else:
            return self.coef[-1,:,:]
        
    def setCoefEdge(self,edge,new_coef):
        '''Set Coef along edge edge'''

        if   edge == 0:
            self.coef[:,0,:] = new_coef
        elif edge == 1:
            self.coef[:,-1,:] = new_coef
        elif edge == 2:
            self.coef[0,:,:] = new_coef
        else:
            self.coef[-1,:,:] = new_coef
        return

    def getNormal(self,u,v):
        '''Get the normal at the surface point u,v'''
        du,dv = self.getDerivative(u,v)

        # Now cross and normalize
        n = zeros(3,'D')
        n[0] = du[1]*dv[2]-du[2]*dv[1]
        n[1] = du[2]*dv[0]-du[0]*dv[2]
        n[2] = du[0]*dv[1]-du[1]*dv[0]
        
        n/= sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
        return n

    def getValue(self,u,v):
        
        '''Get the value of the spline at point u,v'''
        x = zeros([self.nDim],'D')
        for idim in xrange(self.nDim):
            x[idim] = pyspline_cs.b2val(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        
        return x

    def getValueV(self,u,v):
        '''Get the value of a spline at vector of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same length'
        x = zeros((len(u),self.nDim),'D')
        for idim in xrange(self.nDim):
            x[:,idim] = pyspline_cs.b2valv(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        return x

    def getValueM(self,u,v):
        '''Get the value of a spline at matrix of points u,v'''
        assert u.shape == v.shape, 'u and v must be the same shape'
        x = zeros((u.shape[0],u.shape[1],self.nDim),'D')
        for idim in xrange(self.nDim):
            x[:,:,idim] = pyspline_cs.b2valm(u,v,0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return x

    def getDerivative(self,u,v):
        
        '''Get the value of the derivative of spline at point u,v'''
        du = zeros(self.nDim,'D')
        dv = zeros(self.nDim,'D')
        for idim in xrange(self.nDim):
            du[idim] = pyspline_cs.b2val(u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            dv[idim] = pyspline_cs.b2val(u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
        return du,dv


    def getJacobian(self,u,v):
        
        '''Get the jacobian at point u,v'''

        J = zeros((self.nDim,2),'D')
        
        for idim in xrange(self.nDim):
            J[idim,0] = pyspline_cs.b2val(u,v,1,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])
            J[idim,1] = pyspline_cs.b2val(u,v,0,1,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim])

        return J

    def projectPoint(self,x0,u0=0.5,v0=0.5,Niter=25,tol=1e-6):

        '''Project a point x0 onto the surface. i.e. Find the point on
        the surface that minimizes the distance from x0 to surface(u,v)'''

        # We will use a starting point u0,v0 if given

        # Check we have the same dimension:
        assert len(x0) == self.nDim,'Dimension of x0 and the dimension of spline must be the same'
        converged = False
        #print 'x0:',x0
        for i in xrange(Niter):
            D = x0 - self.getValue(u0,v0)
            Ydotu,Ydotv = self.getDerivative(u0,v0)
            updateu = dot(D,Ydotu)/(sqrt(dot(Ydotu,Ydotu)))#/self.length
            updatev = dot(D,Ydotv)/(sqrt(dot(Ydotv,Ydotv)))#/self.length
            # Check to see if update went too far

            if u0+updateu > self.range[1]:
                updateu = self.range[1]-u0
            elif u0+updateu< self.range[0]:
                # Only put the update up to the end
                updateu = self.range[0]-u0
            # end if

            if v0+updatev > self.range[3]:
                updatev = self.range[3]-v0
            elif v0+updatev< self.range[2]:
                # Only put the update up to the end
                updatev = self.range[2]-v0
            # end if
                                
            D2 = x0-self.getValue(u0+updateu,v0+updatev)
            if abs(dot(D2,D2)) > abs(dot(D,D)):
                #print 'half update'
                updateu /= 2
                updatev /= 2
            # end if

            if abs(updateu)<tol and abs(updatev)<tol:
                u0 += updateu
                v0 += updatev
                
                converged=True
                D = x0-self.getValue(u0,v0) # Get the final Difference
                break
            else:
                u0 += updateu
                v0 += updatev
                #print 'u0,v0',u0,v0,'updates:',updateu,updatev,'val:',self.getValue(u0,v0)
            # end if
        # end for
        #print 'Took %d iterations'%(i)
        return u0,v0,D,converged



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

            f = mat(zeros((3,1),'D'))
            f[0] = x[0]-(x0[0]+r[0]*s)
            f[1] = x[1]-(x0[1]+r[1]*s)
            f[2] = x[2]-(x0[2]+r[2]*s)

            J = self.getJacobian(u,v)
            A = mat(zeros((3,3),'D'))
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
                    handle.write('%f %f %f \n'%(real(self.X[i,j,0]),real(self.X[i,j,1]),real(self.X[i,j,2])))
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
            
        u_plot = linspace(self.range[0],self.range[1],25)
        v_plot = linspace(self.range[2],self.range[3],25)

        # Dump re-interpolated surface
        handle.write('Zone T=%s I=%d J = %d\n'%('interpolated',len(u_plot),len(v_plot)))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(len(v_plot)):
            for i in xrange(len(u_plot)):
                for idim in xrange(self.nDim):
                    handle.write('%f '%(real(pyspline_cs.b2val(u_plot[i],v_plot[j],0,0,self.tu,self.tv,self.ku,self.kv,self.coef[:,:,idim]))))
                # end for 
                handle.write('\n')
            # end for
        # end for 

        # Dump Control Points (Always have these :-) ) 
        handle.write('Zone T=%s I=%d J = %d\n'%('control_pts',self.Nctlu,self.Nctlv))
        handle.write('DATAPACKING=POINT\n')
        for j in xrange(self.Nctlv):
            for i in xrange(self.Nctlu):
                handle.write('%f %f %f \n'%(real(self.coef[i,j,0]),real(self.coef[i,j,1]),real(self.coef[i,j,2])))
            # end for
        # end for 

        return


    def writeTecplotEdge(self,handle,edge,*args,**kwargs):
        '''Dump out a linear zone along edge used for visualizing edge connections'''

        N = 25
        if 'name' in kwargs:
            handle.write('Zone T=%s I=%d\n'%(kwargs['name'],N))
        else:
            handle.write('Zone T=%s I=%d\n'%('edge',N))
        handle.write('DATAPACKING=POINT\n')
        s = linspace(0,1,N)
        for i in xrange(N):
            value = self.getValueEdge(edge,s[i])
            handle.write('%f %f %f \n'%(real(value[0]),real(value[1]),real(value[2])))
        # end for
        return


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
            # end if
        # end for

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
                    # end for
                    pos_counter = 0
                    handle.write('%7dP%7d\n'%(Pcount,counter))
                    counter += 1
                # end if
            # end if
        # end for

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
        #print 'pySpline Class Initialization Type: %s'%(task)

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
                    print 'Erorr: For 1D mapping splines, the basis, s, must be given as input'
                    sys.exit(1)
                # end if
            # end if
                
            self.Nctl = self.N
            self.range = array([self.s[0],self.s[-1]])
            self.orig_data = True

            # Sanity check to make sure k is less than N
            if self.N <= self.k:
                self.k = self.N
            # end if

            # Generate the knot vector
            self.t = pyspline_cs.bknot(self.s,self.k)
            if self.nDim > 1:
                self.coef = zeros((self.Nctl,self.nDim),'D')
                for idim in xrange(self.nDim):
                    self.coef[:,idim]= pyspline_cs.bintk(self.s,self.X[:,idim],self.t,self.k)
            else:
                self.coef = pyspline_cs.bintk(self.s,self.X,self.t,self.k)
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
                    print 'Erorr: For 1D mapping splines, the basis, s, must be given as input'
                    sys.exit(1)
                # end if
            # end if
                
            self.Nctl = self.N
            self.range = array([self.s[0],self.s[-1]])
            self.orig_data = True

            # Sanity check to make sure k is less than N
            if self.N <= self.k:
                self.k = self.N
            # end if

            if self.Nctl > self.N:
                self.Nctl  = self.N

            # Generate the knot vector
            self.t = pyspline_cs.bknot(self.s,self.k)
            if self.nDim > 1:
                # Calculate the Jacobian
                self._calcJacobian()
                self.coef = zeros((self.Nctl,self.nDim),'D')
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
        assert len(x0) == self.nDim,'Dimension of x0 and the dimension of spline must be the same'
        converged = False

        for i in xrange(Niter):
            D = x0 - self.getValue(s0)
            Ydot = self.getDerivative(s0)
            update = dot(D,Ydot)/(sqrt(dot(Ydot,Ydot)))/self.length

            # Check to see if update went too far

            if s0+update > self.range[1]:
                update = self.range[1]-s0
            elif s0+update< self.range[0]:
                # Only put the update up to the end
                update = self.range[0]-s0
                                
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
        self.s = zeros(self.N,'D');
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
        s = zeros(self.N,'D');
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
            x = pyspline_cs.bvalu(self.t,self.coef,self.k,0,s)
        else:
            x = zeros([self.nDim],'D')
            for idim in xrange(self.nDim):
                x[idim] = pyspline_cs.bvalu(self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

   
    def getValue(self,s):
        
        '''Get the value of the spline at point u,v'''
        if self.nDim == 1:
            x = pyspline_cs.bvalu(self.t,self.coef,self.k,0,s)
        else:
            x = zeros([self.nDim],'D')
            for idim in xrange(self.nDim):
                x[idim] = pyspline_cs.bvalu(self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

    def getValueV(self,s):
        '''Get the value of a spline at vector of points s'''
        if self.nDim == 1:
            x = pyspline_cs.bvaluv(self.t,self.coef,self.k,0,s)
        else:
            x = zeros((len(s),self.nDim),'D')
            for idim in xrange(self.nDim):
                x[:,idim] = pyspline_cs.bvaluv(self.t,self.coef[:,idim],self.k,0,s)
            # end for
        # end if
        return x

 
    def getDerivative(self,s):
        
        '''Get the value of the derivative of spline at point u,v'''
        x = zeros(self.nDim,'D')
        for idim in xrange(self.nDim):
            x[idim] = pyspline_cs.bvalu(self.t,self.coef[:,idim],self.k,1,s)
            
        return x        


    def _calcJacobian(self):
        
        # Calculate the jacobian J, for fixed t and s
        self.J = zeros([self.N,self.Nctl],'D')
        ctl = zeros([self.Nctl],'D')

        for i in xrange(self.Nctl):
            ctl[i] += 1
            self.J[:,i] = pyspline_cs.bvaluv(self.t,ctl,self.k,0,self.s)
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
                handle.write('%f '%(pyspline_cs.bvalu(self.t,self.coef[:,idim],self.k,0,s_plot[i])))
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
    print 'Testing pySplineCS..\n'
    print 'There is an example in the ./example directory'

