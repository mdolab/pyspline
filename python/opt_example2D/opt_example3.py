# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot, cross
from scipy.io import read_array
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline

# =============================================================================
# Extension modules
# =============================================================================

# pySpline
sys.path.append('../../python')
import pyspline
import pyspline_cs

#pyOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/')
from pyOpt_optimization import Optimization

#pySNOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/pySNOPT')
from pySNOPT import SNOPT

'''This script runs an example of fitting a cubic spline to a set of
3D wing coordinates that were generated using aerosurf. The script is
meant to show how to fit a constrained surface using a non-linear
optimizer'''

# Global Variables

Nctlu =13
Nctlv =7
ku=4
kv=4

Nu = 26
Nv = 25

# ==============================
# Read in the Data 
# ==============================

# X holds all the coordinates
# First index is isurf -> surface id
# Second index is U
# Third index is V
# Fourth index is 0,1,2 for x,y,z

Nsurf = 2
X = zeros([Nsurf,Nu,Nv,3])
X_flat = zeros([Nsurf,Nu*Nv,3])
# Top surface of wing
print 'reading upper surface...'
data = read_array(file("upper_surface.inp"))

for j in xrange(Nv):
    for i in xrange(Nu):
        X[0,i,j,:] = data[j*Nu + i ,:]
    #end for
#end for

# Bottom surface of wing
print 'reading bottom surface...'
data = read_array(file("lower_surface.inp"))
for j in xrange(Nv):
    for i in xrange(Nu):
        X[1,i,j,:] = data[j*Nu + i, :]
    #end for
#end for

# Check the Data

f = open('output.dat','w')
for isurf in xrange(Nsurf):
    f.write ('VARIABLES = "X", "Y","Z"\n')
    f.write('Zone I=%d J = %d\n'%(Nu,Nv))
    f.write('DATAPACKING=POINT\n')
    for j in xrange(Nv):
        for i in xrange(Nu):
            f.write('%f %f %f \n'%(X[isurf,i,j,0],X[isurf,i,j,1],X[isurf,i,j,2]))

# ==============================
# Find the 'u' parameters
# ==============================

u=zeros((Nsurf,Nu))
for isurf in xrange(Nsurf):
    x0=X[isurf,:,0,0];y0=X[isurf,:,0,1];z0=X[isurf,:,0,2]
    for i in xrange(Nu-1):
        u[isurf,i+1] = u[isurf,i] + sqrt((x0[i+1]-x0[i])**2 + \
                                   (y0[i+1]-y0[i])**2+ \
                                   (z0[i+1]-z0[i])**2)
    # end for
    u[isurf,:] = u[isurf,:]/u[isurf,-1] #Normalize u
# end for

# ==============================
# Find the 'v' parameters
# ==============================
v=zeros((Nsurf,Nv))
for isurf in xrange(Nsurf):
    x0=X[isurf,0,:,0];y0=X[isurf,0,:,1];z0=X[isurf,0,:,2]
    for i in xrange(Nv-1):
        v[isurf,i+1] = v[isurf,i] + sqrt((x0[i+1]-x0[i])**2 + \
                                   (y0[i+1]-y0[i])**2+ \
                                   (z0[i+1]-z0[i])**2)
    # end for
    v[isurf,:] = v[isurf,:]/v[isurf,-1] #Normalize u
# end for
# ==============================
# Find the Initial Knot Vectors
# ==============================

# U knots

tu0= zeros(Nctlu + ku)
tu0[ku-1:Nctlu+1]  = 0.5*(1-cos(linspace(0,pi,Nctlu-ku+2)))
tu0[0:ku] = 0
tu0[Nctlu:Nctlu+ku] = 1.0#+0.1*(tu0[Nctlu]-tu0[Nctlu-1])

# V Knots

tv0= zeros(Nctlv + kv)
tv0[kv-1:Nctlv+1]=  linspace(0,1,Nctlv-kv+2)
#tv0[kv-1:Nctlv+1]  = 0.5*(1-cos(linspace(0,pi,Nctlv-kv+2)))
tv0[0:kv] = 0
tv0[Nctlv:Nctlv+kv] = 1.0#+0.1*(tv0[Nctlv]-tv0[Nctlv-1])

# Calculate the jacobian J, for fixed t and s

h = 1.0e-40j
J = zeros([Nsurf,Nu*Nv,Nctlu*Nctlv])
ctl = zeros([Nctlu,Nctlv],'D')

V = zeros((Nsurf,Nu,Nv))
U = zeros((Nsurf,Nu,Nv))

for isurf in xrange(Nsurf):
    [V[isurf], U[isurf]] = meshgrid(v[isurf,:],u[isurf,:])
    for j in xrange(Nctlv):
        for i in xrange(Nctlu):
            ctl[i,j] += h
            val = pyspline_cs.b2valv(U[isurf].flatten(),V[isurf].flatten(),0,0,tu0,tv0,ku,kv,ctl)
            ctl[i,j] -= h    
            J[isurf,:,i*Nctlv + j] = imag(val)/imag(h)
        # end for
    # end for 
# end for

global timeCounter 
timeCounter = 0.0

def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global timeCounter,B
    timeA = time.time()
 
# Unpack the x-values

    ctl = zeros((Nsurf,Nctlu,Nctlv,3))
    for isurf in xrange(Nsurf):
        for idim in xrange(3):
            ctl[isurf,:,:,idim] = reshape(x[isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv : isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv],[Nctlu,Nctlv])
        # end if
    #end if
  
    total = 0.0

    for isurf in xrange(Nsurf):
        for idim in xrange(3):
            total += sum((dot(J[isurf],ctl[isurf,:,:,idim].flatten()) - X[isurf,:,:,idim].flatten())**2)
        # end for
    # end for 
    f = total
    fail = False
    fcon = dot(B,x)

   #  # Calculate the LE constraint
    xA = ctl[0,0,0,:] # Root LE (upper)
    xB = ctl[0,1,0,:] # Root LE upper +1
    xC = ctl[1,-2,0,:]# Root LE lower +1

    v1 = xB-xA #vector from A to B
    v2 = xA-xC #vector from A to C
    #print 'norm',dot(cross(v1,v2),cross(v1,v2))
    # Now take cross product and quadrature sum
    xProduct = (v1[1]*v2[2]-v1[2]*v2[1])**2+ (-v1[0]*v2[2]+v1[2]*v2[0] )**2 +  ( v1[0]*v2[1] - v1[1]*v2[0])**2
    #fcon[-1] = (1.0e6*xProduct)

    
    fcon[-2] = (v1[0]/v2[0])/(v1[1]/v2[1])
    fcon[-1] = (v1[0]/v2[0])/(v1[2]/v2[2])
    #print 'vectors:',v1,v2
#     print 'xProduct',xProduct
#    print 'check',v1[0]/v2[0],v1[1]/v2[1],v1[2]/v2[2]
    timeCounter += time.time()-timeA
    return f,fcon,fail

def sens(x,f_obj,f_con):
    global timeCounter,B
    timeA = time.time()
    ndv = len(x)
    g_obj = zeros(ndv)
     # Unpack the x-values

    ctl = zeros((Nsurf,Nctlu,Nctlv,3))
    for isurf in xrange(Nsurf):
        for idim in xrange(3):
            ctl[isurf,:,:,idim] = reshape(x[isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv : isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv],[Nctlu,Nctlv])
        # end if
    #end if

    for isurf in xrange(Nsurf):
        for idim in xrange(3):
            temp = 2*dot(dot(J[isurf],ctl[isurf,:,:,idim].flatten())-X[isurf,:,:,idim].flatten(),J[isurf])

            g_obj[isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv : isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv] = temp
        # end for
    # end for 
    
    fail = False
    g_con = B
    h = 1.0e-40j
    x = array(x,'D')
    for i in xrange(ndv):
        x[i] += h

        ctl = zeros((Nsurf,Nctlu,Nctlv,3),'D')
        for isurf in xrange(Nsurf):
            for idim in xrange(3):
                ctl[isurf,:,:,idim] = reshape(x[isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv : isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv],[Nctlu,Nctlv])
                
        #  # Calculate the LE constraint
        xA = ctl[0,0,0,:] # Root LE (upper)
        xB = ctl[0,1,0,:] # Root LE upper +1
        xC = ctl[1,-2,0,:]# Root LE lower +1
        
        v1 = xB-xA #vector from A to B
        v2 = xA-xC #vector from A to C
        xProduct = (v1[1]*v2[2]-v1[2]*v2[1])**2+ (-v1[0]*v2[2]+v1[2]*v2[0] )**2 +  ( v1[0]*v2[1] - v1[1]*v2[0])**2
#         # Now take cross product and quadrature sum
#        g_con[-1,i] = imag((1.0e6*xProduct))/imag(h)
        g_con[-2,i] = imag((v1[0]/v2[0])/(v1[1]/v2[1]))/imag(h)
        g_con[-1,i] = imag((v1[0]/v2[0])/(v1[2]/v2[2]))/imag(h)
        x[i] -= h
#     #end for

    timeCounter += time.time()-timeA
    return g_obj,g_con,fail

# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

# ================================================
# Find a good guess of the initial control points 
# ================================================

# ctl = zeros((Nsurf,Nctlu,Nctlv,3))

# #Create the interpolation
# u_interp = 0.5*(1-cos(linspace(0,pi,Nctlu)))
# v_interp = linspace(0,1,Nctlv)

# for isurf in xrange(Nsurf):
#     for idim in xrange(3):
#         I = RectBivariateSpline(u[isurf],v[isurf],X[isurf,:,:,idim],kx=1,ky=1)
#         for i in xrange(Nctlu):
#             for j in xrange(Nctlv):
#                 ctl[isurf,i,j,idim] = I(u_interp[i],v_interp[j])
#             # end for
#         # end for
#     # end for
# # end for 
# ndv = 0
# for isurf in xrange(Nsurf):
#     for idim in xrange(3):
#         if idim == 0: name = 'ctlx'
#         if idim == 1: name = 'ctly'
#         if idim == 2: name = 'ctlz'
#         name +=str(isurf)
#         opt_prob.addVarGroup(name,Nctlu*Nctlv,'c',value=ctl[isurf,:,:,idim].flatten(),lower=-100,upper=100)
#         ndv += Nctlu*Nctlv
#     # end for
# # end for 



ndv = Nctlu*Nctlv*Nsurf*3


# ===================
#  Constraints
# ===================
# First figure out how many constraints we are going to have
ncon = 0

# Each of four corners on each Surf

ncon += 4*Nsurf*3 # three for the dimensions
ncon += 2*Nctlv*3 # Set the LE and TE to be same on upper and lower surfaces
ncon += 2         # one constraint for the continutity constraint
global B
B = zeros([ncon,ndv]) #This is the constraint derivative matrix (constant since these are linear constraints)

# Corner Constraints
counter = 0
for isurf in xrange(Nsurf):
    for idim in xrange(3):
        opt_prob.addCon('coner constr',type= 'i',lower=X[isurf,0,0,idim],upper=X[isurf,0,0,idim])
        opt_prob.addCon('coner constr',type= 'i',lower=X[isurf,0,-1,idim],upper=X[isurf,0,-1,idim])
        opt_prob.addCon('coner constr',type= 'i',lower=X[isurf,-1,0,idim],upper=X[isurf,-1,0,idim])
        opt_prob.addCon('coner constr',type= 'i',lower=X[isurf,-1,-1,idim],upper=X[isurf,-1,-1,idim])

        B[counter  ,isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + 0] = 1
        B[counter+1,isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlv-1] = 1
        B[counter+2,isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv-Nctlv] = 1
        B[counter+3,isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv-1] = 1
        counter += 4
    # end for
# end for 

# Edge Constraints
for idim in xrange(3):
    for j in xrange(Nctlv):
        opt_prob.addCon('edge constraint',type='i',lower=0,upper=0)
        opt_prob.addCon('edge constraint',type='i',lower=0,upper=0)
        B[counter,idim*Nctlv*Nctlu + j] = 1
        B[counter,3*Nctlu*Nctlv + idim*Nctlv*Nctlu + (Nctlu-1)*Nctlv + j] = -1
        B[counter+1,idim*Nctlv*Nctlu + (Nctlu-1)*Nctlv + j] = 1
        B[counter+1,3*Nctlu*Nctlv + idim*Nctlv*Nctlu + j] = -1
        counter +=2
    # end for
# end for
# LE Continutiy Constraint

opt_prob.addCon('LE constraint',type='i',lower=1,upper=1)
opt_prob.addCon('LE constraint',type='i',lower=1,upper=1)

# Lets do a lms fit instead
timeA = time.time()
ctl = pyspline.fit_surf(Nsurf,Nu,Nv,Nctlu,Nctlv,ncon,J,X,B,zeros(ncon))
print 'LMS Fit Time:',time.time()-timeA
for isurf in xrange(Nsurf):
    for idim in xrange(3):
        if idim == 0: name = 'ctlx'
        if idim == 1: name = 'ctly'
        if idim == 2: name = 'ctlz'
        name +=str(isurf)
        opt_prob.addVarGroup(name,Nctlu*Nctlv,'c',value=ctl[isurf,:,:,idim].flatten(),lower=-100,upper=100)
    # end for
# end for 

# ===================
#  Objective
# ===================
opt_prob.addObj('RMSerror')

# ===================
#  Select Optimizer
# ===================

opt = SNOPT()

# ===================
#  SNOPT Options
# ===================
#opt.setOption('Derivative level',0)
#opt.setOption('Verify level',3)
opt.setOption('Major iterations limit',150)
#opt.setOption('Linesearch tolerance',.1)
#opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-5)
opt.setOption('Major feasibility tolerance',1e-5)
opt.setOption('Minor feasibility tolerance',1e-5)

# ===================
#  Run Optimization
# ===================

result = opt(opt_prob,sens)

print opt_prob._solutions[0]
print '#--------------------------------'
print '# RMS Error: ',sqrt(result[0][0]/(2*Nu*Nv))
print '#--------------------------------'


# ===================
#  Print Solution  
# ===================
x = result[1][:]
# Unpack the x-values
ctl = zeros((Nsurf,Nctlu,Nctlv,3))
for isurf in xrange(Nsurf):
    for idim in xrange(3):
        ctl[isurf,:,:,idim] = reshape(x[isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv : isurf*3*Nctlu*Nctlv + idim*Nctlu*Nctlv + Nctlu*Nctlv],[Nctlu,Nctlv])
    # end if
# end if

# Dump re-interpolated surface
for isurf in xrange(Nsurf):
    f.write('Zone I=%d J = %d\n'%(Nu,Nv))
    f.write('DATAPACKING=POINT\n')
    for j in xrange(Nv):
        for i in xrange(Nu):
            f.write('%f %f %f \n'%(pyspline.b2val(u[isurf,i],v[isurf,j],0,0,tu0,tv0,ku,kv,ctl[isurf,:,:,0]),
                                   pyspline.b2val(u[isurf,i],v[isurf,j],0,0,tu0,tv0,ku,kv,ctl[isurf,:,:,1]),
                                   pyspline.b2val(u[isurf,i],v[isurf,j],0,0,tu0,tv0,ku,kv,ctl[isurf,:,:,2])))

# Dump Control Points
for isurf in xrange(Nsurf):
    f.write('Zone I=%d J = %d\n'%(Nctlu,Nctlv))
    f.write('DATAPACKING=POINT\n')
    for j in xrange(Nctlv):
        for i in xrange(Nctlu):
            f.write('%f %f %f \n'%(ctl[isurf,i,j,0],ctl[isurf,i,j,1],ctl[isurf,i,j,2]))



print 'Eval Time:',timeCounter
