# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot
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

Nctlu =10
Nctlv =10
ku=4
kv=4

# ==============================
# Read in the Data 
# ==============================

# Top surface of wing
print 'reading upper surface...'
data = read_array(file("upper_surface.inp"))

Nu = 26
Nv = 25

xu0 = zeros((Nu,Nv))
yu0 = zeros((Nu,Nv))
zu0 = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        xu0[i,j] = data[j*Nu + i ,0]
        yu0[i,j] = data[j*Nu + i ,1]
        zu0[i,j] = data[j*Nu + i ,2]
    #end for
#end for

# Bottom surface of wing
print 'reading bottom surface...'
data = read_array(file("lower_surface.inp"))

xl0 = zeros((Nu,Nv))
yl0 = zeros((Nu,Nv))
zl0 = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        xl0[i,j] = data[j*Nu + i,0]
        yl0[i,j] = data[j*Nu + i,1]
        zl0[i,j] = data[j*Nu + i,2]
    #end for
#end for

#Flatten Data Once
xu0flat = xu0.flatten()
yu0flat = yu0.flatten()
zu0flat = zu0.flatten()
xl0flat = xl0.flatten()
yl0flat = yl0.flatten()
zl0flat = zl0.flatten()

# Check the Data

f = open('output.dat','w')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(xu0[i,j],yu0[i,j],zu0[i,j]))
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(xl0[i,j],yl0[i,j],zl0[i,j]))

# ==============================
# Find the 'u' parameters
# ==============================
x0=xu0[:,0];y0=yu0[:,0];z0=zu0[:,0] 
uu=zeros(Nu)
for i in xrange(Nu-1):
     uu[i+1] = uu[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
uu = uu/uu[-1] #Normalize uu

x0=xu0[:,0];y0=yu0[:,0];z0=zu0[:,0] 
ul=zeros(Nu)
for i in xrange(Nu-1):
     ul[i+1] = ul[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
ul = ul/ul[-1] #Normalize ul

# ==============================
# Find the 'v' parameters
# ==============================
x0=xu0[0,:];y0=yu0[0,:];z0=zu0[0,:] 
vu = zeros(Nv)
for i in xrange(Nv-1):
    vu[i+1] = vu[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
      
                         (z0[i+1]-z0[i])**2)
#end for
vu = vu/vu[-1] #Normalize vu

x0=xl0[0,:];y0=yl0[0,:];z0=zl0[0,:] 
vl = zeros(Nv)
for i in xrange(Nv-1):
    vl[i+1] = vl[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
vl = vl/vl[-1] #Normalize vl

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
Ju = zeros([Nu*Nv,Nctlu*Nctlv])
Jl = zeros([Nu*Nv,Nctlu*Nctlv])
ctl = zeros([Nctlu,Nctlv],'D')

[VU, UU] = meshgrid(vu,uu)
[VL, UL] = meshgrid(vu,uu)

for j in xrange(Nctlv):
    for i in xrange(Nctlu):
        ctl[i,j] += h
        valu = pyspline_cs.b2valv(UU.flatten(),VU.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        vall = pyspline_cs.b2valv(UL.flatten(),VL.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        ctl[i,j] -= h    
        Ju[:,i*Nctlv + j] = imag(valu)/imag(h)
        Jl[:,i*Nctlv + j] = imag(vall)/imag(h)
    # end for 
# end for
global timeCounter 
timeCounter = 0
def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global timeCounter
    timeA = time.time()
 
    ctlxu = x[0             : Nctlu*Nctlv]
    ctlyu = x[1*Nctlu*Nctlv:2*Nctlu*Nctlv]
    ctlzu = x[2*Nctlu*Nctlv:3*Nctlu*Nctlv]
    ctlxl = x[3*Nctlu*Nctlv:4*Nctlu*Nctlv]
    ctlyl = x[4*Nctlu*Nctlv:5*Nctlu*Nctlv]
    ctlzl = x[5*Nctlu*Nctlv:6*Nctlu*Nctlv]

    sum1 = sum((dot(Ju,ctlxu)-xu0flat)**2)
    sum2 = sum((dot(Ju,ctlyu)-yu0flat)**2)
    sum3 = sum((dot(Ju,ctlzu)-zu0flat)**2)
    sum4 = sum((dot(Jl,ctlxl)-xl0flat)**2)
    sum5 = sum((dot(Jl,ctlyl)-yl0flat)**2)
    sum6 = sum((dot(Jl,ctlzl)-zl0flat)**2)

    f= sum1+sum2+sum3+sum4+sum5+sum6
    fcon = array([])
    #rint 'obj_time:',time.time()-timeA
    fail = False
    timeCounter += time.time()-timeA
    return f,fcon,fail

def sens(x,f_obj,f_con):
    global timeCounter
    timeA = time.time()
    ndv = len(x)
    g_obj = zeros(ndv)
   
    ctlxu = x[0             : Nctlu*Nctlv]
    ctlyu = x[1*Nctlu*Nctlv:2*Nctlu*Nctlv]
    ctlzu = x[2*Nctlu*Nctlv:3*Nctlu*Nctlv]
    ctlxl = x[3*Nctlu*Nctlv:4*Nctlu*Nctlv]
    ctlyl = x[4*Nctlu*Nctlv:5*Nctlu*Nctlv]
    ctlzl = x[5*Nctlu*Nctlv:6*Nctlu*Nctlv]
   
    g_obj[            0:  Nctlu*Nctlv] = 2*dot(dot(Ju,ctlxu)-xu0flat,Ju)
    g_obj[  Nctlu*Nctlv:2*Nctlu*Nctlv] = 2*dot(dot(Ju,ctlyu)-yu0flat,Ju)
    g_obj[2*Nctlu*Nctlv:3*Nctlu*Nctlv] = 2*dot(dot(Ju,ctlzu)-zu0flat,Ju)
    g_obj[3*Nctlu*Nctlv:4*Nctlu*Nctlv] = 2*dot(dot(Jl,ctlxl)-xl0flat,Jl)
    g_obj[4*Nctlu*Nctlv:5*Nctlu*Nctlv] = 2*dot(dot(Jl,ctlyl)-yl0flat,Jl)
    g_obj[5*Nctlu*Nctlv:6*Nctlu*Nctlv] = 2*dot(dot(Jl,ctlzl)-zl0flat,Jl)

    g_con = array([])
    fail = False
    #rint 'constr time:',time.time()-timeA
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

ctlxu = zeros((Nctlu,Nctlv))
ctlyu = zeros((Nctlu,Nctlv))
ctlzu = zeros((Nctlu,Nctlv))
ctlxl = zeros((Nctlu,Nctlv))
ctlyl = zeros((Nctlu,Nctlv))
ctlzl = zeros((Nctlu,Nctlv))

#Create the interpolation

Ixu = RectBivariateSpline(uu,vu,xu0,kx=1,ky=1)
Iyu = RectBivariateSpline(uu,vu,yu0,kx=1,ky=1)
Izu = RectBivariateSpline(uu,vu,zu0,kx=1,ky=1)
Ixl = RectBivariateSpline(ul,vl,xl0,kx=1,ky=1)
Iyl = RectBivariateSpline(ul,vl,yl0,kx=1,ky=1)
Izl = RectBivariateSpline(ul,vl,zl0,kx=1,ky=1)

u_interp = 0.5*(1-cos(linspace(0,pi,Nctlu)))
v_interp = linspace(0,1,Nctlv)
for i in xrange(Nctlu):
    for j in xrange(Nctlv):
        ctlxu[i,j] = Ixu(u_interp[i],v_interp[j])
        ctlyu[i,j] = Iyu(u_interp[i],v_interp[j]) 
        ctlzu[i,j] = Izu(u_interp[i],v_interp[j])
        ctlxl[i,j] = Ixl(u_interp[i],v_interp[j])
        ctlyl[i,j] = Iyl(u_interp[i],v_interp[j]) 
        ctlzl[i,j] = Izl(u_interp[i],v_interp[j])
   
opt_prob.addVarGroup('CTLxu',(Nctlu*Nctlv),'c',value=ctlxu.flatten(),lower=-100,upper=100)
opt_prob.addVarGroup('CTLyu',(Nctlu*Nctlv),'c',value=ctlyu.flatten(),lower=-100,upper=100)
opt_prob.addVarGroup('CTLzu',(Nctlu*Nctlv),'c',value=ctlzu.flatten(),lower=-100,upper=100)
opt_prob.addVarGroup('CTLxl',(Nctlu*Nctlv),'c',value=ctlxl.flatten(),lower=-100,upper=100)
opt_prob.addVarGroup('CTLyl',(Nctlu*Nctlv),'c',value=ctlyl.flatten(),lower=-100,upper=100)
opt_prob.addVarGroup('CTLzl',(Nctlu*Nctlv),'c',value=ctlzl.flatten(),lower=-100,upper=100)

# ===================
#  Constraints
# ===================
# Knots are not implemented yet

#opt_prob.addConGroup('LEconstraints',Nu,type = 'i',lower=0.0,upper=0.0)

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
opt.setOption('Major iterations limit',500)
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

ctlxu = reshape(result[1][0             : Nctlu*Nctlv],[Nctlu,Nctlv])
ctlyu = reshape(result[1][1*Nctlu*Nctlv:2*Nctlu*Nctlv],[Nctlu,Nctlv])
ctlzu = reshape(result[1][2*Nctlu*Nctlv:3*Nctlu*Nctlv],[Nctlu,Nctlv])
ctlxl = reshape(result[1][3*Nctlu*Nctlv:4*Nctlu*Nctlv],[Nctlu,Nctlv])
ctlyl = reshape(result[1][4*Nctlu*Nctlv:5*Nctlu*Nctlv],[Nctlu,Nctlv])
ctlzl = reshape(result[1][5*Nctlu*Nctlv:6*Nctlu*Nctlv],[Nctlu,Nctlv])

print 'ctlxu:',ctlxu
print 'ctlyu:',ctlyu
print 'ctlzu:',ctlzu

#Just call the objective

# x = hstack([ctlxu.flatten(),ctlyu.flatten(),ctlzu.flatten()])
# objcon(x)


f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlxu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlyu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlzu)))



f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlxl),
                               pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlyl),
                               pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlzl)))



print 'Eval Time:',timeCounter
