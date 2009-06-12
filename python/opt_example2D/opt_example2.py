# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot
from scipy.io import read_array
import scipy.interpolate

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

Nctlu = 11
Nctlv = 7
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
        xu0[i,j] = data[Nv*i+j,0]
        yu0[i,j] = data[Nv*i+j,1]
        zu0[i,j] = data[Nv*i+j,2]
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
        xl0[i,j] = data[Nv*i+j,0]
        yl0[i,j] = data[Nv*i+j,1]
        zl0[i,j] = data[Nv*i+j,2]
    #end for
#end for

# Check the Data

f = open('output.dat','w')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for i in xrange(Nu):
    for j in xrange(Nv):
        f.write('%f %f %f \n'%(xu0[i,j],yu0[i,j],zu0[i,j]))
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for i in xrange(Nu):
    for j in xrange(Nv):
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
x0=xl0[0,:];y0=yl0[0,:];z0=zl0[0,:] 
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
tu0[0:ku] = 0
tu0[Nctlu:Nctlu+ku] = 1
tu0[ku-1:Nctlu+1]  = 0.5*(1-cos(linspace(0,pi,Nctlu-ku+2)))

# V Knots

tv0= zeros(Nctlv + kv)
tv0[0:kv] = 0
tv0[Nctlv:Nctlv+kv] = 1
tv0[kv-1:Nctlv+1]=  linspace(0,1,Nctlv-kv+2)


# Calculate the jacobian J, for fixed t and s
h = 1.0e-40j
Ju = zeros([Nu*Nv,Nctlu*Nctlv])
Jl = zeros([Nu*Nv,Nctlu*Nctlv])
ctl = zeros([Nctlu,Nctlv],'D')
[Uu,Vu] = meshgrid(uu,vu)
[Ul,Vl] = meshgrid(ul,vl)

for i in xrange(Nctlu):
    for j in xrange(Nctlv):
        ctl[i,j] += h
        valu = pyspline_cs.b2valv(Uu.flatten(),Vu.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        vall = pyspline_cs.b2valv(Ul.flatten(),Vl.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        
        ctl[i,j] -= h    
        Ju[:,j*Nctlu + i] = imag(valu)/imag(h)
        Jl[:,j*Nctlu + i] = imag(vall)/imag(h)
    # end for 
# end for




def objcon(x):
    '''Get the rms error for the given set of design variables'''

#     coefxu = zeros([Nu,Nv])
#     coefyu = zeros([Nu,Nv])
#     coefzu = zeros([Nu,Nv])
#     coefxl = zeros([Nu,Nv])
#     coefyl = zeros([Nu,Nv])
#     coefzl = zeros([Nu,Nv])

#     coefxu = reshape(x[0        :   Nu*Nv- 4],[Nu,Nv])
#     coefyu = reshape(x[1*Nu*Nv-4: 2*Nu*Nv- 8],[Nu,Nv])
#     coefzu = reshape(x[2*Nu*Nv-8: 3*Nu*Nv-12],[Nu,Nv])
#     coefxl = reshape(x[3*Nu*Nv-12:4*Nu*Nv-16],[Nu,Nv])
#     coefyl = reshape(x[4*Nu*Nv-16:5:Nu*Nv-20],[Nu,Nv])
#     coefzl = reshape(x[5*Nu*Nv-20:6*Nu*Nv-24],[Nu,Nv])

#     # Fix all the corners of the grid mesh
#     # Upper
#     ctlxu[0 , 0] = x0u[0 , 0]
#     ctlxu[0 ,-1] = x0u[0 ,-1]
#     ctlxu[-1, 0] = x0u[-1, 0]
#     ctlxu[-1,-1] = x0u[-1,-1]

#     ctlyu[0 , 0] = y0u[0 , 0]
#     ctlyu[0 ,-1] = y0u[0 ,-1]
#     ctlyu[-1, 0] = y0u[-1, 0]
#     ctlyu[-1,-1] = y0u[-1,-1]

#     ctlzu[0 , 0] = z0u[0 , 0]
#     ctlzu[0 ,-1] = z0u[0 ,-1]
#     ctlzu[-1, 0] = z0u[-1, 0]
#     ctlzu[-1,-1] = z0u[-1,-1]

#     # Lower
#     ctlxl[0 , 0] = x0l[0 , 0]
#     ctlxl[0 ,-1] = x0l[0 ,-1]
#     ctlxl[-1, 0] = x0l[-1, 0]
#     ctlxl[-1,-1] = x0l[-1,-1]

#     ctlyl[0 , 0] = y0l[0 , 0]
#     ctlyl[0 ,-1] = y0l[0 ,-1]
#     ctlyl[-1, 0] = y0l[-1, 0]
#     ctlyl[-1,-1] = y0l[-1,-1]
    
#     ctlzl[0 , 0] = z0l[0 , 0]
#     ctlzl[0 ,-1] = z0l[0 ,-1]
#     ctlzl[-1, 0] = z0u[-1, 0]
#     ctlzl[-1,-1] = z0l[-1,-1]

    ctlxu = reshape(x[0             : Nctlu*Nctlv],[Nctlu,Nctlv])
    ctlyu = reshape(x[1*Nctlu*Nctlv:2*Nctlu*Nctlv],[Nctlu,Nctlv])
    ctlzu = reshape(x[2*Nctlu*Nctlv:3*Nctlu*Nctlv],[Nctlu,Nctlv])
#     ctlxl = reshape(x[3*Nctlu*Nctlv:4*Nctlu*Nctlv],[Nctlu,Nctlv])
#     ctlyl = reshape(x[4*Nctlu*Nctlv:5*Nctlu*Nctlv],[Nctlu,Nctlv])
#     ctlzl = reshape(x[5*Nctlu*Nctlv:6*Nctlu*Nctlv],[Nctlu,Nctlv])

#     print 'in objcon'
#     print ctlxu
#     print ctlyu
#     print ctlzu

    sum1 = (dot(Ju,ctlxu.flatten())-xu0.flatten())**2
    sum2 = (dot(Ju,ctlyu.flatten())-yu0.flatten())**2
    sum3 = (dot(Ju,ctlzu.flatten())-zu0.flatten())**2

    sum1 = xu0.flatten()**2
    sum2 = yu0.flatten()**2
    sum3 = zu0.flatten()**2

    print dot(sum1,sum1),dot(sum2,sum2),dot(sum3,sum3)
   
#     f = sum( (    dot(Ju,ctlxu.flatten())-xu0.flatten())**2 + \
#                  (dot(Ju,ctlyu.flatten())-yu0.flatten())**2 + \
#                  (dot(Ju,ctlzu.flatten())-zu0.flatten())**2)
# #                  (dot(Jl,ctlxl.flatten())-xl0.flatten())**2 + \
# #                  (dot(Jl,ctlyl.flatten())-yl0.flatten())**2 + \
# #                  (dot(Jl,ctlzl.flatten())-zl0.flatten())**2   )    
     
    #print 'sums:',sum1,sum2,sum3
    f= sum(sum1+sum2+sum3)
    fcon = array([])
    
    fail = False
    return f,fcon,fail

# def sens(x,f_obj,f_con):
            
#         g_con[:,j] = imag(fcon)/imag(h)
        
#     # end for (dv loop)
#     fail = False
#     return g_obj,g_con,fail

# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

# ==============================
# Find the Control Point Vectors
# ==============================


opt_prob.addVarGroup('CTLux',(Nctlu*Nctlv),'c',value=0,lower=-100,upper=100)
opt_prob.addVarGroup('CTLuy',(Nctlu*Nctlv),'c',value=0,lower=-100,upper=100)
opt_prob.addVarGroup('CTLuz',(Nctlu*Nctlv),'c',value=0,lower=-100,upper=100)
# opt_prob.addVarGroup('CTLlx',(Nctlu*Nctlv),'c',value=ctlxl0.flatten(),lower=-100,upper=100)
# opt_prob.addVarGroup('CTLly',(Nctlu*Nctlv),'c',value=ctlyl0.flatten(),lower=-100,upper=100)
# opt_prob.addVarGroup('CTLlz',(Nctlu*Nctlv),'c',value=ctlzl0.flatten(),lower=-100,upper=100)

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
opt.setOption('Derivative level',0)
opt.setOption('Major iterations limit',1)
opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-6)
opt.setOption('Major feasibility tolerance',1e-6)
opt.setOption('Minor feasibility tolerance',1e-6)

# ===================
#  Run Optimization
# ===================

result = opt(opt_prob)#,sens)

# ===================
#  Print Solution  
# ===================

ctlxu = reshape(result[1][0             : Nctlu*Nctlv],[Nctlu,Nctlv])
ctlyu = reshape(result[1][1*Nctlu*Nctlv:2*Nctlu*Nctlv],[Nctlu,Nctlv])
ctlzu = reshape(result[1][2*Nctlu*Nctlv:3*Nctlu*Nctlv],[Nctlu,Nctlv])
print 'ctlxu:',ctlxu
print 'ctlyu:',ctlyu
print 'ctlzu:',ctlzu

f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for i in xrange(Nu):
    for j in xrange(Nv):
        f.write('%f %f %f \n'%(pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlxu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlyu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlzu)))
