# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,hstack,zeros,ones,sqrt,imag,interp,array,real
from scipy.io import read_array

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

global USE_CTLPTS
global USE_KNOTS
global USE_PARAMS

USE_CTLPTS = True
USE_KNOTS  = False
USE_PARAMS = False

Nctlu = 9
Nctlv = 9
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
        xu0[i,j] = data[Nu*j+i,0]
        yu0[i,j] = data[Nu*j+i,1]
        zu0[i,j] = data[Nu*j+i,2]
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
        xl0[i,j] = data[Nu*j+i,0]
        yl0[i,j] = data[Nu*j+i,1]
        zl0[i,j] = data[Nu*j+i,2]
    #end for
#end for

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
tu0[ku-1:Nctlu+1] = linspace(0,pi,Nctlu-ku-2)

# V Knots

tv0= zeros(Nctlv + kv)
tv0[0:kv] = 0
tv0[Nctlv:Nctlv+kv] = 1
tv0[kv-1:Nctlv+1]=  0.5*(1-cos(linspace(0,pi,Nctlv-kv+2)))


# ==============================
# Find the Control Point Vectors
# ==============================

# Initial Control Pt Vector
ctlxu0 = zeros([Nu,Nv])
ctlyu0 = zeros([Nu,Nv])
ctlzu0 = zeros([Nu,Nv])

ctlxl0 = zeros([Nu,Nv])
ctlyl0 = zeros([Nu,Nv])
ctlzl0 = zeros([Nu,Nv])


print 'Check of inputs:'

print 'uu,ul',uu,ul
print 'vu,vl',vu,vl
print 'tu0,tv0',tu0,tv0

def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global USE_CTLPTS
    global USE_KNOTS
    global USE_PARAMS
    
    dv_counter = 0

    if USE_CTLPTS:
#         coefxu = zeros(Nu,Nv)
#         coefyu = zeros(Nu,Nv)
#         coefzu = zeros(Nu,Nv)
#         coefxl = zeros(Nu,Nv)
#         coefyl = zeros(Nu,Nv)
#         coefzl = zeros(Nu,Nv)

        coefxu = reshape(x[0      :  Nu*NV],[Nu,Nv])
        coefyu = reshape(x[1*Nu*Nv:2*Nu*NV],[Nu,Nv])
        coefzu = reshape(x[2*Nu*Nv:3*Nu*NV],[Nu,Nv])
        coefxl = reshape(x[3*Nu*Nv:4*Nu*NV],[Nu,Nv])
        coefyl = reshape(x[4*Nu*Nv:5:Nu*Nv],[Nu,Nv])
        coefzl = reshape(x[5*Nu*Nv:6*Nu*NV],[Nu,Nv])

        dv_counter += 6*(Nu*Nv)
    else:
        coefxu = coefxu0
        coefyu = coefyu0
        coefzu = coefzu0
        coefxl = coefxl0
        coefyl = coefyl0
        coefzl = coefzl0

        
    if USE_KNOTS:
        t = zeros(Nctl+k)
        t[-4:] = 1
        t[k:Nctl] = x[dv_counter:dv_counter + Nctl-k]

        dv_counter += Nctl-k
    else:
        t = t0

    if USE_PARAMS:
        s = x[dv_counter:dv_counter + N]
        dv_counter += N
    else:
        s = s0

    err = 0.0;

    for i in xrange(u):
        for j in xrange(v):
             x_interp_u = pyspline.b2val(uu[i],vu[j],0,0,tu,tv,ku,kv,coefxu)
             y_interp_u = pyspline.b2val(uu[i],vu[j],0,0,tu,tv,ku,kv,coefyu)
             z_interp_u = pyspline.b2val(uu[i],vu[j],0,0,tu,tv,ku,kv,coefzu)

             x_interp_u = pyspline.b2val(ul[i],vl[j],0,0,tu,tv,ku,kv,coefxl)
             y_interp_u = pyspline.b2val(ul[i],vl[j],0,0,tu,tv,ku,kv,coefyl)
             z_interp_u = pyspline.b2val(ul[i],vl[j],0,0,tu,tv,ku,kv,coefzl)
             
             err += (x_interp_u-xu0[i,j])**2 + (x_interp_l-xl0[i,j])**2 \
                    (y_interp_u-yu0[i,j])**2 + (y_interp_l-yl0[i,j])**2 \
                    (z_interp_u-zu0[i,j])**2 + (z_interp_l-zl0[i,j])**2
        # end for
    # end for

    f = sqrt(err/N)

    fcon = array([])
            
    if USE_KNOTS:
        fcon = hstack([fcon,zeros(Nctl-k)])
        
        for i in xrange(Nctl-k):
            fcon[i] = t[i+k+1]-t[i+k]
        # end for
    # end if 
    fail = False
    return f,fcon,fail

# def sens(x,f_obj,f_con):
#     global USE_CTLPTS
#     global USE_KNOTS
#     global USE_PARAMS
#     use_params_start = 0

#     ndv = len(x)

#     ncon = 0
#     if USE_KNOTS:
#         ncon += Nctl-k
#     ncon += 1 # derivative constraint

#     g_obj = zeros(ndv) #objective derivative vector
#     g_con = zeros([ncon, ndv]) #constraint jacobian
    
#     h = 1.0e-40j #complex step size
    
#     for j in xrange(ndv):

#         x_deriv = zeros([ndv],'complex')
        
#         for i in range(ndv): #copy over design variables
#             x_deriv[i] = complex(x[i],0)

#         x_deriv[j] += h
     
#         dv_counter = 0
    
#         if USE_CTLPTS:
#             ctlx = zeros(Nctl,'complex')
#             ctly = zeros(Nctl,'complex')

#             ctlx[1:Nctl-1] = x_deriv[0:Nctl-2] # These are the control points
#             ctly[1:Nctl-1] = x_deriv[Nctl-2:2*Nctl-4]

#             ctlx[0]  = x0[0]
#             ctlx[-1] = x0[-1]
#             ctly[0]  = y0[0]
#             ctly[-1] = y0[-1]

#             dv_counter += 2*Nctl-4
#         else:
#             ctlx = ctlx0
#             ctly = ctly0

#         if USE_KNOTS:
#             t = zeros(Nctl+k,'complex')
#             t[-4:] = 1
#             t[k:Nctl] = x_deriv[dv_counter:dv_counter + Nctl-k]

#             dv_counter += Nctl-k
#         else:
#             t = t0

#         if USE_PARAMS:
#             use_params_start = dv_counter;
#             s = x_deriv[dv_counter:dv_counter + N]
#             dv_counter += N
#         else:
#             s = s0

#         err = 0.0;

#         if not USE_PARAMS or j<use_params_start: #Do derivatives the 'dumb' way
#             interp_x = pyspline_cs.bvaluv(t,ctlx,k,0,s)
#             interp_y = pyspline_cs.bvaluv(t,ctly,k,0,s)
#             err = sum((interp_x-x0)**2 + (interp_y-y0)**2)
#             f = sqrt(err/N)
#             g_obj[j] = imag(f)/imag(h)
#         else: #Do the s derivatives in a more elegent fashion
#             index = j-use_params_start
#             interp_x = pyspline_cs.bvalu(t,ctlx,k,0,s[index])
#             interp_y = pyspline_cs.bvalu(t,ctly,k,0,s[index])
#             temp = (interp_x-x0[index])**2  + (interp_y-y0[index])**2 
#             err = f_obj**2*N + temp - real(temp)
#             f = sqrt(err/N)
#             g_obj[j] = imag(f)/imag(h)
#         # end if

        
#         fcon = array([],'complex')
            
#         if USE_KNOTS:
#             fcon = hstack([fcon,zeros(Nctl-k,'complex')])
            
#             for i in xrange(Nctl-k):
#                 fcon[i] = t[i+k+1]-t[i+k]
#             # end for
#         # end if 
            
#         # Derivative constraint
#         dx = pyspline_cs.bvalu(t,ctlx,k,1,0)
#         fcon = hstack([fcon,dx])
            
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

#Make sure the end points are what we want:

if USE_CTLPTS:
    opt_prob.addVarGroup('CTLux',(Nctlu*Nctlv),'c',value=ctlxu0,lower=-100,upper=100)
    opt_prob.addVarGroup('CTLuy',(Nctlu*Nctlv),'c',value=ctlyu0,lower=-100,upper=100)
    opt_prob.addVarGroup('CTLuz',(Nctlu*Nctlv),'c',value=ctlzu0,lower=-100,upper=100)
    opt_prob.addVarGroup('CTLlx',(Nctlu*Nctlv),'c',value=ctlxl0,lower=-100,upper=100)
    opt_prob.addVarGroup('CTLly',(Nctlu*Nctlv),'c',value=ctlyl0,lower=-100,upper=100)
    opt_prob.addVarGroup('CTLlz',(Nctlu*Nctlv),'c',value=ctlzl0,lower=-100,upper=100)

# Knots are not implemented yet
if USE_KNOTS:
    opt_prob.addVarGroup('Knot',Nctl-k,'c',value=t0[k:Nctl],lower=0,upper=1)
if USE_PARAMS:
    opt_prob.addVarGroup('s:',N,'c',value=s0,lower=0,upper=1)

# ===================
#  Constraints
# ===================
# Knots are not implemented yet
if USE_KNOTS:
    opt_prob.addConGroup('Knots', Nctl-k,type='i', lower=0.001, upper=.33)
#opt_prob.addCon('LEconstraint',type = 'i',lower=0.0,upper=0.0)

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
opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-6)
opt.setOption('Major feasibility tolerance',1e-6)
opt.setOption('Minor feasibility tolerance',1e-6)

# ===================
#  Run Optimization
# ===================

opt(opt_prob)#,sens)

# ===================
#  Print Solution  
# ===================

print opt_prob._solutions[0]
