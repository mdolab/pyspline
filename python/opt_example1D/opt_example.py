# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,hstack,zeros,ones,sqrt,imag,interp,array,real
from matplotlib.pylab import plot,show

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
airfoil coordinates, specificially the upper surface of a naca 0012
airfoil with a constraint on a vertical leading edge'''


# Global Variables

global ctlx
global ctly

global USE_CTLPTS
global USE_KNOTS
global USE_PARAMS


USE_CTLPTS = True
USE_KNOTS  = True
USE_PARAMS = False

Nctl = 13
k = 4

# Top surface of a naca0012 airfoil
N = 61
x0 = 0.5*(1-cos(linspace(0,pi,N)))
y0= 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)

# Top surface of a WT airfoil
# N = 34
# x0 = array([0.00000,0.00018,0.00255,0.00954,0.02090,0.03650,0.05640,0.08030,0.10801,0.13934,0.17395,0.21146,0.25149,0.29361,0.33736,0.38228,0.42820,0.47526,0.52324,0.57161,0.61980,0.66724,0.71333,0.75749,0.79915,0.83778,0.87287,0.90391,0.93072,0.95355,0.97251,0.98719,0.99668,1.00000])
# y0 = array([0.00000,0.00159,0.00748,0.01640,0.02600,0.03580,0.04560,0.05520,0.06430,0.07290,0.08070,0.08760,0.09340,0.09810,0.10133,0.10294,0.10249,0.10005,0.09610,0.09090,0.08490,0.07820,0.07100,0.06340,0.05570,0.04800,0.04030,0.03260,0.02480,0.01700,0.00982,0.00431,0.00103,0.00000])


# S-Parameter Calculation
s0 = zeros(N)
for i in xrange(N-1):
    s0[i+1] = s0[i] + sqrt((x0[i+1]-x0[i])**2 + (y0[i+1]-y0[i])**2)

s0 = s0/s0[-1] #Normalize s
print 's0:',s0

# Initial Knot Vector Calculation
t0= zeros(Nctl + k)
t0[0:k] = s0[0]
t0[Nctl:Nctl+k] = s0[-1]
t0[k-1:Nctl+1]= 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))

# Initial Control Pt Vector
ctlx0 = 0.5*(1-cos(linspace(0,pi,Nctl)))
ctly0 = interp(ctlx0,x0,y0)

def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global USE_CTLPTS
    global USE_KNOTS
    global USE_PARAMS
    global ctlx
    global ctly

    dv_counter = 0

    if USE_CTLPTS:
        ctlx = zeros(Nctl,'complex')
        ctly = zeros(Nctl,'complex')

        ctlx[1:Nctl-1] = x[0:Nctl-2] # These are the control points
        ctly[1:Nctl-1] = x[Nctl-2:2*Nctl-4]

        ctlx[0]  = x0[0]
        ctlx[-1] = x0[-1]
        ctly[0]  = y0[0]
        ctly[-1] = y0[-1]

        dv_counter += 2*Nctl-4
    else:
        ctlx = ctlx0
        ctly = ctly0
        
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

    interp_x = pyspline.bvaluv(t,ctlx,k,0,s)
    interp_y = pyspline.bvaluv(t,ctly,k,0,s)
    err = sum((interp_x-x0)**2 + (interp_y-y0)**2)
    f = sqrt(err/N)

    fcon = array([])
            
    if USE_KNOTS:
        fcon = hstack([fcon,zeros(Nctl-k)])
        
        for i in xrange(Nctl-k):
            fcon[i] = t[i+k+1]-t[i+k]
        # end for
    # end if 

    # Derivative constraint
    dx = pyspline.bvalu(t,ctlx,k,1,0)
    fcon = hstack([fcon,dx])
    
    fail = False

    return f,fcon,fail

def sens(x,f_obj,f_con):
    global USE_CTLPTS
    global USE_KNOTS
    global USE_PARAMS
    use_params_start = 0

    ndv = len(x)

    ncon = 0
    if USE_KNOTS:
        ncon += Nctl-k
    ncon += 1 # derivative constraint

    g_obj = zeros(ndv) #objective derivative vector
    g_con = zeros([ncon, ndv]) #constraint jacobian
    
    h = 1.0e-40j #complex step size
    
    for j in xrange(ndv):

        x_deriv = zeros([ndv],'complex')
        
        for i in range(ndv): #copy over design variables
            x_deriv[i] = complex(x[i],0)

        x_deriv[j] += h
     
        dv_counter = 0
    
        if USE_CTLPTS:
            ctlx = zeros(Nctl,'complex')
            ctly = zeros(Nctl,'complex')

            ctlx[1:Nctl-1] = x_deriv[0:Nctl-2] # These are the control points
            ctly[1:Nctl-1] = x_deriv[Nctl-2:2*Nctl-4]

            ctlx[0]  = x0[0]
            ctlx[-1] = x0[-1]
            ctly[0]  = y0[0]
            ctly[-1] = y0[-1]

            dv_counter += 2*Nctl-4
        else:
            ctlx = ctlx0
            ctly = ctly0

        if USE_KNOTS:
            t = zeros(Nctl+k,'complex')
            t[-4:] = 1
            t[k:Nctl] = x_deriv[dv_counter:dv_counter + Nctl-k]

            dv_counter += Nctl-k
        else:
            t = t0

        if USE_PARAMS:
            use_params_start = dv_counter;
            s = x_deriv[dv_counter:dv_counter + N]
            dv_counter += N
        else:
            s = s0

        err = 0.0;

        if not USE_PARAMS or j<use_params_start: #Do derivatives the 'dumb' way
            interp_x = pyspline_cs.bvaluv(t,ctlx,k,0,s)
            interp_y = pyspline_cs.bvaluv(t,ctly,k,0,s)
            err = sum((interp_x-x0)**2 + (interp_y-y0)**2)
            f = sqrt(err/N)
            g_obj[j] = imag(f)/imag(h)
        else: #Do the s derivatives in a more elegent fashion
            index = j-use_params_start
            interp_x = pyspline_cs.bvalu(t,ctlx,k,0,s[index])
            interp_y = pyspline_cs.bvalu(t,ctly,k,0,s[index])
            temp = (interp_x-x0[index])**2  + (interp_y-y0[index])**2 
            err = f_obj**2*N + temp - real(temp)
            f = sqrt(err/N)
            g_obj[j] = imag(f)/imag(h)
        # end if

        
        fcon = array([],'complex')
            
        if USE_KNOTS:
            fcon = hstack([fcon,zeros(Nctl-k,'complex')])
            
            for i in xrange(Nctl-k):
                fcon[i] = t[i+k+1]-t[i+k]
            # end for
        # end if 
            
        # Derivative constraint
        dx = pyspline_cs.bvalu(t,ctlx,k,1,0)
        fcon = hstack([fcon,dx])
            
        g_con[:,j] = imag(fcon)/imag(h)
        
    # end for (dv loop)
    fail = False
    return g_obj,g_con,fail

# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

#Make sure the end points are what we want:

if USE_CTLPTS:
    opt_prob.addVarGroup('CTLx',Nctl-2,'c',value=ctlx0,lower=-1,upper=1)
    opt_prob.addVarGroup('CTLy',Nctl-2,'c',value=ctly0,lower=-1,upper=1)
if USE_KNOTS:
    opt_prob.addVarGroup('Knot',Nctl-k,'c',value=t0[k:Nctl],lower=0,upper=1)
if USE_PARAMS:
    opt_prob.addVarGroup('s:',N,'c',value=s0,lower=0,upper=1)

# ===================
#  Constraints
# ===================
if USE_KNOTS:
    opt_prob.addConGroup('Knots', Nctl-k,type='i', lower=0.001, upper=.33)
opt_prob.addCon('LEconstraint',type = 'i',lower=0.0,upper=0.0)

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
opt.setOption('Derivative level',3)
opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-6)
opt.setOption('Major feasibility tolerance',1e-6)
opt.setOption('Minor feasibility tolerance',1e-6)

# ===================
#  Run Optimization
# ===================

opt(opt_prob,sens)

# ===================
#  Print Solution  
# ===================

print opt_prob._solutions[0]

plot(ctlx,ctly,'ro')
plot(x0,y0)
show()
