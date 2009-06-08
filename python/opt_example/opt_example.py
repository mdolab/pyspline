# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,hstack,zeros,ones,sqrt,imag,interp
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

N = 81
Nctl = 13
k = 4
x0 = 0.5*(1-cos(linspace(0,pi,N)))
#Top surface of a naca0012 airfoil
y0= 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)

s0 = zeros(N)
for i in xrange(N-1):
    s0[i+1] = s0[i] + sqrt((x0[i+1]-x0[i])**2 + (y0[i+1]-y0[i])**2)

s0 = s0/s0[-1] #Normalize s

# do a cosine initial knot vector 
t0 = zeros(Nctl + k)
t0[0:k] = s0[0]
t0[Nctl:Nctl+k] = s0[-1]
t0[k-1:Nctl+1]= 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))
#t0[k-1:Nctl+1] = linspace(0,1,Nctl-k+2)
print 'Initialize T:',t0

global counter
global ctlx
global ctly

counter = 0
def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global counter
    global ctlx
    global ctly
    ctlx = zeros(Nctl)
    ctly = zeros(Nctl)

    ctlx[1:Nctl-1] = x[0:Nctl-2] # These are the control points
    ctly[1:Nctl-1] = x[Nctl-2:2*Nctl-4]
    
    #Fix the second ctrl pt to 0
    ctlx[1] = 0

    ctlx[-1] = 1
    ctly[-1] = y0[-1]
#    t = zeros(Nctl+k)
#    t[-4:] = 1
#    t[k:Nctl] = x[2*Nctl-4:2*Nctl-4+Nctl-k]
    t=t0
    #s = x[3*Nctl-k-4:]
    s = s0
    
    err = 0.0;

    for i in xrange(len(s)):
        if s[i] > 1:
            s[i] = 1
        if s[i] < 0:
            s[i] = 0
        
        
        interp_x,inbv = pyspline.bvalu(t,ctlx,k,0,s[i])
        interp_y,inbv = pyspline.bvalu(t,ctly,k,0,s[i])
        err = err + (interp_x-x0[i])**2  + (interp_y-y0[i])**2 

    f = sqrt(err/N)
    fcon = zeros(Nctl-k)
    for i in xrange(Nctl-k):
        fcon[i] = t[i+k+1]-t[i+k]

    fail = False

    print counter,f
    counter += 1

    return f,fcon,fail

def sens(x,f_obj,f_con):
    k=4
    ndv = len(x)
    ncon = Nctl-k
    g_obj = zeros(ndv) #objective derivative vector
    g_con = zeros([ncon, ndv]) #constraint jacobian

    t=t0
    s = s0

    x_passed = copy.deepcopy(x)
    h = 1.0e-40j #complex step size

    for j in xrange(ndv):
        x_deriv = zeros([ndv],'complex')

        for i in range(ndv): #copy over design variables
            x_deriv[i] = complex(x[i],0)

        x_deriv[j] += h

        ctlx = zeros(Nctl,'complex')
        ctly = zeros(Nctl,'complex')

        ctlx[1:Nctl-1] = x_deriv[0:Nctl-2] # These are the control points
        ctly[1:Nctl-1] = x_deriv[Nctl-2:2*Nctl-4]
    
        # Fix the second ctrl pt to 0
        ctlx[1] = 0

        ctlx[-1] = 1
        ctly[-1] = y0[-1]
#        t = zeros(Nctl+k,'complex')
#        t[-4:] = 1
#        t[k:Nctl] = x_deriv[2*Nctl-4:2*Nctl-4+Nctl-k]
#        s = x_deriv[3*Nctl-k-4:]
        err = 0.0;
        
        for i in xrange(len(s)):
            if s[i] > 1:
                s[i] = 1
            if s[i] < 0:
                s[i] = 0
                
            interp_x,inbv = pyspline_cs.bvalu(t,ctlx,k,0,s[i])
            interp_y,inbv = pyspline_cs.bvalu(t,ctly,k,0,s[i])
            err = err + (interp_x-x0[i])**2  + (interp_y-y0[i])**2 
            

        f = sqrt(err/N)
        g_obj[j] = imag(f)/imag(h)

        fcon = zeros(Nctl-k)
        for i in xrange(Nctl-k):
            fcon[i] = t[i+k+1]-t[i+k]
        # end for
            
        g_con[:,j] = imag(fcon)/imag(h)

    #end for (dv loop
    fail = False
    return g_obj,g_con,fail

        


# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

ctlx_i = linspace(0,1,Nctl)
ctly_i = interp(ctlx_i,x0,y0)

opt_prob.addVarGroup('CTLx',Nctl-2,'c',value=ctlx_i[1:Nctl-1],lower=-1,upper=1)
opt_prob.addVarGroup('CTLy',Nctl-2,'c',value=ctly_i[1:Nctl-1],lower=-1,upper=1)
#opt_prob.addVarGroup('Knot',Nctl-k,'c',value=t0[k:Nctl],lower=0,upper=1)
#opt_prob.addVarGroup('s:',N,'c',value=s0,lower=0,upper=1)
# ===================
#  Constraints
# ===================

opt_prob.addConGroup('Knots', Nctl-k,type='i', lower=0.0, upper=10.0)

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

print 'counter:',counter
plot(x0,y0,'k-')
plot(ctlx,ctly,'rs')

show()

