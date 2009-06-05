# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace,cos,pi,hstack,zeros,ones,sqrt
from matplotlib.pylab import plot,show

# =============================================================================
# Extension modules
# =============================================================================

# pySpline
sys.path.append('../../python')
from pyspline import *

#pyOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/')
from pyOpt_optimization import Optimization

#pySNOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/pySNOPT')
from pySNOPT import SNOPT

'''This script runs an example of fitting a cubic spline to a set of
airfoil coordinates, specificially the upper surface of a naca 0012
airfoil with a constraint on a vertical leading edge'''

N = 131
Nctl = 25
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
    t = zeros(Nctl+k)
    t[-4:] = 1
    t[k:Nctl] = x[2*Nctl-4:2*Nctl-4+Nctl-k]
    #t=t0
    s = x[3*Nctl-k-4:]
    #s = s0
    
    err = 0.0;

    for i in xrange(len(s)):
        if s[i] > 1:
            s[i] = 1
        if s[i] < 0:
            s[i] = 0
        
        
        interp_x,inbv = bvalu(t,ctlx,k,0,s[i])
        interp_y,inbv = bvalu(t,ctly,k,0,s[i])
        err = err + (interp_x-x0[i])**2  + (interp_y-y0[i])**2 

    f = sqrt(err/N)
    fcon = zeros(Nctl-k)
    for i in xrange(Nctl-k):
        fcon[i] = t[i+k+1]-t[i+k]

    fail = False
    counter += 1
    return f,fcon,fail




# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

ctlx_i = linspace(0,1,Nctl)
ctly_i = 0.05*ones(Nctl)

opt_prob.addVarGroup('CTLx',Nctl-2,'c',value=ctlx_i[1:Nctl-1], lower=-1, upper=1)
opt_prob.addVarGroup('CTLy',Nctl-2,'c',value=ctly_i[1:Nctl-1], lower=-1, upper=1)
opt_prob.addVarGroup('Knot',Nctl-k,'c',value=t0[k:Nctl],lower=0,upper=1)
opt_prob.addVarGroup('s:',N,'c',value=s0,lower=0,upper=1)
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
#opt.setOption('Verify level',3)
opt.setOption('Derivative level',0)
opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-6)
opt.setOption('Major feasibility tolerance',1e-6)
opt.setOption('Minor feasibility tolerance',1e-6)
#opt.setOption('Major step limit',8e-3)
#opt.setOption('Function precision',1.e-5)

# ===================
#  Run Optimization
# ===================

opt(opt_prob)

# ===================
#  Print Solution  
# ===================

print opt_prob._solutions[0]

print 'counter:',counter
plot(x0,y0,'k-')
plot(ctlx,ctly,'rs')

show()

