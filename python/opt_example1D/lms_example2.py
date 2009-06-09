# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, transpose, dot
import scipy.linalg
from matplotlib.pylab import plot,show

# =============================================================================
# Extension modules
# =============================================================================

# pySpline
sys.path.append('../../python')
import pyspline
import pyspline_cs

'''This script runs an example of lms fitting a cubic spline to a set of
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
N = 51
x0 = 0.5*(1-cos(linspace(0,pi,N)))
y0= 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)

x0 = linspace(0,1,N)
y0 = x0

# Top surface of a WT airfoil
# N = 34
# x0 = array([0.00000,0.00018,0.00255,0.00954,0.02090,0.03650,0.05640,0.08030,0.10801,0.13934,0.17395,0.21146,0.25149,0.29361,0.33736,0.38228,0.42820,0.47526,0.52324,0.57161,0.61980,0.66724,0.71333,0.75749,0.79915,0.83778,0.87287,0.90391,0.93072,0.95355,0.97251,0.98719,0.99668,1.00000])
# y0 = array([0.00000,0.00159,0.00748,0.01640,0.02600,0.03580,0.04560,0.05520,0.06430,0.07290,0.08070,0.08760,0.09340,0.09810,0.10133,0.10294,0.10249,0.10005,0.09610,0.09090,0.08490,0.07820,0.07100,0.06340,0.05570,0.04800,0.04030,0.03260,0.02480,0.01700,0.00982,0.00431,0.00103,0.00000])
#x0 = x0[::-1]
#y0 = y0[::-1]
# S-Parameter Calculation
s = zeros(N)
for i in xrange(N-1):
    s[i+1] = s[i] + sqrt((x0[i+1]-x0[i])**2 + (y0[i+1]-y0[i])**2)

s = s/s[-1] #Normalize s

# Initial Knot Vector Calculation
t= zeros(Nctl + k)
t[0:k] = s[0]
t[Nctl:Nctl+k] = s[-1]
t[k-1:Nctl+1] = 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))

# Initial Control Pt Vector
ctlx0 = 0.5*(1-cos(linspace(0,pi,Nctl)))# good Estimate
ctly0 = interp(ctlx0[0:N],x0[0:N],y0)     # good Estimate

ctlx0 = zeros(Nctl) # Zeros
ctly0 = zeros(Nctl) # Zeros
ctlx = array(ctlx0,'D') # Force to Complex
ctly = array(ctly0,'D') # Force to Complex

ctlx[0]  = x0[0]
ctlx[-1] = x0[-1]
ctly[0]  = y0[0] 
ctly[-1] = y0[-1]

# Determine the Jacobian matrix pf(s)/pC (partial function val wrt coefficient)
h = 1.0e-40j

J = zeros([N,Nctl-2])

timeA = time.time()
Wx = zeros([N,N])
Wy = zeros([N,N])
for i in xrange(N):
    Wx[i,i] = 1
#Wx[0,0] = 10
# Wx[-1,-1] = 10
# Wx[-2,-2] = 10           
# Wx[-3,-3] = 10           
# Wx[-4,-4] = 10           
# Wx[-5,-5] = 10           
# Wx[-6,-6] = 10           
# Wx[-7,-7] = 10           

for i in xrange(N):
    Wy[i,i] = 1
# Wy[0,0] = 10000000
# Wy[-1,-1] = 1000000

for i in xrange(Nctl-2):
    
    ctlx[i+1] += h
    
    interp_x = pyspline_cs.bvaluv(t,ctlx,k,0,s)
    ctlx[i+1] -= h    

    J[:,i] = imag(interp_x)/imag(h)

# end for 

Jt = transpose(J)
# 
#ctlx = dot(dot(scipy.linalg.inv(dot(Jtx,dot(Wx,Jx))),Jtx),dot(Wx,x0))
#ctly = dot(dot(scipy.linalg.inv(dot(Jty,dot(Wy,Jy))),Jty),dot(Wy,y0))

ctlx = dot(dot(scipy.linalg.inv(dot(Jt,J)),Jt),x0)
ctly = dot(dot(scipy.linalg.inv(dot(Jt,J)),Jt),y0)

ctlx = hstack([x0[0],ctlx,x0[-1]])
ctly = hstack([y0[0],ctly,y0[-1]])
ctlx[0]  = x0[0]
ctlx[-1] = x0[-1]
ctly[0]  = y0[0] 
ctly[-1] = y0[-1]
print 'Cx:',ctlx
print 'Cy:',ctly

#Calc the RMS Error
interp_x = pyspline.bvaluv(t,ctlx,k,0,s)
interp_y = pyspline.bvaluv(t,ctly,k,0,s)
err = sum((interp_x-x0)**2) + sum((interp_y-y0)**2)
f = sqrt(err/N)

print 'RMS:',f
print 'time:',time.time()-timeA

print ctlx.shape
print ctly.shape

#plot(interp_x,interp_y)
#plot(ctlx,ctly,'ro-')
#plot(ctlx0,ctly0,'go--')
#plot(x0,y0)

plot(s,interp_x)
plot(s,x0,'o')
plot(s,y0,'x')
show()
