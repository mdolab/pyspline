# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, transpose, dot, sin, cos
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

Nctl = 21  # Number of control points
Ndc  = 0  # Number of derivative constraints
k = 4

# Surface of a naca0012 airfoil
# N = 581
# x0 = 0.5*(1-cos(linspace(0,pi,N)))
# y0 = 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)

# Surface of WT Airfoil
N=32
x0u = array([0.0000000,0.0001800,0.0025500,0.0095400,0.0209000,0.0365000,0.0564000,0.0803000,0.1080100,0.1393400,0.1739500,0.2114600,0.2514900,0.2936100,0.3373600,0.3822800,0.4282000,0.4752600,0.5232400,0.6198000,0.6672400,0.7133300,0.7574900,0.7991500,0.8377800,0.9039100,0.9307200,0.9535500,0.9725100,0.9871900,0.9966800,1.0000000])
y0u = array([0.0000000,0.0015900,0.0074800,0.0164000,0.0260000,0.0358000,0.0456000,0.0552000,0.0643000,0.0729000,0.0807000,0.0876000,0.0934000,0.0981000,0.1013300,0.1029400,0.1024900,0.1000500,0.0961000,0.0849000,0.0782000,0.0710000,0.0634000,0.0557000,0.0480000,0.0326000,0.0248000,0.0170000,0.0098200,0.0043100,0.0010300,0.0000000])
x0l = array([1.0000000,0.9965500,0.9860400,0.9681300,0.9425700,0.9094500,0.8695300,0.8240900,0.7743400,0.7214200,0.6664400,0.6105500,0.5548400,0.5003200,0.4478500,0.3977900,0.3502700,0.3049700,0.2615300,0.2198700,0.1796500,0.1413300,0.1063700,0.0758000,0.0500000,0.0292000,0.0137000,0.0036700,0.0021600,0.0009300,0.0002100,0])
y0l = array([0.0000000,0.0008400,0.0032400,0.0068900,0.0111000,0.0150000,0.0174000,0.0178000,0.0162000,0.0125000,0.0069900,-0.0001500,-0.0084100,-0.0172000,-0.0256000,-0.0328000,-0.0379000,-0.0405000,-0.0406000,-0.0386000,-0.0349000,-0.0309000,-0.0273000,-0.0236000,-0.0196000,-0.0152000,-0.0104000,-0.0052500,-0.0040300,-0.0027400,-0.0014600,0])
x0l = x0l[::-1]
y0l = y0l[::-1]

# theta = 0.0
# x0u = cos(theta)*x0 - sin(theta)*y0
# y0u = sin(theta)*x0 + cos(theta)*y0

# x0l = cos(theta)*x0 - sin(theta)*(-y0)
# y0l = sin(theta)*x0 + cos(theta)*(-y0)

# S-Parameter Calculation (same for upper and lower surfaces)
s = zeros(N)
for i in xrange(N-1):
    s[i+1] = s[i] + sqrt((x0l[i+1]-x0l[i])**2 + (y0l[i+1]-y0l[i])**2)

s = s/s[-1] #Normalize s

# Initial Knot Vector Calculation
t= zeros(Nctl + k)
t[0:k] = s[0]
t[Nctl:Nctl+k] = s[-1]
t[k-1:Nctl+1] = 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))

# Determine the Jacobian matrix pf(s)/pC (partial function val wrt coefficient)
h = 1.0e-40j
J = zeros([N,Nctl+Ndc])
ctl = zeros([Nctl],'D')

timeA = time.time()

for i in xrange(Nctl):
    ctl[i] += h
    val = pyspline_cs.bvaluv(t,ctl,k,0,s)
    ctl[i] -= h    
    J[:,i] = imag(val)/imag(h)
# end for 

for i in xrange(Ndc):
    print 'doine Ndc loop'
    ctl[i] += h
    val = pyspline_cs.bvaluv(t,ctl,k,1,s)
    ctl[i] -= h
    J[:,i+Nctl] = imag(val)/imag(h)
# end for

# Now call the fortran function to fit the splines together

coefs = pyspline.fit_foil(J,x0u,x0l,y0u,y0l)

print 'coefs:',coefs


#Calc the RMS Error
interp_x = pyspline.bvaluv(t,coefs[:,0],k,0,s)
interp_y = pyspline.bvaluv(t,coefs[:,1],k,0,s)
err = sum((interp_x-x0u)**2) + sum((interp_y-y0u)**2)
f = sqrt(err/N)
print 'rms:',f

print 'timeB:',time.time()-timeA


plot(x0u,y0u,'ko-')
plot(x0l,y0l,'ko-')
# plot(coefs[:,0],coefs[:,1],'ro')
# plot(coefs[:,2],coefs[:,3],'go')

s = 0.5*(1-cos(linspace(0,pi,400)))
interp_xu = pyspline.bvaluv(t,coefs[:,0],k,0,s)
interp_yu = pyspline.bvaluv(t,coefs[:,1],k,0,s)
interp_xl = pyspline.bvaluv(t,coefs[:,2],k,0,s)
interp_yl = pyspline.bvaluv(t,coefs[:,3],k,0,s)
plot(interp_xu,interp_yu)
plot(interp_xl,interp_yl)
show()

