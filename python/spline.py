# This is a quick and dirty implementation of a linear spline function
from numpy import *
import pdb,sys,time
import pyspline
import pySpline
from scipy.linalg import lstsq,solve
from mdo_import_helper import *
from scipy.io import savemat
from matplotlib.pylab import plot,show,axis
import copy
from lsqr import *

nacax,nacay = loadtxt('naca0012.dat',unpack=True)

# Get some new data
#n=4
#k=4
#theta = linspace(0,pi/2,10)
#x = cos(theta)
#y = sin(theta)
#upper = vstack([x,y]).T.astype('D')
#top = pySpline.linear_spline(task='lms',k=k,Nctl=n,X=upper)

#knots = top.t
#coef = zeros((n,3))
#coef[:,0:2] = top.coef
#coef[:,2] = 1.0
#print 'knots:',knots
#print 'coef:',coef

# #Test the compute curve
# k = 3
# Nctl = 3
# n = 25
# theta = linspace(0,pi/2,n)
# x = cos(theta)
# y = sin(theta)
# X = vstack([x,y]).T.astype('d')
# s0 = linspace(0,1,n)
# knots = pyspline.knots(s0,Nctl,k)
# s = copy.deepcopy(s0)
# s,coef = pyspline.compute_curve(s,X,knots,k,Nctl)
# print coef

#sys.exit(0)
k = 4
Nctl = 40
n = len(nacax)
X = vstack([nacax,nacay]).T.astype('d')
s0 = linspace(0,1,n)
knots = pyspline.knots(s0,Nctl,k)
s = copy.deepcopy(s0)
#s,coef = pyspline.compute_curve(s,X,knots,k,Nctl)
coef = zeros((Nctl,2))
coef[:,1] = 1.0

Jac = pyspline.curve_jacobian_linear(Nctl,knots,k,s)
global direct_time
global trans_time 
global direct_count
global trans_count
direct_time = 0.0
trans_time = 0.0
direct_count = 0
trans_count = 0
def aprod(mode,m,n,u):
    '''Do the dot product'''
    global direct_time
    global trans_time 
    global direct_count
    global trans_count

    if mode == 1:
        timeA = time.time()
        val = dot(Jac,u)
        direct_time += time.time()-timeA
        direct_count += 1
        return val

#         coef[:,0] = u
#         val = pyspline.eval_curve_v(s,knots,k,coef)
#         return val[0]
    else:
        trans_count += 1
        timeA = time.time()
        val = dot(Jac.T,u)
        trans_time += time.time()-timeA
        return val
    # end if

nrow = n
ncol = Nctl
timeA = time.time()
for i in xrange(100):
    res = lsqr( nrow,ncol, aprod, X, 0, 1e-4,1e-4, 1e8, 100, 0)
print 'lsqr time:',time.time()-timeA
print 'direct_time:',direct_time
print 'trans_time:',trans_time
print 'counts:',direct_count,trans_count
print 'x:',res[0]
print 'res:',res[1:]

timeB = time.time()
for i in xrange(100):
    res2 = lstsq(Jac,X[:,0])
print 'lstsq time:',time.time()-timeB
print 'res2:',res2[0]

diff = res[0]-res2[0]
print 'diff:',sqrt(dot(diff,diff))
sys.exit(0)
# Check the rms
tot =0.0
for i in xrange(n):
    res = pyspline.eval_curve(s[i],knots,k,coef)
    tot = tot + (res[0] - X[i,0])**2 + (res[1] - X[i,1])**2
# end for
print 'rms:',sqrt(tot/n)




sys.exit(0)
# coef = zeros((3,3))
# coef[0,:] = [1,0,1]
# coef[1,:] = [1,1,sqrt(2)/2]
# coef[2,:] = [0,1,1]
# knots = [0,0,0,pi/2,pi/2,pi/2]
# k=3

# # Check Tangent
# val = pyspline.eval_curve_deriv(0.0,knots,k,coef)
# print 'val',val



# Plotting Stuff
n=25
#s = linspace(0,pi/2,n)
#s = linspace(0,1,n)
val = zeros((len(s),2))
for i in xrange(len(s)):
    val[i] = pyspline.eval_curve(s[i],knots,k,coef)
# end for
plot(val[:,0],val[:,1],'ko-')
plot(coef[:,0],coef[:,1],'go')

# Plot an ACTUAL circle
theta = linspace(0,pi/2,n)
x = cos(theta)
y = sin(theta)
plot(x,y,'ks--')
#plot(nacax,nacay,'ks--')
show()



