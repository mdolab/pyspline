# This is a quick and dirty implementation of a linear spline function
from numpy import *
import pdb,sys,time
import pyspline
import pySpline
from scipy.linalg import lstsq,solve
from mdo_import_helper import *
from scipy.io import savemat
from matplotlib.pylab import plot,show,axis

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

# Test the compute curve
# k = 4
# Nctl = 5
# n = 6
# theta = linspace(0,pi/2,n)
# x = cos(theta)
# y = sin(theta)
# X = vstack([x,y]).T.astype('d')
# s = linspace(0,1,n)
# knots = pyspline.knots(s,Nctl,k)
# s,coef = pyspline.compute_curve(s,X,knots,k,Nctl)
# coef[3,2] = 5
# print coef


coef = zeros((3,3))
coef[0,:] = [1,0,1]
coef[1,:] = [1,1,sqrt(2)/2]
coef[2,:] = [0,1,1]
knots = [0,0,0,pi/2,pi/2,pi/2]
k=3

# Check Tangent
val = pyspline.eval_curve_deriv(0.0,knots,k,coef)
print 'val',val



# Plotting Stuff
n=25
s = linspace(0,pi/2,n)
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
show()



