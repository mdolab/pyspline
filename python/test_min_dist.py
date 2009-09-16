from numpy import *
import pySpline
from matplotlib.pylab import *
# Create an linear spline:

x = [1362,1217]
y = [452,452]
z = [128,178]

X = zeros((2,3))
X[:,0] = x
X[:,1] = y
X[:,2] = z

curve1 = pySpline.linear_spline(task='interpolate',X=X,k=2)
x = [1147,1314,1804]
y = [119,427,1156]
z = [176,181,264]

X = zeros((3,3))
X[:,0] = x
X[:,1] = y
X[:,2] = z

curve2 = pySpline.linear_spline(task='interpolate',X=X,k=4)

s = 0.0
t=0.0
s,t,D,converged = curve1.minDistance(curve2,s=s,t=t,Niter=400,tol=1e-10)

print 's,t:',s,t
print 'D:',D
print 'converged:',converged

# s0 = linspace(0,1,50)

# x1 = curve1.getValueV(s0)
# x2 = curve2.getValueV(s0)

# v1 = curve1.getValue(s)
# v2 = curve2.getValue(t)

# s0 = linspace(0,2*pi)
# x = cos(s0)
# y = sin(s0)

# plot(x,y)
# axis('equal')

# plot(x1[:,0],x1[:,1],'ko-')

# plot(x2[:,0],x2[:,1],'gs-')
# plot([v1[0],v2[0]],[v1[1],v2[1]],'rs')

# show()
