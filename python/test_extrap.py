from numpy import *
import pySpline
from matplotlib.pylab import *
# Create an linear spline:

x = [0,1,0,-1,0]
y = [-1,0,1,0,-1]
z = [0,0,0,0,0]

X = zeros((5,3))
X[:,0] = x
X[:,1] = y
X[:,2] = z

curve = pySpline.linear_spline(task='interpolate',X=X,k=4)
N = 120
interp_s = linspace(-.1,1.1,N)
vals = zeros((N,3))
for i in xrange(len(interp_s)):
    vals[i] = curve.getValue(interp_s[i])

plot(vals[:,0],vals[:,1],'ko-')
plot(x,y,'gs')
axis('equal')
show()

