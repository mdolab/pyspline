 # This simple script test some of the volume functionality in pySpline
import numpy
from pyspline import pySpline

X = numpy.zeros((2,2,2,3))
X[0,0,0,:] = [0,0,0]
X[1,0,0,:] = [1.1,-.1,-.1]
X[1,1,0,:] = [0.9,1.05,.2]
X[0,1,0,:] = [-.1,1.1,0]

X[0,0,1,:] = [-.1,.1,1.5]
X[1,0,1,:] = [1.2,-.2,1.8]
X[1,1,1,:] = [1.2,1.0,2]
X[0,1,1,:] = [-.2,1.3,2.1]

vol = pySpline.Volume(X=X, ku=2, kv=2, kw=2, Nctlu=2, Nctlv=2, Nctlw=2)
vol.writeTecplot('vol.dat',orig=True)

# Generate random data 
M = 10000
Y = numpy.zeros((M,3))
for i in range(M):
    Y[i,0] = numpy.random.random()*.5+.25
    Y[i,1] = numpy.random.random()*.5+.25
    Y[i,2] = numpy.random.random()*.5+.25

import time
timeA = time.time()
u,v,w,D = vol.projectPoint(Y)
print('Time to project %d points: %f seconds:'%(M, time.time()-timeA))
