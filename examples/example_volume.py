# This simple script test some of the volume functionality in pySpline
# Standard Python modules
import time

# External modules
import numpy as np

# First party modules
from pyspline import Volume

X = np.zeros((2, 2, 2, 3))
X[0, 0, 0, :] = [0, 0, 0]
X[1, 0, 0, :] = [1.1, -0.1, -0.1]
X[1, 1, 0, :] = [0.9, 1.05, 0.2]
X[0, 1, 0, :] = [-0.1, 1.1, 0]

X[0, 0, 1, :] = [-0.1, 0.1, 1.5]
X[1, 0, 1, :] = [1.2, -0.2, 1.8]
X[1, 1, 1, :] = [1.2, 1.0, 2]
X[0, 1, 1, :] = [-0.2, 1.3, 2.1]

vol = Volume(X=X, ku=2, kv=2, kw=2, Nctlu=2, Nctlv=2, Nctlw=2)
vol.writeTecplot("vol.dat", orig=True)

# Generate random data
M = 10000
Y = np.zeros((M, 3))
for i in range(M):
    Y[i, 0] = np.random.random() * 0.5 + 0.25
    Y[i, 1] = np.random.random() * 0.5 + 0.25
    Y[i, 2] = np.random.random() * 0.5 + 0.25


timeA = time.time()
u, v, w, D = vol.projectPoint(Y)
print("Time to project %d points: %f seconds:" % (M, time.time() - timeA))
