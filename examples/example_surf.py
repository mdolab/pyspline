# This is a test script to test the functionality of the
# pySpline surface
# External modules
import numpy as np

# First party modules
from pyspline import Curve, Surface

# Create a generic surface
nu = 20
nv = 20
u = np.linspace(0, 4, nu)
v = np.linspace(0, 4, nv)
[V, U] = np.meshgrid(v, u)
Z = np.cos(U) * np.sin(V)
surf = Surface(x=U, y=V, z=Z, ku=4, kv=4, Nctlu=5, Nctlv=5)
surf.writeTecplot("surface.dat")

n = 100
theta = np.linspace(0.0000, 2 * np.pi, n)
x = np.cos(theta) - 1
y = np.sin(theta) + 1
z = np.linspace(0, 1, n) + 2
curve = Curve(x=x, y=y, z=z, k=4, Nctl=16, niter=100)
curve.writeTecplot("helix.dat")

u, v, s, D = surf.projectCurve(curve, Niter=100, eps1=1e-10, eps2=1e-10, u=1, v=1, s=1)

print(u, v, s, D)
print(curve(s))
print(surf(u, v))
