# This is a test script to test the functionality of the
# pySpline surface
import numpy
from pyspline import pySpline

# Create a generic surface
nu = 20
nv = 20
u = numpy.linspace(0, 4, nu)
v = numpy.linspace(0, 4, nv)
[V, U] = numpy.meshgrid(v, u)
Z = numpy.cos(U) * numpy.sin(V)
surf = pySpline.Surface(x=U, y=V, z=Z, ku=4, kv=4, Nctlu=5, Nctlv=5)
surf.writeTecplot("surface.dat")

n = 100
theta = numpy.linspace(0.0000, 2 * numpy.pi, n)
x = numpy.cos(theta) - 1
y = numpy.sin(theta) + 1
z = numpy.linspace(0, 1, n) + 2
curve = pySpline.Curve(x=x, y=y, z=z, k=4, Nctl=16, niter=100)
curve.writeTecplot("helix.dat")

u, v, s, D = surf.projectCurve(curve, Niter=100, eps1=1e-10, eps2=1e-10, u=1, v=1, s=1)

print(u, v, s, D)
print(curve(s))
print(surf(u, v))
