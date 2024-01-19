# This is a test script to test the functionality of the
# pySpline surface
# External modules
import numpy as np

# First party modules
from pyspline import Curve, Surface

# Create a generic surface
nu = 20
nv = 20
u = np.linspace(-4, 4, nu)
v = np.linspace(-4, 4, nv)
[V, U] = np.meshgrid(v, u)
Z = np.cos(U) * np.sin(V)
surf = Surface(x=U, y=V, z=Z, ku=4, kv=4, nCtlu=5, nCtlv=5)
surf.writeTecplot("surface.dat")

n = 100
theta = np.linspace(0.0000, 2 * np.pi, n)
x = np.cos(theta) - 1
y = np.sin(theta) + 1
z = np.linspace(0, 1, n) + 2
curve = Curve(x=x, y=y, z=z, k=4, nCtl=16, nIter=100)
curve.writeTecplot("helix.dat")

u, v, s, D = surf.projectCurve(curve, nIter=100, eps=1e-10, u=0.5, v=0.5, s=0.5)
val1 = curve(s)
val2 = surf(u, v)

print(f"{u=}, {v=}, {s=}, {D=}")
print(val1)
print(val2)

# Output the data
f = open("surface_projections.dat", "w")
f.write('VARIABLES = "CoordinateX", "CoordinateY", "CoordinateZ"\n')
f.write("Zone T=surf_curve_proj I=2 \n")
f.write("DATAPACKING=POINT\n")
f.write("%f %f %f\n" % (val1[0], val1[1], val1[2]))
f.write("%f %f %f\n" % (val2[0], val2[1], val2[2]))