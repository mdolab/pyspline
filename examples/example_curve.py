# This is a test script to test the functionality of the
# pySpline curve class

# External modules
import numpy as np

# First party modules
from pyspline import Curve

# Get some Helix-like data

n = 100
theta = np.linspace(0.0000, 2 * np.pi, n)
x = np.cos(theta)
y = np.sin(theta)
z = np.linspace(0, 1, n)
print("Helix Data")
curve = Curve(x=x, y=y, z=z, k=4, nCtl=16, nIter=100)
curve.writeTecplot("helix.dat")

# Load naca0012 data
print("Naca 0012 data")
x, y = np.loadtxt("naca0012", unpack=True)
curve = Curve(x=x, y=y, k=4, nCtl=11, nIter=500)
curve.writeTecplot("naca_data.dat")

# Projection Tests
print("Projection Tests")
x = [0, 2, 3, 5]
y = [-2, 5, 3, 0]
z = [0, 0, 0, 0]
curve1 = Curve(x=x, y=y, z=z, k=4)
curve1.writeTecplot("curve1.dat")

x = [-2, 5, 2, 1]
y = [5, 1, 4, 2]
z = [3, 0, 1, 4]
curve2 = Curve(x=x, y=y, z=z, k=4)
curve2.writeTecplot("curve2.dat")

# Get the minimum distance distance between a point and each curve
x0 = [4, 4, 3]
s1, D1 = curve1.projectPoint(x0, s=0.5)
val1 = curve1(s1)  # Closest point on curve
s2, D2 = curve2.projectPoint(x0, s=1.0)
val2 = curve2(s2)  # Closest point on curve

# Output the data
f = open("projections.dat", "w")
f.write('VARIABLES = "CoordinateX", "CoordinateY", "CoordinateZ"\n')
f.write("Zone T=curve1_proj I=2 \n")
f.write("DATAPACKING=POINT\n")
f.write("%f %f %f\n" % (x0[0], x0[1], x0[2]))
f.write("%f %f %f\n" % (val1[0], val1[1], val1[2]))

f.write("Zone T=curve2_proj I=2 \n")
f.write("DATAPACKING=POINT\n")
f.write("%f %f %f\n" % (x0[0], x0[1], x0[2]))
f.write("%f %f %f\n" % (val2[0], val2[1], val2[2]))

# Get the minimum distance between the two curves
s, t, D = curve1.projectCurve(curve2)
val1 = curve1(s)
val2 = curve2(t)

f.write("Zone T=curve1_curve2 I=2 \n")
f.write("DATAPACKING=POINT\n")
f.write("%f %f %f\n" % (val1[0], val1[1], val1[2]))
f.write("%f %f %f\n" % (val2[0], val2[1], val2[2]))
