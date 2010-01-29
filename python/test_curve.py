from numpy import *
import pySpline

# Get some Helix-like data

n = 100
theta = linspace(0,2*pi,n)
x = cos(theta)
y = sin(theta)
z = linspace(0,1,n)
print 'Helix Data'
curve = pySpline.curve('lms',x=x,y=y,z=z,k=4,Nctl=10)
curve.runParameterCorrection(1000)

# Load naca0012 data
print 'Naca 0012 data'
x,y = loadtxt('naca0012.dat',unpack=True)
curve = pySpline.curve('lms',x=x,y=y,k=4,Nctl=19)
curve.runParameterCorrection(10000)
