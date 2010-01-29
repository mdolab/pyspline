from numpy import *
import pySpline

# Get some Helix-like data

n = 100
theta = linspace(0,2*pi,n)
x = cos(theta)
y = sin(theta)
z = linspace(0,1,n)

curve = pySpline.curve('lms',x=x,y=y,z=z,k=4,Nctl=10)
