from numpy import *
import sys
from matplotlib.pylab import plot,show
from pyspline import *

#Debugging script for testing derivative  constrained 2D b-spline surfaces
Nu = 6
Nv = 12
x0 = array([0,0.000001,0.0325,0.05,0.1,0.2,0.35,0.5,0.65,0.85,0.925,1])
#x0 = linspace(0,1,Nv)
#Top surface of a naca0012 airfoil
y0= 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)

#create v
v = zeros(Nv)

for i in xrange(Nv-1):
    v[i+1] = v[i] + sqrt((x0[i+1]-x0[i])**2 + (y0[i+1]-y0[i])**2)

v = v/v[-1] #Normalize s

#Create u 
print 'v:',v
u = linspace(0,1,Nu)  #Just do a straight wing

#Create the surface points
ku = 4
kv = 4

x = zeros((Nu,Nv))
y = zeros((Nu,Nv))
z = zeros((Nu,Nv))

for i in xrange(Nu):
    x[i,:] = x0
    y[i,:] = y0
    z[i,:] = 5*u[i]

print 'x:',x
print 'y:',y
print 'z:',z

print 'calling b2ink'
tu_x,tv_x,bcoef_x = b2ink(u,v,x,ku,kv)
tu_y,tv_y,bcoef_y = b2ink(u,v,y,ku,kv)
tu_z,tv_z,bcoef_z = b2ink(u,v,z,ku,kv)
print 'bcoef_x:',bcoef_x
bcoef_x[:,1] = 0
print 'bcoef_x:',bcoef_x
#print 'bcoef_y:',bcoef_y
#print 'bcoef_z:',bcoef_z



#val= b2val(0.5,0.5,0,0,tu_x,tv_x,4,4,bcoef_x)

#print 'calling b2ink_mod:'
#tu_x,tv_x,bcoef_x = b2ink_mod(v,u,x.T,ku,kv)

#Very quick tecplot output
Nu = 10
Nv = 1000
u = linspace(0,1,Nu)
v = linspace(0,1,Nv)
f = open('output.dat','w')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(b2val(u[i],v[j],0,0,tu_x,tv_x,4,4,bcoef_x),
                               b2val(u[i],v[j],0,0,tu_y,tv_y,4,4,bcoef_y),
                               b2val(u[i],v[j],0,0,tu_z,tv_z,4,4,bcoef_z)))
        #f.write('%f %f %f \n'%(x[i,j],y[i,j],z[i,j]))
f.close()
