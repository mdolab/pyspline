# This is a test script to test the functionality of the 
# pySpline surface

from numpy import *
import sys
sys.path.append('../')
import pySpline

# Create a generic surface
nu = 20
nv = 20
u = linspace(0,4,nu)
v = linspace(0,4,nv)
[V,U] = meshgrid(v,u)
Z = cos(U)*sin(V)
surf = pySpline.surface(x=U,y=V,z=Z,ku=4,kv=4,Nctlu=5,Nctlv=5)
surf.writeTecplot('surface.dat')

n = 100
theta = linspace(.0000,2*pi,n)
x = (cos(theta)-1)
y = (sin(theta)+1)
z = (linspace(0,1,n) + 2)
curve = pySpline.curve(x=x,y=y,z=z,k=4,Nctl=16,niter=100)
curve.writeTecplot('helix.dat',size=.1)

u,v,s,D = surf.projectCurve(curve,Niter=100,eps1=1e-10,eps2=1e-10,
                            u=1,v=1,s=1)

print u,v,s,D
print curve(s)
print surf(u,v)



# # Test the project point Algorithim
# x0 = [3,3,2]

# u,v,D = surf.projectPoint(x0,Niter=100)

# val = surf(u,v)
# # Output the data
# f = open('projections2.dat','w')
# f.write ('VARIABLES = "X", "Y","Z"\n')
# f.write('Zone T=surf_proj_pt I=2 \n')
# f.write('DATAPACKING=POINT\n')
# f.write('%f %f %f\n'%(x0[0],x0[1],x0[2]))
# f.write('%f %f %f\n'%(val[0],val[1],val[2]))

# # Test the project surface-curve Algorithim

# x = [0,1,2]
# y = [4,3,2]
# z = [-3,1,3]

# curve = pySpline.curve(k=3,x=x,y=y,z=z)
# curve.writeTecplot('curve3.dat',size=.2)
# u,v,s,D = surf.projectCurve(curve)
# val1 = surf(u,v)
# val2 = curve(s)

# f.write('Zone T=surf_proj_curve I=2 \n')
# f.write('DATAPACKING=POINT\n')
# f.write('%f %f %f\n'%(val1[0],val1[1],val1[2]))
# f.write('%f %f %f\n'%(val2[0],val2[1],val2[2]))

# f.close()
