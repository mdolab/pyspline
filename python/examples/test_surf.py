from numpy import *
import pySpline
import sys
nu = 25
nv = 25
u = linspace(0,4,nu)
v = linspace(0,4,nv)
[V,U] = meshgrid(v,u)
Z = cos(U)*sin(V)

surf = pySpline.surface('lms',Nctlu=8,Nctlv=8,x=U,y=V,z=Z,ku=4,kv=4)
surf.writeTecplot('test.dat')
sys.exit(0)
u = 0.5
v = 0.5

u_vec = linspace(0,1,10)
v_vec = linspace(0,1,10)

[v_mat,u_mat] = meshgrid(v_vec,u_vec)

# Check the getvalues
surf(u,v)
surf(u_vec,v_vec)
surf(u_mat,v_mat)

surf.getDerivative(u,v)
surf.getDerivative(u_vec,v_vec)
surf.getDerivative(u_mat,v_mat)

surf.getSecondDerivative(u,v)
surf.getSecondDerivative(u_vec,v_vec)
surf.getSecondDerivative(u_mat,v_mat)


# x0 = surf(u,v)
# x0 = [2,2,0]
# u,v,diff = surf.projectPoint(x0)
# print 'u,v,diff'
# print u,v,diff
