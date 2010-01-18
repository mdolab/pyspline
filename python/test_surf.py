from numpy import *
import pySpline

nu = 20
nv = 20
u = linspace(0,1,nu)
v = linspace(0,1,nv)
[V,U] = meshgrid(v,u)

Z = cos(U)*sin(V)

X = zeros(nu,nv,3)
X[:,:,0] = U
X[:,:,1] = V
X[:,:,2] = Z

#surf = pySpline.surface('lms',Nctlu=5,Nctlv=5,X=X,ku=4,kv=4)
#surf.writeTecplot('test.dat')
