from numpy import *
import sys
from matplotlib.pylab import plot,show
from pyspline import *
from scipy.io import read_array,savemat

data = read_array(file("wing.inp"))

Nu = 26
Nv = 25

x0 = zeros((Nu,Nv))
x1 = zeros((Nu,Nv))
x2 = zeros((Nu,Nv))
w  = ones((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        x0[i,j] = data[Nu*j+i,0]
        x1[i,j] = data[Nu*j+i,1]
        x2[i,j] = data[Nu*j+i,2]
    #end for
#end for

#Find the 'u' parameter
x=x0[:,0];y=x1[:,0];z=x2[:,0] 
u = zeros(Nu)
for i in xrange(Nu-1):
    u[i+1] = u[i] + sqrt((x[i+1]-x[i])**2 + \
                             (y[i+1]-y[i])**2+ \
                             (z[i+1]-z[i])**2)
#end for
u = u/u[-1] #Normalize u

#Find the 'v' parameter

x=x0[0,:];y=x1[0,:];z=x2[0,:] 
v = zeros(Nv)
for i in xrange(Nv-1):
    v[i+1] = v[i] + sqrt((x[i+1]-x[i])**2 + \
                               (y[i+1]-y[i])**2+ \
                               (z[i+1]-z[i])**2)
#end for
v = v/v[-1] #Normalize v

[U,V] = meshgrid(v,u)

#Options
iopt = -1 #least squares surface
ku = 3
kv = 3
nu = 18
nv = 18
smooth = 5e-5

#Knot Stuff

m = len(x0.flatten())
nuest=ku+1+floor(sqrt(m/2))
nvest=kv+1+floor(sqrt(m/2))
# do a uniform knot vector
tu = zeros(nuest)
tu[0:ku+1] = u[0]
tu[nu-ku-1:nu] = u[-1]
tu[ku:nu-ku] = 0.5*(1-cos(linspace(0,pi,nu-2*ku)))

tv= zeros(nvest)
tv[0:kv+1] = v[0]
tv[nv-kv-1:nv] = v[-1]
tv[kv:nv-kv]=0.5*(1-cos(linspace(0,pi,nv-2*kv)))

nx,tx,ny,ty,cx1,fpx,ier = surfit(iopt,U.flatten(),V.flatten(),x0.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
nu,tu,nv,tv,cy1,fpy,ier = surfit(iopt,U.flatten(),V.flatten(),x1.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
nu,tu,nv,tv,cz1,fpz,ier = surfit(iopt,U.flatten(),V.flatten(),x2.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
print 'fp x,y,z:',fpx,fpy,fpz

Nu_plot = 100
Nv_plot = 25
Nu = nu-ku-1
Nv = nv-kv-1
uplot = linspace(-1,1,Nu_plot)
vplot = linspace(-1,1,Nv_plot)
vplot = 0.5*(1-cos(linspace(0,pi,Nv_plot)))
# Start tecplot output

f = open('output3.dat','w')
f.write('Zone I=%d J = %d\n'%(Nu_plot,Nv_plot))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv_plot):
    for i in xrange(Nu_plot):
        xo,ier = bispev(tu[0:nu],tv[0:nv],cx1[0:(nu-kv-1)*(nv-kv-1)],ku,kv,uplot[i],vplot[j])
        yo,ier = bispev(tu[0:nu],tv[0:nv],cy1[0:(nu-kv-1)*(nv-kv-1)],ku,kv,uplot[i],vplot[j])
        zo,ier = bispev(tu[0:nu],tv[0:nv],cz1[0:(nu-kv-1)*(nv-kv-1)],ku,kv,uplot[i],vplot[j])
        f.write('%f %f %f \n'%(xo,yo,zo))
        
 # Also dump out the control points
f.write('Zone I=%d, J=%d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write("%.5g %.5g %.5g \n"%(cx1[i*Nv + j],cy1[i*Nv + j],cz1[i*Nv + j]))

# ------------------------------------
# Now do the FUSE section
# ------------------------------------
data = read_array(file("fuse.inp"))

Nu = 26
Nv = 26

x0 = zeros((Nu,Nv))
x1 = zeros((Nu,Nv))
x2 = zeros((Nu,Nv))
w  = ones((Nu,Nv))
for j in xrange(Nv):
    for i in xrange(Nu):
        x0[i,j] = data[Nu*j+i,0]
        x1[i,j] = data[Nu*j+i,1]
        x2[i,j] = data[Nu*j+i,2]
    #end for
#end for

#Find the u parameter : Must be (and is) Identical to one above

x=x0[:,0];y=x1[:,0];z=x2[:,0] 
u = zeros(Nu)
for i in xrange(Nu-1):
    u[i+1] = u[i] + sqrt((x[i+1]-x[i])**2 + \
                               (y[i+1]-y[i])**2+ \
                               (z[i+1]-z[i])**2)
#end for
u = u/u[-1] #Normalize u

#Find the 'v' parameter

x=x0[0,:];y=x1[0,:];z=x2[0,:] 
v = zeros(Nv)
for i in xrange(Nv-1):
    v[i+1] = v[i] + sqrt((x[i+1]-x[i])**2 + \
                               (y[i+1]-y[i])**2+ \
                               (z[i+1]-z[i])**2)
#end for
v = v/v[-1] #Normalize v

#mesh grid it
[U,V] = meshgrid(v,u)

#Options
m = len(x0.flatten())
nuest=ku+1+floor(sqrt(m/2))
nvest=kv+1+floor(sqrt(m/2))
# do a uniform knot vector
tu = zeros(nuest)
tu[0:ku+1] = u[0]
tu[nu-ku-1:nu] = u[-1]
tu[ku:nu-ku] = 0.5*(1-cos(linspace(0,pi,nu-2*ku)))

tv= zeros(nvest)
tv[0:kv+1] = v[0]
tv[nv-kv-1:nv] = v[-1]
tv[kv:nv-kv]=0.5*(1-cos(linspace(0,pi,nv-2*kv)))

nx,tx,ny,ty,cx2,fpx,ier = surfit(iopt,U.flatten(),V.flatten(),x0.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
nu,tu,nv,tv,cy2,fpy,ier = surfit(iopt,U.flatten(),V.flatten(),x1.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
nu,tu,nv,tv,cz2,fpz,ier = surfit(iopt,U.flatten(),V.flatten(),x2.flatten(),w.flatten()\
                                    ,0,1,0,1,ku,kv,smooth,nuest,nvest,nu,tu,nv,tv)
print 'fp x,y,z:',fpx,fpy,fpz

#Make the control points on edge match
cx2[0:nv-kv-1] = cx1[0:nv-kv-1]
cy2[0:nv-kv-1] = cy1[0:nv-kv-1]
cz2[0:nv-kv-1] = cz1[0:nv-kv-1]

Nu = nu-ku-1 
Nv = nv-kv-1
u = linspace(-1,1,Nu_plot)
v = linspace(-1,1,Nv_plot)
v = 0.5*(1-cos(linspace(0,pi,Nv_plot)))
# Start tecplot output

f.write('Zone I=%d J = %d\n'%(Nu_plot,Nv_plot))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv_plot):
    for i in xrange(Nu_plot):
        xo,ier = bispev(tu[0:nu],tv[0:nv],cx2[0:(nu-kv-1)*(nv-kv-1)],ku,kv,u[i],v[j])
        yo,ier = bispev(tu[0:nu],tv[0:nv],cy2[0:(nu-kv-1)*(nv-kv-1)],ku,kv,u[i],v[j])
        zo,ier = bispev(tu[0:nu],tv[0:nv],cz2[0:(nu-kv-1)*(nv-kv-1)],ku,kv,u[i],v[j])
        f.write('%f %f %f \n'%(xo,yo,zo))
                

 # Also dump out the control points
f.write('Zone I=%d, J=%d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write("%.5g %.5g %.5g \n"%(cx2[i*Nv + j],cy2[i*Nv + j],cz2[i*Nv + j]))

f.close()

