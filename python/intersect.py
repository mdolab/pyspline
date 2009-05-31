from numpy import zeros,sqrt,linspace,cos,sin,pi
import sys
from pyspline import b2ink,b2val
from scipy.io import read_array
import time
data = read_array(file("wing.inp"))

Nu = 26
Nv = 25

x = zeros((Nu,Nv))
y = zeros((Nu,Nv))
z = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        x[i,j] = data[Nu*j+i,0]
        y[i,j] = data[Nu*j+i,1]
        z[i,j] = data[Nu*j+i,2]
    #end for
#end for

#Find the 'u' parameter

x0=x[:,0];y0=y[:,0];z0=z[:,0] 
u = zeros(Nu)
for i in xrange(Nu-1):
    u[i+1] = u[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
u = u/u[-1] #Normalize u

#Find the 'v' parameter

x0=x[0,:];y0=y[0,:];z0=z[0,:] 
v = zeros(Nv)
for i in xrange(Nv-1):
    v[i+1] = v[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
v = v/v[-1] #Normalize v

ku=4
kv=4

tu_x,tv_x,bcoef_x = b2ink(u,v,x,ku,kv)
tu_y,tv_y,bcoef_y = b2ink(u,v,y,ku,kv)
tu_z,tv_z,bcoef_z = b2ink(u,v,z,ku,kv)

# print 'wing coef:'
# print bcoef_x[:,0]
# print bcoef_y[:,0]
# print bcoef_z[:,0]

#Do a linear twist

# xtemp = bcoef_x[:,-1]-7
# ytemp = bcoef_y[:,-1]-1

# theta = 5.0* pi/180;
# print 'xtemp:',xtemp
# print 'ytemp:',ytemp
# for i in xrange(Nu):
#     bcoef_x[i,-1] = cos(theta)*xtemp[i] - sin(theta)*ytemp[i]+7
#     bcoef_y[i,-1] = sin(theta)*xtemp[i] + cos(theta)*ytemp[i]+1




#bcoef_z[5:15,0] += 0.08
#bcoef_y[5:15,0] += 0.03

Nu_plot= 26
Nv_plot = 25

u = linspace(0,1,Nu_plot)
u = 0.5*(1-cos(linspace(0,pi,Nu_plot)))
v = linspace(0,1,Nv_plot)
f = open('output.dat','w')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone I=%d J = %d\n'%(Nu_plot,Nv_plot))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv_plot):
    for i in xrange(Nu_plot):
        f.write('%f %f %f \n'%(b2val(u[i],v[j],0,0,tu_x,tv_x,4,4,bcoef_x),
                               b2val(u[i],v[j],0,0,tu_y,tv_y,4,4,bcoef_y),
                               b2val(u[i],v[j],0,0,tu_z,tv_z,4,4,bcoef_z)))
        
        #f.write('%f %f %f \n'%(x[i,j],y[i,j],z[i,j]))

 # Also dump out the control points
f.write('Zone I=%d, J=%d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write("%.5g %.5g %.5g \n"%(bcoef_x[i,j],\
                                     bcoef_y[i,j],\
                                     bcoef_z[i,j]))

# ------------------------------------
# Now do the FUSE section
# ------------------------------------
data = read_array(file("fuse.inp"))

Nu = 26
Nv = 26

x = zeros((Nu,Nv))
y = zeros((Nu,Nv))
z = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        x[i,j] = data[Nu*j+i,0]
        y[i,j] = data[Nu*j+i,1]
        z[i,j] = data[Nu*j+i,2]
    #end for
#end for

#Find the u parameter

x0=x[:,0];y0=y[:,0];z0=z[:,0] 
u = zeros(Nu)
for i in xrange(Nu-1):
    u[i+1] = u[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
u = u/u[-1] #Normalize u

#Find the 'v' parameter

x0=x[0,:];y0=y[0,:];z0=z[0,:] 
v = zeros(Nv)
for i in xrange(Nv-1):
    v[i+1] = v[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
v = v/v[-1] #Normalize v

tu_x,tv_x,bcoef_x = b2ink(u,v,x,ku,kv)
tu_y,tv_y,bcoef_y = b2ink(u,v,y,ku,kv)
tu_z,tv_z,bcoef_z = b2ink(u,v,z,ku,kv)

# print 'fuse coef:'
# print bcoef_x[:,0]
# print bcoef_y[:,0]
# print bcoef_z[:,0]

# bcoef_z[5:15,0] += 0.08
# bcoef_y[5:15,0] += 0.03

Nu_plot= 26
Nv_plot = 26
u = linspace(0,1,Nu_plot)
u = 0.5*(1-cos(linspace(0,pi,Nu_plot)))
v = linspace(0,1,Nv_plot)

f.write('Zone I=%d J = %d\n'%(Nu_plot,Nv_plot))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv_plot):
    for i in xrange(Nu_plot):
        f.write('%f %f %f \n'%(b2val(u[i],v[j],0,0,tu_x,tv_x,4,4,bcoef_x),
                               b2val(u[i],v[j],0,0,tu_y,tv_y,4,4,bcoef_y),
                               b2val(u[i],v[j],0,0,tu_z,tv_z,4,4,bcoef_z)))
                #        f.write('%f %f %f \n'%(x[i,j],y[i,j],z[i,j]))
                                     

 # Also dump out the control points
f.write('Zone I=%d, J=%d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write("%.5g %.5g %.5g \n"%(bcoef_x[i,j],\
                                         bcoef_y[i,j],\
                                         bcoef_z[i,j]))
        

f.close()
