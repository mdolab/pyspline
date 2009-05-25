#/usr/bin/python

# =============================================================================
# Standard Python modules
# =============================================================================
import sys,os,string,time

# =============================================================================
# Extension modules
# =============================================================================
#Numpy
import numpy
from numpy import array,zeros,sin,cos,array,linspace,pi,meshgrid,ones,where

#Scipy
import scipy
from scipy import interpolate

#pySpline
sys.path.append("../")
import pyspline as spline
#import pyspline_cs as spline_cs
from matplotlib.pylab import *
#Functions to load airfoil file

def load_af(filename):
    ''' Load the airfoil file in precomp format'''
    f = open(filename,'r')

    new_aux = string.split(f.readline())
    naf = int(new_aux[0]) 

    xnodes = zeros(naf,float)
    ynodes = zeros(naf,float)

    f.readline()
    f.readline()
    f.readline()

    for i in xrange(naf):

        new_aux = string.split(f.readline()) # you can use also other separaters,

        xnodes[i] = float(new_aux[0])
        ynodes[i] = float(new_aux[1])
    
    f.close()
    return naf,xnodes, ynodes

# =============================================================================
# Start of Script
# =============================================================================

# =================================
# Setup the blade
# =================================

#Lets start setting things we know we will need
naf = 12
bl_length = 6.15

# Wind Turbine Blade Example
chord = [0.6640,0.6640,0.6640,0.6640,0.6640,0.6640,1.0950,1.6800,\
         1.5390,1.2540,0.9900,0.7900,0.6100,0.4550,0.4540,0.4530]
sloc = [0.0000,0.0236,0.0273,0.0333,0.0393,0.0397,0.1141,\
        0.2184,0.3226,0.4268,0.5310,0.6352,0.7395,0.8437,0.9479,1.0000]
le_loc = [0.5000,0.5000,0.5000,0.5000,0.5000,0.5000,0.3300,\
          0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500]

airfoil_list = ['af1-6.inp','af1-6.inp','af1-6.inp','af1-6.inp','af1-6.inp','af1-6.inp',\
                'af-07.inp','af8-9.inp','af8-9.inp','af-10.inp','af-11.inp','af-12.inp','af-13.inp',\
                'af-14.inp','af15-16.inp','af15-16.inp']

#tw_aero = [-8.0000,0.0000,8.0,16.0,20.3900,16.0200,11.6500,\
         #  6.9600,1.9800,-1.8800,-3.3700,-3.4100,-3.4500,-3.4700]

chord = [0.6640,.6440,1.0950,1.6800,\
         1.5390,1.2540,0.9900,0.7900,0.6100,0.4550,0.4540,0.4530]
chord = 1.0*ones([naf])

sloc = [0.0000,0.0397,0.1141,\
        0.2184,0.3226,0.4268,0.5310,0.6352,0.7395,0.8437,0.9479,1.0000]
le_loc = [0.5000,0.5000,0.3300,\
          0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500]
le_loc = zeros([naf],float)
airfoil_list = ['af1-6.inp','af1-6.inp','af-07.inp',\
                    'af8-9.inp','af8-9.inp','af-10.inp','af-11.inp','af-12.inp',\
                    'af-13.inp', 'af-14.inp','af15-16.inp','af15-16.inp']

airfoil_list = ['af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp',\
                    'af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp','af15-16.inp']
#
max_nodes = 70
x_nodes = zeros([naf,max_nodes,2],float)
y_nodes = zeros([naf,max_nodes,2],float)
n_nodes = zeros([naf,2],int)

N = 10;
x = zeros([naf,N,2],float)
y = zeros([naf,N,2],float)
z = zeros([naf,N,2],float)

#x_interp = linspace(0,1,N)
theta = linspace(0,pi,N)

x_interp= 0.5*(1-cos(theta))
x_interp  = array([0,0.02,0.040,0.125,0.375,0.625,0.875,0.96,0.98,1],float)
print 'x_interp:',x_interp

# Add a custom knot vector (interior knots)
t  = array([0,0.02,0.040,0.125,0.375,0.625,0.875,0.96,0.98,1],float)
t = t[1:-2]
for i in xrange(naf):
    n_temp,temp_x,temp_y = load_af(airfoil_list[i])

    index = where(temp_x == 1)
    te_index = index[0][0]
    # upper Surface
    
    n_nodes[i,0] = te_index+1   # number of nodes on upper surface
    n_nodes[i,1] = int(n_temp-te_index) # number of nodes on lower surface

    x_u = temp_x[0:n_nodes[i,0]]
    y_u = temp_y[0:n_nodes[i,0]]

#    theta = linspace(0,pi,N)
#    t     = 0.5*(1-cos(theta))
   
    tck = interpolate.splrep(x_u,y_u,task=-1,k=3,t=t)
    y_interp   = interpolate.splev(x_interp,tck)

    x[i,:,0] = x_interp*chord[i] - le_loc[i]*chord[i]
    y[i,:,0] = y_interp*chord[i] 
    z[i,:,0] = sloc[i]*bl_length

    # lower Surface
    
    n_nodes[i,0] = te_index+1   # number of nodes on upper surface
    n_nodes[i,1] = n_temp-te_index # number of nodes on lower surface

    x_l = temp_x[te_index:n_temp]
    y_l = temp_y[te_index:n_temp]

    #Reverse coordinates for lower spline
    x_l = x_l[::-1]
    y_l = y_l[::-1]

#    theta = linspace(0,pi,N)
#    t     = 0.5*(1-cos(theta))
#    t = t[1:-2]

    tck = interpolate.splrep(x_l,y_l,task=-1,t=t)
    y_interp   = interpolate.splev(x_interp,tck)

    x[i,:,1] = x_interp*chord[i] - le_loc[i]*chord[i]
    y[i,:,1] = y_interp*chord[i] 
    z[i,:,1] = sloc[i]*bl_length

#end for

u = linspace(-1,1,naf)
#v = linspace(-1,1,N)
#theta = linspace(0,pi,N)
#v     = 1*(1-cos(theta))-1
# Add a custom knot vector (interior knots)
t = array([0,0.02,0.040,0.125,0.375,0.625,0.875,0.96,0.98,1],float)

v = 2*t - 1 
print 'v:',v

ku = 4
kv = 4
print 'Going to do the interpolation'
tu_x =[]; tv_x = []; bcoef_x = []
tu_y =[]; tv_y = []; bcoef_y = []
tu_z =[]; tv_z = []; bcoef_z = []
for isurf in xrange(2):
    tu,tv,bcoef = spline.b2ink(u,v,x[:,:,isurf],ku,kv)
    print 'tv:',tv
    tu_x.append(tu)
    tv_x.append(tv)
    bcoef_x.append(bcoef)

    tu,tv,bcoef = spline.b2ink(u,v,y[:,:,isurf],ku,kv)
    tu_y.append(tu)
    tv_y.append(tv)
    bcoef_y.append(bcoef)

    tu,tv,bcoef = spline.b2ink(u,v,z[:,:,isurf],ku,kv)
    tu_z.append(tu)
    tv_z.append(tv)
    bcoef_z.append(bcoef)

print 'x coef for station 6:'
print bcoef_x[1][6,:]
print 'x:'
print x[6,:,0]
plot(x[6,:,0],y[6,:,0],'ko-')
#plot(bcoef_y[0][6,:],'ko--')
show()

# end for
#  b2val = b2val(xval,yval,idx,idy,tx,ty,kx,ky,bcoef,nx=shape(bcoef,0),ny=shape(bcoef,1),work=)

 #Puturb a coefficient in y

# # Do a "dihedreal" motion

# tip_def = 2 # meters

# local_def = array(sloc)*tip_def;

# for i in xrange(naf):
#     bcoef_y[0][i,:] += local_def[i]
#     bcoef_y[1][i,:] += local_def[i]


# bl_stretch = 1.25;

# bcoef_z[0] *= bl_stretch
# bcoef_z[1] *= bl_stretch

# #increase chord by 20% at station 6

# print 'xcoef:',bcoef_x[0][6,:]
# print 'xcoef:',bcoef_x[0][6,:]

# bcoef_x[0][6,:] *= 1.2
# bcoef_x[1][6,:] *= 1.2

# bcoef_x[0][6,:] += 0.05
# bcoef_x[1][6,:] += 0.05

# bcoef_x[0][7,:] += 0.1
# bcoef_x[1][7,:] += 0.1

# bcoef_x[0][8,:] += 0.2
# bcoef_x[1][8,:] += 0.2

# bcoef_x[0][9,:] += 0.4
# bcoef_x[1][9,:] += 0.4

# bcoef_x[0][10,:] += 0.8
# bcoef_x[1][10,:] += 0.8

# bcoef_x[0][10,:] += 1.6
# bcoef_x[1][10,:] += 1.6

# bcoef_x[0][11,:] += 3.2
# bcoef_x[1][11,:] += 3.2

naf = 60
n_nodes = 20

u = linspace(-1,1,naf)
v = linspace(-1,1,n_nodes)

theta = linspace(0,pi,n_nodes)
v = 1*(1-cos(theta))-1
print v

#Start tecplot output

nodes_total = naf*n_nodes
elements_total = (naf-1)*(n_nodes-1)

f = open("output.dat",'w')
points = zeros([naf,n_nodes,2,3],float) # section, nodes,isurf, [x,y,z]

f.write ('\"Blade Data\"\n')
f.write ('VARIABLES = "X", "Y","Z"\n')

for isurf in xrange(2):
    f.write('Zone N=%d, E=%d\n'%(nodes_total,elements_total))
    f.write('DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL\n')
    for i in xrange(naf):
        for j in xrange(n_nodes):
            points[i,j,isurf,0] = spline.b2val(u[i],v[j],0,0,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])
            points[i,j,isurf,1] = spline.b2val(u[i],v[j],0,0,tu_y[isurf],tv_y[isurf],ku,kv,bcoef_y[isurf])
            points[i,j,isurf,2] = spline.b2val(u[i],v[j],0,0,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])
        # end for
    # end for

    # The next step will be to output all the x-y-z Data
    for i in xrange(naf):
        for j in xrange(n_nodes):
            f.write("%f %f %f\n"%(points[i,j,isurf,0],points[i,j,isurf,1],points[i,j,isurf,2]))
        # end for
    # end for
                                 
           
    # now write out the connectivity
    for i in xrange(naf-1):
        for j in xrange(n_nodes-1):
            f.write( '%d %d %d %d\n'%(i*n_nodes + (j+1), i*n_nodes+(j+2), (i+1)*n_nodes+(j+2),(i+1)*n_nodes + (j+1)))
        # end for
    #end for
# end for

# Also dump out the control points
naf = 12 #original naf
N = 10   #original N

for isurf in xrange(2):
    f.write('Zone I=%d, J=%d\n'%(naf,N))
    f.write('DATAPACKING=POINT\n')
    for j in xrange(N):
        for i in xrange(naf):
            f.write("%f %f %f \n"%(bcoef_x[isurf][i,j],\
                                   bcoef_y[isurf][i,j],\
                                   bcoef_z[isurf][i,j]))
        # end for
    # end for 
# end for
f.close()
