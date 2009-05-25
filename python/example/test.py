#/usr/bin/python

# =============================================================================
# Standard Python modules
# =============================================================================
import sys,os,string,time

# =============================================================================
# Extension modules
# =============================================================================
#Numpy
from numpy import array,zeros,sin,cos,array,linspace,pi,meshgrid,ones,where,\
    interp, append, hstack,sqrt,mat

import scipy.linalg

#pySpline
sys.path.append("../")
import pyspline as spline
import pyspline_cs as spline_cs

#Function to load airfoil file (precomp format)

def load_af(filename):
    ''' Load the airfoil file in precomp format'''
    f = open(filename,'r')

    new_aux = string.split(f.readline())
    naf = int(new_aux[0]) 

    xnodes = zeros(naf)
    ynodes = zeros(naf)

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
naf = 11
bl_length = 21.15

# Wind Turbine Blade Example

tw_aero = [-8.0000,0.0000,8.0,16.0,20.3900,16.0200,11.6500,\
                6.9600,1.9800,-1.8800,-3.3700,-3.4100,-3.4500,-3.4700]

chord = [.6440,1.0950,1.6800,\
         1.5390,1.2540,0.9900,0.7900,0.6100,0.4550,0.4540,0.4530]

sloc = [0.0000,0.1141,\
        0.2184,0.3226,0.4268,0.5310,0.6352,0.7395,0.8437,0.9479,1.0000]
le_loc = [0.5000,0.3300,\
          0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500,0.2500]

airfoil_list = ['af1-6.inp','af-07.inp',\
                    'af8-9.inp','af8-9.inp','af-10.inp','af-11.inp','af-12.inp',\
                    'af-13.inp', 'af-14.inp','af15-16.inp','af15-16.inp']


n_nodes = zeros((naf,2),int)
N = 15;
x = zeros((naf,N*2-1,2),float)
y = zeros((naf,N*2-1,2),float)
z = zeros((naf,N*2-1,2),float)

# This is the standard cosine distribution in x (one sided)
x_interp = 0.5*(1-cos(linspace(0,pi,N)))
#x_interp = linspace(0,1,N)

for i in xrange(naf):
    n_temp,temp_x,temp_y = load_af(airfoil_list[i])
    
    # Find the trailing edge point
    index = where(temp_x == 1)
    te_index = index[0][0]
    n_nodes[i,0] = te_index+1   # number of nodes on upper surface
    n_nodes[i,1] = int(n_temp-te_index) # number of nodes on lower surface
    
    # upper Surface
    x_u = temp_x[0:n_nodes[i,0]]
    y_u = temp_y[0:n_nodes[i,0]]
    
    # linearly interpolate to find the points at the positions we want
    y_interp_u = interp(x_interp,x_u,y_u)
    
    #Reverse upper coordinates to start at te
    y_interp_u = y_interp_u[::-1]
    x_interp_u = x_interp[::-1]
    
    # lower Surface
    
    n_nodes[i,0] = te_index+1   # number of nodes on upper surface
    n_nodes[i,1] = n_temp-te_index # number of nodes on lower surface Add extra 1 for extra le
    
    x_l = temp_x[te_index:n_temp]
    y_l = temp_y[te_index:n_temp]

     # Reverse coordinates for lower spline
    x_l = x_l[::-1]
    y_l = y_l[::-1]

    x_interp_l = x_interp[1:] #DO NOT want 0,0 in lower surface

    # Interpolate
    y_interp_l = interp(x_interp_l,x_l,y_l)


    x_cor_full = hstack([x_interp_u,x_interp_l])
    y_cor_full = hstack([y_interp_u,y_interp_l])
    
    # Finally Set the Coordinates
    x[i,:,0] = x_cor_full*chord[i] - le_loc[i]*chord[i]
    y[i,:,0] = y_cor_full*chord[i] 
    z[i,:,0] = sloc[i]*bl_length

#end for

u = linspace(-1,1,naf)
#v = -1+(1-cos(linspace(0,pi,N)))
v = linspace(-1,1,2*N-1)

ku = 4
kv = 4
print 'Going to do the interpolation'

tu_x =[]; tv_x = []; bcoef_x = []
tu_y =[]; tv_y = []; bcoef_y = []
tu_z =[]; tv_z = []; bcoef_z = []
for isurf in xrange(1):
    tu,tv,bcoef = spline.b2ink(u,v,x[:,:,isurf],ku,kv)
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
#end for 

#Now re-interpolate for the tecplot output

u_plot = 150
v_plot = 31

#v = -1+(1-cos(linspace(0,2*pi,v_plot)))

u = linspace(-1,1,u_plot)
v = linspace(-1,1,v_plot)
#Start tecplot output

nodes_total = u_plot*v_plot
elements_total = (u_plot-1)*(v_plot-1) 

f = open("output.dat",'w')
points = zeros((u_plot,v_plot,2,3),float) # section, nodes,isurf, [x,y,z]

f.write ('\"Blade Data\"\n')
f.write ('VARIABLES = "X", "Y","Z"\n')

for isurf in xrange(1):
    f.write('Zone N=%d, E=%d\n'%(nodes_total,elements_total))
    f.write('DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL\n')
    for i in xrange(u_plot):
        for j in xrange(v_plot):
            points[i,j,isurf,0] = spline.b2val(u[i],v[j],0,0,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])
            points[i,j,isurf,1] = spline.b2val(u[i],v[j],0,0,tu_y[isurf],tv_y[isurf],ku,kv,bcoef_y[isurf])
            points[i,j,isurf,2] = spline.b2val(u[i],v[j],0,0,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])
        # end for
    # end for

    # The next step will be to output all the x-y-z Data
    for i in xrange(u_plot):
        for j in xrange(v_plot):
            f.write("%f %f %f\n"%(points[i,j,isurf,0],points[i,j,isurf,1],points[i,j,isurf,2]))
        # end for
    # end for
                                 
           
    # now write out the connectivity
    for i in xrange(u_plot-1):
        for j in xrange(v_plot-1):
            f.write( '%d %d %d %d\n'%(i*v_plot + (j+1), i*v_plot+(j+2), (i+1)*v_plot+(j+2),(i+1)*v_plot + (j+1)))
        # end for
    #end for
# end for

# Also dump out the control points

for isurf in xrange(1):
    f.write('Zone I=%d, J=%d\n'%(naf,N*2-1))
    f.write('DATAPACKING=POINT\n')
    for j in xrange(2*N-1):
        for i in xrange(naf):
            f.write("%f %f %f \n"%(bcoef_x[isurf][i,j],\
                                   bcoef_y[isurf][i,j],\
                                   bcoef_z[isurf][i,j]))
        # end for
    # end for 
# end for
f.close()
 




# Test of a newton iteration to find a location on the surface given the gloabl x-z coordiante

# For example, find the u-v coordinates for the upper and lower surfaces for x=0,z=10

x_o = 0
z_o = 10
# Note
# te = -1, le =0, te = 1

# We *should* always get a good guess 

u_o = (z_o/bl_length)*2 -1 #This should be excellent
v_o = .8  # should be ok


print 'u_o:',u_o
print 'v_o:',v_o

#eval
isurf = 0
x = spline.b2val(u_o,v_o,0,0,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])
y = spline.b2val(u_o,v_o,0,0,tu_y[isurf],tv_y[isurf],ku,kv,bcoef_y[isurf])
z = spline.b2val(u_o,v_o,0,0,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])

print 'x,y,z:',x,y,z

u = u_o
v = v_o

for iter in xrange(12):


    x_iter = spline.b2val(u,v,0,0,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])
    y_iter = spline.b2val(u,v,0,0,tu_y[isurf],tv_y[isurf],ku,kv,bcoef_y[isurf])
    z_iter = spline.b2val(u,v,0,0,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])

    f = mat([x_iter-x_o,z_iter-z_o])

    #now get the derivative

    dxdu = spline.b2val(u,v,1,0,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])
    dxdv = spline.b2val(u,v,0,1,tu_x[isurf],tv_x[isurf],ku,kv,bcoef_x[isurf])

    dzdu = spline.b2val(u,v,1,0,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])
    dzdv = spline.b2val(u,v,0,1,tu_z[isurf],tv_z[isurf],ku,kv,bcoef_z[isurf])

    
    J = mat([[dxdu,dxdv],[dzdu,dzdv]])
    print 'J:',J
    J_inv = scipy.linalg.inv(mat(J))
    #Now update the u and v
    
    x_up = J_inv*f.T

    u = u - 0.75*x_up[0]
    v = v - 0.75*x_up[1]
    print 'f:',f
    print "Itertion %d, u=%f,v=%f\n"%(iter,u,v)

print x_iter,y_iter,z_iter
