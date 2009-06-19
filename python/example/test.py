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

# =============================================================================
# Extension modules
# =============================================================================
import pyspline 
import pyspline_cs

#pyOPT
sys.path.append(os.path.abspath('../../../../pyACDT/pyACDT/Optimization/pyOpt/'))
from pyOpt_optimization import Optimization

#pySNOPT
sys.path.append(os.path.abspath('../../../../pyACDT/pyACDT/Optimization/pyOpt/pySNOPT'))
from pySNOPT import SNOPT

import pyspline
import pyspline_cs


#pySpline
sys.path.append("../")
import pySpline 

print 'sys.path1:',sys.path

# Define an analytic function

u = linspace(-1,1,10)
v = linspace(-1,1,10)

[U,V] = meshgrid(u,v)

# create the surface object

surface = pySpline.spline(u,v,x,y,z)

# Now re-interpolate for the tecplot output

u_plot = 150
v_plot = 31

u = linspace(-1,1,u_plot)
v = linspace(-1,1,v_plot)
# Start tecplot output

nodes_total = u_plot*v_plot
elements_total = (u_plot-1)*(v_plot-1) 

f = open("output.dat",'w')
points = zeros((u_plot,v_plot,3),float) # u,v, [x,y,z]

f.write ('\"pySpline Data\"\n')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone N=%d, E=%d\n'%(nodes_total,elements_total))
f.write('DATAPACKING=POINT,ZONETYPE=FEQUADRILATERAL\n')
for i in xrange(u_plot):
    for j in xrange(v_plot):
        points[i,j,:] = surface.getValue(u[i],v[j])
    # end for
# end for

# The next step will be to output all the x-y-z Data
for i in xrange(u_plot):
    for j in xrange(v_plot):
        f.write("%5g %5g %5g\n"%(points[i,j,0],points[i,j,1],points[i,j,2]))
    # end for
# end for
                                 
# now write out the connectivity
for i in xrange(u_plot-1):
    for j in xrange(v_plot-1):
        f.write( '%d %d %d %d\n'%(i*v_plot + (j+1), i*v_plot+(j+2), \
                                (i+1)*v_plot+(j+2),(i+1)*v_plot + (j+1)))
    # end for
#end for

# Also dump out the control points
f.write('Zone I=%d, J=%d\n'%(10,10))
f.write('DATAPACKING=POINT\n')
for j in xrange(10):
    for i in xrange(10):
        f.write("%5g %5g %5g \n"%(surface.bcoef_x[i,j],\
                                   surface.bcoef_y[i,j],\
                                   surface.bcoef_z[i,j]))
    # end for
# end for 
f.close()

print 'Now run Tecplot and open the default.lay file"
