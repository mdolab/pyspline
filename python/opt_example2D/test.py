# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot, cross
from scipy.io import read_array
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline

# =============================================================================
# Extension modules
# =============================================================================

# pySpline
sys.path.append('../../python')

#pyOPT
sys.path.append(os.path.abspath('../../../../pyACDT/pyACDT/Optimization/pyOpt/'))

#pySNOPT
sys.path.append(os.path.abspath('../../../../pyACDT/pyACDT/Optimization/pyOpt/pySNOPT'))

#pySpline
import pySpline

'''This script runs an example of fitting a cubic spline to a set of
3D wing coordinates that were generated using aerosurf. The script is
meant to show how to fit a constrained surface using a non-linear
optimizer'''

Nu = 26
Nv = 25

# Global Variables

# ==============================
# Read in the Data 
# ==============================

# X holds all the coordinates
# First index is isurf -> surface id
# Second index is U
# Third index is V
# Fourth index is 0,1,2 for x,y,z

Nsurf = 2
X = zeros([Nsurf,Nu,Nv,3])

# Top surface of wing
print 'reading upper surface...'
data = read_array(file("upper_surface.inp"))

for j in xrange(Nv):
    for i in xrange(Nu):
        X[0,i,j,:] = data[j*Nu + i ,:]
    #end for
#end for

# Bottom surface of wing
print 'reading bottom surface...'
data = read_array(file("lower_surface.inp"))
for j in xrange(Nv):
    for i in xrange(Nu):
        X[1,i,j,:] = data[j*Nu + i, :]
    #end for
#end for

# ==============================
# Find the 'u' parameters
# ==============================

u=zeros((Nsurf,Nu))
for isurf in xrange(Nsurf):
    x0=X[isurf,:,0,0];y0=X[isurf,:,0,1];z0=X[isurf,:,0,2]
    for i in xrange(Nu-1):
        u[isurf,i+1] = u[isurf,i] + sqrt((x0[i+1]-x0[i])**2 + \
                                   (y0[i+1]-y0[i])**2+ \
                                   (z0[i+1]-z0[i])**2)
    # end for
    u[isurf,:] = u[isurf,:]/u[isurf,-1] #Normalize u
# end for

# ==============================
# Find the 'v' parameters
# ==============================
v=zeros((Nsurf,Nv))
for isurf in xrange(Nsurf):
    x0=X[isurf,0,:,0];y0=X[isurf,0,:,1];z0=X[isurf,0,:,2]
    for i in xrange(Nv-1):
        v[isurf,i+1] = v[isurf,i] + sqrt((x0[i+1]-x0[i])**2 + \
                                   (y0[i+1]-y0[i])**2+ \
                                   (z0[i+1]-z0[i])**2)
    # end for
    v[isurf,:] = v[isurf,:]/v[isurf,-1] #Normalize u

spline_test = pySpline.spline(Nsurf,u,v,X,fit_type='lms',Nctlu = 13, Nctlv = 7, ku=4,kv=4)

print 'Done!'

spline_test.writeTecplot('output.dat')

print spline_test.getValue(0,0.5,0.5)
