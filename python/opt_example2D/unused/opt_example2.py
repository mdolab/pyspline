# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy, time

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, pi, hstack, zeros, ones, sqrt, imag, interp, \
    array, real, reshape, meshgrid, dot
from scipy.io import read_array
import scipy.interpolate
from scipy.interpolate import RectBivariateSpline

# =============================================================================
# Extension modules
# =============================================================================

# pySpline
sys.path.append('../../python')
import pyspline
import pyspline_cs

#pyOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/')
from pyOpt_optimization import Optimization

#pySNOPT
sys.path.append('../../../../pyACDT/pyACDT/Optimization/pyOpt/pySNOPT')
from pySNOPT import SNOPT

'''This script runs an example of fitting a cubic spline to a set of
3D wing coordinates that were generated using aerosurf. The script is
meant to show how to fit a constrained surface using a non-linear
optimizer'''

# Global Variables

Nctlu =10
Nctlv =6
ku=4
kv=4

# ==============================
# Read in the Data 
# ==============================

# Top surface of wing
print 'reading upper surface...'
data = read_array(file("upper_surface.inp"))

Nu = 26
Nv = 25

xu0 = zeros((Nu,Nv))
yu0 = zeros((Nu,Nv))
zu0 = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        xu0[i,j] = data[j*Nu + i ,0]
        yu0[i,j] = data[j*Nu + i ,1]
        zu0[i,j] = data[j*Nu + i ,2]
    #end for
#end for

# Bottom surface of wing
print 'reading bottom surface...'
data = read_array(file("lower_surface.inp"))

xl0 = zeros((Nu,Nv))
yl0 = zeros((Nu,Nv))
zl0 = zeros((Nu,Nv))

for j in xrange(Nv):
    for i in xrange(Nu):
        xl0[i,j] = data[j*Nu + i,0]
        yl0[i,j] = data[j*Nu + i,1]
        zl0[i,j] = data[j*Nu + i,2]
    #end for
#end for

#Flatten Data Once
xu0flat = xu0.flatten()
yu0flat = yu0.flatten()
zu0flat = zu0.flatten()
xl0flat = xl0.flatten()
yl0flat = yl0.flatten()
zl0flat = zl0.flatten()

# Check the Data

f = open('output.dat','w')
f.write ('VARIABLES = "X", "Y","Z"\n')
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(xu0[i,j],yu0[i,j],zu0[i,j]))
f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(xl0[i,j],yl0[i,j],zl0[i,j]))

# ==============================
# Find the 'u' parameters
# ==============================
x0=xu0[:,0];y0=yu0[:,0];z0=zu0[:,0] 
uu=zeros(Nu)
for i in xrange(Nu-1):
     uu[i+1] = uu[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
uu = uu/uu[-1] #Normalize uu

x0=xu0[:,0];y0=yu0[:,0];z0=zu0[:,0] 
ul=zeros(Nu)
for i in xrange(Nu-1):
     ul[i+1] = ul[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
ul = ul/ul[-1] #Normalize ul

# ==============================
# Find the 'v' parameters
# ==============================
x0=xu0[0,:];y0=yu0[0,:];z0=zu0[0,:] 
vu = zeros(Nv)
for i in xrange(Nv-1):
    vu[i+1] = vu[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
      
                         (z0[i+1]-z0[i])**2)
#end for
vu = vu/vu[-1] #Normalize vu

x0=xl0[0,:];y0=yl0[0,:];z0=zl0[0,:] 
vl = zeros(Nv)
for i in xrange(Nv-1):
    vl[i+1] = vl[i] + sqrt((x0[i+1]-x0[i])**2 + \
                               (y0[i+1]-y0[i])**2+ \
                               (z0[i+1]-z0[i])**2)
#end for
vl = vl/vl[-1] #Normalize vl

# ==============================
# Find the Initial Knot Vectors
# ==============================

# U knots

tu0= zeros(Nctlu + ku)
tu0[ku-1:Nctlu+1]  = 0.5*(1-cos(linspace(0,pi,Nctlu-ku+2)))
tu0[0:ku] = 0
tu0[Nctlu:Nctlu+ku] = 1.0#+0.1*(tu0[Nctlu]-tu0[Nctlu-1])


# V Knots

tv0= zeros(Nctlv + kv)
tv0[kv-1:Nctlv+1]=  linspace(0,1,Nctlv-kv+2)
#tv0[kv-1:Nctlv+1]  = 0.5*(1-cos(linspace(0,pi,Nctlv-kv+2)))
tv0[0:kv] = 0
tv0[Nctlv:Nctlv+kv] = 1.0#+0.1*(tv0[Nctlv]-tv0[Nctlv-1])

# Calculate the jacobian J, for fixed t and s
# Calculate the jacobian J, for fixed t and s
h = 1.0e-40j
Ju = zeros([Nu*Nv,Nctlu*Nctlv])
Jl = zeros([Nu*Nv,Nctlu*Nctlv])
ctl = zeros([Nctlu,Nctlv],'D')

[VU, UU] = meshgrid(vu,uu)
[VL, UL] = meshgrid(vu,uu)

for j in xrange(Nctlv):
    for i in xrange(Nctlu):
        ctl[i,j] += h
        valu = pyspline_cs.b2valv(UU.flatten(),VU.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        vall = pyspline_cs.b2valv(UL.flatten(),VL.flatten(),0,0,tu0,tv0,ku,kv,ctl)
        ctl[i,j] -= h    
        Ju[:,i*Nctlv + j] = imag(valu)/imag(h)
        Jl[:,i*Nctlv + j] = imag(vall)/imag(h)
    # end for 
# end for

global timeCounter 
timeCounter = 0




def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global timeCounter
    timeA = time.time()
 
    # Unpack the x-values
  
    ctlxu = reshape(hstack([0,x[                  0      :              Nctlv-2 ],0, \
                              x[            Nctlv-2      :Nctlu*Nctlv-Nctlv-2   ],0, \
                              x[Nctlu*Nctlv-Nctlv-2      :Nctlu*Nctlv-4         ],0]),[Nctlu,Nctlv])

    ctlyu = reshape(hstack([0,x[1*Nctlu*Nctlv-4          :1*Nctlu*Nctlv+Nctlv-6 ],0, \
                              x[1*Nctlu*Nctlv+Nctlv-6    :2*Nctlu*Nctlv-Nctlv-6 ],0,
                              x[2*Nctlu*Nctlv-Nctlv-6    :2*Nctlu*Nctlv-8       ],0]),[Nctlu,Nctlv])

    ctlzu = reshape(hstack([0,x[2*Nctlu*Nctlv-8          :2*Nctlu*Nctlv+Nctlv-10],0, \
                              x[2*Nctlu*Nctlv+Nctlv-10   :3*Nctlu*Nctlv-Nctlv-10],0,
                              x[3*Nctlu*Nctlv-Nctlv-10   :3*Nctlu*Nctlv-12      ],0]),[Nctlu,Nctlv])

    ctlxl = reshape(hstack([0,x[3*Nctlu*Nctlv-12         :3*Nctlu*Nctlv+Nctlv-14],0, \
                              x[3*Nctlu*Nctlv+Nctlv-14   :4*Nctlu*Nctlv-Nctlv-14],0,
                              x[4*Nctlu*Nctlv-Nctlv-14   :4*Nctlu*Nctlv-16      ],0]),[Nctlu,Nctlv])

    ctlyl = reshape(hstack([0,x[4*Nctlu*Nctlv-16         :4*Nctlu*Nctlv+Nctlv-18],0, \
                              x[4*Nctlu*Nctlv+Nctlv-18   :5*Nctlu*Nctlv-Nctlv-18],0,
                              x[5*Nctlu*Nctlv-Nctlv-18   :5*Nctlu*Nctlv-20      ],0]),[Nctlu,Nctlv])

    ctlzl = reshape(hstack([0,x[5*Nctlu*Nctlv-20         :5*Nctlu*Nctlv+Nctlv-22],0, \
                              x[5*Nctlu*Nctlv+Nctlv-22   :6*Nctlu*Nctlv-Nctlv-22],0,
                              x[6*Nctlu*Nctlv-Nctlv-22   :6*Nctlu*Nctlv-24      ],0]),[Nctlu,Nctlv])



   #  ctlxl = reshape(hstack([zeros(Nctlv),x[3*Nctlu*Nctlv-12         :4*Nctlu*Nctlv-12-2*Nctlv],zeros(Nctlv)]),[Nctlu,Nctlv])
#     ctlyl = reshape(hstack([zeros(Nctlv),x[4*Nctlu*Nctlv-12-2*Nctlv :5*Nctlu*Nctlv-12-4*Nctlv],zeros(Nctlv)]),[Nctlu,Nctlv])
#     ctlzl = reshape(hstack([zeros(Nctlv),x[5*Nctlu*Nctlv-12-4*Nctlv :6*Nctlu*Nctlv-12-6*Nctlv],zeros(Nctlv)]),[Nctlu,Nctlv])

 
    #Fix in the corners
    ctlxu[0 , 0] = xu0[0 , 0]
    ctlxu[0 ,-1] = xu0[0 ,-1]
    ctlxu[-1, 0] = xu0[-1, 0]
    ctlxu[-1,-1] = xu0[-1,-1]

    ctlyu[0 , 0] = yu0[0 , 0]
    ctlyu[0 ,-1] = yu0[0 ,-1]
    ctlyu[-1, 0] = yu0[-1, 0]
    ctlyu[-1,-1] = yu0[-1,-1]

    ctlzu[0 , 0] = zu0[0 , 0]
    ctlzu[0 ,-1] = zu0[0 ,-1]
    ctlzu[-1, 0] = zu0[-1, 0]
    ctlzu[-1,-1] = zu0[-1,-1]

    ctlxl[0 , 0] = xl0[0 , 0]
    ctlxl[0 ,-1] = xl0[0 ,-1]
    ctlxl[-1, 0] = xl0[-1, 0]
    ctlxl[-1,-1] = xl0[-1,-1]

    ctlyl[0 , 0] = yl0[0 , 0]
    ctlyl[0 ,-1] = yl0[0 ,-1]
    ctlyl[-1, 0] = yl0[-1, 0]
    ctlyl[-1,-1] = yl0[-1,-1]

    ctlzl[0 , 0] = zl0[0 , 0]
    ctlzl[0 ,-1] = zl0[0 ,-1]
    ctlzl[-1, 0] = zl0[-1, 0]
    ctlzl[-1,-1] = zl0[-1,-1]
  #   ctlxu = x[              0:  Nctlu*Nctlv- 4]
#     ctlyu = x[  Nctlu*Nctlv-4:2*Nctlu*Nctlv- 8]
#     ctlzu = x[2*Nctlu*Nctlv-8:3*Nctlu*Nctlv-12]

#     ctlxl = x[3*Nctlu*Nctlv        -12:4*Nctlu*Nctlv-2*Nctlv-12]
#     ctlyl = x[4*Nctlu*Nctlv-2*Nctlv-12:5*Nctlu*Nctlv-4*Nctlv-12]
#     ctlzl = x[5*Nctlu*Nctlv-4*Nctlv-12:6*Nctlu*Nctlv-6*Nctlv-12]
    fcon = hstack([ctlxu[0,:]-ctlxl[-1,:],ctlyu[0,:]-ctlyl[-1,:],ctlyu[0,:]-ctlyl[-1,:]])
    
    ctlxu = ctlxu.flatten()
    ctlyu = ctlyu.flatten()
    ctlzu = ctlzu.flatten()
    ctlxl = ctlxl.flatten()
    ctlyl = ctlyl.flatten()
    ctlzl = ctlzl.flatten()

    sum1 = sum((dot(Ju,ctlxu)-xu0flat)**2)
    sum2 = sum((dot(Ju,ctlyu)-yu0flat)**2)
    sum3 = sum((dot(Ju,ctlzu)-zu0flat)**2)
    sum4 = sum((dot(Jl,ctlxl)-xl0flat)**2)
    sum5 = sum((dot(Jl,ctlyl)-yl0flat)**2)
    sum6 = sum((dot(Jl,ctlzl)-zl0flat)**2)

    f= sum1+sum2+sum3+sum4+sum5+sum6
    fail = False
    timeCounter += time.time()-timeA
    return f,fcon,fail

def sens(x_orig,f_obj,f_con):
    global timeCounter
    timeA = time.time()
    ndv = len(x_orig)
    g_obj = zeros(ndv)
    ncon = 3*Nctlv
    g_con = zeros((ncon,ndv))
    for j in xrange(ndv):
        
        x = zeros([ndv],'complex')
        x[j] += h

    # Unpack the x-values
        ctlxu = reshape(hstack([0,x[                  0      :              Nctlv-2 ],0, \
                                    x[            Nctlv-2      :Nctlu*Nctlv-Nctlv-2   ],0, \
                                    x[Nctlu*Nctlv-Nctlv-2      :Nctlu*Nctlv-4         ],0]),[Nctlu,Nctlv])
        
        ctlyu = reshape(hstack([0,x[1*Nctlu*Nctlv-4          :1*Nctlu*Nctlv+Nctlv-6 ],0, \
                                    x[1*Nctlu*Nctlv+Nctlv-6    :2*Nctlu*Nctlv-Nctlv-6 ],0,
                                x[2*Nctlu*Nctlv-Nctlv-6    :2*Nctlu*Nctlv-8       ],0]),[Nctlu,Nctlv])
        
        ctlzu = reshape(hstack([0,x[2*Nctlu*Nctlv-8          :2*Nctlu*Nctlv+Nctlv-10],0, \
                                    x[2*Nctlu*Nctlv+Nctlv-10   :3*Nctlu*Nctlv-Nctlv-10],0,
                                x[3*Nctlu*Nctlv-Nctlv-10   :3*Nctlu*Nctlv-12      ],0]),[Nctlu,Nctlv])
        
        ctlxl = reshape(hstack([0,x[3*Nctlu*Nctlv-12         :3*Nctlu*Nctlv+Nctlv-14],0, \
                                    x[3*Nctlu*Nctlv+Nctlv-14   :4*Nctlu*Nctlv-Nctlv-14],0,
                                x[4*Nctlu*Nctlv-Nctlv-14   :4*Nctlu*Nctlv-16      ],0]),[Nctlu,Nctlv])
        
        ctlyl = reshape(hstack([0,x[4*Nctlu*Nctlv-16         :4*Nctlu*Nctlv+Nctlv-18],0, \
                                    x[4*Nctlu*Nctlv+Nctlv-18   :5*Nctlu*Nctlv-Nctlv-18],0,
                                x[5*Nctlu*Nctlv-Nctlv-18   :5*Nctlu*Nctlv-20      ],0]),[Nctlu,Nctlv])
        
        ctlzl = reshape(hstack([0,x[5*Nctlu*Nctlv-20         :5*Nctlu*Nctlv+Nctlv-22],0, \
                                    x[5*Nctlu*Nctlv+Nctlv-22   :6*Nctlu*Nctlv-Nctlv-22],0,
                                x[6*Nctlu*Nctlv-Nctlv-22   :6*Nctlu*Nctlv-24      ],0]),[Nctlu,Nctlv])
        

    #Fix in the corners
        ctlxu[0 , 0] = xu0[0 , 0]
        ctlxu[0 ,-1] = xu0[0 ,-1]
        ctlxu[-1, 0] = xu0[-1, 0]
        ctlxu[-1,-1] = xu0[-1,-1]
        
        ctlyu[0 , 0] = yu0[0 , 0]
        ctlyu[0 ,-1] = yu0[0 ,-1]
        ctlyu[-1, 0] = yu0[-1, 0]
        ctlyu[-1,-1] = yu0[-1,-1]

        ctlzu[0 , 0] = zu0[0 , 0]
        ctlzu[0 ,-1] = zu0[0 ,-1]
        ctlzu[-1, 0] = zu0[-1, 0]
        ctlzu[-1,-1] = zu0[-1,-1]
        
        ctlxl[0 , 0] = xl0[0 , 0]
        ctlxl[0 ,-1] = xl0[0 ,-1]
        ctlxl[-1, 0] = xl0[-1, 0]
        ctlxl[-1,-1] = xl0[-1,-1]
        
        ctlyl[0 , 0] = yl0[0 , 0]
        ctlyl[0 ,-1] = yl0[0 ,-1]
        ctlyl[-1, 0] = yl0[-1, 0]
        ctlyl[-1,-1] = yl0[-1,-1]
        
        ctlzl[0 , 0] = zl0[0 , 0]
        ctlzl[0 ,-1] = zl0[0 ,-1]
        ctlzl[-1, 0] = zl0[-1, 0]
        ctlzl[-1,-1] = zl0[-1,-1]
        
       

      #Now for g_con

        fcon = hstack([ctlxu[0,:]-ctlxl[-1,:],ctlyu[0,:]-ctlyl[-1,:],ctlyu[0,:]-ctlyl[-1,:]])
        g_con[:,j] = imag(fcon)/imag(h)
        x[j]-=h
    # end for

    #keep this shit out of the loop
       # Unpack the x-values
    x = x_orig
    ctlxu = reshape(hstack([0,x[                  0      :              Nctlv-2 ],0, \
                                x[            Nctlv-2      :Nctlu*Nctlv-Nctlv-2   ],0, \
                                x[Nctlu*Nctlv-Nctlv-2      :Nctlu*Nctlv-4         ],0]),[Nctlu,Nctlv])
        
    ctlyu = reshape(hstack([0,x[1*Nctlu*Nctlv-4          :1*Nctlu*Nctlv+Nctlv-6 ],0, \
                                x[1*Nctlu*Nctlv+Nctlv-6    :2*Nctlu*Nctlv-Nctlv-6 ],0,
                            x[2*Nctlu*Nctlv-Nctlv-6    :2*Nctlu*Nctlv-8       ],0]),[Nctlu,Nctlv])
    
    ctlzu = reshape(hstack([0,x[2*Nctlu*Nctlv-8          :2*Nctlu*Nctlv+Nctlv-10],0, \
                                x[2*Nctlu*Nctlv+Nctlv-10   :3*Nctlu*Nctlv-Nctlv-10],0,
                            x[3*Nctlu*Nctlv-Nctlv-10   :3*Nctlu*Nctlv-12      ],0]),[Nctlu,Nctlv])
    
    ctlxl = reshape(hstack([0,x[3*Nctlu*Nctlv-12         :3*Nctlu*Nctlv+Nctlv-14],0, \
                                x[3*Nctlu*Nctlv+Nctlv-14   :4*Nctlu*Nctlv-Nctlv-14],0,
                            x[4*Nctlu*Nctlv-Nctlv-14   :4*Nctlu*Nctlv-16      ],0]),[Nctlu,Nctlv])
    
    ctlyl = reshape(hstack([0,x[4*Nctlu*Nctlv-16         :4*Nctlu*Nctlv+Nctlv-18],0, \
                                x[4*Nctlu*Nctlv+Nctlv-18   :5*Nctlu*Nctlv-Nctlv-18],0,
                            x[5*Nctlu*Nctlv-Nctlv-18   :5*Nctlu*Nctlv-20      ],0]),[Nctlu,Nctlv])
    
    ctlzl = reshape(hstack([0,x[5*Nctlu*Nctlv-20         :5*Nctlu*Nctlv+Nctlv-22],0, \
                                x[5*Nctlu*Nctlv+Nctlv-22   :6*Nctlu*Nctlv-Nctlv-22],0,
                            x[6*Nctlu*Nctlv-Nctlv-22   :6*Nctlu*Nctlv-24      ],0]),[Nctlu,Nctlv])
    

    #Fix in the corners
    ctlxu[0 , 0] = xu0[0 , 0]
    ctlxu[0 ,-1] = xu0[0 ,-1]
    ctlxu[-1, 0] = xu0[-1, 0]
    ctlxu[-1,-1] = xu0[-1,-1]
    
    ctlyu[0 , 0] = yu0[0 , 0]
    ctlyu[0 ,-1] = yu0[0 ,-1]
    ctlyu[-1, 0] = yu0[-1, 0]
    ctlyu[-1,-1] = yu0[-1,-1]

    ctlzu[0 , 0] = zu0[0 , 0]
    ctlzu[0 ,-1] = zu0[0 ,-1]
    ctlzu[-1, 0] = zu0[-1, 0]
    ctlzu[-1,-1] = zu0[-1,-1]
    
    ctlxl[0 , 0] = xl0[0 , 0]
    ctlxl[0 ,-1] = xl0[0 ,-1]
    ctlxl[-1, 0] = xl0[-1, 0]
    ctlxl[-1,-1] = xl0[-1,-1]
        
    ctlyl[0 , 0] = yl0[0 , 0]
    ctlyl[0 ,-1] = yl0[0 ,-1]
    ctlyl[-1, 0] = yl0[-1, 0]
    ctlyl[-1,-1] = yl0[-1,-1]
    
    ctlzl[0 , 0] = zl0[0 , 0]
    ctlzl[0 ,-1] = zl0[0 ,-1]
    ctlzl[-1, 0] = zl0[-1, 0]
    ctlzl[-1,-1] = zl0[-1,-1]
        
    ctlxu = real(ctlxu.flatten())
    ctlyu = real(ctlyu.flatten())
    ctlzu = real(ctlzu.flatten())
    ctlxl = real(ctlxl.flatten())
    ctlyl = real(ctlyl.flatten())
    ctlzl = real(ctlzl.flatten())


    temp = 2*dot(dot(Ju,ctlxu)-xu0flat,Ju)
    g_obj[                0: Nctlu*Nctlv- 4] = hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])

    temp = 2*dot(dot(Ju,ctlyu)-yu0flat,Ju)
    g_obj[  Nctlu*Nctlv-4 :2*Nctlu*Nctlv- 8] = hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])

    temp = 2*dot(dot(Ju,ctlzu)-zu0flat,Ju)
    g_obj[2*Nctlu*Nctlv-8 :3*Nctlu*Nctlv-12] = hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])


    temp = 2*dot(dot(Jl,ctlxl)-xl0flat,Jl)
    g_obj[3*Nctlu*Nctlv-12:4*Nctlu*Nctlv-16] =  hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])

    temp = 2*dot(dot(Jl,ctlyl)-yl0flat,Jl)
    g_obj[4*Nctlu*Nctlv-16:5*Nctlu*Nctlv-20] =  hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])

    temp = 2*dot(dot(Jl,ctlzl)-zl0flat,Jl)
    g_obj[5*Nctlu*Nctlv-20:6*Nctlu*Nctlv-24] =  hstack([temp[1:Nctlv-1],temp[Nctlv:Nctlv*Nctlu-Nctlv],temp[Nctlu*Nctlv-Nctlv+1:Nctlu*Nctlv-1]])

#     g_obj[                0: Nctlu*Nctlv- 4] = 2*dot(dot(Ju,ctlxu)-xu0flat,Ju)
#     g_obj[  Nctlu*Nctlv-4 :2*Nctlu*Nctlv- 8] = 2*dot(dot(Ju,ctlyu)-yu0flat,Ju)
#     g_obj[2*Nctlu*Nctlv-8 :3*Nctlu*Nctlv-12] = 2*dot(dot(Ju,ctlzu)-zu0flat,Ju)
#     g_obj[3*Nctlu*Nctlv-12:4*Nctlu*Nctlv-12-2*Nctlv] = 2*dot(dot(Jl,ctlxl)-xl0flat,Jl)
#     g_obj[4*Nctlu*Nctlv-12-2*Nctlv:5*Nctlu*Nctlv-12-4*Nctlv] = 2*dot(dot(Jl,ctlyl)-yl0flat,Jl)
#     g_obj[5*Nctlu*Nctlv-12-4*Nctlv:6*Nctlu*Nctlv-12-6*Nctlv] =2*dot(dot(Jl,ctlzl)-zl0flat,Jl)
    
    fail = False
    #rint 'constr time:',time.time()-timeA
    timeCounter += time.time()-timeA
    return g_obj,g_con,fail

# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

# ================================================
# Find a good guess of the initial control points 
# ================================================

ctlxu = zeros(Nctlu*Nctlv-4)
ctlyu = zeros(Nctlu*Nctlv-4)
ctlzu = zeros(Nctlu*Nctlv-4)
ctlxl = zeros(Nctlu*Nctlv-4)
ctlyl = zeros(Nctlu*Nctlv-4)
ctlzl = zeros(Nctlu*Nctlv-4)

#Create the interpolation

Ixu = RectBivariateSpline(uu,vu,xu0,kx=1,ky=1)
Iyu = RectBivariateSpline(uu,vu,yu0,kx=1,ky=1)
Izu = RectBivariateSpline(uu,vu,zu0,kx=1,ky=1)
Ixl = RectBivariateSpline(ul,vl,xl0,kx=1,ky=1)
Iyl = RectBivariateSpline(ul,vl,yl0,kx=1,ky=1)
Izl = RectBivariateSpline(ul,vl,zl0,kx=1,ky=1)

u_interp = 0.5*(1-cos(linspace(0,pi,Nctlu)))
v_interp = linspace(0,1,Nctlv)

counter = 0
for i in xrange(Nctlu):
    for j in xrange(Nctlv):
        if not(    (i== 0       and j ==0       ) or\
                   (i== 0       and j == Nctlv-1) or\
                   (i== Nctlu-1 and j == 0      ) or\
                   (i== Nctlu-1 and j == Nctlv-1)):
            ctlxu[counter] = Ixu(u_interp[i],v_interp[j])
            ctlyu[counter] = Iyu(u_interp[i],v_interp[j])
            ctlzu[counter] = Izu(u_interp[i],v_interp[j])
            counter +=1
        else:
            print Ixu(u_interp[i],v_interp[j]),Iyu(u_interp[i],v_interp[j]),Izu(u_interp[i],v_interp[j])

counter = 0
for i in xrange(Nctlu):
    for j in xrange(Nctlv):
        if not(    (i== 0       and j ==0       ) or\
                   (i== 0       and j == Nctlv-1) or\
                   (i== Nctlu-1 and j == 0      ) or\
                   (i== Nctlu-1 and j == Nctlv-1)):
            ctlxl[counter] = Ixl(u_interp[i],v_interp[j])
            ctlyl[counter] = Iyl(u_interp[i],v_interp[j])
            ctlzl[counter] = Izl(u_interp[i],v_interp[j])
            counter += 1
print 'init'
print ctlxu
opt_prob.addVarGroup('CTLxu',(Nctlu*Nctlv-4),'c',value=ctlxu,lower=-100,upper=100)
opt_prob.addVarGroup('CTLyu',(Nctlu*Nctlv-4),'c',value=ctlyu,lower=-100,upper=100)
opt_prob.addVarGroup('CTLzu',(Nctlu*Nctlv-4),'c',value=ctlzu,lower=-100,upper=100)
opt_prob.addVarGroup('CTLxl',(Nctlu*Nctlv-4),'c',value=ctlxl,lower=-100,upper=100)
opt_prob.addVarGroup('CTLyl',(Nctlu*Nctlv-4),'c',value=ctlyl,lower=-100,upper=100)
opt_prob.addVarGroup('CTLzl',(Nctlu*Nctlv-4),'c',value=ctlzl,lower=-100,upper=100)

# ===================
#  Constraints
# ===================
opt_prob.addConGroup('match_consts',3*Nctlv,type = 'e',value=0.0)

#opt_prob.addVarGroup('corners',24,'e',value=coners)

# Knots are not implemented yet



# ===================
#  Objective
# ===================
opt_prob.addObj('RMSerror')

# ===================
#  Select Optimizer
# ===================

opt = SNOPT()

# ===================
#  SNOPT Options
# ===================
#opt.setOption('Derivative level',0)
opt.setOption('Verify level',3)
opt.setOption('Major iterations limit',150)
#opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-5)
opt.setOption('Major feasibility tolerance',1e-5)
opt.setOption('Minor feasibility tolerance',1e-5)

# ===================
#  Run Optimization
# ===================

result = opt(opt_prob,sens)

# # make  an x and run obj
# x = array(hstack([ctlxu,ctlyu,ctlzu,ctlxl,ctlyl,ctlzl]),'D')

# fref = objcon(x)
# gobj = sens(x,fref[0],fref[1])
# g_sens = gobj[0][0]

# #x[0] = x[0] + 3.88e-6
# x[0]= x[0]+ 1.0e-40j
# f = objcon(x)
# print 'f:',f
# print 'finite diff:',imag(f[0])/1.0e-40
# print 'sens:',g_sens

# sys.exit(0)

print opt_prob._solutions[0]
print '#--------------------------------'
print '# RMS Error: ',sqrt(result[0][0]/(2*Nu*Nv))
print '#--------------------------------'


# ===================
#  Print Solution  
# ===================
x = result[1][:]
# Unpack the x-values



f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlxu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlyu),
                               pyspline.b2val(uu[i],vu[j],0,0,tu0,tv0,ku,kv,ctlzu)))



f.write('Zone I=%d J = %d\n'%(Nu,Nv))
f.write('DATAPACKING=POINT\n')
for j in xrange(Nv):
    for i in xrange(Nu):
        f.write('%f %f %f \n'%(pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlxl),
                               pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlyl),
                               pyspline.b2val(ul[i],vl[j],0,0,tu0,tv0,ku,kv,ctlzl)))



print 'Eval Time:',timeCounter
