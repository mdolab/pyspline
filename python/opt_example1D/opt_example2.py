# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, string, pdb, copy

# =============================================================================
# External Python modules
# =============================================================================
from numpy import linspace, cos, sin, arccos, pi, hstack, zeros, ones, sqrt,\
    imag, interp, array, real, dot
from matplotlib.pylab import plot,show

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
airfoil coordinates, specificially the upper surface of a naca 0012
airfoil with a constraint on a vertical leading edge'''


# Global Variables

global ctlxu,ctlyu,ctlxl,ctlyl
global su,sl,Ju,Jl

Nctl = 13
k = 4

# Top surface of a naca0012 airfoil
N = 91
x0 = 0.5*(1-cos(linspace(0,pi,N)))
y0 = 0.12/0.2*(0.2969*x0**0.5-0.1260*x0-0.3516*x0**2 + 0.2843*x0**3-0.1015*x0**4)
x0u = x0
x0l = x0
y0u = y0
y0l = -y0
# Surface of WT Airfoil
# N=32
# x0u = array([0.0000000,0.0001800,0.0025500,0.0095400,0.0209000,0.0365000,0.0564000,0.0803000,0.1080100,0.1393400,0.1739500,0.2114600,0.2514900,0.2936100,0.3373600,0.3822800,0.4282000,0.4752600,0.5232400,0.6198000,0.6672400,0.7133300,0.7574900,0.7991500,0.8377800,0.9039100,0.9307200,0.9535500,0.9725100,0.9871900,0.9966800,1.0000000])
# y0u = array([0.0000000,0.0015900,0.0074800,0.0164000,0.0260000,0.0358000,0.0456000,0.0552000,0.0643000,0.0729000,0.0807000,0.0876000,0.0934000,0.0981000,0.1013300,0.1029400,0.1024900,0.1000500,0.0961000,0.0849000,0.0782000,0.0710000,0.0634000,0.0557000,0.0480000,0.0326000,0.0248000,0.0170000,0.0098200,0.0043100,0.0010300,0.0000000])
# x0l = array([1.0000000,0.9965500,0.9860400,0.9681300,0.9425700,0.9094500,0.8695300,0.8240900,0.7743400,0.7214200,0.6664400,0.6105500,0.5548400,0.5003200,0.4478500,0.3977900,0.3502700,0.3049700,0.2615300,0.2198700,0.1796500,0.1413300,0.1063700,0.0758000,0.0500000,0.0292000,0.0137000,0.0036700,0.0021600,0.0009300,0.0002100,0])
# y0l = array([0.0000000,0.0008400,0.0032400,0.0068900,0.0111000,0.0150000,0.0174000,0.0178000,0.0162000,0.0125000,0.0069900,-0.0001500,-0.0084100,-0.0172000,-0.0256000,-0.0328000,-0.0379000,-0.0405000,-0.0406000,-0.0386000,-0.0349000,-0.0309000,-0.0273000,-0.0236000,-0.0196000,-0.0152000,-0.0104000,-0.0052500,-0.0040300,-0.0027400,-0.0014600,0])
# x0l = x0l[::-1]
# y0l = y0l[::-1]

# S-Parameter Calculation
s0u = zeros(N)
s0l = zeros(N)
for i in xrange(N-1):
    s0u[i+1] = s0u[i] + sqrt((x0u[i+1]-x0u[i])**2 + (y0u[i+1]-y0u[i])**2)
    s0l[i+1] = s0l[i] + sqrt((x0l[i+1]-x0l[i])**2 + (y0l[i+1]-y0l[i])**2)
muu = s0u[-1]
mul = s0l[-1]
s0u = s0u/s0u[-1] #Normalize s
s0l = s0l/s0l[-1] #Normalize s

# Initial Knot Vector Calculation
t0u= zeros(Nctl + k)
t0u[0:k] = s0u[0]
t0u[Nctl:Nctl+k] = s0u[-1]
t0u[k-1:Nctl+1]= 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))

t0l= zeros(Nctl + k)
t0l[0:k] = s0l[0]
t0l[Nctl:Nctl+k] = s0l[-1]
t0l[k-1:Nctl+1]= 0.5*(1-cos(linspace(0,pi,Nctl-k+2)))

def objcon(x):
    '''Get the rms error for the given set of design variables'''
    global ctlxu,ctlyu,ctlxl,ctlyl
    global Ju,Jl
    ctlxu = zeros(Nctl)
    ctlyu = zeros(Nctl)
    ctlxl = zeros(Nctl)
    ctlyl = zeros(Nctl)    

    ctlxu[1:Nctl-1] = x[       0:  Nctl-2] # These are the control points
    ctlyu[1:Nctl-1] = x[  Nctl-2:2*Nctl-4]
    ctlxl[1:Nctl-1] = x[2*Nctl-4:3*Nctl-6] # These are the control points
    ctlyl[1:Nctl-1] = x[3*Nctl-6:4*Nctl-8]

    ctlxu[0]  = x0u[0]
    ctlxu[-1] = x0u[-1]
    ctlyu[0]  = y0u[0]
    ctlyu[-1] = y0u[-1]

    ctlxl[0]  = x0l[0]
    ctlxl[-1] = x0l[-1]
    ctlyl[0]  = y0l[0]
    ctlyl[-1] = y0l[-1]


    f = sum(   (dot(Ju,ctlxu)-x0u)**2 + (dot(Ju,ctlyu)-y0u)**2 + (dot(Jl,ctlxl)-x0l)**2 + (dot(Jl,ctlyl)-y0l)**2)
    fail = False
    x1 = ctlxl[1]
    y1 = ctlyl[1]
    x3 = ctlxu[1]
    y3 = ctlyu[1]

    fcon = array([-x1*y3+x3*y1])
    return f,fcon,fail

def sens(x,f_obj,f_con):
    global Ju,Jl

    ndv = len(x)
    ncon = 1
    g_obj = zeros(ndv) #objective derivative vector
    g_con = zeros([ncon, ndv]) #constraint jacobian
    
    ctlxu = zeros(Nctl)
    ctlyu = zeros(Nctl)
    ctlxl = zeros(Nctl)
    ctlyl = zeros(Nctl)    

    ctlxu[1:Nctl-1] = x[       0:  Nctl-2] # These are the control points
    ctlyu[1:Nctl-1] = x[  Nctl-2:2*Nctl-4]
    ctlxl[1:Nctl-1] = x[2*Nctl-4:3*Nctl-6] # These are the control points
    ctlyl[1:Nctl-1] = x[3*Nctl-6:4*Nctl-8]

    ctlxu[0]  = x0u[0]
    ctlxu[-1] = x0u[-1]
    ctlyu[0]  = y0u[0]
    ctlyu[-1] = y0u[-1]

    ctlxl[0]  = x0l[0]
    ctlxl[-1] = x0l[-1]
    ctlyl[0]  = y0l[0]
    ctlyl[-1] = y0l[-1]
 
    g_obj[       0:  Nctl-2] = 2*dot(dot(Ju,ctlxu)-x0u,Ju)[1:Nctl-1]
    g_obj[  Nctl-2:2*Nctl-4] = 2*dot(dot(Ju,ctlyu)-y0u,Ju)[1:Nctl-1]
    g_obj[2*Nctl-4:3*Nctl-6] = 2*dot(dot(Jl,ctlxl)-x0l,Jl)[1:Nctl-1]
    g_obj[3*Nctl-6:4*Nctl-8] = 2*dot(dot(Jl,ctlyl)-y0l,Jl)[1:Nctl-1]

    x1 = ctlxl[1]
    y1 = ctlyl[1]
    x3 = ctlxu[1]
    y3 = ctlyu[1]

    g_con[0,0] = y1   # position of 'x3'
    g_con[0,Nctl-2] = -x1 # position of 'y3'
    g_con[0,2*Nctl-4] = -y3 #position of 'x1'
    g_con[0,3*Nctl-5] = x3

    fcon = array([-x1*y3+x3*y1])

    fail = False
    return g_obj,g_con,fail

# =============================================================================
#  Run Optimization Problem
# =============================================================================

# ===================
#  Variables
# ===================

opt_prob = Optimization('Cubic Spline Optimization Problem',objcon)

#Make sure the end points are what we want:

opt_prob.addVarGroup('CTLxu',Nctl-2,'c',value=0,lower=-1,upper=1)
opt_prob.addVarGroup('CTLyu',Nctl-2,'c',value=0,lower=-1,upper=1)
opt_prob.addVarGroup('CTLxl',Nctl-2,'c',value=0,lower=-1,upper=1)
opt_prob.addVarGroup('CTLyl',Nctl-2,'c',value=0,lower=-1,upper=1)

# ===================
#  Constraints
# ===================
opt_prob.addCon('LEconstraint',type = 'i',lower=0,upper=0)

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
opt.setOption('Derivative level',3)
opt.setOption('Nonderivative linesearch')
opt.setOption('Major optimality tolerance', 1e-6)
opt.setOption('Major feasibility tolerance',1e-6)
opt.setOption('Minor feasibility tolerance',1e-6)

# ===================================================
#  Run Optimization with Hoschek Parameter Correction
# ===================================================
Niter = 50
su = s0u
sl = s0l
#print 'su before:',su
for iter in xrange(Niter):


    # Calculate the jacobian J, for fixed t and s
    h = 1.0e-40j
    Ju = zeros([N,Nctl])
    Jl = zeros([N,Nctl])
    ctl = zeros([Nctl],'D')
    for i in xrange(Nctl):
        ctl[i] += h
        valu = pyspline_cs.bvaluv(t0u,ctl,k,0,su)
        vall = pyspline_cs.bvaluv(t0l,ctl,k,0,sl)
        ctl[i] -= h    
        Ju[:,i] = imag(valu)/imag(h)
        Jl[:,i] = imag(vall)/imag(h)
    # end for 

    # Determine approximation curve using SNOPT
    result = opt(opt_prob,sens)
 
    RMS = sqrt(result[0][0]/(2*N))
    print 'RMS:',RMS

    ctlxu = hstack([x0u[0],result[1][       0:  Nctl-2],x0u[-1]])
    ctlyu = hstack([y0u[0],result[1][  Nctl-2:2*Nctl-4],y0u[-1]])
    ctlxl = hstack([x0l[0],result[1][2*Nctl-4:3*Nctl-6],x0l[-1]])
    ctlyl = hstack([y0l[0],result[1][3*Nctl-6:4*Nctl-8],y0l[-1]])



  #   Du = zeros([N,2])
#     Dl = zeros([N,2])
#     Du[:,0] = dot(Ju,ctlxu)-x0u
#     Du[:,1] = dot(Ju,ctlyu)-y0u
#     Dl[:,0] = dot(Jl,ctlxl)-x0l
#     Dl[:,1] = dot(Jl,ctlyl)-y0l

#     ydotu = zeros([N,2])
#     ydotl = zeros([N,2])
#     ydotu[:,0] = pyspline.bvaluv(t0u,ctlxu,k,1,su)
#     ydotu[:,1] = pyspline.bvaluv(t0u,ctlyu,k,1,su)

#     ydotl[:,0] = pyspline.bvaluv(t0l,ctlxl,k,1,sl)
#     ydotl[:,1] = pyspline.bvaluv(t0l,ctlyl,k,1,sl)



    for i in xrange(N):

        #Upper Surface

        for j in xrange(25):
#            print 'muu:',muu
            if su[i] > 1:
                su[i] = 1
            if su[i] < 0:
                su[i] =0
            
            
            ydotx = pyspline.bvalu(t0u,ctlxu,k,1,su[i])
            ydoty = pyspline.bvalu(t0u,ctlyu,k,1,su[i])
            x=pyspline.bvalu(t0u,ctlxu,k,0,su[i])
            y=pyspline.bvalu(t0u,ctlyu,k,0,su[i])
            
            ydotnorm = sqrt(ydotx**2+ydoty**2)
            delta_ci = (((x0u[i]-x)*ydotx + (y0u[i]-y)*ydoty) /ydotnorm)/(2**j)
            s_tilde = su[i] + delta_ci*(su[-1]-su[0])/muu
           #  print 'j:',j
#             print 'ydotx:',ydotx
#             print 'ydoty:',ydoty
#             print 'ydotnorm:',ydotnorm
#             print 'delta_ci:',delta_ci
#             print 's:',su[i]
#             print 's_tilde:',s_tilde

            if s_tilde > 1:
                s_tilde = 1
            if s_tilde < 0:
                s_tilde = 0

            norm_Di = sqrt((x-x0u[i])**2 + (y-y0u[i])**2)

            x_tilde = pyspline.bvalu(t0u,ctlxu,k,0,s_tilde)
            y_tilde = pyspline.bvalu(t0u,ctlyu,k,0,s_tilde)
            norm_Di_tilde = sqrt((x_tilde-x0u[i])**2 + (y_tilde-y0u[i])**2)

#             print 'norm_Di:',norm_Di
#             print 'norm_Di_tilde:',norm_Di_tilde

            if norm_Di >= norm_Di_tilde:
                su[i] = s_tilde
                break
           # if j == 9:
                #print 'no better top'

        #Lower Surface
    for i in xrange(N):
        for j in xrange(10):

            if sl[i] > 1:
                sl[i] = 1
            if sl[i] < 0:
                sl[i] =0
            ydotx = pyspline.bvalu(t0l,ctlxl,k,1,sl[i])
            ydoty = pyspline.bvalu(t0l,ctlyl,k,1,sl[i])
            x=pyspline.bvalu(t0l,ctlxl,k,0,sl[i])
            y=pyspline.bvalu(t0l,ctlyl,k,0,sl[i])
            
            ydotnorm = sqrt(ydotx**2+ydoty**2)
            delta_ci = (((x-x0l[i])*ydotx + (y-y0l[i])*ydoty) /ydotnorm)/(2**j)
            s_tilde = sl[i] + delta_ci*(sl[-1]-sl[0])/mul

            if s_tilde > 1:
                s_tilde = 1
            if s_tilde < 0:
                s_tilde = 0
            

            norm_Di = sqrt((x-x0l[i])**2 + (y-y0l[i])**2)

            x_tilde = pyspline.bvalu(t0l,ctlxl,k,0,s_tilde)
            y_tilde = pyspline.bvalu(t0l,ctlyl,k,0,s_tilde)
            norm_Di_tilde = sqrt((x_tilde-x0l[i])**2 + (y_tilde-y0l[i])**2)

            if norm_Di >= norm_Di_tilde:
                #print 'delta s:',sl[i]-s_tilde
                sl[i] = s_tilde
                break
          #  if j==9:
                
             #   print 'no better bottom:'
            
        
 



#print 'suafter:',su
# ===================
#  Print Solution  
# ===================

#print help(opt_prob._solutions[0])
# print '#--------------------------------'
# print '# RMS Error: ',sqrt(x_ref/(2*N))
# print '#--------------------------------'

# print 'ctlxu:',ctlxu
# print 'ctlyu:',ctlyu
# print 'ctlxl:',ctlxl
# print 'ctlyl:',ctlyl

#s = 0.5*(1-cos(linspace(0,pi,N00)))
XU = pyspline.bvaluv(t0u,ctlxu,k,0,su)
YU = pyspline.bvaluv(t0u,ctlyu,k,0,su)
XL = pyspline.bvaluv(t0l,ctlxl,k,0,sl)
YL = pyspline.bvaluv(t0l,ctlyl,k,0,sl)

#plot(ctlxu,ctlyu,'ro')
#plot(ctlxl,ctlyl,'go')
plot(x0u,y0u,'ko-')
plot(x0l,y0l,'ko-')
plot(XU,YU,'ro')
plot(XL,YL,'go')
show()