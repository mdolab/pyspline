# This is a quick and dirty implementation of a linear spline function
from numpy import array,zeros,vstack,mat,transpose
import pdb,sys,time
import pyspline
from scipy.linalg import lstsq,solve
from mdo_import_helper import *
exec(import_modules('pySpline','pyOpt_optimization','pySNOPT','geo_utils'))
from scipy.io import savemat

def interval(knots,t):
    '''Get the interval for the value t in knots'''
    right = knots.searchsorted(t,side='right')
    left  = knots.searchsorted(t,side='left')
    return left,right
global call_count
call_count = 0

def basis(t,i,k,knots):

    global call_count
    call_count += 1
    if k == 1:
        if t>=knots[i] and t < knots[i+1]:
            return 1
        else:
            return 0
         
        # end if
    else:
        
        b1 = basis(t,i  ,k-1,knots)
        b2 = basis(t,i+1,k-1,knots)

        if b1 == 0:
            val1 = 0
        else:
            val1 = (t - knots[i])/(knots[i+k-1] - knots[i])*b1
        # end if

        if b2 == 0:
            val2 = 0
        else:
            val2 = (knots[i+k]-t)/(knots[i+k] - knots[i+1])*b2
        # end if

        return val1+val2

    # end if

# def getBlending(t,knots,k):
#     '''Get the full blending matrix'''
#     n = len(knots)-k
#     B = zeros(n)
#     for i in xrange(n):
#         B[i] = basis(t,i,k,knots)
#     # end for

#     return B

def getBasisFortran(t,i,k,knots,deriv):
    n = len(knots)-k
    coef = zeros(n)
    coef[i] = 1
    val = pyspline.bvalu(knots,coef,k,deriv,t)
    return val

def getBlendingFortran(t,knots,k,deriv):
    n = len(knots)-k
    B = zeros(n)
    coef = zeros(n)
    for i in xrange(n):
        coef[i] += 1
        B[i] = pyspline.bvalu(knots,coef,k,deriv,t)
        coef[i] -= 1
    # end for
    
    return B 

# def getValue(t,knots,coef,k,deriv):
#     B = getBlendingFortran(t,knots,k,deriv)
#     val = dot(B,coef)
#     return val

def greville(knots,k):
    ncoef = len(knots)-k
    gpts = zeros(ncoef)
    for i in xrange(ncoef):
        for n in xrange(k-1): #degree n loop
            gpts[i] += knots[i+n+1]
        # end for
        gpts[i] /= (k-1)
    # end for
    return gpts

n = 55
k = 4
N = 80
#start_foil = 'naca0012.dat'
#xu,yu,xl,yl = read_af(start_foil,'xfoil',N)
#upper = vstack([xu,yu]).T.astype('D')
#lower = vstack([xl[::-1],yl[::-1]]).T.astype('D')

# Get some new data
theta = linspace(-pi/6,pi/2,100)
x = cos(theta)
y = sin(theta)
upper = vstack([x,y]).T.astype('D')
top = pySpline.linear_spline(task='lms',k=k,Nctl=n,X=upper,complex=True)

# Open the matfile and save basline curve
savemat('out.mat', {'top':top.getValueV(linspace(0,1,100))})
timeA = time.time()
knots = top.t
coef = top.coef

# Create the k matrix
alpha = 10
beta  = 1
gpts = greville(knots,k)
K0 = pyspline.curve_stiffness(knots,k,alpha,beta,gpts)

# Do an auto-magic peturbation
delta_left = [cos(-pi/6),sin(-pi/6)+.15]

# Get the tangent vector
d0 = top.getDerivative(0)
s0 = dot(delta_left-top(0),d0)/dot(d0,d0)

Constr = zeros((n,4))
for i in xrange(n):
    Constr[i,0] = gpts[i]
    Constr[i,1:3] = top(  (1-gpts[i])*s0+gpts[i]    )-top(gpts[i])
    Constr[i,3] = 1e3
# end for
knew1,f1 = pyspline.apply_constr(knots,k,K0,Constr)

p_delta = s0
gamma = []
gpts = greville(knots,k)
constr= []
for i in xrange(len(gpts)):
    constr.append([gpts[i],top(  (1-gpts[i])*p_delta+gpts[i]    )-top(gpts[i])])
    gamma.append(1e3)
# end for

Constr2 = zeros([2,4])
Constr2[0,0] = 0
Constr2[0,1:3] = delta_left-top(s0)
Constr2[0,3]   = 1e3
Constr2[1,0] = 1
Constr2[1,1:3] = [0,0]
Constr2[1,3]   = 1e3

knew2,f2 = pyspline.apply_constr(knots,k,K0,Constr2)

temp = solve(knew1,f1,sym_pos=1)
temp2 = solve(knew2,f2,sym_pos=1)

coef_x = coef[:,0] +temp[:,0] + temp2[:,0]
coef_y = coef[:,1] +temp[:,1] + temp2[:,1]

top.coef[:,0] = coef_x
top.coef[:,1] = coef_y
timeB = time.time()
savemat('out2.mat', {'d0':d0,'top_mod':top.getValueV(linspace(0,1,100)),'delta_left':delta_left})
print 'Time is:',timeB-timeA



# f = zeros((n,2))
# timeB = time.time()
# K = K0.copy()
# for i in xrange(n):
#     # Get the basis vectors for this constraint
    
#     B = getBlendingFortran(constr[i][0],knots,k,0)
    
#     C = transpose(mat(B))
#     D = mat(B)
#     M = dot(C,D)
#     K += array(M)*gamma[i]
#     f[:,0] += gamma[i]*constr[i][1][0]*(B)
#     f[:,1] += gamma[i]*constr[i][1][1]*(B)

#print 'f:',f
#print 'f1:',f1
#print 'k',K
#print 'knew:',knew1
