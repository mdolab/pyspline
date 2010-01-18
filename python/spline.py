# This is a quick and dirty implementation of a linear spline function
from numpy import *
import pdb,sys,time
import pyspline
import pySpline
from petsc4py import PETSc
from scipy.linalg import lstsq,solve
from mdo_import_helper import *
exec(import_modules('geo_utils'))
from scipy.io import savemat
from matplotlib.pylab import plot,show,axis
import sys
sys.path.append('../../../pyMIL/pyXLIGHT/python')
import copy
import pyXLIGHT



airfoil = pyXLIGHT.xfoilAnalysis('naca0012.dat')
airfoil.re = 100000
airfoil.mach = 0.0
airfoil.iter = 100

nacax,nacay = loadtxt('naca0012.dat',unpack=True)
#nacax,nacay = loadtxt('sc20714.dat',unpack=True)
#nacax,nacay = loadtxt('sg6042.dat',unpack=True)


curve = pySpline.linear_spline(task='lms',k=4,Nctl=10,x=nacax,y=nacay)
x = linspace(-.25,1.25,250)
y = linspace(-.2,.2,110)


vals = curve(curve.s)
plot(vals[:,0],vals[:,1])
axis('equal')
timeA = time.time()
for i in xrange(len(x)):
    for j in xrange(len(y)):
        s = curve.pyspline.project_point_curve([x[i],y[j]],curve.t,curve.k,curve.coef,4,1e-6,1e-6)
        val = curve(s)
        plot([val[0],x[i]],[val[1],y[j]],'k-')
    # end for
# end for
show()
timeB = time.time()
print 'time:',timeB-timeA        
print 'avg time/call:',(timeB-timeA)/(250*110.0)


sys.exit(0)


# Get some new data
#n=4
#k=4
#theta = linspace(0,pi/2,10)
#x = cos(theta)
#y = sin(theta)
#upper = vstack([x,y]).T.astype('D')
#top = pySpline.linear_spline(task='lms',k=k,Nctl=n,X=upper)

#knots = top.t
#coef = zeros((n,3))
#coef[:,0:2] = top.coef
#coef[:,2] = 1.0
#print 'knots:',knots
#print 'coef:',coef

# #Test the compute curve
# k = 3
# Nctl = 3
# n = 25
# theta = linspace(0,pi/2,n)
# x = cos(theta)
# y = sin(theta)
# X = vstack([x,y]).T.astype('d')
# s0 = linspace(0,1,n)
# knots = pyspline.knots(s0,Nctl,k)
# s = copy.deepcopy(s0)
# s,coef = pyspline.compute_curve(s,X,knots,k,Nctl)
# print coef
#sys.exit(0)


# A split section

# Next setup the pySpline Object
# Load the airfoil using geo_utils read_af
Nctl = 13
xu,yu,xl,yl = read_af('naca0012.dat','xfoil',80)
upper = vstack([xu[::-1],yu[::-1]]).T.astype('d')
lower = vstack([xl[::-1],yl[::-1]]).T.astype('d')
top = pySpline.linear_spline(task='lms',k=4,Nctl=Nctl,X=upper)
bottom = pySpline.linear_spline(task='lms',k=4,Nctl=Nctl,X=lower)

print xu,xl
for i in xrange(len(xl)):
    print i,xu[i],xl[::-1][i],yu[i],yl[::-1][i]
print top.t
print bottom.t

f=open('naca0012_mod.dat','w')
for i in xrange(len(top.s)):
    val = top(top.s[i])
    f.write('%f %f\n'%(val[0],val[1]))
for i in xrange(len(bottom.s)):
    val = bottom(bottom.s[i])
    f.write('%f %f\n'%(val[0],val[1]))
f.close()

val_up = top(top.s)
val_bot = bottom(bottom.s)

plot(nacax,nacay)
plot(val_up[:,0],val_up[:,1],'ko-')
plot(val_bot[:,0],-val_bot[:,1],'ks--')
show()

sys.exit(0)

n = len(nacax)
k = 4
Nctl = 25
z = linspace(0,1,n)
X = vstack([nacax,nacay]).T.astype('d')
foil = pySpline.linear_spline(task='lms',k=4,x=nacax,y=nacay,z=z,Nctl=Nctl,
                              niter=1000,rel_tol=5e-4)

val = foil(foil.s)
plot(nacax,nacay)
plot(val[:,0],val[:,1],'ko-')
show()

# Put the points back to to naca0012_mod.dat
f=open('naca0012_mod.dat','w')
for i in xrange(n):
    f.write('%f %f\n'%(val[i,0],val[i,1]))
f.close()

sys.exit(0)



coef = zeros((Nctl,2))

#s =linspace(0,1,n)

s = zeros(n)
for i in xrange(n-1):
    s[i+1] = s[i] + sqrt(dot(X[i+1]-X[i],X[i+1]-X[i]))
s /= s[-1]

s0 = copy.deepcopy(s)
knots = pyspline.knots(s,Nctl,k)
knots0 = copy.deepcopy(knots)

length = pyspline.poly_length(X)
timeA = time.time()
s,knots,coef2 = pyspline.compute_curve(s,X,knots,k,Nctl)
print 'Fortran Time:',time.time()-timeA
# Check the starting rms with lsm
tot =0.0
for i in xrange(n):
    res = pyspline.eval_curve(s[i],knots,k,coef2)
    tot = tot + (res[0] - X[i,0])**2 + (res[1] - X[i,1])**2
# end for
print 'rms_fortran:',sqrt(tot/n)
print


# ------------------- PETSc Implementation -----------------
global old_rnorm
def converge_test(ksp,iter,rnorm):
    # Do a relative norm and a happy breakdown absolut norm
    global old_rnorm
    if iter== 0:
        old_rnorm = rnorm
        print 'rnorm0',rnorm
        if rnorm < 1e-12:
            return 1
        else:
            return 
    else:
        val = abs((rnorm-old_rnorm)/old_rnorm)
        old_rnorm = rnorm
        if val < 1e-12 or rnorm < 1e-12:
            return 1
        # end if
    # end if

J = PETSc.Mat().create(PETSc.COMM_SELF)
J.setSizes([n,Nctl])
J.setType(PETSc.Mat.Type.SEQAIJ)
rhs = PETSc.Vec().createSeq(n)
temp = PETSc.Vec().createSeq(Nctl)

ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)
ksp.getPC().setType('none')
ksp.setType('lsqr')
ksp.setConvergenceTest(converge_test)

s = copy.deepcopy(s0)
knots = copy.deepcopy(knots0)
timeA = time.time()
global old_rms
old_rms = sqrt(tot/n)
for j in xrange(1000):
    J.zeroEntries()
    rows,cols,vals = pyspline.curve_jacobian_linear(knots,k,s,Nctl)
    for i in xrange(len(rows)):
        J.setValue(rows[i],cols[i], vals[i],PETSc.InsertMode.INSERT_VALUES)
    # end for
    J.assemble()
    ksp.setOperators(J)

    for idim in xrange(2):
        rhs[:] = X[:,idim]
        ksp.solve(rhs, temp)
        coef[:,idim] = temp
    # end if
    s,rms = pyspline.curve_para_corr(knots,k,s.copy(),coef,length,X)

    if j == 0:
        old_rms = rms
    else:
        val = abs((rms-old_rms)/old_rms)
        old_rms = rms
        if val < 1e-4 or rms <1e-12:
            break
        # end if
    # en dif
# end for
        
# Check End Rms      
tot =0.0
for i in xrange(n):
    res = pyspline.eval_curve(s[i],knots,k,coef)
    tot = tot + (res[0] - X[i,0])**2 + (res[1] - X[i,1])**2
# end for
print 'actual fucking rms %d'%(j),sqrt(tot/n)


# end for
print 'PETSc iter:',j

print 'PETSc time:',time.time()-timeA




sys.exit(0)

# coef = zeros((3,3))
# coef[0,:] = [1,0,1]
# coef[1,:] = [1,1,sqrt(2)/2]
# coef[2,:] = [0,1,1]
# knots = [0,0,0,pi/2,pi/2,pi/2]
# k=3

# # Check Tangent
# val = pyspline.eval_curve_deriv(0.0,knots,k,coef)
# print 'val',val



# Plotting Stuff
n=25
#s = linspace(0,pi/2,n)
#s = linspace(0,1,n)
val = zeros((len(s),2))
for i in xrange(len(s)):
    val[i] = pyspline.eval_curve(s[i],knots,k,coef)
# end for
plot(val[:,0],val[:,1],'ko-')
plot(coef[:,0],coef[:,1],'go')

# Plot an ACTUAL circle
theta = linspace(0,pi/2,n)
x = cos(theta)
y = sin(theta)
plot(x,y,'ks--')
#plot(nacax,nacay,'ks--')
show()



