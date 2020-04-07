# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_regression_helper import *
from pyspline import pySpline

def eval_test(crv):
    '''Eval fixed points from the curve'''
    # ----------- Evaluation and derivative functions ---------------
    pts = [0.0, 0.5, 0.75, 1.0]
    for pt in pts:
        print('Testing pt %f'%(pt))
        print('Value:')
        reg_write(crv(pt))
        print('Deriv:')
        reg_write(crv.getDerivative(pt))
        print('Second Derivative')
        reg_write(crv.getSecondDerivative(pt),1e-10,1e-10)

def run_curve_test(crv): 
    ''' This function is used to test the functions that are apart of
    the curve class. They operate on the 'crv' that is passed. '''

    # Test the evaluations
    eval_test(crv)

    # Inset a knot at 0.25 with multiplicity of k-1 and retest
    crv.insertKnot(0.25,crv.k-1)
    print('------- Cruve with inserted knots ----------')
    eval_test(crv)
    
    print('Curve length:')
    reg_write(crv.getLength())

    curve_windowed = crv.windowCurve(0.1, 0.90)
    print('These points should be the same (window curve test)')
    reg_write(crv(0.2))
    reg_write(curve_windowed((0.2-0.1)/(0.90-0.1)))

    # Split the curve in two:
    curve1, curve2 = crv.splitCurve(0.5)
    print('These three points should be the same (split curve test)')
    reg_write(curve1(1.0))
    reg_write(curve2(0.0))
    reg_write(crv(0.5))

    # Reverse test:
    print('These two points should be the same (reverse test)')
    reg_write(crv(0.4))
    crv.reverse()
    reg_write(crv(0.6))

def run_project_test(crv):
    # Run point projection and curve projection tests'
    pts = [[0.4,1.5,1.5],[-.1,0.5,1.8], crv(0.75)]
    # Default tolerance is 1e-10. so only check to 1e-9
    s, D = crv.projectPoint(pts)
    for i in range(len(s)):
        print('Project point %f %f %f'%(pts[i][0],pts[i][1],pts[i][2]))
        reg_write(s[i],1e-9,1e-9)
        reg_write(D[i],1e-9,1e-9)


def io_test(crv):
    '''Test the writing functions'''
    crv.writeTecplot('tmp.dat')
    crv.writeTecplot('tmp.dat',coef=False, orig=False)
    os.remove('tmp.dat')

print('+--------------------------------------+')
print('           Create Tests   ')
print('+--------------------------------------+')
print('-----------  2D k=2 test ----------')
k = 2
t = [0,0,0.5,1,1]
coef = numpy.zeros((3,2))
coef[0] = [0,0.6]
coef[1] = [1.1,1.4]
coef[2] = [2.6,5.1]
curve = pySpline.Curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)

print('-----------  2D k=3 test ----------')
k = 3
t = [0,0,0,0.5,1,1,1]
coef = numpy.zeros((4,2))
coef[0] = [0,0.45]
coef[1] = [.71,1.5]
coef[2] = [2.5,5.9]
coef[3] = [4,-2]
curve = pySpline.Curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)

print('-----------  2D k=4 test ----------')
k = 4
t = [0,0,0,0,0.5,1,1,1,1]
coef = numpy.zeros((5,2))
coef[0] = [0,-.60]
coef[1] = [.9,1.6]
coef[2] = [1.6,5.2]
coef[3] = [4.2,-2.24]
coef[4] = [2.9,6.2]
curve = pySpline.Curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)


# Get helix data
n = 100
theta = numpy.linspace(0.0, numpy.pi*2, n)
x = numpy.cos(theta)
y = numpy.sin(theta)
z = numpy.linspace(0,1,n)

print('+--------------------------------------+')
print('              LMS Tests   ')
print('+--------------------------------------+')
for k in [2,3,4]:
    print('--------- Test helix data with k=%d-------'%(k))
    curve = pySpline.Curve(x=x,y=y,z=z,k=k,nCtl=16,niter=50)
    run_curve_test(curve)
    run_project_test(curve)

print('+--------------------------------------+')
print('           Interp Tests   ')
print('+--------------------------------------+')
for k in [2,3,4]:
    print('--------- Test helix data with k=%d-------'%(k))
    curve = pySpline.Curve(x=x,y=y,z=z,k=k)
    run_curve_test(curve)
    run_project_test(curve)
