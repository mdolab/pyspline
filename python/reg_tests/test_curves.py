# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os, sys, argparse

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_regression_helper import *
sys.path.append('../')
import pySpline
import numpy

def eval_test(curve):
    '''Eval fixed points from the curve'''
    # ----------- Evaluation and derivative functions ---------------
    pts = [0.0, 0.5, 0.75, 1.0]
    for pt in pts:
        print 'Testing pt %f'%(pt)
        print 'Value:'
        reg_write(curve(0.0))
        print 'Deriv:'
        reg_write(curve.getDerivative(0))
        print 'Second Derivative'
        reg_write(curve.getSecondDerivative(0))
    # end for


def run_curve_test(curve): 
    ''' This function is used to test the functions that are apart of
    the curve class. They operate on the 'curve' that is passed. '''

    # Test the evaluations
    eval_test(curve)

    # Inset a knot at 0.25 with multiplicity of k-1 and retest
    curve.insertKnot(0.25,curve.k-1)
    print '------- Cruve with inserted knots ----------'
    eval_test(curve)
    
    print 'Curve length:'
    reg_write(curve.getLength())

    curve_windowed = curve.windowCurve(0.1, 0.90)
    print 'These points should be the same (window curve test)'
    reg_write(curve(0.2))
    reg_write(curve_windowed((0.2-0.1)/(0.90-0.1)))

    # Split the curve in two:
    curve1, curve2 = curve.splitCurve(0.5)
    print 'These three points should be the same (split curve test)'
    reg_write(curve1(1.0))
    reg_write(curve2(0.0))
    reg_write(curve(0.5))

    # Reverse test:
    print 'These two points should be the same (reverse test)'
    reg_write(curve(0.4))
    curve.reverse()
    reg_write(curve(0.6))

def run_project_test(curve):
    # Run point projection and curve projection tests'
    pts = [[0.4,1.5,1.5],[-.1,0.5,1.8], curve(0.75)]
    s, D = curve.projectPoint(pts)
    for i in xrange(len(s)):
        print 'Project point %f %f %f'%(pts[i][0],pts[i][1],pts[i][2])
        reg_write(s[i])
        reg_write(D[i])


def io_test(curve):
    '''Test the writing functions'''
    curve.writeTecplot('tmp.dat')
    curve.writeTecplot('tmp.dat',coef=False, orig=False)
    os.remove('tmp.dat')

print '+--------------------------------------+'
print '           Create Tests   '
print '+--------------------------------------+'
print '-----------  2D k=2 test ----------'
k = 2
t = [0,0,0.5,1,1]
coef = numpy.zeros((3,2))
coef[0] = [0,0.6]
coef[1] = [1.1,1.4]
coef[2] = [2.6,5.1]
curve = pySpline.curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)

print '-----------  2D k=3 test ----------'
k = 3
t = [0,0,0,0.5,1,1,1]
coef = numpy.zeros((4,2))
coef[0] = [0,0.45]
coef[1] = [.71,1.5]
coef[2] = [2.5,5.9]
coef[3] = [4,-2]
curve = pySpline.curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)

print '-----------  2D k=4 test ----------'
k = 4
t = [0,0,0,0,0.5,1,1,1,1]
coef = numpy.zeros((5,2))
coef[0] = [0,-.60]
coef[1] = [.9,1.6]
coef[2] = [1.6,5.2]
coef[3] = [4.2,-2.24]
coef[4] = [2.9,6.2]
curve = pySpline.curve(t=t,k=k,coef=coef)
run_curve_test(curve)
io_test(curve)


# Get helix data
n = 100
theta = numpy.linspace(0.0, numpy.pi*2, n)
x = numpy.cos(theta)
y = numpy.sin(theta)
z = numpy.linspace(0,1,n)

print '+--------------------------------------+'
print '              LMS Tests   '
print '+--------------------------------------+'
for k in [2,3,4]:
    print '--------- Test helix data with k=%d-------'%(k)
    curve = pySpline.curve(x=x,y=y,z=z,k=k,Nctl=16,niter=50)
    run_curve_test(curve)
    run_project_test(curve)

print '+--------------------------------------+'
print '           Interp Tests   '
print '+--------------------------------------+'
for k in [2,3,4]:
    print '--------- Test helix data with k=%d-------'%(k)
    curve = pySpline.curve(x=x,y=y,z=z,k=k)
    run_curve_test(curve)
    run_project_test(curve)


