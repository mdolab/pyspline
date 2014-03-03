from __future__ import print_function
# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os, sys

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_regression_helper import *
from pyspline import pySpline

def eval_test(surface):
    '''Eval fixed points from the surface'''
    # Evaluations are only good to about 1e-10 since there is fitting
    # involved

    #----------- Evaluation and derivative functions ---------------
    pts = [[0,0],[1,0],[0,1],[1,1],[.25,.25],[.75,.25]]
    for pt in pts:
        print('Testing pt (%f %f)'%(pt[0],pt[1]))
        print('Value:')
        reg_write(surface(pt[0],pt[1]),1e-10,1e-10)
        print('Deriv:')
        reg_write(surface.getDerivative(pt[0],pt[1]),1e-10,1e-10)
        print('Second Derivative')
        reg_write(surface.getSecondDerivative(pt[0],pt[1]),1e-8,1e-8)

    print('Orig values at each corner')
    if surface.origData:
        reg_write(surface.getOrigValueCorner(0))
        reg_write(surface.getOrigValueCorner(1))
        reg_write(surface.getOrigValueCorner(2))
        reg_write(surface.getOrigValueCorner(3))
        
    print('Orig values on edges')
    if surface.origData:
        reg_write(surface.getOrigValuesEdge(0))
        reg_write(surface.getOrigValuesEdge(1))
        reg_write(surface.getOrigValuesEdge(2))
        reg_write(surface.getOrigValuesEdge(3))

    print('getValueEdge:')
    reg_write(surface.getValueEdge(0, 0.25))
    reg_write(surface.getValueEdge(0, 0.75))
    reg_write(surface.getValueEdge(1, 0.25))
    reg_write(surface.getValueEdge(1, 0.75))
    reg_write(surface.getValueEdge(2, 0.25))
    reg_write(surface.getValueEdge(2, 0.75))
    reg_write(surface.getValueEdge(3, 0.25))
    reg_write(surface.getValueEdge(3, 0.75))


def run_surface_test(surface):
    ''' This function is used to test the functions that are apart of
    the curve class. They operate on the 'curve' that is passed. '''

    # Test the evaluations
    eval_test(surface)

  
    # Test the windowing (same surface)
    surf2 = surface.windowSurface([0,0],[1,1])

    print('Evaluations on surf2 should be same as original:')
    eval_test(surf2)

    surf2 = surface.windowSurface([0.25,.25],[.75,.5])

    print('These points should be the same:')
    reg_write(surface(0.25,.25))
    reg_write(surf2(0,0))

    print('These points should be the same:')
    reg_write(surface(0.75,.5))
    reg_write(surf2(1,1))

    print('Test get bounds')
    reg_write(surface.getBounds())
    

def run_project_test(surface):
    # Run a bunch of point projections: Only try to match to 1e-8
    eps = 1e-8

    print('------------- These points should be fully inside of domain')
    pts= [[0,0,0],[2,3,-1],[3,2.5,-.1]]
    for pt in pts:
        print('Projecting point (%f %f %f)'%(pt[0],pt[1],pt[2]))
        u,v,D = surface.projectPoint(pt,eps=1e-12)
        print('u:')
        reg_write(u, eps, eps)
        print('v:')
        reg_write(v, eps, eps)
        print('D:')
        reg_write(D, eps*10, eps*10)

    print(' ----------- This should be (0,0) corner')
    u,v,D = surface.projectPoint([-1,-1,0],eps=1e-12)
    print('u:')
    reg_write(u, eps, eps)
    print('v:')
    reg_write(v, eps, eps)

    print(' ---------- This should be (0,1) corner')
    u,v,D = surface.projectPoint([-1,5,0],eps=1e-12)
    print('u:')
    reg_write(u, eps, eps)
    print('v:')
    reg_write(v, eps, eps)

    print(' ---------- This should be (1,0) corner')
    u,v,D = surface.projectPoint([6,-1,0],eps=1e-12)
    print('u:')
    reg_write(u, eps, eps)
    print('v:')
    reg_write(v, eps, eps)

    print(' ---------- This should be (1,1) corner')
    u,v,D = surface.projectPoint([6,6,0],eps=1e-12)
    print('u:')
    reg_write(u, eps, eps)
    print('v:')
    reg_write(v, eps, eps)

    print(' ---------- This should be  edge zero (*,0)')
    u,v,D = surface.projectPoint([2.54,-1,0],eps=1e-12)
    print('u:')
    reg_write(u, eps, eps)
    print('v:')
    reg_write(v,eps, eps)


    # Curve projection
    for kc in [2,3,4]:
        x = [0,1,2,0]
        y = [4,3,2,1]
        z = [-3,1,3,5]

        curve = pySpline.Curve(k=kc,x=x,y=y,z=z)
        u,v,s,D = surface.projectCurve(curve)
        print(' ---------- surface-curve projection with kc=%d'%(kc))
        print('u:')
        reg_write(u, eps, eps)
        print('v:')
        reg_write(v, eps, eps)
        print('s:')
        reg_write(s, eps, eps)
        print('D:')
        reg_write(D, eps*10, eps*10)
 
def io_test(surface):
    '''Test the writing functions'''
    surface.writeTecplot('tmp.dat', surf=True, coef=True, orig=True,
                         directions=True)
    f = open('tmp.dat','w')
    # These three calls, are private functions normally only called
    # from pyGeo. We are not checking their output, rather just making
    # sure they run. 
    surface.writeIGES_directory(f, 0, 0)
    surface.writeIGES_directory(f, 0, 0)
    surface.writeTin(f)
    os.remove('tmp.dat')

    return

if __name__ == '__main__':
    # Create a generic surface
    nu = 10
    nv = 10
    u = numpy.linspace(0,4,nu)
    v = numpy.linspace(0,4,nv)
    [V,U] = numpy.meshgrid(v,u)
    Z = numpy.cos(U)*numpy.sin(V)


    for ku in [2, 3, 4]:
        for kv in [2, 3, 4]:
            for nCtlu in [5, 10]:
                for nCtlv in [5, 10]:
                    print('+'+'-'*78+'+')
                    print(' '*20 + 'Testing Surface with ku=%d, kv=%d, nCtlu=%d,\
     nCtlv=%d'%(ku,kv,nCtlu,nCtlv))
                    print('+'+'-'*78+'+')
                    surface = pySpline.Surface(x=U, y=V, z=Z, ku=ku, kv=kv, 
                                               nCtlu=nCtlu, nCtlv=nCtlv)
                    surface.recompute()
                    run_surface_test(surface)
                    run_project_test(surface)
                    io_test(surface)

