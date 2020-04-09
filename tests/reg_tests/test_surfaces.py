# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy
import unittest

# =============================================================================
# Extension modules
# =============================================================================
from pyspline import pySpline
from baseclasses import BaseRegTest


def eval_test(surface, handler):
    '''Eval fixed points from the surface'''
    # Evaluations are only good to about 1e-10 since there is fitting
    # involved

    #----------- Evaluation and derivative functions ---------------
    pts = [[0,0],[1,0],[0,1],[1,1],[.25,.25],[.75,.25]]
    for pt in pts:
        # print('Testing pt (%f %f)'%(pt[0],pt[1]))
        # print('Value:')
        handler.root_add_val(surface(pt[0],pt[1]),1e-10,1e-10)
        # print('Deriv:')
        handler.root_add_val(surface.getDerivative(pt[0],pt[1]),1e-10,1e-10)
        # print('Second Derivative')
        handler.root_add_val(surface.getSecondDerivative(pt[0],pt[1]),1e-8,1e-8)

    # print('Orig values at each corner')
    if surface.origData:
        handler.root_add_val(surface.getOrigValueCorner(0))
        handler.root_add_val(surface.getOrigValueCorner(1))
        handler.root_add_val(surface.getOrigValueCorner(2))
        handler.root_add_val(surface.getOrigValueCorner(3))
        
    # print('Orig values on edges')
    if surface.origData:
        handler.root_add_val(surface.getOrigValuesEdge(0))
        handler.root_add_val(surface.getOrigValuesEdge(1))
        handler.root_add_val(surface.getOrigValuesEdge(2))
        handler.root_add_val(surface.getOrigValuesEdge(3))

    # print('getValueEdge:')
    handler.root_add_val(surface.getValueEdge(0, 0.25))
    handler.root_add_val(surface.getValueEdge(0, 0.75))
    handler.root_add_val(surface.getValueEdge(1, 0.25))
    handler.root_add_val(surface.getValueEdge(1, 0.75))
    handler.root_add_val(surface.getValueEdge(2, 0.25))
    handler.root_add_val(surface.getValueEdge(2, 0.75))
    handler.root_add_val(surface.getValueEdge(3, 0.25))
    handler.root_add_val(surface.getValueEdge(3, 0.75))

def run_surface_test(surface, handler):
    ''' This function is used to test the functions that are apart of
    the curve class. They operate on the 'curve' that is passed. '''

    # Test the evaluations
    eval_test(surface, handler)

  
    # Test the windowing (same surface)
    surf2 = surface.windowSurface([0,0],[1,1])

    # print('Evaluations on surf2 should be same as original:')
    eval_test(surf2, handler)

    surf2 = surface.windowSurface([0.25,.25],[.75,.5])

    # print('These points should be the same:')
    handler.root_add_val(surface(0.25,.25))
    handler.root_add_val(surf2(0,0))

    # print('These points should be the same:')
    handler.root_add_val(surface(0.75,.5))
    handler.root_add_val(surf2(1,1))

    # print('Test get bounds')
    handler.root_add_val(surface.getBounds())

def run_project_test(surface, handler):
    # Run a bunch of point projections: Only try to match to 1e-8
    eps = 1e-8

    # print('------------- These points should be fully inside of domain')
    pts= [[0,0,0],[2,3,-1],[3,2.5,-.1]]
    for pt in pts:
        # print('Projecting point (%f %f %f)'%(pt[0],pt[1],pt[2]))
        u,v,D = surface.projectPoint(pt,eps=1e-12)
        # print('u:')
        handler.root_add_val(u, eps, eps)
        # print('v:')
        handler.root_add_val(v, eps, eps)
        # print('D:')
        handler.root_add_val(D, eps*10, eps*10)

    # print(' ----------- This should be (0,0) corner')
    u,v,D = surface.projectPoint([-1,-1,0],eps=1e-12)
    # print('u:')
    handler.root_add_val(u, eps, eps)
    # print('v:')
    handler.root_add_val(v, eps, eps)

    # print(' ---------- This should be (0,1) corner')
    u,v,D = surface.projectPoint([-1,5,0],eps=1e-12)
    # print('u:')
    handler.root_add_val(u, eps, eps)
    # print('v:')
    handler.root_add_val(v, eps, eps)

    # print(' ---------- This should be (1,0) corner')
    u,v,D = surface.projectPoint([6,-1,0],eps=1e-12)
    # print('u:')
    handler.root_add_val(u, eps, eps)
    # print('v:')
    handler.root_add_val(v, eps, eps)

    # print(' ---------- This should be (1,1) corner')
    u,v,D = surface.projectPoint([6,6,0],eps=1e-12)
    # print('u:')
    handler.root_add_val(u, eps, eps)
    # print('v:')
    handler.root_add_val(v, eps, eps)

    # print(' ---------- This should be  edge zero (*,0)')
    u,v,D = surface.projectPoint([2.54,-1,0],eps=1e-12)
    # print('u:')
    handler.root_add_val(u, eps, eps)
    # print('v:')
    handler.root_add_val(v,eps, eps)


    # Curve projection
    for kc in [2,3,4]:
        x = [0,1,2,0]
        y = [4,3,2,1]
        z = [-3,1,3,5]

        curve = pySpline.Curve(k=kc,x=x,y=y,z=z)
        u,v,s,D = surface.projectCurve(curve)
        # print(' ---------- surface-curve projection with kc=%d'%(kc))
        # print('u:')
        handler.root_add_val(u, eps, eps)
        # print('v:')
        handler.root_add_val(v, eps, eps)
        # print('s:')
        handler.root_add_val(s, eps, eps)
        # print('D:')
        handler.root_add_val(D, eps*10, eps*10)


def io_test(surface, handler):
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


class Test(unittest.TestCase):
    
    def setUp(self):
        self.ref_file = 'ref/test_surfaces.ref'

    def train(self):
        with BaseRegTest(self.ref_file, train=True) as handler:
            self.regression_test(handler)

    def test(self):
        with BaseRegTest(self.ref_file, train=False) as handler:
            self.regression_test(handler)
    


    def regression_test(self, handler, solve=False):
        
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
                        # print('+'+'-'*78+'+')
                        # print(' '*20 + 'Testing Surface with ku=%d, kv=%d, nCtlu=%d,\
         # nCtlv=%d'%(ku,kv,nCtlu,nCtlv))
                        # print('+'+'-'*78+'+')
                        surface = pySpline.Surface(x=U, y=V, z=Z, ku=ku, kv=kv, 
                                                   nCtlu=nCtlu, nCtlv=nCtlv)
                        surface.recompute()
                        run_surface_test(surface, handler)
                        run_project_test(surface, handler)
                        io_test(surface, handler)
