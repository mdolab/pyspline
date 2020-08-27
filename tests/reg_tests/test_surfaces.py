# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from numpy.testing import assert_allclose
import unittest

# =============================================================================
# Extension modules
# =============================================================================
from pyspline import pySpline
from baseclasses import BaseRegTest


def eval_test(surface, handler, test_name):
    '''Eval fixed points from the surface'''
    # Evaluations are only good to about 1e-10 since there is fitting
    # involved

    #----------- Evaluation and derivative functions ---------------
    pts = [[0,0],[1,0],[0,1],[1,1],[.25,.25],[.75,.25]]
    for pt in pts:
        # print('Testing pt (%f %f)'%(pt[0],pt[1]))
        handler.root_add_val(surface(pt[0],pt[1]), '{} evaluate point at {}'.format(test_name, pt), 1e-10,1e-10)
        handler.root_add_val(surface.getDerivative(pt[0], pt[1]), '{} evaluate point deriv at {}'.format(test_name, pt), 1e-10,1e-10)
        handler.root_add_val(surface.getSecondDerivative(pt[0], pt[1]), '{} evaluate point second deriv at {}'.format(test_name, pt), 1e-8,1e-8)

    # Orig values at each corner
    if surface.origData:
        for i in range(4):
            handler.root_add_val(surface.getOrigValueCorner(i), '{} evaluate corner {}'.format(test_name, i))

    # Orig values on edges
    if surface.origData:
        for i in range(4):
            handler.root_add_val(surface.getOrigValuesEdge(i), '{} evaluate edge {}'.format(test_name, i))

    # getValueEdge
    for i in [0, 1, 2, 3]:
        for j in [0.25, 0.75]:
            handler.root_add_val(surface.getValueEdge(i, j), '{} getValueEdge({}, {})'.format(test_name, i, j))

def run_surface_test(surface, handler, test_name):
    ''' This function is used to test the functions that are apart of
    the curve class. They operate on the 'curve' that is passed. '''

    # Test the evaluations
    eval_test(surface, handler, test_name)

    # Test the windowing (same surface)
    surf2 = surface.windowSurface([0,0],[1,1])

    # Evaluations on surf2 should be same as original
    # we use the same test_name to check the same entries
    eval_test(surf2, handler, test_name)

    surf2 = surface.windowSurface([0.25,.25],[.75,.5])
    # these values should be the same
    assert_allclose(surface(0.25,.25), surf2(0,0))
    assert_allclose(surface(0.75,.5), surf2(1,1))

    # Test get bounds
    handler.root_add_val(surface.getBounds(), '{} bounds'.format(test_name))

def run_project_test(surface, handler, test_name):
    # Run a bunch of point projections: Only try to match to 1e-8
    eps = 1e-8

    # print('------------- These points should be fully inside of domain')
    pts= [[0,0,0],[2,3,-1],[3,2.5,-.1]]
    for pt in pts:
        # print('Projecting point (%f %f %f)'%(pt[0],pt[1],pt[2]))
        u,v,D = surface.projectPoint(pt,eps=1e-12)
        handler.root_add_val(u, '{} point {} projection u'.format(test_name, pt), eps, eps)
        handler.root_add_val(v, '{} point {} projection v'.format(test_name, pt), eps, eps)
        handler.root_add_val(D, '{} point {} projection D'.format(test_name, pt), eps*10, eps*10)

    # ----------- This should be (0,0) corner
    u,v,D = surface.projectPoint([-1,-1,0],eps=1e-12)
    handler.root_add_val(u, '{} projected u for (0,0) corner'.format(test_name), eps, eps)
    handler.root_add_val(v, '{} projected v for (0,0) corner'.format(test_name), eps, eps)

    # ---------- This should be (0,1) corner
    u,v,D = surface.projectPoint([-1,5,0],eps=1e-12)
    handler.root_add_val(u, '{} projected u for (0,1) corner'.format(test_name), eps, eps)
    handler.root_add_val(v, '{} projected v for (0,1) corner'.format(test_name), eps, eps)

    # ---------- This should be (1,0) corner
    u,v,D = surface.projectPoint([6,-1,0],eps=1e-12)
    handler.root_add_val(u, '{} projected u for (1,0) corner'.format(test_name), eps, eps)
    handler.root_add_val(v, '{} projected v for (1,0) corner'.format(test_name), eps, eps)

    # ---------- This should be (1,1) corner
    u,v,D = surface.projectPoint([6,6,0],eps=1e-12)
    handler.root_add_val(u, '{} projected u for (1,1) corner'.format(test_name), eps, eps)
    handler.root_add_val(v, '{} projected v for (1,1) corner'.format(test_name), eps, eps)

    # ---------- This should be edge zero (*,0)
    u,v,D = surface.projectPoint([2.54,-1,0],eps=1e-12)
    handler.root_add_val(u, '{} projected u for (*,0) edge'.format(test_name), eps, eps)
    handler.root_add_val(v, '{} projected v for (*,0) edge'.format(test_name), eps, eps)


    # Curve projection
    for kc in [2,3,4]:
        x = [0,1,2,0]
        y = [4,3,2,1]
        z = [-3,1,3,5]

        curve = pySpline.Curve(k=kc,x=x,y=y,z=z)
        u,v,s,D = surface.projectCurve(curve)
        # ---------- surface-curve projection with kc = kc
        handler.root_add_val(u, '{} projected curve u with kc={}'.format(test_name, kc), eps, eps)
        handler.root_add_val(v, '{} projected curve v with kc={}'.format(test_name, kc), eps, eps)
        handler.root_add_val(s, '{} projected curve s with kc={}'.format(test_name, kc), eps, eps)
        handler.root_add_val(D, '{} projected curve D with kc={}'.format(test_name, kc), eps*10, eps*10)


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
            handler.writeRef()

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

        # Testing Surface with ku, kv, nCtlu, nCtlv
        for ku in [2, 3, 4]:
            for kv in [2, 3, 4]:
                for nCtlu in [5, 10]:
                    for nCtlv in [5, 10]:
                        test_name = 'surface with ku={}, kv={}, nCtlu={}, nCtlv={}'.format(ku,kv,nCtlu,nCtlv)
                        surface = pySpline.Surface(x=U, y=V, z=Z, ku=ku, kv=kv, 
                                                   nCtlu=nCtlu, nCtlv=nCtlv)
                        surface.recompute()
                        run_surface_test(surface, handler, test_name)
                        run_project_test(surface, handler, test_name)
                        io_test(surface, handler)
