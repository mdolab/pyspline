# Standard Python modules
import os
import unittest

# External modules
from baseclasses import BaseRegTest
import numpy as np
from numpy.testing import assert_allclose

# First party modules
from pyspline import Curve

baseDir = os.path.dirname(os.path.abspath(__file__))


def eval_test(crv, handler, test_name):
    """Eval fixed points from the curve"""
    # ----------- Evaluation and derivative functions ---------------
    pts = [0.0, 0.5, 0.75, 1.0]
    for pt in pts:
        # print('Testing pt %f'%(pt))
        # print('Value:')
        handler.root_add_val("{} curve evaluated at {}".format(test_name, pt), crv(pt))
        # print('Deriv:')
        handler.root_add_val("{} curve derivative evaluated at {}".format(test_name, pt), crv.getDerivative(pt))
        # print('Second Derivative')
        handler.root_add_val(
            "{} curve second derivative evaluated at {}".format(test_name, pt), crv.getSecondDerivative(pt), tol=1e-10
        )


def run_curve_test(crv, handler, test_name):
    """This function is used to test the functions that are apart of
    the curve class. They operate on the 'crv' that is passed."""

    # Test the evaluations
    eval_test(crv, handler, "{} initial".format(test_name))

    # Inset a knot at 0.25 with multiplicity of k-1 and retest
    crv.insertKnot(0.25, crv.k - 1)
    # print('------- Cruve with inserted knots ----------')
    eval_test(crv, handler, "{} inserted knots".format(test_name))

    handler.root_add_val("{} curve length".format(test_name), crv.getLength())

    a = 0.1
    b = 0.9
    curve_windowed = crv.windowCurve(a, b)
    c = 0.2
    assert_allclose(crv(c), curve_windowed((c - a) / (b - a)), err_msg="points do not match after windowing!")

    # Split the curve in two:
    curve1, curve2 = crv.splitCurve(0.5)
    c1 = curve1(1.0)
    c2 = curve2(0.0)
    c3 = crv(0.5)
    # These three points should be the same (split curve test)
    assert_allclose(c1, c2, err_msg="points do not match after splitting curve")
    assert_allclose(c1, c3, err_msg="points do not match after splitting curve")

    # Reverse test:
    c1 = crv(0.4)
    crv.reverse()
    c2 = crv(0.6)
    assert_allclose(c1, c2, err_msg="points do not match after reversing curve")


def run_project_test(crv, handler, test_name):
    # Run point projection and curve projection tests'
    pts = [[0.4, 1.5, 1.5], [-0.1, 0.5, 1.8], crv(0.75)]
    # Default tolerance is 1e-10. so only check to 1e-9
    s, D = crv.projectPoint(pts)
    for i in range(len(s)):
        # print('Project point %f %f %f'%(pts[i][0],pts[i][1],pts[i][2]))
        handler.root_add_val("{} projection test for point {} solution".format(test_name, i), s[i], tol=1e-9)
        handler.root_add_val("{} projection test for point {} distance".format(test_name, i), D[i], tol=1e-9)


def io_test(crv):
    """Test the writing functions"""
    crv.writeTecplot("tmp_curves.dat")
    crv.writeTecplot("tmp_curves.dat", coef=False, orig=False)
    os.remove("tmp_curves.dat")


class Test(unittest.TestCase):
    def setUp(self):
        self.ref_file = os.path.join(baseDir, "ref/test_curves.ref")

    def train(self):
        with BaseRegTest(self.ref_file, train=True) as handler:
            self.regression_test(handler)
            handler.writeRef()

    def test(self):
        with BaseRegTest(self.ref_file, train=False) as handler:
            self.regression_test(handler)

    def regression_test(self, handler):

        # print('+--------------------------------------+')
        # print('           Create Tests   ')
        # print('+--------------------------------------+')
        # print('-----------  2D k=2 test ----------')
        k = 2
        t = [0, 0, 0.5, 1, 1]
        coef = np.zeros((3, 2))
        coef[0] = [0, 0.6]
        coef[1] = [1.1, 1.4]
        coef[2] = [2.6, 5.1]
        curve = Curve(t=t, k=k, coef=coef)
        run_curve_test(curve, handler, "2D k={}".format(k))
        io_test(curve)

        # print('-----------  2D k=3 test ----------')
        k = 3
        t = [0, 0, 0, 0.5, 1, 1, 1]
        coef = np.zeros((4, 2))
        coef[0] = [0, 0.45]
        coef[1] = [0.71, 1.5]
        coef[2] = [2.5, 5.9]
        coef[3] = [4, -2]
        curve = Curve(t=t, k=k, coef=coef)
        run_curve_test(curve, handler, "2D k={}".format(k))
        io_test(curve)

        # print('-----------  2D k=4 test ----------')
        k = 4
        t = [0, 0, 0, 0, 0.5, 1, 1, 1, 1]
        coef = np.zeros((5, 2))
        coef[0] = [0, -0.60]
        coef[1] = [0.9, 1.6]
        coef[2] = [1.6, 5.2]
        coef[3] = [4.2, -2.24]
        coef[4] = [2.9, 6.2]
        curve = Curve(t=t, k=k, coef=coef)
        run_curve_test(curve, handler, "2D k={}".format(k))
        io_test(curve)

        # Get helix data
        n = 100
        theta = np.linspace(0.0, np.pi * 2, n)
        x = np.cos(theta)
        y = np.sin(theta)
        z = np.linspace(0, 1, n)

        # print('+--------------------------------------+')
        # print('              LMS Tests   ')
        # print('+--------------------------------------+')
        for k in [2, 3, 4]:
            test_name = "LMS test k={}".format(k)
            # print('--------- Test helix data with k=%d-------'%(k))
            curve = Curve(x=x, y=y, z=z, k=k, nCtl=16, niter=50)
            run_curve_test(curve, handler, test_name)
            run_project_test(curve, handler, test_name)

        # print('+--------------------------------------+')
        # print('           Interp Tests   ')
        # print('+--------------------------------------+')
        for k in [2, 3, 4]:
            test_name = "interpolation test k={}".format(k)
            # print('--------- Test helix data with k=%d-------'%(k))
            curve = Curve(x=x, y=y, z=z, k=k)
            run_curve_test(curve, handler, test_name)
            run_project_test(curve, handler, test_name)
