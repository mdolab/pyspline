import sys, os
cur_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(cur_dir + '/python')
import pySpline
__all__ = ['pySpline']
