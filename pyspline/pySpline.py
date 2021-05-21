"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# Local modules
from . import libspline  # noqa: F401
from .pySurface import Surface
from .pyVolume import Volume
from .pyCurve import Curve
from .utils import (
    writeTecplot1D,
    writeTecplot2D,
    writeTecplot3D,
    openTecplot,
    closeTecplot,
    checkInput,
    trilinearVolume,
    bilinearSurface,
    line,
)
