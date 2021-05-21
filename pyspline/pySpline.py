"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# Local modules
from . import libspline  # noqa: F401
from .pyCurve import Curve
from .pySurface import Surface
from .pyVolume import Volume
from .utils import (
    bilinearSurface,
    checkInput,
    closeTecplot,
    line,
    openTecplot,
    trilinearVolume,
    writeTecplot1D,
    writeTecplot2D,
    writeTecplot3D,
)
