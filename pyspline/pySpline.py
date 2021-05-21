"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# Local modules
from . import libspline  # noqa: F401
from .pyCurve import Curve  # noqa: F401
from .pySurface import Surface  # noqa: F401
from .pyVolume import Volume  # noqa: F401
from .utils import (  # noqa: F401
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
