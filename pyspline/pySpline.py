"""
pySpline
--------

Contains classes for working with B-spline :class:`Curve`, :class:`Surface` and
:class:`Volume`
"""

# External modules
import numpy as np

# Local modules
from .pyCurve import Curve
from .pySurface import Surface
from .pyVolume import Volume
from .utils import Error

# ----------------------------------------------------------------------
#                     Misc Helper Functions
# ----------------------------------------------------------------------


def trilinearVolume(*args):
    """This is a short-cut function to create a trilinear b-spline
    volume. It can be created with ``x`` **OR** with ``xmin`` and
    ``xmax``.

    Parameters
    ----------
    x : array of size (2, 2, 2, 3)
        Coordinates of the corners of the box.

    xmin : array of size (3)
        The extreme lower corner of the box
    xmax : array of size (3)
        The extreme upper corner of the box. In this case, by
        construction, the box will be coordinate axis aligned.
    """
    tu = [0, 0, 1, 1]
    tv = [0, 0, 1, 1]
    tw = [0, 0, 1, 1]
    ku = 2
    kv = 2
    kw = 2

    if len(args) == 1:
        return Volume(coef=args[0], tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    elif len(args) == 2:
        xmin = np.array(args[0]).astype("d")
        xmax = np.array(args[1]).astype("d")

        xLow = xmin[0]
        xHigh = xmax[0]
        yLow = xmin[1]
        yHigh = xmax[1]
        zLow = xmin[2]
        zHigh = xmax[2]

        coef = np.zeros((2, 2, 2, 3))
        coef[0, 0, 0, :] = [xLow, yLow, zLow]
        coef[1, 0, 0, :] = [xHigh, yLow, zLow]
        coef[0, 1, 0, :] = [xLow, yHigh, zLow]
        coef[1, 1, 0, :] = [xHigh, yHigh, zLow]
        coef[0, 0, 1, :] = [xLow, yLow, zHigh]
        coef[1, 0, 1, :] = [xHigh, yLow, zHigh]
        coef[0, 1, 1, :] = [xLow, yHigh, zHigh]
        coef[1, 1, 1, :] = [xHigh, yHigh, zHigh]

        return Volume(coef=coef, tu=tu, tv=tv, tw=tw, ku=ku, kv=kv, kw=kw)
    else:
        raise Error("An unknown number of arguments was passed to trilinear  Volume")


def bilinearSurface(*args):
    """This is short-cut function to create a bilinear surface.

    Args can contain:

    1.  ``x`` array of size(4, 3).  The four corners of the array
        arranged in the coordinate direction orientation::

          2          3
          +----------+
          |          |
          |          |
          |          |
          +----------+
          0          1

    2. ``p1``, ``p2``, ``p3``, ``p4``. Individual corner points in CCW Ordering::

          3          2
          +----------+
          |          |
          |          |
          |          |
          +----------+
          0          1
    """
    if len(args) == 1:
        # One argument passed in ... assume its X
        if len(args[0]) != 4:
            raise Error("A single argument passed to bilinear surface must contain 4 points and be of size (4, 3)")
        coef = np.zeros((2, 2, 3))
        coef[0, 0] = args[0][0]
        coef[1, 0] = args[0][1]
        coef[0, 1] = args[0][2]
        coef[1, 1] = args[0][3]
        return Surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)
    else:
        # Assume 4 arguments
        coef = np.zeros([2, 2, 3])
        coef[0, 0] = args[0]
        coef[1, 0] = args[1]
        coef[0, 1] = args[3]
        coef[1, 1] = args[2]
        return Surface(coef=coef, tu=[0, 0, 1, 1], tv=[0, 0, 1, 1], ku=2, kv=2)


def line(*args, **kwargs):
    """This is a short cut function to create a line as a pySpline Curve object.

    Args can contain: (pick one)

    1. ``X`` array of size(2, ndim). The two end points
    2. ``x1``, ``x2`` The two end points (each of size ndim)
    3. ``x```, dir=direction. A point and a displacement vector
    4. ``x1``, dir=direction, length=length. As 3. but with a specific length
    """

    if len(args) == 2:
        # Its a two-point type
        return Curve(coef=[args[0], args[1]], k=2, t=[0, 0, 1, 1])
    elif len(args) == 1:
        if len(args[0]) == 2:  # its X
            return Curve(coef=args[0], k=2, t=[0, 0, 1, 1])
        elif "dir" in kwargs:
            # We have point and direction
            if "length" in kwargs:
                x2 = args[0] + kwargs["dir"] / np.linalg.norm(kwargs["dir"]) * kwargs["length"]
            else:
                x2 = args[0] + kwargs["dir"]

            return Curve(coef=[args[0], x2], k=2, t=[0, 0, 1, 1])
        else:
            Error("Error: dir must be specified if only 1 argument is given")
