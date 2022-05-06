#! /usr/bin/env python
# Standard Python modules
from ast import Import
import sys

print("Testing if module pyspline can be imported...")
try:
    # External modules
    import libspline  # noqa: F401
except ImportError as e:
    print("Error importing libspline.so")
    print(e)
    sys.exit(1)
# end try

print("Module libspline was successfully imported.")
