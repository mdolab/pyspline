#! /usr/bin/env python
# Standard Python modules
import sys

print("Testing if module pyspline can be imported...")
try:
    # External modules
    import libspline  # noqa: F401
except ImportError:
    print("Error importing libspline.so")
    sys.exit(1)
# end try

print("Module libspline was successfully imported.")
