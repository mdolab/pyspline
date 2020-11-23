#! /usr/bin/env python
import sys

print("Testing if module pyspline can be imported...")
try:
    import libspline  # noqa: F401
except ImportError:
    print("Error importing libspline.so")
    sys.exit(1)
# end try

print("Module libspline was successfully imported.")
