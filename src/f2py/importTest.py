#! /usr/bin/env python
from __future__ import print_function
import sys

print('Testing if module pyspline can be imported...')
try:
    import pyspline
except:
    print('Error importing pyspline.so.' )
    sys.exit(1)
# end try

print('Module pyspline was successfully imported.')
