#! /usr/bin/env python
from __future__ import print_function
import sys

print('Testing if module pyspline can be imported...')
try:
    import libspline
except:
    print('Error importing libspline.so')
    sys.exit(1)
# end try

print('Module libspline was successfully imported.')
