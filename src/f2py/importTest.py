#! /usr/bin/env python
# Standard Python modules
import importlib.util
import os
import sys

print("Testing if module pyspline can be imported...")
try:
    # Checks for the libspline module in current working directory
    cwd = os.getcwd()
    spec = importlib.util.spec_from_file_location("libspline", os.path.join(cwd, "libspline.so"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
except ImportError as e:
    print("Error importing libspline.so")
    print(e)
    sys.exit(1)
# end try

print("Module libspline was successfully imported.")
