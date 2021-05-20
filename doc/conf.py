# Standard Python modules
import os
import sys

# External modules
from sphinx_mdolab_theme.config import *

# -- Path setup --------------------------------------------------------------
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.


sys.path.insert(0, os.path.abspath("../"))

# -- Project information -----------------------------------------------------
project = "pySpline"

# mock import for autodoc
autodoc_mock_imports = ["numpy", "scipy", "pyspline.libspline"]
