.. _pySpline_building:

Building
--------

For speed purposes, :ref:`pySpline` uses a small compiled Fortran
library for doing the time consuming computational operators. It is
therefore necessary to build this library before using
:ref:`pySpline`.

To see a list of architectures that :ref:`pySpline` has been known to
compile on run::
   
   make

from the root directory. 

The easiest approach to to try the cloest one to your system and
attempt a build using (for example)::

   make LINUX_INTEL

If everything was successful, the following lines will be printed to
the screen (near the end)::

   Testing if module pyspline can be imported...
   Module libspline was successfully imported.

If you don't see this, it will be necessary to configure the build
manually. To configure manually, first copy a default configuration
file from the defaults folder like this (run this in the root
direcotry)::
  
   cp config/defaults/config.LINUX_INTEL.mk config

Now open ``config/config.LINUX_INTEL.mk`` which should look like::

  # Config File for LINUX and INTEL Compiler
  AR       = ar
  AR_FLAGS = -rvs
  RM       = /bin/rm -rf

  # Fortran compiler and flags
  FF90        = ifort
  FF90_FLAGS  = -r8 -O2 -fPIC

  # C compiler and flags
  CC       = gcc
  CC_FLAGS   = -O2 -fPIC

  # Define potentially different python, python-config and f2py executables:
  PYTHON = python
  PYTHON-CONFIG = python-config
  F2PY = f2py

  # Define additional flags for linking
  LINKER_FLAGS = -nofor_main -lifport
  SO_LINKER_FLAGS =-fPIC -shared

Modify these parameters are required and attempt the build again. If
you have sucessfully compiled the code on a new system, please contact
the developpers such that a new default configuration file can be
added.


