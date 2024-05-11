# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

ICC_EXISTS := $(shell command -v icc)
ifdef ICC_EXISTS
  # icc only exists on older Intel versions
  # Assume that we want to use the old compilers
  FF90 = ifort
  CC = icc
else
  # Use the new compilers
  FF90 = ifx
  CC = icx
endif

# Compiler and flags
FF90_FLAGS  = -r8 -O2 -fPIC -stand f08
CC_FLAGS   = -O2 -fPIC -std=c99

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS =-fPIC -shared
