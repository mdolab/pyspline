# Config File for LINUX and INTEL Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

# Fortran compiler and flags
FF90        = ifort
FF90_FLAGS  = -r8 -O2 -fPIC

# C compiler and flags
CC       = icc
CC_FLAGS   = -O2 -fPIC -std=c99

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS =-fPIC -shared
