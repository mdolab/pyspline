# Config File for LINUX and LLVM (clang/flang) Compilers
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

# Fortran compiler and flags
# Note: flang only accepts -std=f2018 (no f2008), so no -std flag is used here.
FF90        = flang
FF90_FLAGS  = -fdefault-real-8 -O2 -fPIC

# C compiler and flags
CC       = clang
CC_FLAGS   = -O2 -fPIC -std=c99

# Define potentially different python, python-config and f2py executables:
PYTHON = python
PYTHON-CONFIG = python3-config # use python-config for python 2
F2PY = f2py

# Define additional flags for linking
LINKER_FLAGS =
SO_LINKER_FLAGS =-fPIC -shared
