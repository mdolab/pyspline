# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

FF90 = gfortran
FFLAGS    = -fdefault-real-8 -O2 -fPIC

CC       = gcc
CFLAGS   = -O2 -fPIC

COMPILER_NAME = gfortran
EXTRA_LIBS = -lgfortran

# Combine Flags
# ------------------------------------

FF90_FLAGS = $(FFLAGS) $(USE_CGNS) $(USE_TECIO) $(CGNS_INCLUDE) 
CC_FLAGS   = $(CFLAGS)
USER_LIBS = $(CGNS_LIB) $(TEC_LIB) $(LIBSTDCpp) $(EXTRA_LIBS)
