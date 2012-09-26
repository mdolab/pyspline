# Config File for LINUX and INTEL Compiler on SCINET
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf

FF90     = ifort
FFLAGS    =  -r8 -O2 -fPIC

CC       = gcc
CFLAGS   = -O2 -fPIC

COMPILER_NAME = intelem
EXTRA_LIBS = -limf -lifcore

# Optional Things

# Combine Flags
# ------------------------------------

FF90_FLAGS = $(FFLAGS) $(USE_CGNS) $(USE_TECIO) $(CGNS_INCLUDE) 
CC_FLAGS   = $(CFLAGS)
USER_LIBS = $(CGNS_LIB) $(TEC_LIB) $(LIBSTDCpp) $(EXTRA_LIBS)