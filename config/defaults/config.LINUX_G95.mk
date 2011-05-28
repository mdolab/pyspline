# Config File for LINUX and GFORTRAN Compiler
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf
FC       = g95
FLAGS   = -O2 -r8
COMPILER_NAME = g95

CC       = gcc
CFLAGS   = -O2 -fPIC

# Optional Things

# CGNS Functionality
USE_CGNS = -DUSE_CGNS
CGNS_LIB = -lcgns

# Tecplot Binary IO
#USE_TECIO = -DUSE_TECIO 
#TEC_LIB   = ../tecio/tecio.a
#LIBSTDCpp = -lstdc++

# Combine Flags
FFLAGS = $(FLAGS) $(USE_CGNS) ${USE_TECIO} $(CGNS_FLAGS)
# Additional User Libs
USER_LIBS = $(CGNS_LIB) $(TEC_LIB) $(LIBSTDCpp)