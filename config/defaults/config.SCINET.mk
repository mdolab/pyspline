# Config File for SCINET
AR       = ar
AR_FLAGS = -rvs
RM       = /bin/rm -rf
FC       = ifort
FLAGS   =  -r8 -O2 -fPIC
COMPILER_NAME = intelem
EXTRA_LIBS = -limf -lifcore

CC       = gcc
CFLAGS   = -O2 -fPIC

# Optional Things

# CGNS Functionality
# USE_CGNS = -DUSE_CGNS
# CGNS_LIB = -lcgns
# CGNS_INCLUDE = -I/usr/local/include/

# Tecplot Binary IO
# USE_TECIO = -DUSE_TECIO 
# TEC_LIB   = ../tecio/tecio.a
# LIBSTDCpp = -lstdc++

# Combine Flags
FFLAGS = $(FLAGS) $(USE_CGNS) ${USE_TECIO} $(CGNS_INCLUDE)

# Lib flags
USER_LIBS = $(CGNS_LIB) $(TEC_LIB) $(LIBSTDCpp) $(EXTRA_LIBS)