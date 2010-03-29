#*********************************************************
# Makefile for cmlib Tensor B-spline
# G. Kenway May 21, 2009
#*********************************************************

INCLUDE_COMMON = ../common.mk
include ${INCLUDE_COMMON}
RULES = ../config/rules.mk
include ${RULES}

SPLINE_LIB = libspline.a
SPLINE_LIB_SHARED = libspline.so
OBJ =	compute_curve.o\
	compute_surface.o\
	compute_volume.o\
	eval_curve.o\
	eval_surface.o\
	eval_volume.o\
	projections.o\
	bvalu.o\
	intrv.o\
	bspvn.o\
	bspvd.o\
	basis.o\
	knots_lms.o\
	knots_interp.o\
	b2val.o\
	b3val.o\
	insertKnot.o\
	updatecoefficients.o\
	updatecoefficients_c.o\
	getBasisPt.o\
	cgns_functions.o\
	tecplot_functions.o\
	parameterizations.o\
	machcon/i1mach.o\
	machcon/r1mach.o\
	xerror/fdump.o\
	xerror/j4save.o\
	xerror/xerabt.o\
	xerror/xerctl.o\
	xerror/xerprt.o\
	xerror/xerror.o\
	xerror/xerrwv.o\
	xerror/xersav.o\
	xerror/xgetua.o\

###================================================
###  Default compilers and flags

default: libspline

clean:
	$(RM) *.o xerror/*.o machcon/*.o *.mod libspline.a libspline.so

libspline: $(OBJ) 
	$(AR) $(AR_FLAGS) $(SPLINE_LIB) $(OBJ)
	$(FC)  -shared $(OBJ) -o $(SPLINE_LIB_SHARED) 
	f2py  --fcompiler=${COMPILER_NAME} --f90flags='${FFLAGS}' -c -m pyspline pyspline.pyf libspline.a ${USER_LIBS}
	mv pyspline.so ../python	