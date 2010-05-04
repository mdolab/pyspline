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
SPLINE_OBJ =	compute_curve.o\
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
		xerror/xgetua.o

EXTRAS_LIB =    libspline_extras.a
EXTRAS_OBJ = 	cgns_functions.o\
		tecplot_functions.o

###================================================
###  Default compilers and flags

default: libspline_extras libspline

clean:
	$(RM) *.o xerror/*.o machcon/*.o *.mod libspline.a libspline.so libspline_extras.a

libspline_extras: $(EXTRAS_OBJ)	
	$(AR) $(AR_FLAGS) $(EXTRAS_LIB) $(EXTRAS_OBJ)

libspline: $(SPLINE_OBJ) 
	$(AR) $(AR_FLAGS) $(SPLINE_LIB) $(SPLINE_OBJ)
	$(FC)  -shared $(SPLINE_OBJ) -o $(SPLINE_LIB_SHARED) 
	f2py  --fcompiler=${COMPILER_NAME} --f90flags='${FFLAGS}' -c -m pyspline pyspline.pyf libspline.a libspline_extras.a ${USER_LIBS}
	mv pyspline.so ../python	

