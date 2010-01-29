#*********************************************************
# Makefile for cmlib Tensor B-spline
# G. Kenway May 21, 2009
#*********************************************************

RM  	        = /bin/rm -rf

SPLINE_LIB = libspline.a
SPLINE_LIB_SHARED = libspline.so
OBJ =	lsqrDataModule.o\
	lsqrblasInterface.o\
	lsqrCheckModule.o\
	lsqrModule.o\
	lms_jacobian.o\
	compute_curve.o\
	eval_curve.o\
	eval_curve_deriv.o\
	eval_curve_deriv2.o\
	project_point_curve.o\
	compute_surface.o\
	eval_surface.o\
	eval_surface_deriv.o\
	eval_surface_deriv2.o\
	project_point_surface.o\
	bvalu.o\
	intrv.o\
	bspvn.o\
	knots.o\
	b2val.o\
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

### Intel Fortran Compiler
FC = ifort

FFLAGS = -O2 -r8 
#FFLAGS = -O2 -fdefault-real-8
AR       = ar
AR_FLAGS = -rvs

default: libspline

clean:
	@echo "Removing .o files"
	$(RM) *.o
	$(RM) xerror/*.o
	$(RM) machcon/*.o
	@echo "Removing .mod files"
	$(RM) *.mod
	@echo "Removing library"
	$(RM) libspline.a
	$(RM) libspline.so

libspline: $(OBJ) 
	$(AR) $(AR_FLAGS) $(SPLINE_LIB) $(OBJ)
	$(FC)  -shared $(OBJ) -o $(SPLINE_LIB_SHARED) 

.SUFFIXES: .o .f .f90

.f.o:	
	$(FC) $(FFLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f successfully ---"
	@echo	

.f90.o:	
	$(FC) $(FFLAGS) -c $< -o $*.o
	@echo
	@echo "        --- Compiled $*.f90 successfully ---"
	@echo	
