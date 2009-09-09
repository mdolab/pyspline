#*********************************************************
# Makefile for cmlib Tensor B-spline
# G. Kenway May 21, 2009
#*********************************************************

RM  	        = /bin/rm -rf

SPLINE_LIB = libspline_cs.a
SPLINE_LIB_SHARED = libspline_cs.so
OBJ = 	complexify.o \
	bspline/c_bfqad.o\
	bspline/c_bint4.o\
	bspline/c_bintk.o\
	bspline/c_bnfac.o\
	bspline/c_bnslv.o\
	bspline/c_bsgq8.o\
	bspline/c_bspdr.o\
	bspline/c_bspev.o\
	bspline/c_bsppp.o\
	bspline/c_bspvd.o\
	bspline/c_bspvn.o\
	bspline/c_bsqad.o\
	bspline/c_bvalu.o\
	bspline/c_intrv.o\
	bspline/c_pfqad.o\
	bspline/c_ppgq8.o\
	bspline/c_ppqad.o\
	bspline/c_ppval.o\
	machcon/c_i1mach.o\
	machcon/c_r1mach.o\
	tensbs/c_b2ink.o\
	tensbs/c_b2val.o\
	tensbs/c_b3ink.o\
	tensbs/c_b3val.o\
	tensbs/c_bknot.o\
	tensbs/c_btpcf.o\
	xerror/c_fdump.o\
	xerror/c_j4save.o\
	xerror/c_xerabt.o\
	xerror/c_xerctl.o\
	xerror/c_xerprt.o\
	xerror/c_xerror.o\
	xerror/c_xerrwv.o\
	xerror/c_xersav.o\
	xerror/c_xgetua.o\
	misc/c_bvaluV.o\
	misc/c_b2valM.o\
	misc/c_b2valV.o\
	misc/updateSurfacePoints.o\
	misc/getComplexCoef.o\


###================================================
###  Default compilers and flags

### Intel Fortran Compiler
FC = ifort

FFLAGS = -O2 -r8
AR       = ar
AR_FLAGS = -rvs

default: libspline_cs

clean:
	@echo "Removing .o files in ./"
	$(RM) *.o
	@echo "Removing .mod files in ./"
	$(RM) *.mod
	@echo "Removing .o files in bspline"
	$(RM) bspline/*.o
	@echo "Removing .o files in machcon"
	$(RM) machcon/*.o
	@echo "Removing .o files in tensbs"
	$(RM) tensbs/*.o
	@echo "Removing .o files in xerror"
	$(RM) xerror/*.o
	@echo "Removing .o files in misc"
	$(RM) misc/*.o
	@echo "Removing library"
	$(RM) libspline_cs.a
	$(RM) libspline_cs.so
libspline_cs: $(OBJ) 
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
