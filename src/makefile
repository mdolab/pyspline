#*********************************************************
# Makefile for cmlib Tensor B-spline
# G. Kenway May 21, 2009
#*********************************************************

RM  	        = /bin/rm -rf

TENSBS_LIB = libtensbs.a

OBJ =   bspline/bfqad.o\
	bspline/bint4.o\
	bspline/bintk.o\
	bspline/bnfac.o\
	bspline/bnslv.o\
	bspline/bsgq8.o\
	bspline/bspdr.o\
	bspline/bspev.o\
	bspline/bsppp.o\
	bspline/bspvd.o\
	bspline/bspvn.o\
	bspline/bsqad.o\
	bspline/bvalu.o\
	bspline/intrv.o\
	bspline/pfqad.o\
	bspline/ppgq8.o\
	bspline/ppqad.o\
	bspline/ppval.o\
	machcon/i1mach.o\
	machcon/r1mach.o\
	tensbs/b2ink.o\
	tensbs/b2val.o\
	tensbs/b3ink.o\
	tensbs/b3val.o\
	tensbs/bknot.o\
	tensbs/knots.o\
	tensbs/btpcf.o\
	xerror/fdump.o\
	xerror/j4save.o\
	xerror/xerabt.o\
	xerror/xerctl.o\
	xerror/xerprt.o\
	xerror/xerror.o\
	xerror/xerrwv.o\
	xerror/xersav.o\
	xerror/xgetua.o\
	misc/bvaluV.o\
	misc/b2valM.o\
	misc/b2valV.o\
	lmsfit/fit_foil.o\
	lmsfit/fit_surf.o\


###================================================
###  Default compilers and flags

### Intel Fortran Compiler
FC = ifort

FFLAGS = -O2 -r8
#FFLAGS = -O2 -fdefault-real-8
AR       = ar
AR_FLAGS = -rvs

default: libtensbs

clean:
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
	@echo "Removing .o files in fitpack"
	$(RM) fitpack/*.o
	@echo "Removing .o files in lmsfit"
	$(RM) lmsfit/*.o
	@echo "Removing library"
	$(RM) libtensbs.a

libtensbs: $(OBJ) 
	$(AR) $(AR_FLAGS) $(TENSBS_LIB) $(OBJ)


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



# O-files for the fitpack function NO LONGER USED
# 	fitpack/curfit.o\
# 	fitpack/fpbspl.o\
# 	fitpack/fpcurf.o\
# 	fitpack/fpgivs.o\
# 	fitpack/fprati.o\
# 	fitpack/fpback.o\
# 	fitpack/fpchec.o\
# 	fitpack/fpdisc.o\
# 	fitpack/fpknot.o\
# 	fitpack/fprota.o\
# 	fitpack/splev.o\
# 	fitpack/fpbacp.o\
# 	fitpack/fpchep.o\
# 	fitpack/fpgrpa.o\
# 	fitpack/fppasu.o\
# 	fitpack/fptrnp.o\
# 	fitpack/fptrpe.o\
# 	fitpack/parsur.o\
# 	fitpack/fpsuev.o\
# 	fitpack/surev.o\
# 	fitpack/fporde.o\
# 	fitpack/fprank.o\
# 	fitpack/fpsurf.o\
# 	fitpack/surfit.o\
# 	fitpack/bispev.o\
# 	fitpack/fpbisp.o\