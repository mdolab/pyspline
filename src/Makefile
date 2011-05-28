#*********************************************************
# Makefile for pySpline
# G. Kenway May 21, 2009
#*********************************************************

MAIN_DIR = ../
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}
vpath %.o $(OBJDIR)

# Name of files to be compiled:
OBJS = bvalu.o\
	cbvalu.o\
	intrv.o\
	bspvn.o\
	bspvd.o\
	basis.o\
	knots_lms.o\
	knots_interp.o\
	b2val.o\
	b3val.o\
	cb2val.o\
	cb3val.o\
	insertKnot.o\
	getBasisPt.o\
	parameterizations.o\
	i1mach.o\
	r1mach.o\
	fdump.o\
	j4save.o\
	xerabt.o\
	xerctl.o\
	xerprt.o\
	xerror.o\
	xerrwv.o\
	xersav.o\
	xgetua.o\
	compute_curve.o\
	compute_surface.o\
	compute_volume.o\
	eval_curve.o\
	eval_surface.o\
	eval_volume.o\
	projections.o\
	quality.o\
	skew_hexa_b.o\
	scalar_triple_product_b.o\
	adBuffer.o\
	adStack.o\
	cgns_functions.o\
	tecplot_functions.o

default:
	$(MAKE) $(OBJS)

clean:
	@echo "        Making clean in src... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

