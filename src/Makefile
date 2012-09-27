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
MODS = precision.o\

OBJS =  bvalu.o\
	intrv.o\
	bspvn.o\
	bspvd.o\
	b2val.o\
	b3val.o\
	cbvalu.o\
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
	basis.o\
	knots.o\
	insertKnot.o\
	getBasisPt.o\
	parameterizations.o\
	findSpan.o\
	compute_curve.o\
	compute_surface.o\
	compute_volume.o\
	eval_curve.o\
	eval_surface.o\
	eval_volume.o\
	projections.o\
	tfi2d.o\

default:
	$(MAKE) $(MODS)
	cp -f *.mod ../mod
	$(MAKE) $(OBJS)

clean:
	@echo "        Making clean in src... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)

