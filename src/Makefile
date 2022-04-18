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
GIT_VERSION   = $(shell git describe --dirty --always --tags)

# Name of files to be compiled:
MODS = precision.o\
       adtProjections.o

OBJS =  evaluations.o\
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
	version.o

.PHONY: version

default: version
	$(MAKE) $(MODS)
	cp -f *.mod ../mod
	$(MAKE) $(OBJS)

version:
	@echo "Making auto-generated version file"
	@echo GIT_VERSION is $(GIT_VERSION)
	@echo "subroutine version(ver)" > version.f90
	@echo "  ! Automatically generated version function" >> version.f90
	@echo "  !" >> version.f90
	@echo "  ! Desricption of Arguments" >> version.f90
	@echo "  ! Output" >> version.f90
	@echo "  ! ver    - string, contains git tag and commit hash" >> version.f90
	@echo "  use precision" >> version.f90
	@echo "  implicit none" >> version.f90
	@echo "  character(len=256), intent(out) :: ver" >> version.f90
	@echo "  " >> version.f90
	@echo "  ver = '$(GIT_VERSION)'" >> version.f90
	@echo "end subroutine version" >> version.f90

clean:
	@echo "        Making clean in src... "
	rm -f $(MAKE_CLEAN_ARGUMENTS)
