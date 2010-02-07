#*********************************************************
# Makefile for pySpline Library
# G. Kenway Feb 7, 2010
#*********************************************************

default: intel

intel:
	@echo "Linux - Intel"
	cp ./config/config.LINUX_INTEL.mk common.mk
	( cd src && $(MAKE)) || exit 1; 
	f2py  --fcompiler=intel --f90flags=-r8 -c -m pyspline src/pyspline.pyf src/libspline.a -lblas
	mv pyspline.so ./python
	rm common.mk
gfortran:
	@echo "Linux - Gfortran"
	cp ./config/config.LINUX_GFORTRAN.mk common.mk
	( cd src && $(MAKE)) || exit 1; 
	f2py  --fcompiler=gfortran --f90flags=-fdefault-real-8 -c -m pyspline src/pyspline.pyf src/libspline.a -lblas
	mv pyspline.so ./python
	rm common.mk
scinet:
	@echo "Scinet"
	cp ./config/config.SCINET.mk common.mk
	( cd src && $(MAKE)) || exit 1; 
	f2py  --fcompiler=intelem --f90flags=-r8 -c -m pyspline src/pyspline.pyf src/libspline.a -lblas
	mv pyspline.so ./python
	rm common.mk

clean:
	cp ./config/config.LINUX_INTEL.mk common.mk
	(cd src && $(MAKE) $@) || exit 1;
	rm common.mk

# Note we have to copy a dummy config file to common.mk so it doesn't complain
# about the file not existing