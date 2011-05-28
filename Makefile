# Makefile for pySpline

SUBDIR_SRC    = src/

#      ******************************************************************
#      *                                                                *
#      * General targets.                                               *
#      *                                                                *
#      ******************************************************************

default:
	@echo "Usage: make <arch>"
	@echo "Supported architectures: LINUX_INTEL"
	@echo "                         LINUX_GFORTRAN"
	@echo "                         LINUX_INTEL_SCINET"
	@echo "                         BASALT"

all:	 default

clean:
	@echo " Making clean ... "

	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo; \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make $@) || exit 1; \
		done
	rm -f *~ config.mk;
	rm -f lib/lib* mod/* obj/*

#      ******************************************************************
#      *                                                                *
#      * The actual make. This is not a direct target, but is called    *
#      * from the architectures.                                        *
#      *                                                                *
#      ******************************************************************

module:
	@for subdir in $(SUBDIR_SRC) ; \
		do \
			echo "making $@ in $$subdir"; \
			echo; \
			(cd $$subdir && make) || exit 1; \
		done
	(cd lib && make)


#      ******************************************************************
#      *                                                                *
#      * Platform specific targets.                                     *
#      *                                                                *
#      ******************************************************************

LINUX_INTEL:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL.mk" ]; then cp "config/defaults/config.LINUX_INTEL.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL.mk config.mk
	make module
	(cd src/f2py && make f2py)

LINUX_GFORTRAN:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_GFORTRAN.mk" ]; then cp "config/defaults/config.LINUX_GFORTRAN.mk" ./config; fi
	ln -sf config/config.LINUX_GFORTRAN.mk config.mk
	make module
	(cd src/f2py && make f2py)

LINUX_INTEL_SCINET:
	mkdir -p obj
	if [ ! -f "config/config.LINUX_INTEL_SCINET.mk" ]; then cp "config/defaults/config.LINUX_INTEL_SCINET.mk" ./config; fi
	ln -sf config/config.LINUX_INTEL_SCINET.mk config.mk
	make module
	(cd src/f2py && make f2py)

BASALT:
	mkdir -p obj
	if [ ! -f "config/config.BASALT.mk" ]; then cp "config/defaults/config.BASALT.mk" ./config; fi
	ln -sf config/config.BASALT.mk config.mk
	make module
	(cd src/f2py && make f2py)



