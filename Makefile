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




# #*********************************************************
# # Makefile for pySpline Library
# # G. Kenway Feb 7, 2010
# #*********************************************************

# default: 
# 	@echo "Please select one of the following configurations"
# 	@echo "make intel       -> for linux pcs with intel compiler"
# 	@echo "make gfortran    -> for linux pcs with gfortran compiler"
# 	@echo "make scinet      -> for uoft scinet gpc cluster"
# 	@echo "make basalt      -> for sgi altix machine at utias"

# intel:
# 	@echo "Linux - Intel"
# 	-rm common.mk
# 	cp ./config/config.LINUX_INTEL.mk ./common.mk
# 	( cd src && make) || exit 1; 

# gfortran:
# 	@echo "Linux - Gfortran"
# 	-rm common.mk
# 	cp ./config/config.LINUX_GFORTRAN.mk ./common.mk
# 	( cd src && $(MAKE)) || exit 1; 
# 	-rm common.mk

# g95:
# 	@echo "Linux - G95"
# 	-rm common.mk
# 	cp ./config/config.LINUX_G95.mk ./common.mk
# 	( cd src && $(MAKE)) || exit 1; 
# 	-rm common.mk
# scinet:
# 	@echo "Scinet"
# 	-rm common.mk
# 	cp ./config/config.SCINET.mk ./common.mk
# 	( cd src && $(MAKE)) || exit 1; 
# 	-rm common.mk

# basalt:
# 	@echo "Basalt"
# 	-rm common.mk
# 	cp ./config/config.BASALT.mk common.mk
# 	( cd src && $(MAKE)) || exit 1; 
# 	-rm common.mk


# clean:
# 	-rm common.mk
# 	cp ./config/config.LINUX_INTEL.mk common.mk # Doesn't matter which 
# 	(cd src && $(MAKE) $@) || exit 1;
# 	-rm common.mk
# # Note we have to copy a dummy config file to common.mk so it doesn't complain
# # about the file not existing
