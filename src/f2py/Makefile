MAIN_DIR = ../..
COMMON_FILE = $(MAIN_DIR)/Common.mk
RULES_FILE  = $(MAIN_DIR)/Rules.mk
include ${COMMON_FILE}
include ${RULES_FILE}

f2py: 
	f2py --verbose --fcompiler=${COMPILER_NAME} --f90flags='${FF90_ALL_FLAGS}' -c -m pyspline pyspline.pyf $(LIBDIR)/libspline.a ${USER_LIBS}
	python importTest.py
	mv pyspline.so ../../python
