pyexec_LTLIBRARIES = group.la
group_la_SOURCES = pygrplib.c pygrplib.h
group_la_CFLAGS = -I../src $(PYTHON_CPPFLAGS)
group_la_LDFLAGS = -avoid-version -module

#group_la_LIBADD = ../src/libgrp.la


#if CIAO

#MK_TOP = ../../..
#include $(MK_TOP)/Makefile.master
#include $(MK_TOP)/include/Makefile.scidev

#PYMODULE        = group

#PY_MODULE_FILES = $(PYMODULE).so

#INSTALL_PYMODULE= $(PY_MODULE_FILES:%=$(PYTHON_DIR)/%)


#LOCAL_INC       = \
#        -I. -I$(PYTHON)/include/python$(PYTHONVER) \
#        -I../src \
#        -I$(NUMPY)/lib/python$(PYTHONVER)/site-packages/numpy/core/include

#LOCAL_LIB       = #-L$(PYTHON)/lib -lpython$(PYTHONVER)
#MAKECLEAN       = $(PYMODULE).so pygrplib.o

#SLOCPAGE=Libdev
#SLOCEST=
#SLOCITEM=py$(PYMODULE)

#include $(MK_TOP)/Makefile.all

# -------------------------------------
#$(PYMODULE).so:  pygrplib.o
#          $(CCC) $(BUNDLEFLAG) $(DYN_FLAG) -o $@ pygrplib.o $(LFLAGS) $(LOCAL_LIB) $(GRPLIB) $(ERRLIBS)

#announce1:
#	@echo "   /-------------------------------------------------------\ "
#	@echo "   |           Building $(PYMODULE) python module  "
#	@echo "   \-------------------------------------------------------/ "

#else

#endif
