MK_TOP = ../../..
include $(MK_TOP)/Makefile.master
include $(MK_TOP)/include/Makefile.scidev

PYMODULE 	= stk

PY_MODULE_FILES	= $(PYMODULE).so

INSTALL_PYMODULE=$(PY_MODULE_FILES:%=$(PYTHON_DIR)/%)


LOCAL_INC	= \
	-I. -I$(PYTHON)/include/python$(PYTHONVER)

LOCAL_LIB	= $(STKLIB) 
MAKECLEAN	+= $(PYMODULE).so pystk.o

include $(MK_TOP)/Makefile.all

# -------------------------------------

$(PYMODULE).so:  pystk.o
	$(CCC) $(BUNDLEFLAG) $(DYN_FLAG) -o $@ pystk.o \
	$(ASCDMLIBS) $(LFLAGS) $(LOCAL_LIB)

announce1:
	@echo "   /-------------------------------------------------------\ "
	@echo "   |           Building $(PYMODULE) python module          | "
	@echo "   \-------------------------------------------------------/ " 
