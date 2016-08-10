lib_LTLIBRARIES = libstk.la

include_HEADERS = dmstack.h stack.h

libstk_la_SOURCES = logical.h stack.h stack.c

pyexec_LTLIBRARIES = stk.la
stk_la_SOURCES = pystk.c pystk.h logical.h stack.h stack.c
stk_la_LDFLAGS = $(PYTHON_LDFLAGS) -avoid-version -module
stk_la_CFLAGS = $(PYTHON_CPPFLAGS) -I$(NUMPY_INCLUDEDIR)

