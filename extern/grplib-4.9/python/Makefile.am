
pyexec_LTLIBRARIES = group.la
group_la_SOURCES = pygrplib.c pygrplib.h  
DEFAULTCFLAGS = -I$(top_srcdir)/src $(PYTHON_CPPFLAGS) $(GROUP_CFLAGS)
group_la_LDFLAGS = -avoid-version -module
group_la_LIBADD = ../src/libgrp.la

if CIAO
group_la_CFLAGS = $(DEFAULTCFLAGS)
else
group_la_CFLAGS = -DSTUB_ERROR_LIB=1 $(DEFAULTCFLAGS)
endif

