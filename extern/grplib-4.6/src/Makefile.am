
lib_LTLIBRARIES = libgrp.la

include_HEADERS = grplib.h

libgrp_la_SOURCES = grplib.c grp_priv.c stub_dserr.c stub_dserr.h


pyexec_LTLIBRARIES = group.la
group_la_SOURCES = pygrplib.c pygrplib.h grplib.c grp_priv.c stub_dserr.c stub_dserr.h
group_la_LDFLAGS = -avoid-version -module

if CIAO
group_la_CFLAGS = $(PYTHON_CPPFLAGS) $(GROUP_CFLAGS)

else
group_la_CFLAGS = -DSTUB_ERROR_LIB=1 $(PYTHON_CPPFLAGS) $(GROUP_CFLAGS)

endif


EXTRA_DIST = \
dserror_database.h  dserror_global.h      dserror_ptinstrument.h    dserror_structs.h    grp_priv.h  stub_dserr.h\
dserror_funcs.h     dserror_ptanalysis.h  dserror_ptnonsi.h         dserror_trace_fct.h  liberr.h\
dserror_general.h   dserror_ptdmtools.h   dserror_signal_handler.h  grplib.h             pygrplib.h

