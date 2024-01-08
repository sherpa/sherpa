
lib_LTLIBRARIES = libgrp.la

include_HEADERS = grplib.h

libgrp_la_SOURCES = grplib.c grp_priv.c stub_dserr.c stub_dserr.h liberr.h grp_priv.h

if CIAO
libgrp_la_CFLAGS = $(GROUP_CFLAGS)
libgrp_la_LIBADD = -lerr
else
libgrp_la_CFLAGS = -DSTUB_ERROR_LIB=1 $(GROUP_CFLAGS)
libgrp_la_LDFLAGS = -static
endif

