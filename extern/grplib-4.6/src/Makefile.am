
lib_LTLIBRARIES = libgrp.la

include_HEADERS = grplib.h

libgrp_la_SOURCES = grplib.c grp_priv.c stub_dserr.c stub_dserr.h

if CIAO

announce1:
	@echo "   /----------------------------------------------------------\  "
	@echo "   |                    Building grp library                  |  "
	@echo "   \----------------------------------------------------------/  "
else
libgrp_la_CFLAGS = -DSTUB_ERROR_LIB=1
endif


EXTRA_DIST = *.h

