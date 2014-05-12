
SUBDIRS = $(fftw) $(grplib) $(region) $(wcssubs)

WCS = $(wcssubs)

EXTRA_DIST = $(WCS)

#lib_LTLIBRARIES = libwcs.la

#include_HEADERS = $(WCS)/wcslib.h $(WCS)/wcs.h $(WCS)/fitshead.h

#libwcs_la_SOURCES = $(WCS)/cel.c $(WCS)/distort.c $(WCS)/fitsfile.c $(WCS)/hget.c $(WCS)/hput.c \
#	$(WCS)/iget.c $(WCS)/imhfile.c $(WCS)/imio.c $(WCS)/lin.c $(WCS)/platepos.c $(WCS)/poly.c \
#	$(WCS)/proj.c $(WCS)/sph.c $(WCS)/tnxpos.c $(WCS)/zpxpos.c $(WCS)/wcs.c $(WCS)/wcsinit.c \
#	$(WCS)/wcscon.c $(WCS)/wcslib.c $(WCS)/wcstrig.c $(WCS)/worldpos.c \
#	$(WCS)/dateutil.c $(WCS)/fileutil.c $(WCS)/dsspos.c \
#	$(WCS)/wcslib.h $(WCS)/wcs.h $(WCS)/fitshead.h $(WCS)/fitfile.h


#all-local:
#	cd wcssubs-3.8.7 && make && cd -
#
#install-exec-local : wcssubs-3.8.7/libwcs.a
#	cp wcssubs-3.8.7/libwcs.a $(libdir); \
#	cp wcssubs-3.8.7/wcs.h $(includedir); \
#	cp wcssubs-3.8.7/wcslib.h $(includedir); \
#	cp wcssubs-3.8.7/fitshead.h $(includedir);
#
#uninstall-local:
#	rm -f $(libdir)/libwcs.a $(includedir)/wcs.h
#
#clean-local:
#	rm -f wcssubs-3.8.7/*.a wcssubs-3.8.7/*.o
