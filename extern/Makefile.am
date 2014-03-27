SUBDIRS = fftw-3.3.3 grplib-4.6 region-4.6

all-local:
	cd wcssubs-3.8.7 && make && cd -

install-exec-local : wcssubs-3.8.7/libwcs.a
	cp wcssubs-3.8.7/libwcs.a $(libdir); \
	cp wcssubs-3.8.7/wcs.h $(includedir); \
	cp wcssubs-3.8.7/wcslib.h $(includedir); \
	cp wcssubs-3.8.7/fitshead.h $(includedir);

uninstall-local:
	rm -f $(libdir)/libwcs.a $(includedir)/wcs.h

clean-local:
	rm -f wcssubs-3.8.7/*.a wcssubs-3.8.7/*.o
