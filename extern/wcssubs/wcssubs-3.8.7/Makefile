OBJS =  wcsinit.o wcs.o wcscon.o fitsfile.o imhfile.o \
	hget.o hput.o iget.o imio.o worldpos.o platepos.o \
	tnxpos.o zpxpos.o dsspos.o poly.o \
        wcslib.o lin.o cel.o proj.o sph.o wcstrig.o dateutil.o distort.o

libwcs.a:	$(OBJS)
	ar rv $@ $?
	ranlib $@

cel.o:		wcslib.h
distort.o:	wcs.h fitshead.h wcslib.h
fitsfile.o:	fitsfile.h fitshead.h
hget.o:		fitshead.h
hput.o:		fitshead.h
iget.o:		fitshead.h
imhfile.o:	fitsfile.h fitshead.h
imio.o:		fitshead.h
lin.o:		wcslib.h
platepos.o:	wcs.h fitshead.h wcslib.h
poly.o:		wcslib.h
proj.o:		wcslib.h
sph.o:		wcslib.h
tnxpos.o:	wcs.h fitshead.h wcslib.h
zpxpos.o:	wcs.h fitshead.h wcslib.h
wcs.o:		wcs.h fitshead.h wcslib.h
wcsinit.o:	wcs.h fitshead.h wcslib.h
wcscon.o:	wcs.h fitshead.h wcslib.h
wcslib.o:	wcslib.h
wcstrig.o:	wcslib.h
worldpos.o:	wcs.h fitshead.h wcslib.h
dateutil.o:	fitsfile.h fitshead.h
fileutil.o:	fitsfile.h

clean:
	rm -f *.a *.o

