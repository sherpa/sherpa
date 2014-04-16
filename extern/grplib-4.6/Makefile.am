AUTOMAKE_OPTIONS = foreign

if CIAO
SUBDIRS = src python doc
else
SUBDIRS = src
endif

DIST_SUBDIRS = $(SUBDIRS)

forceclean:
	-rm -f Makefile src/Makefile python/Makefile doc/Makefile

