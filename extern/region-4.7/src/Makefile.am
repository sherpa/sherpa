
BUILT_SOURCES = regparser.c reglexer.c

lib_LTLIBRARIES = libregion.la

include_HEADERS = cxcregion.h

libregion_la_SOURCES = \
	region_create.c    \
	region_inside.c    \
	region_util.c      \
	region_extent.c    \
	region_attribute.c \
        regparser.c        \
	reglexer.c         \
        region_priv.h      \
        regpar.h           \
        regparser.h        \
        regparser.tab.h

#INC_FILES = cxcregion.h
#PKG_FILES = cxcregion.pc

regparser.c: regparser.y
	echo $(YACC)
	$(YACC) -bregparser -d $(YFLAGS) $<
	$(SED) 's/yy/regYY/g' regparser.tab.c > regparser.c
	$(SED) 's/yy/regYY/g' regparser.tab.h > regparser.h

resed:
	$(SED) 's/yy/regYY/g' regparser.tab.c > regparser.c
	$(SED) 's/yy/regYY/g' regparser.tab.h > regparser.h
	$(SED) 's/yy/regYY/g' reglexer.c.tmp >reglexer.c

reglexer.c: reglexer.l
	rm -f $@.tmp
	$(LEX) -t $< > $@.tmp
	- chmod +w $@.tmp 
	$(SED) 's/yy/regYY/g' $@.tmp > $@

clean-local:
	-rm -f regparser.c reglexer.c *.o regparser.h regparser.tab.c regparser.tab.h

EXTRA_DIST =  reglexer.c.tmp reglexer.l regparser.y 

