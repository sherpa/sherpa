
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
	reglexer.c

IOSRCS = region_write.c region_read.c region_type.c
IOINC_FILES = regionio.h

INC_FILES = cxcregion.h 
PKG_FILES = cxcregion.pc

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


test1: test.c $(LIB)
	  cc test.c -L. -lregion -lm -g -o test1

test2: test_2.c $(LIB)
	  cc test_2.c -L. -lregion -lm -g -o test2

test3: test3.c $(LIB)
	  cc test3.c -L. -lregion -lm -g -o test3


$(PKG_FILES): $(LIB_FILES)
	@echo
	\rm -f $@
	@echo Building $@
	@echo prefix=$(ASCDS_BLD_INSTALL)       > $@
	@echo exec_prefix=$$\{prefix}	       >> $@
	@echo libdir=$$\{exec_prefix}/lib      >> $@
	@echo includedir=$$\{prefix}/include   >> $@
	@echo				       >> $@     
	@echo Name: cxcregion                  >> $@             
	@echo Description: CXC Region Library  >> $@
	@echo Version:                         >> $@
	@echo Requires:                        >> $@
	@echo Libs: -L$$\{libdir} -lregion $(MATHLIBS)  >> $@
	@echo Cflags: -I$$\{includedir}        >> $@

announce1:
	@echo " /-------------------------------------"
	@echo " | Building region library"
	@echo " \-------------------------------------"
	@echo

