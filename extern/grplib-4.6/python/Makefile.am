
pyexec_LTLIBRARIES = group.la
group_la_SOURCES = pygrplib.c pygrplib.h ../src/grplib.c ../src/grp_priv.c ../src/stub_dserr.c ../src/stub_dserr.h
group_la_LDFLAGS = -avoid-version -module

if CIAO
group_la_CFLAGS = -I../src $(PYTHON_CPPFLAGS) $(GROUP_CFLAGS)

else
group_la_CFLAGS = -DSTUB_ERROR_LIB=1 -I../src $(PYTHON_CPPFLAGS)

endif

