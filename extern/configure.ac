#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT(sherpa-dependencies, 4.7, olaurino@cfa.harvard.edu)
AM_INIT_AUTOMAKE(foreign)
#AC_CONFIG_SRCDIR([wcssubs-3.8.7/cel.c])
#AC_CONFIG_HEADERS([config.h])

LT_INIT

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h sys/file.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor ftruncate gettimeofday memset pow sqrt strcasecmp strchr strncasecmp strrchr])

AC_CONFIG_SUBDIRS([fftw-3.3.3
                   grplib-4.6
                   region-4.6])
AC_OUTPUT([Makefile])
