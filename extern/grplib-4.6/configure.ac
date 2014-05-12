#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(libgrp, 4.6, cxchelp@head.cfa.harvard.edu)

AM_INIT_AUTOMAKE(libgrp, 4.6)

AC_CONFIG_MACRO_DIR([m4])

LT_INIT

AM_PATH_PYTHON([2.7])
AX_PYTHON_DEVEL([>'2.7' <'3.0'])
AX_PYTHON_NUMPY

export GROUP_CFLAGS=$GROUP_CFLAGS
echo $GROUP_CFLAGS

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])

AC_ARG_ENABLE([stuberrorlib],
[ --enable-stuberrorlib Build without liberr support  ],
[ case "${enableval}" in
  yes) ciao=false ;;
  no) ciao=true ;;
  *) AC_MSG_ERROR([bad value ${enableval} for --enable-suberrorlib]) ;;
  esac], [ciao=true])
AM_CONDITIONAL([CIAO], [test x$ciao = xtrue])

AC_ARG_VAR(GROUP_CFLAGS, [Additional flags to be pass to the grouping library])

AC_OUTPUT(Makefile src/Makefile doc/Makefile python/Makefile)

