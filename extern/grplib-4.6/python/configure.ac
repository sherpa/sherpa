#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(group, 4.6, olaurino@cfa.harvard.edu)

AM_INIT_AUTOMAKE(group, 4.6)

AC_CONFIG_MACRO_DIR([m4])

LT_INIT([disable-static])

AM_PATH_PYTHON([2.7])
AX_PYTHON_DEVEL([>'2.7' <'3.0'])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_CHECK_FUNCS([pow sqrt])


AC_OUTPUT(Makefile )
