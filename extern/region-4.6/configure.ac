#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

AC_CONFIG_MACRO_DIR([m4])

AC_INIT([region], [4.6], [cxchelp@head.cfa.harvard.edu])
AC_CONFIG_SRCDIR([src/region_create.c])
AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC

AC_PROG_YACC
if test x"$YACC" = "xyacc"; then
  AC_CHECK_PROG([YACC_EXISTS], [yacc], [yes], [no])
  if test x"$YACC_EXISTS" != xyes; then
    AC_MSG_ERROR([[bison/byacc/yacc not found.
          Please install bison]])
  fi
fi

AC_PROG_LEX
if test x"$LEX" = "x:"; then
  AC_MSG_ERROR([Lexer not found. Please install flex.])
fi

AC_PROG_SED

LT_INIT

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lregion':
AC_CHECK_LIB([region], [main])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([sqrt strcasecmp strchr strpbrk])

AC_OUTPUT(Makefile src/Makefile)
