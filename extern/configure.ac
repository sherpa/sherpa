#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

AC_INIT(sherpa-dependencies, 4.7, cxchelp@head.cfa.harvard.edu)
AM_INIT_AUTOMAKE([foreign])

AC_CONFIG_MACRO_DIR([m4])

AM_MAINTAINER_MODE([enable])

LT_INIT

AC_ARG_ENABLE([fftw],
  AS_HELP_STRING([--enable-fftw], [Build fftw library shipped with Sherpa]))

AS_IF([test "x$enable_fftw" = "xyes"], [
  AC_CONFIG_SUBDIRS([fftw-3.3.4])
  AC_SUBST(fftw, fftw-3.3.4)
  touch fftw-3.3.4/doc/fftw3.info
])

AC_ARG_ENABLE([region],
  AS_HELP_STRING([--enable-region], [Build region library shipped with Sherpa]))

AS_IF([test "x$enable_region" = "xyes"], [
  AC_CONFIG_SUBDIRS([region-4.6])
  AC_SUBST(region, region-4.6)
])

AC_ARG_ENABLE([group],
  AS_HELP_STRING([--enable-group], [Build group module shipped with Sherpa]))

AS_IF([test "x$enable_group" = "xyes"], [
  AC_CONFIG_SUBDIRS([grplib-4.6])
  AC_SUBST(grplib, grplib-4.6)
])

AC_ARG_ENABLE([wcssubs],
  AS_HELP_STRING([--enable-wcssubs], [Build wcssubs library shipped with Sherpa]))

AS_IF([test "x$enable_wcs" = "xyes"], [
  AC_CONFIG_SUBDIRS([wcssubs])
  AC_SUBST(wcssubs, wcssubs)
])

AC_ARG_ENABLE([stk],
  AS_HELP_STRING([--enable-stk], [Build stack library shipped with Sherpa]))

AS_IF([test "x$enable_stk" = "xyes"], [
  AC_CONFIG_SUBDIRS([stklib-4.6])
  AC_SUBST(stklib, stklib-4.6)
])


AC_OUTPUT([Makefile])
