#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])

AC_INIT(sherpa-dependencies, 4.7, olaurino@cfa.harvard.edu)
AM_INIT_AUTOMAKE(foreign)

LT_INIT

AC_ARG_ENABLE([fftw],
  AS_HELP_STRING([--enable-fftw], [Build fftw library shipped with Sherpa]))

AS_IF([test "x$enable_fftw" = "xyes"], [
  AC_CONFIG_SUBDIRS([fftw-3.3.3])
  AC_SUBST(fftw, fftw-3.3.3)
])

AC_ARG_ENABLE([region],
  AS_HELP_STRING([--enable-region], [Build region library shipped with Sherpa]))

AS_IF([test "x$enable_region" = "xyes"], [
  AC_CONFIG_SUBDIRS([region-4.6])
  AC_SUBST(region, region-4.6)
])

AC_CONFIG_SUBDIRS([grplib-4.6])

AC_OUTPUT([Makefile])
