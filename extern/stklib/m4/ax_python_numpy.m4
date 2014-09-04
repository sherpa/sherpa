# ===========================================================================
#
# ===========================================================================
#
# SYNOPSIS
#
#   AX_PYTHON_NUMPY([ACTION_IF_FOUND], [ACTION_IF_NOT_FOUND])
#
# COPYING
#
#   Copyright (c) 2009 David Grundberg
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_PYTHON_NUMPY],[
	NUMPY_INCLUDEDIR=

	AC_MSG_CHECKING([for python])
	AS_IF([test -z "$PYTHON"], [
		AC_MSG_RESULT([unknown])
	],[
		AC_MSG_RESULT([$PYTHON])
	])

	AC_MSG_CHECKING([for numpy includedir])
	AS_IF([test -z "$PYTHON"], [
		AC_MSG_RESULT([no (python unknown)])
	],[
		NUMPY_INCLUDEDIR=`$PYTHON -c '
try:
	from numpy import get_include
	print get_include()
except:
	pass
'`
		AC_SUBST([GROUP_CFLAGS], ["$GROUP_CFLAGS -I$NUMPY_INCLUDEDIR"])

		AC_MSG_RESULT([$NUMPY_INCLUDEDIR])
	])

dnl	AS_IF([test -n "$NUMPY_INCLUDEDIR"], [
dnl
dnl		AC_CHECKING([whether linking to numpy library works], [ax_python_numpy_cv_check],
dnl		[
dnl			ax_python_numpy_cv_check=no
dnl
dnl			AC_LANG_PUSH([C++])
dnl
dnl			ax_python_numpy_cppflags="$CPPFLAGS"
dnl			ax_python_numpy_ldflags="$LDFLAGS"
dnl			CPPFLAGS="$CPPFLAGS -I$NUMPY_INCLUDEDIR $PYTHON_CPPFLAGS"
dnl			LDFLAGS="$LDFLAGS $PYTHON_LDFLAGS"
dnl			
dnl			AC_LANG_ASSERT(C++)
dnl			AC_LINK_IFELSE(
dnl			AC_LANG_PROGRAM(
dnl				[[
dnl#define PY_ARRAY_UNIQUE_SYMBOL my_array_symbol
dnl#include <Python.h>
dnl#include <numpy/oldnumeric.h>
dnl#include <numpy/old_defines.h>
dnl]],
dnl				[[ &PyArray_FromDims; ]]),
dnl				[ax_python_numpy_cv_check=yes],
dnl				[ax_python_numpy_cv_check=no])
dnl			CPPFLAGS="$ax_python_numpy_cppflags"
dnl			LDFLAGS="$ax_python_numpy_ldflags"
dnl			
dnl			AC_LANG_POP([C++])
dnl		])
dnl	])

dnl	AS_IF([test "x$ax_python_numpy_cv_check" != "xyes"], [
dnl		NUMPY_INCLUDEDIR=
dnl
dnl		AC_MSG_WARN([[
dnl========================================================================
dnl Can not link with Numpy.
dnl
dnl Make sure the Numpy development package is installed.
dnl========================================================================]])
dnl	])
dnl
	AC_SUBST([NUMPY_INCLUDEDIR])

	# Execute ACTION_IF_FOUND or ACTION_IF_NOT_FOUND
	if test "x$ax_python_numpy_cv_check" == "xyes" ; then
		m4_ifvaln([$1],[$1],[:])dnl
		m4_ifvaln([$2],[else $2])dnl
	fi

])
