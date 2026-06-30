//
//  Copyright (C) 2026
//  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//

#ifndef __sherpa_python_hh__
#define __sherpa_python_hh__

// Use the same setup for extension modules

// The current (mid 2026) free-threaded build does not support the
// limited API, so only set if not using such a build (although the
// current Sherpa build does not officially support the free-threaded
// build).
//
// This code is from
// https://github.com/python/cpython/blob/main/Modules/xxlimited_35.c
//
#include <pyconfig.h>   // Py_GIL_DISABLED
#ifndef Py_GIL_DISABLED
#  define Py_LIMITED_API 0x030B0000
#endif

#include <Python.h>

#endif /* __sherpa_python_hh__ */
