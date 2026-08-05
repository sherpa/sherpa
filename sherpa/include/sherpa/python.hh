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
		
// Use the same Python setup for extension modules.
//

// Only needed until Python 3.13 is the minimum-supported version:
// https://docs.python.org/3/extending/extending.html
//
#define PY_SSIZE_T_CLEAN

// For now there is no attempt to support
// - limited Python API
// - free-threaded build
//
#include <Python.h>

#endif /* __sherpa_python_hh__ */
