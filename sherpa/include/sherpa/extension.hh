// 
//  Copyright (C) 2007, 2016, 2017  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_extension_hh__
#define __sherpa_extension_hh__

#include <sherpa/array.hh>

typedef unsigned int SherpaUInt;
typedef sherpa::Array< unsigned int, NPY_UINT > SherpaUIntArray;

typedef double SherpaFloat;
typedef sherpa::Array< double, NPY_DOUBLE > DoubleArray;
typedef DoubleArray SherpaFloatArray;

typedef sherpa::Array< int, NPY_INT > IntArray;

typedef int (*converter)( PyObject*, void* );

#define CONVERTME(arg) ((converter) sherpa::convert_to_contig_array<arg>)

#if PY_MAJOR_VERSION >= 3
#define PY3
#endif

#ifdef PY3

#define SHERPAMOD(name, fctlist) \
static struct PyModuleDef module##name = {\
PyModuleDef_HEAD_INIT, \
#name, \
NULL, \
-1, \
fctlist \
}; \
\
PyMODINIT_FUNC PyInit_##name(void) { \
  import_array(); \
  return PyModule_Create(&module##name); \
}

#else

#define SHERPAMOD(name, fctlist) \
PyMODINIT_FUNC init##name(void);\
PyMODINIT_FUNC \
init##name(void) \
{ \
  import_array(); \
  Py_InitModule( (char*)#name, fctlist ); \
}

#endif


#define FCTSPEC(name, func) \
 { (char*)#name, (PyCFunction)func, METH_VARARGS, NULL }


#endif /* __sherpa_extension_hh__ */
