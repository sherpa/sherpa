//
//  Copyright (C) 2007, 2016-2017, 2020, 2022-2023, 2025-2026
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

#ifndef __sherpa_extension_hh__
#define __sherpa_extension_hh__

#include <complex>

#include <sherpa/array.hh>

typedef unsigned int SherpaUInt;
typedef sherpa::Array< unsigned int, NPY_UINT > SherpaUIntArray;

typedef double SherpaFloat;
typedef sherpa::Array< double, NPY_DOUBLE > DoubleArray;
typedef DoubleArray SherpaFloatArray;

typedef sherpa::Array< int, NPY_INT > IntArray;
typedef sherpa::Array< std::complex<double>, NPY_COMPLEX128 > ComplexArray;

typedef int (*converter)( PyObject*, void* );

#define CONVERTME(arg) ((converter) sherpa::convert_to_contig_array<arg>)

// Define the macros for the SHERPAMOD/SHERPAMODDOC macros, which are
// used for extensions with no "state".
//
#define _SHERPAMOD_EXEC(name) \
  static int exec##name(PyObject *Py_UNUSED(m)) { \
    if (PyArray_ImportNumPyAPI() < 0) { return -1; } \
    return 0; \
  }

// The assumption is that with no setup these models have no "state"
// and so can be used with no GIL when building as a free-threaded
// build, which - for Python 3.14 and earlier, means no limited-API.
//
#if defined(Py_GIL_DISABLED) && !defined(Py_LIMITED_API)
#define _SHERPAMOD_SLOTS_GIL {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#else
#define _SHERPAMOD_SLOTS_GIL
#endif

#define _SHERPAMOD_SLOTS(name) \
  static PyModuleDef_Slot slots##name[] = { \
    {Py_mod_exec, (void *) exec##name},    \
    _SHERPAMOD_SLOTS_GIL \
    {0, NULL} \
  };

#define _SHERPAMOD_DEF(name, fctlist, doc)  \
  static struct PyModuleDef module##name = { \
    PyModuleDef_HEAD_INIT, \
    #name, \
    doc, \
    0, \
    fctlist, \
    slots##name, \
    NULL, \
    NULL, \
    NULL \
  };

#define _SHERPAMOD_FUNC(name) \
  PyMODINIT_FUNC PyInit_##name(void) { \
    return PyModuleDef_Init(&module##name); \
  }

// Create the initialization for the module called name, with the list
// of functions in fctlist, and the module documentation set to doc
// (which can be NULL). This is for modules with no state.
//
// The SHERPAMOD and SHERPAMODDOC defines are expected to be used
// rather than this one.
//
#define _SHERPAMOD(name, fctlist, doc) \
  _SHERPAMOD_EXEC(name) \
  _SHERPAMOD_SLOTS(name) \
  _SHERPAMOD_DEF(name, fctlist, doc) \
  _SHERPAMOD_FUNC(name)

// Create the initialization for the module called name, with the list
// of functions in fctlist, and no module documentation.
//
// fctlist should be an array of PyMethodDef elements, which can be
// constructed with FCTSPEC, FCTSPECDOC, and KWSPEC, or created
// directly.
//
#define SHERPAMOD(name, fctlist) \
  _SHERPAMOD(name, fctlist, NULL)

// Similar to SHERMAMOD, but adds the ability to set the module
// documentation (as a string).
//
#define SHERPAMODDOC(name, fctlist, doc) \
  _SHERPAMOD(name, fctlist, PyDoc_STR(doc))

// Create a function specification: name is the name for the
// function (accessible from Python) and func is the function
// to be called. The function will have no documentation.
//
#define FCTSPEC(name, func) \
 { (char*)#name, (PyCFunction)func, METH_VARARGS, NULL }

// Add a docstring to FCTSPEC.
//
#define FCTSPECDOC(name, func, doc) \
  { (char*)#name, (PyCFunction)func, METH_VARARGS, PyDoc_STR(doc) }

// Similar to FCTSPEC but allows the function to have keywords.
//
#define KWSPEC(name, func) \
  { (char*)#name, (PyCFunction)((PyCFunctionWithKeywords)func), \
      METH_VARARGS|METH_KEYWORDS, NULL }

#endif /* __sherpa_extension_hh__ */
