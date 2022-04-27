//
//  Copyright (C) 2009, 2016, 2020, 2021, 2022
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


#include "sherpa/extension.hh"
#include <sstream>
#include <iostream>
#include <string>

extern "C" {
#include "cxcregion.h"

#ifdef USE_CXCDM_PARSER
#include "ascdm.h"
#endif
}

typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  regRegion *region;
} PyRegion;


// Declare symbols
static PyObject* region_union( PyRegion* self, PyObject* args, PyObject *kwargs );
static PyObject* region_subtract( PyRegion* self, PyObject* args, PyObject *kwargs );
static PyObject* region_invert( PyRegion* self, PyObject* args );


// The string handling depends on whether we have access to the
// CIAO data model, in which we case we can use dmRegParse and
// so support FITS region files, or if we do not have it, in
// which case we are restricted to ASCII region files via
// regReadAsciiRegion.
//
static regRegion* parse_string( std::string input, int fileflag ) {

  regRegion *reg = NULL;

#ifdef USE_CXCDM_PARSER

  if( fileflag )
    input = "region(" + input + ")";

  reg = dmRegParse( (char*)input.c_str() );

#else

  if( fileflag ) {
    reg = regReadAsciiRegion( (char*)input.c_str(), 0 ); // Verbosity set to 0
  } else {
    reg = regParse( (char*)input.c_str() );
  }
#endif

  return reg;
}


static PyObject* pyRegion_build(PyTypeObject *type, regRegion *reg) {
  PyRegion *ret = NULL;
  ret = (PyRegion*)type->tp_alloc( type, 0 );
  if( ret )
  ret->region = reg;

  return (PyObject*)ret;
}

// The constructor can be given a string and a flag, which if
// set indicates this is a file name. The file support is
// unnescessary as we can read the data in Python, so this
// API could just deal with strings (it lets the C code deal
// with the length of the region string rather than us).
//
static PyObject* pyRegion_new(PyTypeObject *type, PyObject *args,
			      PyObject *kwargs)
{
  regRegion *reg = NULL;
  char *objS = NULL;
  int fileflag = 0;

  static const char *kwlist[] = {"region", "fileflag", NULL};
  if ( !PyArg_ParseTupleAndKeywords( args, kwargs, "s|p",
                                     const_cast<char**>(kwlist),
				     &objS, &fileflag ) )
    return NULL;

  // The string is not NULL (otherwise the PyArg_ParseTupleAndKeywords
  // call would have failed) so we don't need to check here.
  //
  std::string name(objS);

  reg = parse_string( name, fileflag );
  if (!reg) {
    const char *fmt;
    if (fileflag) {
      fmt = "unable to read region from: %s";
    } else {
      fmt = "unable to parse region string: '%s'";
    }
    return PyErr_Format( PyExc_ValueError, fmt, objS );
  }

  return pyRegion_build( type, reg );
}

// New datatype dealloc method
static void pyRegion_dealloc(PyRegion* reg) {

  if( reg ) {
    if( reg->region )
      regFree( reg->region );

    PyObject_Del( reg );
  }
}

static PyObject* pyRegion_str(PyRegion* reg) {
  regRegion *region = reg->region;
  const char *ret = "";

  if ( region ) {
    ret = regAllocComposeRegion( region );
  }

  return Py_BuildValue("s", ret);
}


// Return the mask indicating whether the array points
// lie inside or outside the region.
//
static PyObject* region_mask( PyRegion* self, PyObject* args, PyObject *kwargs )
{

  DoubleArray xpos, ypos;

  static const char *kwlist[] = {"x0", "x1", NULL};
  if ( !PyArg_ParseTupleAndKeywords( args, kwargs, "O&O&",
                                     const_cast<char**>(kwlist),
				     CONVERTME(DoubleArray), &xpos,
				     CONVERTME(DoubleArray), &ypos))
    return NULL;

  npy_intp size = xpos.get_size();
  if ( size != ypos.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "xpos: " << size << " vs ypos: " << ypos.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  IntArray mask;
  npy_intp dim[1];
  dim[0] = npy_intp(size);

  if ( EXIT_SUCCESS != mask.create( 1, dim ) )
    return NULL;

  // Create the mask
  for ( npy_intp ii = 0; ii < size; ii++)
    mask[ii] = regInsideRegion( self->region, xpos[ii], ypos[ii] );

  return mask.return_new_ref();

}


static PyMethodDef pyRegion_methods[] = {

  { (char*)"mask", (PyCFunction) region_mask,
    METH_VARARGS | METH_KEYWORDS, (char*)"Calculate the mask for a set of points: mask = r.mask(x0, x1)" },

  { (char*)"union", (PyCFunction) region_union,
    METH_VARARGS | METH_KEYWORDS, (char*)"Take the union of two regions: r3 = r1.union(r2)" },

  { (char*)"subtract", (PyCFunction) region_subtract,
    METH_VARARGS | METH_KEYWORDS, (char*)"Remove r2 from r1: r3 = r1.subtract(r2)" },

  { (char*)"invert", (PyCFunction) region_invert,
    METH_NOARGS, (char*)"Invert the region: r.invert()" },

  { NULL, NULL, 0, NULL }

};


// Unfortunately we can not use C99-style designated initializers here
// so we need to keep all the "null" values.
//
static PyTypeObject pyRegion_Type = {
  // Note that there is no semicolon after the PyObject_HEAD_INIT macro;
  // one is included in the macro definition.
  PyVarObject_HEAD_INIT(NULL, 0)    // ob_size
  (char*)"PyRegion",             // tp_name
  sizeof(PyRegion),              // tp_basicsize
  0, 		                 // tp_itemsize
  (destructor)pyRegion_dealloc,  // tp_dealloc
  0,                             // tp_print
  0,                             // tp_getattr
  0,                             // tp_setattr
  0,                             // tp_compare
  (reprfunc)pyRegion_str,        // tp_repr
  0,                             // tp_as_number
  0,                             // tp_as_sequence
  0,                             // tp_as_mapping
  0,                             // tp_hash
  0,                             // tp_call
  (reprfunc)pyRegion_str,        // tp_str
  0,                             // tp_getattro
  0,                             // tp_setattro
  0,                             // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, //tp_flags  DO NOT REMOVE
  (char*)"PyRegion reg = Region(region, fileflag=False)",      // tp_doc, __doc__
  0,		                 // tp_traverse
  0,		                 // tp_clear
  0,		                 // tp_richcompare
  0,		                 // tp_weaklistoffset
  0,		                 // tp_iter, __iter__
  0,		                 // tp_iternext
  pyRegion_methods,              // tp_methods
  0,                             // tp_members
  0,                             // tp_getset
  0,                             // tp_base
  0,                             // tp_dict, __dict__
  0,                             // tp_descr_get
  0,                             // tp_descr_set
  0,                             // tp_dictoffset
  0,                             // tp_init, __init__
  0,                             // tp_alloc
  pyRegion_new,                  // tp_new
};


// Needs to be defined after pyRegion_type is defined.
//
// Take the union of the the argument with the object.
//
static PyObject* region_union( PyRegion* self, PyObject* args, PyObject *kwargs )
{

  PyObject *reg_obj2 = NULL;

  static const char *kwlist[] = {"region", NULL};
  if ( !PyArg_ParseTupleAndKeywords( args, kwargs, "O!",
                                     const_cast<char**>(kwlist),
				     &pyRegion_Type, &reg_obj2))
    return NULL;

  PyRegion *reg2 = (PyRegion*)(reg_obj2);

  regRegion *r1 = self->region;
  regRegion *r2 = reg2->region;

  regRegion *combined = regUnionRegion( r1, r2 );
  if ( NULL == combined ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to union the regions" );
    return NULL;
  }

  // I assume we use "N" rather than "O" because, from the docs,
  // "Useful when the object is created by a call to an object
  // constructor in the argument list." so we don't want to
  // increase the reference count (and if "O" is used it
  // leaks memory...)
  //
  PyRegion *out = (PyRegion*) pyRegion_build( &pyRegion_Type, combined );
  return Py_BuildValue((char*)"N", out);

}

// Remove the argument region from the object.
//
static PyObject* region_subtract( PyRegion* self, PyObject* args, PyObject *kwargs )
{

  PyObject *reg_obj2 = NULL;

  static const char *kwlist[] = {"region", NULL};
  if ( !PyArg_ParseTupleAndKeywords( args, kwargs, "O!",
                                     const_cast<char**>(kwlist),
				     &pyRegion_Type, &reg_obj2))
    return NULL;

  PyRegion *reg2 = (PyRegion*)(reg_obj2);

  regRegion *r1 = self->region;
  regRegion *r2 = regInvert( reg2->region );

  regRegion *combined = regIntersectRegion( r1, r2 );
  regFree( r2 );

  if ( NULL == combined ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to subtract a region" );
    return NULL;
  }

  // I assume we use "N" rather than "O" because, from the docs,
  // "Useful when the object is created by a call to an object
  // constructor in the argument list." so we don't want to
  // increase the reference count (and if "O" is used it
  // leaks memory...)
  //
  PyRegion *out = (PyRegion*) pyRegion_build( &pyRegion_Type, combined );
  return Py_BuildValue((char*)"N", out);

}

// Invert a region. The method takes no arguments, so I assume
// args is NULL.
//
static PyObject* region_invert( PyRegion* self, PyObject* args )
{
  regRegion *rold = self->region;
  regRegion *rnew = regInvert( rold );
  if ( NULL == rnew ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to invert the region" );
    return NULL;
  }

  self->region = rnew;
  regFree( rold );

  Py_RETURN_NONE;
}

#ifndef PyMODINIT_FUNC	// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif

static struct PyModuleDef module_region = {
    PyModuleDef_HEAD_INIT,
    "_region",
    "Defines the Region type.",
    -1,
    NULL
};

PyMODINIT_FUNC
PyInit__region(void)
{

  PyObject *m;

  if (PyType_Ready(&pyRegion_Type) < 0)
    return NULL;

  import_array();

  m = PyModule_Create(&module_region);
  if (m == NULL)
    return NULL;

  Py_INCREF(&pyRegion_Type);
  if (PyModule_AddObject(m, (char*)"Region", (PyObject *)&pyRegion_Type) < 0) {
    Py_DECREF(&pyRegion_Type);
    Py_DECREF(m);
  }
  return m;

}
