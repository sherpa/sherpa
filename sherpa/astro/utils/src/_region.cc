//
//  Copyright (C) 2009, 2016, 2020  Smithsonian Astrophysical Observatory
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

extern "C" {
#include "cxcregion.h"
#include <string>
}

typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  regRegion *region;
} PyRegion;


static regRegion* parse_string( char* str, int fileflag );


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
// API could just deal with strings.
//
static PyObject* pyRegion_new(PyTypeObject *type, PyObject *args,
			      PyObject *kwds)
{
  regRegion *reg = NULL;
  char *objS = NULL;
  int fileflag = 0;

  if ( !PyArg_ParseTuple( args, "|si", &objS, &fileflag ) ) {
    PyErr_SetString( PyExc_TypeError, "Region expects optional string and flag" );
    return NULL;
  }

  if (objS == NULL) {
    reg = regCreateEmptyRegion();
  } else {
    reg = parse_string( objS, fileflag );
    if (!reg) {
      std::ostringstream err;
      if (fileflag) {
	err << "unable to read region from: " << objS;
      } else {
	err << "unable to parse region string: " << objS;
      }
      PyErr_SetString( PyExc_ValueError, err.str().c_str() );
      return NULL;
    }
  }

  return pyRegion_build( type, reg );
}

// New datatype dealloc method
static void pyRegion_dealloc(PyRegion* reg) {

  if( reg ) {
    if( reg->region )
      regFree( reg->region );
  }

  PyObject_Del( reg );
}

static PyObject* pyRegion_str(PyRegion* reg) {
  regRegion *region = reg->region;
  const char *ret = "";

  if ( region ) {
    ret = regAllocComposeRegion( region );
  }

  return Py_BuildValue("s", ret);
}


static regRegion* parse_string( char* str, int fileflag ) {

  regRegion *reg = NULL;
  std::string input(str);
  if( fileflag ) {
  	  reg = regReadAsciiRegion( (char*)input.c_str() , 0 ); // Verbosity set to 0
  } else {
	  reg = regParse( (char*)input.c_str() );
  }

  return reg;
}

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
  (char*)"PyRegion reg = Region('Circle(256,256,25)')",      // tp_doc, __doc__
  0,		                 // tp_traverse
  0,		                 // tp_clear
  0,		                 // tp_richcompare
  0,		                 // tp_weaklistoffset
  0,		                 // tp_iter, __iter__
  0,		                 // tp_iternext
  0,                             // tp_methods
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


// Return the mask indicating whether the array points
// lie inside or outside the region.
//
static PyObject* region_mask( PyObject* self, PyObject* args )
{

  PyObject *reg_obj = NULL;
  DoubleArray xpos, ypos;

  if ( !PyArg_ParseTuple( args, (char*)"O!O&O&",
			  &pyRegion_Type, &reg_obj,
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

  // Extract the region
  PyRegion *region = (PyRegion*)(reg_obj);

  for ( npy_intp ii = 0; ii < size; ii++)
    mask[ii] = regInsideRegion( region->region, xpos[ii], ypos[ii] );

  return mask.return_new_ref();

}


// If reverse is 0 then       reg1&reg2
//                 otherwise  reg1&!reg2
//
static PyObject* region_combine( PyObject* self, PyObject* args )
{

  PyObject *reg_obj1 = NULL;
  PyObject *reg_obj2 = NULL;

  int reverse = 0;

  if ( !PyArg_ParseTuple( args, (char*)"O!O!|i",
			  &pyRegion_Type, &reg_obj1,
			  &pyRegion_Type, &reg_obj2,
			  &reverse))
    return NULL;

  PyRegion *reg1 = (PyRegion*)(reg_obj1);
  PyRegion *reg2 = (PyRegion*)(reg_obj2);

  regRegion *r1 = reg1->region;
  regRegion *r2 = reg2->region;

  if( reverse ) {
    r2 = regInvert( r2 );
  }

  regRegion *combined = regCombineRegion( r1, r2 );
  if ( NULL == combined ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to combine regions successfully" );
    return NULL;
  }

  if( reverse ) {
    regFree( r2 );
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


/*
static PyObject* region_parse( PyObject* self, PyObject* args )
{

  int fileflag = 1;
  char* str = NULL;
  PyObject *cobj = NULL;
  PyRegion *region = NULL;
  regRegion *reg = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s|i", &str, &fileflag))
    return NULL;

  reg = parse_string( str, fileflag );

  if ( NULL == reg) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to parse region string successfully" );
    return NULL;
  }

  // In order to pass in regRegion* as a Python argument, use void* trick
  cobj = PyCObject_FromVoidPtr((void*)reg, NULL);

  region = (PyRegion*) pyRegion_new( &pyRegion_Type,
				     Py_BuildValue((char*)"(O)", cobj ),
				     NULL );

  Py_XINCREF(region);

  return (PyObject*)region;

}
*/

static PyMethodDef RegionFcts[] = {

  { (char*)"region_mask", (PyCFunction) region_mask,
    METH_VARARGS, (char*)"Calculate the mask for a region and set of points" },

  { (char*)"region_combine", (PyCFunction) region_combine,
    METH_VARARGS, (char*)"Combine regions" },

  /*
    { (char*)"region_parse", (PyCFunction)region_parse,
    METH_VARARGS, (char*)"Parse and create a region using CXC region lib" },
  */

  { NULL, NULL, 0, NULL }

};

#ifndef PyMODINIT_FUNC	// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif

static struct PyModuleDef module_region = {
PyModuleDef_HEAD_INIT,
"_region",
NULL,
-1,
RegionFcts
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

  PyModule_AddObject(m, (char*)"Region", (PyObject *)&pyRegion_Type);
  return m;

}
