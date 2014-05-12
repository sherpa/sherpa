//_C++_INSERT_SAO_COPYRIGHT_HERE_(2009)_
//_C++_INSERT_GPL_LICENSE_HERE_

#include "sherpa/extension.hh"
#include <sstream>
#include <iostream>
#include <string>

extern "C" {
  #include "cxcregion.h"
  #include "ascdm.h"
  void init_region();
}

typedef struct {
  // Note that there is no semicolon after the PyObject_HEAD macro;
  // one is included in the macro definition.
  PyObject_HEAD
  regRegion *region;
} PyRegion;

static PyObject* pyRegion_new(PyTypeObject *type, PyObject *args,
			      PyObject *kwds)
{
  PyRegion *self = NULL;
  PyObject *cobj = NULL;
  regRegion *reg = NULL;

  if ( !PyArg_ParseTuple( args, "|O", &cobj ) )
    return NULL;

  // In order to pass in regRegion* as a Python argument, use void* trick
  if(cobj) {
    if( PyString_Check( cobj ) ) {
      reg = regParse( PyString_AS_STRING( cobj ) );
//       if ( !reg ) {
// 	PyErr_SetString( PyExc_TypeError,
// 			 (char*)"unable to parse region string successfully" );
// 	return NULL;
//       }
    }
    else if( PyCObject_Check( cobj ) )
      reg = (regRegion*) PyCObject_AsVoidPtr( cobj );
  }

  self = (PyRegion*)type->tp_alloc( type, 0 );

  if(self)
    self->region = reg;

  return (PyObject*)self;
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
  
  if( !region )
    return PyString_FromString((char*)"");

  return PyString_FromString((char*) regAllocComposeRegion( region ));
}


static regRegion* parse_string( char* str, int fileflag ) {
  
  regRegion *reg = NULL;

  std::string input(str);
  if( fileflag )
    input = "region(" + input + ")";

  reg = dmRegParse( (char*)input.c_str() );
 
  return reg;
}

static PyTypeObject pyRegion_Type = {
  // Note that there is no semicolon after the PyObject_HEAD_INIT macro;
  // one is included in the macro definition.
  PyObject_HEAD_INIT(NULL)
  0,                             // ob_size
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
  (char*)"PyRegion object, reg = Region('Circle(256,256,25)')",      // tp_doc, __doc__
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


static PyObject* region_mask( PyObject* self, PyObject* args )
{

  char* str = NULL;
  int fileflag, reverse;
  DoubleArray xpos, ypos;

  PyObject *reg_obj = NULL;
  PyObject *cobj = NULL;
  PyRegion *region = NULL;
  regRegion *reg = NULL;
  
  if ( !PyArg_ParseTuple( args, (char*)"OsO&O&ii",
			  &reg_obj,
			  &str,
			  CONVERTME(DoubleArray), &xpos,
			  CONVERTME(DoubleArray), &ypos,
			  &fileflag,
			  &reverse))
    return NULL;

  if ( xpos.get_size() != ypos.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "xpos: " << xpos.get_size() << " vs ypos: " << ypos.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  IntArray mask;
  npy_intp dim[1];
  dim[0] = npy_intp(xpos.get_size());

  if ( EXIT_SUCCESS != mask.create( 1, dim ) )
    return NULL;

  reg = parse_string( str, fileflag );

  if ( NULL == reg) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"unable to parse region string successfully" );
    return NULL;
  }

  npy_intp size = xpos.get_size();
  for ( npy_intp ii = 0; ii < size; ii++)
    mask[ii] = regInsideRegion( reg, xpos[ii], ypos[ii] );

  if( reverse ) {
    regRegion *tmp = reg;
    reg = regInvert( tmp );
    if( tmp ) regFree( tmp );
  }

  // If possible, append the new region string to the stack  
  if( reg_obj != Py_None ) {
    region = (PyRegion*)(reg_obj);
    regRegion *tmp = reg;

    reg = regCombineRegion( region->region, tmp );
    
    if ( NULL == reg) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"unable to combine regions successfully" );
      return NULL;
    }

    pyRegion_dealloc( region );
    if( tmp ) regFree( tmp );
  }

  // In order to pass in regRegion* as a Python argument, use void* trick
  cobj = PyCObject_FromVoidPtr((void*)reg, NULL);

  region = (PyRegion*) pyRegion_new( &pyRegion_Type,
				     Py_BuildValue((char*)"(O)", cobj ),
				     NULL );

  Py_XINCREF(region);

  return Py_BuildValue((char*)"ON", region, mask.return_new_ref());
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

  { (char*)"region_mask", (PyCFunction)region_mask,
    METH_VARARGS, (char*)"Create a region using CXC region lib" },

  /*
    { (char*)"region_parse", (PyCFunction)region_parse,
    METH_VARARGS, (char*)"Parse and create a region using CXC region lib" },
  */

  { NULL, NULL, 0, NULL }
  
};

#ifndef PyMODINIT_FUNC	// declarations for DLL import/export
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC 
init_region(void) 
{ 

  PyObject *m;
  
  if (PyType_Ready(&pyRegion_Type) < 0)
    return;

  import_array();
  
  m = Py_InitModule3((char*)"_region", RegionFcts, NULL);
  
  if (m == NULL)
    return;

  Py_INCREF(&pyRegion_Type);

  PyModule_AddObject(m, (char*)"Region", (PyObject *)&pyRegion_Type);

}
