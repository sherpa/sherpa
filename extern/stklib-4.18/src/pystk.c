/*                                                                
**  Copyright (C) 2011,2015-2016,2025  Smithsonian Astrophysical Observatory 
*/                                                                
/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


#include "pystk.h"
#include <stdlib.h>
#include <stack.h>
#include <string.h>

#include <unistd.h>
#include <stdio.h>

static PyMethodDef stkMethods[] =
{
  /* Each entry in the groupMethods array is a PyMethodDef structure containing
   * 1) the Python name,
   * 2) the C-function that implements the function,
   * 3) flags indicating whether or not keywords are accepted for this function,
   * and 4) The docstring for the function.
   */
  { "build", (PyCFunction)_stk_build, METH_VARARGS,
    "Example: \n"
    ">>> foo = build(\"a,b,c,d\")\n"
    ">>> foo\n"
    "['a', 'b', 'c', 'd']\n" },
  { NULL, NULL, 0, NULL }
};


/*
 * Initialize the module
 * */
#if PY_MAJOR_VERSION >= 3
#define MOD_DEF(ob, name, doc, methods) \
        static struct PyModuleDef moduledef = { \
          PyModuleDef_HEAD_INIT, name, doc, -1, methods }; \
        ob = PyModule_Create(&moduledef);
#endif

PyMODINIT_FUNC PyInit_stk(void)
{

  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT, "stk", "stk", -1, stkMethods };

  PyObject *mod = PyModule_Create(&moduledef);
  if (mod == NULL) {
    return NULL;
  }

#ifdef Py_GIL_DISABLED
  PyUnstable_Module_SetGIL(mod, Py_MOD_GIL_NOT_USED);
#endif

  return mod;
}

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Function Definitions
 * * * * * * * * * * * * * * * * * * * * * * * */
/**
 * Function: _stk_build
 *   Wrapper function for stk_build method.
 *
 *   Uses stk library to open and process the provided stack list
 *   generating a List of the stack contents.
 *
 * Returns
 *   Python List containing all stack entries.
 */
static PyObject* _stk_build(PyObject *self, PyObject *args)
{
  char *buff;
  int  status = 0;
  (void)self; /* avoid compiler warning */

  /* Parse arguments, allow str, unicode or bytes input string.  */
#if PY_MAJOR_VERSION >= 3
  if (!PyArg_ParseTuple(args, "s", &buff))
  {
    /* attempt using bytes type before giving up. */
    PyObject *type, *value, *traceback;
    PyErr_Fetch(&type, &value, &traceback);
    PyErr_Clear();
    if (!PyArg_ParseTuple(args, "y", &buff))
    {
      /* did not work either.. restore original error */
      PyErr_Restore(type, value, traceback);
      status = 1;
    }
  }
#else
  if (!PyArg_ParseTuple(args, "s", &buff))
    status = 1;
#endif

  /* handle argument parsing error */
  if ( status == 1 )
  {
    PyErr_SetString(PyExc_Exception, "Could not parse arguments.");
    return NULL;
  }

  if ( ( NULL == buff ) || ( 0 == strlen( buff ) ) )
  {
    PyErr_SetString(PyExc_Exception, "Empty stack string.");
    return NULL;
  }

  Stack *stk;

  /* We trash stk_lib's stdout/stderr and issue our own exception */

  int pipe_fild[2];
  pipe( pipe_fild );
  int orig_err = dup( fileno(stderr ));
  fflush(stderr);
  dup2( pipe_fild[1], fileno(stderr));

  /* Open the stack */
  stk = stk_build( buff );

  /* Return stderr */
  fflush(stderr);
  dup2( orig_err, fileno(stderr));
  close(pipe_fild[0]);
  close(pipe_fild[1]);
  close(orig_err);

  /* Check stack opened OK */
  if ( NULL == stk )
  {
    long maxl = strlen( buff ) + 100;
    char *ibuff = ( char*) calloc( maxl, sizeof(char));
    sprintf( ibuff, "Cannot build stack from string '%s'\n", buff );
    PyErr_SetString(PyExc_IOError, ibuff);
    free(ibuff);
    return NULL;
  }

  /* Get number of stack entries.. require > 0 */
  int num_elem = stk_count(stk);
  if ( 0 == num_elem )
  {
    PyErr_SetString(PyExc_ValueError, "Stack has 0 elements");
    return NULL;
  }

  if ( (1 == num_elem ) && (0 == strlen( stk_read_num( stk, 1 )) ))
  {
    PyErr_SetString(PyExc_ValueError, "Stack has only 1 element and it is blank");
    return NULL;
  }

  /* Create output List */
  PyObject *pylist;
  if ( NULL == (pylist = PyList_New(  (Py_ssize_t) 0 )))
  {
    PyErr_SetString(PyExc_Exception, "Failed to create new list");
    return NULL;
  }

  /* Transfer Stack content to List */
  int ii;
  for ( ii = 1; ii <= num_elem; ii++ )
  {
    /* get stack record */
    char *ibuff = stk_read_num( stk, ii );
    if ( NULL == ibuff )
    {
      PyErr_SetString(PyExc_IndexError, "Invalid stack_read_num");
      return NULL;
    }

    /* Convert to python string object */
    PyObject *pstr = PyUnicode_FromString( ibuff );
    if ( NULL == pstr )
    {
      PyErr_SetString(PyExc_ValueError, "Cannot convert to python string");
      return NULL;
    }

    /* Add to List */
    if ( 0 != PyList_Append( pylist, pstr ) )
    {
      PyErr_SetString( PyExc_Exception, "Failed to append string to list");
      return NULL;
    }

    stk_read_free( ibuff );

  } /* End for ii */

  stk_close( stk );

  return( pylist );
}
