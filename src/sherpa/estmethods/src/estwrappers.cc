// 
//  Copyright (C) 2007, 2016  Smithsonian Astrophysical Observatory
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
#include <cmath>
#include <cstdlib>
#include <float.h>
#include "estutils.hh"

#ifdef __SUNPRO_CC
#include <sunmath.h>
#define NAN quiet_nan(0)
#endif

extern "C" {
  void init_est_funcs();
}

// Keep pointers to the statistic and fitting functions used
// by these methods.

static PyObject* stat_func = NULL;
static PyObject* fit_func = NULL;

// These objects are class objects that are references to various
// estmethod module exceptions.  The idea is that from this C++ code, 
// we want to raise particular types of exceptions--whether a new 
// minimum is found, hard parameter limits are hit, and so on.
// To raise new Python exceptions of these types, we need references
// to the exception classes from which the new exceptions will derive;
// the estmethods exceptions are not pre-defined.  Therefore we need
// to give a reference to the relevant class object, to 
// PyErr_NewException().

static PyObject* est_error = NULL;     // Base error, no more info
static PyObject* est_hardmin = NULL;   // Hard parameter min hit
static PyObject* est_hardmax = NULL;   // Hard parameter max hit
static PyObject* est_newmin = NULL;    // New minimum statistic hit
static PyObject* est_maxiter = NULL;   // Maximum number of iterations
static PyObject* est_hitnan = NULL;    // NaN value encountered


static double statfcn( double* pars, int npars )
{

  if ( NULL == stat_func ) {
    PyErr_SetString( PyExc_SystemError,
		     (char*)"statistic callback is not set (NULL pointer)" );
    return NAN;
  }

  npy_intp dims[1];
  dims[0] = npy_intp( npars );

  DoubleArray pars_obj;
  if ( EXIT_SUCCESS != pars_obj.create( 1, dims, pars ) )
    return DBL_MAX;

  PyObject* rv_obj = NULL;
  if ( NULL == ( rv_obj = PyObject_CallFunction( stat_func, (char*)"N",
						 pars_obj.new_ref() ) ) )
    return NAN;
  
  if ( !PyFloat_Check( rv_obj ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"statistic callback did not return a float" );
    Py_DECREF( rv_obj );
    return NAN;
  }
  
  double rv = PyFloat_AsDouble( rv_obj );
  Py_DECREF( rv_obj );

  return rv;

}


static double fitfcn( double (*dummyfunc)(double*, int),
		      double* pars, double* parmins, 
		      double* parmaxs, int npars, int parnum )
{

  // If all has gone correctly, dummyfunc is just a pointer
  // to the function that came from the stat_func PyObject
  // above.  So instead of wrapping up dummyfunc, just pass
  // a reference to the original stat_func object when
  // calling the function used in this function.

  if ( NULL == stat_func ) {
    PyErr_SetString( PyExc_SystemError,
		     (char*)"statistic callback is not set (NULL pointer)" );
    return NAN;
  }

  if ( NULL == fit_func ) {
    PyErr_SetString( PyExc_SystemError,
		     (char*)"fit callback is not set (NULL pointer)" );
    return NAN;
  }

  npy_intp dims[1];
  dims[0] = npy_intp( npars );

  DoubleArray pars_obj;
  if ( EXIT_SUCCESS != pars_obj.create( 1, dims, pars ) )
    return NAN;

  DoubleArray parmins_obj;
  if ( EXIT_SUCCESS != parmins_obj.create( 1, dims, parmins ) )
    return NAN;

  DoubleArray parmaxs_obj;
  if ( EXIT_SUCCESS != parmaxs_obj.create( 1, dims, parmaxs ) )
    return NAN;
				   
  PyObject* rv_obj = NULL;
  if ( NULL == ( rv_obj = PyObject_CallFunction( fit_func, (char*)"NNNi", 
						 pars_obj.new_ref(),
						 parmins_obj.new_ref(),
						 parmaxs_obj.new_ref(),
						 parnum ) ) )
    return NAN;

  if ( !PyFloat_Check( rv_obj ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"minimize callback did not return a float" );
    Py_DECREF( rv_obj );
    return NAN;
  }
  
  double rv = PyFloat_AsDouble( rv_obj );
  Py_DECREF( rv_obj );

  return rv;
}

static void _get_exception_objects()
{
  // If any of the exception class objects above haven't yet been
  // set, then do so.
  if (NULL == est_error ||
      NULL == est_hardmin ||
      NULL == est_hardmax ||
      NULL == est_newmin ||
      NULL == est_maxiter ||
      NULL == est_hitnan) {
    PyObject* est_module = PyImport_AddModule((char*)"sherpa.estmethods");
    if (NULL == est_module)
      return;
    PyObject* est_dict = PyModule_GetDict(est_module);
    if (NULL == est_dict)
      return;

    // Do I need to call Py_INCREF() on these six PyObjects?
    // Answer:  no, these are only borrowed references, therefore
    // function doesn't own them.
    // http://docs.python.org/api/dictObjects.html says PyDict_GetItemString()
    // returns borrowed reference.
    // http://docs.python.org/api/refcountDetails.html defines the
    // difference between ownership of borrowed and new references.
    // SMD 3/29/07
    if (NULL == est_error)
      est_error = PyDict_GetItemString(est_dict, (char*)"EstMethodError");
    if (NULL == est_hardmin)
      est_hardmin = PyDict_GetItemString(est_dict, (char*)"EstHardMin");
    if (NULL == est_hardmax)
      est_hardmax = PyDict_GetItemString(est_dict, (char*)"EstHardMax");
    if (NULL == est_newmin)
      est_newmin = PyDict_GetItemString(est_dict, (char*)"EstNewMin");
    if (NULL == est_maxiter)
      est_maxiter = PyDict_GetItemString(est_dict, (char*)"EstMaxIter");
    if (NULL == est_hitnan)
      est_hitnan = PyDict_GetItemString(est_dict, (char*)"EstNaN");
  }
  else
    return;
}


static void _raise_python_error(const char* base_message, 
				const est_return_code status)
{
  _get_exception_objects();
  PyObject* est_exception = NULL;
  PyObject* par_number = NULL;
  int set_correct_error = EXIT_FAILURE;

  switch(status.status) {

  case EST_HARDMIN:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstHardMin",
				       est_hardmin,
				       NULL);
    par_number = PyLong_FromLong(status.par_number);
    if (NULL == est_hardmin || NULL == est_exception || NULL == par_number)
      ;
    else {
      PyErr_SetObject(est_exception, par_number);
      set_correct_error = EXIT_SUCCESS;
    }
    break;

  case EST_HARDMAX:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstHardMax",
				       est_hardmax,
				       NULL);
    par_number = PyLong_FromLong(status.par_number);
    if (NULL == est_hardmax || NULL == est_exception || NULL == par_number)
      ;
    else {
      PyErr_SetObject(est_exception, par_number);
      set_correct_error = EXIT_SUCCESS;
    }
    break;

  case EST_NEWMIN:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstNewMin",
				       est_newmin,
				       NULL);
    if (NULL == est_newmin || NULL == est_exception)
      ;
    else {
      PyErr_SetString(est_exception, 
		      (char*)"new minimum found, restarting error method");
      set_correct_error = EXIT_SUCCESS;
    }
    break;
    
  case EST_MAXITER:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstMaxIter",
				       est_maxiter,
				       NULL);
    if (NULL == est_maxiter || NULL == est_exception)
      ;
    else {
      PyErr_SetString(est_exception, 
		      (char*)"maximum number of iterations in scaling function");
      set_correct_error = EXIT_SUCCESS;
    }
    break;

  case EST_HITNAN:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstNaN",
				       est_hitnan,
				       NULL);
    if (NULL == est_hitnan || NULL == est_exception)
      ;
    else {
      PyErr_SetString(est_exception, 
		      (char*)"NaN encountered during error method");
      set_correct_error = EXIT_SUCCESS;
    }
    break;

  default:
    est_exception = PyErr_NewException((char*)"sherpa.estmethods.EstMethodError",
				      est_error,
				      NULL);
    if (NULL == est_error || NULL == est_exception)
      ;
    else {
      PyErr_SetString(est_exception, base_message);
      set_correct_error = EXIT_SUCCESS;
    }
  }

  if (EXIT_FAILURE == set_correct_error) {
    PyErr_SetString(PyExc_RuntimeError, base_message);
  }

  Py_XDECREF(par_number);
  Py_XDECREF(est_exception);
  return;
}

static PyObject* _wrap_info_matrix( PyObject* self, PyObject* args )
{

  DoubleArray pars;
  DoubleArray pars_mins;
  DoubleArray pars_maxs;
  DoubleArray pars_hardmins;
  DoubleArray pars_hardmaxs;
  double sigma;
  double eps;
  int maxiters;
  double remin;

  if ( !PyArg_ParseTuple( args,(char *)"O&O&O&O&O&ddidO",
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_mins,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_maxs,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_hardmins,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_hardmaxs,
			  &sigma,
			  &eps,
			  &maxiters,
			  &remin,
			  &stat_func ) )
    return NULL;

  npy_intp nelem = pars.get_size();

  if ( nelem != pars_mins.get_size() ||
       nelem != pars_maxs.get_size() ||
       nelem != pars_hardmins.get_size() ||
       nelem != pars_hardmaxs.get_size() ) {
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"input array sizes do not match" );
    return NULL;
  }

  npy_intp dims[2];
  dims[0] = nelem;
  dims[1] = nelem;

  PyObject* info_obj = NULL;
  if ( NULL == ( info_obj = PyArray_New( &PyArray_Type, 2, dims,
					 NPY_DOUBLE, NULL, NULL, 0, NPY_ARRAY_CARRAY,
					 NULL ) ) )
    return NULL;

  est_return_code status = info_matrix( &(pars[0]), int( nelem ),
					&(pars_mins[0]), int( nelem ),
					&(pars_maxs[0]), int( nelem ),
					&(pars_hardmins[0]), int( nelem ),
					&(pars_hardmaxs[0]), int( nelem ),
					static_cast< double* >( PyArray_DATA( (PyArrayObject*) info_obj ) ),
					int( nelem ), int( nelem ),
					sigma,
					eps,
					maxiters,
					remin,
					statfcn );

  if ( EST_SUCCESS != status.status ) { 
    if ( NULL == PyErr_Occurred() )
      _raise_python_error((char*)"covariance failed", status);
    Py_DECREF( info_obj );
    return NULL; 
  }
 
  // Use "N" (i.e. don't increment reference count) so that ownership of
  // info_obj passes to caller
  return Py_BuildValue( (char*)"N", info_obj );

}


static PyObject* _wrap_projection( PyObject* self, PyObject* args )
{

  DoubleArray pars;
  DoubleArray pars_mins;
  DoubleArray pars_maxs;
  DoubleArray pars_hardmins;
  DoubleArray pars_hardmaxs;
  IntArray parnums;
  double sigma;
  double eps;
  double tol;
  int maxiters;
  double remin;

  if ( !PyArg_ParseTuple( args,(char *)"O&O&O&O&O&dddidO&OO",
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_mins,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_maxs,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_hardmins,
			  (converter)sherpa::
			  convert_to_contig_array< DoubleArray >,
			  &pars_hardmaxs,
			  &sigma,
			  &eps,
			  &tol,
			  &maxiters,
			  &remin,
			  (converter)sherpa::
			  convert_to_contig_array< IntArray >,
			  &parnums,
			  &stat_func,
			  &fit_func ) )
    return NULL;

  npy_intp nelem = pars.get_size();
  npy_intp parnumsize = parnums.get_size();

  if ( nelem != pars_mins.get_size() ||
       nelem != pars_maxs.get_size() ||
       nelem != pars_hardmins.get_size() ||
       nelem != pars_hardmaxs.get_size() ) {
    PyErr_SetString( PyExc_RuntimeError,
		     (char*)"input array sizes do not match" );
    return NULL;
  }
 
  npy_intp dims[1];
  dims[0] = parnumsize;

  DoubleArray pars_elow;
  if ( EXIT_SUCCESS != pars_elow.create( 1, dims ) )
    return NULL;

  DoubleArray pars_ehi;
  if ( EXIT_SUCCESS != pars_ehi.create( 1, dims ) )
    return NULL;

  IntArray pars_eflags;
  if ( EXIT_SUCCESS != pars_eflags.create( 1, dims ) )
    return NULL;
  
  est_return_code status = projection( &(pars[0]), int( nelem ),
				       &(pars_mins[0]), int( nelem ),
				       &(pars_maxs[0]), int( nelem ),
				       &(pars_hardmins[0]), int( nelem ),
				       &(pars_hardmaxs[0]), int( nelem ),
				       &(pars_elow[0]), int( parnumsize ), 
				       &(pars_ehi[0]), int( parnumsize ),
				       &(pars_eflags[0]), int( parnumsize ),
				       sigma,
				       eps,
				       tol,
				       maxiters,
				       remin,
				       &(parnums[0]), int ( parnumsize ),
				       statfcn,
				       fitfcn );

  if ( EST_SUCCESS != status.status ) { 
    if ( NULL == PyErr_Occurred() )
      _raise_python_error((char*)"projection failed", status);
    return NULL;
  }

  return Py_BuildValue( (char*)"(NNNi)", pars_elow.return_new_ref(),
			pars_ehi.return_new_ref(), 
			pars_eflags.return_new_ref(),
			status.nfits );

}


static PyMethodDef WrapperFcts[] = {

  FCTSPEC( info_matrix, _wrap_info_matrix ),
  FCTSPEC( projection, _wrap_projection ),

  { NULL, NULL, 0, NULL }

};


SHERPAMOD(_est_funcs, WrapperFcts)
