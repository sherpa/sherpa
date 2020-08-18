// 
//  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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


#include <sherpa/extension.hh>
#include <memory>
#include <iostream>
#include "pileup.hh"	     
#include "PyWrapper.hh"

extern "C" {
  void init_pileup();
}

static int pileup_model_func( double* x0, double* x1, double* res, int len,
			      sherpa::PyWrapper* wrapper ) {

  PyObject* py_function = wrapper->get_pileup_function_ptr( );

  if ( NULL == py_function ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"model function pointer is NULL");
    return EXIT_FAILURE;
  }

  if ( !PyCallable_Check( py_function ) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"model function pointer is not callable");
    return EXIT_FAILURE;
  }

  npy_intp dims[1];
  dims[0] = npy_intp( len );

  DoubleArray n_x0, n_x1, n_res;
  if ( EXIT_SUCCESS != n_x0.create( 1, dims, x0 ) )    
    return EXIT_FAILURE;
  
  if ( EXIT_SUCCESS != n_x1.create( 1, dims, x1 ) )
    return EXIT_FAILURE;

  // call python user model
  PyObject* rv_obj = PyObject_CallFunction( py_function, (char*)"NN",
					    n_x0.new_ref(),
					    n_x1.new_ref() );

  if ( rv_obj == NULL || rv_obj == Py_None ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"model evaluation failed\n");
    return EXIT_FAILURE;
  }
  
  // convert pyobject into double array obj
  CONVERTME(DoubleArray)(rv_obj, &n_res);

  // fill res pointer
  for( int ii = 0; ii < len; ii++ )
    res[ii] = n_res[ii];
  
  Py_DECREF( rv_obj );
  
  return EXIT_SUCCESS;
}

static PyObject* _apply_pileup( PyObject* self, PyObject* args )
{
  DoubleArray arf_source, energ_lo, energ_hi, specresp;
  double exposure_time;
  int max_num_terms;
  double fracexpo;
  double frame_time;
  double alpha;
  double g0;
  double num_regions;
  double psf_frac;
  PyObject *py_function=NULL;

  if ( !PyArg_ParseTuple( args, (char*)"O&diO&O&O&ddddddO",
			  CONVERTME(DoubleArray),
			  &arf_source,
			  &exposure_time,
			  &max_num_terms,
			  CONVERTME(DoubleArray),
			  &energ_lo,
			  CONVERTME(DoubleArray),
			  &energ_hi,
			  CONVERTME(DoubleArray),
			  &specresp,
			  &fracexpo,
			  &frame_time,
			  &alpha,
			  &g0,
			  &num_regions,
			  &psf_frac,
			  &py_function) )
    return NULL;

  if ( ( exposure_time <= 0.0 ) ||
       ( max_num_terms < 1 ) || ( max_num_terms > MAX_NUM_TERMS ) ||
       ( fracexpo < 0.0 ) || ( fracexpo > 1.0 ) ||
       ( frame_time <= 0.0 ) ||
       ( alpha < 0.0 ) || ( alpha > 1.0 ) ||
       ( g0 <= 0.0 ) || ( g0 > 1.0 ) ||
       ( num_regions <= 0.0 ) ||
       ( psf_frac < 0.0 ) || ( psf_frac > 1.0 ) ) {
    PyErr_SetString( PyExc_ValueError, (char*)"invalid pileup parameters" );
    return NULL;
  }

  std::auto_ptr< sherpa::PyWrapper >
    wrapper( new sherpa::PyWrapper( py_function ) );

  DoubleArray results;
  if ( EXIT_SUCCESS != results.create( arf_source.get_ndim(),
				       arf_source.get_dims() ) )
    return NULL;

  DoubleArray pileup_fractions;
  npy_intp pf_dims[1];
  pf_dims[0] = (npy_intp)( max_num_terms + 1 );  // FIXME: check for overflow
  if ( EXIT_SUCCESS != pileup_fractions.zeros( 1, pf_dims ) )
    return NULL;

  // FIXME: check for overflow
  unsigned int num_bins = (unsigned int)( arf_source.get_size() );
  double integral_ae;
  unsigned int num_terms = 0;

  if ( EXIT_SUCCESS != apply_pileup( num_bins, &arf_source[0],
				     &results[0], &pileup_fractions[0],
				     &integral_ae, exposure_time,
				     (unsigned int)max_num_terms, &num_terms,
				     &energ_lo[0], &energ_hi[0], &specresp[0],
				     fracexpo, frame_time,
				     alpha, g0, num_regions, psf_frac,
				     pileup_model_func, wrapper.get() ) ) {
    PyErr_SetString( PyExc_ValueError, (char*)"pileup computation failed" );
    return NULL;
  }

  // We really want to use "I" instead of "i", but the former is broken in
  // pre-2.4.3 Python.  It really doesn't matter, since MAX_NUM_TERMS is much
  // smaller than INT_MAX.

  return Py_BuildValue( (char*)"(NNdi)", results.return_new_ref(),
			pileup_fractions.return_new_ref(), integral_ae,
			(int)num_terms );

}


static PyMethodDef PileupFcts[] = {

  FCTSPEC( apply_pileup, _apply_pileup ),
  
  { NULL, NULL, 0, NULL }

};


SHERPAMOD(_pileup, PileupFcts)
