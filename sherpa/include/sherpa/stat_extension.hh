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

#ifndef __sherpa_stat_extension_hh__
#define __sherpa_stat_extension_hh__

#include <sherpa/extension.hh>

namespace sherpa { namespace stats {

  template <typename ArrayType,
	    typename DataType,
	    int (*ErrFunc)( npy_intp num, const ArrayType& yraw,
			    ArrayType& err )>
  PyObject* staterrfct( PyObject* self, PyObject* args )
  {

    ArrayType yraw;

    if ( !PyArg_ParseTuple( args, (char*)"O&",
			    (converter)convert_to_array< ArrayType >, &yraw ) )
      return NULL;

    ArrayType err;
    if ( EXIT_SUCCESS != err.create( yraw.get_ndim(), yraw.get_dims() ) )
      return NULL;

    if ( EXIT_SUCCESS != ErrFunc( yraw.get_size(), yraw, err ) ) {
      PyErr_SetString( PyExc_ValueError,
		       (char*)"calculation of errors has failed using current statistic");
      return NULL;
    }
    
    return err.return_new_ref();

  }


  template <typename ArrayType,
	    typename DataType,
	    int (*StatFunc)( npy_intp num, const ArrayType& yraw,
			     const ArrayType& model,
			     const ArrayType& staterror,
			     const ArrayType& syserror,
			     const ArrayType& weight,
			     ArrayType& dev, DataType& val)>
  PyObject* statfct_noerr( PyObject* self, PyObject* args )
  {

    ArrayType yraw;
    ArrayType model;
    ArrayType staterror;  // Not used in calculation
    ArrayType syserror;   // Not used in calculation
    ArrayType weight;
    PyObject* dummy = NULL;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&OOO&",
			    (converter)convert_to_array< ArrayType >, &yraw,
			    (converter)convert_to_array< ArrayType >, &model,
			    &dummy,
			    &dummy,
			    (converter)array_or_none< ArrayType >, &weight ) )
      return NULL;

    npy_intp nelem = yraw.get_size();

    if ( ( model.get_size() != nelem ) ||
	 ( weight && ( weight.get_size() != nelem ) ) ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"statistic input array sizes do not match" );
      return NULL;
    }

    // dev is needed for temporary storage
    ArrayType dev;
    if ( EXIT_SUCCESS != dev.create( yraw.get_ndim(), yraw.get_dims() ) )
      return NULL;

    DataType val = 0.0;

    if ( EXIT_SUCCESS != StatFunc( nelem, yraw, model, staterror, syserror,
				   weight, dev, val ) ) {
      PyErr_SetString( PyExc_ValueError, (char*)"statistic calculation failed");
      return NULL;
    }

    // Py_None MUST be incremented before being returned!!
    Py_INCREF(Py_None);
    return Py_BuildValue( (char*)"(dO)", val, Py_None );

  }


  template <typename ArrayType,
	    typename DataType,
	    int (*StatFunc)( npy_intp num, const ArrayType& yraw,
			     const ArrayType& model,
			     const ArrayType& staterror,
			     const ArrayType& syserror,
			     const ArrayType& weight,
			     ArrayType& dev, DataType& val, 
			     DataType& trunc_value )>
  PyObject* statfct( PyObject* self, PyObject* args )
  {

    ArrayType yraw;
    ArrayType model;
    ArrayType staterror;
    ArrayType syserror;
    ArrayType weight;
    double trunc_value = -1.0;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&d",
			    (converter)convert_to_array< ArrayType >, &yraw,
			    (converter)convert_to_array< ArrayType >, &model,
			    (converter)convert_to_array< ArrayType >,
			    &staterror,
			    (converter)array_or_none< ArrayType >, &syserror,
			    (converter)array_or_none< ArrayType >, &weight,
			    &trunc_value) )
      return NULL;

    npy_intp nelem = yraw.get_size();

    if ( ( model.get_size() != nelem ) ||
	 ( staterror.get_size() != nelem ) ||
	 ( syserror && ( syserror.get_size() != nelem ) ) ||
	 ( weight && ( weight.get_size() != nelem ) ) ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"statistic input array sizes do not match" );
      return NULL;
    }

    ArrayType dev;
    if ( EXIT_SUCCESS != dev.create( yraw.get_ndim(), yraw.get_dims() ) )
      return NULL;

    DataType val = 0.0;

    if ( EXIT_SUCCESS != StatFunc( nelem, yraw, model, staterror, syserror,
				   weight, dev, val, trunc_value ) ) {
      PyErr_SetString( PyExc_ValueError, (char*)"statistic calculation failed");
      return NULL;
    }

    return Py_BuildValue( (char*)"(dN)", val, dev.return_new_ref() );

  }


}  }  /* namespace stats, namespace sherpa */


#define _STATFCTPTR(name) \
  sherpa::stats::name< SherpaFloatArray, SherpaFloatArray, SherpaFloat, \
                       npy_intp >

#define _STATFCTSPEC(name, ftype) \
  FCTSPEC(name, (sherpa::stats::ftype< SherpaFloatArray, SherpaFloat, \
                                       _STATFCTPTR(name) >))

#define STATERRFCT(name)	_STATFCTSPEC(name, staterrfct)
#define STATFCT(name)		_STATFCTSPEC(name, statfct)
#define STATFCT_NOERR(name)	_STATFCTSPEC(name, statfct_noerr)


#endif /* __sherpa_stat_extension_hh__ */
