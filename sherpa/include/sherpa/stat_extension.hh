// 
//  Copyright (C) 2009, 2015  Smithsonian Astrophysical Observatory
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
			     ArrayType& fvec, DataType& val, 
			     DataType& trunc_value )>
  PyObject* statfct( PyObject* self, PyObject* args )
  {

    ArrayType yraw;
    ArrayType model;
    ArrayType staterror;
    ArrayType syserror;
    ArrayType weight;
    DataType trunc_value = 1.0e-25;

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

    ArrayType fvec;
    if ( EXIT_SUCCESS != fvec.create( yraw.get_ndim(), yraw.get_dims() ) )
      return NULL;

    DataType val = 0.0;

    if ( EXIT_SUCCESS != StatFunc( nelem, yraw, model, staterror, syserror,
				   weight, fvec, val, trunc_value ) ) {
      PyErr_SetString( PyExc_ValueError, (char*)"statistic calculation failed");
      return NULL;
    }

    return Py_BuildValue( (char*)"(dN)", val, fvec.return_new_ref() );

  }

    template <typename ArrayType, typename DataType, typename iArrayType,
              int (*StatFunc)( npy_intp num, const ArrayType& yraw,
                               const ArrayType& model,
                               const iArrayType& data_size,
                               const ArrayType& exposure_time,
                               const ArrayType& bkg,
                               const ArrayType& backscale_ratio,
                               ArrayType& fvec, DataType& val,
                               DataType trunc_value  )>
    PyObject* wstatfct( PyObject* self, PyObject* args ) {

      ArrayType yraw;
      ArrayType model;
      iArrayType data_size;
      ArrayType exposure_time;
      ArrayType bkg;
      ArrayType backscale_ratio;
      DataType trunc_value = 1.0e-25;

      if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&O&d",
                              CONVERTME( ArrayType ), &yraw,
                              CONVERTME( ArrayType ), &model,
                              CONVERTME( iArrayType ), &data_size,
                              CONVERTME( ArrayType ), &exposure_time,
                              CONVERTME( ArrayType ), &bkg,
                              CONVERTME( ArrayType ), &backscale_ratio,
                              &trunc_value ) )
        return NULL;

      const npy_intp nelem = yraw.get_size();

      if ( ( model.get_size( ) != nelem ) || bkg.get_size( ) != nelem ||
           ( 2 * data_size.get_size( ) != exposure_time.get_size( ) ) ||
           data_size.get_size( ) != backscale_ratio.get_size( ) ) {
        PyErr_SetString( PyExc_TypeError,
                         (char*)"statistic input array sizes do not match" );
        return NULL;
      }

      npy_intp sum_data_size = 0;
      for ( npy_intp ii = 0; ii < data_size.get_size( ); ++ii )
        sum_data_size += data_size[ ii ];
      if ( nelem != sum_data_size ) {
        PyErr_SetString( PyExc_TypeError,
                         (char*)"data size do not match" );
        return NULL;
      }

      ArrayType fvec;
      if ( EXIT_SUCCESS != fvec.create( yraw.get_ndim(), yraw.get_dims() ) )
        return NULL;
      DataType val = 0.0;
      if ( EXIT_SUCCESS != StatFunc( nelem, yraw, model, data_size, 
                                     exposure_time, bkg, backscale_ratio,
                                     fvec, val, trunc_value ) ) {
        PyErr_SetString( PyExc_ValueError,
                         (char*)"statistic calculation failed");
        return NULL;
      }

      return Py_BuildValue( (char*)"(dN)", val, fvec.return_new_ref() );

    }


}  }  /* namespace stats, namespace sherpa */


#define _STATFCTPTR(name) \
  sherpa::stats::name< SherpaFloatArray, SherpaFloatArray, SherpaFloat, \
                       npy_intp >
#define _WSTATFCTPTR(name) \
  sherpa::stats::name< SherpaFloatArray, SherpaFloatArray, SherpaFloat, \
                       npy_intp, IntArray >

#define _STATFCTSPEC(name, ftype) \
  FCTSPEC(name, (sherpa::stats::ftype< SherpaFloatArray, SherpaFloat, \
                                       _STATFCTPTR(name) >))
#define _WSTATFCTSPEC(name, ftype) \
  FCTSPEC(name, (sherpa::stats::ftype< SherpaFloatArray, SherpaFloat, \
                 IntArray, _WSTATFCTPTR(name) >))

#define STATERRFCT(name)	_STATFCTSPEC(name, staterrfct)
#define STATFCT(name)		_STATFCTSPEC(name, statfct)
#define WSTATFCT(name)		_WSTATFCTSPEC(name, wstatfct)


#endif /* __sherpa_stat_extension_hh__ */
