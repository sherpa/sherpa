// 
//  Copyright (C) 2007, 2015, 2016, 2018, 2019  Smithsonian Astrophysical Observatory
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

#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <sstream>

#include "sherpa/extension.hh"
#include "sherpa/utils.hh"
#include "sherpa/fcmp.hh"
#include "Faddeeva.hh"

extern "C" {

#include "cephes.h"
  //#include "fcmp.h"

  void init_utils();

}

static PyObject* wofz( PyObject* self, PyObject* args )
{
  
  ComplexArray xxx;
  if( !PyArg_ParseTuple( args, (char*)"O&",
                         CONVERTME(ComplexArray), &xxx) )
    return NULL;

  ComplexArray result;
  if ( EXIT_SUCCESS != result.create( xxx.get_ndim(), xxx.get_dims() ) )
    return NULL;

  npy_intp nelem = xxx.get_size();
  for (npy_intp ii = 0; ii < nelem; ++ii)
    result[ii] = Faddeeva::w(xxx[ii]);

  return result.return_new_ref();

}

static PyObject* ftest( PyObject* self, PyObject* args )
{
  
  DoubleArray dof_1;
  DoubleArray dof_2;
  DoubleArray chisq_1;
  DoubleArray chisq_2;

  if( !PyArg_ParseTuple( args, (char*)"O&O&O&O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &dof_1,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &chisq_1,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &dof_2,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &chisq_2 ) )
    return NULL;

  npy_intp nelem = dof_1.get_size();

  if ( dof_2.get_size() != nelem ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dof_1: " << nelem << " vs dof_2: " << dof_2.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( chisq_1.get_size() != nelem ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dof_1: " << nelem << " vs chisq_1: " << chisq_1.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( chisq_2.get_size() != nelem ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "dof_1: " << nelem << " vs chisq_2: " << chisq_2.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( dof_1.get_ndim(), dof_1.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < nelem; ii++ ) {

    const double delta_dof = dof_1[ii] - dof_2[ii];
    if ( 0.0 == chisq_2[ii] ) {
      std::ostringstream err;
      err << "chisq_2[" << ii << "] cannot be equal to 0: ";
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    if ( 0.0 == dof_2[ii] ) {
      std::ostringstream err;
      err << "dof_2[" << ii << "] cannot be equal to 0: ";
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    if ( 0.0 == delta_dof ) {
      std::ostringstream err;
      err << "dof_1[" << ii << "] cannot be equal to dof_2[" << ii << "]: ";
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    const double delta_chi = chisq_1[ii] - chisq_2[ii];
    const double f = delta_chi / delta_dof / (chisq_2[ii] / dof_2[ii]);
    const double tmp = dof_2[ii] + delta_dof * f;
    if ( 0.0 == tmp ) {
      std::ostringstream err;
      err << "dof_2[" << ii << "] + delta_dof * f cannot be equal to 0,\nwhere f = delta_chi / delta_dof / (chisq_2[" << ii << "] / dof_2[" << ii << " ]: ";
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }
    result[ii] = incbet( dof_2[ii] * 0.5 , delta_dof * 0.5,
                         (dof_2[ii] / ( tmp ) ) );

  }
  
  return result.return_new_ref();

}


static PyObject* mlr( PyObject* self, PyObject* args )
{

  DoubleArray delta_dof;
  DoubleArray delta_chisq;

  if( !PyArg_ParseTuple( args, (char*)"O&O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &delta_dof,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &delta_chisq ) )
    return NULL;
  
  npy_intp nelem = delta_dof.get_size();

  if ( delta_chisq.get_size() != nelem  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "delta_dof: " << nelem
	<< " vs delta_chisq: " << delta_chisq.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( delta_dof.get_ndim(),
				      delta_dof.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < nelem; ii++ )
    result[ii] = igamc( delta_dof[ii] / 2.0, delta_chisq[ii] / 2.0 );  
 
  return result.return_new_ref();

}

extern int sgngam;

static PyObject* gamma( PyObject* self, PyObject* args )
{
  DoubleArray x;

  if( !PyArg_ParseTuple( args, (char*)"O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x ) )
    return NULL;

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( x.get_ndim(), x.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < x.get_size(); ii++ ) {
    result[ii] = Gamma(x[ii]);
    result[ii] *= sgngam;
  }
  
  return result.return_new_ref();

}

static PyObject* igam( PyObject* self, PyObject* args )
{

  DoubleArray a;
  DoubleArray x;

  if( !PyArg_ParseTuple( args, (char*)"O&O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &a,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x ) )
    return NULL;
  
  npy_intp asize = a.get_size();
  npy_intp xsize = x.get_size();

  if ( asize != xsize  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "a: " << asize << " vs x: " << xsize;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( x.get_ndim(),
				      x.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < xsize; ii++ ) {

    if( a[ii] < 0 || x[ii] < 0 ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"igam domain error, a and x must be positive");
      return NULL;
    }

    result[ii] = igam( a[ii], x[ii]);  
 
  }
  return result.return_new_ref();

}

static PyObject* igamc( PyObject* self, PyObject* args )
{

  DoubleArray a;
  DoubleArray x;

  if( !PyArg_ParseTuple( args, (char*)"O&O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &a,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x ) )
    return NULL;
  
  npy_intp asize = a.get_size();
  npy_intp xsize = x.get_size();

  if ( asize != xsize  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "a: " << asize << " vs x: " << xsize;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( x.get_ndim(),
				      x.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < xsize; ii++ ) {

    if( a[ii] < 0 || x[ii] < 0 ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"igamc domain error, a and x must be positive");
      return NULL;
    }

    result[ii] = igamc( a[ii], x[ii]);  
 
  }
  return result.return_new_ref();

}

static PyObject* lgam( PyObject* self, PyObject* args )
{
  DoubleArray x;

  if( !PyArg_ParseTuple( args, (char*)"O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x ) )
    return NULL;

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( x.get_ndim(),
				      x.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < x.get_size(); ii++ ) {
    result[ii] = lgam(x[ii]);
    result[ii] *= sgngam;
  }

  return result.return_new_ref();

}

static PyObject* incbet( PyObject* self, PyObject* args )
{
  DoubleArray x;
  DoubleArray a;
  DoubleArray b;

  if( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &a,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &b,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x) )
    return NULL;

  npy_intp xsize = x.get_size();
  npy_intp asize = a.get_size();
  npy_intp bsize = b.get_size();

  if ( asize != xsize  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "a: " << asize << " vs x: " << xsize;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( asize != bsize  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "a: " << asize << " vs b: " << bsize;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( x.get_ndim(),
				      x.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < xsize; ii++ ) {

    if( x[ii] < 0 || x[ii] > 1 ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"incbeta domain error, 0 <= x <= 1" );
      return NULL;
    }

    if( a[ii] < 0 || b[ii] < 0 ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"incbeta domain error, a and b must be positive");      return NULL;
    }

    result[ii] = incbet( a[ii], b[ii], x[ii]);
  }

  return result.return_new_ref();

}

static PyObject* erf( PyObject* self, PyObject* args )
{

  DoubleArray a;

  if( !PyArg_ParseTuple( args, (char*)"O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &a ) )
    return NULL;

  npy_intp nelem = a.get_size();

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( a.get_ndim(), a.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < nelem; ii++ )
    result[ii] = erf( a[ii] );

  return result.return_new_ref();

}

static PyObject* ndtri( PyObject* self, PyObject* args )
{

  DoubleArray a;

  if( !PyArg_ParseTuple( args, (char*)"O&",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &a ) )
    return NULL;

  npy_intp nelem = a.get_size();

  DoubleArray result;
  if ( EXIT_SUCCESS != result.create( a.get_ndim(), a.get_dims() ) )
    return NULL;

  for ( npy_intp ii = 0; ii < nelem; ii++ )
    result[ii] = ndtri( a[ii] );

  return result.return_new_ref();

}

template <int (*fcmp)( const double x1, const double x2,
		       const double epsilon )>
PyObject* _sherpa_fcmp( PyObject *self, PyObject *args )
{

  DoubleArray x1;
  DoubleArray x2;
  double epsilon;

  if( !PyArg_ParseTuple( args, (char*)"O&O&d",
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x1,
			 (converter)sherpa::convert_to_array< DoubleArray >,
			 &x2,
			 &epsilon ) )
    return NULL;
  
  npy_intp n1 = x1.get_size();
  npy_intp n2 = x2.get_size();

  if ( n1 > 1 && n1 != n2 ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x2: " << n2 << " vs x1: " << n1;
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  IntArray result;
  if ( EXIT_SUCCESS != result.create( x2.get_ndim(), x2.get_dims() ) )
    return NULL;

  // compare x1_0 to x2_i
  if (n1 == 1)
    for ( npy_intp ii = 0; ii < n2; ii++ ) {
      double val = x1[ 0 ];
      result[ ii ] = fcmp( val, x2[ ii ], epsilon );
    }
  // compare x1_i to x2_i
  else
    for ( npy_intp ii = 0; ii < n2; ii++ )
      result[ ii ] = fcmp( x1[ ii ], x2[ ii ], epsilon );

  return result.return_new_ref();

}

template <typename ArrayType, typename DataType>
PyObject* rebin( PyObject* self, PyObject* args )
{
  
  ArrayType x0;
  ArrayType x0lo;
  ArrayType x0hi;
  ArrayType x1;
  ArrayType x1lo;
  ArrayType x1hi;
  
  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&",
			  (converter)sherpa::convert_to_array<ArrayType>,
			  &x0,
			  (converter)sherpa::convert_to_array<ArrayType>,
			  &x0lo,
			  (converter)sherpa::convert_to_array<ArrayType>,
			  &x0hi,
			  (converter)sherpa::convert_to_array<ArrayType>,
			  &x1lo,
			  (converter)sherpa::convert_to_array<ArrayType>,
			  &x1hi ))
    
    return NULL;
  
  if ( x0.get_size() != x0lo.get_size()  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x0: " << x0.get_size() << " vs x0lo: " << x0lo.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( x0hi.get_size() != x0lo.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x0hi: " << x0hi.get_size() << " vs x0lo: " << x0lo.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }	 

  if ( x1hi.get_size() != x1lo.get_size() ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x1hi: " << x1hi.get_size() << " vs x1lo: " << x1lo.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( EXIT_SUCCESS != x1.create( x1lo.get_ndim(), x1lo.get_dims() ) )
    return NULL;

  if ( EXIT_SUCCESS != (sherpa::utils::rebin_histogram <DataType, ArrayType>
			(
			 x0, x0lo, x0hi,
			 x0.get_size(),
			 x1, x1lo, x1hi,
			 x1lo.get_size() ) ) ){
    PyErr_SetString( PyExc_ValueError, (char*)"rebinning data failed" );
    return NULL;
  }
  
  return x1.return_new_ref();
  
}

template <typename ArrayType, typename DataType,
	  typename IndexArrayType, typename IndexType>
PyObject* histogram1d( PyObject* self, PyObject* args )
{
  
  ArrayType x;
  ArrayType x_lo;
  ArrayType x_hi;
  IndexArrayType res;

  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			  CONVERTME(ArrayType),
			  &x,
			  CONVERTME(ArrayType),
			  &x_lo,
			  CONVERTME(ArrayType),
			  &x_hi ))
    return NULL;

  if ( x_lo.get_size() != x_hi.get_size()  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x_lo: " << x_lo.get_size() << " vs x_hi: " << x_hi.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }
  
  if ( (x.get_size() < 1) || (x_lo.get_size() < 1) || (x_hi.get_size() < 1) ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"need at least one element for histogram");
    return NULL;
  }

  if ( EXIT_SUCCESS != res.zeros(  x_lo.get_ndim(), x_lo.get_dims() ) )
    return NULL;

  if ( EXIT_SUCCESS != (sherpa::utils::histogram1d
			<ArrayType, DataType, IndexArrayType, IndexType>
			( x, x_lo, x_hi, res ) ) ) {
    PyErr_SetString( PyExc_ValueError, (char*)"histogram1d failed" );
    return NULL;
  }
  
  return res.return_new_ref();
}

template <typename ArrayType, typename DataType,
	  typename IndexArrayType, typename IndexType>
PyObject* histogram2d( PyObject* self, PyObject* args )
{
  
  ArrayType x;
  ArrayType y;
  ArrayType x_grid;
  ArrayType y_grid;
  IndexArrayType res;

  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&",
			  CONVERTME(ArrayType),
			  &x,
			  CONVERTME(ArrayType),
			  &y,
			  CONVERTME(ArrayType),
			  &x_grid,
			  CONVERTME(ArrayType),
			  &y_grid))
    return NULL;
  
  if ( x.get_size() != y.get_size()  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "x: " << x.get_size() << " vs y: " << y.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( x.get_size() < 1 ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"need at least one element for histogram");
    return NULL;
  }

  npy_intp dims[1];
  dims[0] = npy_intp( x_grid.get_size() * y_grid.get_size() );

  if ( EXIT_SUCCESS != res.zeros( 1, dims ) )
    return NULL;

  if ( EXIT_SUCCESS != (sherpa::utils::histogram2d
			<ArrayType, DataType, IndexArrayType, IndexType>
			( x, y, x_grid, y_grid, res ) ) ) {
    PyErr_SetString( PyExc_ValueError, (char*)"histogram2d failed" );
    return NULL;
  }
  
  return res.return_new_ref();
}


template <typename FloatArrayType, typename IntArrayType,
	  typename DataType, typename IntType, typename IndexType>
PyObject* sum_intervals( PyObject* self, PyObject* args )
{
  
  FloatArrayType src;
  FloatArrayType model;
  IntArrayType indx0;
  IntArrayType indx1;

  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			  CONVERTME(FloatArrayType),
			  &src,
			  CONVERTME(IntArrayType),
			  &indx0,
			  CONVERTME(IntArrayType),
			  &indx1))
    return NULL;
  
  if ( indx0.get_size() != indx1.get_size()  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "indx0: " << indx0.get_size() << " vs indx1: " << indx1.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( EXIT_SUCCESS != model.zeros( indx0.get_ndim(), indx0.get_dims() ) )
    return NULL;
  
  if ( EXIT_SUCCESS != (sherpa::utils::sum_intervals<DataType,
			DataType, IntType, IndexType>
			(&src[0], &indx0[0], &indx1[0], model.get_size(),
			 &model[0] ) ) ){
    PyErr_SetString( PyExc_ValueError, (char*)"sum_intervals" );
    return NULL;
  }
  
  return model.return_new_ref();
}

template <typename ArrayType, typename DataType>
PyObject* neville( PyObject* self, PyObject* args )
{
  
  ArrayType xout;
  ArrayType xin;
  ArrayType yin;
  ArrayType result;

  if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			  CONVERTME(ArrayType),
			  &xout,
			  CONVERTME(ArrayType),
			  &xin,
			  CONVERTME(ArrayType),
			  &yin))
    return NULL;
  
  if ( xin.get_size() != yin.get_size()  ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "xin: " << xin.get_size() << " vs yin: " << yin.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( EXIT_SUCCESS != result.zeros( xout.get_ndim(), xout.get_dims() ) )
    return NULL;
  
  int nin = (int) xin.get_size();
  int nout = (int) xout.get_size();

  for(int ii = 0; ii < nout; ii++) {

    if ( EXIT_SUCCESS != (sherpa::utils::neville<ArrayType, DataType>
			  (nin, xin, yin, xout[ii], result[ii] ) ) ){

      PyErr_SetString( PyExc_ValueError,
		       (char*)"neville interpolation failed" );

      return NULL;
    }
    
  }
  
  return result.return_new_ref();
}


static PyObject* sao_arange( PyObject* self, PyObject* args )
{
  
  double start, stop, step = 1.0;
  DoubleArray result;
  std::vector<double> arr;
  double eps = std::numeric_limits< double >::epsilon();

  if ( !PyArg_ParseTuple( args, (char*)"dd|d", &start, &stop, &step ) )
    return NULL;
  
  int count = 0;
  double bin = start;
  while( sao_fcmp(bin, stop, eps) < 0 ) {
    bin = start + double(count)*step;
    arr.push_back(bin);
    ++count;
  }
  
  npy_intp dim;
  dim = (npy_intp)arr.size();

  if ( EXIT_SUCCESS != result.create( 1, &dim ) )
    return NULL;
  
  for(npy_intp ii = 0; ii < dim; ++ii)
    result[ii] = arr[ii];
   
  return result.return_new_ref();
}

static PyMethodDef UtilsFcts[] = {

  // wofz
  FCTSPEC(calc_wofz, wofz),

  // F-Test
  FCTSPEC(calc_ftest, ftest),
  
  // Maximum likelihood ratio
  FCTSPEC(calc_mlr, mlr),

  // Complement of incomplete gamma function
  FCTSPEC(igamc, igamc),

  // Incomplete gamma function
  FCTSPEC(igam, igam),

  // Incomplete beta function
  FCTSPEC(incbet, incbet),
  
  // Gamma function
  FCTSPEC(gamma, gamma),
  
  // Log gamma function
  FCTSPEC(lgam, lgam),

  // Error Function
  FCTSPEC(erf, erf),

  // Inverse of Normal distribution Function
  FCTSPEC(ndtri, ndtri),

  // This function determines whether x and y are approximately equal
  // to a relative accuracy epsilon.
  FCTSPEC(gsl_fcmp, _sherpa_fcmp< gsl_fcmp >), 

  //Same as gsl_fcmp, but also handles the case where one of the args is 0.
  FCTSPEC(sao_fcmp, _sherpa_fcmp< sao_fcmp >),

  //Rebin to new grid
  FCTSPEC(rebin, (rebin<SherpaFloatArray, SherpaFloat>)),

  //Histogram1d 
  { (char*) "hist1d",(PyCFunction)(histogram1d<SherpaFloatArray, SherpaFloat, IntArray, int>), METH_VARARGS, (char*) " create 1D histogram\n\nExample:\nsherpa> histogram1d( x, xlo, xhi )"},
  
  //Histogram2d 
  { (char*) "hist2d",(PyCFunction)(histogram2d<SherpaFloatArray, SherpaFloat, IntArray, int>), METH_VARARGS, (char*) " create 2D histogram\n\nExample:\nsherpa> histogram2d( x, y, x_grid, y_grid )"},

  FCTSPEC(sum_intervals, (sum_intervals<SherpaFloatArray, IntArray,
			  SherpaFloat, int, npy_intp>)),

  //neville
  FCTSPEC(neville, (neville<SherpaFloatArray, SherpaFloat>)),

  FCTSPEC(sao_arange, sao_arange),
  
  { NULL, NULL, 0, NULL }

};

SHERPAMOD(_utils, UtilsFcts)
