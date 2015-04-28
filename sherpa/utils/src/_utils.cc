// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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
#include "sherpa/utils.hh"
#include "sherpa/fcmp.hh"
#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <sstream>

extern "C" {

#include "cephes.h"
  //#include "fcmp.h"

  void init_utils();

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

  double f, x;

  for ( npy_intp ii = 0; ii < nelem; ii++ ) {

    f = ( chisq_1[ii] / dof_1[ii] ) / ( chisq_2[ii] / dof_2[ii] );
    x = dof_2[ii] / ( dof_2[ii] + dof_1[ii] * f );    
    result[ii] = incbet( dof_2[ii] / 2.0, dof_1[ii] / 2.0, x);

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

// for documentation
#define SEEALSODOC "\nSee also\n--------\n"
#define NOTESDOC "\nNotes\n-----\n"
#define REFERENCESDOC "\nReferences\n----------\n"
#define EXAMPLESDOC "\nExamples\n--------\n"
#define PARAMETERSDOC "\nParameters\n----------\n"
#define RETURNSDOC "\nReturns\n-------\n"

static PyMethodDef UtilsFcts[] = {

  // F-Test
  //FCTSPEC(calc_ftest, ftest),
  { (char*) "calc_ftest", (PyCFunction) ftest, METH_VARARGS,
    (char*) "calc_ftest(dof1, stat1, dof2, stat2)\n\n"
            "Compare two models using the F test.\n\n"
            "The F-test is a model comparison test; that is, it is a test\n"
            "used to select from two competing models which best describes\n"
            "a particular data set. A model comparison test statistic, T,\n"
            "is created from the best-fit statistics of each fit; as with all\n"
            "statistics, it is sampled from a probability distribution p(T).\n"
            "The test significance is defined as the integral of p(T) from the\n"
            "observed value of T to infinity. The significance quantifies the\n"
            "probability that one would select the more complex model when in\n"
            "fact the null hypothesis is correct. See also `calc_mlr`.\n"
            PARAMETERSDOC
            "dof1 : int\n"
            "   degrees of freedom of the simple model\n"
            "stat1 : number\n"
            "   best-fit chi-square statistic value of the simple model\n"
            "dof2 : int\n"
            "   degrees of freedom of the complex model\n"
            "stat2 : number\n"
            "   best-fit chi-square statistic value of the complex model\n"
            RETURNSDOC
            "sig : number\n"
            "   The significance, or p-value. A standard threshold for\n"
            "   selecting the more complex model is significance < 0.05 (the\n"
            "   '95% criterion' of statistics).\n"
            SEEALSODOC
            "calc_mlr : Calculate the maximum-likelihood ratio test.\n"
            "incbet  : Calculate the incomplete Beta function.\n"
            NOTESDOC
            "The F test uses the ratio of the reduced chi2, which follows\n"
            "the F-distribution, (stat1/dof1)/(stat2/dof2). The incomplete\n"
            "Beta function is used to calculate the integral of the tail of\n"
            "the F-distribution.\n\n"
            "The F test should only be used when:\n\n"
            " - the simpler of the two models is nested within the other;\n"
            "   that is, one can obtain the simpler model by setting the extra\n"
            "   parameters of the more complex model (often to zero or one);\n"
            " - those normal distributions are not truncated by parameter space\n"
            "   boundaries;\n"
            " - the best-fit statistics are sampled from the chi-square\n"
            "   distribution.\n"
            EXAMPLESDOC "\n"
            "In this example, the simple model has 5 degrees of freedom and\n"
            "a chi-square statistic of 7.73, while the complex model has 8\n"
            "degrees of freedom and a chi-square statistic of 8.94. The\n"
            "F test does not provide any evidence that the complex model\n"
            "is a better fit to the data than the simple model since the\n"
            "result is much larger than 0.\n\n"
            ">>> calc_ftest(5, 7.73, 8, 8.94)\n0.32480691622314933\n\n"},
  
  // Maximum likelihood ratio
  //FCTSPEC(calc_mlr, mlr),
  { (char*) "calc_mlr", (PyCFunction) mlr, METH_VARARGS,
    (char*) " Calculate the Maximum Likelihood Ratio test\n\nExample:\nsherpa> calc_mlr( <delta_dof>, <delta_stat> )"},
  
  // Complement of incomplete gamma function
  { (char*) "igamc", (PyCFunction) igamc, METH_VARARGS,
    (char*) "igamc(a,x)\n\n"
            "Calculate the complement of the regularized incomplete Gamma\n"
            "function (upper).\n\n"
            "The function is defined using the regularized incomplete Gamma\n"
            "function - igam(a,x) - and the Gamma function - gamma(a) - as\n\n"
            "igamc(a,x) = 1 - igam(a,x)\n"
            "           = 1/gamma(a) Int_x^Inf e^(-t) t^(a-1) dt\n"
            PARAMETERSDOC
            "a : scalar or array\n"
            "   a > 0\n"
            "x : scalar or array\n"
            "   x > 0\n"
            RETURNSDOC
            "val : scalar or array\n"
            SEEALSODOC
            "gamma : The Gamma function.\n"
            "igam  : The incomplete Gamma function (upper).\n"
            NOTESDOC
            "In this implementation, which is provided by the Cephes Math\n"
            "Library [1]_, both arguments must be positive. The integral is\n"
            "evaluated by either a power series or continued fraction expansion,\n"
            "depending on the relative values of a and x. Using IEEE arithmetic,\n"
            "the relative errors are\n\n"
            "========  ======  ========  =======  =======\n"
            " domain   domain  # trials   peak      rms\n"
            "========  ======  ========  =======  =======\n"
            "0.5,100   0,100   200000    1.9e-14  1.7e-15\n"
            "0.01,0.5  0,100   200000    1.4e-13  1.6e-15\n"
            "========  ======  ========  =======  =======\n"
            REFERENCESDOC "\n"
            ".. [1] Cephes Math Library Release 2.0:  April, 1987.\n"
            "       Copyright 1985, 1987 by Stephen L. Moshier.\n"
            "       Direct inquiries to 30 Frost Street, Cambridge, MA 02140.\n"
            EXAMPLESDOC "\n"
            ">>> igamc(1, 2)\n0.1353352832366127\n\n"
            ">>> igamc([1,1], [2,3])\narray([ 0.13533528,  0.04978707])\n\n"},

  // Incomplete gamma function
  { (char*) "igam", (PyCFunction) igam, METH_VARARGS,
    (char*) "igam(a,x)\n\n"
            "Calculate the regularized incomplete Gamma function (lower).\n\n"
            "The function is defined using the complete Gamma function -\n"
            "gamma(a) - as\n\n"
            "igam(a,x) = 1/gamma(a) Int_0^x e^(-t) t^(a^-1) dt\n"
            PARAMETERSDOC
            "a : scalar or array\n"
            "   a > 0\n"
            "x : scalar or array\n"
            "   x > 0\n"
            RETURNSDOC
            "val : scalar or array\n"
            SEEALSODOC
            "gamma : The Gamma function.\n"
            "igamc : The complement of the incomplete Gamma function (upper).\n"
            NOTESDOC
            "In this implementation, which is provided by the Cephes Math\n"
            "Library [1]_, both arguments must be positive. The integral is\n"
            "evaluated by either a power series or continued fraction expansion,\n"
            "depending on the relative values of a and x. Using IEEE arithmetic,\n"
            "the relative errors are\n\n"
            "======  ========  =======  =======\n"
            "domain  # trials   peak      rms\n"
            "======  ========  =======  =======\n"
            "0,30    200000    3.6e-14  2.9e-15\n"
            "0,100   300000    9.9e-14  1.5e-14\n"
            "======  ========  =======  =======\n"
            REFERENCESDOC "\n"
            ".. [1] Cephes Math Library Release 2.0:  April, 1987.\n"
            "       Copyright 1985, 1987 by Stephen L. Moshier.\n"
            "       Direct inquiries to 30 Frost Street, Cambridge, MA 02140.\n"
            EXAMPLESDOC "\n"
            ">>> igam(1, 2)\n0.8646647167633873\n\n"
            ">>> igam([1,1], [2,3])\narray([ 0.86466472,  0.95021293])\n\n"},

  // Incomplete beta function
  
  { (char*) "incbet", (PyCFunction) incbet, METH_VARARGS,
    (char*) "incbet(a,b,x)\n\n"
            "Calculate the incomplete Beta function\n\n"
            "The function is defined as\n"
            "sqrt(a+b)/(sqrt(a) sqrt(b)) Int_0^x t^(a-1) (1-t)^(b-1) dt\n"
            "and the integral from x to 1 can be obtained using the relation\n"
            "1 - incbet(a, b, x) = incbet(b, a, 1-x)\n"
            PARAMETERSDOC
            "a : scalar or array\n"
            "   a > 0\n"
            "b : scalar or array\n"
            "   b > 0\n"
            "x : scalar or array\n"
            "   0 <= x <= 1\n"
            RETURNSDOC
            "val : scalar or array\n"
            SEEALSODOC
            "calc_ftest : Compare two models using the F test.\n"
            NOTESDOC
            "In this implementation, which is provided by the Cephes Math\n"
            "Library [1]_, the integral is evaluated by a continued fraction\n"
            "expansion or, when b*x is small, by a power series.\n\n"
            "Using IEEE arithmetic, the relative errors are (tested uniformly\n"
            "distributed random points (a,b,x) with a and b in 'domain' and\n"
            "x between 0 and 1):\n\n"
            "========  ========  =======  =======\n"
            " domain   # trials   peak      rms\n"
            "========  ========  =======  =======\n"
            "0,5       10000     6.9e-15  4.5e-16\n"
            "0,85      250000    2.2e-13  1.7e-14\n" // DOC-TODO: is 25000 correct here?
            "0,1000    30000     5.3e-12  6.3e-13\n"
            "0,1000    250000    9.3e-11  7.1e-12\n" // DOC-TODO: is 25000 correct here?
            "0,100000  10000     8.7e-10  4.8e-11\n"
            "========  ========  =======  =======\n\n"
            "Outputs smaller than the IEEE gradual underflow threshold were\n"
            "excluded from these statistics.\n"
            REFERENCESDOC "\n"
            ".. [1] Cephes Math Library Release 2.0:  April, 1987.\n"
            "       Copyright 1985, 1987 by Stephen L. Moshier.\n"
            "       Direct inquiries to 30 Frost Street, Cambridge, MA 02140.\n"
            EXAMPLESDOC "\n"
            ">>> incbet(0.3, 0.6, 0.5)\n0.68786273145845922\n\n"
            ">>> incbet([0.3,0.3], [0.6,0.7], [0.5,0.4])\n"
            "array([ 0.68786273,  0.67356524])\n\n"},
  
  // Gamma function
  { (char*) "gamma", (PyCFunction) gamma, METH_VARARGS,
    (char*) "gamma(z)\n\n"
            "Calculate the Gamma function.\n"
            PARAMETERSDOC
            "z : scalar or array\n"
            "   -171 <= z <- 171.6\n"
            RETURNSDOC
            "val : scalar or array\n"
            SEEALSODOC
            "lgam : The log of the Gamma function.\n"
            "igam : The incomplete Gamma function.\n"
            NOTESDOC
            "This implementation is provided by the Cephes Math Library [1]_.\n"
            "Arguments |x| >= 34 are reduced by recurrence and the function\n"
            "approximated by a rational function of degree 6/7 in the interval\n"
            "(2,3). Large arguments are handled by Stirling's formula. Large\n"
            "negative arguments are made positive using a reflection formula.\n"
            "Relative errors are\n\n"
            "========  ========  =======  =======\n"
            " domain   # trials   peak      rms\n"
            "========  ========  =======  =======\n"
            "-170,33   20000     2.3e-15  3.3e-16\n"
            "-33,33    20000     9.4e-16  2.2e-16\n"
            "33,171.6  20000     2.3e-15  3.2e-16\n"
            "========  ========  =======  =======\n\n"
            "Errors for arguments outside the test range will be larger owing\n"
            "to amplification by the exponential function.\n"
            REFERENCESDOC "\n"
            ".. [1] Cephes Math Library Release 2.0:  April, 1987.\n"
            "       Copyright 1985, 1987 by Stephen L. Moshier.\n"
            "       Direct inquiries to 30 Frost Street, Cambridge, MA 02140.\n"
            EXAMPLESDOC "\n"
            ">>> gamma(2.3)\n1.1667119051981603\n\n"
            ">>> gamma([2.3,1.9])\narray([ 1.16671191,  0.96176583])\n\n"},
  
  // Log gamma function
  { (char*) "lgam", (PyCFunction) lgam, METH_VARARGS,
    (char*) "lgam(z)\n\n"
            "Calculate the log (base e) of the Gamma function.\n"
            PARAMETERSDOC
            "z : scalar or array\n"
            "   0 <= z <= 2.556348e305\n"
            RETURNSDOC
            "val : scalar or array\n"
            SEEALSODOC
            "gamma : The Gamma function.\n"
            "igam : The incomplete Gamma function.\n"
            NOTESDOC
            "This implementation is provided by the Cephes Math Library [1]_.\n"
            "For arguments greater than 13, the logarithm of the Gamma function\n"
            "is approximated by the logarithmic version of Stirling's formula\n"
            "using a polynomial approximation of degree 4. Arguments\n"
            "between -33 and +33 are reduced by recurrence to the interval [2,3]\n"
            "of a rational approximation. The cosecant reflection formula is\n"
            "employed for arguments less than -33.\n\n"
            "Relative errors are\n\n"
            "===============  ========  =======  =======\n"
            "    domain       # trials   peak      rms\n"
            "===============  ========  =======  =======\n"
            "0,3              28000     5.4e-16  1.1e-16\n"
            "2.718,2.556e305  40000     3.5e-16  8.3e-17\n"
            "===============  ========  =======  =======\n\n"
            "The error criterion was relative when the function magnitude was\n"
            "greater than one but absolute when it was less than one.\n\n"
            "The following test used the relative error criterion, though at\n"
            "certain points the relative error could be much higher than\n"
            "indicated.\n\n"
            "=======  ========  =======  =======\n"
            "domain   # trials   peak      rms\n"
            "=======  ========  =======  =======\n"
            "-200,-4  10000     4.8e-16  1.3e-16\n"
            "=======  ========  =======  =======\n"
            REFERENCESDOC "\n"
            ".. [1] Cephes Math Library Release 2.0:  April, 1987.\n"
            "       Copyright 1985, 1987 by Stephen L. Moshier.\n"
            "       Direct inquiries to 30 Frost Street, Cambridge, MA 02140.\n"
            EXAMPLESDOC "\n"
            ">>> lgam(104.56)\n380.21387239435785\n\n"
            ">>> lgam([104.56,2823.4])\narray([   380.21387239,  19607.42734396])\n\n"},

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
  { (char*) "rebin", (PyCFunction)(rebin<SherpaFloatArray, SherpaFloat>), METH_VARARGS,
    (char*) "integrate from: (x0lo, x0hi, x0) to: (x1lo, x1hi) yielding => ty\n\nExample:\nsherpa> rebin( x0, x0lo, x0hi, x1lo, x1hi )"},

  //Histogram1d 
  { (char*) "hist1d",(PyCFunction)(histogram1d<SherpaFloatArray, SherpaFloat, IntArray, int>), METH_VARARGS, (char*) " create 1D histogram\n\nExample:\nsherpa> histogram1d( x, xlo, xhi )"},
  
  //Histogram2d 
  { (char*) "hist2d",(PyCFunction)(histogram2d<SherpaFloatArray, SherpaFloat, IntArray, int>), METH_VARARGS, (char*) " create 2D histogram\n\nExample:\nsherpa> histogram2d( x, y, x_grid, y_grid )"},

  FCTSPEC(sum_intervals, (sum_intervals<SherpaFloatArray, IntArray,
			  SherpaFloat, int, npy_intp>)),

  //neville
  { (char*) "neville",(PyCFunction)(neville<SherpaFloatArray, SherpaFloat>), METH_VARARGS, (char*) " neville interpolation\n\nExample:\nsherpa> yout = neville(xout, xin, yin)"},

  FCTSPEC(sao_arange, sao_arange),
  
  { NULL, NULL, 0, NULL }

};

SHERPAMOD(_utils, UtilsFcts)
