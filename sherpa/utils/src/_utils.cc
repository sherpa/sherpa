//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#include "sherpa/extension.hh"
#include "sherpa/utils.hh"
#include <cmath>
#include <vector>
#include <limits>
#include <iostream>
#include <sstream>

extern "C" {

#include "cephes.h"
#include "fcmp.h"

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


PyObject* sao_arange( PyObject* self, PyObject* args )
{
  
  double start, stop, step = 1.0;
  DoubleArray result;
  std::vector<double> arr;
  double eps = std::numeric_limits< double >::epsilon();

  if ( !PyArg_ParseTuple( args, (char*)"dd|d", &start, &stop, &step ) )
    return NULL;
  
  int count = 0;
  double bin = start;
  while( _sao_fcmp(bin, stop, eps) < 0 ) {
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

  // F-Test
  //FCTSPEC(calc_ftest, ftest),
  { (char*) "calc_ftest", (PyCFunction) ftest, METH_VARARGS,
    (char*) " Calculate the significance using the F test\n\nExample:\nsherpa> calc_ftest( <dof_1>, <stat_1>, <dof_2>, <stat_2> )"},
  
  // Maximum likelihood ratio
  //FCTSPEC(calc_mlr, mlr),
  { (char*) "calc_mlr", (PyCFunction) mlr, METH_VARARGS,
    (char*) " Calculate the Maximum Likelihood Ratio test\n\nExample:\nsherpa> calc_mlr( <delta_dof>, <delta_stat> )"},
  
  // Compliment of incomplete gamma function
  { (char*) "igamc", (PyCFunction) igamc, METH_VARARGS,
    (char*) " Calculate the compliment of the incomplete Gamma (upper) function [a>0; x>0]\n\nExample:\nsherpa> igamc( a, x )"},

  // Incomplete gamma function
  { (char*) "igam", (PyCFunction) igam, METH_VARARGS,
    (char*) " Calculate the incomplete Gamma (lower) function [a>0; x>0]\n\nExample:\nsherpa> igam( a, x )"},

  // Incomplete beta function
  
  { (char*) "incbet", (PyCFunction) incbet, METH_VARARGS,
    (char*) " Calculate the incomplete Beta function [a>0; b>0; 0<=x<=1]\n\nExample:\nsherpa> incbet( a, b, x )"},
  
  // Gamma function
  { (char*) "gamma", (PyCFunction) gamma, METH_VARARGS,
    (char*) " Calculate the Gamma function [-170<=x<=171]\n\nExample:\nsherpa> gamma( x )"},
  
  // Log gamma function
  { (char*) "lgam", (PyCFunction) lgam, METH_VARARGS,
    (char*) " Calculate the log of the Gamma function\n\nExample:\nsherpa> lgam( x )"},

  // Error Function
  FCTSPEC(erf, erf),

  // Inverse of Normal distribution Function
  FCTSPEC(ndtri, ndtri),

  // This function determines whether x and y are approximately equal
  // to a relative accuracy epsilon.
  FCTSPEC(gsl_fcmp, _sherpa_fcmp< _gsl_fcmp >), 

  //Same as gsl_fcmp, but also handles the case where one of the args is 0.
  FCTSPEC(sao_fcmp, _sherpa_fcmp< _sao_fcmp >),

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
