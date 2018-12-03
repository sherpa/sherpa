//
//  Copyright (C) 2007, 2018  Smithsonian Astrophysical Observatory
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
#include <sherpa/functor.hh>

#include <iostream>
#include <stdexcept>

#include "DifEvo.hh"
#include "minim.hh"
#include "Opt.hh"
#include "NelderMead.hh"

#include "minpack/LevMar.hh"


//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//
static void lmdif_callback_fcn( int mfct, int npar, double* xpars,
                             double* fvec, int& ierr, PyObject* py_fcn ) {

  DoubleArray pars_array;
  npy_intp dims[1];

  dims[0] = npar;
  if ( EXIT_SUCCESS != pars_array.create( 1, dims, xpars ) ) {
    ierr = EXIT_FAILURE;
    return;
  }

  PyObject* rv = PyObject_CallFunction( py_fcn, (char*)"N",
                                        pars_array.new_ref() );
  if ( NULL == rv ) {
    ierr = EXIT_FAILURE;
    return;
  }

  DoubleArray vals_array;
  int stat = vals_array.from_obj( rv );
  Py_DECREF( rv );
  if ( EXIT_SUCCESS != stat ) {
    ierr = EXIT_FAILURE;
    return;
  }

  if ( vals_array.get_size() != mfct ) {
    PyErr_SetString( PyExc_TypeError,
		     "callback function returned wrong number of values" );
    ierr = EXIT_FAILURE;
    return;
  }

  std::copy( &vals_array[0], &vals_array[0] + mfct, fvec );

  return;

}

//*****************************************************************************
//
// py_cpp_lmdif:  Python wrapper function for C++ function lmdif
//
//*****************************************************************************
template< typename Func >
static PyObject* py_cpp_lmdif( PyObject* self, PyObject* args, Func func ) {

  PyObject* py_function=NULL;
  DoubleArray par, lb, ub, fjac;
  int mfct, maxnfev, nfev, info, verbose;
  double fval, ftol, xtol, gtol, epsfcn, factor;

  if ( !PyArg_ParseTuple( args, (char*) "OiO&dddiddiO&O&O&",
			  &py_function,
			  &mfct,
			  CONVERTME(DoubleArray), &par,
			  &ftol, &xtol, &gtol, &maxnfev,
			  &epsfcn, &factor, &verbose,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &fjac
                          ) ) {
    return NULL;
  }

  const int npar = par.get_size( );
  const int mn = mfct * npar;
  std::vector<double> covarerr( npar );
  std::vector<double> jacobian( mn );

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  if ( mn != fjac.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(fjac)=%d != m * n =%d",
		  static_cast<int>( ub.get_size( ) ), mn );
    return NULL;
  }


  try {

    minpack::LevMarDif< Func, PyObject*, double > levmar( func, py_function,
                                                          mfct );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    sherpa::Bounds<double> bounds( mylb, myub );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    info = levmar( npar, ftol, xtol, gtol, maxnfev, epsfcn, factor, verbose,
		   mypar, nfev, fval, bounds, jacobian, covarerr );

    // info > 0 means par needs to be updated
    if ( info > 0 )
      for ( int ii = 0; ii < npar; ++ii )
        par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) "Unknown exception caught" );
    return NULL;
  }

  // info == 0 is a no-op.

  if ( info < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) "function call failed" );
    return NULL;
  }

  std::copy( &jacobian[0], &jacobian[0] + ( mn ), &fjac[0] );
  return Py_BuildValue( (char*)"(NdiiN)", par.return_new_ref(), fval, nfev,
			info, fjac.return_new_ref() );
}
static PyObject* py_lmdif( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_cpp_lmdif( self, args, sherpa::fct_ptr( lmdif_callback_fcn ) );

}
//*****************************************************************************
//
// py_cpp_lmdif:  Python wrapper function for C++ function lmdif
//
//*****************************************************************************
//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//lmdif//

//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//
static void lmder_callback_fcn( int mfct, int npar, double* xpars,
                                double* fvec, int& iflag, PyObject* py_fcn ) {

  DoubleArray pars_array;
  npy_intp dims[1];

  dims[0] = npar;
  if ( EXIT_SUCCESS != pars_array.create( 1, dims, xpars ) ) {
    iflag = EXIT_FAILURE;
    return;
  }

  PyObject* rv = PyObject_CallFunction( py_fcn, (char*)"Ni",
                                        pars_array.new_ref(),
                                        iflag );
  if ( NULL == rv ) {
    iflag = -1;
    return;
  }

  DoubleArray vals_array;
  int stat = vals_array.from_obj( rv );
  Py_DECREF( rv );
  if ( EXIT_SUCCESS != stat ) {
    iflag = -1;
    return;
  }

  const int num = iflag == 1 ? mfct : mfct * npar;

  if ( vals_array.get_size() != num ) {
    PyErr_SetString( PyExc_TypeError,
		     "callback function returned wrong number of values" );
    iflag = -1;
    return;
  }

  std::copy( &vals_array[0], &vals_array[0] + num, fvec );

  return;

}

//*****************************************************************************
//
// py_cpp_lmder:  Python wrapper function for C++ function lmder
//
//*****************************************************************************
template< typename Func >
static PyObject* py_cpp_lmder( PyObject* self, PyObject* args, Func func ) {

  PyObject* py_function=NULL;
  DoubleArray par, lb, ub, fjac;
  int mfct, maxnfev, nfev, njev, info, rank, verbose;
  double fval, ftol, xtol, gtol, factor;

  if ( !PyArg_ParseTuple( args, (char*) "OiO&dddidiO&O&O&",
			  &py_function,
			  &mfct,
			  CONVERTME(DoubleArray), &par,
			  &ftol, &xtol, &gtol, &maxnfev,
			  &factor, &verbose,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &fjac
                          ) ) {
    return NULL;
  }

  const int npar = par.get_size( );
  const int mn = mfct * npar;
  std::vector<double> covarerr( npar );
  std::vector<double> jacobian( mn );

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  if ( mn != fjac.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(fjac)=%d != m * n =%d",
		  static_cast<int>( ub.get_size( ) ), mn );
    return NULL;
  }


  try {

    minpack::LevMarDer< Func, PyObject*, double > levmar( func, py_function,
                                                          mfct );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    sherpa::Bounds<double> bounds( mylb, myub );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    info = levmar( npar, ftol, xtol, gtol, maxnfev, factor, verbose,
		   mypar, nfev, njev, fval, jacobian, bounds, rank );
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) "Unknown exception caught" );
    return NULL;
  }

  if ( info < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) "function call failed" );
    return NULL;
  }

  std::copy( &jacobian[0], &jacobian[0] + ( mn ), &fjac[0] );
  return Py_BuildValue( (char*)"(NdiiiN)", par.return_new_ref(), fval, nfev,
                        njev, info, fjac.return_new_ref() );
}
static PyObject* py_lmder( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_cpp_lmder( self, args, sherpa::fct_ptr( lmder_callback_fcn ) );

}
//*****************************************************************************
//
// py_cpp_lmder:  Python wrapper function for C++ function lmder
//
//*****************************************************************************
//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//lmder//


//*****************************************************************************
//
// py_difevo_lm:  Python wrapper function for C++ function difevo
//
//*****************************************************************************
template< typename Func >
static PyObject* py_difevo_levmar( PyObject* self, PyObject* args,
				   Func callback_func ) {

  PyObject* py_function=NULL;
  DoubleArray par, step, lb, ub;
  int verbose, maxnfev, seed, population_size, mfcts, nfev, ierr;
  double fval, tol, xprob, weighting_factor;

  if ( !PyArg_ParseTuple( args, (char*) "iiiidddO&O&O&Oi",
			  &verbose,
			  &maxnfev,
			  &seed,
			  &population_size,
			  &tol,
			  &xprob,
			  &weighting_factor,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &par,
			  &py_function, &mfcts ) ) {
    return NULL;
  }

  const int npar = par.get_size( );

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  try {

    sherpa::DifEvo< Func, PyObject*, minpack::LevMarDif< Func, PyObject*, double >, double >
      difevo( callback_func, py_function, mfcts );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    ierr = difevo( verbose, maxnfev, tol, population_size, seed, xprob,
		  weighting_factor, npar, mylb, myub, mypar, nfev, fval );
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {

    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"Unknown exception caught" );
    return NULL;
  }


  if ( ierr < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"function call failed" );
    return NULL;
  }

  return Py_BuildValue( (char*)"(Ndii)", par.return_new_ref(), fval, nfev,
			ierr );
}
static PyObject* py_difevo_lm( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_difevo_levmar( self, args, sherpa::fct_ptr( lmdif_callback_fcn ) );

}
//*****************************************************************************
//
// py_difevo_lm:  Python wrapper function for C++ function difevo
//
//*****************************************************************************


static void sao_callback_func( int npar, double* xpars, double& fval,
			       int& ierr, PyObject* py_function ) {

  DoubleArray py_xpars;
  npy_intp dim[1];

  dim[0] = npar;
  if ( EXIT_SUCCESS != py_xpars.create( 1, dim, xpars ) ) {
    ierr = EXIT_FAILURE;
    return;
  }

  PyObject* return_val = PyObject_CallFunction( py_function, (char*)"O",
						py_xpars.borrowed_ref() );

  if ( NULL == return_val || Py_None == return_val ) {
    ierr = EXIT_FAILURE;
    return;
  }
  if ( !PyFloat_Check( return_val ) ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"callback did not return a float" );
      Py_DECREF( return_val );
      ierr = -1;
      return;
  }

  fval = PyFloat_AsDouble( return_val );

  Py_DECREF( return_val );

}

//*****************************************************************************
//
// py_difevo_nm:  Python wrapper function for C++ function difevo
//
//*****************************************************************************
template< typename Func >
static PyObject* py_difevo_neldermead( PyObject* self, PyObject* args,
				       Func func ) {

  PyObject* py_function=NULL;
  DoubleArray par, step, lb, ub;
  int verbose, maxnfev, seed, population_size, nfev, ierr;
  double fval, tol, xprob, weighting_factor;

  if ( !PyArg_ParseTuple( args, (char*) "iiiidddO&O&O&O",
			  &verbose,
			  &maxnfev,
			  &seed,
			  &population_size,
			  &tol,
			  &xprob,
			  &weighting_factor,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &par,
			  &py_function ) ) {
    return NULL;
  }

  const int npar = par.get_size( );

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  try {

    sherpa::DifEvo< Func, PyObject*, sherpa::NelderMead< Func, PyObject*, double >, double > difevo( func, py_function );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    ierr = difevo( verbose, maxnfev, tol, population_size, seed, xprob,
		   weighting_factor, npar, mylb, myub, mypar, nfev, fval );
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"Unknown exception caught" );
    return NULL;
  }

  if ( ierr < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"function call failed" );
    return NULL;
  }

  return Py_BuildValue( (char*)"(Ndii)", par.return_new_ref(), fval, nfev,
			ierr );
}
static PyObject* py_difevo_nm( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_difevo_neldermead( self, args,
			       sherpa::fct_ptr( sao_callback_func ) );

}
//*****************************************************************************
//
// py_difevo_nm:  Python wrapper function for C++ function difevo
//
//*****************************************************************************

//*****************************************************************************
//
// py_difevo:  Python wrapper function for C++ function difevo
//
//*****************************************************************************
template< typename Func >
static PyObject* py_difevo( PyObject* self, PyObject* args, Func func ) {

  PyObject* py_function=NULL;
  DoubleArray par, step, lb, ub;
  int verbose, maxnfev, seed, population_size, nfev, ierr;
  double fval, tol, xprob, weighting_factor;

  if ( !PyArg_ParseTuple( args, (char*) "iiiidddO&O&O&O",
			  &verbose,
			  &maxnfev,
			  &seed,
			  &population_size,
			  &tol,
			  &xprob,
			  &weighting_factor,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &par,
			  &py_function ) ) {
    return NULL;
  }

  const int npar = par.get_size( );

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*) "len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  try {

    sherpa::DifEvo< Func, PyObject*, sherpa::OptFunc< Func, PyObject*, double >,double >
      difevo( func, py_function );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    ierr = difevo( verbose, maxnfev, tol, population_size, seed, xprob,
		   weighting_factor, npar, mylb, myub, mypar, nfev, fval );
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"Unknown exception caught" );
    return NULL;
  }

  if ( ierr < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"function call failed" );
    return NULL;
  }

  return Py_BuildValue( (char*)"(Ndii)", par.return_new_ref(), fval, nfev,
			ierr );
}
static PyObject* py_difevo( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_difevo( self, args, sherpa::fct_ptr( sao_callback_func ) );

}
//*****************************************************************************
//
// py_difevo:  Python wrapper function for C++ function difevo
//
//*****************************************************************************


//*****************************************************************************
//
// py_nm: Python wrapper function for C++ function neldermead
//
//*****************************************************************************
template< typename Func >
static PyObject* py_neldermead( PyObject* self, PyObject* args,
				Func callback_func ) {

  PyObject* py_function=NULL;
  DoubleArray par, step, lb, ub;
  IntArray finalsimplex;
  int verbose, maxnfev, nfev, initsimplex, ierr;
  double fval, tol;

  if ( !PyArg_ParseTuple( args, (char*) "iiiO&dO&O&O&O&O",
			  &verbose,
			  &maxnfev,
			  &initsimplex,
			  CONVERTME(IntArray), &finalsimplex,
			  &tol,
			  CONVERTME(DoubleArray), &step,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &par,
			  &py_function ) ) {
    return NULL;
  }

  const int npar = par.get_size( );

  if ( npar != step.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(step)=%d != len(par)=%d",
		  static_cast<int>( step.get_size( ) ), npar );
    return NULL;
  }

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  try {

    sherpa::NelderMead< Func, PyObject*, double > nm( callback_func,
                                                      py_function );
    std::vector<int> myfinalsimplex( &finalsimplex[0], &finalsimplex[0] +
				     finalsimplex.get_size( ) );
    std::vector<double> mystep( &step[0], &step[0] + step.get_size( ) );
    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    ierr = nm( verbose, maxnfev, tol, npar, initsimplex, myfinalsimplex, mylb,
	       myub, mystep, mypar, nfev, fval );
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;

  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"Unknown exception caught" );
    return NULL;
  }

  if ( ierr < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"function call failed" );
    return NULL;
  }

  return Py_BuildValue( (char*)"(Ndii)", par.return_new_ref(), fval, nfev,
			ierr );

}
static PyObject* py_nm( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_neldermead( self, args, sherpa::fct_ptr( sao_callback_func ) );

}
//*****************************************************************************
//
// py_nm: Python wrapper function for C++ function neldermead
//
//*****************************************************************************

//??
//*****************************************************************************
//
// py_minim: Python wrapper function for C++ function minim
//
//*****************************************************************************
template< typename Func >
static PyObject* py_minim( PyObject* self, PyObject* args,
				Func callback_func ) {

  PyObject* py_function=NULL;
  DoubleArray par, step, lb, ub;
  IntArray finalsimplex;
  int verbose, maxnfev, nfev, iquad, ierr, initsimplex;
  double fval, ftol, simp;

  if ( !PyArg_ParseTuple( args, (char*) "iiiiddO&O&O&O&O",
			  &verbose,
			  &maxnfev,
			  &initsimplex,
                          &iquad,
                          &simp,
			  &ftol,
			  CONVERTME(DoubleArray), &step,
			  CONVERTME(DoubleArray), &lb,
			  CONVERTME(DoubleArray), &ub,
			  CONVERTME(DoubleArray), &par,
			  &py_function ) ) {
    return NULL;
  }

  const int npar = par.get_size( );

  if ( npar != step.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(step)=%d != len(par)=%d",
		  static_cast<int>( step.get_size( ) ), npar );
    return NULL;
  }

  if ( npar != lb.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(lb)=%d != len(par)=%d",
		  static_cast<int>( lb.get_size( ) ), npar);
    return NULL;
  }

  if ( npar != ub.get_size( ) ) {
    PyErr_Format( PyExc_ValueError, (char*)"len(ub)=%d != len(par)=%d",
		  static_cast<int>( ub.get_size( ) ), npar );
    return NULL;
  }

  try {

    std::vector<double> mylb( &lb[0], &lb[0] + npar );
    std::vector<double> myub( &ub[0], &ub[0] + npar );
    std::vector<double> mypar( &par[0], &par[0] + npar );
    std::vector<double> mystep( &step[0], &step[0] + step.get_size( ) );

    const sherpa::Bounds<double> bounds(mylb, myub);
    sherpa::Minim<Func, PyObject*, double> minim(callback_func, py_function );

    std::vector<double> vc(npar*(npar+1)/2);

    minim.minim( mypar, mystep, npar, fval, maxnfev, verbose, ftol, iquad,
                 simp, vc, ierr, nfev, bounds);
    for ( int ii = 0; ii < npar; ++ii )
      par[ ii ] = mypar[ ii ];

    // {
    //   int counter = 0;
    //   for ( int ii = 0; ii < npar; ++ii )
    //     for ( int jj = 0; jj < ii + 1; +jj ) {
    //       covar[ii][jj] = vc[ counter ];
    //       covar[jj][ii] = vc[ counter ];
    //       counter += 1;
    //     }
    // }

  } catch( sherpa::OptErr& oe ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError,
		       (char*) "The parameters are out of bounds\n" );
    return NULL;
  } catch( std::runtime_error& re ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*) re.what() );
    return NULL;
  } catch ( ... ) {
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"Unknown exception caught" );
    return NULL;
  }

  if ( ierr < 0 ) {
    // Make sure an exception is set
    if ( NULL == PyErr_Occurred() )
      PyErr_SetString( PyExc_RuntimeError, (char*)"function call failed" );
    return NULL;
  }

  return Py_BuildValue( (char*)"(Ndii)", par.return_new_ref(), fval, nfev,
			ierr );

}
static PyObject* py_nm_minim( PyObject* self, PyObject* args ) {

  //
  // it looks like an extra indirection but fct_ptr does the nasty
  // work so I do not have to worry about the function prototype.
  //
  return py_minim( self, args, sherpa::fct_ptr( sao_callback_func ) );

}
//*****************************************************************************
//
// py_minim: Python wrapper function for C++ function neldermead
//
//*****************************************************************************
//??

//*****************************************************************************
//
// Module initialization
//
//*****************************************************************************
static PyMethodDef WrapperFcts[] = {


  FCTSPEC(difevo, py_difevo),
  FCTSPEC(nm_difevo, py_difevo_nm),
  FCTSPEC(lm_difevo, py_difevo_lm),
  FCTSPEC(cpp_lmder, py_lmder),
  FCTSPEC(cpp_lmdif, py_lmdif),
  FCTSPEC(neldermead, py_nm),
  FCTSPEC(minim, py_nm_minim),
  { NULL, NULL, 0, NULL }

};

SHERPAMOD(_saoopt, WrapperFcts)
//*****************************************************************************
//
// Module initialization
//
//*****************************************************************************
