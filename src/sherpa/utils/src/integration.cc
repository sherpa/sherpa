// 
//  Copyright (C) 2007, 2016, 2020  Smithsonian Astrophysical Observatory
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

#define _INTEGRATIONMODULE
#include <Python.h>
#include "sherpa/integration.hh"
#include <iostream>
#include <limits>
#include "gsl_errno.h"
#include "gsl_integration.h"
#include "adapt_integrate.h"
#include "qng.h"

extern "C" {
  void initintegration();
}

namespace sherpa { namespace integration {

  static int gsl_int_flag = 1;
  static int sao_int_flag = 1;
  
  int integrate_1d( integrand_1d fct, void* params,
		    double xlo, double xhi,
		    unsigned int maxeval, double epsabs, double epsrel,
		    double& result, double& abserr )
  {

    if ( NULL == fct )
      return EXIT_FAILURE;

    int retval;
    gsl_function func;
    func.function = fct;
    func.params = params;

    size_t neval = size_t( maxeval );

    gsl_set_error_handler_off ();
 
    retval = gsl_integration_qng( &func, xlo, xhi, epsabs, epsrel, &result,
				  &abserr, &neval );
    
    if (retval != EXIT_SUCCESS) {
      if( gsl_int_flag ) {
	std::cerr << "WARNING: Gauss-Kronrod integration failed "
		  << "with tolerance " << epsabs << ", trying lower tolerance..."
		  << std::endl;
	
	double tol = std::numeric_limits< float >::epsilon();
	retval = gsl_integration_qng( &func, xlo, xhi, tol, epsrel, &result,
				      &abserr, &neval );
	
	if (retval != EXIT_SUCCESS) {
	  std::cerr << "integration failed with tolerance " << tol
		    << ", resorting to trapezoid method"
		    << std::endl;

	  result = 0.5 * ( xhi - xlo ) * ( fct( xlo, params) + fct( xhi, params) );
	}
	else {
	  std::cerr << "integration succeeded with tolerance " << tol
		    << std::endl;
	}
	gsl_int_flag = 0;
      }      
    }
    
    return EXIT_SUCCESS;

  }


  int integrate_Nd( integrand_Nd fct, void* params,
                    unsigned int ndim, const double* xlo, const double* xhi,
                    unsigned int maxeval, double epsabs, double epsrel,
                    double& result, double& abserr )
  {

    if ( NULL == fct || NULL == xlo || NULL == xhi )
      return EXIT_FAILURE;

    if ( 0 != adapt_integrate( fct, params, ndim, xlo, xhi, maxeval,
			       epsabs, epsrel, &result, &abserr ) )
      return EXIT_FAILURE;

    return EXIT_SUCCESS;

  }

  int py_integrate_1d( integrand_1d_vec fct, void* params,
		       const double xlo, const double xhi,
		       unsigned int maxeval, double epsabs, double epsrel,
		       double &result, double& abserr, int errflag,
		       std::ostringstream& err)
  {

    if ( NULL == fct )
      return EXIT_FAILURE;

    int retval;

    size_t neval = size_t( maxeval );

    gsl_set_error_handler_off ();
 
    retval = sao_integration_qng( fct, xlo, xhi, params, epsabs, epsrel,
				  &result, &abserr, &neval );
    
    if (retval == -1) {
      return EXIT_FAILURE;
    }
    
    if (retval != EXIT_SUCCESS) {
      if (sao_int_flag) {
	err << "Gauss-Kronrod integration failed "
	    << "with tolerance " << epsabs << ", trying lower tolerance...";
	
	double tol = std::numeric_limits< float >::epsilon();
	retval = sao_integration_qng( fct, xlo, xhi, params, tol, epsrel,
				      &result, &abserr, &neval );
	
	if (retval != EXIT_SUCCESS) {
	  err << std::endl << "integration failed with tolerance " << tol
	      << ", resorting to trapezoid method";

	  
	  double loval[1], hival[1];
	  loval[0] = xlo;
	  hival[0] = xhi;
	  if (-1 == fct(loval, 1, params) )
	    return EXIT_FAILURE;
	  
	  if (-1 == fct(hival, 1, params) )
	    return EXIT_FAILURE;
	  
	  result = 0.5 * ( xhi - xlo ) * ( loval[0] + hival[0] );
	}
	else {
	  err << std::endl << "integration succeeded with tolerance " << tol;
	}
      }
      sao_int_flag=0;
    }
    return EXIT_SUCCESS;
    
  }
    
}  }  /* namespace integration, namespace sherpa */

static struct PyModuleDef integration = {
    PyModuleDef_HEAD_INIT,
    "integration",
    NULL,
    -1,
    NULL
};

PyMODINIT_FUNC PyInit_integration(void) {

  static void *Integration_API[3];
  Integration_API[0] = (void*)sherpa::integration::integrate_1d;
  Integration_API[1] = (void*)sherpa::integration::integrate_Nd;
  Integration_API[2] = (void*)sherpa::integration::py_integrate_1d;

  PyObject *m;
  PyObject *api_cobject;

  if ( NULL == ( m = PyModule_Create( &integration ) ) )
    return NULL;

  if ( NULL == ( api_cobject = PyCapsule_New( (void*)Integration_API,
                              NULL,
						      NULL) ) )
    return NULL;

  // Since the actual data is static, we can let PyModule_AddObject()
  // steal the reference
  PyModule_AddObject( m, (char*)"_C_API", api_cobject );

  return m;
}
