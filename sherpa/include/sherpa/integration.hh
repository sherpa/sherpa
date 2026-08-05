//
//  Copyright (C) 2007, 2016, 2020, 2026  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_integration_hh__
#define __sherpa_integration_hh__

#include <sstream>

extern "C" {
  typedef double (*integrand_1d)( double x, void* params );
  typedef double (*integrand_Nd)( unsigned int ndim, const double* x,
				  void* params );
  typedef int (*integrand_1d_vec)( double* x, int len, void* params );
}

#define INTEGRATION_CAPSULE_NAME "sherpa.utils.integration._C_API"


#if defined(_INTEGRATIONMODULE) || !defined(Py_PYTHON_H)

namespace sherpa { namespace integration {

int integrate_1d( integrand_1d fct, void* params,
		  double xlo, double xhi,
		  unsigned int maxeval, double epsabs, double epsrel,
		  double& result, double& abserr );


int integrate_Nd( integrand_Nd fct, void* params,
		  unsigned int ndim, const double* xlo, const double* xhi,
		  unsigned int maxeval, double epsabs, double epsrel,
		  double& result, double& abserr );


int py_integrate_1d( integrand_1d_vec fct, void* params,
		     const double xlo, const double xhi,
		     unsigned int maxeval, double epsabs, double epsrel,
		     double &result, double& abserr, int errflag,
		     std::ostringstream& err);


}  }  /* namespace integration, namespace sherpa */

#endif  /* defined(_INTEGRATIONMODULE) || !defined(Py_PYTHON_H) */


#if !defined(_INTEGRATIONMODULE) && defined(Py_PYTHON_H)


typedef int (*_integrate_1d)( integrand_1d fct,
			      void* params, double xlo, double xhi,
			      unsigned int maxeval,
			      double epsabs, double epsrel,
			      double& result, double& abserr );

typedef int (*_integrate_Nd)( integrand_Nd fct,
			      void* params, unsigned int ndim,
			      const double* xlo, const double* xhi,
			      unsigned int maxeval,
			      double epsabs, double epsrel,
			      double& result, double& abserr );

typedef int (*_py_integrate_1d)( integrand_1d_vec fct, void* params,
				 const double xlo, const double xhi,
				 unsigned int maxeval,
				 double epsabs, double epsrel,
				 double &result, double& abserr, int errflag,
				 std::ostringstream& err );

static void **Integration_API;

#define integrate_1d ((_integrate_1d)Integration_API[0])
#define integrate_Nd ((_integrate_Nd)Integration_API[1])
#define py_integrate_1d ((_py_integrate_1d)Integration_API[2])


static int
import_integration(void)
{

  Integration_API = (void **) PyCapsule_Import(INTEGRATION_CAPSULE_NAME, 0);
  return (Integration_API != NULL) ? 0 : -1;
}


#endif  /* !defined(_INTEGRATIONMODULE) && defined(Py_PYTHON_H) */


#endif  /* __sherpa_integration_hh__ */
