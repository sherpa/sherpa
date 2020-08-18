// 
//  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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


#ifndef __sherpa_constants_hh__
#define __sherpa_constants_hh__

#include <cmath>
#include <limits>

//
// The following constants come from GNU's math.h.  They provide
// enough digits for the 128-bit IEEE quad.
//
#ifndef M_LN2l
# define M_LN2l		0.6931471805599453094172321214581766L  /* log_e 2 */
#endif
#ifndef M_LN10l
# define M_LN10l	2.3025850929940456840179914546843642L  /* log_e 10 */
#endif
#ifndef M_PIl
# define M_PIl		3.1415926535897932384626433832795029L  /* pi */
#endif


namespace sherpa { namespace constants {


  //
  // Min/max values to be used in place of DBL_MIN/DBL_MAX.  ("smp" stands
  // for "Sherpa model parameter", but obviously they can be used other
  // places.)
  //

  template <typename Type>
  inline Type smp_min() { return 1.0E-120; }

  template <typename Type>
  inline Type smp_max() { return 1.0E+120; }

  //
  // Mathematical constants
  //

  template <typename Type>
  inline Type logten() { return M_LN10l; }

  template <typename Type>
  inline Type pi() { return M_PIl; }

  template <typename Type>
  inline Type twopi() { return 2.0L * M_PIl; }

  template <typename Type>
  inline Type sqrt_pi() { return std::sqrt( M_PIl ); }

  // The formula for gfactor is derived towards the end of the file
  template <typename Type>
  inline Type gfactor() { return 4.0L * M_LN2l; }

  template <typename Type>
  inline Type sqrt_gfactor()
  {
    return std::sqrt( gfactor< long double >() );
  }

  //
  // Physical constants
  //

  // nist Electron volt [J]
  template <typename Type>
  inline Type e_ev() { return 1.60217653e-9L; }

  // nist Planck constant [erg * s]
  template <typename Type>
  inline Type h_erg() { return 6.6260693e-27L; }

  // nist Planck constant [KeV * s]
  template <typename Type>
  inline Type h_kev() { return h_erg< Type >() / e_ev< Type >(); }

  // gnu  Speed of light vacuum [km/s]
  template <typename Type>
  inline Type c_km() { return 2.99792458e+5L; }

  // gnu  Speed of light vacuum [cm/s]
  template <typename Type>
  inline Type c_cm() { return 2.99792458e+10L; }

  // gnu  Speed of light vacuum [ang/s]
  template <typename Type>
  inline Type c_ang() { return 2.99792458e+18L; }

  // nist Boltzmann constant [erg / K]
  template <typename Type>
  inline Type k() { return 1.3806505e-16L; }

  template <typename Type>
  inline Type two_h_over_c_squared()
  {
    return 2.0L * h_erg< Type >() / c_cm< Type >() / c_cm< Type >();
  }

  template <typename Type>
  inline Type h_over_k() { return h_erg< Type >() / k< Type >(); }


}  }  /* namespace constants, namespace sherpa */


#ifndef _SHERPA_CONST
#define _SHERPA_CONST(name)	sherpa::constants::name< DataType >()
#endif

#define SMP_MIN			_SHERPA_CONST(smp_min)
#define SMP_MAX			_SHERPA_CONST(smp_max)
#define LOGTEN			_SHERPA_CONST(logten)
#define PI			_SHERPA_CONST(pi)
#define TWOPI			_SHERPA_CONST(twopi)
#define H_ERG			_SHERPA_CONST(h_erg)
#define E_EV			_SHERPA_CONST(e_ev)
#define H_KEV			_SHERPA_CONST(h_kev)
#define C_KM			_SHERPA_CONST(c_km)
#define C_CM			_SHERPA_CONST(c_cm)
#define C_ANG			_SHERPA_CONST(c_ang)
#define K			_SHERPA_CONST(k)
#define TWO_H_OVER_C_SQUARED	_SHERPA_CONST(two_h_over_c_squared)
#define H_OVER_K		_SHERPA_CONST(h_over_k)
#define GFACTOR			_SHERPA_CONST(gfactor)
#define SQRT_GFACTOR		_SHERPA_CONST(sqrt_gfactor)
#define SQRT_PI			_SHERPA_CONST(sqrt_pi)


/*
  dumbo-357: maple
  |\^/|     Maple V Release 2 (Harvard University)
  ._|\|   |/|_. Copyright (c) 1981-1992 by the University of Waterloo.
  \  MAPLE  /  All rights reserved. MAPLE is a registered trademark of
  <____ ____>  Waterloo Maple Software.
  |       Type ? for help.
  > f:= x-> exp( - 1/2 * ( (x-mu)/sigma)^2 );

                              2
                      (x - mu)
  f := x -> exp(- 1/2 ---------)
                             2
                       sigma

  > f( mu - fwhm/2) = f(mu)/2;

                 2
             fwhm
  exp(- 1/8 ------) = 1/2
                 2
            sigma

  > solve(",fwhm);

  1/2      1/2             1/2      1/2
  2 sigma 2    ln(2)   , - 2 sigma 2    ln(2)


  gfactor is defined to be :
  fwhm = 2 * sqrt( 2 * ln( 2 ) )
  gfactor = fwhm^2 / 2
  gfactor = ( 4 * 2 * ln( 2 ) ) / 2 
  gfactor =  4 * ln( 2 )
*/


#endif /* __sherpa_constants_hh__ */
