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

#ifndef __sherpa_models_hh__
#define __sherpa_models_hh__

#include <algorithm>
#include <sherpa/utils.hh>
#include <sherpa/constants.hh>


namespace sherpa { namespace models {


  template <typename DataType, typename ConstDataType>
  inline int lfactorial( ConstDataType arg, DataType& answer )
  {

    if ( arg < 0.0 ) {
      //*answer = NAN;
      return EXIT_FAILURE;
    }

    answer = exp ( LGAMMA( arg + 1 ) );
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int box1d_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( x < p[0] || x > p[1] ) 
      val = 0.0;
    else
      val = p[2];

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int box1d_integrated( const ConstArrayType& p, 
			       DataType xlo, DataType xhi, DataType& val )
  {

    if ( p[1] <= xlo || p[0] >= xhi ) 
      val = 0.0;
    else {
      val = p[2] * ( std::min(xhi,p[1]) - std::max(xlo,p[0]) );
    }

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int const1d_point( const ConstArrayType& p, DataType x,
			    DataType& val )
  {

    val = p[0];
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int const1d_integrated( const ConstArrayType& p,
				 DataType xlo, DataType xhi, DataType& val )
  {

    val = p[0]*(xhi-xlo);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int cos_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp = TWOPI*(x-p[1])/p[0];
    val = p[2]*COS(tmp);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int cos_integrated( const ConstArrayType& p,
			     DataType xlo, DataType xhi, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp1 = TWOPI*(xlo-p[1])/p[0];
    register DataType tmp2 = TWOPI*(xhi-p[1])/p[0];
    val = p[2]*p[0]*(SIN(tmp2)-SIN(tmp1))/TWOPI;
    return EXIT_SUCCESS;
  
  }


  template <typename DataType, typename ConstArrayType>
  inline int delta1d_point( const ConstArrayType& p, DataType x,
			    DataType& val )
  {

    if ( x == p[0] )
      val = p[1];
    else
      val = 0.0;

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int delta1d_integrated( const ConstArrayType& p,
				 DataType xlo, DataType xhi,  DataType& val )
  {

    if ( p[0] >= xlo && p[0] < xhi )
      val = p[1];
    else
      val = 0.0;

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int erf_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( 0.0 != p[2] ) {

      val = ERF( ( x - p[1] ) / p[2] );

    } else {

      if ( x == p[1] ) {  // 0/0 is undefined
	//val = std::numeric_limits< DataType >::quiet_NaN();
	return EXIT_FAILURE;
      }

      if ( x > p[1] )
	val = 1.0;
      else
	val = -1.0;

    }

    val *= p[0];
  
    return EXIT_SUCCESS;

  }


  //
  // The formula for the integral of erf() comes from
  // http://mathworld.wolfram.com/Erf.html
  //

  template <typename DataType, typename ConstArrayType>
  inline DataType _erf_sub1( const ConstArrayType& p, DataType x )
  {
    register DataType arg = ( x - p[1] ) / p[2];
    return arg * ERF( arg ) + EXP( -arg*arg )
      / constants::sqrt_pi< DataType >();
  }


  template <typename DataType, typename ConstArrayType>
  inline DataType _erf_sub2( const ConstArrayType& p, DataType x )
  {
    register DataType arg = ( x - p[1] ) / p[2];
    if ( x < p[1] )  arg *= -1.0;
    return arg;
  }


  template <typename DataType, typename ConstArrayType>
  inline int erf_integrated( const ConstArrayType& p,
			     DataType xlo, DataType xhi, DataType& val )
  {

    if ( 0.0 != p[2] ) {

      val = _erf_sub1( p, xhi ) - _erf_sub1( p, xlo );

    } else {

      if ( ( xlo == p[1] ) || ( xhi == p[1] ) ) {  // 0/0 is undefined
	//val = std::numeric_limits< DataType >::quiet_NaN();
	return EXIT_FAILURE;
      }

      val = _erf_sub2( p, xhi ) - _erf_sub2( p, xlo );

    }

    val *= p[0] * p[2];

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int erfc_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( 0.0 != p[2] ) {

      val = ERFC( ( x - p[1] ) / p[2] );

    } else {

      if ( x == p[1] ) {  // 0/0 is undefined
	//val = std::numeric_limits< DataType >::quiet_NaN();
	return EXIT_FAILURE;
      }

      if ( x > p[1] )
	val = 0.0;
      else
	val = 2.0;

    }

    val *= p[0];

    return EXIT_SUCCESS;

  }


  //
  // The formula for the integral of erfc() comes from
  // http://mathworld.wolfram.com/Erfc.html
  //

  template <typename DataType, typename ConstArrayType>
  inline DataType _erfc_sub1( const ConstArrayType& p, DataType x )
  {
    register DataType arg = ( x - p[1] ) / p[2];
    return arg * ERFC( arg ) - EXP( -arg*arg )
      / constants::sqrt_pi< DataType >();
  }


  template <typename DataType, typename ConstArrayType>
  inline DataType _erfc_sub2( const ConstArrayType& p, DataType x )
  {
    if ( x > p[1] )  return 0.0;
    return 2.0 * ( x - p[1] ) / p[2];
  }


  template <typename DataType, typename ConstArrayType>
  inline int erfc_integrated( const ConstArrayType& p,
			      DataType xlo, DataType xhi, DataType& val )
  {

    if ( 0.0 != p[2] ) {

      val = _erfc_sub1( p, xhi ) - _erfc_sub1( p, xlo );

    } else {

      if ( ( xlo == p[1] ) || ( xhi == p[1] ) ) {  // 0/0 is undefined
	//val = std::numeric_limits< DataType >::quiet_NaN();
	return EXIT_FAILURE;
      }

      val = _erfc_sub2( p, xhi ) - _erfc_sub2( p, xlo );

    }

    val *= p[0] * p[2];

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int exp_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    val = p[2] * EXP( p[1] * (x - p[0]));
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int exp_integrated( const ConstArrayType& p,
			     DataType xlo, DataType xhi, DataType& val )
  {
  
    if ( p[1] != 0.0 ) {
      register DataType y2 = p[1]*(xhi-p[0]);
      register DataType y1 = p[1]*(xlo-p[0]);
      val = (p[2]/p[1])*(EXP(y2)-EXP(y1));
      return EXIT_SUCCESS;
    } else {
      val = p[2]*(xhi-xlo);
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int exp10_point( const ConstArrayType& p, DataType x,
			  DataType& val )
  {

    val = p[2] * POW(10.0, p[1] * (x - p[0]));
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int exp10_integrated( const ConstArrayType& p,
			       DataType xlo, DataType xhi, DataType& val )
  {
  
    if ( p[1] != 0.0 ) {
      register DataType y2 = LOGTEN*p[1]*(xhi-p[0]);
      register DataType y1 = LOGTEN*p[1]*(xlo-p[0]);
      val = (p[2]/p[1]/LOGTEN)*(EXP(y2)-EXP(y1));
      return EXIT_SUCCESS;
    } else {
      val = p[2]*(xhi-xlo);
      return EXIT_SUCCESS;
    }

  }


  //                                                      2
  //                                    GFACTOR (x - p[1])
  //                         p[2] exp(- -------------------)
  //                                               2
  //                                           p[0]
  // 
  template <typename DataType, typename ConstArrayType>
  inline int gauss1d_point( const ConstArrayType& p, DataType x,
			    DataType& val )
  {

    if ( p[0] == 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
  
    // the compiler will take care of optimization.
    val = p[2] * EXP( - GFACTOR * ( x - p[1] ) * ( x - p[1] ) / p[0] / p[0] );

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int gauss1d_integrated( const ConstArrayType& p,
				 DataType xlo, DataType xhi, DataType& val )
  {

    if ( p[0] == 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }

    register DataType z2 = SQRT_GFACTOR * ( xhi - p[1] ) / p[0];
    register DataType z1 = SQRT_GFACTOR * ( xlo - p[1] ) / p[0];

    val =
      p[2] * p[0] * SQRT_PI * ( ERF(z2) - ERF(z1) ) / ( 2. * SQRT_GFACTOR );

    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int log_point( const ConstArrayType& p, DataType x, DataType& val )
  {
  
    if ( (p[1] * (x - p[0])) > 0.0  ) {
      val = p[2] * LOG( p[1] * (x - p[0]));
      return EXIT_SUCCESS;
    } else {
      // val = NAN;
      return EXIT_FAILURE;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int log_integrated( const ConstArrayType& p,
			      DataType xlo, DataType xhi, DataType& val )
  {

    if ( p[1] == 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType y1 = p[1]*(xlo-p[0]);
    register DataType y2 = p[1]*(xhi-p[0]);

    if ( y1 > 0.0 && y2 > 0.0 ) {
      val = p[2]*(y2*LOG(y2)-y1*LOG(y1)-y2+y1)/p[1];
      return EXIT_SUCCESS;
    } else {
      // val = NAN;
      return EXIT_FAILURE;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int log10_point( const ConstArrayType& p, DataType x,
			  DataType& val )
  {

    if ( (p[1] * (x - p[0])) > 0.0  ) {
      val = p[2] * LOG10( p[1] * (x - p[0]));
      return EXIT_SUCCESS;
    } else {
      // val = NAN;
      return EXIT_FAILURE;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int log10_integrated( const ConstArrayType& p,
			       DataType xlo, DataType xhi, DataType& val )
  {

    if ( p[1] == 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType y1 = p[1]*(xlo-p[0]);
    register DataType y2 = p[1]*(xhi-p[0]);

    if ( y1 > 0.0 && y2 > 0.0 ) {
      val = p[2]*(y2*LOG(y2)-y1*LOG(y1)-y2+y1)/p[1]/LOG(10.0);
      return EXIT_SUCCESS;
    }

    else {
      // val = NAN;
      return EXIT_FAILURE;
    }  

  }


  template <typename DataType, typename ConstArrayType>
  inline int ngauss1d_point( const ConstArrayType& p, DataType x,
			     DataType& val )
  {
  
    if( p[0] == 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }

    register DataType norm = SQRT(PI/GFACTOR)*p[0];
    val = (p[2]/norm)*EXP(-GFACTOR*(x-p[1])*(x-p[1])/p[0]/p[0]);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int ngauss1d_integrated( const ConstArrayType& p,
				  DataType xlo, DataType xhi, DataType& val )
  {

    if( p[0] == 0.0 )
      // val = NAN
      return EXIT_FAILURE;
    else {
      register DataType z2 = SQRT_GFACTOR*((xhi-p[1])/p[0]);
      register DataType z1 = SQRT_GFACTOR*((xlo-p[1])/p[0]);
      val = (p[2]*(ERF(z2)-ERF(z1))/2.0);
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int poisson_point( const ConstArrayType& p, DataType x,
			    DataType& val )
  {

    register DataType p_zero_fact;
    register DataType x_fact;

    if( EXIT_SUCCESS != lfactorial(p[0], p_zero_fact)) {
      return EXIT_FAILURE;
    }
  
    if( EXIT_SUCCESS != lfactorial(x, x_fact)) {
      return EXIT_FAILURE;
    }

    if( p[0] > 0.0 ) {
      val = p[1]*EXP((x-p[0])*LOG(p[0])+ p_zero_fact - x_fact);
      return EXIT_SUCCESS;
    }
    else {
      // val = NAN;
      return EXIT_FAILURE;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int poly1d_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    register DataType xtemp = x - p[9];
    register DataType retval = p[8];
    int ii;
    for ( ii = 7; ii >= 0; ii--) {
      retval = retval*xtemp + p[ii];
    }
    val = retval;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int poly1d_integrated( const ConstArrayType& p,
				DataType xlo, DataType xhi, DataType& val )
  {

    register DataType xtemp1 = xlo - p[9];
    register DataType xtemp2 = xhi - p[9];
    register DataType retval = 0.0;
    int ii;
    for( ii = 0; ii <= 8; ii++) {
      register DataType pexp = (DataType)(ii+1);
      retval += p[ii]*(POW(xtemp2,pexp)-POW(xtemp1,pexp))/pexp;
    }
    val = retval;
    return EXIT_SUCCESS;


  }


  //
  //                                    /  x \(- p[0])
  //                               p[2] |----|
  //                                    \p[1]/
  // 
  template <typename DataType, typename ConstArrayType>
  inline int powlaw_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    // p[1] is an always frozen parameter, currently initialized to 1
    // x must be zero or greater
    if ( !(x < 0.0) ) {
      val = p[2] * POW( x / p[ 1 ], - p[ 0 ] );
      return EXIT_SUCCESS;
    }
    
    val = 0.0;
    return EXIT_FAILURE;

  }


  template <typename DataType, typename ConstArrayType>
  inline int powlaw_integrated( const ConstArrayType& p,
				DataType xlo, DataType xhi, DataType& val )
  {

    // xlo must be 0 or greater
    if ( !(xlo < 0.0) ) {
      if ( p[0] == 1.0 ) {
	if ( xlo > 0.0 ) {
	  ;
	}
	else {
	  // If we get here, xlo == 0.0
	  // Stub in Sherpa minimum value for xlo,
	  // so we can take its log
	  xlo = SMP_MIN;
	}	
	val = p[2] * p[1] * ( LOG(xhi) - LOG(xlo) );
	return EXIT_SUCCESS;
      } else {
	register DataType p1 = POW( xlo, 1.0-p[0] ) / (1.0-p[0]);
	register DataType p2 = POW( xhi, 1.0-p[0] ) / (1.0-p[0]);
	val = p[2] / POW( p[1], -p[0] ) * ( p2 - p1 );
	return EXIT_SUCCESS;
      }
    }

    val = 0.0;
    return EXIT_FAILURE;

  }


  //
  //                                    /  x \ - p[1] - p[2] * log10(x/p[0])
  //                               p[3] |----|
  //                                    \p[0]/
  // 
  template <typename DataType, typename ConstArrayType>
  inline int logparabola_point( const ConstArrayType& p, DataType x,
			      DataType& val )
  {

    if ( p[0] != 0 ) {
      
      register DataType frac = (x / p[0]);
      
      if ( frac > 0.0 ) {
	val = p[3] * POW( frac, - p[1] - p[2] * LOG10( frac ) );
	return EXIT_SUCCESS;
      }

    }
    
    val = 0.0;
    return EXIT_FAILURE;
  }


  template <typename DataType, typename ConstArrayType>
  inline int sin_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp = TWOPI*(x-p[1])/p[0];
    val = p[2]*SIN(tmp);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int sin_integrated( const ConstArrayType& p,
			     DataType xlo, DataType xhi, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp1 = TWOPI*(xlo-p[1])/p[0];
    register DataType tmp2 = TWOPI*(xhi-p[1])/p[0];
    val = -1.0*p[2]*p[0]*(COS(tmp2)-COS(tmp1))/TWOPI;
    return EXIT_SUCCESS;
  
  }


  template <typename DataType, typename ConstArrayType>
  inline int sqrt_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( (x-p[0]) < 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
  
    else {
      val = p[1]*SQRT(x-p[0]);
      return EXIT_SUCCESS;
    }
  
  }


  template <typename DataType, typename ConstArrayType>
  inline int sqrt_integrated( const ConstArrayType& p,
			      DataType xlo, DataType xhi, DataType& val )
  {

    if( (xlo-p[0]) < 0.0 || (xhi-p[0]) < 0.0 ) {
      // val = NAN;
      return EXIT_FAILURE;
    }

    register DataType tmp1 = POW(xlo-p[0], 1.50);
    register DataType tmp2 = POW(xhi-p[0], 1.50);
    val = 2.0*p[1]*(tmp2-tmp1)/3.0;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int stephi1d_point( const ConstArrayType& p, DataType x,
			     DataType& val )
  {

    if( x >= p[0] ) {
      val = p[1];
      return EXIT_SUCCESS;
    }
    else {
      val = 0.0;
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int stephi1d_integrated( const ConstArrayType& p,
				  DataType xlo, DataType xhi, DataType& val )
  {

    if( xlo <= p[0] && xhi >= p[0] ) {
      val = (xhi - p[0]) * p[1];
      return EXIT_SUCCESS;
    }
    else if( xlo > p[0] ) {
      val = (xhi - xlo) * p[1];
      return EXIT_SUCCESS;
    }
    else {
      val = 0.0;
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int steplo1d_point( const ConstArrayType& p, DataType x,
			     DataType& val )
  {

    if( x <= p[0] ) {
      val = p[1];
      return EXIT_SUCCESS;
    }
    else {
      val = 0.0;
      return EXIT_SUCCESS;
    }

  }

  
  template <typename DataType, typename ConstArrayType>
  inline int steplo1d_integrated( const ConstArrayType& p,
				  DataType xlo, DataType xhi, DataType& val )
  {

    if( xlo <= p[0] && xhi >= p[0] ) {
      val = (p[0] - xlo) * p[1];
      return EXIT_SUCCESS;
    }
    else if( xhi < p[0] ) {
      val = (xhi - xlo) * p[1];
      return EXIT_SUCCESS;
    }
    else {
      val = 0.0;
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int tan_point( const ConstArrayType& p, DataType x, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp = TWOPI*(x-p[1])/p[0];
    val = p[2]*TAN(tmp);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int tan_integrated( const ConstArrayType& p,
			     DataType xlo, DataType xhi, DataType& val )
  {

    if ( 0.0 == p[0] ) {
      // val = NAN;
      return EXIT_FAILURE;
    }
    register DataType tmp1 = TWOPI*(xlo-p[1])/p[0];
    register DataType tmp2 = TWOPI*(xhi-p[1])/p[0];
    val = -1.0*p[2]*p[0]*(LOG(COS(tmp2))-LOG(COS(tmp1)))/TWOPI;
    return EXIT_SUCCESS;
  
  }


  // =============================================================


  template <typename DataType, typename ConstArrayType>
  inline int box2d_point( const ConstArrayType& p,
			  DataType x0, DataType x1, DataType& val )
  {
  
    if ( p[1] <= x0 || p[0] >= x0 || p[3] <= x1 || p[2] >= x1 )
      val = 0.0;
    else
      val = p[4];
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int box2d_integrated( const ConstArrayType& p,
			       DataType x0lo, DataType x0hi,
			       DataType x1lo, DataType x1hi, DataType& val )
  {
 
    if ( p[1] <= x0lo || p[0] >= x0hi || p[3] <= x1lo || p[2] >= x1hi ) {
      val = 0.0;
      return EXIT_SUCCESS;
    } else {
      register DataType x0frac = (std::min(x0hi,p[1])-std::max(x0lo,p[0]))/(x0hi-x0lo);
      register DataType x1frac = (std::min(x1hi,p[3])-std::max(x1lo,p[2]))/(x1hi-x1lo);
      val = p[4]*x0frac*x1frac;
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int const2d_point( const ConstArrayType& p,
			    DataType x0, DataType x1, DataType& val )
  {

    (void)x0;
    (void)x1; 
    val = p[0];
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int const2d_integrated( const ConstArrayType& p,
				 DataType x0lo, DataType x0hi,
				 DataType x1lo, DataType x1hi, DataType& val )
  {
  
    val = p[0]*(x1hi-x1lo)*(x0hi-x0lo);
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int delta2d_point( const ConstArrayType& p,
			    DataType x0, DataType x1, DataType& val )
  {
  
    if ( p[0] == x0 && p[1] == x1 ) {
      val = p[2];
      return EXIT_SUCCESS;
    }
    val = 0.0;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int delta2d_integrated( const ConstArrayType& p,
				 DataType x0lo, DataType x0hi,
				 DataType x1lo, DataType x1hi, DataType& val )
  {

    if ( p[0] >= x0lo && p[0] < x0hi && p[1] >= x1lo && p[1] < x1hi ) {
      val = p[2];
      return EXIT_SUCCESS;
    }
    val = 0.0;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int gauss2d_point( const ConstArrayType& p,
			    DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r = 0.0;
  
    if( EXIT_SUCCESS != sherpa::utils::radius2(p, x0, x1, r)) {
      return EXIT_FAILURE;
    }
    if( p[0] == 0.0 )
      return EXIT_FAILURE;
    else {
      val = p[5]*EXP(-r/(p[0]*p[0])*GFACTOR);
      return EXIT_SUCCESS;
    }

  }

  template <typename DataType, typename ConstArrayType>
  inline int sigmagauss2d_point( const ConstArrayType& p,
				 DataType x0, DataType x1, DataType& val ) {
  
    register DataType r = 0.0;
    
    if ( 0 == p[0] || 0 == p[1] )
      return EXIT_FAILURE;
  
    if( EXIT_SUCCESS != sherpa::utils::sigmaradius2(p, x0, x1, r)) {
      return EXIT_FAILURE;
    }
    val = p[5]*EXP(-r/2.0);
    return EXIT_SUCCESS;

  }


  // FIXME: work the analytic part of this back in
  /*
  template <typename DataType, typename ConstArrayType>
  inline int gauss2d_integrated( const ConstArrayType& p,
				 DataType x0lo, DataType x0hi,
				 DataType x1lo, DataType x1hi, DataType& val )
  {

    if ( p[3] == 0 ) {
      if ( 0 == p[0] ) {
	return EXIT_FAILURE;
      } else {
	DataType prefix = p[5]*p[0]*p[0]*PI/(4.0*GFACTOR);
	DataType zx2 = SQRT_GFACTOR*((x0hi-p[1])/p[0]);
	DataType zx1 = SQRT_GFACTOR*((x0lo-p[1])/p[0]);
	DataType zy2 = SQRT_GFACTOR*((x1hi-p[2])/p[0]);
	DataType zy1 = SQRT_GFACTOR*((x1lo-p[2])/p[0]);
	val = prefix*(ERF(zx2)-ERF(zx1))*(ERF(zy2)-ERF(zy1));
	return EXIT_SUCCESS;
      }
    } else {
      DataType lolft;
      DataType lorgt;
      DataType hilft;
      DataType hirgt;

      if( EXIT_SUCCESS != gauss2d_point(p,x0lo,x1lo, lolft)) {
	return EXIT_FAILURE;
      }
    
      if( EXIT_SUCCESS != gauss2d_point(p,x0hi,x1lo, lorgt)) {
	return EXIT_FAILURE;
      }
    
      if( EXIT_SUCCESS != gauss2d_point(p,x0lo,x1hi, hilft)) {
	return EXIT_FAILURE;
      }
    
      if( EXIT_SUCCESS != gauss2d_point(p,x0hi,x1hi, hirgt)) {
	return EXIT_FAILURE;
      }
      if( EXIT_SUCCESS !=
	  sherpa::utils::numerical_integration2d
	  ( gauss2d_point< DataType, ConstArrayType >,
	    p, x0lo, x0hi, x1lo, x1hi, lolft, lorgt,
	    hilft, hirgt, val )) {
	return EXIT_FAILURE;
      }
      else
	return EXIT_SUCCESS;
    }

    val = 0.0;
    return EXIT_SUCCESS;

  }
  */


  template <typename DataType, typename ConstArrayType>
  inline int ngauss2d_point( const ConstArrayType& p,
			     DataType x0, DataType x1, DataType& val )
  {
  
    register DataType r = 0.0;
  
    if( EXIT_SUCCESS != sherpa::utils::radius2(p, x0, x1, r)) {
      return EXIT_FAILURE;
    }
    if( p[0] == 0.0 )
      return EXIT_FAILURE;
    else {
      register DataType norm = (PI/GFACTOR)*p[0]*p[0]*SQRT(1.0 - (p[3]*p[3]));
      val = (p[5]/norm)*EXP(-r/(p[0]*p[0])*GFACTOR);
      return EXIT_SUCCESS;
    }

  }


  template <typename DataType, typename ConstArrayType>
  inline int poly2d_point( const ConstArrayType& p,
			   DataType x0, DataType x1, DataType& val )
  {
  
    register int ix;
    register int iy;
    register DataType retval = 0.0;
    for ( ix = 0 ; ix < 3 ; ix++ ) {
      for ( iy = 0 ; iy < 3 ; iy++ ) {
	retval += POW(x0,(DataType)ix)*POW(x1,(DataType)iy)*p[3*ix+iy];
      }
    }
    val = retval;
    return EXIT_SUCCESS;

  }


  template <typename DataType, typename ConstArrayType>
  inline int poly2d_integrated( const ConstArrayType& p,
				DataType x0lo, DataType x0hi,
				DataType x1lo, DataType x1hi, DataType& val )
  {
  
    register int ix;
    register int iy;
    register DataType retval = 0.0;
    register DataType u[3];
    register DataType v[3];
    u[0] = x0hi - x0lo;
    u[1] = (POW(x0hi,2.0)/2.0) - (POW(x0lo,2.0)/2.0);
    u[2] = (POW(x0hi,3.0)/3.0) - (POW(x0lo,3.0)/3.0);
    v[0] = x1hi - x1lo;
    v[1] = (POW(x1hi,2.0)/2.0) - (POW(x1lo,2.0)/2.0);
    v[2] = (POW(x1hi,3.0)/3.0) - (POW(x1lo,3.0)/3.0);
    for ( ix = 0 ; ix < 3 ; ix++ ) {
      for ( iy = 0 ; iy < 3 ; iy++ ) {
	retval += u[ix]*v[iy]*p[3*ix+iy];
      }
    }
    val = retval;
    return EXIT_SUCCESS;

  }


}  }  /* namespace models, namespace sherpa */


#endif /* __sherpa_models_hh__ */
