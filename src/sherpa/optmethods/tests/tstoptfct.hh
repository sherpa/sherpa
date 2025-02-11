#ifndef tstoptfct_hh
#define tstoptfct_hh

//
//  Copyright (C) 2007, 2017, 2019, 2020, 2021
//     Smithsonian Astrophysical Observatory
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


//
// J. MORE', B. GARBOW & K. HILLSTROM,
// "Algorithm 566: Fortran Subroutines for Testing Unconstrained
// Optimization Software.", ACM TOMS, VOL. 7, PAGES 14-41 AND 136-140, 1981
//
// See also https://www.osti.gov/biblio/6650344
//
// The comments on the solutions were taken from the site:
//
// https://web.archive.org/web/20060221201940/http://www.uni-graz.at/imawww/kuntsevich/solvopt/results/moreset.html
//
// Jan 2008 D. T. Nguyen
//
//

#include <cstdio>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

/*
 **********Sample of user supplied function****************
 * m = number of functions
 * n = number of variables
 * x = vector of function arguments
 * fvec = vector of function values
 * iflag = error return variable
 */

namespace tstoptfct {

  template<typename Real>
  Real mysqr( Real arg ) {
    return arg * arg;
  }

  template<typename Real, typename Type>
  void Ackley( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    const Real a = 20.0, b = 0.2, c = 2 * pi;
    Real s1 = 0.0, s2 = 0.0;
    for ( int ii = 0; ii < npar; ++ii ) {
      s1 += mysqr( x[ ii ] );
      s2 += cos( c * x[ ii ] );
    }
    fval = -a * exp( -b * sqrt( 1.0 / npar * s1 ) ) -
      exp( 1.0 / npar * s2 ) + a + exp( 1.0 );
  }
  template<typename Real>
  void AckleyInit( int npar, int& mfct, Real& answer, Real* x,
		   Real* lo, Real* hi ) {
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = - 15.0;
      hi[ ii ] =   30.0;
      x[ ii ] = 17.0;
    }
    answer = 0.0;
  }
  //
  // f( 0, 0, ...., 0 ) = 0
  //

  template<typename Real, typename Type>
  void Bard( int mfct, int npar, Real* x, Real* fvec, int& ierr, Type xptr ) {

    Real yi[] = { 0.14, 0.18, 0.22, 0.25, 0.29, 0.32, 0.35, 0.39, 0.37, 0.58,
		  0.73, 0.96, 1.34, 2.10, 4.39 };

    for ( int ii = 0; ii < npar; ii += 3 )
      for ( int jj = 0; jj < 15; jj++ ) {
	Real wi = 0.0;
	Real ui = jj + 1;
	Real vi = 15.0 - jj;
	if ( ui < vi )
	  wi = ui;
	else
	  wi = vi;
	fvec[ 15 * ii/3 + jj ] =
	  yi[ jj ] - ( x[ ii ] + ui / ( vi * x[ ii + 1 ] + wi * x[ ii + 2 ] ) );
      }

  }
  template<typename Real, typename Type>
  void Bard( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 15 * npar / 3;
    std::vector< Real> fvec( mfct );
    Bard( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BardInit( int npar, int& mfct, Real& answer, Real* x,
		 Real* lo, Real* hi ) {

    if ( npar % 3 )
      throw std::runtime_error( "npar for the Bard func must be multiple of 3\n" );

    mfct = 15 * npar / 3;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e12;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = (npar/3*8.21487e-3);

  }
  //
  // This function has the only minimum
  // f( 0.0824106, 1.13304, 2.3437) = 0.00821487
  // However, at x(1)=0.84066666669332, it becomes flatter and flatter
  // as the other two variables (x(2) and x(3)) decrease.  Starting at
  // x(2),x(3)=-1.e-10, the gradient equals zero to the extend of a
  // machine epsilon.
  //

  template<typename Real, typename Type>
  void Beale( int mfct, int npar, Real* x, Real* fvec, int& ierr,
	      Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 2 ) {
      fvec[ 3*ii/2 ] = 1.5 - x[ ii ] * ( 1.0 - x[ ii + 1 ] );
      fvec[ 3*ii/2 + 1 ] = 2.25 - x[ ii ] *
	( 1.0 - x[ ii + 1 ] * x[ ii + 1 ] );
      fvec[ 3*ii/2 + 2 ] = 2.625  - x[ ii ] *
	( 1.0 - x[ ii + 1 ] * x[ ii + 1 ] * x[ ii + 1 ]);
    }

  }
  template<typename Real, typename Type>
  void Beale( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 15 * npar / 3;
    std::vector< Real> fvec( mfct );
    Beale( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BealeInit( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the Beale func must be even\n" );

    mfct = 3 * npar/2;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // f( 3, 0.5 ) = 0.0
  //

  template<typename Real, typename Type>
  void Booth( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 2 != npar ) {
      if ( npar % 2 ) {
        throw std::runtime_error( "npar for the Booth func must be 2\n" );
      }
      return;
    }

    fval = pow( x[ 0 ] + 2 * x[ 1 ] - 7, 2.0 ) +
      pow( 2 * x[ 0 ] + x[ 1 ] - 5, 2.0 );

  }
  template<typename Real>
  void BoothInit( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the Beale func must be even\n" );

    mfct = 0;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = -5.0 + ii;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e1;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e1;

    answer = 0.0;

  }
  //
  // f( 1, 3 ) = 0.0
  //

  template<typename Real, typename Type>
  void Biggs( int mfct, int npar, Real* x, Real* fvec, int& ierr,
	      Type xptr ) {

    for ( int i = 0; i < mfct; i++ ) {
      Real ti = 0.1 * i;
      Real yi = exp( - ti * x[0] ) - 5.0 * exp( -10.0 * ti ) +
	3.0 * exp( -4.0 * ti );
      fvec[ i ] = x[2] * exp( - ti * x[0] ) - x[3] * exp( - ti * x[1] ) +
	x[5] *exp( - ti * x[4] ) - yi;
    }

  }
  template<typename Real, typename Type>
  void Biggs( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 6;
    std::vector< Real> fvec( mfct );
    Biggs( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BiggsInit( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi ) {

    if ( 6 != npar )
      throw std::runtime_error( "npar for the Biggs func must be 6\n" );

    mfct = 6;
    x[ 0 ] = 1.0;
    x[ 1 ] = 10.0;
    x[ 2 ] = 1.0;
    x[ 3 ] = 5.0;
    x[ 4 ] = 4.0;
    x[ 5 ] = 3.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e3;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e3;

    answer = 0.0;

  }
  //
  // This function has the global minimum f( 1, 10, 1, 5, 4, 3 ) = 0.0
  // if m = 13 and local minimum
  // f( 1.711416, 17.6831975, 1.16314365, 5.18656, 1.711416,
  //    1.16314365 ) = 0.00565565
  // at the point , which are known. In addition, the function takes a local
  // minimum 0.30636677262479 on the plane
  // {x1=1.22755594752403; x2>>0; x3=0.83270306333466; x4<<0; x5=x1; x6=x3}
  //

  template<typename Real, typename Type>
  void Bohachevsky1( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = mysqr( x[ 0 ] ) + 2 * mysqr( x[ 1 ] ) -
      0.3 * cos( 3 * pi * x[ 0 ] ) - 0.4 * cos(4 * pi * x[ 1 ] ) + 0.7;
  }
  template < typename Real>
  void Bohachevsky1Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Bohachevsky func must be 2\n" );

    for ( int ii = 0; ii < npar; ++ii ) {

      lo[ ii ] = - 100.0;
      hi[ ii ] =   100.0;

      x[ ii ] = 35.0;

    }
    answer = 0.0;

  }
  //
  // f( 0, 0 ) = 0
  //

  template<typename Real, typename Type>
  void Bohachevsky2( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = mysqr( x[ 0 ] ) + 2 * mysqr( x[ 1 ] ) -
      0.3*cos( 3*pi*x[ 0 ] ) * cos( 4*pi*x[ 1 ] ) + 0.3;
  }
  template < typename Real>
  void Bohachevsky2Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Bohachevsky func must be 2\n" );

    for ( int ii = 0; ii < npar; ++ii ) {

      lo[ ii ] = - 100.0;
      hi[ ii ] =   100.0;

      x[ ii ] = 35.0;

    }
    answer = 0.0;

  }

  template<typename Real, typename Type>
  void Bohachevsky3( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = mysqr( x[ 0 ] ) + 2 * mysqr( x[ 1 ] ) -
      0.3*cos( 3*pi*x[ 0 ] + 4*pi*x[ 1 ] ) + 0.3;
  }
  template < typename Real>
  void Bohachevsky3Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Bohachevsky func must be 2\n" );

    for ( int ii = 0; ii < npar; ++ii ) {

      lo[ ii ] = - 100.0;
      hi[ ii ] =   100.0;

      x[ ii ] = 35.0;

    }
    answer = 0.0;

  }

  template<typename Real, typename Type>
  void BoxBetts( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    fval = 0.0;
    for ( int ii = 1; ii <= 10; ++ii ) {
      Real e0 = exp( -0.1 * ii * x[ 0 ] );
      Real e1 = exp( -0.1 * ii * x[ 1 ] );
      Real e2 = exp( -0.1 * ii ) - exp( static_cast< Real>( - ii ) );
      Real tmp = e0 - e1 - e2 * x[ 2 ];
      fval += tmp * tmp;
    }

  }
  template < typename Real>
  void BoxBettsInit( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 3 != npar )
      throw std::runtime_error( "npar for the BoxBetts func must be 3\n" );

    x[ 0 ] = 1.1; x[ 1 ] = 9.5; x[ 2 ] = 0.99;
    lo[ 0 ] = 0.9; lo[ 1 ] = 9.0; lo[ 2 ] = 0.9;
    hi[ 0 ] = 1.2; hi[ 1 ] = 11.2; hi[ 2 ] = 1.2;
    answer = 0.0;

  }
  //
  // f( 1, 10, 1 ) = 0.0
  //

  template<typename Real, typename Type>
  void Box3d( int mfct, int npar, Real* x, Real* fvec, int& ierr,
	      Type xptr ) {

    for ( int ii = 0; ii < mfct; ++ii ) {
      Real ti = ( ii + 1 ) * 0.1;
      fvec[ ii ] =  exp( -ti * x[ 0 ] ) - exp( -ti * x[ 1 ] ) -
	x[ 2 ] * ( exp( - ti ) - exp( - 10.0 * ti ) );
    }

  }
  template<typename Real, typename Type>
  void Box3d( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 16;
    std::vector< Real> fvec( mfct );
    Box3d( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template < typename Real>
  void Box3dInit( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi ) {

    if ( 3 != npar )
      throw std::runtime_error( "npar for the Box3d func must be 3\n" );

    mfct = 16;
    x[ 0 ] = 0.0;
    x[ 1 ] = 10.0;
    x[ 2 ] = 20.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e2;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e2;

    answer = 0.0;

  }
  //
  // f=0 at (1,10,1), (10,1,-1) and wherever x1=x2 and x3-0
  //

  template<typename Real, typename Type>
  void Branin( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Branin func must be 2\n" );

    double x0 = x[0];
    double s  = x[1] - ( 5.1/(4.0*M_PI*M_PI) * x0 - 5.0/M_PI) * x0 - 6.0;
    fval = s*s + 10*(1.0 - 1.0/(8.0*M_PI)) * cos(x0) + 10.0;

  }
  template<typename Real>
  void BraninInit( int npar, int& mfct, Real& answer, Real* x,
			      Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Branin func must be 2\n" );
    lo[ 0 ] = -5.0; lo[ 1 ] = 0.0;
    hi[ 0 ] = 10.0; hi[ 1 ] = 15.0;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lo[ ii ] + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.397889;
  }
  //
  // f( -3.142, 12.275 ) = f( 3.142, 2.275 ) = f( 9.425, 2.425 ) = 0.397889
  //

  template<typename Real, typename Type>
  void Branin2( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Branin2 func must be 2\n" );

    fval = pow( 1.0 - 2.0 * x[ 1 ] + sin( 4.0 * M_PI * x[ 1 ] ) / 20.0 -
		x[ 0 ], 2.0 ) +
      pow( x[ 1 ] - sin( 2.0 * M_PI * x[ 0 ] ) / 2.0, 2.0 );
  }
  template<typename Real>
  void Branin2Init( int npar, int& mfct, Real& answer, Real* x,
			      Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Branin2 func must be 2\n" );
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -10.0;
      hi[ ii ] = 10.0;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lo[ ii ] + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0;
  }

  template<typename Real, typename Type>
  void BrownAlmostLinear( int mfct, int npar, Real* x, Real* fvec,
			  int& ierr, Type xptr ) {

    Real sum = 0.0;
    Real prd = 1.0;
    for ( int ii = 1; ii <= npar; ++ii ) {
      Real xi  = x[ ii - 1 ];
      sum += xi;
      prd *= xi;
      fvec[ ii - 1 ] = xi - ( npar + 1.0 );
    }

    int nvm1 = npar - 1;
    for ( int ii = 1; ii <= nvm1; ++ii )
      fvec[ ii - 1 ] += sum;

    fvec[ npar - 1 ] = prd - 1.0;

  }
  template<typename Real, typename Type>
  void BrownAlmostLinear( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    BrownAlmostLinear( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BrownAlmostLinearInit( int npar, int& mfct, Real& answer, Real* x,
			      Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 0.5;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 1.0;

  }

  template<typename Real, typename Type>
  void BrownBadlyScaled( int mfct, int npar, Real* x, Real* fvec,
			 int& ierr, Type xptr ) {
    for ( int ii = 0; ii < npar; ii += 2 ) {
      fvec[ ii ] = x[ ii ] - 1.0e6;
      fvec[ ii + 1 ] = x[ ii + 1 ] - 2.0e-6;
      fvec[ ii + 2 ] = x[ ii ] * x[ ii + 1 ] - 2.0;
    }

  }
  template<typename Real, typename Type>
  void BrownBadlyScaled( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar + npar/2;
    std::vector< Real> fvec( mfct );
    BrownBadlyScaled( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BrownBadlyScaledInit( int npar, int& mfct, Real& answer, Real* x,
			     Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the BrownBadlyScaled func must be even\n" );

    mfct = npar + npar/2;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e2;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e9;

    answer = 0.0;

  }
  //
  // f( 1e+6, 2e-6 ) = 0.0
  //

  template<typename Real, typename Type>
  void BrownDennis( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		    Type xptr ) {

    for ( int ii = 0; ii < mfct; ++ii ) {
      Real ti = (ii+1)/5.0;
      fvec[ ii ] = pow( x[0] + ti * x[1] - exp( ti ), 2.0 ) +
	pow( x[2] + x[3]*sin(ti) - cos(ti), 2.0 );
    }

  }
  template<typename Real, typename Type>
  void BrownDennis( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 20;
    std::vector< Real> fvec( mfct );
    BrownDennis( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BrownDennisInit( int npar, int& mfct, Real& answer, Real* x,
			Real* lo, Real* hi ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the BrownDennis func must be 4\n" );

    mfct = 20;
    x[ 0 ] = 25.0;
    x[ 1 ] = 5.0;
    x[ 2 ] = -5.0;
    x[ 3 ] = -1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 85822.2;
  }

  template<typename Real, typename Type>
  void BroydenBanded( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		      Type xptr ) {

    const int ilbw   = 5;
    const int iubw   = 1;
    for ( int i = 1; i <= npar; i++ ) {
      Real xi   = x[ i - 1 ];
      fvec[ i - 1 ] = xi*(2.0 + 5.0*xi*xi) + 1.0;
      int j1 = 1 > i - ilbw ? 1 : i - ilbw;
      int j2 = npar < i + iubw ? npar : i + iubw;
      for ( int j = j1; j <= j2; j++ ) {
	Real xj = x[ j - 1 ];
	if ( j != i )
	  fvec[ i - 1 ] -= xj*(1.0 + xj);
      }
    }

  }
  template<typename Real, typename Type>
  void BroydenBanded( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    BroydenBanded( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BroydenBandedInit( int npar, int& mfct, Real& answer, Real* x,
			  Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = -1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // In addition to the global minimum f0=0, this function has at least
  // two more local minima:
  //
  // f( 0.210713858177301, 0.260976685957722, 0.227947672085427,
  //    0.197025588320616, 0.203229694083309, 0.221915978634404,
  //    0.16281223895881, 0.0483025286174326, -0.0793565823288942,
  //    -0.153125186088696] = 3.05727843
  //
  // f( 0.21717225508353, 0.30646717454530, 0.35678613574215,
  //    0.44839494272506, 0.57763800138617, 0.71310355796924,
  //    0.82395652593316, 0.92190937955671, 0.99699281917885,
  //    0.96346230718557 ) = 2.68021992072616
  //

  template<typename Real, typename Type>
  void BroydenTridiagonal( int mfct, int npar, Real* x, Real* fvec,
			   int& ierr, Type xptr ) {

    for ( int ii = 1; ii <= npar; ++ii ) {
      Real xiiminus1 = ( ii == 1 ) ? 0.0 :  x[ ii - 1 ];
      Real xiiplus1  = ( ii == npar ) ? 0.0 :  x[ ii ];
      fvec[ ii - 1 ] = ( 3.0 - 2 * x[ ii - 1 ] ) * x[ ii - 1 ] - xiiminus1 -
	2 * xiiplus1 + 1.0;
    }

  }
  template<typename Real, typename Type>
  void BroydenTridiagonal( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    BroydenTridiagonal( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void BroydenTridiagonalInit( int npar, int& mfct, Real& answer, Real* x,
			       Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = -1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // In addition to the global minimum f0=0, this function has a number of
  // local minima. They are:
  //
  // f( 0.89571637066795, 1.43694322050304, -0.08787533663029;
  //    -0.51565133590034, -0.56252955081589, -0.42678887777818;
  //    -0.03488113362957, 0.73773668071208, 1.31192604432832;
  //    0.23349654053705 ) = 1.36025590473840
  //
  // f( 1.57293556968858, 0.37536914085555, 0.16829821000091,
  //    0.60512754208275, 1.05196056858588, 0.58619628989910,
  //    0.42021181170450, 0.76713494858584, 1.17153397159970,
  //    0.26636104677652 ) = 1.02865203567795
  //
  // f( 1.76423409318896, 0.03091840088999, -0.32488644889735,
  //    -0.06855209483793, 0.70236806246591, 1.49772106403168,
  //    -0.06585213920739, -0.50614689954560, -0.54944142333092,
  //    -0.40108999846985 ) = 1.05122618838356
  //
  // f( 1.73905838100313, 0.08286257408083, -0.24655906721381,
  //    0.04553098676945, 0.74854483640949, 1.17092063399475,
  //    0.39373320226816, 0.28526714046312, 0.79188653014133,
  //    1.31220034435141 ) = 1.05122618838356
  //
  // f( 1.83141075408277, -0.10695150598288, -0.58756136817458,
  //    -0.67305617728853, -0.66862208919252, -0.61247812756592,
  //    -0.45433311772938, -0.05529093386451, 0.75788912482266,
  //    1.41326822790186 ) = 0.71260601731262
  //
  // f( 1.82910567383287, -0.10197809329322, -0.57788049656156,
  //    -0.64962654169425, -0.60679005996046, -0.44981389867908,
  //    -0.05474342424961, 0.72306839424618, 1.31937148707469,
  //    0.23503639006362 ) = 0.39737346895853.
  //

  template<typename Real, typename Type>
  void Chichinadze( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Chichinadze func must be 2\n" );

    fval = x[ 0 ] * x[ 0 ] - 12 * x[ 0 ] + 11 + 10 * cos( M_PI / 2.0 * x[ 0 ] )
      + 8 * sin( 5 * M_PI * x[ 0 ] ) - exp( -(x[ 1 ] - 0.5) * 0.5) / sqrt(5.0);
  }
  template<typename Real>
  void ChichinadzeInit( int npar, int& mfct, Real& answer, Real* x,
			Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Chichinadze func must be 2\n" );
    lo[ 0 ] = -30;
    lo[ 1 ] = -30;
    hi[ 0 ] = 30.0;
    hi[ 1 ] = 30.0;
    x[ 0 ] = 25.0;
    x[ 1 ] = 9.0;
    answer = - 43.3159;

  }
  //
  // f( 5.90133, 0.5 ) = - 43.3159
  //

  template<typename Real, typename Type>
  void Chebyquad( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		  Type xptr ) {

    for ( int  ii = 0; ii < mfct; ++ii )
      fvec[ ii ] = 0.0;

    for ( int ii = 1; ii <= npar; ++ii ) {
      Real tmp1 = 1.0;
      Real tmp2 = 2.0 * x[ ii - 1 ] - 1.0;
      Real temp  = 2.0 *tmp2;
      for ( int jj = 1; jj <= mfct; ++jj ) {
	fvec[ jj - 1 ] += tmp2;
	Real ti = temp*tmp2 - tmp1;
	tmp1 = tmp2;
	tmp2 = ti;
      }
    }

    Real dx = 1.0 / npar;
    int iev = -1;
    for ( int ii = 1; ii <= mfct; ++ii ) {
      fvec[ ii - 1 ] *= dx;
      if ( iev > 0 )
	fvec[ ii - 1 ] += 1.0 / ( ii * ii - 1.0);
      iev = -iev;
    }

  }
  template<typename Real, typename Type>
  void Chebyquad( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    Chebyquad( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void ChebyquadInit( int npar, int& mfct, Real& answer, Real* x,
		      Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = ( ii + 1.0 ) / ( npar + 1.0 );

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }

  template<typename Real, typename Type>
  void Cola( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    static Real d_mds[] = { 0.0, 1.27, 1.69, 1.43, 2.04, 2.35, 2.43, 3.09,
			    3.18, 3.26, 2.85, 3.20, 3.22, 3.27, 2.88, 1.55,
			    2.86, 2.56, 2.58, 2.59, 3.12, 3.06, 3.17, 3.18,
			    3.18, 3.12, 1.31, 1.64, 3.00, 3.21, 3.18, 3.18,
			    3.17, 1.70, 1.36, 2.95, 1.32, 2.38, 2.31, 2.42,
			    1.94, 2.85, 2.81, 2.56, 2.91, 2.97 };

    fval = 0.0;
    std::vector< Real> mt( 20, 0.0 );
    for( int i = 4; i < 20; i++ )
      mt[ i ] = x[ i-3 ];
    for( int i = 1, k = 1; i < 10; i++ )
      for( int j = 0; j < i; j++ ) {
	Real temp = 0.0;
	for( int t = 0; t < 2; t++ )
	  temp += mysqr( mt[ i*2+t ] - mt[ j*2+t ] );
	fval += mysqr( d_mds[ k ] - sqrt( temp ) );
	k++;
      }

  }
  template<typename Real>
  void ColaInit( int npar, int& mfct, Real& answer, Real* x,
		 Real* lo, Real* hi ) {

    if ( 17 != npar )
      throw std::runtime_error( "npar for the Cola func must be 17\n" );

    lo[ 0 ] = 0.0;
    hi[ 0 ] = 4.0;
    x[ 0 ] = 2.0;
    for ( int ii = 1; ii < npar; ++ii ) {
      lo[ ii ] = -4.0;
      x[ ii ] = 0.0;
      hi[ ii ] = 4.0;
    }

    answer = 12.8154;

  }

  template<typename Real, typename Type>
  void Colville( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the Colville func must be 4\n" );

    fval = 100 * pow(x[0]-x[1]*x[1], 2) + pow(1-x[0], 2)
      + 90 * pow(x[3] -x[2]*x[2], 2) + pow(1-x[2], 2)
	   + 10.1 * (pow(x[1] -1, 2) +  pow(x[3] - 1, 2))
	   + 19.8 * (x[1] - 1)*(x[3] - 1);

  }
  template<typename Real>
  void ColvilleInit( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the Colville func must be 4\n" );

    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 0.0;
      hi[ ii ] = 10.0;
      x[ ii ] = 8.5;
    }

    answer = 0.0;

  }
  //
  // f( 1, 1, 1, 1 ) = 0.0
  //

  //
  //     Deflected Corrugated Spring function
  //     dcs(5, 5, ..., 5) = -1 for any k and alpha=5 and npar
  //
  template<typename Real, typename Type>
  void dcs( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

  int k=5;
  Real alpha=5.0;
  std::vector< Real> c( npar, alpha );
  Real r2 = 0.0;
  for ( int ii = 0; ii < npar; ++ii )
    r2 += pow( x[ ii ] - c[ ii ], 2.0 );
  Real r = sqrt( r2 );
  fval = - cos( k * r ) + 0.1 * r2;

  }
  template<typename Real>
  void dcsInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    srand(1);
    Real lolo = -1.0e-3, hihi = 1.0e3;
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = lolo;
      hi[ ii ] = hihi;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lolo + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = -1.0;
  }

  //
  // decanom( 2, -3, ) = 0.0
  //
  template<typename Real, typename Type>
  void decanom( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    double x0 = x[ 0 ];
    double x1 = x[ 1 ];

    double f1 = fabs( (((((((((1.0*x0-20)*x0+180)*x0-960)*x0+3360)*x0-8064)*x0+
			  13340)*x0-15360)*x0+11520)*x0-5120)*x0+2624 );
    double f2 = fabs( (((1.0*x1+12.0)*x1+54)*x1+108)*x1+81.0 );

    fval = 0.0010 * pow( f1 + f2, 2.0 );

  }
  template<typename Real>
  void decanomInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {

    if ( npar != 2 )
      throw std::runtime_error( "npar for the dodeccal func must be 2\n" );
    Real low=-1.0e3, high=1.0e3;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }


  //
  // dodecal(1,2,3) = 0.
  //
  template<typename Real, typename Type>
  void dodecal( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 3 != npar ) {
      ierr = EXIT_FAILURE;
      return;
    }
    double f1 = 2*pow(x[0],3)+5*x[0]*x[1]+4*x[2]-2*x[0]*x[0]*x[2]-18;
    double f2 = x[0]+x[1]*x[2]+x[0]*x[1]*x[1]+x[0]*x[2]*x[2]-22.;
    double f3 = 8*x[0]*x[0]+2*x[1]*x[2]+2*x[1]*x[1]+3*pow(x[1],3)-52.;
    fval = pow(f1*f3*f2*f2+f1*f2*f3*f3+f2*f2+pow(x[0]+x[1]-x[2],2), 2 );
  }
  template<typename Real>
  void dodecalInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
      if ( npar != 3 )
      throw std::runtime_error( "npar for the dodeccal func must be 3\n" );
    Real low=-5, high=5;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }

  template<typename Real, typename Type>
  void DiscreteBoundary( int mfct, int npar, Real* x, Real* fvec,
			 int& ierr, Type xptr ) {

    Real h = 1.0 / ( npar + 1 );
    for ( int ii = 1; ii <= npar; ++ii ) {
      Real ti = ii * h;
      Real xiiminus1 = ( ii == 1 ) ? 0.0 :  x[ ii - 1 ];
      Real xiiplus1  = ( ii == npar ) ? 0.0 :  x[ ii ];

      fvec[ ii - 1 ] = 2.0 * x[ ii - 1 ]  - xiiminus1 - xiiplus1 +
	h*h * pow( x[ ii - 1 ] + ti + 1.0, 3.0 ) / 2.0;
    }

  }
  template<typename Real, typename Type>
  void DiscreteBoundary( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    DiscreteBoundary( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void DiscreteBoundaryInit( int npar, int& mfct, Real& answer,
			     Real* x, Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii ) {
      Real tj = ii / ( npar + 1.0 );
      x[ ii ] = tj * ( tj - 1.0 );
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }


  template<typename Real, typename Type>
  void DiscreteIntegral( int mfct, int npar, Real* x, Real* fvec,
			 int& ierr, Type xptr ) {

    const Real h = 1.0 / ( npar + 1 );
    const Real h2 = 0.5 * h;

    for ( int j = 1; j <= npar; j++ )
      fvec[ j - 1 ] = x[ j - 1 ];

    for ( int j = 1; j <= npar; j++ ) {
      Real tj = j * h;
      Real oj = 1.0 - tj;
      Real cj = h2 * pow( x[ j - 1 ] + tj + 1.0, 3.0 );
      for ( int i = 1; i <= npar; i++ ) {
	Real ti = i * h;
	Real oi = 1.0 - ti;
	if ( j <= i )
	  fvec[ i - 1 ] += oi * tj * cj;
	if ( j > i )
	  fvec[ i - 1 ] += ti * oj * cj;
      }
    }

  }
  template<typename Real, typename Type>
  void DiscreteIntegral( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    DiscreteIntegral( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void DiscreteIntegralInit( int npar, int& mfct, Real& answer, Real* x,
			     Real* lo, Real* hi ) {

    mfct = npar;

    for ( int ii = 0; ii < npar; ++ii ) {
      Real tj = ii / ( npar + 1.0 );
      x[ ii ] = tj * ( tj - 1.0 );
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }

  template<typename Real, typename Type>
  void DixonPrice( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    fval = 0.0;
    for ( int ii = 1; ii < npar; ++ii )
      fval += ( ii + 1 ) * mysqr( 2.0 * mysqr( x[ ii ] ) - x[ ii - 1 ] );

    fval += mysqr( x[ 0 ] - 1.0 );

  }
  template<typename Real>
  void DixonPriceInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = - 10.0;
      hi[ ii ] =   10.0;
      x[ ii ] = 2.5;
    }
    answer = 0.0;
  }

  template<typename Real, typename Type>
  void Easom( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = - cos(x[ 0 ]) * cos(x[ 1 ]) *
      exp( - mysqr(x[ 0 ]-pi) - mysqr( x[ 1 ] - pi ) );
  }
  template<typename Real>
  void EasomInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = - 100.0;
      hi[ ii ] =   100.0;
      x[ ii ] = 25.0;
    }
    answer = -1.0;
  }
  //
  // f( pi, pi ) = -1.0
  //

  //
  // factor1( 1, 2, 3, ..., npar ) = 0.0
  //
  template<typename Real, typename Type>
  void factor( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    fval = 0.0;
    Real p=1.0, fact=1.0;
    for ( int ii = 0; ii < npar; ++ii ) {
      p *= x[ ii ];
      fact *= ii + 1;
      fval += std::fabs( std::pow( p - fact, 2.0 ) );
    }

  }
  template<typename Real>
  void factorInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    Real low=-1.0e3, high=1.0e3;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }


  template<typename Real, typename Type>
  void FreudensteinRoth( int mfct, int npar, Real* x, Real* fvec,
			 int& ierr, Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 2 ) {
      fvec[ ii ] = - 13.0 + x[ ii ] +
	x[ ii + 1 ] * (( 5.0 - x[ ii + 1 ] ) * x[ ii + 1 ] - 2.0 );
      fvec[ ii + 1 ] = - 29.0 + x[ ii ] +
	 x[ ii + 1 ] * (( x[ ii + 1 ] + 1.0 ) * x[ ii + 1 ] - 14.0 );
    }

  }
  template<typename Real, typename Type>
  void FreudensteinRoth( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    FreudensteinRoth( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void FreudensteinRothInit( int npar, int& mfct, Real& answer, Real* x,
			     Real* lo, Real* hi ) {


    mfct = npar;
    for ( int ii = 0; ii < npar; ii += 2 ) {
      x[ ii ] = 0.5;
      x[ ii + 1 ] = -2.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e3;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e3;

    // The global min is 0.0, there is a local min at (npar/2 * 48.9842);
    answer = 0.0;

  }
  //
  // For n = 2,
  // this function has the global minimum f( 4, 5 ) = 0
  // and a local minimum at f( 11.4125; -0.89685 ) = 48.98425
  //    answer = (npar/2 * 48.9842);
  //

  template<typename Real, typename Type>
  void Func1( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    fval = sqrt( fabs( cos( sqrt( fabs(x[0]*x[0] + x[1] ) ) ) ) ) +
      0.01*( x[0] + x[1] );
  }
  template<typename Real>
  void Func1Init( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi, Real low=-1.0e1, Real high=1.0e1 ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Func1 must be == 2\n" );
    mfct = npar;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = low + (high - low) * rand( ) / ( RAND_MAX + 1.0 );
    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = low;
    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = high;
    //     FUNCTION #1 MIN AT -0.18467 APPROX AT (-8.4666, -10) APPROX
    answer = -0.18467;
  }

  template<typename Real, typename Type>
  void Gaussian( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		 Type xptr ) {

    Real yi[] = { 0.0009, 0.0044, 0.0175, 0.0540, 0.1295, 0.2420, 0.3521,
		  0.3989, 0.3521, 0.2420, 0.1295, 0.0540, 0.0175, 0.0044,
		  0.0009 };

    for ( int i = 0; i < 15; i++ ) {
      Real tmp2 = pow( ( 7.0 - i ) / 2.0 - x[ 2 ], 2.0 );
      fvec[ i ] = x[ 0 ] * exp( - x[ 1 ] * 0.5 * tmp2 ) - yi[ i ];
    }

  }
  template<typename Real, typename Type>
  void Gaussian( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 15;
    std::vector< Real> fvec( mfct );
    Gaussian( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void GaussianInit( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 3 != npar )
      throw std::runtime_error( "npar for the Gaussian func must be 3\n" );

    mfct = 15;
    x[ 0 ] = 0.4;
    x[ 1 ] = 1.0;
    x[ 2 ] = 0.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e5;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e5;

    answer = 1.12793e-8;

  }

  template<typename Real, typename Type>
  void GoldsteinPrice( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    Real x0 = x[0];
    Real x1 = x[1];
    Real x0x1 = x0 * x1;
    Real x1x1 = x1 * x1;
    Real a1 = x0 + x1 + 1;
    Real a12 = mysqr( a1 );

    Real a2 = 19 - 14 * x0 + 3 * x0 * x0 - 14 * x1 + 6 * x0x1 + 3 * x1x1;
    Real b1 = 2 * x0 - 3 * x1;
    Real b12 = mysqr( b1 );
    Real b2 = 18 - 32 * x0 + 12 * x0 * x0 + 48 * x1 - 36 * x0x1 +
      27 * x1x1;
    fval = (1 + a12 * a2) * (30 + b12 * b2);
  }
  template<typename Real>
  void GoldsteinPriceInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the GoldsteinPricefunc must be equal to 2\n" );

    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -2.0;
      hi[ ii ] = 2.0;
    }
    x[ 0 ] = -1.5;
    x[ 1 ] = 1.4;
    answer = 3.0;

  }
  //
  // f( ?, ? ) = 3.0
  //

  template<typename Real, typename Type>
  void Griewank( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    Real sum=0.0,prod=1.0;
    for (int ii = 0 ; ii < npar; ++ii ) {
      sum += x[ ii ] * x[ ii ] ;
      prod*=cos( x[ ii ] / sqrt( Real( ii + 1 ) ) );
    }
    fval = sum / 4000.0 - prod + 1.0 ;
  }
  template<typename Real>
  void GriewankInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    mfct = 0;
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -600.0;
      hi[ ii ] =  600.0;
      x[ ii ] = 100.0 * pow( -1.0, ii );
    }
    // The global minima f( 0, 0, ...., 0 ) = 0
    answer = 0.0;
  }

  template<typename Real, typename Type>
  void GulfResearchDevelopment( int mfct, int npar, Real* x, Real* fvec,
				int& ierr, Type xptr ) {

    Real two3rd = 2.0 / 3.0;
    for ( int i = 0; i < mfct; i++ ) {
      Real ti = ( i + 1 ) * 0.01;
      Real yi = 25.0 + pow( -50.0 * log( ti ), two3rd );
      Real ai = yi - x[ 1 ];
      fvec[ i ] = exp( -(pow( fabs(ai), x[ 2 ] ) / x[ 0 ]) ) - ti;
    }

  }
  template<typename Real, typename Type>
  void GulfResearchDevelopment( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    GulfResearchDevelopment( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void GulfResearchDevelopmentInit( int npar, int& mfct, Real& answer,
				    Real* x, Real* lo, Real* hi ) {

    if ( 3 != npar )
      throw std::runtime_error( "npar for the GulfResearchDevelopment func must be equal to 3\n" );

    mfct = npar;
    x[ 0 ] = 5.0;
    x[ 1 ] = 2.5;
    x[ 2 ] = 0.15;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e3;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e3;

    answer = 0.0;
  }
  //
  // This function has the global minimum f( 50, 25, 1.5 ) = 0.0
  // and a local minimum at
  // f( 99.89537833835886, 60.61453903025014, 9.16124389236433 ) = 0.038
  // Additionally, this local minimum is very flat and the minimizer is
  // surrounded by a "plateau", where the gradient is zero everywhere and
  // the function equals 0.0385.
  // Another local minimum with exactly the same function value (0.0380) is
  // reached at the point
  // f( 201.662589489426; 60.616331504682; 10.224891158489 ) = 0.0380
  //


  template<typename Real>
  Real mytheta( Real x1, Real x2 ) {

    if ( 0.0 == x1 )
      return 1.0e128;

    Real piinv = 0.25 / atan( 1.0 );
    Real quo   = x2 / x1;
    Real theta = 0.5 * piinv * atan( quo );


    if ( x1 < 0.0 )
      theta += 0.5;

    return theta;

  }

  //
  // compute sqrt( a^2 + b^2 ) without destructive underflow or overflow
  //
  template<typename Real>
  Real mypythag( Real a, Real b ) {

    Real absa = fabs( a );
    Real absb = fabs( b );
    if ( absa > absb )
      return absa * sqrt( 1.0 + mysqr( absb /absa ) );
    else {
      if ( absb == 0.0 )
	return 0.0;
      else
	absb * sqrt( 1.0 + mysqr( absa / absb ) );
    }

  }

  template<typename Real, typename Type>
  void Hansen( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 2 != npar )
      throw std::runtime_error( "npar for the Hansen func must be 2\n" );

    fval = (cos(1.0)+2.0*cos(x[0]+2.0)+3.0*cos(2.0*x[0]+3.0)+
	     4.0*cos(3.0*x[0]+4.0)+5.0*cos(4.0*x[0]+5.0))*
      (cos(2.0*x[1]+1.0)+2.0*cos(3.0*x[1]+2.0)+
       3.0*cos(4.0*x[1]+3.0)+4.0*cos(5.0*x[1]+4.0)+5.0*cos(6.0*x[1]+5.0));

  }
  template<typename Real>
  void HansenInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Hansen func must be 2\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -10.0;
      hi[ ii ] = 10.0;
      x[ ii ] = 3.5;
    }
    answer = -176.54;
  }

  template<typename Real, typename Type>
  void Hartman6( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    if ( 6 != npar )
      throw std::runtime_error( "npar for the Hartman6 func must be 6\n" );

    static double a[4][6] = {
      {10,  3,  17,   3.5,  1.7,  8.0},
      {0.05,  10, 17, 0.1, 8,  14},
      {3, 3.5, 1.7, 10, 17, 8},
      {17, 8, 0.05, 10, 0.1, 14}
    };

    static double c[] = {1, 1.2, 3, 3.2};
    static double p[4][6] = {
      {0.1312, 0.1696, 0.5569, 0.0124, 0.8283, 0.5886},
      {0.2329, 0.4135, 0.8307, 0.3736, 0.1004, 0.9991},
      {0.2348, 0.1415, 0.3522, 0.2883, 0.3047, 0.6650},
      {0.4047, 0.8828, 0.8732, 0.5743, 0.1091, 0.0381}
    };

    double s = 0.0;
    for( int ii = 0; ii < 4; ++ii ) {
      double t = 0.0;
      for( int jj = 0; jj < npar; ++jj ) {
	double t1 = x[ jj ]- p[ ii ][ jj ];
	t += a[ ii ][ jj ] * ( t1 * t1 );
      }
      s += c[ ii ] * exp( - t  );
    };
    fval = -s;
  }
  template<typename Real>
  void Hartman6Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {
    if ( 6 != npar )
      throw std::runtime_error( "npar for the Hartman6 func must be 6\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 0.0;
      hi[ ii ] = 1.0;
      x[ ii ] = 0.5;
    }
    answer = -3.32;
  }
  //
  // f( 0.201, 0.150, 0.477, 0.275, 0.311, 0.657 ) = -3.32
  //

  template<typename Real, typename Type>
  void HelicalValley( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		      Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 3 ) {
      fvec[ ii ] = 10.0 * ( x[ ii + 2 ] -
			    10 * mytheta( x[ ii ], x[ ii + 1 ] ) );
      fvec[ ii + 1 ] = 10.0 * ( sqrt( x[ ii ] * x[ ii ] +
				      x[ ii + 1 ] * x[ ii + 1 ] ) - 1.0 );
      fvec[ ii + 2 ] = x[ ii + 2 ];
    }

  }
  template<typename Real, typename Type>
  void HelicalValley( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    fval = 0.0;
    for ( int ii = 0; ii < npar; ii += 3 ) {
      Real a = 10.0 * ( x[ ii + 2 ] -
			    10 * mytheta( x[ ii ], x[ ii + 1 ] ) );
      Real b = 10.0 * ( sqrt( x[ ii ] * x[ ii ] +
				      x[ ii + 1 ] * x[ ii + 1 ] ) - 1.0 );
      Real c = x[ ii + 2 ];

      fval += a*a + b*b + c*c;
    }

  }
  template<typename Real>
  void HelicalValleyInit( int npar, int& mfct, Real& answer, Real* x,
			  Real* lo, Real* hi ) {

    if ( npar % 3 )
      throw std::runtime_error( "npar for the HelicalValley func must be multiple of 3\n" );

    mfct = npar;

    for ( int ii = 0; ii < npar; ii += 3 ) {
      x[ ii ] = -1.0;
      x[ ii + 1 ] = 0.0;
      x[ ii + 2 ] = 0.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // f( 1, 0, 0 ) = 0.0
  //


  template<typename Real, typename Type>
  void Himmelblau( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Himmelblau func must be 2\n" );
    fval = mysqr( x[ 0 ] * x[ 0 ]+x[ 1 ] - 11.0 ) +
      mysqr( x[ 0 ]+x[ 1 ]* x[ 1 ] - 7 );
  }
  template<typename Real>
  void HimmelblauInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Himmelblau func must be 2\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -5.0;
      hi[ ii ] = 5.0;
      x[ ii ] = 0.0;
    }
    answer = -0.0;
  }
  //
  // f( 3, 2 ) = 0.0
  //

  template<typename Real, typename Type>
  void Holzman1( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 3 != npar )
      throw std::runtime_error( "npar for the Holzman1 func must be 3\n" );
    fval = 0.0;
    Real twothird = 2.0 / 3.0;
    for ( int ii = 0; ii < 99; ++ii ) {
      double ui = 25 + pow( -50.0 * log( 0.01 *( ii + 1 ) ), twothird );
      fval += -0.1 * ( ii + 1 ) + exp( pow( ui - x[ 1 ], x[ 2 ] ) / x[ 0 ] );
    }
    return;
  }
  template<typename Real>
  void Holzman1Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {
    if ( 3 != npar )
      throw std::runtime_error( "npar for the Holzman1 func must be 3\n" );
    lo[ 0 ] = 0.1; hi[ 0 ] = 100.0;
    lo[ 1 ] = 0.0; hi[ 1 ] = 25.6;
    lo[ 2 ] = 0.0; hi[ 2 ] = 5.0;
    // for ( int ii = 0; ii < npar; ++ii )
    //   x[ ii ] = 5.0;
    x[0] = 50; x[1] = 25; x[2] = 1.5;
    answer = -395.005;
  }
  //
  // f( 50, 25, 1.5 ) = 0.0;
  //

  template<typename Real, typename Type>
  void Holzman2( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    fval = 0.0;
    for ( int ii = 0; ii < npar; ++ii ) {
      fval += ii * pow( x[ ii ] , 4);
    }
  }
  template<typename Real>
  void Holzman2Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {
    for ( int ii = 0; ii < npar; ++ii ) {
      x[ ii ] = 5.0;
      lo[ ii ] = -10.0;
      hi[ ii ] = 10.0;
    }
    answer = 0.0;
  }
  //
  // f( 0, 0, ..., 0 ) = 0.0
  //

  template<typename Real, typename Type>
  void JennrichSampson( int mfct, int npar, Real* x, Real* fvec, int& ierr,
			Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 2 )
      for ( int i = 1; i <= 10; i++ )
	fvec[ 10*ii/2 + i - 1 ] = 2 * ( 1 + i ) -
	  ( exp( i * x[ ii ] ) + exp( i * x[ ii + 1 ] ) );

  }
  template<typename Real, typename Type>
  void JennrichSampson( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    fval = 0.0;
    for ( int ii = 0; ii < npar; ii += 2 )
      for ( int i = 1; i <= 10; i++ ) {
	Real foo = 2 * ( 1 + i ) -
	  ( exp( i * x[ ii ] ) + exp( i * x[ ii + 1 ] ) );
	fval += foo*foo;
      }

  }
  template<typename Real>
  void JennrichSampsonInit( int npar, int& mfct, Real& answer, Real* x,
			    Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the JennrichSampson func must be even\n" );

    mfct = 10 * npar/2;
    for ( int ii = 0; ii < npar; ii += 2 ) {
      x[ ii ] = 0.3;
      x[ ii + 1 ] = 0.4;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e5;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e5;

    answer = (npar/2 * 124.362);
  }
  //
  // This function has the global minimum f( 0.257825, 0.257825 ) = 124.362
  // for m = 10, and a local minimum f=259.58019078047 on the plane
  // x1=0.33148432
  //


  //
  // From The Theory And Practice Of Econometrics, 2nd ed., pp. 956-7.
  // Judge et al: judge(0.86479,1.2357)=16.0817307, the global minimum
  // judge(2.35,-0.319)=20.9805, the  local minimum
  //
  template<typename Real, typename Type>
  void Judge( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    static const double y[] = {4.284,4.149,3.877,0.533,2.211,2.389,2.145,
			       3.231,1.998,1.379,2.106,1.428,1.011,2.179,
			       2.858,1.388,1.651,1.593,1.046,2.152 };
    static const double x2[] = {.286,.973,.384,.276,.973,.543,.957,.948,.543,
				.797,.936,.889,.006,.828,.399,.617,.939,.784,
				.072,.889};
    static const double x3[] = {.645,.585,.310,.058,.455,.779,.259,.202,.028,
				.099,.142,.296,.175,.180,.842,.039,.103,.620,
				.158,.704};
    const int size = sizeof(x3)/sizeof(double);
    fval = 0.0;
    for ( int ii = 0; ii < size; ++ii )
      fval += pow(x[0] + x[1]*x2[ii] + (x[1]*x[1])*x3[ii] - y[ii], 2.0);

  }
  template<typename Real>
  void JudgeInit( int npar, int& mfct, Real& answer, Real* x, Real* lo,
		Real* hi, Real low=-1.0e3, Real high=1.0e3 ) {
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 16.0817307;
  }

  // template<typename Real, typename Type>
  // void Katsuuras( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
  //   int d = 32;
  //   fval = 1.0;
  //   for ( int ii = 0; ii < npar; ++ii ) {
  //     double s = 0.0;
  //     for ( int kk = 1; kk <= d; ++kk ) {
  //       double pow2 = std::pow( Real(2), Real(kk) );
  //       s += round(pow2 * x[ ii ]) / pow2;
  //     }
  //     fval *= 1.0 + ( ii + 1) * s;
  //   }
  // }
  // template<typename Real>
  // void KatsuurasInit( int npar, int& mfct, Real& answer, Real* x, Real* lo,
  //       	      Real* hi, Real low=-1.0e3, Real high=1.0e3 ) {
  //   srand(1);
  //   for ( int ii = 0; ii < npar; ++ii ) {
  //     lo[ ii ] = low;
  //     hi[ ii ] = high;
  //     double tmp = rand( ) / (RAND_MAX + 1.0);
  //     x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
  //   }
  //   answer = 1.0;
  // }
  // //
  // // f( 0, 0, ..., 0 ) = 1.0
  // //

  template<typename Real, typename Type>
  void KowalikOsborne( int mfct, int npar, Real* x, Real* fvec,
		       int& ierr, Type xptr ) {

    Real yi[] = { 0.1957, 0.1947, 0.1735, 0.1600, 0.0844, 0.0627, 0.0456,
		  0.0342, 0.0323, 0.0235, 0.0246 };
    Real ui[] = { 4.0000, 2.0000, 1.0000, 0.5000, 0.2500, 0.1670, 0.1250,
		  0.1000, 0.0833, 0.0714, 0.0625 };

    for ( int i = 0; i < mfct; i++ )
      fvec[ i ] = yi[ i ] - x[ 0 ] * ui[ i ] *
	( ui[ i ] + x[ 1 ] ) / ( ui[ i ] * ( ui[ i ] + x[ 2 ] ) + x[ 3 ] );

  }
  template<typename Real, typename Type>
  void KowalikOsborne( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 11;
    std::vector< Real> fvec( mfct );
    KowalikOsborne( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void KowalikOsborneInit( int npar, int& mfct, Real& answer, Real* x,
			   Real* lo, Real* hi ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the KowalikOsborne func must be 4\n" );

    mfct = 11;
    x[ 0 ] = 0.25;
    x[ 1 ] = 0.39;
    x[ 2 ] = 0.415;
    x[ 3 ] = 0.39;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 3.07505e-4;

  }
  //
  // This function has the global minimum
  // f( .192807, 0.191282, 0.123056, 0.136062 ) = 0.000307506
  // and a local minimum at
  // f ( inf, -14.07..., -inf, -inf ) = 1.02734...10^-3
  // Additionally, it takes a local minimum f=0.00179453906640 on a big flat
  // plateau. The following are the minimum points found on this plateau by
  // SolvOpt:
  // [0.816961, 2045080438, 7189803671, 2927176299],
  // [0.734588, 2648629959, 8372807699, 3408733052],
  // [0.537138, 1852002979, 4280874360, 1742904647].
  // These particular points appear to be "false" local minima on a flat
  // plateau, the fact caused by the finite accuracy of calculations.
  //

  template<typename Real, typename Type>
  void Levy( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = 0.0;
    for (int ii = 0 ; ii <= npar - 2 ; ++ii )
      fval += pow( x[ ii ] - 1.0, 2.0) *
	( 1.0 + pow( sin( 3.0 * pi * x[ ii + 1 ]), 2.0 ) );
    fval += pow( sin( 3 * pi * x[ 0 ] ), 2.0 ) +
      ( x[ npar-1 ] - 1.0 ) *
      ( 1.0 + pow( sin( 2.0 * pi * x[ npar - 1 ]), 2.0 ) );

  }
  template<typename Real>
  void LevyInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 4 == npar )
      for ( int ii = 0; ii < npar; ++ii ) {
	lo[ ii ] = - 10.0;
	hi[ ii ] =   10.0;
	x[ ii ] = 7.0;
      }
    else
      for ( int ii = 0; ii < npar; ++ii ) {
	lo[ ii ] = - 5.0;
	hi[ ii ] =   5.0;
	x[ ii ] = 4.0;
      }

    switch( npar ) {

    case 4:
      answer = -21.502;
      break;
    case 5:
    case 6:
    case 7:
      answer = -11.504;
      break;
    default:
      answer = std::numeric_limits< double >::max( );
      return;

    }

  }

  template<typename Real, typename Type>
  void LinearFullRank( int mfct, int npar, Real* x, Real* fvec,
		       int& ierr, Type xptr ) {

    Real sum = 0.0;
    for ( int j = 1; j <= npar; j++ )
      sum += x[ j - 1 ];

    const Real dnfi2 = 2.0 / mfct;
    for ( int i = 1; i <= mfct; i++ ) {
      if (i <= npar)
	fvec[ i - 1 ] = x[ i - 1 ] - dnfi2*sum - 1.0;
      if (i > npar)
	fvec[ i - 1 ] =      - dnfi2*sum - 1.0;
    }

  }
  template<typename Real, typename Type>
  void LinearFullRank( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    LinearFullRank( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void LinearFullRankInit( int npar, int& mfct, Real& answer, Real* x,
			   Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = static_cast< Real>( mfct - npar );
  }
  //
  // f( -1, ... , -1) = m - n
  //

  template<typename Real, typename Type>
  void LinearFullRank1( int mfct, int npar, Real* x, Real* fvec,
			int& ierr, Type xptr ) {

    Real sum = 0.0;
    for ( int jj = 1; jj <= npar; ++jj )
      sum += x[ jj - 1 ] * jj;

    for ( int ii = 1; ii <= mfct; ++ii )
      fvec[ ii - 1 ] = ii * sum - 1.0;

  }
  template<typename Real, typename Type>
  void LinearFullRank1( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    LinearFullRank1( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void LinearFullRank1Init( int npar, int& mfct, Real& answer, Real* x,
			    Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer =
      static_cast< Real>( mfct * ( mfct - 1.0 ) / 2.0 / ( 2.0 * mfct + 1.0 ) );

  }

  template<typename Real, typename Type>
  void LinearFullRank0cols0rows( int mfct, int npar, Real* x, Real* fvec,
				 int& ierr, Type xptr ) {

    const int nxm1   = npar - 1;
    const int nfm1   = mfct - 1;
    Real sum = 0.0;
    for ( int j = 2; j <= nxm1; j++ )
      sum += x[ j - 1 ] * j;

    fvec[ 0 ] = - 1.0;
    for ( int i = 2; i <= nfm1; i++ )
      fvec[ i - 1 ] = ( i - 1 ) * sum - 1.0;

    fvec[ mfct - 1 ] = - 1.0;

  }
  template<typename Real, typename Type>
  void LinearFullRank0cols0rows( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    LinearFullRank0cols0rows( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void LinearFullRank0cols0rowsInit( int npar, int& mfct, Real& answer,
				     Real* x, Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = static_cast< Real>( ( mfct * ( mfct + 3 ) - 6.0 ) / 2.0 /
				  ( 2.0 * mfct - 3.0 ) );

  }

  template<typename Real, typename Type>
  void McCormick( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    Real a = x[ 0 ] + x[ 1 ], b = x[ 0 ] - x[ 1 ];
    fval = sin( a ) + b * b - 1.5 * x[ 0 ] + 2.5 * x[ 1 ] + 1.0;
  }
  template<typename Real>
  void McCormickInit( int npar, int& mfct, Real& answer, Real* x,
		      Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the McCormick func must be 2\n" );
    x[ 0 ] = 0.0; x[ 1 ] = 0.0;
    lo[ 0 ] = -1.5; lo[ 1 ] = -3.0;
    hi[ 0 ] =  4.0; hi[ 1 ] =  4.0;
    answer = -1.9132;
  }
  //
  // f ( -0.547197553, -1.54719756 ) = -1.91
  //

  template<typename Real, typename Type>
  void McKinnon( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real tau = 3.0, theta = 6.0, phi = 400.0;
    const Real tmp = 1.0 + x[ 1 ];
    if ( x[ 0 ] <= 0.0 )
        fval = theta * phi * pow( std::fabs( x[ 0 ] ), tau ) * tmp;
    else
	fval = theta * pow( x[ 0 ], tau ) + tmp;
  }
  template<typename Real>
  void McKinnonInit( int npar, int& mfct, Real& answer, Real* x,
		      Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the McKinnon func must be 2\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      x[ ii ] = 1.0;
      lo[ ii ] = -1.0e2;
      hi[ ii ] = 1.0e2;
    }
    answer = - 0.25;
  }
  //
  // f ( 0.0, -0.5 ) = - 0.25
  //

  template<typename Real, typename Type>
  void Meyer( int mfct, int npar, Real* x, Real* fvec, int& ierr,
	      Type xptr ) {

    Real yi[] = { 34780.0, 28610.0, 23650.0, 19630.0, 16370.0, 13720.0,
		  11540.0, 9744.0, 8261.0, 7030.0, 6005.0, 5147.0, 4427.0,
		  3820.0, 3307.0, 2872.0 };

    int ii = 0;
    for ( int jj = 0; jj < 16; jj++ )
      fvec[ 16 * ii/3 + jj ] = x[ ii ] *
	exp( x[ ii + 1 ] / ( 45.0 + 5.0 * (jj+1) + x[ ii + 2 ] ) ) - yi[ jj ];

  }
  template<typename Real, typename Type>
  void Meyer( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 16;
    std::vector< Real> fvec( mfct );
    Meyer( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void MeyerInit( int npar, int& mfct, Real& answer, Real* x,
		  Real* lo, Real* hi ) {

    if ( 3 != npar )
      throw std::runtime_error( "npar for the Meyer func must be 3\n" );

    mfct = 16 * npar / 3;
    for ( int ii = 0; ii < npar; ii += 3 ) {
      x[ ii ] = 0.02;
      x[ ii + 1 ] = 4000.0;
      x[ ii + 2 ] = 250.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e3;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e5;

    answer = 87.9458;
  }
  template<typename Real, typename Type>
  void Michalewicz( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = 0.0;
    for ( int ii = 0; ii < npar; ++ii )
      fval -= sin( x[ ii ] ) *
	pow( sin( (ii+1) * x[ ii ] * x[ ii ] / pi ), 20.0 );
  }
  template<typename Real>
  void MichalewiczInit( int npar, int& mfct, Real& answer, Real* x,
			Real* lo, Real* hi ) {
    switch( npar ) {
    case 2:
      answer = - 1.8013;
      break;
    case 5:
      answer = - 4.687658;
      break;
    case 10:
      answer = - 9.66015;
      break;
    default:
      printf( "MichalewiczInit(): npar=%d\t is not an option\n", npar );
      return;
    }
    const Real pi = std::fabs( acos( - 1.0 ) );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 0.0;
      hi[ ii ] = pi;
      x[ ii ] = pi / 2.0;
    }
  }

  template<typename Real, typename Type>
  void Osborne1( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		 Type xptr ) {

    Real yi[] = { 0.844, 0.908, 0.932, 0.936, 0.925, 0.908, 0.881, 0.850,
		  0.818, 0.784, 0.751, 0.718, 0.685, 0.658, 0.628, 0.603,
		  0.580, 0.558, 0.538, 0.522, 0.506, 0.490, 0.478, 0.467,
		  0.457, 0.448, 0.438, 0.431, 0.424, 0.420, 0.414, 0.411,
		  0.406 };

    for ( int i = 0; i < mfct; i++ ) {
      Real ti = 10.0 * i;
      fvec[ i ] = yi[ i ] - ( x[0] + x[1] * exp( -ti*x[3] ) +
			      x[2] * exp( -ti*x[4] ) );
    }

  }
  template<typename Real, typename Type>
  void Osborne1( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 33;
    std::vector< Real> fvec( mfct );
    Osborne1( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void Osborne1Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 5 != npar )
      throw std::runtime_error( "npar for the Osborne1Init func must be 5\n" );
    mfct = 33;
    x[ 0 ] = 0.5;
    x[ 1 ] = 1.5;
    x[ 2 ] = -1.0;
    x[ 3 ] = 0.01;
    x[ 4 ] = 0.02;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 5.46489e-5;

  }
  //
  // This function has the global minimum
  // f( 0.3754; 1.9358; -1.4647; 0.01287; 0.02212 ) = 5.46489e-5
  // It also has a number of "false" minima which are points of flatness
  // and look like a minimum.
  // In particular, gradient is approximately zero
  // ( [1.55e-15; 4.2e-11; 4.2e-11; 0; 0] ) everywhere on a big "plateau" at
  // x1=0.62415625, x2~125*10i, x3~-125*10i, where i=0,1,..., and x4>10, x5>10.
  // The function equals 1.10603621875 everywhere on this "plateau"
  // (for example, try the point
  // [0.62415625; 12500.1402308134; -12499.9203870634; 100; 200]).
  // As soon as the starting point is far from the global minimum, this flat
  // area is possibly the most attractive for a gradient method.
  //

  template<typename Real, typename Type>
  void Osborne2( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		 Type xptr ) {

    Real yi[] = { 1.366, 1.191, 1.112, 1.013, 0.991, 0.885, 0.831, 0.847,
		  0.786, 0.725, 0.746, 0.679, 0.608, 0.655, 0.615, 0.606,
		  0.602, 0.626, 0.651, 0.724, 0.649, 0.649, 0.694, 0.644,
		  0.624, 0.661, 0.612, 0.558, 0.533, 0.495, 0.500, 0.423,
		  0.395, 0.375, 0.372, 0.391, 0.396, 0.405, 0.428, 0.429,
		  0.523, 0.562, 0.607, 0.653, 0.672, 0.708, 0.633, 0.668,
		  0.645, 0.632, 0.591, 0.559, 0.597, 0.625, 0.739, 0.710,
		  0.729, 0.720, 0.636, 0.581, 0.428, 0.292, 0.162, 0.098,
		  0.054 };

    for ( int ii = 0; ii < mfct; ++ii ) {
      Real ti = 0.1 * ii;
      Real a = x[0] * exp( - ti * x[4] );
      Real b = x[1] * exp( - x[5] * (ti-x[8])*(ti-x[8]) );
      Real c = x[2] * exp( - x[6] * (ti-x[9])*(ti-x[9]) );
      Real d = x[3] * exp( - x[7] * (ti-x[10])*(ti-x[10]) );
      fvec[ ii ] = yi[ ii ] - ( a + b + c + d );
    }

  }
  template<typename Real, typename Type>
  void Osborne2( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 65;
    std::vector< Real> fvec( mfct );
    Osborne2( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void Osborne2Init( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 11 != npar )
      throw std::runtime_error( "npar for the Osborne2Init func must be 11\n" );

    mfct = 65;
    x[ 0 ] = 1.3;
    x[ 1 ] = 0.65;
    x[ 2 ] = 0.65;
    x[ 3 ] = 0.7;
    x[ 4 ] = 0.6;
    x[ 5 ] = 3.0;
    x[ 6 ] = 5.0;
    x[ 7 ] = 7.0;
    x[ 8 ] = 2.0;
    x[ 9 ] = 4.5;
    x[ 10 ] = 5.5;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e3;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e3;

    answer = 4.01377e-2;
  }
  //
  // f( 1.30997, 0.431459, 0.633631, 0.599303, 0.753915, 0.905575,
  //    1.36504, 4.82482, 2.39882, 4.56887, 5.67537 ) = 0.0401377
  //

  template<typename Real, typename Type>
  void Paviani( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    Real mul = 1.0;
    fval = 0.0;
    for ( int ii=0 ; ii <npar ; ++ii ) {
      Real a = log( x[ ii ] - 2.0 );
      Real b = log( 10.0 - x[ ii ] );
      fval += a * a +  b * b;
      mul *= x[ ii ];
    }
    fval -= pow( mul, 0.2 );
  }
  template<typename Real>
  void PavianiInit( int npar, int& mfct, Real& answer, Real* x,
		    Real* lo, Real* hi ) {
    if ( 10 != npar )
      throw std::runtime_error( "npar for the PavianiInit func must be 10\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 2.001;
      x[ ii ] = 5.0;
      hi[ ii ] = 9.99;
    }
    answer = -45.778;
  }
  //
  // f( 9.35026583, 9.35026583, ..., 9.35026583, 9.35026583 ) = -45.7;
  //

  template<typename Real, typename Type>
  void PenaltyI( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		 Type xptr ) {

    Real qtr = 1.0/4.0;
    Real sum = - qtr;
    Real a2  = sqrt( 1.0e-5 );
    for ( int ii = 1; ii <= npar; ++ii ) {
      Real xi   = x[ii-1];
      sum  += xi*xi;
      fvec[ii-1] = a2*(xi - 1.0);
    }

    fvec[npar] = sum;

  }
  template<typename Real, typename Type>
  void PenaltyI( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar + 1;
    std::vector< Real> fvec( mfct );
    PenaltyI( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void PenaltyIInit( int npar, int& mfct, Real& answer, Real* x,
		     Real* lo, Real* hi ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the PenaltyI func must be 4\n" );

    mfct = npar + 1;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = ii + 1;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 2.24997e-5;

  }


  template<typename Real, typename Type>
  void PenaltyII( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		  Type xptr ) {

    Real sum = 0.0;
    Real sqrta = sqrt( 1.0e-5 );
    Real exp01 = exp( - 0.1 );

    fvec[ 0 ] = x[ 0 ] - 0.2;

    for ( int ii = 2; ii <= npar; ++ii ) {
      Real yi = exp( ii * 0.1 ) + exp( (ii-1) * 0.1 );
      fvec[ ii - 1 ] = sqrta * ( exp( x[ ii - 1 ] * 0.1 ) +
				exp( x[ ii - 2 ] * 0.1 ) - yi );

    }

    for ( int ii = npar + 1; ii <= 2 * npar - 1; ++ii ) {
      fvec[ ii - 1 ] = sqrta * ( exp( 0.1 * x[ ii - npar ] ) - exp01 );
    }

    for ( int ii = 1; ii <= npar; ++ii )
      sum += ( npar - ii + 1 ) * x[ ii - 1 ] * x[ ii - 1 ];

    fvec[ 2 * npar - 1 ] = sum - 1.0;


  }
  template<typename Real, typename Type>
  void PenaltyII( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 2 * npar;
    std::vector< Real> fvec( mfct );
    PenaltyII( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void PenaltyIIInit( int npar, int& mfct, Real& answer, Real* x,
		      Real* lo, Real* hi ) {

    if ( 4 != npar )
      throw std::runtime_error( "npar for the PenaltyII func must be 4\n" );

    mfct = 2 * npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 0.5;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 9.37629e-6;

  }
  //
  // f=9.37629...10^(-6)  if n=4
  // f=2.93660...10^(-4)  if n=10
  //

  template<typename Real, typename Type>
  void PowellBadlyScaled( int mfct, int npar, Real* x, Real* fvec,
			  int& ierr, Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 2 ) {
      fvec[ ii ] = 1.0e4 * x[ ii ] * x[ ii + 1 ] - 1;
      fvec[ ii + 1 ] = exp( - x[ ii  ] ) + exp( - x[ ii + 1 ] ) - 1.0001;
    }

  }
  template<typename Real, typename Type>
  void PowellBadlyScaled( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    PowellBadlyScaled( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void PowellBadlyScaledInit( int npar, int& mfct, Real& answer, Real* x,
			      Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the PowellBadlyScaled func must be even\n" );

    mfct = npar;
    for ( int ii = 0; ii < npar; ii += 2 ) {
      x[ ii ] = 0.0;
      x[ ii + 1 ] = 1.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // For n = 2, this function has a global minimum at
  // f( 1.09815933e-5, 9.106146738 ) = 0.0
  // t is very flat near the minimum point and a solution x with
  // f(x) <= 10-8 is not a satisfactory one in view of the argument.
  // A satisfactory solution providing sufficiently accurate estimate
  // for the minimizer is less then f(x)~10-13.
  //

  template<typename Real, typename Type>
  void PowellSingular( int mfct, int npar, Real* x, Real* fvec,
		       int& ierr, Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 4 ) {
      fvec[ ii ] = x[ ii ] + 10.0 * x[ ii + 1 ];
      fvec[ ii + 1 ] = sqrt( 5.0 ) * ( x[ ii + 2 ] - x[ ii + 3 ] );
      fvec[ ii + 2 ] = pow( x[ ii + 1 ] - 2.0 * x[ ii + 2 ], 2.0 );
      fvec[ ii + 3 ] = sqrt( 10.0 ) * pow( x[ ii ] - x[ ii + 3 ], 2.0 );
    }

  }
  template<typename Real, typename Type>
  void PowellSingular( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    PowellSingular( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void PowellSingularInit( int npar, int& mfct, Real& answer, Real* x,
			   Real* lo, Real* hi ) {

    if ( npar % 4 )
      throw std::runtime_error( "npar for the PowellSingular func must be multiple of 4\n" );

    mfct = npar;
    for ( int ii = 0; ii < npar; ii += 4 ) {
      x[ ii ] = 3.0;
      x[ ii + 1 ] = -1.0;
      x[ ii + 2 ] = 0.0;
      x[ ii + 3 ] = 1.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // f( origin ) = 0.0
  //


/*
  //
  // quintic( 2, 2,..., 2, -0.402627941 ) = 0.0
  // quintic( -1, -1,..., -1, -0.402627941 ) = 0.0
  // or any permutation of (2, 2,..., 2, -1, -1, ..., -1, -0.402627941)
  static void quintic( int npar, double* x, double& fval, int& ierr ) {

    fval = 0.0;
    for ( int ii = 0; ii < npar; ++ii )
      fval += std::fabs( ((((x[ii]-3)*x[ii]+4)*x[ii]+2)*x[ii]-10)*x[ii] - 4 );

  }

  template<typename Real>
  void quinticInit( int npar, Real& answer, Real* x, Real* lo,
		    Real* hi, Real low=-1.0e3, Real high=1.0e3 ) {
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }
*/


  template<typename Real, typename Type>
  void Rastrigin( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    const Real pi = std::fabs( acos( - 1.0 ) );
    fval = 0.0;
    for ( int ii = 0; ii < npar; ++ii )
      fval += ( mysqr( x[ ii ] ) - 10.0 * cos( 2 * pi * x[ ii ] ) );
    fval += 10.0 * npar;
  }
  template<typename Real>
  void RastriginInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = - 5.12;
      hi[ ii ] =   5.12;
      x[ ii ] = 3.14;
    }
    answer = 0.0;
  }
  //
  // f( 0, 0, ..., 0 ) = 0
  //

  template<typename Real, typename Type>
  void Rosenbrock( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		   Type xptr ) {

    for ( int ii = 0; ii < npar; ii += 2 ) {

      fvec[ ii ] = 1.0 - x[ ii ];
      fvec[ ii + 1 ] = 10.0 * ( x[ ii + 1 ] - x[ ii ] * x[ ii ] );

    }

  }
  template<typename Real, typename Type>
  void Rosenbrock( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    Rosenbrock( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void RosenbrockInit( int npar, int& mfct, Real& answer, Real* x,
		       Real* lo, Real* hi ) {

    if ( npar % 2 )
      throw std::runtime_error( "npar for the Rosenbrock func must be even\n" );

    mfct = npar;
    for ( int ii = 0; ii < npar; ii += 2 ) {
      x[ ii ] = -1.2;
      x[ ii + 1 ] = 1.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e2;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e2;

    answer = 0.0;

  }
  //
  // This well-known test function becomes very nasty as soon as a point
  // is far from the minimum. The function has a long ravine with very steep
  // walls and almost flat bottom. If x1>108 and x2=x12, the two successive
  // gradients are exactly opposite to each other at a very small possible
  // step size. This means, one cannot reach the bottom because of the finite
  // accuracy of computer calculations and one cannot calculate a right
  // direction to go along the ravine. This is why a computer application of
  // any gradient method may fail to minimize the function. The more so, if
  // gradients are approximated by the finite differences.
  //

  //
  // seqp(0,0) = seqp(2,2) = 0, find the values such that x0+x1 = x0*x1
  //
  template<typename Real, typename Type>
  void seqp( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    double x0 = x[0];
    double x1 = x[1];
    fval = std::fabs( (x0+x1) - (x0*x1) );
  }
  template<typename Real>
  void seqpInit( int npar, int& mfct, Real& answer, Real* x, Real* lo,
		 Real* hi, Real low=-1.0e3, Real high=1.0e3 ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the seqp func must be 2\n" );
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = low;
      hi[ ii ] = high;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = low + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }

  template<typename Real, typename Type>
  void Shekel5( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel5 func must be 4\n" );
    static const Real a[10][4] = { {4,4,4,4}, {1,1,1,1}, {8,8,8,8},
				   {6,6,6,6}, {3,7,3,7}, {2,9,2,9},
				   {5,5,3,3}, {8,1,8,1}, {6,2,6,2},
				   {7,3.6,7,3.6} };
    static const double c[10] = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3,
				  0.7, 0.5, 0.5 };
    fval = 0.0;
    for( int i=0; i < 5; i++ ) {
      double tmp=0;
      for( int j=0; j < npar; j++ )
	tmp +=pow( x[j]-a[i][j], 2.0 );
      fval -= 1.0 / ( tmp + c[i] );
    }
  }
  template<typename Real, typename Type>
  void Shekel7( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel7 func must be 4\n" );
    static const Real a[10][4] = { {4,4,4,4}, {1,1,1,1}, {8,8,8,8},
				   {6,6,6,6}, {3,7,3,7}, {2,9,2,9},
				   {5,5,3,3}, {8,1,8,1}, {6,2,6,2},
				   {7,3.6,7,3.6} };
    static const double c[10] = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3,
				  0.7, 0.5, 0.5 };
    fval = 0.0;
    for( int i=0; i < 7; i++ ) {
      double tmp=0;
      for( int j=0; j < npar; j++ )
	tmp +=pow( x[j]-a[i][j], 2.0 );
      fval -= 1.0 / ( tmp + c[i] );
    }
  }
  template<typename Real, typename Type>
  void Shekel10( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel10 func must be 4\n" );
    static const Real a[10][4] = { {4,4,4,4}, {1,1,1,1}, {8,8,8,8},
				   {6,6,6,6}, {3,7,3,7}, {2,9,2,9},
				   {5,5,3,3}, {8,1,8,1}, {6,2,6,2},
				   {7,3.6,7,3.6} };
    static const double c[10] = { 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3,
				  0.7, 0.5, 0.5 };
    fval = 0.0;
    for( int i=0; i < 10; i++ ) {
      double tmp=0;
      for( int j=0; j < npar; j++ )
	tmp +=pow( x[j]-a[i][j], 2.0 );
      fval -= 1.0 / ( tmp + c[i] );
    }
  }
  template<typename Real>
  void ShekelInit( int npar, int& mfct, Real& answer, Real* x,
		   Real* lo, Real* hi, int m=5) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel func must 4\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 0.0;
      x[ ii ] = 8.0;
      hi[ ii ] = 10.0;
    }

    switch( m ) {
    case 5:
      answer = - 10.1532;
      break;
    case 7:
      answer = - 10.4029;
      break;
    case 10:
      answer = - 10.5364;
      break;
    default:
      printf( "ShekelInit(): m=%d\t is not an option\n", m );
      return;
    }

  }
  template<typename Real>
  void Shekel5Init( int npar, int& mfct, Real& answer, Real* x,
		    Real* lo, Real* hi ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel func must be 4\n" );
    ShekelInit( npar, mfct, answer, x, lo, hi, 5 );
  }
  template<typename Real>
  void Shekel7Init( int npar, int& mfct, Real& answer, Real* x,
		    Real* lo, Real* hi ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel func must be 7\n" );
    ShekelInit( npar, mfct, answer, x, lo, hi, 7 );
  }
  template<typename Real>
  void Shekel10Init( int npar, int& mfct, Real& answer, Real* x,
		    Real* lo, Real* hi ) {
    if ( 4 != npar )
      throw std::runtime_error( "npar for the Shekel func must be 10\n" );
    ShekelInit( npar, mfct, answer, x, lo, hi, 10 );
  }

  template<typename Real, typename Type>
  void ShekelModified( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    static const double a[30][10] = {
      {9.681, 0.667, 4.783, 9.095, 3.517, 9.325, 6.544, 0.211, 5.122, 2.020},
      {9.400, 2.041, 3.788, 7.931, 2.882, 2.672, 3.568, 1.284, 7.033, 7.374},
      {8.025, 9.152, 5.114, 7.621, 4.564, 4.711, 2.996, 6.126, 0.734, 4.982},
      {2.196, 0.415, 5.649, 6.979, 9.510, 9.166, 6.304, 6.054, 9.377, 1.426},
      {8.074, 8.777, 3.467, 1.863, 6.708, 6.349, 4.534, 0.276, 7.633, 1.567},
      {7.650, 5.658, 0.720, 2.764, 3.278, 5.283, 7.474, 6.274, 1.409, 8.208},
      {1.256, 3.605, 8.623, 6.905, 0.584, 8.133, 6.071, 6.888, 4.187, 5.448},
      {8.314, 2.261, 4.224, 1.781, 4.124, 0.932, 8.129, 8.658, 1.208, 5.762},
      {0.226, 8.858, 1.420, 0.945, 1.622, 4.698, 6.228, 9.096, 0.972, 7.637},
      {7.305, 2.228, 1.242, 5.928, 9.133, 1.826, 4.060, 5.204, 8.713, 8.247},
      {0.652, 7.027, 0.508, 4.876, 8.807, 4.632, 5.808, 6.937, 3.291, 7.016},
      {2.699, 3.516, 5.874, 4.119, 4.461, 7.496, 8.817, 0.690, 6.593, 9.789},
      {8.327, 3.897, 2.017, 9.570, 9.825, 1.150, 1.395, 3.885, 6.354, 0.109},
      {2.132, 7.006, 7.136, 2.641, 1.882, 5.943, 7.273, 7.691, 2.880, 0.564},
      {4.707, 5.579, 4.080, 0.581, 9.698, 8.542, 8.077, 8.515, 9.231, 4.670},
      {8.304, 7.559, 8.567, 0.322, 7.128, 8.392, 1.472, 8.524, 2.277, 7.826},
      {8.632, 4.409, 4.832, 5.768, 7.050, 6.715, 1.711, 4.323, 4.405, 4.591},
      {4.887, 9.112, 0.170, 8.967, 9.693, 9.867, 7.508, 7.770, 8.382, 6.740},
      {2.440, 6.686, 4.299, 1.007, 7.008, 1.427, 9.398, 8.480, 9.950, 1.675},
      {6.306, 8.583, 6.084, 1.138, 4.350, 3.134, 7.853, 6.061, 7.457, 2.258},
      {0.652, 2.343, 1.370, 0.821, 1.310, 1.063, 0.689, 8.819, 8.833, 9.070},
      {5.558, 1.272, 5.756, 9.857, 2.279, 2.764, 1.284, 1.677, 1.244, 1.234},
      {3.352, 7.549, 9.817, 9.437, 8.687, 4.167, 2.570, 6.540, 0.228, 0.027},
      {8.798, 0.880, 2.370, 0.168, 1.701, 3.680, 1.231, 2.390, 2.499, 0.064},
      {1.460, 8.057, 1.336, 7.217, 7.914, 3.615, 9.981, 9.198, 5.292, 1.224},
      {0.432, 8.645, 8.774, 0.249, 8.081, 7.461, 4.416, 0.652, 4.002, 4.644},
      {0.679, 2.800, 5.523, 3.049, 2.968, 7.225, 6.730, 4.199, 9.614, 9.229},
      {4.263, 1.074, 7.286, 5.599, 8.291, 5.200, 9.214, 8.272, 4.398, 4.506},
      {9.496, 4.830, 3.150, 8.270, 5.079, 1.231, 5.731, 9.494, 1.883, 9.732},
      {4.138, 2.562, 2.532, 9.661, 5.611, 5.500, 6.886, 2.341, 9.699, 6.500}};
    static double c[] = { 0.806, 0.517, 0.1, 0.908, 0.965, 0.669, 0.524,
			  0.902, 0.531, 0.876, 0.462, 0.491, 0.463, 0.714,
			  0.352, 0.869, 0.813, 0.811, 0.828, 0.964, 0.789,
			  0.360, 0.369, 0.992, 0.332, 0.817, 0.632, 0.883,
			  0.608, 0.326};
    for( int ii = 0; ii < 30; ++ii ) {
      double sp = 0.0;
      for ( int jj = 0; jj < npar; ++jj ) {
	double tmp = x[ jj ] - a[ ii ][ jj ];
	sp += tmp * tmp;
      }
      fval += 1.0 / ( sp + c[ ii ] );
    }
    fval = - fval;
  }
  template<typename Real>
  void ShekelModifiedInit( int npar, int& mfct, Real& answer, Real* x,
			   Real* lo, Real* hi ) {
    if ( 2 == npar )
      // f( 8.0240653, 9.1465340 ) = -12.11900837975
      answer = -12.11900837975;
    else if ( 5 == npar )
      // f( 8.02491489, 9.15172576, 5.11392781, 7.62086096, 4.56408839 ) = -10.40561723899245
      answer = -10.40561723899245;
    else
      throw std::runtime_error( "npar for the Shekel func must either 2 or 5\n" );
    mfct = 0;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = 0.0;
      hi[ ii ] = 10.0;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lo[ ii ] + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
  }

  template<typename Real, typename Type>
  void Shubert( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Shubert func must be 2\n" );
    fval = -sin(2.0*x[0]+1.0)-2.0*sin(3.0*x[0]+2.0)-3.0*sin(4.0*x[0]+3.0)
      -4.0*sin(5.0*x[0]+4.0)-5.0*sin(6.0*x[0]+5.0)-sin(2.0*x[1]+1.0)
      -2.0*sin(3.0*x[1]+2.0)-3.0*sin(4.0*x[1]+3.0)-4.0*sin(5.0*x[1]+4.0)
      -5.0*sin(6.0*x[1]+5.0);
  }
  template<typename Real>
  void ShubertInit( int npar, int& mfct, Real& answer, Real* x,
		    Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Shubert func must be 2\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -10.0;
      hi[ ii ] = 10.0;
      x[ ii ] = 7.0;
    }
    answer = -24.06;
  }

  template<typename Real, typename Type>
  void SixHumpCamel( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the SixHumpCamel func must be 2\n" );
    Real x1=x[ 0 ], x2 = x[ 1 ];
    fval = 4.0*x1*x1-0.21E1*pow(x1,4.0)+pow(x1,6.0)/3+x1*x2-4.0*x2*x2
      + 4.0*pow(x2,4.0);
  }
  template<typename Real>
  void SixHumpCamelInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the SixHumpCamel func must be 2\n" );
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -5.0;
      hi[ ii ] = 5.0;
      x[ ii ] = 4.0;
    }
    answer = -1.03;
  }

  template<typename Real, typename Type>
  void Trecanni( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Trecanni func must be 2\n" );
    fval = ( ( x[ 0 ]  +  4.0 ) * x[ 0 ] + 4.0 )  * mysqr( x[ 0 ]  ) +
      mysqr( x[ 1 ] );
  }
  template<typename Real>
  void TrecanniInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Trecanni func must be 2\n" );
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      lo[ ii ] = -5.0;
      hi[ ii ] = 5.0;
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lo[ ii ] + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = 0.0;
  }
  //
  // f( 0, 0 ) = f( -2, 0 ) = 0.0
  //

  template<typename Real, typename Type>
  void Trefethen4( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Trefethen4 func must be 2\n" );
    fval = ( ( x[ 0 ]  +  4.0 ) * x[ 0 ] + 4.0 )  * mysqr( x[ 0 ]  ) +
      mysqr( x[ 1 ] );
  }
  template<typename Real>
  void Trefethen4Init( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {
    if ( 2 != npar )
      throw std::runtime_error( "npar for the Trefethen4 func must be 2\n" );
    lo[ 0 ] = -6.5;
    lo[ 1 ] = -4.5;
    hi[ 0 ] = 6.5;
    hi[ 1 ] = 4.5;
    srand(1);
    for ( int ii = 0; ii < npar; ++ii ) {
      double tmp = rand( ) / (RAND_MAX + 1.0);
      x[ ii ] = lo[ ii ] + ( hi[ ii ] - lo[ ii ] ) * tmp;
    }
    answer = -3.30686865;
  }
  //
  // f( -0.0244031, 0.2106124 ) = -3.30686865
  //

  template<typename Real, typename Type>
  void Trigonometric( int mfct, int npar, Real* x, Real* fvec, int& ierr,
		      Type xptr ) {

    Real sum = cos( x[ 0 ] );
    for ( int ii = 1; ii < npar; ++ii )
      sum += cos( x[ ii ] );

    for ( int ii = 0; ii < npar; ++ii )
      fvec[ ii ] = npar - sum -
	npar * ii * ( 1.0 - cos( x[ ii ] ) ) - sin( x[ ii ] );

  }
  template<typename Real, typename Type>
  void Trigonometric( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = npar;
    std::vector< Real> fvec( mfct );
    Trigonometric( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void TrigonometricInit( int npar, int& mfct, Real& answer, Real* x,
			  Real* lo, Real* hi ) {

    mfct = npar;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0 / npar;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // This function is known to have the global minimum
  // f( 0.0429645656464827, 0.043976286943042, 0.0450933996565844;
  //    0.0463389157816542, 0.0477443839560646, 0.0493547352444078;
  //    0.0512373449259557, 0.195209463715277, 0.164977664720328;
  //    0.0601485854398078 ) = 0.0
  // and a local minimum
  // f(  0.0551509, 0.0568408, 0.0587639, 0.0609906, 0.0636263, 0.0668432;
  //     0.208162, 0.164363, 0.0850068, 0.0914314 ) = 2.79506e-5
  //

  template<typename Real, typename Type>
  void VariablyDimensioned( int mfct, int npar, Real* x, Real* fvec,
			    int& ierr, Type xptr ) {

    Real sum = 0.0;
    for ( int ii = 1; ii <= npar; ++ii ) {
      Real xim1 = x[ii-1] - 1.0;
      sum  += ii*xim1;
      fvec[ ii - 1 ] = xim1;
    }
    fvec[ npar ] = sum;
    fvec[ npar+1 ] = sum*sum;

  }
  template<typename Real, typename Type>
  void VariablyDimensioned( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 2 + npar;
    std::vector< Real> fvec( mfct );
    VariablyDimensioned( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void VariablyDimensionedInit( int npar, int& mfct, Real& answer,
				Real* x, Real* lo, Real* hi ) {

    mfct = npar + 2;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 1.0 - ( ii + 1 ) / npar;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // f( 1, ... , 1 ) = 0.0
  //

  template<typename Real, typename Type>
  void Watson( int mfct, int npar, Real* x, Real* fvec, int& ierr,
	       Type xptr ) {

    for ( int ii = 1; ii <= 29; ++ii ) {
      Real div = ii / 29.0;
      Real s1 = 0.0;
      Real dx = 1.0;
      for ( int  jj = 2; jj <= npar; ++jj ) {
	s1 += ( jj - 1 ) * dx * x[ jj - 1 ];
	dx *= div;
      }
      Real s2 = 0.0;
      dx = 1.0;
      for ( int jj = 1; jj <= npar; ++jj ) {
	s2 += dx * x[ jj - 1 ];
	dx *= div;
      }
      fvec[ ii - 1 ]  = s1 - s2 * s2 - 1.0;
    }
    fvec[ 29 ] = x[ 0 ];
    fvec[ 30 ] = x[ 1 ] - x[ 0 ] * x[ 0 ] - 1.0;

  }
  template<typename Real, typename Type>
  void Watson( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 31;
    std::vector< Real> fvec( mfct );
    Watson( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void WatsonInit( int npar, int& mfct, Real& answer, Real* x,
		   Real* lo, Real* hi ) {

    if ( 6 != npar )
      throw std::runtime_error( "npar for the Watson func must be 6\n" );

    mfct = 31;
    for ( int ii = 0; ii < npar; ++ii )
      x[ ii ] = 0.0;

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 2.28767e-3;
  }
  //
  // f=2.28767...10^(-3)    if n=6
  // f=1.39976...10^(-6)    if n=9
  // f=4.72238...10^(-10)   if n=12
  // f=2.48631...10^(-20)   if n=20
  //

  template<typename Real, typename Type>
  void Wood( int mfct, int npar, Real* x, Real* fvec, int& ierr, Type xptr ) {

    Real sqrt10 = sqrt( 10.0 );
    Real sqrt90 = sqrt( 90.0 );
    Real invsqrt10 = 1.0 / sqrt10;

    fvec[ 0 ] = 10.0 * ( x[ 1 ] - x[ 0 ] * x[ 0 ] );
    fvec[ 1 ] = 1.0 - x[ 0 ];
    fvec[ 2 ] = sqrt90 * ( x[ 3 ] - x[ 2 ] * x[ 2 ] );
    fvec[ 3 ] = 1.0 - x[ 2 ];
    fvec[ 4 ] = sqrt10 * ( x[ 1 ] + x[ 3 ] - 2.0 );
    fvec[ 5 ] = invsqrt10 * ( x[ 1 ] - x[ 3 ] );

  }
  template<typename Real, typename Type>
  void Wood( int npar, Real* x, Real& fval, int& ierr, Type xptr ) {

    int mfct = 6;
    std::vector< Real> fvec( mfct );
    Wood( mfct, npar, x, &fvec[0], ierr, 0 );
    fval = 0.0;
    for ( int ii = mfct - 1; ii >= 0; --ii )
      fval += fvec[ii]*fvec[ii];

  }
  template<typename Real>
  void WoodInit( int npar, int& mfct, Real& answer, Real* x, Real* lo, Real* hi ) {

    if ( npar % 4 )
      throw std::runtime_error( "npar for the Wood func must be multiple of4\n" );

    mfct = 6;
    for ( int ii = 0; ii < npar; ii += 4 ) {
      x[ ii ] = -3.0;
      x[ ii + 1 ] = -1.0;
      x[ ii + 2 ] = -3.0;
      x[ ii + 3 ] = -1.0;
    }

    for ( int ii = 0; ii < npar; ++ii )
      lo[ ii ] = -1.0e6;

    for ( int ii = 0; ii < npar; ++ii )
      hi[ ii ] = 1.0e6;

    answer = 0.0;

  }
  //
  // f( 1, 1, 1, 1 ) = 0.0
  //

} // namespace

#endif // #ifndef uncoptfunc_hh
