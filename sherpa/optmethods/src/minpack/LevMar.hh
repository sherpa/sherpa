#ifndef LevMar_hh
#define LevMar_hh

//
//
// Minpack Copyright Notice (1999) University of Chicago.  All rights reserved
// 
// Redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the
// following conditions are met:
// 
// 1. Redistributions of source code must retain the above
// copyright notice, this list of conditions and the following
// disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following
// disclaimer in the documentation and/or other materials
// provided with the distribution.
// 
// 3. The end-user documentation included with the
// redistribution, if any, must include the following
// acknowledgment:
// 
//    "This product includes software developed by the
//    University of Chicago, as Operator of Argonne National
//    Laboratory.
// 
// Alternately, this acknowledgment may appear in the software
// itself, if and wherever such third-party acknowledgments
// normally appear.
// 
// 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
// WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
// UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
// THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
// OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
// OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
// USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
// THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
// DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
// UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
// BE CORRECTED.
// 
// 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
// HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
// ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
// INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
// ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
// PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
// SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
// (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
// EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
// POSSIBILITY OF SUCH LOSS OR DAMAGES.
//
//

/* lmdif.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <cmath>
#include "../Opt.hh"
namespace minpack {

  /*
  // from size_t to int
  template < typename Type0, typename Type1 >
  Type1 convertme( Type0 arg0 ) {

    if ( arg0 < std::numeric_limits< Type1 >::max() )
	 return static_cast< Type1 >( arg0 );
  }
  */

  template < typename Func, typename Data >
  class LevMar : public sherpa::Opt {
    
  public:

    LevMar( Func func, Data xdata, int mfct )
      : sherpa::Opt( ), usr_func( func ), usr_data( xdata ), myfvec( mfct ) { }

    int operator( )( int n, double ftol, double xtol,
		     double gtol, int maxfev, double epsfcn,
		     double factor, int nprint, 
		     const std::vector<double>& low,
		     const std::vector<double>& high, std::vector<double>& x,
		     int& nfev, double& fmin, std::vector<double>& covarerr ) {

      int info = 0;

      try {

	const sherpa::Opt::mypair limits( low, high );
	if ( sherpa::Opt::are_pars_outside_limits( n, limits, x ) )
	  throw sherpa::OptErr( sherpa::OptErr::OutOfBound );

	int m = static_cast<int>( myfvec.size( ) );

	std::vector<double> diag( n ), qtf( n ), wa1( n ), wa2( n ), wa3( n );
	std::vector<double> wa4( m ), fjac( m * n);
	std::vector<int> ipvt( n );

	const int mode = 1;
	const int ldfjac = m;

	info = lmdif( usr_func, usr_data, m, n, &x[0], &myfvec[0], ftol, xtol,
		      gtol, maxfev, epsfcn, &diag[0], mode, factor, nprint,
		      nfev, &fjac[0], ldfjac, &ipvt[0], &qtf[0], &wa1[ 0 ],
		      &wa2[0], &wa3[0], &wa4[0], low, high );

	covar( n, &fjac[ 0 ], ldfjac, &ipvt[0], ftol, &wa1[0] );

	for ( int ii = 0; ii < n; ++ii )
	  if ( fjac[ ii + ldfjac * ii ] > 0.0 )
	    covarerr[ ii ] = sqrt( fjac[ ii + ldfjac * ii ] );
	  else
	    covarerr[ ii ] = 0.0;

      } catch( sherpa::OptErr& oe ) {

	if ( nprint )
	  std::cerr << oe << '\n';
	info = oe.err;

      } catch( std::runtime_error& re ) {

	if ( nprint )
	  std::cerr << re.what( ) << '\n';
	info = sherpa::OptErr::Unknown;

      } catch( std::exception& e ) {

	if ( nprint )
	  std::cerr << e.what( ) << '\n';
	info = sherpa::OptErr::Unknown;

      }

      fmin = std::pow( enorm( myfvec.size( ), &myfvec[0] ), 2.0 );
      return info;

    }

    double eval_func( int maxnfev, const sherpa::Opt::mypair& limits,
		      int npar, sherpa::Opt::myvec& par, int& nfev ) {


      double fval = std::numeric_limits< double >::max( );
      if ( sherpa::Opt::are_pars_outside_limits( npar, limits, par ) )
	return fval;

      ++nfev;
      int ierr=EXIT_SUCCESS;

      int mymfct = static_cast<int>( myfvec.size( ) );

      usr_func( mymfct, npar, &par[0], &myfvec[0], ierr, usr_data );

      fval = pow( enorm( myfvec.size( ), &myfvec[0] ), 2.0 );
      if ( EXIT_SUCCESS != ierr )
	throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
      if ( nfev >= maxnfev )
	throw sherpa::OptErr( sherpa::OptErr::MaxFev );

      return fval;

    }


    // de
    int minimize( int maxnfev, const sherpa::Opt::mypair& limits,
		  double tol, int npar, sherpa::Opt::myvec& par, double& fmin,
		  int& nfev ) {
      int nprint = 0;
      double epsfcn =
	std::sqrt( std::numeric_limits< double >::epsilon( ) );
      double factor = 100.0;
      
      const Opt::myvec& low = limits.first;
      const Opt::myvec& high = limits.second;
      std::vector<double> diag( npar );
      std::vector<double> covarerr( npar );
      return this->operator( )( npar, tol, tol, tol, maxnfev, epsfcn, factor,
				nprint, low, high, par, nfev, fmin, covarerr );
    }
    // de

    // c     **********
    // c
    // c     subroutine covar
    // c
    // c     given an m by n matrix a, the problem is to determine
    // c     the covariance matrix corresponding to a, defined as
    // c
    // c                    t
    // c           inverse(a *a) .
    // c
    // c     this subroutine completes the solution of the problem
    // c     if it is provided with the necessary information from the
    // c     qr factorization, with column pivoting, of a. that is, if
    // c     a*p = q*r, where p is a permutation matrix, q has orthogonal
    // c     columns, and r is an upper triangular matrix with diagonal
    // c     elements of nonincreasing magnitude, then covar expects
    // c     the full upper triangle of r and the permutation matrix p.
    // c     the covariance matrix is then computed as
    // c
    // c                      t     t
    // c           p*inverse(r *r)*p  .
    // c
    // c     if a is nearly rank deficient, it may be desirable to compute
    // c     the covariance matrix corresponding to the linearly independent
    // c     columns of a. to define the numerical rank of a, covar uses
    // c     the tolerance tol. if l is the largest integer such that
    // c
    // c           abs(r(l,l)) .gt. tol*abs(r(1,1)) ,
    // c
    // c     then covar computes the covariance matrix corresponding to
    // c     the first l columns of r. for k greater than l, column
    // c     and row ipvt(k) of the covariance matrix are set to zero.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine covar(n,r,ldr,ipvt,tol,wa)
    // c
    // c     where
    // c
    // c       n is a positive integer input variable set to the order of r.
    // c
    // c       r is an n by n array. on input the full upper triangle must
    // c         contain the full upper triangle of the matrix r. on output
    // c         r contains the square symmetric covariance matrix.
    // c
    // c       ldr is a positive integer input variable not less than n
    // c         which specifies the leading dimension of the array r.
    // c
    // c       ipvt is an integer input array of length n which defines the
    // c         permutation matrix p such that a*p = q*r. column j of p
    // c         is column ipvt(j) of the identity matrix.
    // c
    // c       tol is a nonnegative input variable used to define the
    // c         numerical rank of a in the manner described above.
    // c
    // c       wa is a work array of length n.
    // c
    // c     subprograms called
    // c
    // c       fortran-supplied ... dabs
    // c
    // c     argonne national laboratory. minpack project. august 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    void covar( int n, double *r__, int ldr, const int *ipvt, double tol,
		double *wa ) {
      // System generated locals
      int r_dim1, r_offset, i__1, i__2, i__3;

      // Local variables
      int i__, j, k, l, ii, jj, km1;
      int sing;
      double temp, tolr;
      // Parameter adjustments
      --wa;
      --ipvt;
      tolr = tol * fabs(r__[0]);
      r_dim1 = ldr;
      r_offset = 1 + r_dim1;
      r__ -= r_offset;

      // Function Body

      //     form the inverse of r in the full upper triangle of r.

      l = 0;
      i__1 = n;
      for (k = 1; k <= i__1; ++k) {
	if (fabs(r__[k + k * r_dim1]) <= tolr) {
	  goto L50;
	}
	r__[k + k * r_dim1] = 1. / r__[k + k * r_dim1];
	km1 = k - 1;
	if (km1 < 1) {
	  goto L30;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	  temp = r__[k + k * r_dim1] * r__[j + k * r_dim1];
	  r__[j + k * r_dim1] = 0.;
	  i__3 = j;
	  for (i__ = 1; i__ <= i__3; ++i__) {
	    r__[i__ + k * r_dim1] -= temp * r__[i__ + j * r_dim1];
	    // L10:
	  }
	  // L20:
	}
      L30:
	l = k;
	// L40:
      }
    L50:

      //     form the full upper triangle of the inverse of (r transpose)*r
      //     in the full upper triangle of r.

      if (l < 1) {
	goto L110;
      }
      i__1 = l;
      for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	if (km1 < 1) {
	  goto L80;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	  temp = r__[j + k * r_dim1];
	  i__3 = j;
	  for (i__ = 1; i__ <= i__3; ++i__) {
	    r__[i__ + j * r_dim1] += temp * r__[i__ + k * r_dim1];
	    // L60:
	  }
	  // L70:
	}
      L80:
	temp = r__[k + k * r_dim1];
	i__2 = k;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  r__[i__ + k * r_dim1] = temp * r__[i__ + k * r_dim1];
	  // L90:
	}
	// L100:
      }
    L110:

      //     form the full lower triangle of the covariance matrix
      //     in the strict lower triangle of r and in wa.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	jj = ipvt[j];
	sing = j > l;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  if (sing) {
	    r__[i__ + j * r_dim1] = 0.;
	  }
	  ii = ipvt[i__];
	  if (ii > jj) {
	    r__[ii + jj * r_dim1] = r__[i__ + j * r_dim1];
	  }
	  if (ii < jj) {
	    r__[jj + ii * r_dim1] = r__[i__ + j * r_dim1];
	  }
	  // L120:
	}
	wa[jj] = r__[j + j * r_dim1];
	// L130:
      }

      //     symmetrize the covariance matrix in r.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
	  // L140:
	}
	r__[j + j * r_dim1] = wa[j];
	// L150:
      }
      //return 0

      //     last card of subroutine covar.

    } // covar

    //
    // c     **********
    // c
    // c     function enorm
    // c
    // c     given an n-vector x, this function calculates the
    // c     euclidean norm of x.
    // c
    // c     the euclidean norm is computed by accumulating the sum of
    // c     squares in three different sums. the sums of squares for the
    // c     small and large components are scaled so that no overflows
    // c     occur. non-destructive underflows are permitted. underflows
    // c     and overflows do not occur in the computation of the unscaled
    // c     sum of squares for the intermediate components.
    // c     the definitions of small, intermediate and large components
    // c     depend on two constants, rdwarf and rgiant. the main
    // c     restrictions on these constants are that rdwarf**2 not
    // c     underflow and rgiant**2 not overflow. the constants
    // c     given here are suitable for every known computer.
    // c
    // c     the function statement is
    // c
    // c       double precision function enorm(n,x)
    // c
    // c     where
    // c
    // c       n is a positive integer input variable.
    // c
    // c       x is an input array of length n.
    // c
    // c     subprograms called
    // c
    // c       fortran-supplied ... dabs,dsqrt
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    //

    template < typename Type >
    double enorm( Type n, const double *x) {
      
      // Initialized data
      const double rdwarf=3.834e-20;
      const double rgiant= 1.304e19;

      // System generated locals
      Type i__1;
      double ret_val=0.0, d__1;

      // Local variables
      Type i__;
      double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

      // Parameter adjustments
      --x;

      // Function Body
      s1 = 0.;
      s2 = 0.;
      s3 = 0.;
      x1max = 0.;
      x3max = 0.;
      floatn = (double) (n);
      agiant = rgiant / floatn;
      i__1 = n;
      for (i__ = 1; i__ <= i__1; ++i__) {
	xabs = fabs(x[i__]);
	if (xabs > rdwarf && xabs < agiant) {
	  goto L70;
	}
	if (xabs <= rdwarf) {
	  goto L30;
	}

	//              sum for large components.

	if (xabs <= x1max) {
	  goto L10;
	}
	// Computing 2nd power
	d__1 = x1max / xabs;
	s1 = 1. + s1 * (d__1 * d__1);
	x1max = xabs;
	goto L20;
      L10:
	// Computing 2nd power
	d__1 = xabs / x1max;
	s1 += d__1 * d__1;
      L20:
	goto L60;
      L30:

	//              sum for small components.

	if (xabs <= x3max) {
	  goto L40;
	}
	// Computing 2nd power
	d__1 = x3max / xabs;
	s3 = 1. + s3 * (d__1 * d__1);
	x3max = xabs;
	goto L50;
      L40:
	if (xabs != 0.) {
	  // Computing 2nd power
	  d__1 = xabs / x3max;
	  s3 += d__1 * d__1;
	}
      L50:
      L60:
	goto L80;
      L70:

	//           sum for intermediate components.

	// Computing 2nd power
	d__1 = xabs;
	s2 += d__1 * d__1;
      L80:
	// L90:
	;
      }

      //     calculation of norm.

      if (s1 == 0.) {
	goto L100;
      }
      ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
      goto L130;
    L100:
      if (s2 == 0.) {
	goto L110;
      }
      if (s2 >= x3max) {
	ret_val = sqrt(s2 * (1. + x3max / s2 * (x3max * s3)));
      }
      if (s2 < x3max) {
	ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
      }
      goto L120;
    L110:
      ret_val = x3max * sqrt(s3);
    L120:
    L130:
      return ret_val;

      //     last card of function enorm.

    } // enorm

  private:

    Func usr_func;
    Data usr_data;
    std::vector< double > myfvec;

    //
    // c     **********
    // c
    // c     subroutine fdjac2
    // c
    // c     this subroutine computes a forward-difference approximation
    // c     to the m by n jacobian matrix associated with a specified
    // c     problem of m functions in n variables.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
    // c
    // c     where
    // c
    // c       fcn is the name of the user-supplied subroutine which
    // c         calculates the functions. fcn must be declared
    // c         in an external statement in the user calling
    // c         program, and should be written as follows.
    // c
    // c         subroutine fcn(m,n,x,fvec,iflag)
    // c         integer m,n,iflag
    // c         double precision x(n),fvec(m)
    // c         ----------
    // c         calculate the functions at x and
    // c         return this vector in fvec.
    // c         ----------
    // c         return
    // c         end
    // c
    // c         the value of iflag should not be changed by fcn unless
    // c         the user wants to terminate execution of fdjac2.
    // c         in this case set iflag to a negative integer.
    // c
    // c       m is a positive integer input variable set to the number
    // c         of functions.
    // c
    // c       n is a positive integer input variable set to the number
    // c         of variables. n must not exceed m.
    // c
    // c       x is an input array of length n.
    // c
    // c       fvec is an input array of length m which must contain the
    // c         functions evaluated at x.
    // c
    // c       fjac is an output m by n array which contains the
    // c         approximation to the jacobian matrix evaluated at x.
    // c
    // c       ldfjac is a positive integer input variable not less than m
    // c         which specifies the leading dimension of the array fjac.
    // c
    // c       iflag is an integer variable which can be used to terminate
    // c         the execution of fdjac2. see description of fcn.
    // c
    // c       epsfcn is an input variable used in determining a suitable
    // c         step length for the forward-difference approximation. this
    // c         approximation assumes that the relative errors in the
    // c         functions are of the order of epsfcn. if epsfcn is less
    // c         than the machine precision, it is assumed that the relative
    // c         errors in the functions are of the order of the machine
    // c         precision.
    // c
    // c       wa is a work array of length m.
    // c
    // c     subprograms called
    // c
    // c       user-supplied ...... fcn
    // c
    // c       minpack-supplied ... dpmpar
    // c
    // c       fortran-supplied ... dabs,dmax1,dsqrt
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    // 
    int fdjac2( Func fcn, int m, int n, double *x, 
		const double *fvec, double *fjac, int ldfjac,
		double epsfcn, double *wa, Data xptr,
		const std::vector<double>& high ) {

    // System generated locals
    int fjac_dim1, fjac_offset, i__1, i__2;

    // Local variables
    double h__;
    int i__, j;
    double eps, temp, epsmch;
    // dtn
    int iflag=0;
    // dtn

    // Parameter adjustments
    --wa;
    --fvec;
    --x;
    fjac_dim1 = ldfjac;
    fjac_offset = 1 + fjac_dim1 * 1;
    fjac -= fjac_offset;

    // Function Body

    //     epsmch is the machine precision.

    epsmch = std::numeric_limits< double >::epsilon( );

    eps = sqrt((std::max(epsfcn,epsmch)));
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
      temp = x[j];
      h__ = eps * fabs(temp);
      if (h__ == 0.) {
	h__ = eps;
      }
      // dtn
      // if the parameter is beyond the upper boundary then
      // perform backwards-difference approximation
      if ( x[ j ] + h__ > high[ j - 1 ] )
	h__ = - h__;
      // dtn
      x[j] = temp + h__;
      fcn(m, n, &x[1], &wa[1], iflag, xptr);
      if (iflag < 0) {
	// goto L30;
	return iflag;
      }
      x[j] = temp;
      i__2 = m;
      for (i__ = 1; i__ <= i__2; ++i__) {
	fjac[i__ + j * fjac_dim1] = (wa[i__] - fvec[i__]) / h__;
	// L10:
      }
      // L20:
    }
    // L30:
    return iflag;

    //     last card of subroutine fdjac2.

    } // fdjac2

    //
    // c     **********
    // c
    // c     subroutine lmdif
    // c
    // c     the purpose of lmdif is to minimize the sum of the squares of
    // c     m nonlinear functions in n variables by a modification of
    // c     the levenberg-marquardt algorithm. the user must provide a
    // c     subroutine which calculates the functions. the jacobian is
    // c     then calculated by a forward-difference approximation.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,
    // c                        diag,mode,factor,nprint,info,nfev,fjac,
    // c                        ldfjac,ipvt,qtf,wa1,wa2,wa3,wa4)
    // c
    // c     where
    // c
    // c       fcn is the name of the user-supplied subroutine which
    // c         calculates the functions. fcn must be declared
    // c         in an external statement in the user calling
    // c         program, and should be written as follows.
    // c
    // c         subroutine fcn(m,n,x,fvec,iflag)
    // c         integer m,n,iflag
    // c         double precision x(n),fvec(m)
    // c         ----------
    // c         calculate the functions at x and
    // c         return this vector in fvec.
    // c         ----------
    // c         return
    // c         end
    // c
    // c         the value of iflag should not be changed by fcn unless
    // c         the user wants to terminate execution of lmdif.
    // c         in this case set iflag to a negative integer.
    // c
    // c       m is a positive integer input variable set to the number
    // c         of functions.
    // c
    // c       n is a positive integer input variable set to the number
    // c         of variables. n must not exceed m.
    // c
    // c       x is an array of length n. on input x must contain
    // c         an initial estimate of the solution vector. on output x
    // c         contains the final estimate of the solution vector.
    // c
    // c       fvec is an output array of length m which contains
    // c         the functions evaluated at the output x.
    // c
    // c       ftol is a nonnegative input variable. termination
    // c         occurs when both the actual and predicted relative
    // c         reductions in the sum of squares are at most ftol.
    // c         therefore, ftol measures the relative error desired
    // c         in the sum of squares.
    // c
    // c       xtol is a nonnegative input variable. termination
    // c         occurs when the relative error between two consecutive
    // c         iterates is at most xtol. therefore, xtol measures the
    // c         relative error desired in the approximate solution.
    // c
    // c       gtol is a nonnegative input variable. termination
    // c         occurs when the cosine of the angle between fvec and
    // c         any column of the jacobian is at most gtol in absolute
    // c         value. therefore, gtol measures the orthogonality
    // c         desired between the function vector and the columns
    // c         of the jacobian.
    // c
    // c       maxfev is a positive integer input variable. termination
    // c         occurs when the number of calls to fcn is at least
    // c         maxfev by the end of an iteration.
    // c
    // c       epsfcn is an input variable used in determining a suitable
    // c         step length for the forward-difference approximation. this
    // c         approximation assumes that the relative errors in the
    // c         functions are of the order of epsfcn. if epsfcn is less
    // c         than the machine precision, it is assumed that the relative
    // c         errors in the functions are of the order of the machine
    // c         precision.
    // c
    // c       diag is an array of length n. if mode = 1 (see
    // c         below), diag is internally set. if mode = 2, diag
    // c         must contain positive entries that serve as
    // c         multiplicative scale factors for the variables.
    // c
    // c       mode is an integer input variable. if mode = 1, the
    // c         variables will be scaled internally. if mode = 2,
    // c         the scaling is specified by the input diag. other
    // c         values of mode are equivalent to mode = 1.
    // c
    // c       factor is a positive input variable used in determining the
    // c         initial step bound. this bound is set to the product of
    // c         factor and the euclidean norm of diag*x if nonzero, or else
    // c         to factor itself. in most cases factor should lie in the
    // c         interval (.1,100.). 100. is a generally recommended value.
    // c
    // c       nprint is an integer input variable that enables controlled
    // c         printing of iterates if it is positive. in this case,
    // c         fcn is called with iflag = 0 at the beginning of the first
    // c         iteration and every nprint iterations thereafter and
    // c         immediately prior to return, with x and fvec available
    // c         for printing. if nprint is not positive, no special calls
    // c         of fcn with iflag = 0 are made.
    // c
    // c       info is an integer output variable. if the user has
    // c         terminated execution, info is set to the (negative)
    // c         value of iflag. see description of fcn. otherwise,
    // c         info is set as follows.
    // c
    // c         info = 0  improper input parameters.
    // c
    // c         info = 1  both actual and predicted relative reductions
    // c                   in the sum of squares are at most ftol.
    // c
    // c         info = 2  relative error between two consecutive iterates
    // c                   is at most xtol.
    // c
    // c         info = 3  conditions for info = 1 and info = 2 both hold.
    // c
    // c         info = 4  the cosine of the angle between fvec and any
    // c                   column of the jacobian is at most gtol in
    // c                   absolute value.
    // c
    // c         info = 5  number of calls to fcn has reached or
    // c                   exceeded maxfev.
    // c
    // c         info = 6  ftol is too small. no further reduction in
    // c                   the sum of squares is possible.
    // c
    // c         info = 7  xtol is too small. no further improvement in
    // c                   the approximate solution x is possible.
    // c
    // c         info = 8  gtol is too small. fvec is orthogonal to the
    // c                   columns of the jacobian to machine precision.
    // c
    // c       nfev is an integer output variable set to the number of
    // c         calls to fcn.
    // c
    // c       fjac is an output m by n array. the upper n by n submatrix
    // c         of fjac contains an upper triangular matrix r with
    // c         diagonal elements of nonincreasing magnitude such that
    // c
    // c                t     t           t
    // c               p *(jac *jac)*p = r *r,
    // c
    // c         where p is a permutation matrix and jac is the final
    // c         calculated jacobian. column j of p is column ipvt(j)
    // c         (see below) of the identity matrix. the lower trapezoidal
    // c         part of fjac contains information generated during
    // c         the computation of r.
    // c
    // c       ldfjac is a positive integer input variable not less than m
    // c         which specifies the leading dimension of the array fjac.
    // c
    // c       ipvt is an integer output array of length n. ipvt
    // c         defines a permutation matrix p such that jac*p = q*r,
    // c         where jac is the final calculated jacobian, q is
    // c         orthogonal (not stored), and r is upper triangular
    // c         with diagonal elements of nonincreasing magnitude.
    // c         column j of p is column ipvt(j) of the identity matrix.
    // c
    // c       qtf is an output array of length n which contains
    // c         the first n elements of the vector (q transpose)*fvec.
    // c
    // c       wa1, wa2, and wa3 are work arrays of length n.
    // c
    // c       wa4 is a work array of length m.
    // c
    // c     subprograms called
    // c
    // c       user-supplied ...... fcn
    // c
    // c       minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac
    // c
    // c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    // 
    int lmdif( Func fcn, Data xptr, int m, int n, double *x, 
	       double *fvec, double ftol, double xtol, double gtol,
	       int maxfev, double epsfcn, double *diag, int mode,
	       double factor, int nprint, int& nfev, double *fjac,
	       int ldfjac, int *ipvt, double * qtf, double *wa1,
	       double *wa2, double *wa3, double *wa4,
	       const std::vector<double>& low,
	       const std::vector<double>& high ) {

      // Initialized data

      const double p1=0.1;
      const double p5=0.5;
      const double p25=0.25;
      const double p75=0.75;
      const double p0001=1.0e-4;

      // System generated locals
      int fjac_dim1, fjac_offset, i__1, i__2;
      double d__1, d__2, d__3;

      // Local variables
      int i__, j, l;
      double par, sum;
      int iter;
      double temp=0.0, temp1, temp2;
      int iflag=0;
      double delta=0.0;
      double ratio;
      double fnorm, gnorm;
      double pnorm, xnorm=0.0, fnorm1, actred, dirder, epsmch, prered;
      int info=0;

      // Parameter adjustments
      --wa4;
      --fvec;
      --wa3;
      --wa2;
      --wa1;
      --qtf;
      --ipvt;
      --diag;
      --x;
      fjac_dim1 = ldfjac;
      fjac_offset = 1 + fjac_dim1 * 1;
      fjac -= fjac_offset;

      // Function Body

      //     epsmch is the machine precision.

      epsmch = std::numeric_limits< double >::epsilon( );

      nfev = 0;

      //     check the input parameters for errors.

      if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. || 
	  gtol < 0. || maxfev <= 0 || factor <= 0.) {
	goto L300;
      }
      if (mode != 2) {
	goto L20;
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	if (diag[j] <= 0.) {
	  goto L300;
	}
	// L10:
      }
    L20:

      //     evaluate the function at the starting point
      //     and calculate its norm.

      fcn(m, n, &x[1], &fvec[1], iflag, xptr);
      nfev = 1;
      if (iflag < 0) {
	goto L300;
      }
      fnorm = enorm(m, &fvec[1]);

      //     initialize levenberg-marquardt parameter and iteration counter.

      par = 0.;
      iter = 1;

      //     beginning of the outer loop.

    L30:

      //        calculate the jacobian matrix.

      iflag = fdjac2(fcn, m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac,
		     epsfcn, &wa4[1], xptr, high );
      nfev += n;
      if (iflag < 0) {
	goto L300;
      }

      //        if requested, call fcn to enable printing of iterates.

      if (nprint <= 0) {
	goto L40;
      }
      iflag = 0;
      if ((iter - 1) % nprint == 0) {
	fcn(m, n, &x[1], &fvec[1], iflag, xptr);
      }
      if (iflag < 0) {
	goto L300;
      }
    L40:

      //        compute the qr factorization of the jacobian.

      qrfac(m, n, &fjac[fjac_offset], ldfjac, 1, &ipvt[1], n, &wa1[1], &
	    wa2[1], &wa3[1]);

      //        on the first iteration and if mode is 1, scale according
      //        to the norms of the columns of the initial jacobian.

      if (iter != 1) {
	goto L80;
      }
      if (mode == 2) {
	goto L60;
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	diag[j] = wa2[j];
	if (wa2[j] == 0.) {
	  diag[j] = 1.;
	}
	// L50:
      }
    L60:

      //        on the first iteration, calculate the norm of the scaled x
      //        and initialize the step bound delta.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa3[j] = diag[j] * x[j];
	// L70:
      }
      xnorm = enorm(n, &wa3[1]);
      delta = factor * xnorm;
      if (delta == 0.) {
	delta = factor;
      }
    L80:

      //        form (q transpose)*fvec and store the first n components in
      //        qtf.

      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
	wa4[i__] = fvec[i__];
	// L90:
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	if (fjac[j + j * fjac_dim1] == 0.) {
	  goto L120;
	}
	sum = 0.;
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	  sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
	  // L100:
	}
	temp = -sum / fjac[j + j * fjac_dim1];
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	  wa4[i__] += fjac[i__ + j * fjac_dim1] * temp;
	  // L110:
	}
      L120:
	fjac[j + j * fjac_dim1] = wa1[j];
	qtf[j] = wa4[j];
	// L130:
      }

      //        compute the norm of the scaled gradient.

      gnorm = 0.;
      if (fnorm == 0.) {
	goto L170;
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	if (wa2[l] == 0.) {
	  goto L150;
	}
	sum = 0.;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  sum += fjac[i__ + j * fjac_dim1] * (qtf[i__] / fnorm);
	  // L140:
	}
	// Computing MAX
	d__2 = gnorm, d__3 = fabs(sum / wa2[l]);
	gnorm = std::max(d__2,d__3);
      L150:
	// L160:
	;
      }
    L170:

      //        test for convergence of the gradient norm.

      if (gnorm <= gtol) {
	info = 4;
      }
      if (info != 0) {
	goto L300;
      }

      //        rescale if necessary.

      if (mode == 2) {
	goto L190;
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	// Computing MAX
	d__1 = diag[j], d__2 = wa2[j];
	diag[j] = std::max(d__1,d__2);
	// L180:
      }
    L190:

      //        beginning of the inner loop.

    L200:

      //           determine the levenberg-marquardt parameter.

      lmpar(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], delta,
	    &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

      //           store the direction p and x + p. calculate the norm of p.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa1[j] = -wa1[j];
	wa2[j] = x[j] + wa1[j];
	wa3[j] = diag[j] * wa1[j];
	// L210:
      }
      pnorm = enorm(n, &wa3[1]);

      //           on the first iteration, adjust the initial step bound.

      // dtn
      // If any of the parameter, wa2,
      // is outside the open interval deal with it
      for ( j = 1; j <= n; ++j )
	wa2[ j ] = std::max( low[ j - 1 ], std::min( wa2[ j ], high[ j - 1 ] ) );
      // dtn

      if (iter == 1) {
	delta = std::min(delta,pnorm);
      }

      //           evaluate the function at x + p and calculate its norm.

      fcn(m, n, &wa2[1], &wa4[1], iflag, xptr);
      ++nfev;
      if (iflag < 0) {
	goto L300;
      }
      fnorm1 = enorm(m, &wa4[1]);

      //           compute the scaled actual reduction.

      actred = -1.;
      if (p1 * fnorm1 < fnorm) {
	// Computing 2nd power
	d__1 = fnorm1 / fnorm;
	actred = 1. - d__1 * d__1;
      }

      //           compute the scaled predicted reduction and
      //           the scaled directional derivative.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa3[j] = 0.;
	l = ipvt[j];
	temp = wa1[l];
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  wa3[i__] += fjac[i__ + j * fjac_dim1] * temp;
	  // L220:
	}
	// L230:
      }
      temp1 = enorm(n, &wa3[1]) / fnorm;
      temp2 = sqrt(par) * pnorm / fnorm;
      // Computing 2nd power
      d__1 = temp1;
      // Computing 2nd power
      d__2 = temp2;
      prered = d__1 * d__1 + d__2 * d__2 / p5;
      // Computing 2nd power
      d__1 = temp1;
      // Computing 2nd power
      d__2 = temp2;
      dirder = -(d__1 * d__1 + d__2 * d__2);

      //           compute the ratio of the actual to the predicted
      //           reduction.

      ratio = 0.;
      if (prered != 0.) {
	ratio = actred / prered;
      }

      //           update the step bound.

      if (ratio > p25) {
	goto L240;
      }
      if (actred >= 0.) {
	temp = p5;
      }
      if (actred < 0.) {
	temp = p5 * dirder / (dirder + p5 * actred);
      }
      if (p1 * fnorm1 >= fnorm || temp < p1) {
	temp = p1;
      }
      // Computing MIN
      d__1 = delta, d__2 = pnorm / p1;
      delta = temp * std::min(d__1,d__2);
      par /= temp;
      goto L260;
    L240:
      if (par != 0. && ratio < p75) {
	goto L250;
      }
      delta = pnorm / p5;
      par = p5 * par;
    L250:
    L260:

      //           test for successful iteration.

      if (ratio < p0001) {
	goto L290;
      }

      //           successful iteration. update x, fvec, and their norms.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	x[j] = wa2[j];
	wa2[j] = diag[j] * x[j];
	// L270:
      }
      i__1 = m;
      for (i__ = 1; i__ <= i__1; ++i__) {
	fvec[i__] = wa4[i__];
	// L280:
      }
      xnorm = enorm(n, &wa2[1]);
      fnorm = fnorm1;
      ++iter;
    L290:

      //           tests for convergence.

      if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
	info = 1;
      }
      if (delta <= xtol * xnorm) {
	info = 2;
      }
      if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info 
	  == 2) {
	info = 3;
      }
      if (info != 0) {
	goto L300;
      }

      //           tests for termination and stringent tolerances.

      if (nfev >= maxfev) {
	info = 5;
      }
      if (fabs(actred) <= epsmch && prered <= epsmch && p5 * ratio <= 1.) {
	info = 6;
      }
      if (delta <= epsmch * xnorm) {
	info = 7;
      }
      if (gnorm <= epsmch) {
	info = 8;
      }
      if (info != 0) {
	goto L300;
      }

      //           end of the inner loop. repeat if iteration unsuccessful.

      if (ratio < p0001) {
	goto L200;
      }

      //        end of the outer loop.

      goto L30;
    L300:

      //     termination, either normal or user imposed.

      if (iflag < 0) {
	info = iflag;
      }
      iflag = 0;
      if (nprint > 0) {
	fcn(m, n, &x[1], &fvec[1], iflag, xptr);
      }
      return info;

      //     last card of subroutine lmdif.

    } // lmdif


    //
    // c     **********
    // c
    // c     subroutine lmpar
    // c
    // c     given an m by n matrix a, an n by n nonsingular diagonal
    // c     matrix d, an m-vector b, and a positive number delta,
    // c     the problem is to determine a value for the parameter
    // c     par such that if x solves the system
    // c
    // c           a*x = b ,     sqrt(par)*d*x = 0 ,
    // c
    // c     in the least squares sense, and dxnorm is the euclidean
    // c     norm of d*x, then either par is zero and
    // c
    // c           (dxnorm-delta) .le. 0.1*delta ,
    // c
    // c     or par is positive and
    // c
    // c           abs(dxnorm-delta) .le. 0.1*delta .
    // c
    // c     this subroutine completes the solution of the problem
    // c     if it is provided with the necessary information from the
    // c     qr factorization, with column pivoting, of a. that is, if
    // c     a*p = q*r, where p is a permutation matrix, q has orthogonal
    // c     columns, and r is an upper triangular matrix with diagonal
    // c     elements of nonincreasing magnitude, then lmpar expects
    // c     the full upper triangle of r, the permutation matrix p,
    // c     and the first n components of (q transpose)*b. on output
    // c     lmpar also provides an upper triangular matrix s such that
    // c
    // c            t   t                   t
    // c           p *(a *a + par*d*d)*p = s *s .
    // c
    // c     s is employed within lmpar and may be of separate interest.
    // c
    // c     only a few iterations are generally needed for convergence
    // c     of the algorithm. if, however, the limit of 10 iterations
    // c     is reached, then the output par will contain the best
    // c     value obtained so far.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
    // c                        wa1,wa2)
    // c
    // c     where
    // c
    // c       n is a positive integer input variable set to the order of r.
    // c
    // c       r is an n by n array. on input the full upper triangle
    // c         must contain the full upper triangle of the matrix r.
    // c         on output the full upper triangle is unaltered, and the
    // c         strict lower triangle contains the strict upper triangle
    // c         (transposed) of the upper triangular matrix s.
    // c
    // c       ldr is a positive integer input variable not less than n
    // c         which specifies the leading dimension of the array r.
    // c
    // c       ipvt is an integer input array of length n which defines the
    // c         permutation matrix p such that a*p = q*r. column j of p
    // c         is column ipvt(j) of the identity matrix.
    // c
    // c       diag is an input array of length n which must contain the
    // c         diagonal elements of the matrix d.
    // c
    // c       qtb is an input array of length n which must contain the first
    // c         n elements of the vector (q transpose)*b.
    // c
    // c       delta is a positive input variable which specifies an upper
    // c         bound on the euclidean norm of d*x.
    // c
    // c       par is a nonnegative variable. on input par contains an
    // c         initial estimate of the levenberg-marquardt parameter.
    // c         on output par contains the final estimate.
    // c
    // c       x is an output array of length n which contains the least
    // c         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
    // c         for the output par.
    // c
    // c       sdiag is an output array of length n which contains the
    // c         diagonal elements of the upper triangular matrix s.
    // c
    // c       wa1 and wa2 are work arrays of length n.
    // c
    // c     subprograms called
    // c
    // c       minpack-supplied ... dpmpar,enorm,qrsolv
    // c
    // c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    //
    void lmpar(int n, double *r__, int ldr, const int *ipvt,
	       const double *diag, const double *qtb, double delta, 
	       double *par, double *x, double *sdiag, double *wa1, 
	       double *wa2) {

      // Initialized data

      const double p1=0.1;
      const double p001=0.001;

      // System generated locals
      int r_dim1, r_offset, i__1, i__2;
      double d__1, d__2;

      // Local variables
      int i__, j, k, l;
      double fp;
      int jm1, jp1;
      double sum, parc, parl;
      int iter;
      double temp, paru, dwarf;
      int nsing;
      double gnorm;
      double dxnorm;

      // Parameter adjustments
      --wa2;
      --wa1;
      --sdiag;
      --x;
      --qtb;
      --diag;
      --ipvt;
      r_dim1 = ldr;
      r_offset = 1 + r_dim1 * 1;
      r__ -= r_offset;

      // Function Body

      //     dwarf is the smallest positive magnitude.

      dwarf = std::numeric_limits< double >::min();

      //     compute and store in x the gauss-newton direction. if the
      //     jacobian is rank-deficient, obtain a least squares solution.

      nsing = n;
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa1[j] = qtb[j];
	if (r__[j + j * r_dim1] == 0. && nsing == n) {
	  nsing = j - 1;
	}
	if (nsing < n) {
	  wa1[j] = 0.;
	}
	// L10:
      }
      if (nsing < 1) {
	goto L50;
      }
      i__1 = nsing;
      for (k = 1; k <= i__1; ++k) {
	j = nsing - k + 1;
	wa1[j] /= r__[j + j * r_dim1];
	temp = wa1[j];
	jm1 = j - 1;
	if (jm1 < 1) {
	  goto L30;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	  // L20:
	}
      L30:
	// L40:
	;
      }
    L50:
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
	// L60:
      }

      //     initialize the iteration counter.
      //     evaluate the function at the origin, and test
      //     for acceptance of the gauss-newton direction.

      iter = 0;
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa2[j] = diag[j] * x[j];
	// L70:
      }
      dxnorm = enorm(n, &wa2[1]);
      fp = dxnorm - delta;
      if (fp <= p1 * delta) {
	goto L220;
      }

      //     if the jacobian is not rank deficient, the newton
      //     step provides a lower bound, parl, for the zero of
      //     the function. otherwise set this bound to zero.

      parl = 0.;
      if (nsing < n) {
	goto L120;
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
	// L80:
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	  goto L100;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  sum += r__[i__ + j * r_dim1] * wa1[i__];
	  // L90:
	}
      L100:
	wa1[j] = (wa1[j] - sum) / r__[j + j * r_dim1];
	// L110:
      }
      temp = enorm(n, &wa1[1]);
      parl = fp / delta / temp / temp;
    L120:

      //     calculate an upper bound, paru, for the zero of the function.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	sum = 0.;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  sum += r__[i__ + j * r_dim1] * qtb[i__];
	  // L130:
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
	// L140:
      }
      gnorm = enorm(n, &wa1[1]);
      paru = gnorm / delta;
      if (paru == 0.) {
	paru = dwarf / std::min(delta,p1);
      }

      //     if the input par lies outside of the interval (parl,paru),
      //     set par to the closer endpoint.

      *par = std::max(*par,parl);
      *par = std::min(*par,paru);
      if (*par == 0.) {
	*par = gnorm / dxnorm;
      }

      //     beginning of an iteration.

    L150:
      ++iter;

      //        evaluate the function at the current value of par.

      if (*par == 0.) {
	// Computing MAX
	d__1 = dwarf, d__2 = p001 * paru;
	*par = std::max(d__1,d__2);
      }
      temp = sqrt(*par);
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa1[j] = temp * diag[j];
	// L160:
      }
      qrsolv(n, &r__[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], &sdiag[
									       1], &wa2[1]);
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa2[j] = diag[j] * x[j];
	// L170:
      }
      dxnorm = enorm(n, &wa2[1]);
      temp = fp;
      fp = dxnorm - delta;

      //        if the function is small enough, accept the current value
      //        of par. also test for the exceptional cases where parl
      //        is zero or the number of iterations has reached 10.

      if (fabs(fp) <= p1 * delta || (parl == 0. && fp <= temp && temp < 0.) ||
	  iter == 10) {
	goto L220;
      }

      //        compute the newton correction.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
	// L180:
      }
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	wa1[j] /= sdiag[j];
	temp = wa1[j];
	jp1 = j + 1;
	if (n < jp1) {
	  goto L200;
	}
	i__2 = n;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	  wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	  // L190:
	}
      L200:
	// L210:
	;
      }
      temp = enorm(n, &wa1[1]);
      parc = fp / delta / temp / temp;

      //        depending on the sign of the function, update parl or paru.

      if (fp > 0.) {
	parl = std::max(parl,*par);
      }
      if (fp < 0.) {
	paru = std::min(paru,*par);
      }

      //        compute an improved estimate for par.

      // Computing MAX
      d__1 = parl, d__2 = *par + parc;
      *par = std::max(d__1,d__2);

      //        end of an iteration.

      goto L150;
    L220:

      //     termination.

      if (iter == 0) {
	*par = 0.;
      }
      return;

      //     last card of subroutine lmpar.

    } // lmpar

    //
    // c     **********
    // c
    // c     subroutine qrfac
    // c
    // c     this subroutine uses householder transformations with column
    // c     pivoting (optional) to compute a qr factorization of the
    // c     m by n matrix a. that is, qrfac determines an orthogonal
    // c     matrix q, a permutation matrix p, and an upper trapezoidal
    // c     matrix r with diagonal elements of nonincreasing magnitude,
    // c     such that a*p = q*r. the householder transformation for
    // c     column k, k = 1,2,...,min(m,n), is of the form
    // c
    // c                           t
    // c           i - (1/u(k))*u*u
    // c
    // c     where u has zeros in the first k-1 positions. the form of
    // c     this transformation and the method of pivoting first
    // c     appeared in the corresponding linpack subroutine.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
    // c
    // c     where
    // c
    // c       m is a positive integer input variable set to the number
    // c         of rows of a.
    // c
    // c       n is a positive integer input variable set to the number
    // c         of columns of a.
    // c
    // c       a is an m by n array. on input a contains the matrix for
    // c         which the qr factorization is to be computed. on output
    // c         the strict upper trapezoidal part of a contains the strict
    // c         upper trapezoidal part of r, and the lower trapezoidal
    // c         part of a contains a factored form of q (the non-trivial
    // c         elements of the u vectors described above).
    // c
    // c       lda is a positive integer input variable not less than m
    // c         which specifies the leading dimension of the array a.
    // c
    // c       pivot is a logical input variable. if pivot is set true,
    // c         then column pivoting is enforced. if pivot is set false,
    // c         then no column pivoting is done.
    // c
    // c       ipvt is an integer output array of length lipvt. ipvt
    // c         defines the permutation matrix p such that a*p = q*r.
    // c         column j of p is column ipvt(j) of the identity matrix.
    // c         if pivot is false, ipvt is not referenced.
    // c
    // c       lipvt is a positive integer input variable. if pivot is false,
    // c         then lipvt may be as small as 1. if pivot is true, then
    // c         lipvt must be at least n.
    // c
    // c       rdiag is an output array of length n which contains the
    // c         diagonal elements of r.
    // c
    // c       acnorm is an output array of length n which contains the
    // c         norms of the corresponding columns of the input matrix a.
    // c         if this information is not needed, then acnorm can coincide
    // c         with rdiag.
    // c
    // c       wa is a work array of length n. if pivot is false, then wa
    // c         can coincide with rdiag.
    // c
    // c     subprograms called
    // c
    // c       minpack-supplied ... dpmpar,enorm
    // c
    // c       fortran-supplied ... dmax1,dsqrt,min0
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    // 
    void qrfac(int m, int n, double *a, int lda, int pivot,
	       int *ipvt, int lipvt, double *rdiag,
	       double *acnorm, double *wa) {

      // Initialized data

      const double p05=0.05;

      // System generated locals
      int a_dim1, a_offset, i__1, i__2, i__3;
      double d__1, d__2, d__3;

      // Local variables
      int i__, j, k, jp1;
      double sum;
      int kmax;
      double temp;
      int minmn;
      double epsmch;
      double ajnorm;

      // Parameter adjustments
      --wa;
      --acnorm;
      --rdiag;
      a_dim1 = lda;
      a_offset = 1 + a_dim1 * 1;
      a -= a_offset;
      --ipvt;

      // Function Body

      //     epsmch is the machine precision.

      epsmch = std::numeric_limits< double >::epsilon( );

      //     compute the initial column norms and initialize several arrays.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	acnorm[j] = enorm(m, &a[j * a_dim1 + 1]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (pivot) {
	  ipvt[j] = j;
	}
	// L10:
      }

      //     reduce a to r with householder transformations.

      minmn = std::min(m,n);
      i__1 = minmn;
      for (j = 1; j <= i__1; ++j) {
	if (! (pivot)) {
	  goto L40;
	}

	//        bring the column of largest norm into the pivot position.

	kmax = j;
	i__2 = n;
	for (k = j; k <= i__2; ++k) {
	  if (rdiag[k] > rdiag[kmax]) {
	    kmax = k;
	  }
	  // L20:
	}
	if (kmax == j) {
	  goto L40;
	}
	i__2 = m;
	for (i__ = 1; i__ <= i__2; ++i__) {
	  temp = a[i__ + j * a_dim1];
	  a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
	  a[i__ + kmax * a_dim1] = temp;
	  // L30:
	}
	rdiag[kmax] = rdiag[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;
      L40:

	//        compute the householder transformation to reduce the
	//        j-th column of a to a multiple of the j-th unit vector.

	i__2 = m - j + 1;
	ajnorm = enorm(i__2, &a[j + j * a_dim1]);
	if (ajnorm == 0.) {
	  goto L100;
	}
	if (a[j + j * a_dim1] < 0.) {
	  ajnorm = -ajnorm;
	}
	i__2 = m;
	for (i__ = j; i__ <= i__2; ++i__) {
	  a[i__ + j * a_dim1] /= ajnorm;
	  // L50:
	}
	a[j + j * a_dim1] += 1.;

	//        apply the transformation to the remaining columns
	//        and update the norms.

	jp1 = j + 1;
	if (n < jp1) {
	  goto L100;
	}
	i__2 = n;
	for (k = jp1; k <= i__2; ++k) {
	  sum = 0.;
	  i__3 = m;
	  for (i__ = j; i__ <= i__3; ++i__) {
	    sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
	    // L60:
	  }
	  temp = sum / a[j + j * a_dim1];
	  i__3 = m;
	  for (i__ = j; i__ <= i__3; ++i__) {
	    a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
	    // L70:
	  }
	  if (! (pivot) || rdiag[k] == 0.) {
	    goto L80;
	  }
	  temp = a[j + k * a_dim1] / rdiag[k];
	  // Computing MAX
	  // Computing 2nd power
	  d__3 = temp;
	  d__1 = 0., d__2 = 1. - d__3 * d__3;
	  rdiag[k] *= sqrt((std::max(d__1,d__2)));
	  // Computing 2nd power
	  d__1 = rdiag[k] / wa[k];
	  if (p05 * (d__1 * d__1) > epsmch) {
	    goto L80;
	  }
	  i__3 = m - j;
	  rdiag[k] = enorm(i__3, &a[jp1 + k * a_dim1]);
	  wa[k] = rdiag[k];
	L80:
	  // L90:
	  ;
	}
      L100:
	rdiag[j] = -ajnorm;
	// L110:
      }
      return;

      //     last card of subroutine qrfac.

    } // qrfac

    //
    // c     **********
    // c
    // c     subroutine qrsolv
    // c
    // c     given an m by n matrix a, an n by n diagonal matrix d,
    // c     and an m-vector b, the problem is to determine an x which
    // c     solves the system
    // c
    // c           a*x = b ,     d*x = 0 ,
    // c
    // c     in the least squares sense.
    // c
    // c     this subroutine completes the solution of the problem
    // c     if it is provided with the necessary information from the
    // c     qr factorization, with column pivoting, of a. that is, if
    // c     a*p = q*r, where p is a permutation matrix, q has orthogonal
    // c     columns, and r is an upper triangular matrix with diagonal
    // c     elements of nonincreasing magnitude, then qrsolv expects
    // c     the full upper triangle of r, the permutation matrix p,
    // c     and the first n components of (q transpose)*b. the system
    // c     a*x = b, d*x = 0, is then equivalent to
    // c
    // c                  t       t
    // c           r*z = q *b ,  p *d*p*z = 0 ,
    // c
    // c     where x = p*z. if this system does not have full rank,
    // c     then a least squares solution is obtained. on output qrsolv
    // c     also provides an upper triangular matrix s such that
    // c
    // c            t   t               t
    // c           p *(a *a + d*d)*p = s *s .
    // c
    // c     s is computed within qrsolv and may be of separate interest.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)
    // c
    // c     where
    // c
    // c       n is a positive integer input variable set to the order of r.
    // c
    // c       r is an n by n array. on input the full upper triangle
    // c         must contain the full upper triangle of the matrix r.
    // c         on output the full upper triangle is unaltered, and the
    // c         strict lower triangle contains the strict upper triangle
    // c         (transposed) of the upper triangular matrix s.
    // c
    // c       ldr is a positive integer input variable not less than n
    // c         which specifies the leading dimension of the array r.
    // c
    // c       ipvt is an integer input array of length n which defines the
    // c         permutation matrix p such that a*p = q*r. column j of p
    // c         is column ipvt(j) of the identity matrix.
    // c
    // c       diag is an input array of length n which must contain the
    // c         diagonal elements of the matrix d.
    // c
    // c       qtb is an input array of length n which must contain the first
    // c         n elements of the vector (q transpose)*b.
    // c
    // c       x is an output array of length n which contains the least
    // c         squares solution of the system a*x = b, d*x = 0.
    // c
    // c       sdiag is an output array of length n which contains the
    // c         diagonal elements of the upper triangular matrix s.
    // c
    // c       wa is a work array of length n.
    // c
    // c     subprograms called
    // c
    // c       fortran-supplied ... dabs,dsqrt
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********
    //
    void qrsolv(int n, double *r__, int ldr, 
		const int *ipvt, const double *diag, const double *qtb, double *x, 
		double *sdiag, double *wa) {
      // Initialized data

      const double p5=0.5;
      const double p25=0.25;

      // System generated locals
      int r_dim1, r_offset, i__1, i__2, i__3;
      double d__1, d__2;

      // Local variables
      int i__, j, k, l, jp1, kp1;
      double tan__, cos__, sin__, sum, temp, cotan;
      int nsing;
      double qtbpj;

      // Parameter adjustments
      --wa;
      --sdiag;
      --x;
      --qtb;
      --diag;
      --ipvt;
      r_dim1 = ldr;
      r_offset = 1 + r_dim1 * 1;
      r__ -= r_offset;

      // Function Body

      //     copy r and (q transpose)*b to preserve input and initialize s.
      //     in particular, save the diagonal elements of r in x.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	i__2 = n;
	for (i__ = j; i__ <= i__2; ++i__) {
	  r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
	  // L10:
	}
	x[j] = r__[j + j * r_dim1];
	wa[j] = qtb[j];
	// L20:
      }

      //     eliminate the diagonal matrix d using a givens rotation.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {

	//        prepare the row of d to be eliminated, locating the
	//        diagonal element using p from the qr factorization.

	l = ipvt[j];
	if (diag[l] == 0.) {
	  goto L90;
	}
	i__2 = n;
	for (k = j; k <= i__2; ++k) {
	  sdiag[k] = 0.;
	  // L30:
	}
	sdiag[j] = diag[l];

	//        the transformations to eliminate the row of d
	//        modify only a single element of (q transpose)*b
	//        beyond the first n, which is initially zero.

	qtbpj = 0.;
	i__2 = n;
	for (k = j; k <= i__2; ++k) {

	  //           determine a givens rotation which eliminates the
	  //           appropriate element in the current row of d.

	  if (sdiag[k] == 0.) {
	    goto L70;
	  }
	  if ((d__1 = r__[k + k * r_dim1], fabs(d__1)) >= (d__2 = sdiag[k], 
							  fabs(d__2))) {
	    goto L40;
	  }
	  cotan = r__[k + k * r_dim1] / sdiag[k];
	  // Computing 2nd power
	  d__1 = cotan;
	  sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	  cos__ = sin__ * cotan;
	  goto L50;
	L40:
	  tan__ = sdiag[k] / r__[k + k * r_dim1];
	  // Computing 2nd power
	  d__1 = tan__;
	  cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	  sin__ = cos__ * tan__;
	L50:

	  //           compute the modified diagonal element of r and
	  //           the modified element of ((q transpose)*b,0).

	  r__[k + k * r_dim1] = cos__ * r__[k + k * r_dim1] + sin__ * sdiag[
									    k];
	  temp = cos__ * wa[k] + sin__ * qtbpj;
	  qtbpj = -sin__ * wa[k] + cos__ * qtbpj;
	  wa[k] = temp;

	  //           accumulate the tranformation in the row of s.

	  kp1 = k + 1;
	  if (n < kp1) {
	    goto L70;
	  }
	  i__3 = n;
	  for (i__ = kp1; i__ <= i__3; ++i__) {
	    temp = cos__ * r__[i__ + k * r_dim1] + sin__ * sdiag[i__];
	    sdiag[i__] = -sin__ * r__[i__ + k * r_dim1] + cos__ * sdiag[
									i__];
	    r__[i__ + k * r_dim1] = temp;
	    // L60:
	  }
	L70:
	  // L80:
	  ;
	}
      L90:

	//        store the diagonal element of s and restore
	//        the corresponding diagonal element of r.

	sdiag[j] = r__[j + j * r_dim1];
	r__[j + j * r_dim1] = x[j];
	// L100:
      }

      //     solve the triangular system for z. if the system is
      //     singular, then obtain a least squares solution.

      nsing = n;
      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	if (sdiag[j] == 0. && nsing == n) {
	  nsing = j - 1;
	}
	if (nsing < n) {
	  wa[j] = 0.;
	}
	// L110:
      }
      if (nsing < 1) {
	goto L150;
      }
      i__1 = nsing;
      for (k = 1; k <= i__1; ++k) {
	j = nsing - k + 1;
	sum = 0.;
	jp1 = j + 1;
	if (nsing < jp1) {
	  goto L130;
	}
	i__2 = nsing;
	for (i__ = jp1; i__ <= i__2; ++i__) {
	  sum += r__[i__ + j * r_dim1] * wa[i__];
	  // L120:
	}
      L130:
	wa[j] = (wa[j] - sum) / sdiag[j];
	// L140:
      }
    L150:

      //     permute the components of z back to components of x.

      i__1 = n;
      for (j = 1; j <= i__1; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
	// L160:
      }
      return;

      //     last card of subroutine qrsolv.

    } // qrsolv

  }; // class

  
} // namespace

#endif


