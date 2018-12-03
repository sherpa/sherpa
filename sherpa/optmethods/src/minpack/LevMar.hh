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
//
// @misc{cminpack,
//   title={C/C++ Minpack},
//   author={Devernay, Fr{\'e}d{\'e}ric},
//   year={2007},
//   howpublished = "\url{http://devernay.free.fr/hacks/cminpack/}",
// }
//

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


  template <typename Func, typename Data, typename real >
  class LevMar : public sherpa::Opt<Data, real> {

  public:
    LevMar( Func func, Data xdata ) : sherpa::Opt<Data, real>( xdata ),
                               usr_func( func ) { }

  protected:

    Func usr_func;
    Func get_usr_func( ) { return usr_func; }

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
    // c       real precision function enorm(n,x)
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


    void print_progress( int m, int n, const real* x, const real* fvec )
      const {
      const real fval = pow( this->enorm(m, fvec), 2.0 );
      std::cout << "f( " << x[ 0 ];
      for ( int iii = 1; iii < n; ++iii )
        std::cout << ", " << x[ iii ];
      std::cout << " ) = " << fval << '\n';
    }

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
    void covar( int n, real *r, int ldr, const int *ipvt, real tol,
		real *wa ) {

      /* Local variables */
      int i, j, k, l, ii, jj;
      int sing;
      real temp, tolr;

      tolr = tol * fabs(r[0]);

      /*     form the inverse of r in the full upper triangle of r. */

      l = -1;
      for (k = 0; k < n; ++k) {
        if (fabs(r[k + k * ldr]) <= tolr) {
          break;
        }
        r[k + k * ldr] = 1. / r[k + k * ldr];
        if (k > 0) {
          for (j = 0; j < k; ++j) {
            // coverity[copy_paste_error]
            temp = r[k + k * ldr] * r[j + k * ldr];
            r[j + k * ldr] = 0.;
            for (i = 0; i <= j; ++i) {
              r[i + k * ldr] -= temp * r[i + j * ldr];
            }
          }
        }
        l = k;
      }

      /*     form the full upper triangle of the inverse of (r transpose)*r */
      /*     in the full upper triangle of r. */

      if (l >= 0) {
        for (k = 0; k <= l; ++k) {
          if (k > 0) {
            for (j = 0; j < k; ++j) {
              temp = r[j + k * ldr];
              for (i = 0; i <= j; ++i) {
                r[i + j * ldr] += temp * r[i + k * ldr];
              }
            }
          }
          temp = r[k + k * ldr];
          for (i = 0; i <= k; ++i) {
            r[i + k * ldr] *= temp;
          }
        }
      }

      /*     form the full lower triangle of the covariance matrix */
      /*     in the strict lower triangle of r and in wa. */

      for (j = 0; j < n; ++j) {
        jj = ipvt[j]-1;
        sing = j > l;
        for (i = 0; i <= j; ++i) {
          if (sing) {
            r[i + j * ldr] = 0.;
          }
          ii = ipvt[i]-1;
          if (ii > jj) {
            r[ii + jj * ldr] = r[i + j * ldr];
          }
          else if (ii < jj) {
            r[jj + ii * ldr] = r[i + j * ldr];
          }
        }
        wa[jj] = r[j + j * ldr];
      }

      /*     symmetrize the covariance matrix in r. */

      for (j = 0; j < n; ++j) {
        for (i = 0; i < j; ++i) {
          r[i + j * ldr] = r[j + i * ldr];
        }
        r[j + j * ldr] = wa[j];

      }      //     last card of subroutine covar.

    } // covar

    real enorm( int n, const real *x) const {

      // Initialized data
/*
  About the values for rdwarf and rgiant.
  The original values, both in single-precision FORTRAN source code and in double-precision code were:
#define rdwarf 3.834e-20
#define rgiant 1.304e19
  See for example:
    http://www.netlib.org/slatec/src/denorm.f
    http://www.netlib.org/slatec/src/enorm.f
  However, rdwarf is smaller than sqrt(FLT_MIN) = 1.0842021724855044e-19, so that rdwarf**2 will
  underflow. This contradicts the constraints expressed in the comments below.
  We changed these constants to those proposed by the
  implementation found in MPFIT http://cow.physics.wisc.edu/~craigm/idl/fitting.html
 cmpfit-1.2 proposes the following definitions:
  rdwarf = sqrt(dpmpar(2)*1.5) * 10
  rgiant = sqrt(dpmpar(3)) * 0.1
 The half version does not really worked that way, so we use for half:
  rdwarf = sqrt(dpmpar(2)) * 2
  rgiant = sqrt(dpmpar(3)) * 0.5
 Any suggestion is welcome. Half CMINPACK is really only a
 proof-of-concept anyway.
 See the example/tenorm*c, which computes these values
*/
      // const real rdwarf=3.834e-20;
      // const real rgiant= 1.304e19;

      const real rdwarf = sqrt(std::numeric_limits< real >::min() * 1.5) * 10;
      const real rgiant = sqrt(std::numeric_limits< real >::max()) / 10;
#ifdef USE_CBLAS
    return cblas_dnrm2(n, x, 1);
#else /* !USE_CBLAS */
    /* System generated locals */
      real ret_val, d1;

      /* Local variables */
      real s1, s2, s3, xabs, x1max, x3max, agiant;

      /*     ********** */

      /*     function enorm */

      /*     given an n-vector x, this function calculates the */
      /*     euclidean norm of x. */

      /*     the euclidean norm is computed by accumulating the sum of */
      /*     squares in three different sums. the sums of squares for the */
      /*     small and large components are scaled so that no overflows */
      /*     occur. non-destructive underflows are permitted. underflows */
      /*     and overflows do not occur in the computation of the unscaled */
      /*     sum of squares for the intermediate components. */
      /*     the definitions of small, intermediate and large components */
      /*     depend on two constants, rdwarf and rgiant. the main */
      /*     restrictions on these constants are that rdwarf**2 not */
      /*     underflow and rgiant**2 not overflow. the constants */
      /*     given here are suitable for every known computer. */

      /*     the function statement is */

      /*       real precision function enorm(n,x) */

      /*     where */

      /*       n is a positive integer input variable. */

      /*       x is an input array of length n. */

      /*     subprograms called */

      /*       fortran-supplied ... dabs,dsqrt */

      /*     argonne national laboratory. minpack project. march 1980. */
      /*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

      /*     ********** */

      s1 = 0.;
      s2 = 0.;
      s3 = 0.;
      x1max = 0.;
      x3max = 0.;
      agiant = rgiant / (real)n;
      for (int i = 0; i < n; ++i) {
	xabs = fabs(x[i]);
        if (xabs >= agiant) {
          /*              sum for large components. */
          if (xabs > x1max) {
            /* Computing 2nd power */
            d1 = x1max / xabs;
            s1 = 1. + s1 * (d1 * d1);
            x1max = xabs;
          } else {
            /* Computing 2nd power */
            d1 = xabs / x1max;
            s1 += d1 * d1;
          }
        } else if (xabs <= rdwarf) {
          /*              sum for small components. */
          if (xabs > x3max) {
            /* Computing 2nd power */
            d1 = x3max / xabs;
            s3 = 1. + s3 * (d1 * d1);
            x3max = xabs;
          } else if (xabs != 0.) {
            /* Computing 2nd power */
            d1 = xabs / x3max;
            s3 += d1 * d1;
          }
	} else {
          /*           sum for intermediate components. */
          /* Computing 2nd power */
          s2 += xabs * xabs;
        }
      }

      /*     calculation of norm. */

      if (s1 != 0.) {
        ret_val = x1max * sqrt(s1 + (s2 / x1max) / x1max);
      } else if (s2 != 0.) {
        if (s2 >= x3max) {
          ret_val = sqrt(s2 * (1. + (x3max / s2) * (x3max * s3)));
        } else {
          ret_val = sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
        }
      } else {
        ret_val = x3max * sqrt(s3);
      }
      return ret_val;

      //     last card of function enorm.
#endif /* !USE_CBLAS */
    } // enorm


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
    void lmpar(int n, real *r, int ldr, const int *ipvt,
	       const real *diag, const real *qtb, real delta,
	       real *par, real *x, real *sdiag, real *wa1,
	       real *wa2) {

      // Initialized data

      const real p1=0.1;
      const real p001=0.001;

      real dwarf = std::numeric_limits< real >::min();
      /* System generated locals */
      real d1, d2;

      /* Local variables */
      int j, l;
      real fp;
      real parc, parl;
      int iter;
      real temp, paru;
      int nsing;
      real gnorm;
      real dxnorm;


      /*     compute and store in x the gauss-newton direction. if the */
      /*     jacobian is rank-deficient, obtain a least squares solution. */

      nsing = n;
      for (j = 0; j < n; ++j) {
	wa1[j] = qtb[j];
	if (r[j + j * ldr] == 0. && nsing == n) {
          nsing = j;
	}
	if (nsing < n) {
          wa1[j] = 0.;
	}
      }
# ifdef USE_CBLAS
      cblas_dtrsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, nsing, r, ldr, wa1, 1);
# else
      if (nsing >= 1) {
        int k;
        for (k = 1; k <= nsing; ++k) {
          j = nsing - k;
          wa1[j] /= r[j + j * ldr];
          temp = wa1[j];
          if (j >= 1) {
            int i;
            for (i = 0; i < j; ++i) {
              wa1[i] -= r[i + j * ldr] * temp;
            }
          }
        }
      }
# endif
      for (j = 0; j < n; ++j) {
	l = ipvt[j]-1;
	x[l] = wa1[j];
      }

      /*     initialize the iteration counter. */
      /*     evaluate the function at the origin, and test */
      /*     for acceptance of the gauss-newton direction. */

      iter = 0;
      for (j = 0; j < n; ++j) {
	wa2[j] = diag[j] * x[j];
      }
      dxnorm = this->enorm(n, wa2);
      fp = dxnorm - delta;
      if (fp <= p1 * delta) {
	goto TERMINATE;
      }

      /*     if the jacobian is not rank deficient, the newton */
      /*     step provides a lower bound, parl, for the zero of */
      /*     the function. otherwise set this bound to zero. */

      parl = 0.;
      if (nsing >= n) {
        for (j = 0; j < n; ++j) {
          l = ipvt[j]-1;
          wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
#     ifdef USE_CBLAS
        cblas_dtrsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, r, ldr, wa1, 1);
#     else
        for (j = 0; j < n; ++j) {
          real sum = 0.;
          if (j >= 1) {
            int i;
            for (i = 0; i < j; ++i) {
              sum += r[i + j * ldr] * wa1[i];
            }
          }
          wa1[j] = (wa1[j] - sum) / r[j + j * ldr];
        }
#     endif
        temp = this->enorm(n, wa1);
        parl = fp / delta / temp / temp;
      }

      /*     calculate an upper bound, paru, for the zero of the function. */

      for (j = 0; j < n; ++j) {
        real sum;
#     ifdef USE_CBLAS
        sum = cblas_ddot(j+1, &r[j*ldr], 1, qtb, 1);
#     else
        int i;
        sum = 0.;
        for (i = 0; i <= j; ++i) {
          sum += r[i + j * ldr] * qtb[i];
        }
#     endif
        l = ipvt[j]-1;
        wa1[j] = sum / diag[l];
      }
      gnorm = this->enorm(n, wa1);
      paru = gnorm / delta;
      if (paru == 0.) {
        paru = dwarf / std::min(delta,(real)p1) /* / p001 ??? */;
      }

      /*     if the input par lies outside of the interval (parl,paru), */
      /*     set par to the closer endpoint. */

      *par = std::max(*par,parl);
      *par = std::min(*par,paru);
      if (*par == 0.) {
        *par = gnorm / dxnorm;
      }

      /*     beginning of an iteration. */

      for (;;) {
        ++iter;

        /*        evaluate the function at the current value of par. */

        if (*par == 0.) {
          /* Computing MAX */
          d1 = dwarf, d2 = p001 * paru;
          *par = std::max(d1,d2);
        }
        temp = sqrt(*par);
        for (j = 0; j < n; ++j) {
          wa1[j] = temp * diag[j];
        }
        qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
        for (j = 0; j < n; ++j) {
          wa2[j] = diag[j] * x[j];
        }
        dxnorm = this->enorm(n, wa2);
        temp = fp;
        fp = dxnorm - delta;

        /*        if the function is small enough, accept the current value */
        /*        of par. also test for the exceptional cases where parl */
        /*        is zero or the number of iterations has reached 10. */

        if (fabs(fp) <= p1 * delta || (parl == 0. && fp <= temp && temp < 0.) || iter == 10) {
          goto TERMINATE;
        }

        /*        compute the newton correction. */

#     ifdef USE_CBLAS
        for (j = 0; j < nsing; ++j) {
          l = ipvt[j]-1;
          wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for (j = nsing; j < n; ++j) {
          wa1[j] = 0.;
        }
        /* exchange the diagonal of r with sdiag */
        cblas_dswap(n, r, ldr+1, sdiag, 1);
        /* solve lower(r).x = wa1, result id put in wa1 */
        cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, nsing, r, ldr, wa1, 1);
        /* exchange the diagonal of r with sdiag */
        cblas_dswap( n, r, ldr+1, sdiag, 1);
#     else /* !USE_CBLAS */
        for (j = 0; j < n; ++j) {
          l = ipvt[j]-1;
          wa1[j] = diag[l] * (wa2[l] / dxnorm);
        }
        for (j = 0; j < n; ++j) {
          wa1[j] /= sdiag[j];
          temp = wa1[j];
          if (n > j+1) {
            int i;
            for (i = j+1; i < n; ++i) {
              wa1[i] -= r[i + j * ldr] * temp;
            }
          }
        }
#     endif /* !USE_CBLAS */
        temp = this->enorm(n, wa1);
        parc = fp / delta / temp / temp;

        /*        depending on the sign of the function, update parl or paru. */

        if (fp > 0.) {
          parl = std::max(parl,*par);
        }
        if (fp < 0.) {
          paru = std::min(paru,*par);
        }

        /*        compute an improved estimate for par. */

        /* Computing MAX */
        d1 = parl, d2 = *par + parc;
        *par = std::max(d1,d2);

        /*        end of an iteration. */

      }
    TERMINATE:

      /*     termination. */

      if (iter == 0) {
	*par = 0.;
      }
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
    void qrfac(int m, int n, real *a, int lda, int pivot,
	       int *ipvt, int lipvt, real *rdiag,
	       real *acnorm, real *wa) {

      // Initialized data

      const real p05=0.05;

#ifdef USE_LAPACK
      __CLPK_integer m_ = m;
      __CLPK_integer n_ = n;
      __CLPK_integer lda_ = lda;
      __CLPK_integer *jpvt;

      int i, j, k;
      real t;
      real* tau = wa;
      const __CLPK_integer ltau = m > n ? n : m;
      __CLPK_integer lwork = -1;
      __CLPK_integer info = 0;
      real* work;

      if (pivot) {
        assert( lipvt >= n );
        if (sizeof(__CLPK_integer) != sizeof(ipvt[0])) {
          jpvt = malloc(n*sizeof(__CLPK_integer));
        } else {
          /* __CLPK_integer is actually an int, just do a cast */
          jpvt = (__CLPK_integer *)ipvt;
        }
        /* set all columns free */
        memset(jpvt, 0, sizeof(int)*n);
      }

      /* query optimal size of work */
      lwork = -1;
      if (pivot) {
        dgeqp3_(&m_,&n_,a,&lda_,jpvt,tau,tau,&lwork,&info);
        lwork = (int)tau[0];
        assert( lwork >= 3*n+1  );
      } else {
        dgeqrf_(&m_,&n_,a,&lda_,tau,tau,&lwork,&info);
        lwork = (int)tau[0];
        assert( lwork >= 1 && lwork >= n );
      }

      assert( info == 0 );

      /* alloc work area */
      work = (real *)malloc(sizeof(real)*lwork);
      assert(work != NULL);

      /* set acnorm first (from the doc of qrfac, acnorm may point to the same area as rdiag) */
      if (acnorm != rdiag) {
        for (j = 0; j < n; ++j) {
          acnorm[j] = this->enorm(m, &a[j * lda]);
        }
      }

      /* QR decomposition */
      if (pivot) {
        dgeqp3_(&m_,&n_,a,&lda_,jpvt,tau,work,&lwork,&info);
      } else {
        dgeqrf_(&m_,&n_,a,&lda_,tau,work,&lwork,&info);
      }
      assert(info == 0);

      /* set rdiag, before the diagonal is replaced */
      memset(rdiag, 0, sizeof(real)*n);
      for(i=0 ; i<n ; ++i) {
        rdiag[i] = a[i*lda+i];
      }

      /* modify lower trinagular part to look like qrfac's output */
      for(i=0 ; i<ltau ; ++i) {
        k = i*lda+i;
        t = tau[i];
        a[k] = t;
        for(j=i+1 ; j<m ; j++) {
          k++;
          a[k] *= t;
        }
      }

      free(work);
      if (pivot) {
        /* convert back jpvt to ipvt */
        if (sizeof(__CLPK_integer) != sizeof(ipvt[0])) {
          for(i=0; i<n; ++i) {
            ipvt[i] = jpvt[i];
          }
          free(jpvt);
        }
      }
#else /* !USE_LAPACK */
      /* Initialized data */

      // #define p05 .05

        /* System generated locals */
      real d1;

      /* Local variables */
      int i, j, k, jp1;
      real sum;
      real temp;
      int minmn;
      real epsmch;
      real ajnorm;

      /*     ********** */

      /*     subroutine qrfac */

      /*     this subroutine uses householder transformations with column */
      /*     pivoting (optional) to compute a qr factorization of the */
      /*     m by n matrix a. that is, qrfac determines an orthogonal */
      /*     matrix q, a permutation matrix p, and an upper trapezoidal */
      /*     matrix r with diagonal elements of nonincreasing magnitude, */
      /*     such that a*p = q*r. the householder transformation for */
      /*     column k, k = 1,2,...,min(m,n), is of the form */

      /*                           t */
      /*           i - (1/u(k))*u*u */

      /*     where u has zeros in the first k-1 positions. the form of */
      /*     this transformation and the method of pivoting first */
      /*     appeared in the corresponding linpack subroutine. */

      /*     the subroutine statement is */

      /*       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) */

      /*     where */

      /*       m is a positive integer input variable set to the number */
      /*         of rows of a. */

      /*       n is a positive integer input variable set to the number */
      /*         of columns of a. */

      /*       a is an m by n array. on input a contains the matrix for */
      /*         which the qr factorization is to be computed. on output */
      /*         the strict upper trapezoidal part of a contains the strict */
      /*         upper trapezoidal part of r, and the lower trapezoidal */
      /*         part of a contains a factored form of q (the non-trivial */
      /*         elements of the u vectors described above). */

      /*       lda is a positive integer input variable not less than m */
      /*         which specifies the leading dimension of the array a. */

      /*       pivot is a logical input variable. if pivot is set true, */
      /*         then column pivoting is enforced. if pivot is set false, */
      /*         then no column pivoting is done. */

      /*       ipvt is an integer output array of length lipvt. ipvt */
      /*         defines the permutation matrix p such that a*p = q*r. */
      /*         column j of p is column ipvt(j) of the identity matrix. */
      /*         if pivot is false, ipvt is not referenced. */

      /*       lipvt is a positive integer input variable. if pivot is false, */
      /*         then lipvt may be as small as 1. if pivot is true, then */
      /*         lipvt must be at least n. */

      /*       rdiag is an output array of length n which contains the */
      /*         diagonal elements of r. */

      /*       acnorm is an output array of length n which contains the */
      /*         norms of the corresponding columns of the input matrix a. */
      /*         if this information is not needed, then acnorm can coincide */
      /*         with rdiag. */

      /*       wa is a work array of length n. if pivot is false, then wa */
      /*         can coincide with rdiag. */

      /*     subprograms called */

      /*       minpack-supplied ... dpmpar,enorm */

      /*       fortran-supplied ... dmax1,dsqrt,min0 */

      /*     argonne national laboratory. minpack project. march 1980. */
      /*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

      /*     ********** */
      (void)lipvt;

      /*     epsmch is the machine precision. */
      epsmch = std::numeric_limits< real >::epsilon( );

      /*     compute the initial column norms and initialize several arrays. */

      for (j = 0; j < n; ++j) {
	acnorm[j] = this->enorm(m, &a[j * lda + 0]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (pivot) {
          ipvt[j] = j+1;
	}
      }

      /*     reduce a to r with householder transformations. */

      minmn = std::min(m,n);
      for (j = 0; j < minmn; ++j) {
	if (pivot) {

          /*        bring the column of largest norm into the pivot position. */

          int kmax = j;
          for (k = j; k < n; ++k) {
            if (rdiag[k] > rdiag[kmax]) {
              kmax = k;
            }
          }
          if (kmax != j) {
            for (i = 0; i < m; ++i) {
              temp = a[i + j * lda];
              a[i + j * lda] = a[i + kmax * lda];
              a[i + kmax * lda] = temp;
            }
            rdiag[kmax] = rdiag[j];
            wa[kmax] = wa[j];
            k = ipvt[j];
            ipvt[j] = ipvt[kmax];
            ipvt[kmax] = k;
          }
        }

        /*        compute the householder transformation to reduce the */
        /*        j-th column of a to a multiple of the j-th unit vector. */

	ajnorm = this->enorm(m - (j+1) + 1, &a[j + j * lda]);
	if (ajnorm != 0.) {
          if (a[j + j * lda] < 0.) {
            ajnorm = -ajnorm;
          }
          for (i = j; i < m; ++i) {
            a[i + j * lda] /= ajnorm;
          }
          a[j + j * lda] += 1.;

          /*        apply the transformation to the remaining columns */
          /*        and update the norms. */

          jp1 = j + 1;
          if (n > jp1) {
            for (k = jp1; k < n; ++k) {
              sum = 0.;
              for (i = j; i < m; ++i) {
                sum += a[i + j * lda] * a[i + k * lda];
              }
              temp = sum / a[j + j * lda];
              for (i = j; i < m; ++i) {
                a[i + k * lda] -= temp * a[i + j * lda];
              }
              if (pivot && rdiag[k] != 0.) {
                temp = a[j + k * lda] / rdiag[k];
                /* Computing MAX */
                d1 = 1. - temp * temp;
                rdiag[k] *= sqrt((std::max((real)0.,d1)));
                /* Computing 2nd power */
                d1 = rdiag[k] / wa[k];
                if (p05 * (d1 * d1) <= epsmch) {
                  rdiag[k] = this->enorm(m - (j+1), &a[jp1 + k * lda]);
                  wa[k] = rdiag[k];
                }
              }
            }
          }
        }
	rdiag[j] = -ajnorm;
      }

#endif /* !USE_LAPACK */

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
    void qrsolv(int n, real *r, int ldr, const int *ipvt,
                const real *diag, const real *qtb, real *x,
		real *sdiag, real *wa) {
      // Initialized data

      const real p5=0.5;
      const real p25=0.25;

      /* Local variables */
      int i, j, k, l;
      real cos, sin, sum, temp;
      int nsing;
      real qtbpj;

      /*     ********** */

      /*     subroutine qrsolv */

      /*     given an m by n matrix a, an n by n diagonal matrix d, */
      /*     and an m-vector b, the problem is to determine an x which */
      /*     solves the system */

      /*           a*x = b ,     d*x = 0 , */

      /*     in the least squares sense. */

      /*     this subroutine completes the solution of the problem */
      /*     if it is provided with the necessary information from the */
      /*     qr factorization, with column pivoting, of a. that is, if */
      /*     a*p = q*r, where p is a permutation matrix, q has orthogonal */
      /*     columns, and r is an upper triangular matrix with diagonal */
      /*     elements of nonincreasing magnitude, then qrsolv expects */
      /*     the full upper triangle of r, the permutation matrix p, */
      /*     and the first n components of (q transpose)*b. the system */
      /*     a*x = b, d*x = 0, is then equivalent to */

      /*                  t       t */
      /*           r*z = q *b ,  p *d*p*z = 0 , */

      /*     where x = p*z. if this system does not have full rank, */
      /*     then a least squares solution is obtained. on output qrsolv */
      /*     also provides an upper triangular matrix s such that */

      /*            t   t               t */
      /*           p *(a *a + d*d)*p = s *s . */

      /*     s is computed within qrsolv and may be of separate interest. */

      /*     the subroutine statement is */

      /*       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa) */

      /*     where */

      /*       n is a positive integer input variable set to the order of r. */

      /*       r is an n by n array. on input the full upper triangle */
      /*         must contain the full upper triangle of the matrix r. */
      /*         on output the full upper triangle is unaltered, and the */
      /*         strict lower triangle contains the strict upper triangle */
      /*         (transposed) of the upper triangular matrix s. */

      /*       ldr is a positive integer input variable not less than n */
      /*         which specifies the leading dimension of the array r. */

      /*       ipvt is an integer input array of length n which defines the */
      /*         permutation matrix p such that a*p = q*r. column j of p */
      /*         is column ipvt(j) of the identity matrix. */

      /*       diag is an input array of length n which must contain the */
      /*         diagonal elements of the matrix d. */

      /*       qtb is an input array of length n which must contain the first */
      /*         n elements of the vector (q transpose)*b. */

      /*       x is an output array of length n which contains the least */
      /*         squares solution of the system a*x = b, d*x = 0. */

      /*       sdiag is an output array of length n which contains the */
      /*         diagonal elements of the upper triangular matrix s. */

      /*       wa is a work array of length n. */

      /*     subprograms called */

      /*       fortran-supplied ... dabs,dsqrt */

      /*     argonne national laboratory. minpack project. march 1980. */
      /*     burton s. garbow, kenneth e. hillstrom, jorge j. more */

      /*     ********** */

      /*     copy r and (q transpose)*b to preserve input and initialize s. */
      /*     in particular, save the diagonal elements of r in x. */

      for (j = 0; j < n; ++j) {
	for (i = j; i < n; ++i) {
          r[i + j * ldr] = r[j + i * ldr];
	}
	x[j] = r[j + j * ldr];
	wa[j] = qtb[j];
      }

      /*     eliminate the diagonal matrix d using a givens rotation. */

      for (j = 0; j < n; ++j) {

        /*        prepare the row of d to be eliminated, locating the */
        /*        diagonal element using p from the qr factorization. */

	l = ipvt[j]-1;
	if (diag[l] != 0.) {
          for (k = j; k < n; ++k) {
            sdiag[k] = 0.;
          }
          sdiag[j] = diag[l];

          /*        the transformations to eliminate the row of d */
          /*        modify only a single element of (q transpose)*b */
          /*        beyond the first n, which is initially zero. */

          qtbpj = 0.;
          for (k = j; k < n; ++k) {

            /*           determine a givens rotation which eliminates the */
            /*           appropriate element in the current row of d. */

            if (sdiag[k] != 0.) {
#                 ifdef USE_LAPACK
              dlartg_( &r[k + k * ldr], &sdiag[k], &cos, &sin, &temp );
#                 else /* !USE_LAPACK */
              if (fabs(r[k + k * ldr]) < fabs(sdiag[k])) {
                real cotan;
                cotan = r[k + k * ldr] / sdiag[k];
                sin = p5 / sqrt(p25 + p25 * (cotan * cotan));
                cos = sin * cotan;
              } else {
                real tan;
                tan = sdiag[k] / r[k + k * ldr];
                cos = p5 / sqrt(p25 + p25 * (tan * tan));
                sin = cos * tan;
              }

              /*           compute the modified diagonal element of r and */
              /*           the modified element of ((q transpose)*b,0). */

#                 endif /* !USE_LAPACK */
              temp = cos * wa[k] + sin * qtbpj;
              qtbpj = -sin * wa[k] + cos * qtbpj;
              wa[k] = temp;

              /*           accumulate the tranformation in the row of s. */
#                 ifdef USE_CBLAS
              cblas_drot( n-k, &r[k + k * ldr], 1, &sdiag[k], 1, cos, sin );
#                 else /* !USE_CBLAS */
              r[k + k * ldr] = cos * r[k + k * ldr] + sin * sdiag[k];
              if (n > k+1) {
                for (i = k+1; i < n; ++i) {
                  temp = cos * r[i + k * ldr] + sin * sdiag[i];
                  sdiag[i] = -sin * r[i + k * ldr] + cos * sdiag[i];
                  r[i + k * ldr] = temp;
                }
              }
#                 endif /* !USE_CBLAS */
            }
          }
        }

        /*        store the diagonal element of s and restore */
        /*        the corresponding diagonal element of r. */

	sdiag[j] = r[j + j * ldr];
	r[j + j * ldr] = x[j];
      }

      /*     solve the triangular system for z. if the system is */
      /*     singular, then obtain a least squares solution. */

      nsing = n;
      for (j = 0; j < n; ++j) {
	if (sdiag[j] == 0. && nsing == n) {
          nsing = j;
	}
	if (nsing < n) {
          wa[j] = 0.;
	}
      }
      if (nsing >= 1) {
        for (k = 1; k <= nsing; ++k) {
          j = nsing - k;
          sum = 0.;
          if (nsing > j+1) {
            for (i = j+1; i < nsing; ++i) {
              sum += r[i + j * ldr] * wa[i];
            }
          }
          wa[j] = (wa[j] - sum) / sdiag[j];
        }
      }

      /*     permute the components of z back to components of x. */

      for (j = 0; j < n; ++j) {
	l = ipvt[j]-1;
	x[l] = wa[j];
      }
      return;

    } // qrsolv


  };

  template < typename Func, typename Data, typename real >
  class LevMarDif : public LevMar<Func, Data, real> {

  public:

    LevMarDif( Func func, Data xdata, int mfct )
      : LevMar<Func, Data, real>( func, xdata ), myfvec( mfct ) { }

    int operator( )( int n, real ftol, real xtol,
                     real gtol, int maxfev, real epsfcn,
                     real factor, int nprint,
                     std::vector<real>& x,
                     int& nfev, real& fmin,
                     const sherpa::Bounds<real>& bounds,
                     std::vector<real>& fjac,
                     std::vector<real>& covarerr ) {

      int info = 0;

      try {

        if ( sherpa::Opt<Data, real>::are_pars_outside_limits( n, x, bounds ) )
	  throw sherpa::OptErr( sherpa::OptErr::OutOfBound );

	int m = static_cast<int>( myfvec.size( ) );

	std::vector<real> diag( n ), qtf( n ), wa1( n ), wa2( n ), wa3( n );
	std::vector<real> wa4( m );
	std::vector<int> ipvt( n );

	const int mode = 1;
	const int ldfjac = m;

        Func usr_func = LevMar<Func, Data, real>::get_usr_func( );
        Data lm_usr_data = sherpa::Opt<Data, real>::get_usr_data();
	info = lmdif( usr_func, lm_usr_data, m, n, &x[0], &myfvec[0], ftol,
                      xtol, gtol, maxfev, epsfcn, &diag[0], mode, factor,
                      nprint, nfev, &fjac[0], ldfjac, &ipvt[0], &qtf[0],
                      &wa1[ 0 ], &wa2[0], &wa3[0], &wa4[0], bounds );

        if ( info > 0 ) {
          // The covariance matrix does not need to be calculated if there was no fit.
          this->covar( n, &fjac[ 0 ], ldfjac, &ipvt[0], ftol, &wa1[0] );

          for ( int ii = 0; ii < n; ++ii )
            if ( fjac[ ii + ldfjac * ii ] > 0.0 )
              covarerr[ ii ] = sqrt( fjac[ ii + ldfjac * ii ] );
            else
              covarerr[ ii ] = 0.0;
        }

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

      fmin = std::pow( this->enorm( myfvec.size( ), &myfvec[0] ), 2.0 );
      return info;

    }

    virtual real eval_func( int maxnfev, const sherpa::Bounds<real>& limits,
                            int npar, std::vector<real>& par, int& nfev ) {

      if ( sherpa::Opt<Data, real>::are_pars_outside_limits( npar, par,
                                                             limits ) ) {
        return std::numeric_limits< real >::max( );
      }

      ++nfev;
      int ierr=EXIT_SUCCESS;

      int mymfct = static_cast<int>( myfvec.size( ) );

      Data my_usr_data = sherpa::Opt<Data, real>::get_usr_data();
      this->usr_func( mymfct, npar, &par[0], &myfvec[0], ierr, my_usr_data );

      real fval = pow( this->enorm( mymfct, &myfvec[0] ), 2.0 );
      if ( EXIT_SUCCESS != ierr )
        throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
      if ( nfev >= maxnfev )
        throw sherpa::OptErr( sherpa::OptErr::MaxFev );

      return fval;

    }


    // de
    int minimize( int maxnfev, const sherpa::Bounds<real>& limits, real tol,
                  int npar, std::vector<real>& par, real& fmin, int& nfev ) {
      int nprint = 0;
      real epsfcn =
	std::sqrt( std::numeric_limits< real >::epsilon( ) );
      real factor = 100.0;

      std::vector<real> diag( npar );
      std::vector<real> covarerr( npar );
      std::vector<real> fjac( myfvec.size() * npar );
      return this->operator( )( npar, tol, tol, tol, maxnfev, epsfcn, factor,
				nprint, par, nfev, fmin, limits, fjac,
                                covarerr );
    }
    // de

  protected:

    std::vector< real > myfvec;
    // std::vector<real>& get_myfec( ) { return myfvec; }


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
    // c         real precision x(n),fvec(m)
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
    int fdjac2( Func fcn, int m, int n, real *x,
		const real *fvec, real *fjac, int ldfjac,
		real epsfcn, real *wa, Data xptr,
		const std::vector<real>& high ) {

      /* Local variables */
      real h;
      int i, j;
      real eps, temp, epsmch;
      int iflag=0;
      epsmch = std::numeric_limits< real >::epsilon( );
      eps = sqrt((std::max(epsfcn,epsmch)));
      for (j = 0; j < n; ++j) {
	temp = x[j];
	h = eps * fabs(temp);
	if (h == 0.) {
          h = eps;
        }
        // dtn
        // if the parameter is beyond the upper boundary then
        // perform backwards-difference approximation
        if ( x[ j ] + h > high[ j ] )
          h = - h;
        // dtn
 	x[j] = temp + h;
        /* the last parameter of fcn_mn() is set to 2 to differentiate
           calls made to compute the function from calls made to compute
           the Jacobian (see fcn() in examples/lmfdrv.c, and how njev
           is used to compute the number of Jacobian evaluations) */
        fcn(m, n, x, wa, iflag, xptr);
	if (iflag < 0) {
          return iflag;
	}
	x[j] = temp;
	for (i = 0; i < m; ++i) {
          fjac[i + j * ldfjac] = (wa[i] - fvec[i]) / h;
	}
      }

      // L30:
      return iflag;

      //     last card of subroutine fdjac2.

    } // fdjac2


  private:

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
    // c         real precision x(n),fvec(m)
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
    int lmdif( Func fcn, Data xptr, int m, int n, real *x,
	       real *fvec, real ftol, real xtol, real gtol,
	       int maxfev, real epsfcn, real *diag, int mode,
	       real factor, int nprint, int& nfev, real *fjac,
	       int ldfjac, int *ipvt, real * qtf, real *wa1,
	       real *wa2, real *wa3, real *wa4,
               const sherpa::Bounds<real>& bounds ) {

      const std::vector<real>& low = bounds.get_lb();
      const std::vector<real>& high = bounds.get_ub();

      // Initialized data

      const real p1=0.1;
      const real p5=0.5;
      const real p25=0.25;
      const real p75=0.75;
      const real p0001=1.0e-4;

      /* System generated locals */
      real d1, d2;

      /* Local variables */
      int i, j, l;
      real par, sum;
      int iter;
      real temp, temp1, temp2;
      int iflag;
      real delta = 0.;
      real ratio;
      real fnorm, gnorm;
      real pnorm, xnorm = 0., fnorm1, actred, dirder, epsmch, prered;
      int info;

      epsmch = std::numeric_limits< real >::epsilon( );

      info = 0;
      iflag = 0;
      nfev = 0;

      /*     check the input parameters for errors. */

      if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. ||
          gtol < 0. || maxfev <= 0 || factor <= 0.) {
	goto TERMINATE;
      }
      if (mode == 2) {
        for (j = 0; j < n; ++j) {
          if (diag[j] <= 0.) {
            goto TERMINATE;
          }
        }
      }

      /*     evaluate the function at the starting point */
      /*     and calculate its norm. */
      fcn(m, n, x, fvec, iflag, xptr);
      nfev = 1;
      if (iflag < 0) {
	goto TERMINATE;
      }
      fnorm = this->enorm(m, fvec);

      /*     initialize levenberg-marquardt parameter and iteration counter. */

      par = 0.;
      iter = 1;

      /*     beginning of the outer loop. */

      for (;;) {

        /*        calculate the jacobian matrix. */

        iflag = fdjac2(fcn, m, n, x, fvec, fjac, ldfjac, epsfcn, wa4, xptr, high);
        nfev += n;
        if (iflag < 0) {
          goto TERMINATE;
        }

        /*        if requested, call fcn to enable printing of iterates. */

        if (nprint > 0) {
          iflag = 0;
          if ((iter - 1) % nprint == 0) {
            // fcn(m, n, x, fvec, iflag, xptr);
            this->print_progress(m, n, x, fvec);
          }
          if (iflag < 0) {
            goto TERMINATE;
          }
        }

        /*        compute the qr factorization of the jacobian. */

        this->qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

        /*        on the first iteration and if mode is 1, scale according */
        /*        to the norms of the columns of the initial jacobian. */

        if (iter == 1) {
          if (mode != 2) {
            for (j = 0; j < n; ++j) {
              diag[j] = wa2[j];
              if (wa2[j] == 0.) {
                diag[j] = 1.;
              }
            }
          }

          /*        on the first iteration, calculate the norm of the scaled x */
          /*        and initialize the step bound delta. */

          for (j = 0; j < n; ++j) {
            wa3[j] = diag[j] * x[j];
          }
          xnorm = this->enorm(n, wa3);
          delta = factor * xnorm;
          if (delta == 0.) {
            delta = factor;
          }
        }

        /*        form (q transpose)*fvec and store the first n components in */
        /*        qtf. */

        for (i = 0; i < m; ++i) {
          wa4[i] = fvec[i];
        }
        for (j = 0; j < n; ++j) {
          if (fjac[j + j * ldfjac] != 0.) {
            sum = 0.;
            for (i = j; i < m; ++i) {
              sum += fjac[i + j * ldfjac] * wa4[i];
            }
            temp = -sum / fjac[j + j * ldfjac];
            for (i = j; i < m; ++i) {
              wa4[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          fjac[j + j * ldfjac] = wa1[j];
          qtf[j] = wa4[j];
        }

        /*        compute the norm of the scaled gradient. */

        gnorm = 0.;
        if (fnorm != 0.) {
          for (j = 0; j < n; ++j) {
            l = ipvt[j]-1;
            if (wa2[l] != 0.) {
              sum = 0.;
              for (i = 0; i <= j; ++i) {
                sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
              }
              /* Computing MAX */
              d1 = fabs(sum / wa2[l]);
              gnorm = std::max(gnorm,d1);
            }
          }
        }

        /*        test for convergence of the gradient norm. */

        if (gnorm <= gtol) {
          info = 4;
        }
        if (info != 0) {
          goto TERMINATE;
        }

        /*        rescale if necessary. */

        if (mode != 2) {
          for (j = 0; j < n; ++j) {
            /* Computing MAX */
            d1 = diag[j], d2 = wa2[j];
            diag[j] = std::max(d1,d2);
          }
        }

        /*        beginning of the inner loop. */

        do {

          /*           determine the levenberg-marquardt parameter. */

          this->lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par, wa1, wa2, wa3,
                wa4);

          /*           store the direction p and x + p. calculate the norm of p. */

          for (j = 0; j < n; ++j) {
            wa1[j] = -wa1[j];
            wa2[j] = x[j] + wa1[j];
            wa3[j] = diag[j] * wa1[j];
          }
          pnorm = this->enorm(n, wa3);

          /*           on the first iteration, adjust the initial step bound. */

          // dtn
          // If any of the parameter, wa2,
          // is outside the open interval deal with it
          for ( j = 0; j < n; ++j)
            wa2[ j ] = std::max( low[ j ], std::min( wa2[ j ], high[ j ] ) );
          // dtn

          if (iter == 1) {
            delta = std::min(delta,pnorm);
          }

          /*           evaluate the function at x + p and calculate its norm. */

          fcn(m, n, wa2, wa4, iflag, xptr);
          ++(nfev);
          if (iflag < 0) {
            goto TERMINATE;
          }
          fnorm1 = this->enorm(m, wa4);

          /*           compute the scaled actual reduction. */

          actred = -1.;
          if (p1 * fnorm1 < fnorm) {
            /* Computing 2nd power */
            d1 = fnorm1 / fnorm;
            actred = 1. - d1 * d1;
          }

          /*           compute the scaled predicted reduction and */
          /*           the scaled directional derivative. */

          for (j = 0; j < n; ++j) {
            wa3[j] = 0.;
            l = ipvt[j]-1;
            temp = wa1[l];
            for (i = 0; i <= j; ++i) {
              wa3[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          temp1 = this->enorm(n, wa3) / fnorm;
          temp2 = (sqrt(par) * pnorm) / fnorm;
          prered = temp1 * temp1 + temp2 * temp2 / p5;
          dirder = -(temp1 * temp1 + temp2 * temp2);

          /*           compute the ratio of the actual to the predicted */
          /*           reduction. */

          ratio = 0.;
          if (prered != 0.) {
            ratio = actred / prered;
          }

          /*           update the step bound. */

          if (ratio <= p25) {
            if (actred >= 0.) {
              temp = p5;
            } else {
              temp = p5 * dirder / (dirder + p5 * actred);
            }
            if (p1 * fnorm1 >= fnorm || temp < p1) {
              temp = p1;
            }
            /* Computing MIN */
            d1 = pnorm / p1;
            delta = temp * std::min(delta,d1);
            par /= temp;
          } else {
            if (par == 0. || ratio >= p75) {
              delta = pnorm / p5;
              par = p5 * par;
            }
          }

          /*           test for successful iteration. */

          if (ratio >= p0001) {

            /*           successful iteration. update x, fvec, and their norms. */

            for (j = 0; j < n; ++j) {
              x[j] = wa2[j];
              wa2[j] = diag[j] * x[j];
            }
            for (i = 0; i < m; ++i) {
              fvec[i] = wa4[i];
            }
            xnorm = this->enorm(n, wa2);
            fnorm = fnorm1;
            ++iter;
          }

          /*           tests for convergence. */

          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
            info = 1;
          }
          if (delta <= xtol * xnorm) {
            info = 2;
          }
          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info == 2) {
            info = 3;
          }
          if (info != 0) {
            goto TERMINATE;
          }

          /*           tests for termination and stringent tolerances. */

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
            goto TERMINATE;
          }

          /*           end of the inner loop. repeat if iteration unsuccessful. */

        } while (ratio < p0001);

        /*        end of the outer loop. */

      }
    TERMINATE:

      /*     termination, either normal or user imposed. */

      if (iflag < 0) {
	info = iflag;
      }
      if (nprint > 0) {
        // fcn(m, n, x, fvec, iflag, xptr);
        this->print_progress(m, n, x, fvec);
      }
      return info;
      //     last card of subroutine lmdif.

    } // lmdif

  }; // class


  template < typename Func, typename Data, typename real >
  class LevMarDer : public LevMar<Func, Data, real> {

  public:

    int operator( )( int n, real ftol, real xtol, real gtol, int maxfev,
                     real factor, int nprint, std::vector<real>& x,
                     int& nfev, int& njev, real& fmin, std::vector<real>& fjac,
                     const sherpa::Bounds<real>& bounds, int& rank ) {

      int m = static_cast<int>( myfvec.size( ) );

      std::vector<real> diag( n ), qtf( n ), wa1( n ), wa2( n ), wa3( n );
      std::vector<real> wa4( m );
      std::vector<int> ipvt( n );

      const int mode = 1;
      const int ldfjac = m;

      Data usrdata = sherpa::Opt<Data, real>::get_usr_data();
      const std::vector<real>& low = bounds.get_lb();
      const std::vector<real>& high = bounds.get_ub();
      int info = lmder_2( usr_fcn, usrdata, m, n, &x[0],
                          &myfvec[0], &fjac[0], ldfjac, ftol, xtol, gtol,
                          maxfev, &diag[0], mode, factor, nprint,
                          nfev, njev, &ipvt[0], &qtf[0], &wa1[0], &wa2[0],
                          &wa3[0], &wa4[0], low, high);
        real fnorm = this->enorm(m, &myfvec[0]);
        this->covar( n, &fjac[ 0 ], ldfjac, &ipvt[0], ftol, &wa1[0] );
        // rank = covar1( m, n, fnorm * fnorm, &fjac[0], ldfjac, &ipvt[0],
        //                ftol, &wa1[0] );
        fmin = std::pow( fnorm, 2.0 );

        return info;

    }


    LevMarDer( Func fcn, Data xdata, int mfct )
      : LevMar<Func, Data, real>( fcn, xdata ), usr_fcn( fcn ),
        myfvec( mfct ) { }

    int fitme( int n, real ftol, real xtol,
               real gtol, int maxfev, real epsfcn,
               real factor, int nprint, std::vector<real>& x,
               int& nfev, int& njev, real& fmin, std::vector<real>& fjac,
               const sherpa::Bounds<real>& bounds, int& rank ) {

      int info = 0;

      try {

        if ( sherpa::Opt<Data, real>::are_pars_outside_limits( n, x,
                                                               bounds ) )
          throw sherpa::OptErr( sherpa::OptErr::OutOfBound );

	int m = static_cast<int>( myfvec.size( ) );

	std::vector<real> diag( n ), qtf( n ), wa1( n ), wa2( n ), wa3( n );
	std::vector<real> wa4( m );
	std::vector<int> ipvt( n );

	const int mode = 1;
	const int ldfjac = m;

        Data usrdata = sherpa::Opt<Data, real>::get_usr_data();
        const std::vector<real>& low = bounds.get_lb();
        const std::vector<real>& high = bounds.get_ub();
	info = lmder( usr_fcn, usrdata, m, n, &x[0], &myfvec[0], &fjac[0],
                      ldfjac, ftol, xtol, gtol, maxfev, &diag[0], mode, factor,
                      nprint, nfev, njev, &ipvt[0], &qtf[0], &wa1[ 0 ],
                      &wa2[0], &wa3[0], &wa4[0], low, high);

        real fnorm = this->enorm(m, &myfvec[0]);
        rank = covar1( m, n, fnorm * fnorm, &fjac[0], ldfjac, &ipvt[0],
                       ftol, &wa1[0] );
        fmin = std::pow( fnorm, 2.0 );

        return info;

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

      // fmin = std::pow( this->enorm( myfvec.size( ), &myfvec[0] ), 2.0 );
      return info;

    }

    // real eval_func( int maxnfev, const sherpa::Bounds<real>& limits,
    //                 int npar, std::vector<real>& par, int& nfev ) {

    //   if ( sherpa::Opt<Data, real>::are_pars_outside_limits( npar, par,
    //                                                          limits ) ) {
    //     return std::numeric_limits< real >::max( );
    //   }

    //   ++nfev;
    //   int ierr=EXIT_SUCCESS;

    //   const int m = static_cast<int>( myfvec.size( ) );

    //   Data usrdata = sherpa::Opt<Data, real>::get_usr_data();
    //   usr_func( m, npar, &par[0], &myfvec[0], &myfjac[0], m, ierr, usrdata );

    //   real fval = pow( this->enorm( m, &myfvec[0] ), 2.0 );
    //   if ( EXIT_SUCCESS != ierr )
    //     throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
    //   if ( nfev >= maxnfev )
    //     throw sherpa::OptErr( sherpa::OptErr::MaxFev );

    //   return fval;

    // }

  private:

    Func usr_fcn;
    std::vector< real > myfvec;

    int covar1(int m, int n, real fsumsq, real *r, int ldr,
               const int *ipvt, real tol, real *wa) {

      /* Local variables */
      int i, j, k, l, ii, jj;
      int sing;
      real temp, tolr;
      tolr = tol * fabs(r[0]);

      /*     form the inverse of r in the full upper triangle of r. */

      l = -1;
      for (k = 0; k < n; ++k) {
	if (fabs(r[k + k * ldr]) <= tolr) {
          break;
	}
	r[k + k * ldr] = 1. / r[k + k * ldr];
	if (k > 0) {
          for (j = 0; j < k; ++j) {
            // coverity[copy_paste_error]
            temp = r[k + k * ldr] * r[j + k * ldr];
            r[j + k * ldr] = 0.;
            for (i = 0; i <= j; ++i) {
              r[i + k * ldr] -= temp * r[i + j * ldr];
            }
          }
        }
	l = k;
      }

      /*     form the full upper triangle of the inverse of (r transpose)*r */
      /*     in the full upper triangle of r. */

      if (l >= 0) {
        for (k = 0; k <= l; ++k) {
          if (k > 0) {
            for (j = 0; j < k; ++j) {
              temp = r[j + k * ldr];
              for (i = 0; i <= j; ++i) {
                r[i + j * ldr] += temp * r[i + k * ldr];
              }
            }
          }
          temp = r[k + k * ldr];
          for (i = 0; i <= k; ++i) {
            r[i + k * ldr] *= temp;
          }
        }
      }

      /*     form the full lower triangle of the covariance matrix */
      /*     in the strict lower triangle of r and in wa. */

      for (j = 0; j < n; ++j) {
	jj = ipvt[j]-1;
	sing = j > l;
	for (i = 0; i <= j; ++i) {
          if (sing) {
            r[i + j * ldr] = 0.;
          }
          ii = ipvt[i]-1;
          if (ii > jj) {
            r[ii + jj * ldr] = r[i + j * ldr];
          }
          else if (ii < jj) {
            r[jj + ii * ldr] = r[i + j * ldr];
          }
	}
	wa[jj] = r[j + j * ldr];
      }

      /*     symmetrize the covariance matrix in r. */

      temp = fsumsq / (m - (l + 1));
      for (j = 0; j < n; ++j) {
	for (i = 0; i < j; ++i) {
          r[j + i * ldr] *= temp;
          r[i + j * ldr] = r[j + i * ldr];
	}
	r[j + j * ldr] = temp * wa[j];
      }

      /*     last card of subroutine covar. */
      if (l == (n - 1)) {
        return 0;
      }
      return l + 1;

    }

    // c     **********
    // c
    // c     subroutine lmder
    // c
    // c     the purpose of lmder is to minimize the sum of the squares of
    // c     m nonlinear functions in n variables by a modification of
    // c     the levenberg-marquardt algorithm. the user must provide a
    // c     subroutine which calculates the functions and the jacobian.
    // c
    // c     the subroutine statement is
    // c
    // c       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
    // c                        maxfev,diag,mode,factor,nprint,info,nfev,
    // c                        njev,ipvt,qtf,wa1,wa2,wa3,wa4)
    // c
    // c     where
    // c
    // c       fcn is the name of the user-supplied subroutine which
    // c         calculates the functions and the jacobian. fcn must
    // c         be declared in an external statement in the user
    // c         calling program, and should be written as follows.
    // c
    // c         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
    // c         integer m,n,ldfjac,iflag
    // c         double precision x(n),fvec(m),fjac(ldfjac,n)
    // c         ----------
    // c         if iflag = 1 calculate the functions at x and
    // c         return this vector in fvec. do not alter fjac.
    // c         if iflag = 2 calculate the jacobian at x and
    // c         return this matrix in fjac. do not alter fvec.
    // c         ----------
    // c         return
    // c         end
    // c
    // c         the value of iflag should not be changed by fcn unless
    // c         the user wants to terminate execution of lmder.
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
    // c         occurs when the number of calls to fcn with iflag = 1
    // c         has reached maxfev.
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
    // c         interval (.1,100.).100. is a generally recommended value.
    // c
    // c       nprint is an integer input variable that enables controlled
    // c         printing of iterates if it is positive. in this case,
    // c         fcn is called with iflag = 0 at the beginning of the first
    // c         iteration and every nprint iterations thereafter and
    // c         immediately prior to return, with x, fvec, and fjac
    // c         available for printing. fvec and fjac should not be
    // c         altered. if nprint is not positive, no special calls
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
    // c         info = 5  number of calls to fcn with iflag = 1 has
    // c                   reached maxfev.
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
    // c         calls to fcn with iflag = 1.
    // c
    // c       njev is an integer output variable set to the number of
    // c         calls to fcn with iflag = 2.
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
    // c       minpack-supplied ... dpmpar,enorm,lmpar,qrfac
    // c
    // c       fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod
    // c
    // c     argonne national laboratory. minpack project. march 1980.
    // c     burton s. garbow, kenneth e. hillstrom, jorge j. more
    // c
    // c     **********

    int lmder( Func fcn, Data xptr, int m, int n, real *x,
               real *fvec, real *fjac, int ldfjac, real ftol,
               real xtol, real gtol, int maxfev, real *diag,
               int mode, real factor, int nprint,
               int& nfev, int& njev, int *ipvt, real *qtf,
               real *wa1, real *wa2, real *wa3, real *wa4,
               const std::vector<real>& low,
	       const std::vector<real>& high ) {

      // Initialized data

      const real p1=0.1;
      const real p5=0.5;
      const real p25=0.25;
      const real p75=0.75;
      const real p0001=1.0e-4;
      /* System generated locals */
      real d1, d2;

      /* Local variables */
      int i, j, l;
      real par, sum;
      int iter;
      real temp, temp1, temp2;
      int iflag;
      real delta = 0.;
      real ratio;
      real fnorm, gnorm, pnorm, xnorm = 0., fnorm1, actred, dirder,
        epsmch, prered;
      int info;


      /*     epsmch is the machine precision. */

      epsmch = std::numeric_limits< real >::epsilon( );

      info = 0;
      iflag = 0;
      nfev = 0;
      njev = 0;

      /*     check the input parameters for errors. */

      if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. ||
          gtol < 0. || maxfev <= 0 || factor <= 0.) {
	goto TERMINATE;
      }
      if (mode == 2) {
        for (j = 0; j < n; ++j) {
          if (diag[j] <= 0.) {
            goto TERMINATE;
          }
        }
      }

      /*     evaluate the function at the starting point */
      /*     and calculate its norm. */

      iflag = fcn(m, n, x, fvec, fjac, ldfjac, 1, xptr);
      nfev = 1;
      if (iflag < 0) {
	goto TERMINATE;
      }
      fnorm = this->enorm(m, fvec);

      /*     initialize levenberg-marquardt parameter and iteration counter. */

      par = 0.;
      iter = 1;

      /*     beginning of the outer loop. */

      for (;;) {

        /*        calculate the jacobian matrix. */

        iflag = fcn(m, n, x, fvec, fjac, ldfjac, 2, xptr);
        ++njev;
        if (iflag < 0) {
          goto TERMINATE;
        }

        /*        if requested, call fcn to enable printing of iterates. */

        if (nprint > 0) {
          iflag = 0;
          if ((iter - 1) % nprint == 0) {
            // iflag = fcn(m, n, x, fvec, fjac, ldfjac, 0, xptr);
            this->print_progress(m, n, x, fvec);
          }
          if (iflag < 0) {
            goto TERMINATE;
          }
        }

        /*        compute the qr factorization of the jacobian. */

        this->qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

        /*        on the first iteration and if mode is 1, scale according */
        /*        to the norms of the columns of the initial jacobian. */

        if (iter == 1) {
          if (mode != 2) {
            for (j = 0; j < n; ++j) {
              diag[j] = wa2[j];
              if (wa2[j] == 0.) {
                diag[j] = 1.;
              }
            }
          }

          /*        on the first iteration, calculate the norm of the scaled x */
          /*        and initialize the step bound delta. */

          for (j = 0; j < n; ++j) {
            wa3[j] = diag[j] * x[j];
          }
          xnorm = this->enorm(n, wa3);
          delta = factor * xnorm;
          if (delta == 0.) {
            delta = factor;
          }
        }

        /*        form (q transpose)*fvec and store the first n components in */
        /*        qtf. */

        for (i = 0; i < m; ++i) {
          wa4[i] = fvec[i];
        }
        for (j = 0; j < n; ++j) {
          if (fjac[j + j * ldfjac] != 0.) {
            sum = 0.;
            for (i = j; i < m; ++i) {
              sum += fjac[i + j * ldfjac] * wa4[i];
            }
            temp = -sum / fjac[j + j * ldfjac];
            for (i = j; i < m; ++i) {
              wa4[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          fjac[j + j * ldfjac] = wa1[j];
          qtf[j] = wa4[j];
        }

        /*        compute the norm of the scaled gradient. */

        gnorm = 0.;
        if (fnorm != 0.) {
          for (j = 0; j < n; ++j) {
            l = ipvt[j]-1;
            if (wa2[l] != 0.) {
              sum = 0.;
              for (i = 0; i <= j; ++i) {
                sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
              }
              /* Computing MAX */
              d1 = fabs(sum / wa2[l]);
              gnorm = std::max(gnorm,d1);
            }
          }
        }

        /*        test for convergence of the gradient norm. */

        if (gnorm <= gtol) {
          info = 4;
        }
        if (info != 0) {
          goto TERMINATE;
        }

        /*        rescale if necessary. */

        if (mode != 2) {
          for (j = 0; j < n; ++j) {
            /* Computing MAX */
            d1 = diag[j], d2 = wa2[j];
            diag[j] = std::max(d1,d2);
          }
        }

        /*        beginning of the inner loop. */

        do {

          /*           determine the levenberg-marquardt parameter. */

          this->lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta,
                &par, wa1, wa2, wa3, wa4);

          /*           store the direction p and x + p. calculate the norm of p. */

          for (j = 0; j < n; ++j) {
            wa1[j] = -wa1[j];
            wa2[j] = x[j] + wa1[j];
            wa3[j] = diag[j] * wa1[j];
          }
          pnorm = this->enorm(n, wa3);

          /*           on the first iteration, adjust the initial step bound. */

          // dtn
          // If any of the parameter, wa2,
          // is outside the open interval deal with it
          for ( j = 0; j < n; ++j)
            wa2[ j ] = std::max( low[ j ], std::min( wa2[ j ], high[ j ] ) );
          // dtn

          if (iter == 1) {
            delta = std::min(delta,pnorm);
          }

          /*           evaluate the function at x + p and calculate its norm. */

          iflag = fcn(m, n, wa2, wa4, fjac, ldfjac, 1, xptr);
          ++nfev;
          if (iflag < 0) {
            goto TERMINATE;
          }
          fnorm1 = this->enorm(m, wa4);

          /*           compute the scaled actual reduction. */

          actred = -1.;
          if (p1 * fnorm1 < fnorm) {
            /* Computing 2nd power */
            d1 = fnorm1 / fnorm;
            actred = 1. - d1 * d1;
          }

          /*           compute the scaled predicted reduction and */
          /*           the scaled directional derivative. */

          for (j = 0; j < n; ++j) {
            wa3[j] = 0.;
            l = ipvt[j]-1;
            temp = wa1[l];
            for (i = 0; i <= j; ++i) {
              wa3[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          temp1 = this->enorm(n, wa3) / fnorm;
          temp2 = (sqrt(par) * pnorm) / fnorm;
          prered = temp1 * temp1 + temp2 * temp2 / p5;
          dirder = -(temp1 * temp1 + temp2 * temp2);

          /*           compute the ratio of the actual to the predicted */
          /*           reduction. */

          ratio = 0.;
          if (prered != 0.) {
            ratio = actred / prered;
          }

          /*           update the step bound. */

          if (ratio <= p25) {
            if (actred >= 0.) {
              temp = p5;
            } else {
              temp = p5 * dirder / (dirder + p5 * actred);
            }
            if (p1 * fnorm1 >= fnorm || temp < p1) {
              temp = p1;
            }
            /* Computing MIN */
            d1 = pnorm / p1;
            delta = temp * std::min(delta,d1);
            par /= temp;
          } else {
            if (par == 0. || ratio >= p75) {
              delta = pnorm / p5;
              par = p5 * par;
            }
          }

          /*           test for successful iteration. */

          if (ratio >= p0001) {

            /*           successful iteration. update x, fvec, and their norms. */

            for (j = 0; j < n; ++j) {
              x[j] = wa2[j];
              wa2[j] = diag[j] * x[j];
            }
            for (i = 0; i < m; ++i) {
              fvec[i] = wa4[i];
            }
            xnorm = this->enorm(n, wa2);
            fnorm = fnorm1;
            ++iter;
          }

          /*           tests for convergence. */

          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
            info = 1;
          }
          if (delta <= xtol * xnorm) {
            info = 2;
          }
          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info == 2) {
            info = 3;
          }
          if (info != 0) {
            goto TERMINATE;
          }

          /*           tests for termination and stringent tolerances. */

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
            goto TERMINATE;
          }

          /*           end of the inner loop. repeat if iteration unsuccessful. */

        } while (ratio < p0001);

        /*        end of the outer loop. */

      }
    TERMINATE:

      /*     termination, either normal or user imposed. */

      if (iflag < 0) {
	info = iflag;
      }
      if (nprint > 0) {
	// fcn(m, n, x, fvec, fjac, ldfjac, 0, xptr);
        this->print_progress(m, n, x, fvec);
      }
      return info;

      /*     last card of subroutine lmder. */

    }

    int lmder_2( Func fcn, Data xptr, int m, int n, real *x, real *fvec,
                 real *fjac, int ldfjac, real ftol, real xtol, real gtol,
                 int maxfev, real *diag, int mode, real factor, int nprint,
                 int& nfev, int& njev, int *ipvt, real *qtf,
                 real *wa1, real *wa2, real *wa3, real *wa4,
                 const std::vector<real>& low,
                 const std::vector<real>& high ) {

      // Initialized data

      const real p1=0.1;
      const real p5=0.5;
      const real p25=0.25;
      const real p75=0.75;
      const real p0001=1.0e-4;
      /* System generated locals */
      real d1, d2;

      /* Local variables */
      int i, j, l;
      real par, sum;
      int iter;
      real temp, temp1, temp2;
      int iflag;
      real delta = 0.;
      real ratio;
      real fnorm, gnorm, pnorm, xnorm = 0., fnorm1, actred, dirder,
        epsmch, prered;
      int info;


      /*     epsmch is the machine precision. */

      epsmch = std::numeric_limits< real >::epsilon( );

      info = 0;
      iflag = 0;
      nfev = 0;
      njev = 0;

      /*     check the input parameters for errors. */

      if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. ||
          gtol < 0. || maxfev <= 0 || factor <= 0.) {
	goto TERMINATE;
      }
      if (mode == 2) {
        for (j = 0; j < n; ++j) {
          if (diag[j] <= 0.) {
            goto TERMINATE;
          }
        }
      }

      /*     evaluate the function at the starting point */
      /*     and calculate its norm. */

      // iflag = fcn(m, n, x, fvec, fjac, ldfjac, 1, xptr);
      iflag = 1;
      fcn(m, n, x, fvec, iflag, xptr);
      nfev = 1;
      if (iflag < 0) {
	goto TERMINATE;
      }
      fnorm = this->enorm(m, fvec);

      /*     initialize levenberg-marquardt parameter and iteration counter. */

      par = 0.;
      iter = 1;

      /*     beginning of the outer loop. */

      for (;;) {

        /*        calculate the jacobian matrix. */

        // iflag = fcn(m, n, x, fvec, fjac, ldfjac, 2, xptr);
        iflag = 2;
        fcn(m, n, x, fjac, iflag, xptr);
        ++njev;
        if (iflag < 0) {
          goto TERMINATE;
        }

        /*        if requested, call fcn to enable printing of iterates. */

        if (nprint > 0) {
          iflag = 0;
          if ((iter - 1) % nprint == 0) {
            // iflag = fcn(m, n, x, fvec, fjac, ldfjac, 0, xptr);
            this->print_progress(m, n, x, fvec);
          }
          if (iflag < 0) {
            goto TERMINATE;
          }
        }

        /*        compute the qr factorization of the jacobian. */

        this->qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

        /*        on the first iteration and if mode is 1, scale according */
        /*        to the norms of the columns of the initial jacobian. */

        if (iter == 1) {
          if (mode != 2) {
            for (j = 0; j < n; ++j) {
              diag[j] = wa2[j];
              if (wa2[j] == 0.) {
                diag[j] = 1.;
              }
            }
          }

          /*        on the first iteration, calculate the norm of the scaled x */
          /*        and initialize the step bound delta. */

          for (j = 0; j < n; ++j) {
            wa3[j] = diag[j] * x[j];
          }
          xnorm = this->enorm(n, wa3);
          delta = factor * xnorm;
          if (delta == 0.) {
            delta = factor;
          }
        }

        /*        form (q transpose)*fvec and store the first n components in */
        /*        qtf. */

        for (i = 0; i < m; ++i) {
          wa4[i] = fvec[i];
        }
        for (j = 0; j < n; ++j) {
          if (fjac[j + j * ldfjac] != 0.) {
            sum = 0.;
            for (i = j; i < m; ++i) {
              sum += fjac[i + j * ldfjac] * wa4[i];
            }
            temp = -sum / fjac[j + j * ldfjac];
            for (i = j; i < m; ++i) {
              wa4[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          fjac[j + j * ldfjac] = wa1[j];
          qtf[j] = wa4[j];
        }

        /*        compute the norm of the scaled gradient. */

        gnorm = 0.;
        if (fnorm != 0.) {
          for (j = 0; j < n; ++j) {
            l = ipvt[j]-1;
            if (wa2[l] != 0.) {
              sum = 0.;
              for (i = 0; i <= j; ++i) {
                sum += fjac[i + j * ldfjac] * (qtf[i] / fnorm);
              }
              /* Computing MAX */
              d1 = fabs(sum / wa2[l]);
              gnorm = std::max(gnorm,d1);
            }
          }
        }

        /*        test for convergence of the gradient norm. */

        if (gnorm <= gtol) {
          info = 4;
        }
        if (info != 0) {
          goto TERMINATE;
        }

        /*        rescale if necessary. */

        if (mode != 2) {
          for (j = 0; j < n; ++j) {
            /* Computing MAX */
            d1 = diag[j], d2 = wa2[j];
            diag[j] = std::max(d1,d2);
          }
        }

        /*        beginning of the inner loop. */

        do {

          /*           determine the levenberg-marquardt parameter. */

          this->lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta,
                &par, wa1, wa2, wa3, wa4);

          /*           store the direction p and x + p. calculate the norm of p. */

          for (j = 0; j < n; ++j) {
            wa1[j] = -wa1[j];
            wa2[j] = x[j] + wa1[j];
            wa3[j] = diag[j] * wa1[j];
          }
          pnorm = this->enorm(n, wa3);

          /*           on the first iteration, adjust the initial step bound. */

          // dtn
          // If any of the parameter, wa2,
          // is outside the open interval deal with it
          for ( j = 0; j < n; ++j)
            wa2[ j ] = std::max( low[ j ], std::min( wa2[ j ], high[ j ] ) );
          // dtn

          if (iter == 1) {
            delta = std::min(delta,pnorm);
          }

          /*           evaluate the function at x + p and calculate its norm. */

          // iflag = fcn(m, n, wa2, wa4, fjac, ldfjac, 1, xptr);
          iflag = 1;
          fcn(m, n, wa2, wa4, iflag, xptr);
          ++nfev;
          if (iflag < 0) {
            goto TERMINATE;
          }
          fnorm1 = this->enorm(m, wa4);

          /*           compute the scaled actual reduction. */

          actred = -1.;
          if (p1 * fnorm1 < fnorm) {
            /* Computing 2nd power */
            d1 = fnorm1 / fnorm;
            actred = 1. - d1 * d1;
          }

          /*           compute the scaled predicted reduction and */
          /*           the scaled directional derivative. */

          for (j = 0; j < n; ++j) {
            wa3[j] = 0.;
            l = ipvt[j]-1;
            temp = wa1[l];
            for (i = 0; i <= j; ++i) {
              wa3[i] += fjac[i + j * ldfjac] * temp;
            }
          }
          temp1 = this->enorm(n, wa3) / fnorm;
          temp2 = (sqrt(par) * pnorm) / fnorm;
          prered = temp1 * temp1 + temp2 * temp2 / p5;
          dirder = -(temp1 * temp1 + temp2 * temp2);

          /*           compute the ratio of the actual to the predicted */
          /*           reduction. */

          ratio = 0.;
          if (prered != 0.) {
            ratio = actred / prered;
          }

          /*           update the step bound. */

          if (ratio <= p25) {
            if (actred >= 0.) {
              temp = p5;
            } else {
              temp = p5 * dirder / (dirder + p5 * actred);
            }
            if (p1 * fnorm1 >= fnorm || temp < p1) {
              temp = p1;
            }
            /* Computing MIN */
            d1 = pnorm / p1;
            delta = temp * std::min(delta,d1);
            par /= temp;
          } else {
            if (par == 0. || ratio >= p75) {
              delta = pnorm / p5;
              par = p5 * par;
            }
          }

          /*           test for successful iteration. */

          if (ratio >= p0001) {

            /*           successful iteration. update x, fvec, and their norms. */

            for (j = 0; j < n; ++j) {
              x[j] = wa2[j];
              wa2[j] = diag[j] * x[j];
            }
            for (i = 0; i < m; ++i) {
              fvec[i] = wa4[i];
            }
            xnorm = this->enorm(n, wa2);
            fnorm = fnorm1;
            ++iter;
          }

          /*           tests for convergence. */

          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) {
            info = 1;
          }
          if (delta <= xtol * xnorm) {
            info = 2;
          }
          if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1. && info == 2) {
            info = 3;
          }
          if (info != 0) {
            goto TERMINATE;
          }

          /*           tests for termination and stringent tolerances. */

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
            goto TERMINATE;
          }

          /*           end of the inner loop. repeat if iteration unsuccessful. */

        } while (ratio < p0001);

        /*        end of the outer loop. */

      }
    TERMINATE:

      /*     termination, either normal or user imposed. */

      if (iflag < 0) {
	info = iflag;
      }
      if (nprint > 0) {
	// fcn(m, n, x, fvec, fjac, ldfjac, 0, xptr);
        this->print_progress(m, n, x, fvec);
      }
      return info;

      /*     last card of subroutine lmder. */

    }

  };

} // namespace

#endif
