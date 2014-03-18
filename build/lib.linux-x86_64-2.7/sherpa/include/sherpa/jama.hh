#ifndef jama_hh
#define jama_hh


//
//     Adapted from the TNT project
//     http://math.nist.gov/tnt/download.html
//
//     This software was developed at the National Institute of Standards
//     and Technology (NIST) by employees of the Federal Government in the
//     course of their official duties. Pursuant to title 17 Section 105
//     of the United States Code this software is not subject to copyright
//     protection and is in the public domain. NIST assumes no responsibility
//     whatsoever for its use by other parties, and makes no guarantees,
//     expressed or implied, about its quality, reliability, or any other
//     characteristic.
//
//     We would appreciate acknowledgement if the software is incorporated in
//     redistributable libraries or applications.


#include <algorithm>
#include <assert.h>
#include <cmath>

#include "myArray.hh"
#include "utils.hh"

namespace jama {

  //
  // result = lhs * rhs
  // 
  template< typename T >
  inline void matrix_multiply( const sherpa::Array2d<T>& lhs,
			       const sherpa::Array2d<T>& rhs,
			       sherpa::Array2d<T>& result ) {


    //
    // col size of lhs must == row size of rhs for multiplication to take place
    //
    const int lhs_col = lhs.ncols( );
    assert( lhs_col == rhs.nrows( ) );

    //
    // the result row size must be correct
    //
    const int lhs_row = lhs.nrows( );
    assert ( result.nrows( ) == lhs_row );

    //
    // the result col size must be correct
    //
    const int rhs_col = rhs.ncols( );
    assert( result.ncols( ) == rhs_col );

    for( int rr = 0;  rr < lhs_row;  rr++  ) {
      for( int cc = 0;  cc < rhs_col;  cc++  ) {
	T result_ij = lhs[ rr ][ 0 ] * rhs[ 0 ][ cc ];
	for( int k = 1;  k < lhs_col;  k++  ) {
	  result_ij += lhs[ rr ][ k ] * rhs[ k ][ cc ];
	}
	result[ rr ][ cc ] = result_ij;
      }
    }

  }

  //
  // result = lhs * rhs
  // 
  template< typename T >
  inline void matrix_multiply( const sherpa::Array2d<T>& lhs,
			       const std::vector<T>& rhs,
			       sherpa::Array2d<T>& result ) {

    const int rhs_nrow = rhs.size( );
    assert( lhs.ncols( ) == rhs_nrow );
    for ( int ii = 0; ii < rhs_nrow; ++ii )
      for ( int jj = 0; jj < rhs_nrow; ++jj )
	result[ ii ][ jj ] = lhs[ ii ][ jj ] * rhs[ jj ];

  }

  template< typename T >
  inline void transpose( sherpa::Array2d<T>& arg ) {

    const int nrows = arg.nrows( );
    assert( nrows == arg.ncols( ) );
    for( int ii = 0; ii < nrows; ++ii )
      for ( int jj = 0; jj < ii; ++jj )
	std::swap( arg[ ii ][ jj ], arg[ jj ][ ii ] );

  }

  template< typename Real >
  Real hypot( const Real& a, const Real& b ) {
    if ( 0 == a ) {
      return std::fabs( b );
    } else {
      Real c = b / a;
      return std::fabs( a ) * std::sqrt( 1.0 + c * c );
    }
  }


  /** Singular Value Decomposition.
      <P>
      For an m-by-n matrix A with m >= n, the singular value decomposition is
      an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
      an n-by-n orthogonal matrix V so that A = U*S*V'.
      <P>
      The singular values, sigma[k] = S[k][k], are ordered so that
      sigma[0] >= sigma[1] >= ... >= sigma[n-1].
      <P>
      The singular value decompostion always exists, so the constructor will
      never fail.  The matrix condition number and the effective numerical
      rank can be computed from this decomposition.

      <p>
      (Adapted from JAMA, a Java Matrix Library, developed by jointly 
      by the Mathworks and NIST; see http://math.nist.gov/javanumerics/jama).
  */

  template< typename Real >
  class SVD {

  private:
    const int m, n;
    std::vector< Real > s; // one-dimensional array of singular values
    sherpa::Array2d< Real > U;
    sherpa::Array2d< Real > V; // the right singular vectors

  public:

    SVD( const sherpa::Array2d< Real >& Arg ) :
      m( Arg.nrows( ) ), n( Arg.ncols( ) ), s( std::min( m + 1, n ) ),
      U( m, std::min( m, n ) ), V( n, n ) {

      int nu = std::min(m,n);

      // dtn
      //       m = Arg.dim1();
      //       n = Arg.dim2();
      // int nu = std::min(m,n);
      //       s = std::vector< Real >(std::min(m+1,n)); 
      //       U = sherpa::Array2d< Real >(m, nu, Real(0));
      //       V = sherpa::Array2d< Real >(n,n);
      std::vector< Real > e(n);
      std::vector< Real > work(m);
      //      sherpa::Array2d< Real > A(Arg.copy());
      sherpa::Array2d< Real > A( m, n );
      for ( int ii = 0; ii < m; ++ii )
	for ( int jj = 0; jj < n; ++jj )
	  A[ ii ][ jj ] = Arg[ ii ][ jj ];
      // dtn
      int wantu = 1;  					/* boolean */
      int wantv = 1;  					/* boolean */
      int i=0, j=0, k=0;

      // Reduce A to bidiagonal form, storing the diagonal elements
      // in s and the super-diagonal elements in e.

      int nct = std::min(m-1,n);
      int nrt = std::max(0,std::min(n-2,m));
      for (k = 0; k < std::max(nct,nrt); k++) {
	if (k < nct) {

	  // Compute the transformation for the k-th column and
	  // place the k-th diagonal in s[k].
	  // Compute 2-norm of k-th column without under/overflow.
	  s[k] = 0;
	  for (i = k; i < m; i++) {
	    s[k] = hypot(s[k],A[i][k]);
	  }
	  if (s[k] != 0.0) {
	    if (A[k][k] < 0.0) {
	      s[k] = -s[k];
	    }
	    for (i = k; i < m; i++) {
	      A[i][k] /= s[k];
	    }
	    A[k][k] += 1.0;
	  }
	  s[k] = -s[k];
	}
	for (j = k+1; j < n; j++) {
	  if ((k < nct) && (s[k] != 0.0))  {

	    // Apply the transformation.

	    Real t(0.0);
	    for (i = k; i < m; i++) {
	      t += A[i][k]*A[i][j];
	    }
	    t = -t/A[k][k];
	    for (i = k; i < m; i++) {
	      A[i][j] += t*A[i][k];
	    }
	  }

	  // Place the k-th row of A into e for the
	  // subsequent calculation of the row transformation.

	  e[j] = A[k][j];
	}
	if (wantu & (k < nct)) {

	  // Place the transformation in U for subsequent back
	  // multiplication.

	  for (i = k; i < m; i++) {
	    U[i][k] = A[i][k];
	  }
	}
	if (k < nrt) {

	  // Compute the k-th row transformation and place the
	  // k-th super-diagonal in e[k].
	  // Compute 2-norm without under/overflow.
	  e[k] = 0;
	  for (i = k+1; i < n; i++) {
	    e[k] = hypot(e[k],e[i]);
	  }
	  if (e[k] != 0.0) {
	    if (e[k+1] < 0.0) {
	      e[k] = -e[k];
	    }
	    for (i = k+1; i < n; i++) {
	      e[i] /= e[k];
	    }
	    e[k+1] += 1.0;
	  }
	  e[k] = -e[k];
	  if ((k+1 < m) & (e[k] != 0.0)) {

	    // Apply the transformation.

	    for (i = k+1; i < m; i++) {
	      work[i] = 0.0;
	    }
	    for (j = k+1; j < n; j++) {
	      for (i = k+1; i < m; i++) {
		work[i] += e[j]*A[i][j];
	      }
	    }
	    for (j = k+1; j < n; j++) {
	      Real t(-e[j]/e[k+1]);
	      for (i = k+1; i < m; i++) {
		A[i][j] += t*work[i];
	      }
	    }
	  }
	  if (wantv) {

	    // Place the transformation in V for subsequent
	    // back multiplication.

	    for (i = k+1; i < n; i++) {
	      V[i][k] = e[i];
	    }
	  }
	}
      }

      // Set up the final bidiagonal matrix or order p.

      int p = std::min(n,m+1);
      if (nct < n) {
	s[nct] = A[nct][nct];
      }
      if (m < p) {
	s[p-1] = 0.0;
      }
      if (nrt+1 < p) {
	e[nrt] = A[nrt][p-1];
      }
      e[p-1] = 0.0;

      // If required, generate U.

      if (wantu) {
	for (j = nct; j < nu; j++) {
	  for (i = 0; i < m; i++) {
	    U[i][j] = 0.0;
	  }
	  U[j][j] = 1.0;
	}
	for (k = nct-1; k >= 0; k--) {
	  if (s[k] != 0.0) {
	    for (j = k+1; j < nu; j++) {
	      Real t(0.0);
	      for (i = k; i < m; i++) {
		t += U[i][k]*U[i][j];
	      }
	      t = -t/U[k][k];
	      for (i = k; i < m; i++) {
		U[i][j] += t*U[i][k];
	      }
	    }
	    for (i = k; i < m; i++ ) {
	      U[i][k] = -U[i][k];
	    }
	    U[k][k] = 1.0 + U[k][k];
	    for (i = 0; i < k-1; i++) {
	      U[i][k] = 0.0;
	    }
	  } else {
	    for (i = 0; i < m; i++) {
	      U[i][k] = 0.0;
	    }
	    U[k][k] = 1.0;
	  }
	}
      }

      // If required, generate V.

      if (wantv) {
	for (k = n-1; k >= 0; k--) {
	  if ((k < nrt) & (e[k] != 0.0)) {
	    for (j = k+1; j < nu; j++) {
	      Real t(0.0);
	      for (i = k+1; i < n; i++) {
		t += V[i][k]*V[i][j];
	      }
	      t = -t/V[k+1][k];
	      for (i = k+1; i < n; i++) {
		V[i][j] += t*V[i][k];
	      }
	    }
	  }
	  for (i = 0; i < n; i++) {
	    V[i][k] = 0.0;
	  }
	  V[k][k] = 1.0;
	}
      }

      // Main iteration loop for the singular values.

      int pp = p-1;
      int iter = 0;
      Real eps(std::pow(2.0,-52.0));
      while (p > 0) {
	int k=0;
	int kase=0;

	// Here is where a test for too many iterations would go.

	// This section of the program inspects for
	// negligible elements in the s and e arrays.  On
	// completion the variables kase and k are set as follows.

	// kase = 1     if s(p) and e[k-1] are negligible and k<p
	// kase = 2     if s(k) is negligible and k<p
	// kase = 3     if e[k-1] is negligible, k<p, and
	//              s(k), ..., s(p) are not negligible (qr step).
	// kase = 4     if e(p-1) is negligible (convergence).

	for (k = p-2; k >= -1; k--) {
	  if (k == -1) {
	    break;
	  }
	  if (std::fabs(e[k]) <= eps*(std::fabs(s[k]) + std::fabs(s[k+1]))) {
	    e[k] = 0.0;
	    break;
	  }
	}
	if (k == p-2) {
	  kase = 4;
	} else {
	  int ks;
	  for (ks = p-1; ks >= k; ks--) {
	    if (ks == k) {
	      break;
	    }
	    Real t( (ks != p ? std::fabs(e[ks]) : 0.) + 
		    (ks != k+1 ? std::fabs(e[ks-1]) : 0.));
	    if (std::fabs(s[ks]) <= eps*t)  {
	      s[ks] = 0.0;
	      break;
	    }
	  }
	  if (ks == k) {
	    kase = 3;
	  } else if (ks == p-1) {
	    kase = 1;
	  } else {
	    kase = 2;
	    k = ks;
	  }
	}
	k++;

	// Perform the task indicated by kase.

	switch (kase) {

	  // Deflate negligible s(p).

	case 1: {
	  Real f(e[p-2]);
	  e[p-2] = 0.0;
	  for (j = p-2; j >= k; j--) {
	    Real t( hypot(s[j],f));
	    Real cs(s[j]/t);
	    Real sn(f/t);
	    s[j] = t;
	    if (j != k) {
	      f = -sn*e[j-1];
	      e[j-1] = cs*e[j-1];
	    }
	    if (wantv) {
	      for (i = 0; i < n; i++) {
		t = cs*V[i][j] + sn*V[i][p-1];
		V[i][p-1] = -sn*V[i][j] + cs*V[i][p-1];
		V[i][j] = t;
	      }
	    }
	  }
	}
	  break;

	  // Split at negligible s(k).

	case 2: {
	  Real f(e[k-1]);
	  e[k-1] = 0.0;
	  for (j = k; j < p; j++) {
	    Real t(hypot(s[j],f));
	    Real cs( s[j]/t);
	    Real sn(f/t);
	    s[j] = t;
	    f = -sn*e[j];
	    e[j] = cs*e[j];
	    if (wantu) {
	      for (i = 0; i < m; i++) {
		t = cs*U[i][j] + sn*U[i][k-1];
		U[i][k-1] = -sn*U[i][j] + cs*U[i][k-1];
		U[i][j] = t;
	      }
	    }
	  }
	}
	  break;

	  // Perform one qr step.

	case 3: {

	  // Calculate the shift.
   
	  Real scale = std::max(std::max(std::max(std::max(std::fabs(s[p-1]),std::fabs(s[p-2])),std::fabs(e[p-2])), std::fabs(s[k])),std::fabs(e[k]));
	  Real sp = s[p-1]/scale;
	  Real spm1 = s[p-2]/scale;
	  Real epm1 = e[p-2]/scale;
	  Real sk = s[k]/scale;
	  Real ek = e[k]/scale;
	  Real b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
	  Real c = (sp*epm1)*(sp*epm1);
	  Real shift = 0.0;
	  if ((b != 0.0) || (c != 0.0)) {
	    shift = sqrt(b*b + c);
	    if (b < 0.0) {
	      shift = -shift;
	    }
	    shift = c/(b + shift);
	  }
	  Real f = (sk + sp)*(sk - sp) + shift;
	  Real g = sk*ek;
   
	  // Chase zeros.
   
	  for (j = k; j < p-1; j++) {
	    Real t = hypot(f,g);
	    Real cs = f/t;
	    Real sn = g/t;
	    if (j != k) {
	      e[j-1] = t;
	    }
	    f = cs*s[j] + sn*e[j];
	    e[j] = cs*e[j] - sn*s[j];
	    g = sn*s[j+1];
	    s[j+1] = cs*s[j+1];
	    if (wantv) {
	      for (i = 0; i < n; i++) {
		t = cs*V[i][j] + sn*V[i][j+1];
		V[i][j+1] = -sn*V[i][j] + cs*V[i][j+1];
		V[i][j] = t;
	      }
	    }
	    t = hypot(f,g);
	    cs = f/t;
	    sn = g/t;
	    s[j] = t;
	    f = cs*e[j] + sn*s[j+1];
	    s[j+1] = -sn*e[j] + cs*s[j+1];
	    g = sn*e[j+1];
	    e[j+1] = cs*e[j+1];
	    if (wantu && (j < m-1)) {
	      for (i = 0; i < m; i++) {
		t = cs*U[i][j] + sn*U[i][j+1];
		U[i][j+1] = -sn*U[i][j] + cs*U[i][j+1];
		U[i][j] = t;
	      }
	    }
	  }
	  e[p-2] = f;
	  iter = iter + 1;
	}
	  break;

	  // Convergence.

	case 4: {

	  // Make the singular values positive.
   
	  if (s[k] <= 0.0) {
	    s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
	    if (wantv) {
	      for (i = 0; i <= pp; i++) {
		V[i][k] = -V[i][k];
	      }
	    }
	  }
   
	  // Order the singular values.
   
	  while (k < pp) {
	    if (s[k] >= s[k+1]) {
	      break;
	    }
	    Real t = s[k];
	    s[k] = s[k+1];
	    s[k+1] = t;
	    if (wantv && (k < n-1)) {
	      for (i = 0; i < n; i++) {
		t = V[i][k+1]; V[i][k+1] = V[i][k]; V[i][k] = t;
	      }
	    }
	    if (wantu && (k < m-1)) {
	      for (i = 0; i < m; i++) {
		t = U[i][k+1]; U[i][k+1] = U[i][k]; U[i][k] = t;
	      }
	    }
	    k++;
	  }
	  iter = 0;
	  p--;
	}
	  break;
	}
      }
    }

    /** Two norm  (max(S)) */

    Real norm2 () {
      return s[0];
    }

    /** Two norm of condition number (max(S)/min(S)) */

    Real cond () {
      return s[0]/s[std::min(m,n)-1];
    }

    /** Effective numerical matrix rank
	@return     Number of nonnegligible singular values.
    */
    int rank( ) {
      Real eps = std::pow(2.0,-52.0);
      Real tol = std::max(m,n)*s[0]*eps;
      int r = 0;
      for (int i = 0; i < n; i++) {
	if (s[i] > tol) {
	  r++;
	}
      }
      return r;
    }

    //
    //   -1              -1
    //  A    = ( U S V' )
    //
    //               -1    -1   -1
    //       == ( V' )    S    U
    //
    //             -1  t
    //       == V S   U
    //
    void inverse( sherpa::Array2d< Real >& inv ) {

      assert( U.nrows( ) == U.ncols( ) );
      assert( U.nrows( ) == V.nrows( ) );
      assert( U.ncols( ) == V.ncols( ) );

      sherpa::Array2d< Real > tmp( n, n );

      for ( int ii = 0; ii < n; ++ii )
	if ( 0.0 != s[ ii ] )
	  s[ ii ] = 1.0 / s[ ii ];

      matrix_multiply( V, s, tmp );
      transpose( U );
      matrix_multiply( tmp, U, inv );

      //
      // restore to its original state? (a bit time consuming).
      //
      transpose( U );
      for ( int ii = 0; ii < n; ++ii )
	if ( 0.0 != s[ ii ] )
	  s[ ii ] = 1.0 / s[ ii ];

    }

  };


  /** LU Decomposition.
      <P>
      For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
      unit lower triangular matrix L, an n-by-n upper triangular matrix U,
      and a permutation vector piv of length m so that A(piv,:) = L*U.
      If m < n, then L is m-by-m and U is m-by-n.
      <P>
      The LU decompostion with pivoting always exists, even if the matrix is
      singular, so the constructor will never fail.  The primary use of the
      LU decomposition is in the solution of square systems of simultaneous
      linear equations.  This will fail if isNonsingular() returns false.
  */

  template< typename Real >
  class LU {

  public:

    /** LU Decomposition
	@param  A   Rectangular matrix
	@return     LU Decomposition object to access L, U and piv.
    */
    LU( const sherpa::Array2d<Real>& A ) : 
      m( A.nrows( ) ), n( A.ncols( ) ), pivsign( 1 ), piv( m ), LU_( m, n ) {

      for ( int ii = 0; ii < m; ++ii )
	for ( int jj = 0; jj <n; ++jj )
	  LU_[ ii ][ jj ] = A[ ii ][ jj ];

      // Use a "left-looking", dot-product, Crout/Doolittle algorithm.


      for (int i = 0; i < m; i++) {
	piv[i] = i;
      }
      pivsign = 1;
      //Real *LUrowi = 0;;
      // const std::vector<Real>& LUrowi = LU_[ 0 ];
      std::vector<Real> LUcolj(m);

      // Outer loop.

      for (int j = 0; j < n; j++) {

	// Make a copy of the j-th column to localize references.

	for (int i = 0; i < m; i++) {
	  LUcolj[i] = LU_[i][j];
	}

	// Apply previous transformations.

	for (int i = 0; i < m; i++) {
	  // LUrowi = LU_[i];

	  // Most of the time is spent in the following dot product.

	  int kmax = std::min(i,j);
	  double s = 0.0;
	  for (int k = 0; k < kmax; k++) {
	    // s += LUrowi[k]*LUcolj[k];
	    s += LU_[i][k]*LUcolj[k];
	  }

	  //LUrowi[j] = LUcolj[i] -= s;
	  LU_[i][j] = LUcolj[i] -= s;
	}
   
	// Find pivot and exchange if necessary.

	int p = j;
	for (int i = j+1; i < m; i++) {
	  if (std::fabs(LUcolj[i]) > std::fabs(LUcolj[p])) {
	    p = i;
	  }
	}
	if (p != j) {
	  int k=0;
	  for (k = 0; k < n; k++) {
	    double t = LU_[p][k]; 
	    LU_[p][k] = LU_[j][k]; 
	    LU_[j][k] = t;
	  }
	  k = piv[p]; 
	  piv[p] = piv[j]; 
	  piv[j] = k;
	  pivsign = -pivsign;
	}

	// Compute multipliers.
         
	if ((j < m) && (LU_[j][j] != 0.0)) {
	  for (int i = j+1; i < m; i++) {
	    LU_[i][j] /= LU_[j][j];
	  }
	}
      }
    }



    /** Is the matrix nonsingular?
	@return     1 (true)  if upper triangular factor U (and hence A) 
	is nonsingular, 0 otherwise.
    */
    int isNonsingular( ) {
      for (int j = 0; j < n; j++) {
	if (LU_[j][j] == 0)
	  return 0;
      }
      return 1;
    }

    /** Compute determinant using LU factors.
	@return     determinant of A, or 0 if A is not square.
    */

    Real det( ) {
      assert( LU_.nrows( ) == LU_.ncols( ) );
      const int n = LU_.nrows( );
      Real d = Real(pivsign);
      for (int j = 0; j < n; j++) {
	d *= LU_[j][j];
      }
      return d;
    }

    /** Solve A*X = B
	@param  B   A Matrix with as many rows as A and any number of columns.
	@return     X so that L*U*X = B(piv,:), if B is nonconformant, returns
	0x0 (null) array.
    */
    void solve( const sherpa::Array2d<Real>& B, sherpa::Array2d<Real>& X ) {

      /* Dimensions: A is mxn, X is nxk, B is mxk */

      const int m = LU_.nrows( );
      const int n = LU_.ncols( );
      assert( m == B.nrows( ) );

      if (!isNonsingular( ))
	throw std::runtime_error( "singular matrix" );

      // Copy right hand side with pivoting
      int nx = B.ncols();


      permute_copy(B, piv, 0, nx-1, X);

      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++) {
	for (int i = k+1; i < n; i++) {
	  for (int j = 0; j < nx; j++) {
	    X[i][j] -= X[k][j]*LU_[i][k];
	  }
	}
      }
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
	for (int j = 0; j < nx; j++) {
	  X[k][j] /= LU_[k][k];
	}
	for (int i = 0; i < k; i++) {
	  for (int j = 0; j < nx; j++) {
	    X[i][j] -= X[k][j]*LU_[i][k];
	  }
	}
      }
      return;
    }


    /** Solve A*x = b, where x and b are vectors of length equal	
	to the number of rows in A.

	@param  b   a vector (std::vector> of length equal to the first dimension
	of A.
	@return x a vector (std::vector> so that L*U*x = b(piv), if B is nonconformant,
	returns 0x0 (null) array.
    */
    void solve( const std::vector<Real>& b, std::vector<Real>& x ) {

      /* Dimensions: A is mxn, X is nxk, B is mxk */
      const int m = LU_.nrows( );
      const int n = LU_.ncols( );

      assert( b.size( ) == (size_t) m );
      if ( !isNonsingular( ) )
	throw std::runtime_error( "singular matrix" );

      permute_copy(b, piv, x);

      // Solve L*Y = B(piv)
      for (int k = 0; k < n; k++) {
	for (int i = k+1; i < n; i++) {
	  x[i] -= x[k]*LU_[i][k];
	}
      }
      
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
	x[k] /= LU_[k][k];
	for (int i = 0; i < k; i++) 
	  x[i] -= x[k]*LU_[i][k];
      }
     
      return;
    }

    void inverse( sherpa::Array2d<Real>& inv ) {
      
      assert( m == n );

      std::vector<Real> x( n );
      for ( int ii = 0; ii < n; ++ii ) {
	std::vector<double> b( n, 0.0 );
	b[ ii ] = 1.0;
	solve( b, x );
	for ( int jj = 0; jj < n; ++jj )
	  inv[ jj ][ ii ] = x[ jj ];
      }

    }

  private:

    const int m, n;
    int pivsign;
    std::vector<int> piv;
    sherpa::Array2d< Real > LU_;
    
    void permute_copy( const sherpa::Array2d<Real>& A,
		       const std::vector<int>& piv, int j0, int j1,
		       sherpa::Array2d<Real>& X ) {
      
      int piv_length = piv.size();
      
      X.resize( piv_length, j1-j0+1 );
      
      for (int i = 0; i < piv_length; i++) 
	for (int j = j0; j <= j1; j++) 
	  X[i][j-j0] = A[piv[i]][j];
    
      return;
    }

    void permute_copy( const std::vector<Real>& A, const std::vector<int>& piv,
		       std::vector<Real>& x ) {
      size_t piv_length = piv.size();
      assert( piv_length == A.size( ) );
    
      x.resize( piv_length );

      for (size_t i = 0; i < piv_length; i++) 
	x[i] = A[piv[i]];
    
      return;
    }

  };

  /** 
      <P>
      For a symmetric, positive definite matrix A, this function
      computes the Cholesky factorization, i.e. it computes a lower 
      triangular matrix L such that A = L*L'.
      If the matrix is not symmetric or positive definite, the function
      computes only a partial decomposition.  This can be tested with
      the is_spd() flag.

      <p>Typical usage looks like:
      <pre>
      Array2D<double> A(n,n);
      Array2D<double> L;

      ... 

      Cholesky<double> chol(A);

      if (chol.is_spd())
      L = chol.getL();
		
      else
      cout << "factorization was not complete.\n";

      </pre>


      <p>
      (Adapted from JAMA, a Java Matrix Library, developed by jointly 
      by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).

  */

  template< typename Real >
  class Cholesky {

  public:

    Cholesky( const sherpa::Array2d<Real>& A ) :
      L_( A.nrows( ), A.ncols( ) ) {

      int m = A.nrows();
      int n = A.ncols();
	
      isspd = (m == n);

      if (m != n) {
	// L_ = sherpa::Array2d<Real>(0,0);
	return;
      }

      // L_ = sherpa::Array2d<Real>(n,n);

      // Main loop.
      for (int j = 0; j < n; j++) {
	Real d(0.0);
	for (int k = 0; k < j; k++) {
	  Real s(0.0);
	  for (int i = 0; i < k; i++) {
	    s += L_[k][i]*L_[j][i];
	  }
	  L_[j][k] = s = (A[j][k] - s)/L_[k][k];
	  d = d + s*s;
	  isspd = isspd && (A[k][j] == A[j][k]); 
	}
	d = A[j][j] - d;
	isspd = isspd && (d > 0.0);
	L_[j][j] = sqrt(d > 0.0 ? d : 0.0);
	for (int k = j+1; k < n; k++) {
	  L_[j][k] = 0.0;
	}
      }
    }

    /**

    Solve a linear system A*x = b, using the previously computed
    cholesky factorization of A: L*L'.

    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     x so that L*L'*x = b.  If b is nonconformat, or if A
    was not symmetric posidtive definite, a null (0x0)
    array is returned.
    */
    void solve(const std::vector<Real> &b, std::vector<Real>& x )
    {
      int n = L_.nrows();
      if (b.dim1() != n)
	return;

      x = b;

      // Solve L*y = b;
      for (int k = 0; k < n; k++) 
	{
	  for (int i = 0; i < k; i++) 
	    x[k] -= x[i]*L_[k][i];
	  x[k] /= L_[k][k];
		
	}

      // Solve L'*X = Y;
      for (int k = n-1; k >= 0; k--) 
	{
	  for (int i = k+1; i < n; i++) 
	    x[k] -= x[i]*L_[i][k];
	  x[k] /= L_[k][k];
	}

      return;
    }


    /**

    Solve a linear system A*X = B, using the previously computed
    cholesky factorization of A: L*L'.

    @param  B   A Matrix with as many rows as A and any number of columns.
    @return     X so that L*L'*X = B.  If B is nonconformat, or if A
    was not symmetric posidtive definite, a null (0x0)
    array is returned.
    */
    void solve( const sherpa::Array2d<Real>& B,
		sherpa::Array2d<Real>& X ) {

      int n = L_.nrows();
      if (B.dim1() != n)
	return;

      for ( int ii = 0; ii < B.nrows( ); ++ii )
	for ( int jj = 0; jj < B.ncols( ); ++jj )
	  X[ ii ][ jj ] = B[ ii ][ jj ];

      int nx = B.ncols();

      // Cleve's original code
#if 0
      // Solve L*Y = B;
      for (int k = 0; k < n; k++) {
	for (int i = k+1; i < n; i++) {
	  for (int j = 0; j < nx; j++) {
	    X[i][j] -= X[k][j]*L_[k][i];
	  }
	}
	for (int j = 0; j < nx; j++) {
	  X[k][j] /= L_[k][k];
	}
      }

      // Solve L'*X = Y;
      for (int k = n-1; k >= 0; k--) {
	for (int j = 0; j < nx; j++) {
	  X[k][j] /= L_[k][k];
	}
	for (int i = 0; i < k; i++) {
	  for (int j = 0; j < nx; j++) {
	    X[i][j] -= X[k][j]*L_[k][i];
	  }
	}
      }
#endif

      // Solve L*y = b;
      for (int j=0; j< nx; j++)
	{
	  for (int k = 0; k < n; k++) 
	    {
	      for (int i = 0; i < k; i++) 
		X[k][j] -= X[i][j]*L_[k][i];
	      X[k][j] /= L_[k][k];
	    }
	}

      // Solve L'*X = Y;
      for (int j=0; j<nx; j++)
	{
	  for (int k = n-1; k >= 0; k--) 
	    {
	      for (int i = k+1; i < n; i++) 
		X[k][j] -= X[i][j]*L_[i][k];
	      X[k][j] /= L_[k][k];
	    }
	}

      return;

    }

  private:

    sherpa::Array2d<Real>& L_;
    int isspd;

  }; 

  /** 
      <p>
      Classical QR Decompisition:
      for an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
      orthogonal matrix Q and an n-by-n upper triangular matrix R so that
      A = Q*R.
      <P>
      The QR decompostion always exists, even if the matrix does not have
      full rank, so the constructor will never fail.  The primary use of the
      QR decomposition is in the least squares solution of nonsquare systems
      of simultaneous linear equations.  This will fail if isFullRank()
      returns 0 (false).

      <p>
      The Q and R factors can be retrived via the getQ() and getR()
      methods. Furthermore, a solve() method is provided to find the
      least squares solution of Ax=b using the QR factors.  

      <p>
      (Adapted from JAMA, a Java Matrix Library, developed by jointly 
      by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
  */
  template< typename Real >
  class QR {

  private:

    /** Array for internal storage of decomposition.
	@serial internal array storage.
    */
    sherpa::Array2d<Real> QR_;
    
    /** Row and column dimensions.
	@serial column dimension.
	@serial row dimension.
    */
    const int m, n;
    
    /** Array for internal storage of diagonal of R.
	@serial diagonal of R.
    */
    std::vector<Real> Rdiag;
    
  public:

    /**
       Create a QR factorization object for A.

       @param A rectangular (m>=n) matrix.
    */
    QR(const sherpa::Array2d<Real> &A) {
      QR_ = A.copy();
      m = A.dim1();
      n = A.dim2();
      Rdiag = std::vector<Real>(n);
      int i=0, j=0, k=0;

      // Main loop.
      for (k = 0; k < n; k++) {
	// Compute 2-norm of k-th column without under/overflow.
	Real nrm = 0;
	for (i = k; i < m; i++) {
	  nrm = hypot(nrm,QR_[i][k]);
	}

	if (nrm != 0.0) {
	  // Form k-th Householder vector.
	  if (QR_[k][k] < 0) {
	    nrm = -nrm;
	  }
	  for (i = k; i < m; i++) {
	    QR_[i][k] /= nrm;
	  }
	  QR_[k][k] += 1.0;

	  // Apply transformation to remaining columns.
	  for (j = k+1; j < n; j++) {
	    Real s = 0.0; 
	    for (i = k; i < m; i++) {
	      s += QR_[i][k]*QR_[i][j];
	    }
	    s = -s/QR_[k][k];
	    for (i = k; i < m; i++) {
	      QR_[i][j] += s*QR_[i][k];
	    }
	  }
	}
	Rdiag[k] = -nrm;
      }
    }


    /**
       Flag to denote the matrix is of full rank.

       @return 1 if matrix is full rank, 0 otherwise.
    */
    int isFullRank() const		
    {
      for (int j = 0; j < n; j++) 
	{
	  if (Rdiag[j] == 0)
            return 0;
	}
      return 1;
    }
	
	


    /** 
   
    Retreive the Householder vectors from QR factorization
    @returns lower trapezoidal matrix whose columns define the reflections
    */

    sherpa::Array2d<Real> getHouseholder (void)  const
    {
      sherpa::Array2d<Real> H(m,n);

      /* note: H is completely filled in by algorithm, so
	 initializaiton of H is not necessary.
      */
      for (int i = 0; i < m; i++) 
	{
	  for (int j = 0; j < n; j++) 
	    {
	      if (i >= j) {
		H[i][j] = QR_[i][j];
	      } else {
		H[i][j] = 0.0;
	      }
	    }
	}
      return H;
    }



    /** Return the upper triangular factor, R, of the QR factorization
	@return     R
    */

    sherpa::Array2d<Real> getR() const
    {
      sherpa::Array2d<Real> R(n,n);
      for (int i = 0; i < n; i++) {
	for (int j = 0; j < n; j++) {
	  if (i < j) {
	    R[i][j] = QR_[i][j];
	  } else if (i == j) {
	    R[i][j] = Rdiag[i];
	  } else {
	    R[i][j] = 0.0;
	  }
	}
      }
      return R;
    }
	
	



    /** 
   	Generate and return the (economy-sized) orthogonal factor
	@param     Q the (ecnomy-sized) orthogonal factor (Q*R=A).
    */

    sherpa::Array2d<Real> getQ() const
    {
      int i=0, j=0, k=0;

      sherpa::Array2d<Real> Q(m,n);
      for (k = n-1; k >= 0; k--) {
	for (i = 0; i < m; i++) {
	  Q[i][k] = 0.0;
	}
	Q[k][k] = 1.0;
	for (j = k; j < n; j++) {
	  if (QR_[k][k] != 0) {
	    Real s = 0.0;
	    for (i = k; i < m; i++) {
	      s += QR_[i][k]*Q[i][j];
	    }
	    s = -s/QR_[k][k];
	    for (i = k; i < m; i++) {
	      Q[i][j] += s*QR_[i][k];
	    }
	  }
	}
      }
      return Q;
    }


    /** Least squares solution of A*x = b
	@param B     m-length array (vector).
	@return x    n-length array (vector) that minimizes the two norm of Q*R*X-B.
	If B is non-conformant, or if QR.isFullRank() is false,
	the routine returns a null (0-length) vector.
    */

    std::vector<Real> solve(const std::vector<Real> &b) const
    {
      if (b.dim1() != m)		/* arrays must be conformant */
	return std::vector<Real>();

      if ( !isFullRank() )		/* matrix is rank deficient */
	{
	  return std::vector<Real>();
	}

      std::vector<Real> x = b.copy();

      // Compute Y = transpose(Q)*b
      for (int k = 0; k < n; k++) 
	{
	  Real s = 0.0; 
	  for (int i = k; i < m; i++) 
	    {
	      s += QR_[i][k]*x[i];
            }
	  s = -s/QR_[k][k];
	  for (int i = k; i < m; i++) 
	    {
	      x[i] += s*QR_[i][k];
            }
	}
      // Solve R*X = Y;
      for (int k = n-1; k >= 0; k--) 
	{
	  x[k] /= Rdiag[k];
	  for (int i = 0; i < k; i++) {
	    x[i] -= x[k]*QR_[i][k];
	  }
	}


      /* return n x nx portion of X */
      std::vector<Real> x_(n);
      for (int i=0; i<n; i++)
	x_[i] = x[i];

      return x_;
    }

    /** Least squares solution of A*X = B
	@param B     m x k Array (must conform).
	@return X     n x k Array that minimizes the two norm of Q*R*X-B. If
	B is non-conformant, or if QR.isFullRank() is false,
	the routine returns a null (0x0) array.
    */

    sherpa::Array2d<Real> solve(const sherpa::Array2d<Real> &B) const
    {
      if (B.dim1() != m)		/* arrays must be conformant */
	return sherpa::Array2d<Real>(0,0);

      if ( !isFullRank() )		/* matrix is rank deficient */
	{
	  return sherpa::Array2d<Real>(0,0);
	}

      int nx = B.dim2(); 
      sherpa::Array2d<Real> X = B.copy();
      int i=0, j=0, k=0;

      // Compute Y = transpose(Q)*B
      for (k = 0; k < n; k++) {
	for (j = 0; j < nx; j++) {
	  Real s = 0.0; 
	  for (i = k; i < m; i++) {
	    s += QR_[i][k]*X[i][j];
	  }
	  s = -s/QR_[k][k];
	  for (i = k; i < m; i++) {
	    X[i][j] += s*QR_[i][k];
	  }
	}
      }
      // Solve R*X = Y;
      for (k = n-1; k >= 0; k--) {
	for (j = 0; j < nx; j++) {
	  X[k][j] /= Rdiag[k];
	}
	for (i = 0; i < k; i++) {
	  for (j = 0; j < nx; j++) {
	    X[i][j] -= X[k][j]*QR_[i][k];
	  }
	}
      }


      /* return n x nx portion of X */
      sherpa::Array2d<Real> X_(n,nx);
      for (i=0; i<n; i++)
	for (j=0; j<nx; j++)
	  X_[i][j] = X[i][j];

      return X_;
    }

  };

  namespace svd {

    //
    // A inv = I
    //
    template< typename T >
    int inverse( const sherpa::Array2d<T>& A, sherpa::Array2d<T>& inv ) {

      SVD<T> mysvd( A );
      mysvd.inverse( inv );
      return EXIT_SUCCESS;

    }

  }                                                            // namespace svd

  namespace lu {

    template< typename T >
    int inverse( const sherpa::Array2d<T>& A, sherpa::Array2d<T>& inv ) {

      LU<T> mylu( A );
      mylu.inverse( inv );
      return 0;

    }

    //
    // A X = B
    //
    template< typename T >
    int solve( const sherpa::Array2d<T>& A, const sherpa::Array2d<T>& B,
	       sherpa::Array2d<T>& X ) {

      LU<T> mylu( A );
      mylu.solve( B, X );

      return EXIT_SUCCESS;
    }

    //
    // A x = b
    //
    template< typename T >
    int solve( const sherpa::Array2d<T>& A, const std::vector<T>& b,
	       std::vector<T>& x ) {
      LU<T> mylu( A );
      mylu.solve( b, x );

      return EXIT_SUCCESS;
    }

  }                                                             // namespace lu

}                                                             // namespace jama

#endif

#ifdef testJama

#include "functor.hh"
#include "StopWatch.hh"


typedef int (*invfuncproto)( const sherpa::Array2d< double >& arg, sherpa::Array2d< double >& arg_inv );
void tst_inverse( invfuncproto func, int numiter, int n, double tol,
		  bool printit, const char* msg ) {

  StopWatch stopwatch( msg );

  srand( numiter * n + 1 );
  sherpa::Array2d< double > a( n, n ), a_inv( n, n ), result( n, n );

  for ( int kk = 0; kk < numiter; ++kk ) {

    for ( int ii = 0; ii < n; ++ii )
      for ( int jj = 0; jj < n; ++jj )
	a[ ii ][ jj ] = rand( ) % 4097;

    func( a, a_inv );
    if ( printit ) {
      jama::matrix_multiply( a, a_inv, result );
      std::cout << "result =\n" << result << '\n';
    }
    func( a_inv, result );
    if ( printit )
      std::cout << "a = \n" << a << "\na_inv_inv = \n" << result<< "\n\n";

    // should be same as the original matrix
    for ( int ii = 0; ii < n; ++ii )
      for ( int jj = 0; jj < ii; ++jj ) {
	bool tmp = sherpa::utils::Knuth_close( a[ ii ][ jj ], result[ ii ][ jj ],
					     tol );
	if ( !tmp ) {
	  fprintf( stdout, "error at ii=%d, jj=%d: %.14e\tvs\t%.14e\n",
		   ii, jj, a[ ii ][ jj ], result[ ii ][ jj ] );
	  assert( tmp );
	}
      }

  }

}

typedef int (*solvefuncproto)( const sherpa::Array2d<double>& A, const sherpa::Array2d<double>& B, sherpa::Array2d<double>& X );
void tst_solve( solvefuncproto func, int numiter, int m, 
		double tol, bool printit, const char* msg ) {

  StopWatch stopwatch( msg );

  const int n = m, k = m;

  srand( numiter * n + 1 );
  sherpa::Array2d< double > A( m, n ), X( n, k ), B( m, k ), result( m, k );

  try {

    for ( int kk = 0; kk < numiter; ++kk ) {

      for ( int ii = 0; ii < m; ++ii ) {
	for ( int jj = 0; jj < n; ++jj )
	  A[ ii ][ jj ] = rand( ) % 4097;
	for ( int kk = 0; kk < k; ++kk )
	  B[ ii ][ kk ] = rand( ) % 4097;
      }

      jama::lu::solve( A, B, X );

      jama::matrix_multiply( A, X, result );

      for ( int ii = 0; ii < m; ++ii )
	for ( int jj = 0; jj < k; ++jj ) {
	  bool tmp = sherpa::utils::Knuth_close( B[ ii ][ jj ], result[ ii ][ jj ], tol );
	  if ( !tmp ) {
	    fprintf( stdout, "error at ii=%d, jj=%d: %.14e\tvs\t%.14e\n",
		     ii, jj, B[ ii ][ jj ], result[ ii ][ jj ] );
	    assert( tmp );
	  }
	}

    }

  } catch( std::runtime_error& re ) {

    std::cerr << "ouch\n";
    std::cerr << re.what() << std::endl;

  }

}

/*
void tst_cholesky( int num, double tol ) {

  sherpa::Array2d< double > A( num, num ), L_t( num, num ), result( num, num ),
    L_;
  for ( int ii = 0; ii < num - 1; ++ii ) {
    A[ ii ][ ii ] = 2;
    A[ ii ][ ii + 1 ] = -1;
    A[ ii + 1 ][ ii ] = -1;
  }
  A[ num - 1 ][ num - 1 ] = 2;

  int isspd;
  jama::cholesky::cholesky( A, L_, isspd );
  for ( int ii = 0; ii < num; ++ii )
    for ( int jj = 0; jj < num; ++jj )
      L_t[ ii ][ jj ] = L_[ ii ][ jj ];
  jama::transpose( L_t );
  

  jama::matrix_multiply( L_, L_t, result );
  for ( int ii = 0; ii < num; ++ii )
    for ( int jj = 0; jj < num; ++jj )
      assert( sherpa::utils::Knuth_close_enough( A[ ii ][ jj ], result[ ii ][ jj ],
					tol ) );
}
*/

int main( int argc, char* argv[] ) {

  bool debug = false, printit=false;
  int num = 32, numiter = 10;
  double tol = 1.0e-8;

  tst_inverse( jama::lu::inverse<double>, numiter, num, tol, debug, "lu inv" );
  tst_inverse( jama::svd::inverse<double>, numiter, num, tol, debug,
	       "svd inv" );
  tst_solve( jama::lu::solve, numiter, num, tol, printit, "lu solve" );
  // tst_cholesky( num, tol );

  return 0;

}


#endif

//
// cp jama.hh tmp.cc; g++ -I.. -Wall -pedantic -ansi -O3 -DtestJama tmp.cc; rm -f tmp.cc; a.out
//
