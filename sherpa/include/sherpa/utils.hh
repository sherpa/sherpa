//_C++_INSERT_SAO_COPYRIGHT_HERE_(2008)_
//_C++_INSERT_GPL_LICENSE_HERE_

#ifndef __sherpa_utils_hh__
#define __sherpa_utils_hh__

#include <cstdlib>
#include <cmath>

#ifdef __SUNPRO_CC
#include <sunmath.h>
#define NAN quiet_nan(0)
#endif

#include <vector>
#include <limits>

#include <sherpa/constants.hh>
#include <sherpa/myArray.hh>
#include "sherpa/fcmp.hh" 

#define ACOS(x)			std::acos(x)
#define ASIN(x)			std::asin(x)
#define ATAN(x)			std::atan(x)
#define COS(x)			std::cos(x)
#define ERF(x)			sherpa::utils::erf(x)
#define ERFC(x)			sherpa::utils::erfc(x)
#define EXP(x)			std::exp(x)
#define FABS(x)			std::fabs(x)
#define LGAMMA(x)      		sherpa::utils::lgamma(x)
#define LOG(x)			std::log(x)
#define LOG10(x)       		std::log10(x)
#define POW(x, y)		std::pow(x, y)
#define SIN(x)			std::sin(x)
#define SQRT(x)			std::sqrt(x)
#define TAN(x)			std::tan(x)


namespace sherpa { namespace utils {


  inline double erf( double x )
  {
    return ::erf( x );
  }


  inline double erfc( double x )
  {
    return ::erfc( x );
  }


  inline double lgamma( double x )
  {
    return ::lgamma( x );
  }


  /*
   *     **********
   *
   *     function enorm
   *
   *     given an n-vector x, this function calculates the
   *     euclidean norm of x.
   *
   *     the euclidean norm is computed by accumulating the sum of
   *     squares in three different sums. the sums of squares for the
   *     small and large components are scaled so that no overflows
   *     occur. non-destructive underflows are permitted. underflows
   *     and overflows do not occur in the computation of the unscaled
   *     sum of squares for the intermediate components.
   *     the definitions of small, intermediate and large components
   *     depend on two constants, rdwarf and rgiant. the main
   *     restrictions on these constants are that rdwarf**2 not
   *     underflow and rgiant**2 not overflow. the constants
   *     given here are suitable for every known computer.
   *
   *     the function statement is
   *
   *	Real precision function enorm(n,x)
   *
   *     where
   *
   *	n is a positive integer input variable.
   *
   *	x is an input array of length n.
   *
   *     subprograms called
   *
   *	fortran-supplied ... dabs,dsqrt
   *
   *     argonne national laboratory. minpack project. march 1980.
   *     burton s. garbow, kenneth e. hillstrom, jorge j. more
   *
   *     **********
   */
  template <typename ConstArrayType, typename DataType, typename IndexType>
  inline DataType enorm2( IndexType n, const ConstArrayType& x ) {

    // int i;
    DataType agiant, floatn, s1, s2, s3, xabs, x1max, x3max;
    DataType ans, temp;
    static DataType rdwarf = 3.834e-20;
    static DataType rgiant = 1.304e19;
    static DataType zero = 0.0;
    static DataType one = 1.0;
    /*double fabs(), sqrt();*/

    s1 = zero;
    s2 = zero;
    s3 = zero;
    x1max = zero;
    x3max = zero;
    floatn = n;
    agiant = rgiant / floatn;

    for (IndexType i = 0; i < n; i++) {
      xabs = std::fabs(x[i]);
      if ((xabs > rdwarf) && (xabs < agiant)) {
	/*
	 *	    sum for intermediate components.
	 */
	s2 += xabs * xabs;
	continue;
      }

      if (xabs > rdwarf) {
	/*
	 *	       sum for large components.
	 */
	if (xabs > x1max) {
	  temp = x1max / xabs;
	  s1 = one + s1 * temp * temp;
	  x1max = xabs;
	} else {
	  temp = xabs / x1max;
	  s1 += temp * temp;
	}
	continue;
      }
      /*
       *	       sum for small components.
       */
      if (xabs > x3max) {
	temp = x3max / xabs;
	s3 = one + s3 * temp * temp;
	x3max = xabs;
      } else {
	if (xabs != zero) {
	  temp = xabs / x3max;
	  s3 += temp * temp;
	}
      }
    }
    /*
     *     calculation of norm.
     */
    if (s1 != zero) {
      temp = s1 + (s2 / x1max) / x1max;
      // ans = x1max * sqrt(temp);
      ans = x1max * temp;
      return (ans);
    }
    if (s2 != zero) {
      if (s2 >= x3max)
	temp = s2 * (one + (x3max / s2) * (x3max * s3));
      else
	temp = x3max * ((s2 / x3max) + (x3max * s3));
      // ans = sqrt(temp);
      ans = temp;
    } else {
      // ans = x3max * sqrt(s3);
      ans = x3max * s3;
    }
    return (ans);
    /*
     *     last card of function enorm2.
     */

  }

  template <typename ConstArrayType, typename DataType, typename IndexType>
  inline DataType simple_sum( IndexType num, const ConstArrayType& vals ) {

    DataType sum = vals[0];

    for ( IndexType ii = num - 1; 0 < ii; ii-- )
      sum += vals[ ii ];

    return sum;

  }


  template <typename ConstArrayType, typename DataType, typename IndexType>
  inline DataType kahan_sum( IndexType num, const ConstArrayType& vals ) {

    DataType sum = vals[ 0 ];
    DataType correction = 0.0;

    for ( IndexType ii = 1; ii < num; ii++ ) {

      DataType y = vals[ ii ] - correction;
      DataType t = sum + y;
      correction = t - sum - y;
      sum = t;

    }

    return sum;

  }


  template <typename DataType, typename ConstArrayType>
  inline int radius2( const ConstArrayType& p,
		      DataType x0, DataType x1, DataType& val )
  {
  
    register DataType deltaX = x0 - p[1];
    register DataType deltaY = x1 - p[2];
    register DataType p_four = p[4];
    if( p[3] != 0 ) {
      while( p_four >= 2*PI ) {
	p_four -= 2*PI;
      }
      while( p_four < 0.0 ) {
	p_four += 2*PI;
      }
      register DataType cosTheta = COS(p_four);
      register DataType sinTheta = SIN(p_four);
      register DataType newX = deltaX * cosTheta + deltaY * sinTheta;
      register DataType newY = deltaY * cosTheta - deltaX * sinTheta;
      if( p[3] == 1 )
	return EXIT_FAILURE;
      else {
	register DataType ellip2 = (1. - p[3]) * (1. - p[3]);
	val = (newX * newX * ellip2 + newY * newY) / ellip2;
	return EXIT_SUCCESS;
      }
    } else {
      val = deltaX * deltaX + deltaY * deltaY;
      return EXIT_SUCCESS;
    }
  
    val = 0.0;
    return EXIT_SUCCESS;
  
  }

  template <typename DataType, typename ConstArrayType>
  inline int sigmaradius2( const ConstArrayType& p,
			   DataType x0, DataType x1, DataType& val )
  {
    DataType sigmaa = p[0];
    DataType sigmab = p[1];
    register DataType deltaX = x0 - p[2];
    register DataType deltaY = x1 - p[3];
    register DataType theta = p[4];

    if ( 0 == sigmaa || 0 == sigmab ) {
      val = 0.0;
      return EXIT_FAILURE;
    }

    while( theta >= 2*PI ) {
      theta -= 2*PI;
    }
    while( theta < 0.0 ) {
      theta += 2*PI;
    }
    register DataType cosTheta = COS(theta);
    register DataType sinTheta = SIN(theta);
    register DataType newX = deltaX * cosTheta + deltaY * sinTheta;
    register DataType newY = deltaY * cosTheta - deltaX * sinTheta;
    DataType a = newX / sigmaa;
    DataType b = newY / sigmab;
    val = a*a + b*b;
    return EXIT_SUCCESS;
  }


  template <typename DataType, typename ConstArrayType>
  inline int radius( const ConstArrayType& p,
		     DataType x0, DataType x1, DataType& val )
  {

    if( EXIT_SUCCESS != radius2( p, x0, x1, val ) ) {
      return EXIT_FAILURE;
    }

    val = SQRT( val );

    return EXIT_SUCCESS;
    
  }

  // {
  
  //   register DataType deltaX = x0 - p[1];
  //   register DataType deltaY = x1 - p[2];
  //   register DataType p_four = p[4];
  //   if( p[3] != 0 ) {
  //     while( p_four >= 2*PI ) {
  // 	p_four -= 2*PI;
  //     }
  //     while( p_four < 0.0 ) {
  // 	p_four += 2*PI;
  //     }
  //     register DataType cosTheta = COS(p_four);
  //     register DataType sinTheta = SIN(p_four);
  //     register DataType newX = deltaX * cosTheta + deltaY * sinTheta;
  //     register DataType newY = deltaY * cosTheta - deltaX * sinTheta;
  //     if( p[3] == 1 )
  // 	return EXIT_FAILURE;
  //     else {
  // 	val = SQRT(newX*newX*(1-p[3])*(1-p[3])+newY*newY)/(1-p[3]);
  // 	return EXIT_SUCCESS;
  //     }
  //   } else {
  //     val = SQRT( deltaX * deltaX + deltaY * deltaY );
  //     return EXIT_SUCCESS;
  //   }
  
  //   val = 0.0;
  //   return EXIT_SUCCESS;
  
  // }



//     Copyright (C) 2000 Massachusetts Institute of Technology 
 
//     Author:  John E. Davis <davis@space.mit.edu>

//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details. 

//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

//     2007: Modified to utilize templates, B Refsdal


  // integrate from: (flo, fhi, fy) to: (tlo, thi) yielding => ty
  template <typename DataType, typename ConstArrayType,
	    typename ArrayType, typename IndexType>
  inline int rebin_histogram (const ConstArrayType& fy,
			      const ConstArrayType& flo,
			      const ConstArrayType& fhi,
			      IndexType nf,
			      ArrayType& ty,
			      const ConstArrayType& tlo,
			      const ConstArrayType& thi,
			      IndexType nt)
  {
    IndexType f, t;
    DataType tol = std::numeric_limits< DataType >::epsilon();

    f = 0;
    for (t = 0; t < nt; t++)
      {
        DataType t0, t1, s;

        t0 = tlo[t];
        t1 = thi[t];
	
        /* Accumulate sum over 'from' bins
         * which overlap 'to' bin t
         */
        s = 0.0;
        for ( ;f < nf; f++)
          {
	    DataType f0, f1;
	    DataType min_max, max_min;
	    
	    f0 = flo[f];
	    f1 = fhi[f];
	    
	    // t0 > f1
	    if ( sao_fcmp( t0, f1, tol ) == 1 )
	      continue;

	    // f0 > t1
	    if ( sao_fcmp( f0, t1, tol ) == 1 )
	      break;
	    
	    // t0 > f0
	    if ( sao_fcmp( t0, f0, tol ) == 1 )
	      max_min = t0;
	    else
	      max_min = f0;
	    
	    // t1 < f1
	    if ( sao_fcmp( t1, f1, tol ) == -1 )
	      min_max = t1;
	    else
	      min_max = f1;
	    
	    /* prevent division by zero */
	    // f0 == f1
	    if ( sao_fcmp( f0, f1, tol ) == 0 )
	      return EXIT_FAILURE;
	    
	    s += fy[f] * (min_max - max_min) / (f1 - f0);
	    
	    /* s += fy[f] * overlap_frac (f0, f1, t0, t1); */
	    
	    if (f1 > t1)
	      break;
          }
	
        ty[t] = s;
      }
    
    return EXIT_SUCCESS;
  }


  template<typename ConstFloatArrayType, typename FloatType>
  inline int neville( int n,
	       const ConstFloatArrayType& x, const ConstFloatArrayType& y,
	       FloatType xinterp, FloatType& answer ) {

    // ConstFloatArrayType may not be aligned in memory!
    // Sherpa Arrays overload [] to make the correct jump.

    //std::vector<FloatType> p( y, y+n );

    std::vector<FloatType> p( n );
    for ( int ii = 0; ii < n; ii++ )
      p[ii] = y[ii];

    for ( int jj = 1; jj < n; jj++ )
      for ( int ii = n-1; ii >= jj; ii-- ) {
	FloatType denom = x[ii] - x[ii-jj];
	if ( 0.0 == denom )
	  return EXIT_FAILURE;
	p[ii] = ( (xinterp-x[ii-jj])*p[ii] - (xinterp-x[ii])*p[ii-1] ) / denom;
      }
    
    answer = p[n-1];
    return EXIT_SUCCESS;
    
  }
    
  template<typename ConstFloatArrayType, typename FloatType>
  inline int neville2d( int m, int n, const ConstFloatArrayType& x, 
			const ConstFloatArrayType& y, 
			const sherpa::Array2d< FloatType >& fxy,
			FloatType xinterp, FloatType yinterp,
			FloatType& answer ) {

    std::vector< FloatType > tmp( m );
    for ( int ii = 0; ii < m; ++ii )
      if ( EXIT_SUCCESS != neville( n, y, fxy[ ii ], yinterp, tmp[ ii ] ) )
	return EXIT_FAILURE;
    return neville( m, x, tmp, xinterp, answer );

  }

  //
  // The following text was taken verbatim from:
  //        
  // http://www.boost.org/doc/libs/1_35_0/libs/test/doc/components/test_tools/floating_point_comparison.html#Introduction
  //
  // In most cases it is unreasonable to use an operator==(...)
  // for a floating-point values equality check. The simple solution
  // like abs(f1-f2) <= e does not work for very small or very big values.
  // This floating-point comparison algorithm is based on the more
  // confident solution presented by D. E. Knuth in 'The art of computer
  // programming (vol II)'. For a given floating point values u and v and
  // a tolerance e:
  //    
  // | u - v | <= e * |u| and | u - v | <= e * |v|                    (1)
  // defines a "very close with tolerance e" relationship between u and v
  //    
  // | u - v | <= e * |u| or   | u - v | <= e * |v|                   (2)
  // defines a "close enough with tolerance e" relationship between
  // u and v. Both relationships are commutative but are not transitive.
  // The relationship defined by inequations (1) is stronger that the
  // relationship defined by inequations (2) (i.e. (1) => (2) ).
  //
  template< typename T >
  T safe_div( T num, T denom ) {

    T my_max = std::numeric_limits< T >::max( );
    T my_min = std::numeric_limits< T >::min( );

    // avoid overflow
    if ( denom < 1 && num > denom * my_max )
      return my_max;

    // avoid underflow
    if ( 0.0 == num || denom > 1 && num < denom * my_min )
      return T(0);
    
    return num / denom;

  }

  //
  // more strict then ||
  //
  template< typename T >
  inline bool Knuth_more_strict_cmp( T diff_x, T diff_y, T tol ) {
    return diff_x <= tol && diff_y <= tol;
  }
  //
  // more lenient then &&
  //
  template< typename T >
  inline bool Knuth_less_strict_cmp( T diff_x, T diff_y, T tol ) {
    return diff_x <= tol || diff_y <= tol;
  }

  template< typename T >
  inline bool Knuth_close( T x, T y, T tol,
			   bool (*cmp)( T, T, T )=Knuth_less_strict_cmp ) {
    T diff = std::fabs( x - y );
    if ( 0.0 == x || 0.0 == y )
      return diff <= tol;
    T diff_x = safe_div( diff, x );
    T diff_y = safe_div( diff, y );
    return cmp( diff_x, diff_y, tol );
  }

  template <typename ConstFloatType, typename FloatType,
	    typename ConstIntType, typename IndexType>
  inline int sum_intervals (const ConstFloatType *fy,
			    const ConstIntType *indx0,
			    const ConstIntType *indx1,
			    IndexType size,
			    FloatType *ty)
  {

    for(int ii = 0; ii < size; ++ii) {
      if ( indx0[ ii ] > indx1[ ii ] ) {
	//continue;
	return EXIT_FAILURE;
      }
      for(int jj = indx0[ ii ]; jj <= indx1[ ii ]; ++jj)
	ty[ ii ] += fy[ jj ];
    }
    return EXIT_SUCCESS;

  }

  // histogram bin containing x:
  template <typename FloatType, typename ConstArrayType, typename IntType>
  inline IntType find_bin (FloatType x,
			     const ConstArrayType& lo,
			     const ConstArrayType& hi,
			     IntType n)
  {
    IntType n0, n1, n2;
    FloatType tol = std::numeric_limits< FloatType >::epsilon();

    if (!lo || !hi || n <= 0)
      return -1;

#ifdef __SUNPRO_CC
    if (isnan(x))
#else
    if (std::isnan(x))
#endif 
     return -1;
    
    if (sao_fcmp(x, lo[0], tol) == -1 ||
    	sao_fcmp(hi[n-1], x, tol) == -1 )
      return -1;
    
    n0 = 0;
    n1 = n;
    
    while (n1 > n0 + 1) {
    	n2 = (n0 + n1) / 2;
    	if (sao_fcmp(x, hi[n2], tol) == -1) {
    	    if (sao_fcmp(lo[n2], x, tol) <= 0)
    	      return n2;
    	    n1 = n2;
    	  }
    	else n0 = n2;
      }
    
    return n0;
  }

  template <typename ConstArrayType, typename FloatType,
	    typename IntArrayType, typename IntType>
  inline int histogram1d (const ConstArrayType& bv,
			  const ConstArrayType& x_lo,
			  const ConstArrayType& x_hi,
			  IntArrayType& res) {
    IntType n, nbins;
    
    n = bv.get_size();
    nbins = x_lo.get_size();

    for (IntType i = 0; i < n; i++)
      {
        FloatType t = bv[i];
        IntType k = find_bin (t, x_lo, x_hi, nbins);
	if (k >= n)
	  return EXIT_FAILURE;
	if (k >= 0)
	  res[k] += 1;
      }
    
    return EXIT_SUCCESS;
  }

  template <typename ConstArrayType, typename FloatType,
	    typename IntArrayType, typename IntType>
  inline int histogram2d (const ConstArrayType& py_x,
			  const ConstArrayType& py_y,
			  const ConstArrayType& grid_x,
			  const ConstArrayType& grid_y,
			  IntArrayType& res) {
    IntType i, n, nx, ny;
    FloatType xmax, ymax, *x, *y, *bx, *by;
    
    n = (IntType)py_x.get_size();
    nx = (IntType)grid_x.get_size();
    ny = (IntType)grid_y.get_size();

    bx = (FloatType*) &py_x[0];
    by = (FloatType*) &py_y[0];
    x = (FloatType*) &grid_x[0];
    y = (FloatType*) &grid_y[0];
    
    xmax = x[nx-1];
    ymax = y[ny-1];

    for (i = 0; i < n; i++)
      {
        FloatType b_x = bx[i];
        FloatType b_y = by[i];
        IntType ix, iy, k;
	
        if (b_x >= xmax)
          ix = nx-1;
        else if ((ix = find_bin (b_x, x, x+1, nx-1)) < 0)
          continue;
	
        if (b_y >= ymax)
          iy = ny-1;
        else if ((iy = find_bin (b_y, y, y+1, ny-1)) < 0)
          continue;
	
        k = iy + ny * ix;
	
        res[k] += 1;
      }
    
    return EXIT_SUCCESS;
  }
    
  }  }  // namespace utils, namespace sherpa


#endif // __sherpa_utils_hh__
