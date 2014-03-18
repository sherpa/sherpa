#ifndef fcmp_hh
#define fcmp_hh

#include <cmath>

/*
 fcmp
 Copyright (c) 1998-2000 Theodore C. Belding
 University of Michigan Center for the Study of Complex Systems
 <mailto:Ted.Belding@umich.edu>
 <http://www-personal.umich.edu/~streak/>		

 This file is part of the fcmp distribution. fcmp is free software;
 you can redistribute and modify it under the terms of the GNU Library
 General Public License (LGPL), version 2 or later.  This software
 comes with absolutely no warranty. See the file COPYING for details
 and terms of copying.

 File: fcmp.c

 Description: see fcmp.h and README files.
*/
template< typename Real >
int fcmp( Real x1, Real x2, Real epsilon ) {
  int exponent;
  Real delta;
  Real difference;
  
  /* Get exponent(max(fabs(x1), fabs(x2))) and store it in exponent. */

  /* If neither x1 nor x2 is 0, */
  /* this is equivalent to max(exponent(x1), exponent(x2)). */

  /* If either x1 or x2 is 0, its exponent returned by frexp would be 0, */
  /* which is much larger than the exponents of numbers close to 0 in */
  /* magnitude. But the exponent of 0 should be less than any number */
  /* whose magnitude is greater than 0. */
  
  /* So we only want to set exponent to 0 if both x1 and */
  /* x2 are 0. Hence, the following works for all x1 and x2. */

  frexp(fabs(x1) > fabs(x2) ? x1 : x2, &exponent);

  /* Do the comparison. */

  /* delta = epsilon * pow(2, exponent) */

  /* Form a neighborhood around x2 of size delta in either direction. */
  /* If x1 is within this delta neighborhood of x2, x1 == x2. */
  /* Otherwise x1 > x2 or x1 < x2, depending on which side of */
  /* the neighborhood x1 is on. */
  
  delta = ldexp(epsilon, exponent); 
  
  difference = x1 - x2;

  if (difference > delta)
    return 1; /* x1 > x2 */
  else if (difference < -delta) 
    return -1;  /* x1 < x2 */
  else /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
}

template< typename Real >
int sao_fcmp( const Real x1, const Real x2, const Real epsilon ) {

  if ( x1 == x2 ) {
    // might as well start at the simplest comparison... 
    return 0;
  } else if ( static_cast< Real>( 0 ) == x1 || 
	      static_cast< Real >( 0 ) == x2 ) {

    // x1 || x2 is 0 so the following substraction is kosher
    if ( fabs( x1 - x2 ) < epsilon )
      return 0;
    else {
      if ( x1 > x2 )
	return 1;
      else
	return -1;
    }

  } else {

    return fcmp( x1, x2, epsilon );

  }

}


// template< typename Real >
// bool approximately_equal( Real a, Real b, Real epsilon ) {
//   return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// template< typename Real >
// bool essentially_equal( Real a, Real b, Real epsilon ) {
//   return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon );
// }

// template< typename Real >
// bool essentially_greater_than( Real a, Real b, Real epsilon ) {
//   return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

// template< typename Real >
// bool essentially_less_than( Real a, Real b, Real epsilon ) {
//   return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
// }

#endif
