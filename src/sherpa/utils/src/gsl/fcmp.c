/* sys/gsl_compare.c
 * 
 * Copyright (C) 2002 Gert Van den Eynde
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 * 
 * Based on fcmp 1.2.2 Copyright (c) 1998-2000 Theodore C. Belding
 * University of Michigan Center for the Study of Complex Systems
 * Ted.Belding@umich.edu
 *
 */

/*
#include <config.h>
#include <gsl/gsl_sys.h>
*/
#include <math.h>
#include "fcmp.h"

int _gsl_fcmp (const double x1, const double x2, const double epsilon) 
{

  int exponent, result;
  double delta, difference;

  /* Find exponent of largest absolute value */
  {
    double max = (fabs (x1) > fabs (x2)) ? x1 : x2;

    frexp (max, &exponent);
  }

  /* Form a neighborhood of size  2 * delta */

  delta = ldexp (epsilon, exponent);

  difference = x1 - x2;

  if (difference > delta)       /* x1 > x2 */
    {
      result = 1;
    }
  else if (difference < -delta) /* x1 < x2 */
    {
      result = -1;
    }
  else                          /* -delta <= difference <= delta */
    {
      result = 0;                 /* x1 ~=~ x2 */
    }
  return result;

}

int _sao_fcmp(const double x1, const double x2, const double epsilon) {

  int result;

  if ( x1 == x2 ) {
    /* might as well start at the simplest comparison... */
    result = 0;
  } else if ( 0 == x1 || 0 == x2 ) {

    /* x1 || x2 is 0 so the following substraction is kosher */
    if ( fabs( x1 - x2 ) < epsilon )
      result = 0;
    else {
      if ( x1 > x2 )
	result = 1;
      else
	result = -1;
    }

  } else {

    int exponent;
    double delta, difference;

    /* Find exponent of largest absolute value */
    {
      double max = (fabs (x1) > fabs (x2)) ? x1 : x2;

      frexp (max, &exponent);
    }

    /* Form a neighborhood of size  2 * delta */

    delta = ldexp (epsilon, exponent);

    difference = x1 - x2;

    if (difference > delta) {
      /* x1 > x2 */
      result = 1;
    } else if (difference < -delta) {
      /* x1 < x2 */
      result = -1;
    } else {
      /* -delta <= difference <= delta */
      result = 0;                 /* x1 ~=~ x2 */
    }

  }

  return result;

}
