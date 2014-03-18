/* mode: C; mode: fold */

/* This file is a subset of routines from the JDMath library by
 * John E. Davis, modified by John C. Houck (5/4/2000)
 * to work independently of JDMath
 *
 * Copyright (c) 1994, 1998 John E. Davis 
 *
 * You may distribute this file under the terms the GNU General Public
 * License.  See the file COPYING for more information.
 */
/*--------------------------------*-C-*---------------------------------*
 * File:
 *	fftn.c
 *
 *	multivariate complex Fourier transform, computed in place
 *	using mixed-radix Fast Fourier Transform algorithm.
 *
 *	Fortran code by:
 *	    RC Singleton, Stanford Research Institute, Sept. 1968
 *		NIST Guide to Available Math Software.
 *		Source for module FFT from package GO.
 *		Retrieved from NETLIB on Wed Jul  5 11:50:07 1995.
 *	translated by f2c (version 19950721) and with lots of cleanup
 *	to make it resemble C by:
 *	    MJ Olesen, Queen's University at Kingston, 1995-97
 */
/*{{{ Copyright: */
/*
 * Copyright(c)1995,97 Mark Olesen <olesen@me.QueensU.CA>
 *		Queen's Univ at Kingston (Canada)
 *
 * Permission to use, copy, modify, and distribute this software for
 * any purpose without fee is hereby granted, provided that this
 * entire notice is included in all copies of any software which is
 * or includes a copy or modification of this software and in all
 * copies of the supporting documentation for such software.
 *
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE AUTHOR NOR QUEEN'S
 * UNIVERSITY AT KINGSTON MAKES ANY REPRESENTATION OR WARRANTY OF ANY
 * KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS
 * FITNESS FOR ANY PARTICULAR PURPOSE.
 *
 * All of which is to say that you can do what you like with this
 * source code provided you don't try to sell it as your own and you
 * include an unaltered copy of this message (including the
 * copyright).
 *
 * It is also implicitly understood that bug fixes and improvements
 * should make their way back to the general Internet community so
 * that everyone benefits.
 *----------------------------------------------------------------------*/
/*}}}*/
/*{{{ notes: */
/*
 * Public:
 *	fft_free
 *	fftn / fftnf
 *	(these are documented in the header file)
 *
 * Private:
 *	fftradix / fftradixf
 *
 * ----------------------------------------------------------------------*
 * int fftradix (REAL Re[], REAL Im[], size_t nTotal, size_t nPass,
 *		 size_t nSpan, int iSign, size_t maxFactors,
 *		 size_t maxPerm);
 *
 * RE and IM hold the real and imaginary components of the data, and
 * return the resulting real and imaginary Fourier coefficients.
 * Multidimensional data *must* be allocated contiguously.  There is
 * no limit on the number of dimensions.
 *
 *
 * Although there is no limit on the number of dimensions, fftradix()
 * must be called once for each dimension, but the calls may be in
 * any order.
 *
 * NTOTAL = the total number of complex data values
 * NPASS  = the dimension of the current variable
 * NSPAN/NPASS = the spacing of consecutive data values while indexing
 *	the current variable
 * ISIGN - see above documentation
 *
 * example:
 * tri-variate transform with Re[n1][n2][n3], Im[n1][n2][n3]
 *
 *	fftradix (Re, Im, n1*n2*n3, n1,       n1, 1, maxf, maxp);
 *	fftradix (Re, Im, n1*n2*n3, n2,    n1*n2, 1, maxf, maxp);
 *	fftradix (Re, Im, n1*n2*n3, n3, n1*n2*n3, 1, maxf, maxp);
 *
 * single-variate transform,
 *    NTOTAL = N = NSPAN = (number of complex data values),
 *
 *	fftradix (Re, Im, n, n, n, 1, maxf, maxp);
 *
 * The data can also be stored in a single array with alternating
 * real and imaginary parts, the magnitude of ISIGN is changed to 2
 * to give correct indexing increment, and data [0] and data [1] used
 * to pass the initial addresses for the sequences of real and
 * imaginary values,
 *
 * example:
 *	REAL data [2*NTOTAL];
 *	fftradix (&data[0], &data[1], NTOTAL, nPass, nSpan, 2, maxf, maxp);
 *
 * for temporary allocation:
 *
 * MAXFACTORS	>= the maximum prime factor of NPASS
 * MAXPERM	>= the number of prime factors of NPASS.  In addition,
 *	if the square-free portion K of NPASS has two or more prime
 *	factors, then MAXPERM >= (K-1)
 *
 * storage in FACTOR for a maximum of 15 prime factors of NPASS. if
 * NPASS has more than one square-free factor, the product of the
 * square-free factors must be <= 210 array storage for maximum prime
 * factor of 23 the following two constants should agree with the
 * array dimensions.
 * ----------------------------------------------------------------------*/
/*}}}*/
/*{{{ Revisions: */
/*
 * 26 July 95	John Beale
 *	- added maxf and maxp as parameters to fftradix()
 *
 * 28 July 95	Mark Olesen <olesen@me.QueensU.CA>
 *	- cleaned-up the Fortran 66 goto spaghetti, only 3 labels remain.
 *
 *	- added fft_free() to provide some measure of control over
 *	  allocation/deallocation.
 *
 *	- added fftn() wrapper for multidimensional FFTs
 *
 *	- use -DFFT_NOFLOAT or -DFFT_NODOUBLE to avoid compiling that
 *	  precision. Note suffix `f' on the function names indicates
 *	  float precision.
 *
 *	- revised documentation
 *
 * 31 July 95	Mark Olesen <olesen@me.QueensU.CA>
 *	- added GNU Public License
 *	- more cleanup
 *	- define SUN_BROKEN_REALLOC to use malloc() instead of realloc()
 *	  on the first pass through, apparently needed for old libc
 *	- removed #error directive in favour of some code that simply
 *	  won't compile (generate an error that way)
 *
 * 1 Aug 95	Mark Olesen <olesen@me.QueensU.CA>
 *	- define FFT_RADIX4 to only have radix 2, radix 4 transforms
 *	- made fftradix /fftradixf () static scope, just use fftn()
 *	  instead.  If you have good ideas about fixing the factors
 *	  in fftn() please do so.
 *
 * 8 Jan 95	mj olesen <olesen@me.QueensU.CA>
 *	- fixed typo's, including one that broke scaling for scaling by
 *	  total number of matrix elements or the square root of same
 *	- removed unnecessary casts from allocations
 *
 * 10 Dec 96	mj olesen <olesen@me.QueensU.CA>
 *	- changes defines to compile *without* float support by default,
 *	  use -DFFT_FLOAT to enable.
 *	- shifted some variables to local scope	(better hints for optimizer)
 *	- added Michael Steffens <Michael.Steffens@mbox.muk.uni-hannover.de>
 *	  Fortran 90 module
 *	- made it simpler to pass dimensions for 1D FFT.
 *
 * 23 Feb 97	Mark Olesen <olesen@me.QueensU.CA>
 *	- removed the GNU Public License (see 21 July 1995 entry),
 *	  which should make it clear why I have the right to do so.
 *	- Added copyright notice and submitted to netlib
 *	- Moved documentation for the public functions to the header
 *	  files where is will always be available.
 *
 * 13 Mar 98    John E. Davis <davis@space.mit.edu>
 *      - Made changes to incorporate the code into the JDMath library.  In
 *        particular, the names of the public variables were changed to be
 *        consistent with the library.
 *
 * 5 May 2000   John C. Houck <houck@space.mit.edu>
 *      - 'borrowed' from the JDMath library;  renamed JDM memory allocation
 *        subroutine calls to use the versions in isismath.
 * 
 * ----------------------------------------------------------------------*/
/*}}}*/
/*#include "config.h"*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*#include "isis.h"*/
/*#include "isismath.h"*/
#define ISIS_FREE(p) do {if (p) free ((char *)p); p = NULL;} while (0)

/*{{{ defines/constants */
#ifndef M_PI
# define M_PI	3.14159265358979323846264338327950288
#endif

#ifndef SIN60
# define SIN60	0.86602540378443865	/* sin(60 deg) */
# define COS72	0.30901699437494742	/* cos(72 deg) */
# define SIN72	0.95105651629515357	/* sin(72 deg) */
#endif
/*}}}*/

/*{{{ static parameters - for memory management */
static size_t SpaceAlloced = 0;
static size_t MaxPermAlloced = 0;

/* temp space, (void *) since both float and double routines use it */
static void * Tmp0 = NULL;	/* temp space for real part */
static void * Tmp1 = NULL;	/* temp space for imaginary part */
static void * Tmp2 = NULL;	/* temp space for Cosine values */
static void * Tmp3 = NULL;	/* temp space for Sine values */
static int  * Perm = NULL;	/* Permutation vector */

#define NFACTOR	11
static int factor [NFACTOR];
/*}}}*/

/*{{{ fft_free() */
void
JDMfft_free (void)
{
   SpaceAlloced = MaxPermAlloced = 0;
   ISIS_FREE (Tmp0);
   ISIS_FREE (Tmp1);
   ISIS_FREE (Tmp2);
   ISIS_FREE (Tmp3);
   ISIS_FREE (Perm);
   Tmp0 = Tmp1 = Tmp2 = Tmp3 = NULL;
   Perm = NULL;
}
/*}}}*/

/*{{{ static factorize */

/* return the number of factors */
static int
factorize (int nPass, int * kt)
{
   int nFactor = 0;
   int j, jj;

   *kt = 0;
   /* determine the factors of n */
   while ((nPass % 16) == 0)	/* factors of 4 */
     {
	factor [nFactor++] = 4;
	nPass /= 16;
     }
   j = 3; jj = 9;		/* factors of 3, 5, 7, ... */
   do
     {
	while ((nPass % jj) == 0)
	  {
	     factor [nFactor++] = j;
	     nPass /= jj;
	  }
	j += 2;
	jj = j * j;
     }
   while (jj <= nPass);
   if (nPass <= 4)
     {
	*kt = nFactor;
	factor [nFactor] = nPass;
	if (nPass != 1)
	  nFactor++;
     }
   else
     {
	if (nPass - (nPass / 4 << 2) == 0)
	  {
	     factor [nFactor++] = 2;
	     nPass /= 4;
	  }
	*kt = nFactor;
	j = 2;
	do
	  {
	     if ((nPass % j) == 0)
	       {
		  factor [nFactor++] = j;
		  nPass /= j;
	       }
	     j = ((j + 1) / 2 << 1) + 1;
	  }
	while (j <= nPass);
     }
   if (*kt)
     {
	j = *kt;
	do
	  factor [nFactor++] = factor [--j];
	while (j);
     }

   return nFactor;
}

/*}}}*/

/* defines for double */
#define REAL		double
#define FFTN		JDMfftn
#define FFTNS		"JDMfftn"
#define FFTRADIX	fftradix
#define FFTRADIXS	"fftradix"
/* double precision routine */
static int
fftradix (double Re[], double Im[],
	  size_t nTotal, size_t nPass, size_t nSpan, int isign,
	  int maxFactors, int maxPerm);
#include "fftn.inc"

#undef REAL
#undef FFTN
#undef FFTNS
#undef FFTRADIX
#undef FFTRADIXS

/* defines for float */
#define REAL		float
#define FFTN		JDMfftnf	/* trailing 'f' for float */
#define FFTNS		"JDMfftnf"	/* name for error message */
#define FFTRADIX	fftradixf	/* trailing 'f' for float */
#define FFTRADIXS	"fftradixf"	/* name for error message */
/* float precision routine */
static int
fftradixf (float Re[], float Im[],
	   size_t nTotal, size_t nPass, size_t nSpan, int isign,
	   int maxFactors, int maxPerm);
#include "fftn.inc"
