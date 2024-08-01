/*                                                                
**  Copyright (C) 1998-2007  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


/*H*****************************************************************
 *
 * FILE NAME:  convolve/misc/tcdInitConvolveOut.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * this file contains routines needed allocate and free memory used
 * needed by the convolution routines 
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.        Date
 ----------       -----
 preClearCase     30March1998


 *
 H***************************************************************** */


#include "tcd.h"
#include "tcd_private.h"

/*
  +----------------------------------------------------
  +
  + Allocate memory for convolve output 
  +
  +----------------------------------------------------
  */
int tcdInitConvolveOutput(
			  long    nAxes, /* i: number of axes */
			  long   *lAxes, /* i: length of axes */
			  float **output /* o: output array   */
			  )
{
     
  long ii;
  long ntotal=1;
  int status;

  /* check input */

  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS) return(status);

  /* determine size */

  for ( ii=0;ii<nAxes; ii++) ntotal *= lAxes[ii];

  /* allocate memory */

  *output = ( float *)calloc( ntotal, sizeof(float ));
  if ( *output == NULL ) return( tcdERROR_ALLOC);

  return( tcdSUCCESS);

}

/* double precision */
int tcdInitConvolveOutputD(
			  long    nAxes, /* i: number of axes */
			  long   *lAxes, /* i: length of axes */
			  double **output /* o: output array   */
			  )
{
     
  long ii;
  long ntotal=1;
  int status;

  /* check input */

  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS) return(status);

  /* determine size */

  for ( ii=0;ii<nAxes; ii++) ntotal *= lAxes[ii];

  /* allocate memory */

  *output = ( double *)calloc( ntotal, sizeof(double ));
  if ( *output == NULL ) return( tcdERROR_ALLOC);

  return( tcdSUCCESS);

}



/*
  +--------------------------------------------
  +
  + Free memory for output array 
  +
  +--------------------------------------------
  */
int tcdFreeConvolveOutput(
			  float **output /* i/o: array to free */
			  )
{

  if (*output) free( *output );
  *output = NULL;

  return( tcdSUCCESS );

}

/* double precision */
int tcdFreeConvolveOutputD(
			  double **output /* i/o: array to free */
			  )
{

  if (*output) free( *output );
  *output = NULL;

  return( tcdSUCCESS );

}
