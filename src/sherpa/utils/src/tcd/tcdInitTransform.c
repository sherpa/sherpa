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
 * FILE NAME:  transform/tcdInitTransform.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * Allocate and free memory needed by transform routines
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
  +---------------------------------------------
  +
  + Allocate and copy data to complex array from 
  + real and imaginary parts
  +
  +---------------------------------------------
  */   
int tcdInitTransform(
		     tcdDATATYPE dType, /* i: input data type        */
		     void *real,        /* i: real part of data      */
		     void *imag,        /* i: imaginary part of data */
		     long  nAxes,       /* i: number of axes         */
		     long *lAxes,       /* i: length of axes         */
		     tcdComplex **data  /* o: init. output array     */
		     )
{
  long ii;
  long nTotal = 1;
  int status;
  
  /* error checking */
  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status);

  /* determine size */
  for (ii=0;ii<nAxes; ii++) nTotal *= lAxes[ii];

  /* alloc memory */
  *data = ( tcdComplex *)calloc(nTotal, sizeof( tcdComplex));
  if (*data == NULL ) return( tcdERROR_ALLOC );

  /* copy data */
  return( tcdCastToComplex( dType, real, imag, nAxes, lAxes, *data ) );

}

/* double precision version */
int tcdInitTransformD(
		     tcdDATATYPE dType, /* i: input data type        */
		     void *real,        /* i: real part of data      */
		     void *imag,        /* i: imaginary part of data */
		     long  nAxes,       /* i: number of axes         */
		     long *lAxes,       /* i: length of axes         */
		     tcdDComplex **data  /* o: init. output array     */
		     )
{
  long ii;
  long nTotal = 1;
  int status;
  
  /* error checking */
  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status);

  /* determine size */
  for (ii=0;ii<nAxes; ii++) nTotal *= lAxes[ii];

  /* alloc memory */
  *data = ( tcdDComplex *)calloc(nTotal, sizeof( tcdDComplex));
  if (*data == NULL ) return( tcdERROR_ALLOC );

  /* copy data */
  return( tcdCastToDComplex( dType, real, imag, nAxes, lAxes, *data ) );

}



/*
  +---------------------------------------------
  +
  + free memory for complex data set
  +
  +---------------------------------------------
  */
int tcdFreeTransform(
		     tcdComplex **data  /* i/o: data array */
		     )
{
  if ( *data ) free(*data);

  *data = NULL;
  
  return( tcdSUCCESS );
}


/* double precision */
int tcdFreeTransformD(
		     tcdDComplex **data  /* i/o: data array */
		     )
{
  if ( *data ) free(*data);

  *data = NULL;
  
  return( tcdSUCCESS );
}
