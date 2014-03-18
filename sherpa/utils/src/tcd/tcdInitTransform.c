/*_C_INSERT_SAO_COPYRIGHT_HERE_(1998-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

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
