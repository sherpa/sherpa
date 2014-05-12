/*_C_INSERT_SAO_COPYRIGHT_HERE_(1998-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

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
