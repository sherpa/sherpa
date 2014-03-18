/*_C_INSERT_SAO_COPYRIGHT_HERE_(1998-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

/*H*****************************************************************
 * FILE NAME:  misc/tcdError.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains routines to check and report error conditions.
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
  +-------------------------------------
  +
  + check that the data input is properly defined 
  +
  +---------------------------------------
  */
int tcdCheckData(
		 void *data,  /* i: data array      */
		 long  nAxes, /* i: number of axes  */
		 long *lAxes  /* i: length of axes  */
		 )
{
  if ( data  == NULL ) return( tcdERROR_NULLPTR );

  return( tcdCheckAxes( nAxes, lAxes ));

}



/*
  +---------------------------------------
  +
  + check that the number of axes is not 0
  +
  +---------------------------------------
  */
int tcdCheckAxes(
		 long nAxes, /* i: number of axes */
		 long *lAxes /* i: length of axes */
		 )
{
  long ii;

  if ( nAxes <= 0    ) return( tcdERROR_NAXES0  );
  if ( lAxes == NULL ) return( tcdERROR_NULLPTR );
  for (ii=0;ii<nAxes; ii++)
    {
      if (lAxes[ii] <= 0) return( tcdERROR_LAXES0 );
    }

  return( tcdSUCCESS );
}
