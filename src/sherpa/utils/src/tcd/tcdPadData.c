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
 * FILE NAME:  misc/tcdPadData.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains files that are needed to pad a data array to a
 * specified size.
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
  +-------------------------------------------------
  +
  + free memory allocated by tcdPadData routines
  +
  +-------------------------------------------------
  */
int tcdFreePadData(
		   void *data,   /* i/o: pointer to data array    */
		   long **lAxes  /* i/o: pointer to legths array  */
		   )
{
  void **p;
  p = data;

  if (*p )    free (*p);
  if (*lAxes) free(*lAxes);

  *p     = NULL;
  *lAxes = NULL;

  return(tcdSUCCESS);

}



/*
  +-----------------------------------------
  +
  + Pad data to the next power of 2^(N+param)
  +
  + This is a wrapper for the tcdPadData routine.
  +
  +-----------------------------------------
  */
int tcdPadData2N(
		 tcdDATATYPE dtype, /* i: input data type                */
		 void  *data,       /* i: input data array               */
		 long   nAxes,      /* i: number of axes                 */
		 long  *lAxes,      /* i: length of axes                 */
		 long   param,      /* i: increase to which power of 2   */
		 void  *output,     /* o: output array same type as input*/
		 long **nlAxes      /* o: new axes output lengths        */
		 )
{

  int status;
  /* std error checking */
  
  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  /* do padding */

  return( tcdPadData( tcdPAD2N, &param, dtype, data, nAxes,
		      lAxes, output, nlAxes ) );
}



/*
  +----------------------------------------------
  +
  + Pad data to specified length.
  +
  + This is a wrapper for the tcdPadData routine
  +
  +----------------------------------------------
  */

int tcdPadDataSpec(
		   tcdDATATYPE dtype, /* i: input data type               */
		   void  *data,       /* i: pointer to data               */
		   long   nAxes,      /* i: number of axes                */
		   long  *lAxes,      /* i: lengths of axes               */ 
		   long  *newAxes,    /* i: new lengths of axes           */
		   void  *output      /* o: output data array             */
		   )
{
  long *nlAxes;
  int status;
  
  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );


  status = tcdPadData( tcdPAD2SPEC, newAxes, dtype, data, nAxes,
		       lAxes, output, &nlAxes ); 

  if (nlAxes ) free( nlAxes);

  return(status);

}   


/*
  +----------------------------------------------------
  +
  + Pad data with the additional amount of space specified
  +
  + This is a wrapper for tcdPadData
  +
  +----------------------------------------------------
  */
int tcdPadDataWith(
		   tcdDATATYPE dtype,  /* i: input data type              */
		   void  *data,        /* i: pointer to data array        */
		   long   nAxes,       /* i: number of axes               */
		   long  *lAxes,       /* i: lengths of axes              */
		   long  *addAxes,     /* i: amount to increase axes by   */
		   void  *output,      /* o: pointer to output array      */
		   long **nlAxes       /* o: new lengths array            */
		   )
{

  int status;

  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  return( tcdPadData( tcdPAD2KERNEL, addAxes, dtype, data, nAxes,
		      lAxes, output, nlAxes ) );

}



/*
  +--------------------------------------------------
  +
  + Generic pad data routine.
  +
  +--------------------------------------------------
  */  
int tcdPadData(
	       tcdPADTYPE   ptype, /* i: how to pad data               */
	       long        *param, /* i: parameters needed to pad data */
	       tcdDATATYPE  dtype, /* i: input data type               */
	       void        *data,  /* i: pointer to input data array   */
	       long         nAxes, /* i: number of axes                */
	       long        *lAxes, /* i: lengths of axes               */
	       void       *output, /* o: output data array             */
	       long       **nlAxes /* o: new axes lengths              */
	       )
{
  long ii;  /* generic loop variable */

  double index;
  long   newIndex; /* index in padded data array */

  long newTotal;
  long oldTotal;
  long *oldPixel; /* pixel location in old array */



  short  *data_s, **out_s;  /* pointers to supported data types */
  long   *data_l, **out_l;
  float  *data_f, **out_f;
  double *data_d, **out_d;
  char   *data_b, **out_b;

  int status;

  data_s = NULL;
  out_s = NULL;
  data_l = NULL;
  out_l = NULL;
  data_f = NULL;
  out_f = NULL;
  data_d = NULL;
  out_d = NULL;
  data_b = NULL;
  out_b = NULL;

  /* error checking */

  status = tcdCheckData( data, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  *nlAxes = ( long *)calloc( nAxes, sizeof(long));
  if ( *nlAxes == NULL) return ( tcdERROR_ALLOC );

  if ( param == NULL ) return( tcdERROR_NULLPTR );


  newTotal = 1;
  oldTotal = 1;
  
  /* determine new size of data array */

  switch ( ptype )
    {

    case tcdPAD2KERNEL:
      for (ii=0;ii<nAxes; ii++)
	{
	  if ( param[ii] <= 0 ) return( tcdERROR_PADLTOLD );
	  *(*nlAxes+ii)  = param[ii] + lAxes[ii];
	  newTotal      *= *(*nlAxes+ii);
	  oldTotal      *= lAxes[ii];
	}
      break; /* end tcdPAD2KERNEL */

    case tcdPAD2N:
      for (ii=0;ii<nAxes; ii++)
	{
	  index = log((double) lAxes[ii])/log((double )2.0);
	  index = (double) ((int ) index == index ? index: (int)index+1);
	  
	  (*nlAxes)[ii] = (long )pow((double )2.0, index + (*param));
	  newTotal *=  (*nlAxes)[ii];
	  oldTotal *= lAxes[ii];
	}
      break; /* end tcdPAD2N */

    case tcdPAD2SPEC:
      for (ii=0; ii<nAxes; ii++)
	{
	  if ( param[ii] < lAxes[ii] ) return( tcdERROR_PADLTOLD );

	  (*nlAxes)[ii] = param[ii];
	  newTotal *=  (*nlAxes)[ii];
	  oldTotal *= lAxes[ii];
	} /* end for ii */
      break; /* end tcdPAD2SPEC */

    default:
      return( tcdERROR_UNKWNPAD );

    } /* end switch */

  


  /* allocate memory */


  switch ( dtype )
    {

    case tcdFLOAT:
      data_f = (float *) data;
      out_f = output;
      *out_f = (float *)calloc(newTotal, sizeof(float ));
      if ( *out_f == NULL) return(tcdERROR_ALLOC);
      break;

    case tcdBYTE:
      data_b = (char *) data;
      out_b = output;
      *out_b = (char *)calloc(newTotal, sizeof(char ));
      if ( *out_b == NULL) return(tcdERROR_ALLOC);
      break;

    case tcdDOUBLE:
      data_d = (double *) data;
      out_d = output;
      *out_d = (double *)calloc(newTotal, sizeof(double ));
      if ( *out_d == NULL) return(tcdERROR_ALLOC);
      break;

    case tcdSHORT:
      data_s = (short *) data;
      out_s = output;
      *out_s = (short *) calloc(newTotal, sizeof(short ));
      if ( *out_s == NULL) return(tcdERROR_ALLOC);
      break;

    case tcdLONG:
      data_l = (long *) data;
      out_l = output;
      *out_l = (long *)calloc(newTotal, sizeof(long ));
      if ( *out_l == NULL) return(tcdERROR_ALLOC);
      break;

    default:
      return( tcdERROR_UNSUPORTTYPE );
    }

  
  oldPixel = (long *)calloc(nAxes, sizeof(long));
  if ( oldPixel == NULL ) return(tcdERROR_ALLOC);

  /* copy data */

  for ( ii=0;ii< oldTotal; ii++)
    {

      status = tcdOffsetToPixel( nAxes, lAxes, NULL, ii, oldPixel );
      if ( status != tcdSUCCESS ) return(status);

      status = tcdPixelToOffset( nAxes, *nlAxes, NULL, oldPixel, &newIndex );
      if ( status != tcdSUCCESS ) return(status);

      switch ( dtype )
	{

	case tcdBYTE:
	  (*out_b)[ newIndex] = data_b[ ii ];
	  break;

	case tcdSHORT:
	  (*out_s)[ newIndex] = data_s[ ii ];
	  break;

	case tcdLONG:
	  (*out_l)[ newIndex] = data_l[ ii ];
	  break;

	case tcdDOUBLE:
	  (*out_d)[ newIndex] = data_d[ ii ];
	  break;

	case tcdFLOAT:
	  (*out_f)[ newIndex] = data_f[ ii ];
	  break;


	default:
	  return( tcdERROR_UNSUPORTTYPE );
	} /* end switch */

    } /* end for ii */
  
  free( oldPixel );
  
  return( tcdSUCCESS );

}

