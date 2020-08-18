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
 * FILE NAME:  misc/tcdCastArray.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * This file contains routines needed cast arrays from 1 type to another
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
  +-----------------------------------------------------
  +
  + Cast an array from 1 type to another.
  +
  +-----------------------------------------------------
  */
int tcdCastArray( tcdDATATYPE inType,  /* i: input data type           */
		  void *inData,        /* i: input data                */
		  long  nAxes,         /* i: number of axes            */
		  long *lAxes,         /* i: length of axes            */
		  tcdDATATYPE outType, /* i: output data type          */
		  void *outData        /* o: output data array         */
		  )
{

  long nTotal;  /* total length of array */
  long ii;      /* loop variable */
  int status;


  char        **outData_b, *inData_b;
  short       **outData_s, *inData_s;  /* pointers to alloced data types */
  long        **outData_l, *inData_l;
  float       **outData_f, *inData_f;
  double      **outData_d, *inData_d;
  tcdComplex  **outData_c, *inData_c;
  tcdDComplex **outData_D, *inData_D;

  outData_b = NULL;
  inData_b = NULL;
  outData_s = NULL;
  inData_s = NULL;
  outData_l = NULL;
  inData_l = NULL;
  outData_f = NULL;
  inData_f = NULL;
  outData_d = NULL;
  inData_d = NULL;
  outData_c = NULL;
  inData_c = NULL;
  outData_D = NULL;
  inData_D = NULL;



  status = tcdCheckData( inData, nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  nTotal = 1;
  for (ii=0; ii<nAxes; ii++) nTotal *= lAxes[ii];


  /* allocate memory */

  switch ( outType )
    {
    case tcdSHORT:
      outData_s = outData;
      *outData_s = (short *)calloc( nTotal , sizeof( short ));
      if (*outData_s == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdBYTE:
      outData_b = outData;
      *outData_b = (char *)calloc( nTotal , sizeof( char ));
      if (*outData_b == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdLONG:
      outData_l = outData;
      *outData_l = (long *)calloc( nTotal , sizeof( long ));
      if (*outData_l == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdFLOAT:
      outData_f = outData;
      *outData_f = (float *)calloc( nTotal , sizeof( float ));
      if (*outData_f == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdDOUBLE:
      outData_d = outData;
      *outData_d = (double *)calloc( nTotal , sizeof( double ));
      if (*outData_d == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdCOMPLEX:
      outData_c = outData;
      *outData_c = (tcdComplex *)calloc( nTotal , sizeof( tcdComplex ));
      if (*outData_c == NULL) return( tcdERROR_ALLOC );
      break;

    case tcdDCOMPLEX:
      outData_D = outData;
      *outData_D = (tcdDComplex *)calloc( nTotal , sizeof( tcdDComplex ));
      if (*outData_D == NULL) return( tcdERROR_ALLOC );
      break;

    default:
      return( tcdERROR_UNSUPORTTYPE );

    } /* end switch outType */



  /* determine input */

  switch ( inType )
    {
    case tcdBYTE:
      inData_b = inData; break;
    case tcdSHORT:
      inData_s = inData; break;
    case tcdLONG:
      inData_l = inData; break;
    case tcdFLOAT:
      inData_f = inData; break;
    case tcdDOUBLE:
      inData_d = inData; break;
    case tcdCOMPLEX:
      inData_c = inData; break;
    case tcdDCOMPLEX:
      inData_D = inData; break;
    default:
      return( tcdERROR_UNSUPORTTYPE );
    }


  /* now cast data */


  for ( ii=0; ii<nTotal; ii++ )
    {

      switch ( outType )
	{

	case tcdSHORT:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_s)[ii] = (short )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_s)[ii] = (short )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_s)[ii] = (short )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_s)[ii] = (short )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_s)[ii] = (short )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_s)[ii] = (short )inData_c[ii].r; break;
	    case tcdDCOMPLEX:
	      (*outData_s)[ii] = (short )inData_D[ii].r; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output SHORT */
	  break;

	case tcdBYTE:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_b)[ii] = (char )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_b)[ii] = (char )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_b)[ii] = (char )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_b)[ii] = (char )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_b)[ii] = (char )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_b)[ii] = (char )inData_c[ii].r; break;
	    case tcdDCOMPLEX:
	      (*outData_b)[ii] = (char )inData_D[ii].r; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output SHORT */
	  break;

	case tcdLONG:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_l)[ii] = (long )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_l)[ii] = (long )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_l)[ii] = (long )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_l)[ii] = (long )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_l)[ii] = (long )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_l)[ii] = (long )inData_c[ii].r; break;
	    case tcdDCOMPLEX:
	      (*outData_l)[ii] = (long )inData_D[ii].r; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output LONG */
	  break;


	case tcdFLOAT:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_f)[ii] = (float )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_f)[ii] = (float )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_f)[ii] = (float )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_f)[ii] = (float )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_f)[ii] = (float )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_f)[ii] = (float )inData_c[ii].r; break;
	    case tcdDCOMPLEX:
	      (*outData_f)[ii] = (float )inData_D[ii].r; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output FLOAT */
	  break;


	case tcdDOUBLE:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_d)[ii] = (double )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_d)[ii] = (double )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_d)[ii] = (double )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_d)[ii] = (double )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_d)[ii] = (double )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_d)[ii] = (double )inData_c[ii].r; break;
	    case tcdDCOMPLEX:
	      (*outData_d)[ii] = (double )inData_D[ii].r; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output DOUBLE */
	  break;


	case tcdCOMPLEX:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_c)[ii].r = (float )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_c)[ii].r = (float )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_c)[ii].r = (float )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_c)[ii].r = (float )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_c)[ii].r = (float )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_c)[ii].r = (float )inData_c[ii].r;
	      (*outData_c)[ii].i = (float )inData_c[ii].i; break;
	    case tcdDCOMPLEX:
	      (*outData_c)[ii].r = (float )inData_D[ii].r;
	      (*outData_c)[ii].i = (float )inData_D[ii].i; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output COMPLEX*/
	  break;



	case tcdDCOMPLEX:
	  switch( inType )
	    {
	    case tcdBYTE:
	      (*outData_D)[ii].r = (double )inData_b[ii]; break;
	    case tcdSHORT:
	      (*outData_D)[ii].r = (double )inData_s[ii]; break;
	    case tcdLONG:
	      (*outData_D)[ii].r = (double )inData_l[ii]; break;
	    case tcdFLOAT:
	      (*outData_D)[ii].r = (double )inData_f[ii]; break;
	    case tcdDOUBLE:
	      (*outData_D)[ii].r = (double )inData_d[ii]; break;
	    case tcdCOMPLEX:
	      (*outData_D)[ii].r = (double )inData_c[ii].r;
	      (*outData_D)[ii].i = (double )inData_c[ii].i; break;
	    case tcdDCOMPLEX:
	      (*outData_D)[ii].r = (double )inData_D[ii].r;
	      (*outData_D)[ii].i = (double )inData_D[ii].i; break;
	    default:
	      return(tcdERROR_UNSUPORTTYPE);
	    } /* end output DOUBLE*/
	  break;
	  
	default:
	  return(tcdERROR_UNSUPORTTYPE);
	} /* end switch outType */


    } /* end for ii */


  return(tcdSUCCESS);


}



/*
  +------------------------------------------------
  +
  + Take two arrays of floats and merge into a complex
  + data array (essentally casting data).  Memory
  + for output array is already allocated by user.
  +
  +------------------------------------------------
  */  
int tcdCastToComplex(
		      tcdDATATYPE dtype,/*i: input data type               */
		      void    *real,    /*i: pointer to real part of array */
		      void    *imag,    /*i: pointer to imag part of array */
		      long      nAxes,  /*i: number of axes                */
		      long     *lAxes,  /*i: length of axes                */
		      tcdComplex *data  /*o: complex data returned         */
		      )
{
  long ii;  /* loop variable */
  int status;

  long length;

  char   *real_b, *imag_b; /* pointers to allocate data types */
  short  *real_s, *imag_s;
  long   *real_l, *imag_l;
  float  *real_f, *imag_f;
  double *real_d, *imag_d;
  
  length = 1;

  real_b = imag_b = NULL;
  real_s = imag_s = NULL;
  real_l = imag_l = NULL;
  real_f = imag_f = NULL;
  real_d = imag_d = NULL;


  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  /* User must allocate memory for data */
  if ( data == NULL )  return( tcdERROR_NULLPTR );


  switch ( dtype )
    {
    case tcdBYTE:
      real_b = real;
      imag_b = imag;
      break;
      
    case tcdSHORT:
      real_s = real;
      imag_s = imag;
      break;

    case tcdLONG:
      real_l = real;
      imag_l = imag;
      break;

    case tcdFLOAT:
      real_f = real;
      imag_f = imag;
      break;

    case tcdDOUBLE:
      real_d = real;
      imag_d = imag;
      break;

    default:
      return( tcdERROR_UNSUPORTTYPE );
    }


  for ( ii=0; ii<nAxes; ii++ ) length *= lAxes[ii];

  for (ii=0; ii<length; ii++ )
    {

      switch ( dtype )
	{
	case tcdBYTE:
	  data[ii].r = (real_b) ? (float )real_b[ii] : 0.0;
	  data[ii].i = (imag_b) ? (float )imag_b[ii] : 0.0;
	  break;


	case tcdSHORT:
	  data[ii].r = (real_s) ? (float )real_s[ii] : 0.0;
	  data[ii].i = (imag_s) ? (float )imag_s[ii] : 0.0;
	  break;

	case tcdLONG:
	  data[ii].r = (real_l) ? (float )real_l[ii] : 0.0;
	  data[ii].i = (imag_l) ? (float )imag_l[ii] : 0.0;
	  break;

	case tcdFLOAT:
	  data[ii].r = (real_f) ? (float )real_f[ii] : 0.0;
	  data[ii].i = (imag_f) ? (float )imag_f[ii] : 0.0;
	  break;

	case tcdDOUBLE:
	  data[ii].r = (real_d) ? (float )real_d[ii] : 0.0;
	  data[ii].i = (imag_d) ? (float )imag_d[ii] : 0.0;
	  break;

	default:
	  return( tcdERROR_UNSUPORTTYPE );
	}


    } /* end for ii */

  return( tcdSUCCESS );

}



/*
  +------------------------------------------------
  +
  + Take two arrays of floats and merge into a complex
  + data array (essentally casting data).  Memory
  + for output array is already allocated by user.
  +
  +------------------------------------------------
  */  
int tcdCastToDComplex(
		      tcdDATATYPE dtype,/*i: input data type               */
		      void    *real,    /*i: pointer to real part of array */
		      void    *imag,    /*i: pointer to imag part of array */
		      long      nAxes,  /*i: number of axes                */
		      long     *lAxes,  /*i: length of axes                */
		      tcdDComplex *data  /*o: complex data returned         */
		      )
{
  long ii;  /* loop variable */
  int status;

  long length;

  char   *real_b, *imag_b; /* pointers to allocate data types */
  short  *real_s, *imag_s;
  long   *real_l, *imag_l;
  float  *real_f, *imag_f;
  double *real_d, *imag_d;
  
  real_b = imag_b = NULL;
  real_s = imag_s = NULL;
  real_l = imag_l = NULL;
  real_f = imag_f = NULL;
  real_d = imag_d = NULL;

  length = 1;

  status = tcdCheckAxes( nAxes, lAxes );
  if ( status != tcdSUCCESS ) return( status );

  /* User must allocate memory for data */
  if ( data == NULL )  return( tcdERROR_NULLPTR );


  switch ( dtype )
    {
    case tcdBYTE:
      real_b = real;
      imag_b = imag;
      break;
      
    case tcdSHORT:
      real_s = real;
      imag_s = imag;
      break;

    case tcdLONG:
      real_l = real;
      imag_l = imag;
      break;

    case tcdFLOAT:
      real_f = real;
      imag_f = imag;
      break;

    case tcdDOUBLE:
      real_d = real;
      imag_d = imag;
      break;

    default:
      return( tcdERROR_UNSUPORTTYPE );
    }


  for ( ii=0; ii<nAxes; ii++ ) length *= lAxes[ii];

  for (ii=0; ii<length; ii++ )
    {

      switch ( dtype )
	{
	case tcdBYTE:
	  data[ii].r = (real_b) ? (double )real_b[ii] : 0.0;
	  data[ii].i = (imag_b) ? (double )imag_b[ii] : 0.0;
	  break;


	case tcdSHORT:
	  data[ii].r = (real_s) ? (double )real_s[ii] : 0.0;
	  data[ii].i = (imag_s) ? (double )imag_s[ii] : 0.0;
	  break;

	case tcdLONG:
	  data[ii].r = (real_l) ? (double )real_l[ii] : 0.0;
	  data[ii].i = (imag_l) ? (double )imag_l[ii] : 0.0;
	  break;

	case tcdFLOAT:
	  data[ii].r = (real_f) ? (double )real_f[ii] : 0.0;
	  data[ii].i = (imag_f) ? (double )imag_f[ii] : 0.0;
	  break;

	case tcdDOUBLE:
	  data[ii].r = (real_d) ? (double )real_d[ii] : 0.0;
	  data[ii].i = (imag_d) ? (double )imag_d[ii] : 0.0;
	  break;

	default:
	  return( tcdERROR_UNSUPORTTYPE );
	}


    } /* end for ii */

  return( tcdSUCCESS );

}

