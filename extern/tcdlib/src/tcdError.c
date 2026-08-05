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
