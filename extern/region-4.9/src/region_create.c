/*                                                                
**  Copyright (C) 2007,2013  Smithsonian Astrophysical Observatory 
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

#include "cxcregion.h"
#include "region_priv.h"
#include <float.h>


/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
regRegion* regCreateEmptyRegion( void )
{
  return regCreateRegion( NULL,NULL );
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

regRegion* regCreateRegion( void *xcol, void*ycol )
{
  regRegion *region = ( regRegion *) calloc( 1, sizeof( regRegion ) );  
  return( region );  
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

regRegion* regCopyRegion( regRegion* inRegion )
{
  double fx[2] ={ -DBL_MAX, DBL_MAX };
  double fy[2] ={ -DBL_MAX, DBL_MAX };
  regRegion* Region;
  regShape* Shape;
  regShape* inShape;
  regMath glue;
  
  long lastComponent = 1;
  
  if ( !inRegion )
  return NULL;
  
  Region = regCreateRegion(NULL, NULL);
  
  inShape = inRegion->shape;
  
  while (inShape != NULL ) {
    Shape = regCopyShape( inShape );
    if ( inShape->component == lastComponent ) {
      glue = regAND;
    } else {
      glue = regOR;
    }

    lastComponent = inShape->component;
    regAddShape( Region, glue, Shape );
    
    inShape = inShape->next;
  }

  regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);

  return Region;
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

void regFree( regRegion* region )
{
  regShape *atShape;
  regShape* shape;
  /* Free shapes */
  
  if ( !region )
    return;
  
  atShape = region->shape;

  /* Free shape attributes */
  while (atShape != NULL )
    {
      shape = atShape;
      atShape = atShape->next;
      regFreeShape( region, shape );
      shape = NULL;
      
    }

  free( region );
  region = NULL;
}

