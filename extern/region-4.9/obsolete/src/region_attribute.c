/*                                                                
**  Copyright (C) 2007  Smithsonian Astrophysical Observatory 
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

#include "region_priv.h"

void regFreeShape( regRegion* region, regShape* atShape )
{
  if( !atShape )
    return;
  
  if(atShape->user)
    atShape->user->free( atShape);
  else {
    if ( atShape->xpos ) free( atShape->xpos );
    if ( atShape->ypos ) free( atShape->ypos );
    if ( atShape->angle ) free( atShape->angle );
    if ( atShape->radius ) free( atShape->radius );
    if ( atShape->sin_theta) free( atShape->sin_theta);
    if ( atShape->cos_theta) free( atShape->cos_theta);
  }
  free( atShape );    
  atShape=NULL;

}


/* Attribute stuff */

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
    }
 free( region );
}

