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

#include <math.h>
#include <float.h>

#include "region_priv.h"


int regInsideRegion(
		    regRegion *region, 
		    double xpos,
		    double ypos
		    )
{
    regShape *atShape;
    int retval=0;
    int tmpval;
    int state;
  
    if ( !region ) {
        return 0;
    }

    if ( ( xpos < region->xregbounds[0] ) || ( xpos > region->xregbounds[1] ) ||
         ( ypos < region->yregbounds[0] ) || ( ypos > region->yregbounds[1] ) )
    {
        return(0);
    }
  
    atShape = region->shape;
  
    while( atShape != NULL) {
      
        tmpval = 1;

        do {
	        tmpval &= atShape->isInside( atShape, xpos, ypos );

	        state = 1;
	        if ( atShape->next == NULL ) {
	            state = 0;
	        }
	        else if ( atShape->next->component != atShape->component ) {
	            state = 0;
	        }

	        atShape = atShape->next;
        } while ( state );

        retval |= tmpval;
    }
  
    return( retval );
}


/* --------------------------------------------------------------------------- */
/* regCompareRegion                                                            */
/*    Return 1 if regions are identical                                        */
/* --------------------------------------------------------------------------- */
int regCompareRegion( regRegion* Region1, regRegion* Region2 )
{
    regShape* Shape1;
    regShape* Shape2;
    int true  = 1;
    int false = 0;
    Shape1 = Region1->shape;
    Shape2 = Region2->shape;

    while (Shape1 != NULL )
    {
        if ( !Shape2 )
            return false;

        if ( Shape1->component != Shape2->component )
            return false;

        if( !Shape1->isEqual(Shape1, Shape2) )
            return false;

        Shape1 = Shape1->next;
        Shape2 = Shape2->next;
    }

    if ( Shape2 )
        return false;

    return true;
}


/*
 *  return the maximum number of points of any
 *  shape in the region. Usually 2, but if polygons
 *  are present may be arb.
 */
long regGetMaxPoints(const regRegion * region)
{
    regShape *Shape;
    long n = 0;
    
    if (!region)
	    return 0;

    Shape = region->shape;
    while (Shape) {
	    if (Shape->nPoints > n) {
	        n = Shape->nPoints;
        }

	    Shape = Shape->next;
    }
    return n;
}


regShape *regGetShapeNo(const regRegion * region, long shapeNo)
{
    regShape *Shape;
    long no = 1;

    if (!region)
	    return NULL;
    
    Shape = region->shape;
    while (no < shapeNo) {
	    if (!Shape)
	        return NULL;
	    Shape = Shape->next;
	    no++;
    }
    return Shape;
}


long regGetNoShapes(const regRegion * region)
{
    regShape *Shape;
    long n = 0;
    
    if (!region)
	    return 0;

    Shape = region->shape;
    while (Shape) {
	    n++;
	    Shape = Shape->next;
    }
    return n;
}

/**
 * Return 1 if the regions bounding rectangles overlap. 0 otherwise.
 */
int regOverlapRegion(regRegion* region1, regRegion* region2) {
  
  // Null regions do not overlap
  if (!region1 || !region2) {
    return 0;
  }

  return reg_rectangle_overlap(
          region1->xregbounds,
          region1->yregbounds,
          region2->xregbounds,
          region2->yregbounds);
}

int reg_compare_shape( regShape* Shape1, regShape* Shape2 )
{
    // NEW CODE
    if ((!Shape1) && (!Shape2)) {
        return 1;
    }

    if ((!Shape1) || (!Shape2)) {
        return 0;
    }

    return Shape1->isEqual(Shape1, Shape2);
}

