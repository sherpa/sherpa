/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
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

