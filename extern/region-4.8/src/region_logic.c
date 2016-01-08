/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007,2013)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "cxcregion.h"
#include "region_priv.h"
#include <float.h>

int reg_rectangle_inside( double* xpos1, double* ypos1, double* xpos2, double* ypos2 );
int reg_is_rect( regShape* shape );

regRegion* regCombineRegion( regRegion* Region1, regRegion* Region2 )
{
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };
    regRegion* Region;
    regShape* Shape;
    regShape* inShape;
    regMath glue;

    long lastComponent = 1;

    /* Copy shapes */
    if ( !Region1 )
    {
        if ( !Region2 ) {
            return NULL;
        }
        return regCopyRegion( Region2 );
    }
 
    Region = regCopyRegion( Region1 );
    inShape = Region2->shape;

    while (inShape != NULL )
    {
        Shape = regCopyShape(inShape);
        if ( inShape->component == lastComponent )
        {
            glue = regAND;
        } 
        else {
            glue = regOR;
        }
        lastComponent = inShape->component;
        regAddShape( Region, glue, Shape );

        inShape = inShape->next;
    }
    regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
    return Region;
}


regRegion* regUnionRegion( regRegion* Region1, regRegion* Region2 )
{
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };

    regRegion* region;
    regShape*  shape;
    regShape*  inShape;
    regMath    glue;
  
    long lastComponent;
    int  haveMore;
  
    /* Copy shapes */
    if ( !Region1 )
    {
        if ( !Region2 ) {
            return NULL;
        }
        return regCopyRegion( Region2 );
    }
  
    /* If regions are equal just return copy of one */
    if ( regCompareRegion( Region1, Region2 )) {
         return regCopyRegion( Region1 );
    }
  
    /* Make a new region with all the combined components of the input regions */
    /*  - Put USER defined shapes first in the series.                         */
    region = regCreateRegion(NULL, NULL); 

    /* Transfer Region components with USER shapes.                            */
    /*   NOTE: expectation is that USER shapes are first on the component      */
    haveMore = 1;
    inShape = Region1->shape;
    
    while ( inShape != NULL )
    {
        if ( inShape->type == regMASK )
        {
            glue = regOR;
            lastComponent = inShape->component;
            while ( inShape && inShape->component == lastComponent )
            {
                shape = regCopyShape(inShape);
                regAddShape( region, glue, shape );
                glue = regAND;
                inShape = inShape->next;
            }
        }
        else
        {
            inShape = reg_next_component( inShape );
        }
        
        if ( (inShape == NULL) && haveMore )
        {
            /* scan second region components */
            inShape = Region2->shape;
            haveMore = 0;
        }
    }
  
    /* Transfer Region components Non-USER shapes.*/
    haveMore = 1;
    inShape = Region1->shape;
    while ( inShape != NULL )
    {
        if ( inShape->type == regMASK )
        {
            inShape = reg_next_component( inShape );
        } 
        else 
        {
            glue = regOR;
            lastComponent = inShape->component;
            while ( inShape && inShape->component == lastComponent )
            {
                shape = regCopyShape(inShape);
                regAddShape( region, glue, shape );
                glue = regAND;
                inShape = inShape->next;
            }
        }

        if ( (inShape == NULL) && haveMore )
        {
            /* scan second region components */
            inShape = Region2->shape;
            haveMore = 0;
        }
    }

    /* re-calculate the region extent */
    regExtent(region, fx, fy, region->xregbounds, region->yregbounds);

    return region;
}


regRegion* regIntersectRegion( regRegion* Region1, regRegion* Region2 )
{
    double fx[2] = { -DBL_MAX, DBL_MAX };
    double fy[2] = { -DBL_MAX, DBL_MAX };
    regRegion *Region = NULL;
    regShape  *Shape1 = NULL;
    regShape  *Shape2 = NULL;


    /* Copy shapes */
    if ( !Region1 ) {
        return regCopyRegion( Region2 );
    }
    if (!Region2 ) {
        return regCopyRegion( Region1 );
    }
    

    /* If regions are equal, just return copy of one */

    if ( regCompareRegion( Region1, Region2 )) {
        return regCopyRegion( Region1 );
    }
 
    Region = regCreateEmptyRegion();
    Shape1 = Region1->shape;


    while (Shape1 != NULL ) /* Loop over components of Region 1 */
    {
        Shape2 = Region2->shape; 


        while( Shape2 != NULL ) /* Loop over components of Region 2 */
        {
	    reg_intersect_component( Region, Shape1, Shape2 );
            Shape2 = reg_next_component( Shape2 );
        }

        Shape1 = reg_next_component( Shape1 );
    }
 
    regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);

    return Region;
}


/*
 *  Do two components possibly have nonzero overlap?
 */
int reg_intersect_component( regRegion* region, regShape* Shape1, regShape* Shape2 )
{
    long Cpt1, Cpt2;
    int ok = 1;
    regShape* shape1;
    regShape* shape2;
    regShape** nshape1;
    regShape** nshape2;
    regMath glue = regOR;
    long n1=0;
    long n2=0;
    long i;
    long j;
    int bIntersectOK = 1;
    long* index1;
    long* index2;
    long* isuser1;
    long* isuser2;

    if ( !Shape1 || !Shape2 ) {
        return 0;
    }

    Cpt1 = Shape1->component;
    Cpt2 = Shape2->component;

    // Count shapes
    shape1 = Shape1;
    shape2 = Shape2;

    while (shape1 != NULL && shape1->component == Cpt1 && ok )
    {
        ++n1;
        shape1 = shape1->next;
    }
 
    while (shape2 != NULL && shape2->component == Cpt2 && ok )
    {
        ++n2;
        shape2 = shape2->next;
    }
    
    index1  = calloc( n1, sizeof( long ));
    index2  = calloc( n2, sizeof( long ));
    isuser1 = calloc( n1, sizeof( long ));
    isuser2 = calloc( n2, sizeof( long ));
    nshape1 = calloc( n1, sizeof( regShape* ));
    nshape2 = calloc( n2, sizeof( regShape* ));
    shape1  = Shape1;
    shape2  = Shape2;

    // Make copies of Set 1 shapes.. set flag for all shapes 'ON'
    for ( i = 0; i < n1 && shape1; i++ )
    {
        nshape1[i] = regCopyShape(shape1);
        index1[i] = 1;
        isuser1[i] = 0;
        if (nshape1[i]->type==regMASK) {
            isuser1[i] = 1;
        }
        shape1 = shape1->next;
    }

    // Make copies of Set 2 shapes.. set flag for all shapes 'ON'
    for ( i = 0; i < n2 && shape2; i++ )
    {
        nshape2[i] = regCopyShape(shape2);
        index2[i] = 1;
        isuser2[i] = 0;
        if (nshape2[i]->type==regMASK) {
            isuser2[i] = 1;
        }
        shape2 = shape2->next;
    }

    // Intersect Set 1 with Set 2 - reducing Shapes as possible 
    // Turns flag 'OFF' for absorbed shapes..                   
    //   Code Change  Bug 13500  January 31, 2013   
           
    //   Stop looping if any intersection fails                 
    for ( i = 0; bIntersectOK && (i < n1); i++ )
    {
        for ( j = 0; bIntersectOK && (j < n2); j++ )
        {
            bIntersectOK = reg_shape_intersect( nshape1[i], nshape2[j], &index1[i], &index2[j] );
        }
    }

    // Now patch things together with positive regions first 
    //  only if intersection is successful                   
    if ( bIntersectOK )
    {
        // User defined shapes from Set 1 which are still 'ON'  
        // NOTE: This will do INCLUDE and EXCLUDE, which is     
        //       contrary to the other shape types.. OK?        
        for ( i = 0; i < n1; i++ )
        {
            if ( index1[i] && isuser1[i] )
            {
                regAddShape( region, glue, nshape1[i] );
                glue = regAND;
            }
        }

        // User defined shapes from Set 2 which are still 'ON'  
        // NOTE: This will do INCLUDE and EXCLUDE, which is     
        //       contrary to the other shape types.. OK?        
        for ( i = 0; i < n2; i++ )
        {
            if ( index2[i] && isuser2[i] )
            {
                regAddShape( region, glue, nshape2[i] );
                glue = regAND;
            }
        }

        // Include shapes from Set 1 which are still 'ON' 
        for ( i = 0; i < n1; i++ )
        { 
            if ( index1[i] && nshape1[i]->include==regInclude && !isuser1[i] )
            {
                regAddShape( region, glue, nshape1[i] );
                glue = regAND;
            }
        }

        // Include shapes from Set 2 which are still 'ON' 
        for ( i = 0; i < n2; i++ )
        {
            if ( index2[i] && nshape2[i]->include == regInclude && !isuser2[i] )
            {
                regAddShape( region, glue, nshape2[i] );
                glue = regAND;
            }
        }

        // Exclude shapes from Set 1 which are still 'ON' 
        for ( i = 0; i < n1; i++ )
        {
            if ( index1[i] && nshape1[i]->include != regInclude && !isuser1[i] )
            {
                if ( glue == regAND ) {
		    regAddShape( region, glue, nshape1[i] );
                }
            }
        }

        // Exclude shapes from Set 2 which are still 'ON' 
        for ( i = 0; i < n2; i++ )
        {
            if ( index2[i] && nshape2[i]->include != regInclude && !isuser2[i] )
            {
                if ( glue == regAND ) {
	            regAddShape( region, glue, nshape2[i] );
                }
            }
        }

    } // end bIntetsectOK 
 

    // clean up 
    for ( i = 0; i < n1; i++ )
    {
        if ( !index1[i] || (0 == bIntersectOK) ) {
            regFreeShape( NULL, nshape1[i] );
        }
    }
 
    for ( i = 0; i < n2; i++ )
    {
        if ( !index2[i] || (0 == bIntersectOK) ) {
            regFreeShape( NULL, nshape2[i] );
        }
    }
 
    free( index1 );
    free( index2 );
    free( nshape1 );
    free( nshape2 );
    free( isuser1 );
    free( isuser2 );

    return bIntersectOK;
}


/*
 * Determines if two shapes intersect, and whether or not one (or both) should be 
 * discarded in the interesction of those shapes. (e.g., if one rectangle entirely
 * contains the other, then there is no reason to put the contained rectangle inside
 * the other.)
 */
int reg_shape_intersect( regShape* shape1, regShape* shape2, long* index1, long* index2 )
{
    double fx[2] = { -DBL_MAX, DBL_MAX };
    double fy[2] = { -DBL_MAX, DBL_MAX };
    double xpos1[2];
    double ypos1[2];
    double xpos2[2];
    double ypos2[2];
    int ok = 1;

    if ( !shape1 || !shape2 ) {
        return 0;
    }

    reg_extent_shape_raw( shape1, fx, fy, xpos1, ypos1 );
    reg_extent_shape_raw( shape2, fx, fy, xpos2, ypos2 );

    ok = reg_rectangle_overlap( xpos1, ypos1, xpos2, ypos2 );
    if ( ok )
    {
        /* If shapes are identical, just return one of them */
        if ( regCompareShape( shape1, shape2, 0 ))
        {
            *index2 = 0;
            return ok;
        }
	else if ( shape1->include != shape2->include )
	{
	    /* The regions are exact opposites, they cancel out             
	     *  This often happens with bounding areas. */
	  if ( regCompareShape( shape1, shape2, 1 ))
	      {
		*index1 = 0;
		*index2 = 0;
		ok = 0;
		return ok;
	      }
	}

        /* Exclusions always overlap */
        if ( shape1->include != regInclude || shape2->include != regInclude ) {
            return ok;
        }
  
        /* The temporary hack to support only rectangles */
        /* If shape1 is entirely within rectangle 2, ignore rectangle 2. */
        if ( reg_is_rect( shape2 ))
        {
            if ( reg_rectangle_inside( xpos1, ypos1, xpos2, ypos2 ))
            {
                *index2 = 0;
            }
            else if ( reg_is_rect( shape1 ))
            {
                /* Replace shape1 with intersected rectangle */
                *index2 = 0;
                if ( shape2->xpos[0] > shape1->xpos[0] ) shape1->xpos[0] = shape2->xpos[0];
                if ( shape2->xpos[1] < shape1->xpos[1] ) shape1->xpos[1] = shape2->xpos[1];
                if ( shape2->ypos[0] > shape1->ypos[0] ) shape1->ypos[0] = shape2->ypos[0];
                if ( shape2->ypos[1] < shape1->ypos[1] ) shape1->ypos[1] = shape2->ypos[1];    
            }
        }
        else if (  reg_is_rect( shape1 ) && reg_rectangle_inside( xpos2, ypos2, xpos1, ypos1 ))
        {
            *index1 = 0;
        }

        return ok;
    }

    ok = 1;
    if ( shape1->include != regInclude )
    {
        if ( shape2->include != regInclude ) {
            return 1;
        }
        /* Exclude region outside included area of interest - ignore */
        *index1 = 0;
    }
    else if (  shape2->include != regInclude )
    {
        /* Exclude region outside included area of interest - ignore */
        *index2 = 0;
    }
    else
    {
        *index1 = 0;
        *index2 = 0; 
            /*
             * In this case, two included regions have no overlap. 
             * This is the case in which we throw out the whole intersection as null.
             */
        ok = 0;
    }

    return ok;
}


int reg_is_rect( regShape* shape )
{
    int ok;
    ok = ( shape->type == regRECTANGLE && shape->angle[0] == 0.0 );
    return ok;
}


int reg_rectangle_inside( double* xpos1, double* ypos1, double* xpos2, double* ypos2 )
{
    int ok
        = ( xpos1[0] >= xpos2[0] && xpos1[0] <= xpos2[1] &&
            xpos1[1] >= xpos2[0] && xpos1[1] <= xpos2[1] &&
            ypos1[0] >= ypos2[0] && ypos1[0] <= ypos2[1] &&
            ypos1[1] >= ypos2[0] && ypos1[1] <= ypos2[1] );
    return ok;
}


regShape* reg_next_component( regShape* shape )
{
    long cpt;
    if ( shape ) {
        cpt = shape->component;
    }
    while ( shape && shape->component == cpt ) {
        shape = shape->next;
    }
    return shape;
}
