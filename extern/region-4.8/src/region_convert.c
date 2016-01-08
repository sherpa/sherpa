/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "region_priv.h"
#include <float.h>


/*
 *  Convert a region which may contain world coord shapes to
 *  one with only pixel coord shapes. Requires a function to convert
 *  from world to pixel, e.g.
 *       void CoordInvert( void* data, double* world, double* pixel )
 * The first argument is provided for users to pass in control info.
 */
void regConvertWorldRegion(
    regRegion* Region,         // (i) Region which may have shapes in world units
    double scale,              // (i) degrees per pixel
    regInvertFunction invert)  // (i) Function to convert world to pixel coords
{
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };
    regShape* Shape;
    int force = 0;

    if (!Region) {
        return;
    }

    Shape = Region->shape;

    while (Shape != NULL )
    {
        reg_convert_world_shape( Shape, scale, invert, force );
        Shape = Shape->next;
    }

    regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
}


void regConvertRegion(
    regRegion* Region,         // (i) Region which may have shapes in world units 
    double scale,              // (i) degrees per pixel 
    regInvertFunction invert,  // (i) Function to convert world to pixel coords
    int force)
{
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };
    regShape* Shape;
    
    if (!Region) {
        return;
    }

    Shape = Region->shape;

    while (Shape != NULL )
    {
        reg_convert_world_shape( Shape, scale, invert, force );
        Shape = Shape->next;
    }

    regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
}


/*
 *
 */

void regResolveField( regRegion* region, double* xx, double* yy )
{
    regShape *newShape;
    if ( !region ) return;
    
    if ( region->shape->type == regFIELD ) {

        newShape = regCreateRectangle(regInclude, xx, yy, 0, RC_PHYSICAL, RC_PHYSICAL);
        newShape->component = 1;

        // Replace head of shapes list with the new rectangle.
        newShape->next = region->shape->next;
        regFreeShape(region, region->shape);
        region->shape = newShape;
        newShape->region = region;

    } else {
        return;
    }

    region->shape = newShape;

}

/*
 *  Convert a shape expressed in degrees to one in pixels.
 *  Positions are converted using the invert function
 *  Radii are converted by dividing by the CDELT scale value
 *  We support the possibility that positions are expressed in pixels
 *  but radii in angular measure, and vice versa.
 */
void reg_convert_world_shape( regShape* Shape, double scale, regInvertFunction invert, int force )
{
    long ii;
    long nradii = 0;
    long npoints = Shape->nPoints;
    double world[2];
    double pixel[2];

    if ( Shape->flag_coord == RC_WORLD || force ) {
        for ( ii=0; ii<npoints;ii++) {
            world[0] = Shape->xpos[ii];
            world[1] = Shape->ypos[ii];
            invert( world, pixel );
            // printf("\n%s (xpos,ypos)=(%10.3f, %10.3f) bef", __FUNCTION__, Shape->xpos[ii],Shape->ypos[ii]);
            Shape->xpos[ii] = pixel[0];       
            Shape->ypos[ii] = pixel[1];       
            // printf("\n%s (xpos,ypos)=(%10.3f, %10.3f) aft", __FUNCTION__, Shape->xpos[ii],Shape->ypos[ii]);
        }
	    
        // Down-cast the shape flag_coord
	    if ( Shape->flag_coord == RC_WORLD )
	        Shape->flag_coord = RC_PHYSICAL;
	    else if ( Shape->flag_coord == RC_PHYSICAL )
	        Shape->flag_coord = RC_LOGICAL;
    }

    if (Shape->flag_radius == RC_WORLD || force) {
        nradii = reg_shape_radii( Shape );
        for ( ii = 0; ii < nradii; ii++ ) {
            Shape->radius[ii] = Shape->radius[ii] / scale;
        }

        // Down-cast the shape flag_radius
	    if ( Shape->flag_radius == RC_WORLD )
	        Shape->flag_radius = RC_PHYSICAL;
	    else if ( Shape->flag_radius == RC_PHYSICAL )
	        Shape->flag_radius = RC_LOGICAL;
    }
}


/*
 *  Returns allocated short array *mask
 */
int regRegionToMask(regRegion * region,
		    double xmin,
		    double xmax,
		    double ymin,
		    double ymax,
		    double bin, short **mask, long *xlen, long *ylen)
{
    long ii, jj;

    *xlen = (xmax - xmin) / bin + 1;
    *ylen = (ymax - ymin) / bin + 1;

    *mask = (short *) calloc((*xlen) * (*ylen), sizeof(short));
    if (*mask == NULL) {
	    return (-1);
    }

    for (ii = 0; ii < *xlen; ii++) {
	    for (jj = 0; jj < *ylen; jj++) {
	        *(*mask + (ii + jj * (*xlen))) =
		    regInsideRegion(region, xmin + bin * ii, ymin + bin * jj);
	    }
    }

    return (0);
}


/*
 *  Returns allocated double arrays *xat and *yat, with the number of points
 *  in *nat.
 */
int regRegionToList(regRegion * region,
		    double xmin,
		    double xmax,
		    double ymin,
		    double ymax,
		    double bin, double **xat, double **yat, long *nat)
{
    long ii, jj;
    double xpos, ypos;
    double xlen = (xmax - xmin) / bin + 1;
    double ylen = (ymax - ymin) / bin + 1;
    int stat;
    long Buffer = 200;

    *nat = 0;
    
    if (!region) {
        return 1;
    }
    *xat = (double *) calloc(Buffer, sizeof(double));
    *yat = (double *) calloc(Buffer, sizeof(double));


    for (ii = 0; ii < xlen; ii++) {
	    xpos = xmin + bin * ii;
	    for (jj = 0; jj < ylen; jj++) {
	        ypos = ymin + bin * jj;

	        stat = regInsideRegion(region, xpos, ypos);

	        if (stat) {
                
                // If we're out of space increase the buffer size.
		        *nat += 1;
		        if ((*nat >= Buffer)) {
		            Buffer = 2 * (*nat / Buffer + 1) * Buffer;
		            *xat =
			            (double *) realloc(*xat, Buffer * sizeof(double));
		            *yat =
			            (double *) realloc(*yat, Buffer * sizeof(double));
		        }

		        *(*xat + (*nat - 1)) = xpos;
		        *(*yat + (*nat - 1)) = ypos;
	        }
        }
    }

    return 0;
}


/* regTerm structure
   Used to hold information about regions we are inverting 
*/
struct regTerm 
{
    struct regTerm *next;  // next OR term
    struct regTerm *prev;	// previous OR term
    regShape *first;	// first shape in term 
    regShape *current;	// current shape in term
    regShape *last;	// last shape in term
};


/* Invert a region

   We come here with a region structure containing a linked list of "shapes"

   Region->shapeA->shapeB->shapeC->shapeD etc

   Each shape has a flag which indicates its relation to the previous shape,
   whether AND or OR.  AND operations take precedence over OR, so we can
   view the shape list as a series of AND terms separated by OR operations:

	P*Q+R*S*T+W*Z = (P*Q)+(R*S*T)+(W*Z)

   This inversion function has two stages.  First it creates a linked
   list of terms.  A term is a set of shapes grouped by AND.  Successive
   entries in the term list are joined by OR.


   From DeMorgan's laws:

   ___________________     _ _ _   _ _ _   _ _ _   _ _ _   _ _ _   _ _ _
   (P*Q)+(R*S*T)+(W*Z)  =  P*R*W + P*R*Z + P*S*W + P*S*Z + P*T*W + P*T*Z +
                           _ _ _   _ _ _   _ _ _   _ _ _   _ _ _   _ _ _
   			   Q*R*W + Q*R*Z + Q*S*W + Q*S*Z + Q*T*W + Q*T*Z 

   Stated more simply, the left hand side is composed of all combinations
   of one shape chosen from each term, with each shape negated.  This is NOT
   the most compact way to write the results of DeMorgan's transformation,
   but is a representation which requires no explicit grouping syntax: AND 
   over OR precedence being adequate.

   The second stage of this function cycles through the term list, pulling
   shapes out, negating them, and adding them to the inverted region shape
   list.
*/
regRegion* regInvert( regRegion* inRegion )
{
    double fx[2] ={ -DBL_MAX, DBL_MAX };
    double fy[2] ={ -DBL_MAX, DBL_MAX };
    regRegion* Region;
    regShape* Shape;
    regShape* inShape;
    struct regTerm *term;
    struct regTerm *nextTerm;
    struct regTerm *firstTerm;
    struct regTerm *safeTerm;
    #ifndef TRUE
    #define TRUE 1
    #define FALSE 0
    typedef int boolean;
    #endif

    boolean done;

    // Copy shapes
    if ( !inRegion ) {
        return NULL;
    }

    Region = regCreateRegion(NULL, NULL);

    if (inRegion->shape == NULL) {
        return Region;
    }

    inShape = inRegion->shape;

    // parse the region shapes in "terms" connected by OR 
    term = (struct regTerm *) (malloc(sizeof(struct regTerm)));
    term->next = NULL;
    term->prev = NULL;
			  
    firstTerm = term;

    // initialize first term
    term->first = inShape;
    term->current = inShape;
    term->last = inShape;

    // add shapes to the current term
    while (inShape->next != NULL) 
    {
        // look for OR glue 
        if (inShape->component != inShape->next->component) {
	        // found it - finish off current term, and start new one
	        term->last = inShape;  // remember this shape as "last'

	        // malloc new term
	        nextTerm = (struct regTerm *) (malloc(sizeof(struct regTerm)));

	        // init new term
	        nextTerm->first = inShape->next;
	        nextTerm->current = inShape->next;
	        term->next = nextTerm;
	        nextTerm->prev = term;
	        nextTerm->next = NULL;

	        // new term is current term */
	        term = nextTerm;
        }
        inShape = inShape->next;
    }
    term->last = inShape;

    // now do the inversion math
    done = FALSE;
    while (!done)
    {
        // get a shape from each term 
        term = firstTerm;
        do {
	        // remember our current term (we will advance past it)
	        safeTerm = term;
	        Shape = regCopyShape(term->current );
	        Shape->include = Shape->include ? regExclude : regInclude;

	        // first shape in term is "joined" to predecessor (nothing) by OR
	        if (term == firstTerm) {
	            regAddShape( Region, regOR, Shape );
            }
	        else {
	            // all other terms join shapes by AND
	            regAddShape( Region, regAND, Shape );
            }
	 
            term = term->next;
        
        } while (term != NULL);

        // we have now got a shape from each term
        term = safeTerm;

        // advance to next shape in last term */
        if (term->current != term->last) {
            term->current = term->current->next;
        }
        else {
	        // already at last shape, reset current term to first
	        if (term != firstTerm) {
	            term->current = term->first;
            }
	        else {
	            done = TRUE;
            }

	        while (term != firstTerm)
	        {
	            // go to previous term and advance or reset as appropriate
	            term = term->prev;

	            // any shapes left?
	            if (term->current != term->last)
	            {
		            // yes advance and break out of loop 
		            term->current = term->current->next;
		            break;
	            }
	            else {
	                // no - reset to first shape...*/
	                if (term != firstTerm) {
		                term->current = term->first;
                    }
	                else {
		                // ...unless this is the first term,
		                //   in which case we are done
		                done = TRUE;
                    }
                }
	        }
       }
   }


    // dispose of our parsing structures
    term = firstTerm;
    do {
        nextTerm = term->next;
        free(term);
        term = nextTerm;
    } while (term != NULL);

    regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
 
    return Region;
}
