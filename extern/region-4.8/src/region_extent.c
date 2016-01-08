/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include <float.h>
#include <math.h>
#include "region_priv.h"

#define REGION_MAX_PIXEL_SIZE 0.5
// Externally accessed functions

// Utility functions
int reg_zero_bounds(double *, double *);
void reg_set_bounds(double *xpos, double *ypos, double *fieldx, double *fieldy);


// Computes the bounding rectangle corners of the region w.r.t. the field 
// and places the coordinates into xpos and ypos. Returns 1 if the shape is 
// completely contained within the field. 0 if the coordinates had to be 
// trimmed.
int regExtent(regRegion * region,
	          double *fieldx, 
              double *fieldy, 
              double *xpos, 
              double *ypos) 
{
    regShape *atShape;
    int retval = 1;
    int start;
    int cstart;
    int state = 0;
    
    // Shape extents
    double sxpos[2];
    double sypos[2];

    // Component extents
    double cxpos[2];
    double cypos[2];

    // The null region is taken to be the field
    if (!region) {
    	xpos[0] = fieldx[0];
	    xpos[1] = fieldx[1];
	    ypos[0] = fieldy[0];
	    ypos[1] = fieldy[1];
	    return retval;
    }

    // Set arrays to 0 and start/cstart to 1
    start = reg_zero_bounds(xpos, ypos);
    cstart = reg_zero_bounds(cxpos, cypos);

    atShape = region->shape;
    while (atShape != NULL) {
	    do {
	        reg_extent_shape(atShape, fieldx, fieldy, sxpos, sypos);
	        reg_trim_extent(cxpos, cypos, sxpos, sypos, cstart);
	        state = 1;
	        cstart = 0;
	        if (atShape->next == NULL) {
		        state = 0;
	        } 
            else if (atShape->next->component != atShape->component) {
		        state = 0;
	        }

	        atShape = atShape->next;
        
        } while (state);

        /* End of component */
	    reg_union_extent(xpos, ypos, cxpos, cypos, start);
	    start = 0;
	    cstart = reg_zero_bounds(cxpos, cypos);
    }

    /* Trim to field */
    retval = reg_trim_extent(xpos, ypos, fieldx, fieldy, 0);
    
    return (retval);
}


// Increase cx to include sx
int reg_union_extent(double *cxpos, 
                   double *cypos, 
                   double *sxpos, 
                   double *sypos,
	               int cstart)
{
    /* Here we are doing the union, so the ranges get bigger */
    int retval = 0;
    
    if (cstart || sxpos[0] < cxpos[0]) {
	    cxpos[0] = sxpos[0];
	    retval = 1;
    }
    
    if (cstart || sxpos[1] > cxpos[1]) {
	    cxpos[1] = sxpos[1];
	    retval = 1;
    }
    
    if (cstart || sypos[0] < cypos[0]) {
	    cypos[0] = sypos[0];
	    retval = 1;
    }
    
    if (cstart || sypos[1] > cypos[1]) {
	    cypos[1] = sypos[1];
	    retval = 1;
    }
    
    if (retval) {
	    if (cxpos[0] > cxpos[1]) {
	        cxpos[0] = cxpos[1];
        }
	    if (cypos[0] > cypos[1]) {
	        cypos[0] = cypos[1];
        }
    }
    return retval;
}


// Determines the extent of a shape within the boundaries of the field.
int reg_extent_shape(regShape * shape, double *fieldx, double *fieldy, double *xpos, double *ypos)
{
    int retval = 0;

    // All exclusions imply the whole field
    // except MASK which is defined to be bound by the mask field.
    if (shape->include == regExclude && shape->type != regMASK) {
        reg_set_bounds(xpos, ypos, fieldx, fieldy);  
        retval = 1;
    }
    else {
        retval = reg_extent_shape_raw(shape, fieldx, fieldy, xpos, ypos);
    }

    return retval;
}


// Determines the extent of the shape using it's internal shape function, then trims
int reg_extent_shape_raw(regShape * shape, double *fieldx, double *fieldy, double *xpos, double *ypos)
{
    /* init to field */
    reg_set_bounds(xpos, ypos, fieldx, fieldy);
    shape->calcExtent(shape, xpos, ypos);
    
    return 1;
}


// Initializes arrays to 0.
int reg_zero_bounds(double *xpos, double *ypos)
{
    int retval = 1;
    xpos[0] = 0.0;
    xpos[1] = 0.0;
    ypos[0] = 0.0;
    ypos[1] = 0.0;
    return retval;
}

// Copies field bounds into local bounds.
void reg_set_bounds(double *xpos, double *ypos, double *fieldx, double *fieldy)
{
    xpos[0] = fieldx[0];
    xpos[1] = fieldx[1];
    ypos[0] = fieldy[0];
    ypos[1] = fieldy[1];
}
