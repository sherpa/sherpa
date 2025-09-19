/*
 * Includes operations for regAnnulus type shapes.
 */
#include "region_priv.h"

regShape* regCopyAnnulus( regShape * );
int regIsEqualAnnulus( regShape *, regShape * );
double regCalcAreaAnnulus( regShape * );
int regCalcExtentAnnulus( regShape *, double *, double * );
int regInsideAnnulus( regShape *, double, double );
void regToStringAnnulus( regShape *, char *, long );

regShape* regCreateAnnulus(regFlavor include,
                           double *xpos, 
                           double *ypos,
                           double *radius,
                           int wcoord,
                           int wsize)
{
    // Input validation
    if ( !xpos || !ypos || !radius ) {
        fprintf( stderr, "ERROR: Null input for regCreateAnnulus\n");
        return (NULL);
    }

    if ( radius[0] < 0 ) {
	    fprintf( stderr, "ERROR: inner radius of annulus must be positive\n");
	    return(NULL);
    }
    if ( radius[1] < 0 ) {
	    fprintf( stderr, "ERROR: outer radius of annulus must be positive\n");
	    return(NULL);
    }
    if ( radius[1] < radius[0] ) {
	    fprintf( stderr, "ERROR: annulus outer radius must be larger than "
		         "inner radius\n");
	    return(NULL);
    } 

    // create new shape
    regShape *newShape;
    newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );

    // shape type
    newShape->type = regANNULUS;
    newShape->name = "Annulus";

    // World coords, point number, and inclusion
    newShape->include = include;
    newShape->nPoints = 1;
    newShape->flag_coord = wcoord;
    newShape->flag_radius = wsize;

    // Allocate memory for point, radii, angles
    newShape->xpos = (double *) calloc(1, sizeof(double));
    newShape->ypos = (double *) calloc(1, sizeof(double));
    newShape->radius = (double *) calloc(2, sizeof(double));
    
    // Fill in values
    newShape->angle = NULL;
    newShape->sin_theta = NULL;
    newShape->cos_theta = NULL;

    newShape->xpos[0] = xpos[0];
    newShape->ypos[0] = ypos[0];

    newShape->radius[0] = radius[0];
    newShape->radius[1] = radius[1];
    
    // Add relevant methods
    newShape->calcArea = regCalcAreaAnnulus;
    newShape->calcExtent = regCalcExtentAnnulus;
    newShape->copy = regCopyAnnulus;
    newShape->isEqual = regIsEqualAnnulus;
    newShape->isInside = regInsideAnnulus;
    newShape->toString = regToStringAnnulus;
    newShape->free     = NULL;

    newShape->region = NULL;
    newShape->next   = NULL;
    
    // return the shape
    return newShape;
}

// copy
regShape* regCopyAnnulus( regShape* shape ) {
    if (shape->type != regANNULUS) {
	    fprintf( stderr, "ERROR: Attempting to copy %s as an Annulus\n", shape->name);
	    return(NULL);
    }

    return regCreateAnnulus(shape->include,
                            shape->xpos,
                            shape->ypos,
                            shape->radius,
                            shape->flag_coord,
                            shape->flag_radius);
}

// equals
int regIsEqualAnnulus( regShape* thisShape, regShape* otherShape ) {
    if ( !thisShape && !otherShape ) {
        return 1;
    }

    if ( !thisShape || !otherShape) {
        return 0;
    }

    if (thisShape->type != regANNULUS) {
	    fprintf( stderr, "ERROR: not comparing an Annulus shape\n");
    }

    if (otherShape->type != regANNULUS) {
        return 0;
    }

    if (thisShape->include != otherShape->include ||
        thisShape->xpos[0] != otherShape->xpos[0] ||
        thisShape->ypos[0] != otherShape->ypos[0] ||
        thisShape->radius[0] != otherShape->radius[0] ||
        thisShape->radius[1] != otherShape->radius[1] ||
        thisShape->flag_coord != otherShape->flag_coord ||
        thisShape->flag_radius != otherShape->flag_radius)
    {
        return 0;
    }

    return 1;
}

// calcArea
double regCalcAreaAnnulus( regShape* shape ) {
    double area = 0.0;

	area = PI * (shape->radius[1] * shape->radius[1] -
		         shape->radius[0] * shape->radius[0]);

    return area;
}

// calcExtent
int regCalcExtentAnnulus( regShape* shape, double* xpos, double* ypos ) {
	
    xpos[0] = shape->xpos[0] - shape->radius[1];
	xpos[1] = shape->xpos[0] + shape->radius[1];
	ypos[0] = shape->ypos[0] - shape->radius[1];
	ypos[1] = shape->ypos[0] + shape->radius[1];

    return 1;
}

// inside
int regInsideAnnulus( regShape* shape, double xpos, double ypos ) {
    int retval;
    double d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		  ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
 	
 
    retval = (d <= shape->radius[1] && d >= shape->radius[0] ) ? 1 : 0;
    
    if ( shape->include == regInclude ) 
        return( retval );
    else
        return( 1 - retval );
}

// toString
void regToStringAnnulus( regShape* shape, char *ptr, long maxlength ) {
    if (!shape) return;
    
    if (shape->include == regExclude) {
        *ptr++ = '!';
    }

    long bsize = 80;

    char *x = calloc(bsize, sizeof(char));
    char *y = calloc(bsize, sizeof(char));
    reg_print_pos_pair(shape->xpos[0], shape->ypos[0], 
            shape->flag_coord, x, y);

    char *r1 = calloc(bsize, sizeof(char));
    char *r2 = calloc(bsize, sizeof(char));
	reg_print_radius(shape->radius[0], shape->flag_radius, r1);
	reg_print_radius(shape->radius[1], shape->flag_radius, r2);

    snprintf(ptr, maxlength, "Annulus(%s,%s,%s,%s)", x, y, r1, r2);
    
    free(x);
    free(y);
    free(r1);
    free(r2);
}

