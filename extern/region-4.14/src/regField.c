/* 
 * Includes operations for regField type shapes.
 */
#include "region_priv.h"

regShape* regCopyField( regShape * );
int regIsEqualField( regShape *, regShape * );
double regCalcAreaField( regShape * );
int regCalcExtentField( regShape *, double *, double * );
int regInsideField( regShape *, double, double );
void regToStringField( regShape *, char *, long );

regShape* regCreateField( regFlavor include, int wcoord, int wsize )
{
    // create new shape
    regShape* newShape;
    newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
    
    // Shape type
    newShape->type = regFIELD;
    newShape->name = "Field";

    // World coords, point number, and inclusion
    newShape->include = include;
    newShape->nPoints = 0;
    newShape->flag_coord = wcoord;
    newShape->flag_radius = wsize;
    
    // Fill in values
    newShape->xpos = NULL;
    newShape->ypos = NULL;
    newShape->angle = NULL;
    newShape->radius = NULL;
    newShape->sin_theta = NULL;
    newShape->cos_theta = NULL;
    
    // Add relevant methods
    newShape->calcArea = regCalcAreaField;
    newShape->calcExtent = regCalcExtentField;
    newShape->copy = regCopyField;
    newShape->isEqual = regIsEqualField;
    newShape->isInside = regInsideField;
    newShape->toString = regToStringField;
    newShape->free     = NULL;

    newShape->region = NULL;
    newShape->next   = NULL;

    return newShape;
}

// copy
regShape* regCopyField( regShape* shape ) {
    if (shape->type != regFIELD) {
	    fprintf( stderr, "ERROR: Attempting to copy %s as a Field\n", shape->name);
	    return(NULL);
    }

    return regCreateField(shape->include, 
                          shape->flag_coord,
                          shape->flag_radius);
}

// equals
int regIsEqualField( regShape* thisShape, regShape* otherShape ) {
    if ( !thisShape && !otherShape ) {
        return 1;
    }

    if ( !thisShape || !otherShape) {
        return 0;
    }

    if (thisShape->type != regFIELD) {
	    fprintf( stderr, "ERROR: not comparing a Field shape\n");
    }

    if (otherShape->type != regFIELD) {
        return 0;
    }

    if (thisShape->include != otherShape->include)
    {
        return 0;
    }

    return 1;
}

// calcArea
// this is currently calculated in the regBoundsArea function and
// passed through to area calculation. A field shape has no knowledge 
// the field it is in, so in some sense this function is unecessary.
// That being said, it should* return field_area as defined in
// region_extent.
double regCalcAreaField( regShape* shape ) {
    return 0;
}

// calcExtent
// Does nothing for fields
int regCalcExtentField( regShape* shape, double* xpos, double* ypos ) {
    return 1;
}

// inside
int regInsideField( regShape* shape, double x, double y ) {
    
    if (shape->include == regInclude) {
        return 1;
    }

    return 0;
}

// toString
void regToStringField( regShape* shape, char* ptr, long maxlength) {
    if (!shape) return;

    if (shape->include == regExclude) {
        *ptr++ = '!';
    }

    snprintf(ptr, maxlength, "Field()");
}
