/*
 * Includes operations for regMask type shapes.
 *
 * AFAICT regMASK internals are never accessed by utility functions in
 * the region libray. We are including the implementation of skeleton 
 * functions here, but they will all be dummy functions that spit out an
 * alert if they are actually called.
 *
 */
#define MEMCHECK(xx) if ((xx)==NULL) return(NULL)
#include "region_priv.h"

regShape* regCopyMask( regShape * );
int regIsEqualMask( regShape *, regShape * );
double regCalcAreaMask( regShape * );
int regCalcExtentMask( regShape *, double *, double * );
int regInsideMask( regShape *, double, double );
void regToStringMask( regShape *, char *, long );

regShape* regCreateMask( regFlavor include, int wcoord, int wsize )
{
    // NOTE:
    //  Only weird thing I see going on is in regConvertWorldShape. DM_FLAG_LOGICAL?
	fprintf( stderr, "ERROR: Cannot create a regMask type shape\n" ); 

    // create new shape
    regShape *newShape = NULL;
    /**
    newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
    MEMCHECK(newShape);
    
    // Shape type
    newShape->type = regMASK;
    newShape->name = "Mask";

    // World coords, point number, and inclusion
    newShape->include = include;

    // Add relevant methods
    newShape->calcArea = regCalcAreaMask;
    newShape->calcExtent = regCalcExtentMask;
    newShape->copy = regCopyMask;
    newShape->isEqual = regIsEqualMask;
    newShape->isInside = regInsideMask;
    newShape->toString = regToStringMask;
    */
    return newShape;
}

// copy
regShape* regCopyMask( regShape* shape ) {
    // Ditto, not sure what's going on with the user field in regMASK objects.
	fprintf( stderr, "ERROR: Cannot copy regMask type shape\n" );
	return(NULL);
}

// equals
int regIsEqualMask( regShape* thisShape, regShape* otherShape ) {
	fprintf( stderr, "ERROR: Cannot test equality for a regMask type shape\n" ); 
    return 1;
}

// calcArea
double regCalcAreaMask( regShape* shape ) {
	fprintf( stderr, "ERROR: Cannot calculate area for a regMask type shape\n" ); 
    return 0;
}

// calcExtent
// TODO:
//   Based on what's in region_extent.c there's probably more to do for this.
int regCalcExtentMask( regShape* shape, double* xpos, double* ypos ) {
	fprintf( stderr, "ERROR: Cannot calculate extent for a regMask type shape\n" ); 
    return 0;
}

// inside
int regInsideMask( regShape* shape, double x, double y ) {
	fprintf( stderr, "ERROR: Cannot compute inside for a regMask type shape\n" ); 
    return 0;
}

// toString
void regToStringMask( regShape* shape, char* ptr, long maxlength) {
	fprintf( stderr, "ERROR: Cannot create string for a regMask type shape\n" ); 
}
