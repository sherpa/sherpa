/*                                                                
**  Copyright (C) -2025 2026  Smithsonian Astrophysical Observatory 
*/      


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
/*                                                                
**  Copyright (C) -2025 2026  Smithsonian Astrophysical Observatory 
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

/* calcArea
 * this is currently calculated in the regBoundsArea function and
 * passed through to area calculation. A field shape has no knowledge 
 * the field it is in, so in some sense this function is unecessary.
 * That being said, it should* return field_area as defined in
 * region_extent.
 */
double regCalcAreaField( regShape* shape ) {
    (void)shape;  /* suppress 'warning: unused parameter' */
    return 0;
}

/* calcExtent
 * Does nothing for fields
 */
int regCalcExtentField( regShape* shape, double* xpos, double* ypos ) {
    (void)xpos;   /* suppress 'warning: unused parameter' */
    (void)ypos;   /* suppress 'warning: unused parameter' */
    (void)shape;  /* suppress 'warning: unused parameter' */

  return 1;
}

/* inside */
int regInsideField( regShape* shape, double x, double y ) {
    (void)x;  /* suppress 'warning: unused parameter' */
    (void)y;  /* suppress 'warning: unused parameter' */

    if (shape->include == regInclude) {
        return 1;
    }

    return 0;
}

/* toString */
void regToStringField( regShape* shape, char* ptr, long maxlength) {
    if (!shape) return;

    if (shape->include == regExclude) {
        *ptr++ = '!';
    }

    snprintf(ptr, maxlength, "Field()");
}
