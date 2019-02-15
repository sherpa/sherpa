/*
 * Includes operations for regSector type shapes.
 */
#include <math.h>
#include <float.h>
#include "region_priv.h"

regShape* regCopySector( regShape * );
int regIsEqualSector( regShape *, regShape * );
double regCalcAreaSector( regShape * );
int regCalcExtentSector( regShape *, double *, double * );
int regInsideSector( regShape *, double, double );
void regToStringSector( regShape *, char *, long );

regShape* regCreateSector(regFlavor include,
                          double *xpos, 
                          double *ypos,
                          double *angle,
                          int wcoord,
                          int wsize)
{
    // Input validation
    if (!xpos || !ypos || !angle) {
        fprintf( stderr, "ERROR: Null input for regCreateSector\n");
        return (NULL);
    }

    double theta;

    // create new shape
    regShape *newShape;
    newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
    
    // Shape type
    newShape->type = regSECTOR;
    newShape->name = "Sector";

    // World coords, point number, and inclusion
    newShape->include = include;
    newShape->nPoints = 1;
    newShape->flag_coord = wcoord;
    newShape->flag_radius = wsize;

    // Allocate memory for point, radii, angles
    newShape->xpos = (double *) calloc(1, sizeof(double));
    newShape->ypos = (double *) calloc(1, sizeof(double));
    newShape->angle = (double *) calloc(2,sizeof(double));
    newShape->sin_theta = (double *) calloc(1,sizeof(double));
    newShape->cos_theta = (double *) calloc(1,sizeof(double));
    
    // Fill in values
    newShape->xpos[0] = xpos[0];
    newShape->ypos[0] = ypos[0];

    newShape->angle[0] = angle[0];
    newShape->angle[1] = angle[1];

    theta = newShape->angle[0]*PI/180.0;
    *newShape->sin_theta = sin(theta);
    *newShape->cos_theta = cos(theta);

    newShape->radius = NULL;
    
    // Add relevant methods
    newShape->calcArea = regCalcAreaSector;
    newShape->calcExtent = regCalcExtentSector;
    newShape->copy = regCopySector;
    newShape->isEqual = regIsEqualSector;
    newShape->isInside = regInsideSector;
    newShape->toString = regToStringSector;
    newShape->free     = NULL;

    newShape->region = NULL;
    newShape->next   = NULL;

    return newShape;
}

// copy
regShape* regCopySector( regShape* shape ) {
    if (shape->type != regSECTOR) {
	    fprintf( stderr, "ERROR: Attempting to copy %s as a Sector\n", shape->name);
	    return(NULL);
    }

    return regCreateSector(shape->include,
                           shape->xpos,
                           shape->ypos,
                           shape->angle,
                           shape->flag_coord,
                           shape->flag_radius);
}

// equals
int regIsEqualSector( regShape* thisShape, regShape* otherShape ) {
    if ( !thisShape && !otherShape ) {
        return 1;
    }

    if ( !thisShape || !otherShape) {
        return 0;
    }
    
    if (thisShape->type != regSECTOR) {
	    fprintf( stderr, "ERROR: not comparing a Sector shape\n");
    }

    if (otherShape->type != regSECTOR) {
        return 0;
    }

    if (thisShape->include != otherShape->include ||
        thisShape->xpos[0] != otherShape->xpos[0] ||
        thisShape->ypos[0] != otherShape->ypos[0] ||
        thisShape->angle[0] != otherShape->angle[0] ||
        thisShape->angle[1] != otherShape->angle[1])
    {
        return 0;
    }

    return 1;
}

// calcArea
double regCalcAreaSector( regShape* shape ) {
    // This is meaningless without information about the field. AFAICT we
    // are always returning 0 now, so maintaining that behavior.
    return 1.0;
}

// calcExtent
int regCalcExtentSector( regShape* shape, double* xpos, double* ypos ) {
    // Nothing needs to be calculated here.  Maintaining previous behavior.
    return 1;
}

// inside
int regInsideSector( regShape* shape, double xpos, double ypos ) {
  int retval;
  double d;
  double ang1 = reg_mod_angle(shape->angle[0]);
  double ang2 = reg_mod_angle(shape->angle[1]);
  double angat;

  angat = atan2( (ypos - shape->ypos[0]) , ( xpos  - shape->xpos[0] ) );
  angat = reg_mod_angle( angat * 180.0/PI );

  if ( ang1 < ang2 )
    {
      if ( ( angat >= ang1 ) && ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }
  else
    {
      if ( ( angat >= ang1 ) || ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }


  if (( retval ) && ( shape->radius) ) /* sectors have NULL radius */
  {
   d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		  ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
 	
 
   retval = (d <= shape->radius[1] && d >= shape->radius[0] ) ? 1 : 0;
  } 

  if ( ( xpos == shape->xpos[0] ) && ( ypos == shape->ypos[0] )) 
    {
      if (shape->radius) 
	if ( shape->radius[0] == 0.0 ) retval=1;
    }
  
  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );


}

// toString
void regToStringSector( regShape* shape, char * ptr, long maxlength ) {

    if (!shape) return;
    
    if (shape->include == regExclude) {
        *ptr++ = '!';
    }

    long bsize = 80;

    char *x = calloc(bsize, sizeof(char));
    char *y = calloc(bsize, sizeof(char));
    reg_print_pos_pair(shape->xpos[0], shape->ypos[0], 
            shape->flag_coord, x, y);
  
    char *a1 = calloc(bsize, sizeof(char));
    char *a2 = calloc(bsize, sizeof(char));
    reg_print_val(shape->angle[0], a1);
    reg_print_val(shape->angle[1], a2);
  
    snprintf(ptr, maxlength, "Sector(%s,%s,%s,%s)", x, y, a1, a2);
  
    free(x);
    free(y);
    free(a1);
    free(a2);
}

