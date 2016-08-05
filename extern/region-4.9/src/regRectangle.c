#include "region_priv.h"


regShape* regCopyRectangle( regShape * );
int regIsEqualRectangle( regShape *, regShape * );
double regCalcAreaRectangle( regShape * );
int regCalcExtentRectangle( regShape *, double *, double * );
int regInsideRectangle( regShape *, double, double );
void regToStringRectangle( regShape *, char*, long );


regShape* regCreateRectangle(regFlavor include, 
			     double *xpos, 
			     double *ypos,
			     double *angle, 
			     int wcoord,
                 int wsize) {

  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreateRectangle() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  regShape* newShape;
  long npoints = 2;

  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Rectangle";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( npoints, sizeof(double ));
  newShape->ypos = ( double *)calloc( npoints, sizeof(double) );

  long ii;
  for ( ii=0; ii<npoints;ii++)
    {
      newShape->xpos[ii] = xpos[ii]; 
      newShape->ypos[ii] = ypos[ii];
    }

  newShape->include = include;
  newShape->nPoints = npoints;
 
  newShape->angle = (double *)calloc(1,sizeof(double ));
  newShape->sin_theta = (double *)calloc(1,sizeof(double ));
  newShape->cos_theta = (double *)calloc(1,sizeof(double ));

  newShape->angle[0] = angle ? angle[0] : 0.0;

  double theta = newShape->angle[0]*PI/180.0;
  *newShape->sin_theta = sin(theta);
  *newShape->cos_theta = cos(theta);

  newShape->type = angle ? regROTRECTANGLE : regRECTANGLE;

  newShape->radius = NULL;

  newShape->calcArea   = regCalcAreaRectangle;
  newShape->calcExtent = regCalcExtentRectangle;
  newShape->copy       = regCopyRectangle;
  newShape->isEqual    = regIsEqualRectangle;
  newShape->isInside   = regInsideRectangle;
  newShape->toString   = regToStringRectangle;

  return newShape;

}

// copy
regShape* regCopyRectangle( regShape* inShape ) {

  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyRectangle() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regRECTANGLE) {
    fprintf( stderr, "ERROR: regCopyRectangle() incorrect regShape type");
    return NULL;
  }

  return regCreateRectangle( inShape->include, inShape->xpos, 
			     inShape->ypos,
			     inShape->angle, inShape->flag_coord, inShape->flag_radius);
}


// equals
int regIsEqualRectangle( regShape* thisShape, regShape* otherShape ) {
    
  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }

  if ( (thisShape->type != regRECTANGLE) ) { 
    fprintf( stderr, "ERROR: regIsEqualRectangle() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regRECTANGLE ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ||
       thisShape->angle == NULL ||  otherShape->angle == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualRectangle() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->xpos[1] != otherShape->xpos[1]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->ypos[1] != otherShape->ypos[1]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) ||
       (thisShape->angle[0] != otherShape->angle[0]) )
    return 0;

  return 1;
}


// calcArea
double regCalcAreaRectangle( regShape* shape ) {

    if (shape == NULL) {
      fprintf( stderr, "ERROR: regCalcAreaRectangle() requires a regShape as input");
      return 0;
    }

    if (shape->type != regRECTANGLE) {
      fprintf( stderr, "ERROR: regCalcAreaRectangle() incorrect regShape type");
      return 0;
    }

    double area;
    double xr, yr;
    reg_rectangle_sides(shape, &xr, &yr);
    area = xr * yr;
    return area;

}


// calcExtent
int regCalcExtentRectangle( regShape* shape, double* xpos, double* ypos ) {
    if (shape == NULL) {
      fprintf( stderr, "ERROR: regCalcExtentRectangle() requires a regShape as input");
      return 0;
    }
    
    if (shape->type != regRECTANGLE) {
      fprintf( stderr, "ERROR: regCalcExtentRectangle() incorrect regShape type");
      return 0;
    }
    
    if ((xpos == NULL) || (ypos == NULL)) {
      fprintf( stderr, "ERROR: regCalcExtentRectangle() requires pre-allocated memory for xpos and ypos");
      return 0;
    }
    
    double xcor[4], ycor[4];
    reg_rectangle_corners(shape, xcor, ycor);
    reg_corner_bounds(xcor, ycor, xpos, ypos);

    return 0;

}


// inside
int regInsideRectangle( regShape* shape, double xpos, double ypos ) {
  
  if (shape == NULL) {
    fprintf( stderr, "ERROR: regInsideRectangle() requires a regShape as input");
    return 0;
  }

  if (shape->type != regRECTANGLE) {
    fprintf( stderr, "ERROR: regInsideRectangle() incorrect regShape type");
    return 0;
  }

  int retval;

  if ( shape->angle[0] == 0.0 )
    {
      if ( ( xpos >= shape->xpos[0] ) &&
	   ( xpos <= shape->xpos[1] ) &&
	   ( ypos >= shape->ypos[0] ) &&
	   ( ypos <= shape->ypos[1] ) )
	retval = 1;
      else
	retval = 0;
    }
  else
    {
      double xr[2], yr[2];
      double xpr, ypr;
      double xcen = ( shape->xpos[1] + shape->xpos[0] )/2.0;
      double ycen = ( shape->ypos[1] + shape->ypos[0] )/2.0;

      reg_rotated_coords( shape, xpos, ypos, xcen, ycen, &xpr, &ypr );
      reg_rotated_coords( shape, shape->xpos[0], shape->ypos[0], xcen, ycen, &xr[0], &yr[0] );
      reg_rotated_coords( shape, shape->xpos[1], shape->ypos[1], xcen, ycen, &xr[1], &yr[1] );

      if ( ( xpr >= xr[0] ) &&
	   ( xpr <= xr[1] ) &&
	   ( ypr >= yr[0] ) &&
	   ( ypr <= yr[1] ) )
	retval = 1;
      else
	retval = 0;
    }

  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}


// toString
void regToStringRectangle( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringRectangle() requires a regShape as input");
    return;
  }

  if (shape->type != regRECTANGLE) {
    fprintf( stderr, "ERROR: regToStringRectangle() incorrect regShape type");
    return;
  }

  if (shape->include == regExclude) {
    *buf++ = '!';
  }
  
  long bsize = 80;

  char *x1 = calloc(bsize, sizeof(char));
  char *y1 = calloc(bsize, sizeof(char));
  reg_print_pos_pair(shape->xpos[0], shape->ypos[0], 
            shape->flag_coord, x1, y1);

  char *x2 = calloc(bsize, sizeof(char));
  char *y2 = calloc(bsize, sizeof(char));
  reg_print_pos_pair(shape->xpos[1], shape->ypos[1], 
            shape->flag_coord, x2, y2);
  
  char *a = calloc(bsize, sizeof(char));
  reg_print_val(shape->angle[0], a);
  
  
  if (shape->angle[0] == 0.0)
    snprintf(buf, maxlen, "Rectangle(%s,%s,%s,%s)", x1, y1, x2, y2);
  else
    snprintf(buf, maxlen, "RotRectangle(%s,%s,%s,%s,%s)", x1, y1, x2, y2, a);
  
  free(x1);
  free(y1);
  free(x2);
  free(y2);
  free(a);
}
