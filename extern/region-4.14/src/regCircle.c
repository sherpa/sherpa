#include "region_priv.h"

regShape* regCopyCircle( regShape * );
int regIsEqualCircle( regShape *, regShape * );
double regCalcAreaCircle( regShape * );
int regCalcExtentCircle( regShape *, double *, double * );
int regInsideCircle( regShape *, double, double );
void regToStringCircle( regShape *, char*, long );



regShape* regCreateCircle(regFlavor include, 
			  double *xpos, 
			  double *ypos,
			  double *radius, 
			  int wcoord,
              int wsize) {

  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreateCircle() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  if (radius == NULL) {
     fprintf( stderr, "ERROR: regCreateCircle() requires an input radius.");
     return NULL;
  }

  if ( radius[0] < 0 ) {
    fprintf( stderr, "ERROR: regCreateCircle() radius of circle must be positive\n");
    return(NULL);
  }

  regShape* newShape;
   
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Circle";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( 1, sizeof(double) );
  newShape->ypos = ( double *)calloc( 1, sizeof(double) );

  newShape->xpos[0] = xpos[0]; 
  newShape->ypos[0] = ypos[0];

  newShape->type = regCIRCLE;
  newShape->include = include;
  newShape->nPoints = 1;

  newShape->radius = ( double *)calloc(1,sizeof(double));
  newShape->radius[0] = radius[0];

  newShape->angle = NULL;
  newShape->sin_theta = NULL;
  newShape->cos_theta = NULL;

  newShape->calcArea   = regCalcAreaCircle;
  newShape->calcExtent = regCalcExtentCircle;
  newShape->copy       = regCopyCircle;
  newShape->isEqual    = regIsEqualCircle;
  newShape->isInside   = regInsideCircle;
  newShape->toString   = regToStringCircle;
  newShape->free       = NULL;

  newShape->region = NULL;
  newShape->next   = NULL;

  return newShape;
}


// copy
regShape* regCopyCircle( regShape* inShape ) {
  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyCircle() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regCIRCLE) {
    fprintf( stderr, "ERROR: regCopyCircle() incorrect regShape type");
    return NULL;
  }

  return regCreateCircle( inShape->include, inShape->xpos, 
			  inShape->ypos, inShape->radius,
			  inShape->flag_coord, inShape->flag_radius );
}


// equals
int regIsEqualCircle( regShape* thisShape, regShape* otherShape ) {

  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }
  
  if ( (thisShape->type != regCIRCLE) ) { 
    fprintf( stderr, "ERROR: regIsEqualCircle() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regCIRCLE ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ||
       thisShape->radius == NULL ||  otherShape->radius == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualCircle() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->radius[0] != otherShape->radius[0]) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) )
    return 0;
       
  return 1;
}


// calcArea
double regCalcAreaCircle( regShape* shape ) {
    if (shape == NULL) {
      fprintf( stderr, "ERROR: regCalcAreaCircle() requires a regShape as input");
      return 0;
    }

    if (shape->type != regCIRCLE) {
      fprintf( stderr, "ERROR: regCalcAreaCircle() incorrect regShape type");
      return 0;
    }
  
    double area;
    area = PI * shape->radius[0] * shape->radius[0];
    return area;
}


// calcExtent
int regCalcExtentCircle( regShape* shape, double* xpos, double* ypos ) {
  if (shape == NULL) {
    fprintf( stderr, "ERROR: regCalcExtentCircle() requires a regShape as input");
    return 0;
  }

  if (shape->type != regCIRCLE) {
    fprintf( stderr, "ERROR: regCalcExtentCircle() incorrect regShape type");
    return 0;
  }
  
  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCalcExtentCircle() requires pre-allocated memory for xpos and ypos");
    return 0;
  }

  xpos[0] = shape->xpos[0] - shape->radius[0];
  xpos[1] = shape->xpos[0] + shape->radius[0];
  ypos[0] = shape->ypos[0] - shape->radius[0];
  ypos[1] = shape->ypos[0] + shape->radius[0];

  return 0;
}

// inside
int regInsideCircle( regShape* shape, double xpos, double ypos ) {
 
  if (shape == NULL) {
    fprintf( stderr, "ERROR: regInsideCircle() requires a regShape as input");
    return 0;
  }

  if (shape->type != regCIRCLE) {
      fprintf( stderr, "ERROR: regInsideCircle() incorrect regShape type");
    return 0;
  }

  int retval;
  double d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		   ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
  
  retval = (d <= shape->radius[0]) ? 1 : 0;

  if ( shape->include == regInclude ) 
    return( retval );
  else
    return( 1 - retval );

}

// toString
void regToStringCircle( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringCircle() requires a regShape as input");
    return;
  }

  if (shape->type != regCIRCLE) {
    fprintf( stderr, "ERROR: regToStringCircle() incorrect regShape type");
    return;
  }

  if (shape->include == regExclude) {
    *buf++ = '!';
  }
    
  long bsize = 80;

  char *x = calloc(bsize, sizeof(char));
  char *y = calloc(bsize, sizeof(char));
  reg_print_pos_pair(shape->xpos[0], shape->ypos[0], 
            shape->flag_coord, x, y);
  
  char *r = calloc(bsize, sizeof(char));
  reg_print_radius(shape->radius[0], shape->flag_radius, r);
  
  snprintf(buf, maxlen, "Circle(%s,%s,%s)", x, y, r);

  free(x);
  free(y);
  free(r);
}
