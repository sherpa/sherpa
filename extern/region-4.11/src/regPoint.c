#include "region_priv.h"

regShape* regCopyPoint( regShape * );
int regIsEqualPoint( regShape *, regShape * );
double regCalcAreaPoint( regShape * );
int regCalcExtentPoint( regShape *, double *, double * );
int regInsidePoint( regShape *, double, double );
void regToStringPoint( regShape *, char*, long );



regShape* regCreatePoint(regFlavor include, 
			 double *xpos, 
			 double *ypos,
			 int wcoord,
             int wsize) {

  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreatePoint() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  regShape* newShape;
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Point";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( 1, sizeof(double) );
  newShape->ypos = ( double *)calloc( 1, sizeof(double) );

  newShape->xpos[0] = xpos[0]; 
  newShape->ypos[0] = ypos[0];

  newShape->type = regPOINT;
  newShape->include = include;
  newShape->nPoints = 1;

  newShape->angle = NULL;
  newShape->radius = NULL;
  newShape->sin_theta = NULL;
  newShape->cos_theta = NULL;

  newShape->calcArea   = regCalcAreaPoint;
  newShape->calcExtent = regCalcExtentPoint;
  newShape->copy       = regCopyPoint;
  newShape->isEqual    = regIsEqualPoint;
  newShape->isInside   = regInsidePoint;
  newShape->toString   = regToStringPoint;
  newShape->free       = NULL;

  newShape->region = NULL;
  newShape->next   = NULL;
  
  return newShape;
}

// copy
regShape* regCopyPoint( regShape* inShape ) {

  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyPoint() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regPOINT) {
    fprintf( stderr, "ERROR: regCopyPoint() incorrect regShape type");
    return NULL;
  }

  return regCreatePoint( inShape->include, inShape->xpos, 
			 inShape->ypos, inShape->flag_coord, inShape->flag_radius );

}

// equals
int regIsEqualPoint( regShape* thisShape, regShape* otherShape ) {

  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }
  
  if ( (thisShape->type != regPOINT) ) { 
    fprintf( stderr, "ERROR: regIsEqualPoint() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regPOINT ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualPoint() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) )
    return 0;
       
  return 1;
}


// calcArea
double regCalcAreaPoint( regShape* shape ) {
  return 0;  // always zero
}


// calcExtent
int regCalcExtentPoint( regShape* shape, double* xpos, double* ypos ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regCalcExtentPoint() requires a regShape as input");
    return 0;
  }

  if (shape->type != regPOINT) {
    fprintf( stderr, "ERROR: regCalcExtentPoint() incorrect regShape type");
    return 0;
  }
  
  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCalcExtentPoint() requires pre-allocated memory for xpos and ypos");
    return 0;
  }


  xpos[0] = shape->xpos[0];
  xpos[1] = shape->xpos[0];
  ypos[0] = shape->ypos[0];
  ypos[1] = shape->ypos[0];
  return 0;
}

// inside
int regInsidePoint( regShape* shape, double xpos, double ypos ) {
  int retval;
  
  if (shape == NULL) {
    fprintf( stderr, "ERROR: regInsidePoint() requires a regShape as input");
    return 0;
  }

  if (shape->type != regPOINT) {
      fprintf( stderr, "ERROR: regInsidePoint() incorrect regShape type");
    return 0;
  }

  if (( xpos == shape->xpos[0] ) && ( ypos == shape->ypos[0] ))
    retval = 1;
  else
    retval = 0;

  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );


}

// toString
void regToStringPoint( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringPoint() requires a regShape as input");
    return;
  }

  if (shape->type != regPOINT) {
    fprintf( stderr, "ERROR: regToStringPoint() incorrect regShape type");
    return;
  }

  if (shape->include == regExclude) {
    *buf++ = '!';
  }

  long bsize = 80;
  char *x = calloc(bsize, sizeof(char));
  char *y = calloc(bsize, sizeof(char));
  reg_print_pos_pair(shape->xpos[0], shape->ypos[0], 
                     shape->flag_coord,
                     x, y);

  snprintf(buf, maxlen, "Point(%s,%s)", x, y);

  free(x);
  free(y);
}
