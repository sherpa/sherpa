#include "region_priv.h"

regShape* regCopyEllipse( regShape * );
int regIsEqualEllipse( regShape *, regShape * );
double regCalcAreaEllipse( regShape * );
int regCalcExtentEllipse( regShape *, double *, double * );
int regInsideEllipse( regShape *, double, double );
void regToStringEllipse( regShape *, char*, long );


regShape* regCreateEllipse(regFlavor include, 
			   double *xpos,
			   double *ypos,
			   double *radius, 
			   double *angle, 
			   int wcoord,
               int wsize) {

  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreateEllipse() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  if (radius == NULL) {
     fprintf( stderr, "ERROR: regCreateEllipse() requires two (2) radii.");
     return NULL;
  }

  if ( radius[0] < 0 ) {
    fprintf( stderr, "ERROR: regCreateEllipse() 1st radius of ellipse must be positive\n");
    return(NULL);
  }

  if ( radius[1] < 0 ) {
    fprintf( stderr, "ERROR: regCreateEllipse() 2nd radius of ellipse must be positive\n");
    return(NULL);
  }

  regShape* newShape;
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Ellipse";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( 1, sizeof(double ));
  newShape->ypos = ( double *)calloc( 1, sizeof(double) );

  newShape->xpos[0] = xpos[0]; 
  newShape->ypos[0] = ypos[0];

  newShape->type = regELLIPSE;
  newShape->include = include;
  newShape->nPoints = 1;
  
  newShape->angle = (double *)calloc(1,sizeof(double ));
  newShape->sin_theta = (double *)calloc(1,sizeof(double ));
  newShape->cos_theta = (double *)calloc(1,sizeof(double ));
  
  newShape->angle[0] = angle ? angle[0] : 0.0;
  double theta = newShape->angle[0]*PI/180.0;
  *newShape->sin_theta = sin(theta);
  *newShape->cos_theta = cos(theta);

  newShape->radius = (double *)calloc(2,sizeof(double ));
  newShape->radius[0] = radius[0];      
  newShape->radius[1] = radius[1];
  
  newShape->calcArea   = regCalcAreaEllipse;
  newShape->calcExtent = regCalcExtentEllipse;
  newShape->copy       = regCopyEllipse;
  newShape->isEqual    = regIsEqualEllipse;
  newShape->isInside   = regInsideEllipse;
  newShape->toString   = regToStringEllipse;
  newShape->free       = NULL;

  newShape->region = NULL;
  newShape->next   = NULL;

  return newShape;
  
}


// copy
regShape* regCopyEllipse( regShape* inShape ) {

  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyEllipse() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regELLIPSE) {
    fprintf( stderr, "ERROR: regCopyEllipse() incorrect regShape type");
    return NULL;
  }

  return regCreateEllipse( inShape->include, inShape->xpos, 
			   inShape->ypos, inShape->radius,
			   inShape->angle, inShape->flag_coord, inShape->flag_radius);
}


// equals
int regIsEqualEllipse( regShape* thisShape, regShape* otherShape ) {
  
  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }
  
  if ( (thisShape->type != regELLIPSE) ) { 
    fprintf( stderr, "ERROR: regIsEqualEllipse() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regELLIPSE ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ||
       thisShape->radius == NULL ||  otherShape->radius == NULL ||
       thisShape->angle == NULL ||  otherShape->angle == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualEllipse() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->radius[0] != otherShape->radius[0]) ||
       (thisShape->radius[1] != otherShape->radius[1]) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) ||
       (thisShape->angle[0] != otherShape->angle[0]) )
    return 0;

  return 1;
}


// calcArea
double regCalcAreaEllipse( regShape* shape ) {

    if (shape == NULL) {
      fprintf( stderr, "ERROR: regCalcAreaEllipse() requires a regShape as input");
      return 0;
    }

    if (shape->type != regELLIPSE) {
      fprintf( stderr, "ERROR: regCalcAreaEllipse() incorrect regShape type");
      return 0;
    }

    double area;
    area = PI * shape->radius[0] * shape->radius[1];
    return area;
}


// calcExtent
int regCalcExtentEllipse( regShape* shape, double* xpos, double* ypos ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regCalcExtentEllipse() requires a regShape as input");
    return 0;
  }

  if (shape->type != regELLIPSE) {
    fprintf( stderr, "ERROR: regCalcExtentEllipse() incorrect regShape type");
    return 0;
  }
  
  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCalcExtentEllipse() requires pre-allocated memory for xpos and ypos");
    return 0;
  }
  
  double xcor[4], ycor[4];
  reg_box_corners(shape, xcor, ycor);
  reg_corner_bounds(xcor, ycor, xpos, ypos);
  return 0;

}


// inside
int regInsideEllipse( regShape* shape, double xpos, double ypos ) {
 
  if (shape == NULL) {
    fprintf( stderr, "ERROR: regInsideEllipse() requires a regShape as input");
    return 0;
  }

  if (shape->type != regELLIPSE) {
    fprintf( stderr, "ERROR: regInsideEllipse() incorrect regShape type");
    return 0;
  }

  int retval=0;
  double xr, yr;

  if ( shape->angle[0] == 0 )
    {
      xr = ( xpos - shape->xpos[0] ) / shape->radius[0] ;
      yr = ( ypos - shape->ypos[0] ) / shape->radius[1] ;
    }
  else
    {
      double xm, ym;
      
      /*double theta = shape->angle[0] * PI / 180.0;*/
      
      xm = xpos - shape->xpos[0];
      ym = ypos - shape->ypos[0];
      
      xr = ( *shape->cos_theta * xm + *shape->sin_theta * ym)/shape->radius[0];
      yr = (-*shape->sin_theta * xm + *shape->cos_theta * ym)/shape->radius[1];
    }  
  
  if ( ( xr * xr + yr * yr ) <= 1.0 ) retval = 1;
  
  if ( shape->include == regInclude )
    return( retval );
  else
    return( 1 - retval );
  
}


// toString
void regToStringEllipse( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringEllipse() requires a regShape as input");
    return;
  }

  if (shape->type != regELLIPSE) {
    fprintf( stderr, "ERROR: regToStringEllipse() incorrect regShape type");
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
  
  char *r1 = calloc(bsize, sizeof(char));
  char *r2 = calloc(bsize, sizeof(char));
  reg_print_radius(shape->radius[0], shape->flag_radius, r1);
  reg_print_radius(shape->radius[1], shape->flag_radius, r2);
  
  char *a = calloc(bsize, sizeof(char));
  reg_print_val(shape->angle[0], a);
  
  snprintf(buf, maxlen, "Ellipse(%s,%s,%s,%s,%s)", x, y, r1, r2, a);

  free(x);
  free(y);
  free(r1);
  free(r2);
  free(a);
}
