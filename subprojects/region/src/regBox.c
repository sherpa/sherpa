#include "region_priv.h"


regShape* regCopyBox( regShape * );
int regIsEqualBox( regShape *, regShape * );
double regCalcAreaBox( regShape * );
int regCalcExtentBox( regShape *, double *, double * );
int regInsideBox( regShape *, double, double );
void regToStringBox( regShape *, char*, long );



regShape* regCreateBox(regFlavor include, 
		       double *xpos, 
		       double *ypos,
		       double *radius, 
		       double *angle, 
		       int wcoord,
               int wsize) {

  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreateBox() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  if (radius == NULL) {
     fprintf( stderr, "ERROR: regCreateBox() requires two (2) radii.");
     return NULL;
  }

  if ( radius[0] < 0 ) {
    fprintf( stderr, "ERROR: 1st radius of box must be positive\n");
    return(NULL);
  }

  if ( radius[1] < 0 ) {
    fprintf( stderr, "ERROR: 2nd radius of box must be positive\n");
    return(NULL);
  }

  regShape* newShape;
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Box";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( 1, sizeof(double ));
  newShape->ypos = ( double *)calloc( 1, sizeof(double) );

  newShape->xpos[0] = xpos[0]; 
  newShape->ypos[0] = ypos[0];

  newShape->include = include;
  newShape->nPoints = 1;

  newShape->angle = (double *)calloc(1,sizeof(double ));
  newShape->sin_theta = (double *)calloc(1,sizeof(double ));
  newShape->cos_theta = (double *)calloc(1,sizeof(double ));

  newShape->angle[0] = angle ? angle[0] : 0.0;

  double theta = newShape->angle[0]*PI/180.0;
  *newShape->sin_theta = sin(theta);
  *newShape->cos_theta = cos(theta);

  newShape->type = angle ? regROTBOX : regBOX;

  newShape->radius = (double *)calloc(2,sizeof(double ));
  newShape->radius[0] = radius[0];      
  newShape->radius[1] = radius[1];
 
  newShape->calcArea   = regCalcAreaBox;
  newShape->calcExtent = regCalcExtentBox;
  newShape->copy       = regCopyBox;
  newShape->isEqual    = regIsEqualBox;
  newShape->isInside   = regInsideBox;
  newShape->toString   = regToStringBox;
  newShape->free       = NULL;

  newShape->region = NULL;
  newShape->next   = NULL;

  return newShape;

}


// copy
regShape* regCopyBox( regShape* inShape ) {

  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyBox() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regBOX) {
    fprintf( stderr, "ERROR: regCopyBox() incorrect regShape type");
    return NULL;
  }

  return regCreateBox( inShape->include, inShape->xpos, 
		       inShape->ypos, inShape->radius,
		       inShape->angle, inShape->flag_coord, inShape->flag_radius);
}


// equals
int regIsEqualBox( regShape* thisShape, regShape* otherShape ) {

  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }
  
  if ( (thisShape->type != regBOX) ) { 
    fprintf( stderr, "ERROR: regIsEqualBox() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regBOX ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ||
       thisShape->radius == NULL ||  otherShape->radius == NULL ||
       thisShape->angle == NULL ||  otherShape->angle == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualBox() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) ||
       (thisShape->angle[0] != otherShape->angle[0])   ||
       (thisShape->radius[0] != otherShape->radius[0]) ||
       (thisShape->radius[1] != otherShape->radius[1]))
    return 0;
  
  return 1;
}


// calcArea
double regCalcAreaBox( regShape* shape ) {
    if (!shape) {
      fprintf( stderr, "ERROR: regCalcAreaBox() requires a regShape as input");
      return 0;
    }

    if (shape->type != regBOX) {
      fprintf( stderr, "ERROR: regCalcAreaBox() incorrect regShape type");
      return 0;
    }

    double area;
    area = shape->radius[1] * shape->radius[0];
    return area;
}


// calcExtent
int regCalcExtentBox( regShape* shape, double* xpos, double* ypos ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regCalcExtentBox() requires a regShape as input");
    return 0;
  }

  if (shape->type != regBOX) {
    fprintf( stderr, "ERROR: regCalcExtentBox() incorrect regShape type");
    return 0;
  }
  
  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCalcExtentBox() requires pre-allocated memory for xpos and ypos");
    return 0;
  }

  double xcor[4], ycor[4];
  reg_box_corners(shape, xcor, ycor);
  reg_corner_bounds(xcor, ycor, xpos, ypos);
  return 0;

}


// inside
int regInsideBox( regShape* shape, double xpos, double ypos ) {

    if (!shape) {
      fprintf( stderr, "ERROR: regInsideBox() requires a regShape as input");
      return 0;
    }
    
    if (shape->type != regBOX) {
      fprintf( stderr, "ERROR: regInsideBox() incorrect regShape type");
      return 0;
    }

    int retval;
    double dx;
    double dy;
    
    dx = shape->radius[0] / 2.0;
    dy = shape->radius[1] / 2.0;
    
    if ( shape->angle[0] != 0.0 )
    {
	    /*double theta = shape->angle[0]*PI/180.0;*/
	    double xm = xpos - shape->xpos[0];
	    double ym = ypos - shape->ypos[0];
	
	    double xp =  *shape->cos_theta * xm + *shape->sin_theta * ym;
	    double yp = -*shape->sin_theta * xm + *shape->cos_theta * ym;

	    if ( ( fabs( xp ) <= dx )  &&
	         ( fabs( yp ) <= dy )  )
	    {
	        retval = 1;
	    }
	    else
	        retval = 0;
    }
    else
    {
	    if  ( ( xpos >= shape->xpos[0] - dx ) &&
	          ( xpos <= shape->xpos[0] + dx ) &&
	          ( ypos >= shape->ypos[0] - dy ) &&
	          ( ypos <= shape->ypos[0] + dy ) )
	    {
	        retval = 1;
	    }
	    else
	        retval = 0;
    }
   
    if ( shape->include == regInclude )
      return(retval );
    else
      return( 1 - retval );
}


// toString
void regToStringBox( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringBox() requires a regShape as input");
    return;
  }

  if (shape->type != regBOX) {
    fprintf( stderr, "ERROR: regToStringBox() incorrect regShape type");
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
  
  if (shape->angle[0] == 0.0)
    snprintf(buf, maxlen, "Box(%s,%s,%s,%s)", x, y, r1, r2);
  else
    snprintf(buf, maxlen, "RotBox(%s,%s,%s,%s,%s)", x, y, r1, r2, a);
  
  free(x); 
  free(y);
  free(r1);
  free(r2);
  free(a);
}
