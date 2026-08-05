#include "region_priv.h"

regShape* regCopyPie( regShape * );
int regIsEqualPie( regShape *, regShape * );
double regCalcAreaPie( regShape * );
int regCalcExtentPie( regShape *, double *, double * );
int regInsidePie( regShape *, double, double );
void regToStringPie( regShape *, char*, long );

void reg_pie_bounds(double ang1, double ang2, double r2, double r1,
		    double *xoff, double *yoff);
long reg_quadrant(double ang1);


regShape* regCreatePie(regFlavor include, 
		       double *xpos, 
		       double *ypos, 
		       double *radius,
		       double *angle, 
		       int wcoord,
               int wsize) {


  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCreatePie() requires [xpos, ypos] coordinate pair.");
    return NULL;
  }

  if (radius == NULL) {
     fprintf( stderr, "ERROR: regCreatePie() requires two (2) radii.");
     return NULL;
  }

  if ( radius[0] < 0 ) {
    fprintf( stderr, "ERROR: regCreatePie() inner radius of pie must be positive\n");
    return NULL;
  }

  if ( radius[1] < 0 ) {
    fprintf( stderr, "ERROR: regCreatePie() outer radius of pie must be positive\n");
    return NULL;
  }

  if ( radius[1] < radius[0] ) {
    fprintf( stderr, "ERROR: regCreatePie() pie outer radius must be larger than inner radius\n");
    return NULL;
  }


  regShape* newShape;
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  newShape->name = "Pie";
  newShape->flag_coord = wcoord;
  newShape->flag_radius = wsize;

  newShape->xpos = ( double *)calloc( 1, sizeof(double) );
  newShape->ypos = ( double *)calloc( 1, sizeof(double) );

  newShape->xpos[0] = xpos[0]; 
  newShape->ypos[0] = ypos[0];

  newShape->type = regPIE;
  newShape->include = include;
  newShape->nPoints = 1;

  newShape->angle = (double *)calloc(2,sizeof(double ));
  newShape->sin_theta = (double *)calloc(1,sizeof(double ));
  newShape->cos_theta = (double *)calloc(1,sizeof(double ));

  newShape->angle[0] = angle ? angle[0] : 0.0;
  newShape->angle[1] = angle ? angle[1] : 0.0;
  double theta = newShape->angle[0]*PI/180.0;
  *newShape->sin_theta = sin(theta);
  *newShape->cos_theta = cos(theta);

  newShape->radius = (double *)calloc(2,sizeof(double ));  
  newShape->radius[0] = radius[0];      
  newShape->radius[1] = radius[1];
  
  newShape->calcArea   = regCalcAreaPie;
  newShape->calcExtent = regCalcExtentPie;
  newShape->copy       = regCopyPie;
  newShape->isEqual    = regIsEqualPie;
  newShape->isInside   = regInsidePie;
  newShape->toString   = regToStringPie;
  newShape->free       = NULL;

  newShape->region = NULL;
  newShape->next   = NULL;

  return newShape;
}


// copy
regShape* regCopyPie( regShape* inShape ) {

  if (inShape == NULL) {
    fprintf( stderr, "ERROR: regCopyPie() requires a regShape as input");
    return NULL;
  }
  
  if (inShape->type != regPIE) {
    fprintf( stderr, "ERROR: regCopyPie() incorrect regShape type");
    return NULL;
  }

  return regCreatePie( inShape->include, inShape->xpos, 
			   inShape->ypos, inShape->radius,
			   inShape->angle, inShape->flag_coord, inShape->flag_radius );
  
}


// equals
int regIsEqualPie( regShape* thisShape, regShape* otherShape ) {

  if (!thisShape && !otherShape) {
      return 1;
  }

  if (!thisShape || !otherShape) {
      return 0;
  }
  
  if ( (thisShape->type != regPIE) ) { 
    fprintf( stderr, "ERROR: regIsEqualPie() unable to compare shapes of different types");
    return 0;
  }

  if ( otherShape->type != regPIE ) {
    return 0;
  }

  if ( thisShape->xpos == NULL ||  otherShape->xpos == NULL ||
       thisShape->ypos == NULL ||  otherShape->ypos == NULL ||
       thisShape->radius == NULL ||  otherShape->radius == NULL ||
       thisShape->angle == NULL ||  otherShape->angle == NULL ) {
    fprintf( stderr, "ERROR: regIsEqualPie() unable to compare incomplete shapes");
    return 0;
  }

  if ( (thisShape->xpos[0] != otherShape->xpos[0]) ||
       (thisShape->ypos[0] != otherShape->ypos[0]) ||
       (thisShape->include != otherShape->include) ||
       (thisShape->radius[0] != otherShape->radius[0]) ||
       (thisShape->radius[1] != otherShape->radius[1]) ||
       (thisShape->angle[0] != otherShape->angle[0]) ||
       (thisShape->angle[1] != otherShape->angle[1]) ||
       (thisShape->flag_coord != otherShape->flag_coord) ||
       (thisShape->flag_radius != otherShape->flag_radius) )
    return 0;
  
  return 1;
}


// calcArea
double regCalcAreaPie( regShape* shape ) {

    if (shape == NULL) {
      fprintf( stderr, "ERROR: regCalcAreaPie() requires a regShape as input");
      return 0;
    }

    if (shape->type != regPIE) {
      fprintf( stderr, "ERROR: regCalcAreaPie() incorrect regShape type");
      return 0;
    }

    double area;
    double ang1, ang2, theta;
    ang1 = reg_mod_angle(shape->angle[0]);
    ang2 = reg_mod_angle(shape->angle[1]);

    if (ang1 < ang2)
	theta = ang2 - ang1;
    else
	theta = 360.0 - (ang1 - ang2);

    area = PI * (theta / 360.0) * (shape->radius[1] * shape->radius[1] -
				   shape->radius[0] * shape->radius[0]);
    return area;
}


// calcExtent
int regCalcExtentPie( regShape* shape, double* xpos, double* ypos ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regCalcExtentPie() requires a regShape as input");
    return 0;
  }

  if (shape->type != regPIE) {
    fprintf( stderr, "ERROR: regCalcExtentPie() incorrect regShape type");
    return 0;
  }
  
  if ((xpos == NULL) || (ypos == NULL)) {
    fprintf( stderr, "ERROR: regCalcExtentPie() requires pre-allocated memory for xpos and ypos");
    return 0;
  }

  double txpos[2];
  double typos[2];

  double ang1 = shape->angle[0];
  double ang2 = shape->angle[1];

  reg_pie_bounds(ang1, ang2, shape->radius[0], shape->radius[1],
		 txpos, typos);
  
  xpos[0] = txpos[0] + shape->xpos[0];
  xpos[1] = txpos[1] + shape->xpos[0];
  ypos[0] = typos[0] + shape->ypos[0];
  ypos[1] = typos[1] + shape->ypos[0];
  
  return 0;
}


// inside
int regInsidePie( regShape* shape, double xpos, double ypos ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regInsidePie() requires a regShape as input");
    return 0;
  }
  
  if (shape->type != regPIE) {
    fprintf( stderr, "ERROR: regInsidePie() incorrect regShape type");
    return 0;
  }
  
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
void regToStringPie( regShape* shape, char* buf, long maxlen ) {

  if (shape == NULL) {
    fprintf( stderr, "ERROR: regToStringPie() requires a regShape as input");
    return;
  }

  if (shape->type != regPIE) {
    fprintf( stderr, "ERROR: regToStringPie() incorrect regShape type");
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
  
  char *a1 = calloc(bsize, sizeof(char));
  char *a2 = calloc(bsize, sizeof(char));
  reg_print_val(shape->angle[0], a1);
  reg_print_val(shape->angle[1], a2);
  
  snprintf(buf, maxlen, "Pie(%s,%s,%s,%s,%s,%s)", x, y, r1, r2, a1, a2);
  
  free(x);
  free(y);
  free(r1);
  free(r2);
  free(a1);
  free(a2);
}


/***  Helper Functions ***/


void reg_pie_bounds(double ang1, double ang2, double r2, double r1,
		    double *xoff, double *yoff)
{
    long q1;
    long q2;
    double c1, c2;
    double s1, s2;
    double t1, t2;

    ang1 = reg_mod_angle(ang1);
    ang2 = reg_mod_angle(ang2);

    q1 = reg_quadrant(ang1);
    q2 = reg_quadrant(ang2);

    xoff[0] = -r1;
    xoff[1] = r1;
    yoff[0] = -r1;
    yoff[1] = r1;

/* Case of default bounds */
    if (q1 == q2 && ang1 >= ang2)
	return;

    c1 = cos(PI * ang1 / 180.0);
    c2 = cos(PI * ang2 / 180.0);
    s1 = sin(PI * ang1 / 180.0);
    s2 = sin(PI * ang2 / 180.0);
    if (q1 == 1) {
	if (q2 == 1) {
	    xoff[0] = r2 * c2;
	    xoff[1] = r1 * c1;
	    yoff[0] = r2 * s1;
	    yoff[1] = r1 * s2;
	} else if (q2 == 2) {
	    xoff[0] = r1 * c2;
	    xoff[1] = r1 * c1;
	    t1 = r2 * s1;
	    t2 = r2 * s2;
	    yoff[0] = t1 < t2 ? t1 : t2;
	} else if (q2 == 3) {
	    xoff[1] = r1 * c1;
	    yoff[0] = r1 * s2;
	} else if (q2 == 4) {
	    t1 = r1 * c1;
	    t2 = r1 * c2;
	    xoff[1] = t1 > t2 ? t1 : t2;
	}

    } else if (q1 == 2) {
	if (q2 == 1) {
	    t1 = r1 * s1;
	    t2 = r1 * s2;
	    yoff[1] = t1 > t2 ? t1 : t2;
	} else if (q2 == 2) {
	    xoff[0] = r1 * c2;
	    xoff[1] = r2 * c1;
	    yoff[0] = r2 * s2;
	    yoff[1] = r1 * s1;
	} else if (q2 == 3) {
	    t1 = r2 * c1;
	    t2 = r2 * c2;
	    xoff[1] = t1 > t2 ? t1 : t2;
	    yoff[0] = r1 * s2;
	    yoff[1] = r1 * s1;
	} else if (q2 == 4) {
	    xoff[1] = r1 * c2;
	    yoff[1] = r1 * s1;
	}

    } else if (q1 == 3) {
	if (q2 == 1) {
	    xoff[0] = r1 * c1;
	    yoff[1] = r1 * s2;
	} else if (q2 == 2) {
	    t1 = r1 * c1;
	    t2 = r1 * c2;
	    xoff[0] = t1 < t2 ? t1 : t2;
	} else if (q2 == 3) {
	    xoff[0] = r1 * c1;
	    xoff[1] = r2 * c2;
	    yoff[0] = r1 * s2;
	    yoff[1] = r2 * s1;
	} else if (q2 == 4) {
	    xoff[0] = r1 * c1;
	    xoff[1] = r1 * c2;
	    t1 = r2 * s1;
	    t2 = r2 * s2;
	    yoff[1] = (t1 > t2) ? t1 : t2;
	}
    } else if (q1 == 4) {
	if (q2 == 1) {
	    t1 = r2 * c1;
	    t2 = r2 * c2;
	    xoff[0] = (t1 < t2) ? t1 : t2;
	    yoff[0] = r1 * s1;
	    yoff[1] = r1 * s2;
	} else if (q2 == 2) {
	    xoff[0] = r1 * c2;
	    yoff[0] = r1 * s1;
	} else if (q2 == 3) {
	    t1 = r1 * s1;
	    t2 = r1 * s2;
	    yoff[0] = (t1 < t2) ? t1 : t2;
	} else if (q2 == 4) {
	    xoff[0] = r2 * c1;
	    xoff[1] = r1 * c2;
	    yoff[0] = r1 * s1;
	    yoff[1] = r2 * s2;
	}
    }
}



long reg_quadrant(double ang1)
{
    double ang = ang1;
    long quad = 4;

    
    while (ang < 0.0)
      ang += 360.0;
    if (ang == 360.0)
      return quad;
    ang = fmod( ang, 360.0 );
    quad = (ang / 90.0) + 1;
    return quad;
}
