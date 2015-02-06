/*                                                                
**  Copyright (C) 2007,2013  Smithsonian Astrophysical Observatory 
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

#include "region_priv.h"
#include <float.h>
#define MEMCHECK(xx) if ((xx)==NULL) return(NULL)
int regIntersectComponent( regRegion* region, regShape* Shape1, regShape* Shape2 );
int regShapeIntersect( regShape* shape1, regShape* shape2, long* index1, long* index2 );
int regExtentShapeRaw( regShape* shape, double* fieldx, double* fieldy, double* xpos, double* ypos );
int regRectangleInside( double* xpos1, double* ypos1, double* xpos2, double* ypos2 );
int regIsRect( regShape* shape );



/* regTerm structure
   Used to hold information about regions we are inverting 
*/
struct regTerm 
{
  struct regTerm *next;  /* next OR term */
  struct regTerm *prev;	/* previous OR term */
  regShape *first;	/* first shape in term */
  regShape *current;	/* current shape in term */
  regShape *last;	/* last shape in term */
};


/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
regRegion* regCreateEmptyRegion( void )
{
 return regCreateRegion( NULL,NULL );
}

regRegion* regCreateRegion( void *xcol, void*ycol )
{
  regRegion *region = ( regRegion *) calloc( 1, sizeof( regRegion ) );
  
  region->xcol[0] = xcol; 
  region->xcol[1] = ycol;
  return( region );

}


/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

long regAddShape(
		 regRegion *region,
		 regMath glue,
		 regShape *shape
		 )
{
  double fx[2] ={ -DBL_MAX, DBL_MAX };
  double fy[2] ={ -DBL_MAX, DBL_MAX };
  regShape *atShape;

  if ( region->shape == NULL)
    {
      shape->component =  1;
      region->shape = shape;
    }
  else
    {
      atShape = region->shape;
      while ( atShape->next != NULL)
	{
	  atShape = atShape->next;
	}
      atShape->next = shape;
      
      if ( glue == regAND )
	{
	  shape->component = atShape->component;
	}
      else
	{
	  shape->component = atShape->component+1;
	}

    }

  shape->region = region;
  regExtent(region, fx, fy, region->xregbounds, region->yregbounds);

  return( shape->component );

}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
regShape *regCreateWorldShape(
			 regRegion *region,
			 regMath glue,
			 regGeometry type,
			 regFlavor include,
			 double *xpos,
			 double *ypos,
			 long   npoints,
			 double *radius,
			 double *angle,
                         int wcoord,
                         int wsize
			 )
{
  regShape *newShape = NULL;
  long ii, jj;
  double theta = 0;

  /* Region library can not create user-defined shapes */
  if ( type == regMASK || type == regUSER){ return newShape; };

  /* copy all inputs */
  newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );
  MEMCHECK(newShape);

  newShape->world_coord = wcoord;
  newShape->world_size  = wsize;

  newShape->xpos = ( double *)calloc( npoints, sizeof(double ));
  newShape->ypos = ( double *)calloc( npoints, sizeof(double) );

  newShape->type = type;
  newShape->include = include;
  newShape->nPoints = npoints;

  /* just do a memcopy? */
  for ( ii=0; ii<npoints;ii++)
    {
      newShape->xpos[ii] = xpos[ii]; 
      newShape->ypos[ii] = ypos[ii];
    }

  if ( ( type == regPOLYGON ) &&  (
       ( newShape->xpos[0] != newShape->xpos[npoints-1] ) ||
       ( newShape->ypos[0] != newShape->ypos[npoints-1] )   ) )
    {
      npoints += 1;
      newShape->nPoints = npoints;      
      newShape->xpos = realloc( newShape->xpos, npoints*sizeof(double));
      newShape->ypos = realloc( newShape->ypos, npoints*sizeof(double));
      newShape->xpos[npoints-1]=xpos[0];
      newShape->ypos[npoints-1]=ypos[0];
    }



  switch ( type )
    {
    case regPOINT:
      newShape->angle = NULL;
      newShape->radius = NULL;
      newShape->inside = regInsidePoint;
      break;
    
    case regCIRCLE:
      newShape->angle = NULL;
      MEMCHECK( radius );
      newShape->radius = ( double *)calloc(1,sizeof(double));
      MEMCHECK(newShape->radius);
      if ( radius[0] < 0 ) {
	fprintf( stderr, "ERROR: radius of circle must be positive\n");
	return(NULL);
      }
      newShape->radius[0] = radius[0];
      newShape->inside = regInsideCircle;
      break;
 
    case regANNULUS:
      newShape->angle = NULL;
      MEMCHECK( radius );
      newShape->radius = ( double *)calloc(2,sizeof(double));
      MEMCHECK(newShape->radius);
      if ( radius[0] < 0 ) {
	fprintf( stderr, "ERROR: inner radius of annulus must be positive\n");
	return(NULL);
      }
      if ( radius[1] < 0 ) {
	fprintf( stderr, "ERROR: outer radius of annulus must be positive\n");
	return(NULL);
      }
      if ( radius[1] < radius[0] ) {
	fprintf( stderr, "ERROR: annulus outer radius must be larger than "
		 "inner radius\n");
	return(NULL);
      }
      newShape->radius[0] = radius[0];
      newShape->radius[1] = radius[1];
      newShape->inside = regInsideAnnulus;
      break;
 
    case regELLIPSE:
      MEMCHECK( radius );
      MEMCHECK( angle );
      newShape->angle = (double *)calloc(1,sizeof(double ));
      newShape->sin_theta = (double *)calloc(1,sizeof(double ));
      newShape->cos_theta = (double *)calloc(1,sizeof(double ));
      MEMCHECK(newShape->angle);
      MEMCHECK(newShape->sin_theta);
      MEMCHECK(newShape->cos_theta);
      newShape->angle[0] = angle[0];
      theta = newShape->angle[0]*PI/180.0;
      *newShape->sin_theta = sin(theta);
      *newShape->cos_theta = cos(theta);
      newShape->radius = (double *)calloc(2,sizeof(double ));
      MEMCHECK(newShape->radius);
      if ( radius[0] < 0 ) {
	fprintf( stderr, "ERROR: 1st radius of ellipse must be positive\n");
	return(NULL);
      }
      if ( radius[1] < 0 ) {
	fprintf( stderr, "ERROR: 2nd radius of ellipse must be positive\n");
	return(NULL);
      }
      newShape->radius[0] = radius[0];      
      newShape->radius[1] = radius[1];
      newShape->inside = regInsideEllipse;
      break;
         
    case regBOX:
      MEMCHECK( radius );
      newShape->angle = (double *)calloc(1,sizeof(double ));
      newShape->sin_theta = (double *)calloc(1,sizeof(double ));
      newShape->cos_theta = (double *)calloc(1,sizeof(double ));
      MEMCHECK(newShape->angle);
      MEMCHECK(newShape->sin_theta);
      MEMCHECK(newShape->cos_theta);
      newShape->angle[0] = angle ? angle[0] : 0.0;
      theta = newShape->angle[0]*PI/180.0;
      *newShape->sin_theta = sin(theta);
      *newShape->cos_theta = cos(theta);
      newShape->radius = (double *)calloc(2,sizeof(double ));
      MEMCHECK(newShape->radius);
      if ( radius[0] < 0 ) {
	fprintf( stderr, "ERROR: 1st radius of box must be positive\n");
	return(NULL);
      }
      if ( radius[1] < 0 ) {
	fprintf( stderr, "ERROR: 2nd radius of box must be positive\n");
	return(NULL);
      }
      newShape->radius[0] = radius[0];      
      newShape->radius[1] = radius[1];
      newShape->inside = regInsideBox;
      break;

    case regRECTANGLE:
      newShape->angle = (double *)calloc(1,sizeof(double ));
      newShape->sin_theta = (double *)calloc(1,sizeof(double ));
      newShape->cos_theta = (double *)calloc(1,sizeof(double ));
      MEMCHECK(newShape->angle);
      MEMCHECK(newShape->sin_theta);
      MEMCHECK(newShape->cos_theta);
      newShape->angle[0] = angle ? angle[0] : 0.0;
      theta = newShape->angle[0]*PI/180.0;
      *newShape->sin_theta = sin(theta);
      *newShape->cos_theta = cos(theta);
      newShape->radius = NULL;
      newShape->inside = regInsideRectangle;
      break;


    case regPOLYGON:
      newShape->angle = NULL;
      newShape->radius = NULL;
      newShape->inside = regInsidePolygon;
      
      if ( newShape->nPoints < 4 ) {
	fprintf(stderr, "ERROR: Polygons must have at least 3 vertexes\n");
	return(NULL);
      }
      
      for (ii=0;ii<npoints-2;ii++) {
	if ((newShape->xpos[ii] == newShape->xpos[ii+2]) &&
	    (newShape->ypos[ii] == newShape->ypos[ii+2]) &&
	    (ii+2 != (npoints-1)) ) {
	  fprintf(stderr, "WARNING: Polgyon must have finite width; adjacent line segments with ends at (%g,%g) "
		  "overlap completely\n",
		  newShape->xpos[ii],newShape->ypos[ii] );
	  //return(NULL);
	}
      }
      
      for (ii=0;ii<npoints-2;ii++) {
	if ((newShape->xpos[ii] == newShape->xpos[ii+1]) &&
	    (newShape->ypos[ii] == newShape->ypos[ii+1]) ) {
	  fprintf(stderr, "WARNING: Zero length polygon line segment at (%g,%g).\n",
		  newShape->xpos[ii], newShape->ypos[ii]);
	  
	  break;
	}
      }
      

      /* This does a check to see if any of the lines of the
	 polgyons over lap (thus making a complex polygon).

	 Make a warning for now
      */
      
      for (ii=0;ii<npoints-3;ii++) {
	double x1, x2, y1, y2;

	x1=newShape->xpos[ii];
	x2=newShape->xpos[ii+1];
	y1=newShape->ypos[ii];
	y2=newShape->ypos[ii+1];

	for (jj=ii+2;jj<=npoints-2;jj++) {
	  double x3,x4, y3,y4;
	  
	  /* last point is == 1st by definition */
	  if ( ( ii == 0 ) && ( jj == npoints-2 )) continue;

	  
	  x3=newShape->xpos[jj];
	  x4=newShape->xpos[jj+1];
	  y3=newShape->ypos[jj];
	  y4=newShape->ypos[jj+1];


	  double ua;
	  double ub;
	  double denom;

	  denom=(y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
	  if ( fabs(denom) < FLT_MIN) {
	    // nothing; good, lines are pallel
	  } else {
	    ua=(x4-x3)*(y1-y3)-(y4-y3)*(x1-x3);
	    ub=(x2-x1)*(y1-y3)-(y2-y1)*(x1-x3);
	    
	    ua/= denom;
	    ub/= denom;
	    
	    if ( (0 <= ua) && ( ua <= 1 ) &&
		 (0 <= ub) && ( ub <= 1 )    ) {
	      fprintf(stderr, "WARNING: The lines of the polygon cannot "
		      "intersect.  Area may not be correct\n");

	      /* ---- good for debugging
		 fprintf( stderr, "  linesegment #1:  (xlo,ylo) = (%g,%g)\t (xhi,yhi) = (%g,%g)\n",x1,y1,x2,y2);
		 fprintf( stderr, "  linesegment #2:  (xlo,ylo) = (%g,%g)\t (xhi,yhi) = (%g,%g)\n",x3,y3,x4,y4);
		 fprintf( stderr, "  ii=%ld \t jj=%ld\n", ii, jj);
	      */

	      /* Go ahead and break out; no need to test the remaining data points. */
	      ii=npoints;
	      jj=npoints;
	    }

	  } /* end else */
	} /* end jj */
      } /* end ii */
      
      
      break;

    case regPIE:
      MEMCHECK( angle );
      newShape->angle = (double *)calloc(2,sizeof(double ));
      newShape->sin_theta = (double *)calloc(1,sizeof(double ));
      newShape->cos_theta = (double *)calloc(1,sizeof(double ));
      MEMCHECK(newShape->angle);
      MEMCHECK(newShape->sin_theta);
      MEMCHECK(newShape->cos_theta);
      newShape->angle[0] = angle[0];
      newShape->angle[1] = angle[1];
      theta = newShape->angle[0]*PI/180.0;
      *newShape->sin_theta = sin(theta);
      *newShape->cos_theta = cos(theta);
      newShape->radius = (double *)calloc(2,sizeof(double ));
      MEMCHECK(newShape->radius);
      if ( radius[0] < 0 ) {
	fprintf( stderr, "ERROR: inner radius of pie must be positive\n");
	return(NULL);
      }
      if ( radius[1] < 0 ) {
	fprintf( stderr, "ERROR: outer radius of pie must be positive\n");
	return(NULL);
      }
      if ( radius[1] < radius[0] ) {
	fprintf( stderr, "ERROR: pie outer radius must be larger than "
		 "inner radius\n");
	return(NULL);
      }
      newShape->radius[0] = radius[0];      
      newShape->radius[1] = radius[1];
      newShape->inside = regInsidePie;
      break;

    case regSECTOR:
      MEMCHECK( angle );
      newShape->angle = (double *)calloc(2,sizeof(double ));
      newShape->sin_theta = (double *)calloc(1,sizeof(double ));
      newShape->cos_theta = (double *)calloc(1,sizeof(double ));
      MEMCHECK(newShape->angle);
      MEMCHECK(newShape->sin_theta);
      MEMCHECK(newShape->cos_theta);
      newShape->angle[0] = angle[0];
      newShape->angle[1] = angle[1];
      theta = newShape->angle[0]*PI/180.0;
      *newShape->sin_theta = sin(theta);
      *newShape->cos_theta = cos(theta);
      newShape->radius = NULL;
      newShape->inside = regInsidePie;
      break;

    case regFIELD:
      newShape->angle = NULL;
      newShape->radius = NULL;
      newShape->inside = regInsideField;
      break;

   default:
      return(NULL);
    }
  /* attributes get assinged/init when glued to a region */

  if ( region != NULL )
    regAddShape( region, glue, newShape );
 
  return( newShape );
}


regRegion* regCreateRectangle( double* x, double* y )
{
 regRegion*  field = regCreateRegion( NULL, NULL );
 regCreateShape( field, regAND, regRECTANGLE, regInclude, x, y, 2, NULL, NULL );
 return field;
}

void regAppendShape( regRegion* region, 
		     char* shapeName, int includeFlag, int orFlag, double* xpos, double* ypos,
		     long nmax, double* radius, double* angle, int world_coord, int world_size )
{
 regMath glue;
 regGeometry type;
 regFlavor include;
 long npoints =0;
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };

 glue    = orFlag ? regOR : regAND;
 include = includeFlag? regInclude : regExclude;
 if ( !strcmp( shapeName, "npolygon" ) || !strcmp( shapeName, "NPOLYGON" ))
 {
  type = regPOLYGON;
  npoints = nmax;
 } else {
   /* not able for regUSER case */
   type = regShapeNameToGeometry( shapeName ); 
   npoints = regShapeFindNPoints( type, xpos, ypos, nmax ); 
 }
       
 if ( type == regMASK || type== regUSER){ return; /* add some warning */}

 regCreateWorldShape( region, glue, type, include, xpos, ypos, npoints,
			radius, angle, world_coord, world_size );

 regExtent(region, fx, fy, region->xregbounds, region->yregbounds);
}


void regResolveField( regRegion* region, double* x, double* y )
{
  regShape* Shape;
  if ( !region ) return;

  Shape = region->shape;

  if ( Shape->type == regFIELD )
  {
   Shape->type = regRECTANGLE;
   Shape->xpos = calloc( 2, sizeof( double ));
   Shape->ypos = calloc( 2, sizeof( double ));
   if ( !Shape->angle ) Shape->angle = calloc( 1, sizeof( double ));
   Shape->nPoints = 2;
   Shape->xpos[0] = x[0];
   Shape->xpos[1] = x[1];
   Shape->ypos[0] = y[0];
   Shape->ypos[1] = y[1];
  }
}

/* --------------------------------------------------------------------------- */
/* regCreateNewWorldShape                                                      */
/*   Create and populate Shape object with provided inputs.                    */
/*                                                                             */
/*   NOTE: not appropriate for user defined shapes ( regUSER,regMASK)          */
/* --------------------------------------------------------------------------- */
regShape *regCreateNewWorldShape(
			 regGeometry type,
			 regFlavor include,
			 double *xpos,
			 double *ypos,
			 long   npoints,
			 double *radius,
			 double *angle,
                         int world_coord,    /* xpos ypos in degrees or pixels? */
                         int world_size      /* radius in degrees or pixels? */
			 )
{
  return( regCreateWorldShape( NULL, regOR, type, include, xpos,
			  ypos, npoints, radius, angle, world_coord, world_size ) );
}


regShape* regNextComponent( regShape* shape )
{
 long cpt;
 if ( shape ) cpt = shape->component;
 while ( shape && shape->component == cpt )
  shape = shape->next;
 return shape;
}

/* --------------------------------------------------------------------------- */
/* regCopyShape                                                                */
/*   Create a new shape with same definition as the old one.                   */
/*   NOTE:                                                                     */
/*     component, region, next attributes are NOT copied                       */
/* --------------------------------------------------------------------------- */
regShape* regCopyShape( regShape* inShape )
{ 
  regShape *shp=NULL;
  if ( !inShape ) return NULL;
  
  if(inShape->user)
  {
    shp=inShape->user->copy(inShape);
  } else { 
    shp = regCreateNewWorldShape( inShape->type, inShape->include, 
				 inShape->xpos, inShape->ypos, inShape->nPoints, inShape->radius,
				 inShape->angle, inShape->world_coord, inShape->world_size );
  }

  return shp;
}

/*
 *  Reverse the polarity
 */
void regNegate(regShape* Shape )
{
 Shape->include = (Shape->include)? regExclude : regInclude; 
}

/*
 *  Convert a region which may contain world coord shapes to
 *  one with only pixel coord shapes. Requires a function to convert
 *  from world to pixel, e.g.
 *       void CoordInvert( void* data, double* world, double* pixel )
 * The first argument is provided for users to pass in control info.
 */
void regConvertWorldRegion(
    regRegion* Region,         /* (i) Region which may have shapes in world units */
    double scale,              /* (i) degrees per pixel */
    regInvertFunction invert      /* (i) Function to convert world to pixel coords */
 )
{
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };
 regShape* Shape;
 int force = 0;
 Shape = Region->shape;

 while (Shape != NULL )
   {
     regConvertWorldShape( Shape, scale, invert, force );
     Shape = Shape->next;
   }

 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
 
}

void regConvertRegion(
    regRegion* Region,         /* (i) Region which may have shapes in world units */
    double scale,              /* (i) degrees per pixel */
    regInvertFunction invert,      /* (i) Function to convert world to pixel coords */
    int force
 )
{
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };
 regShape* Shape;

 Shape = Region->shape;

 while (Shape != NULL )
    {
      regConvertWorldShape( Shape, scale, invert, force );
      Shape = Shape->next;
    }

 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);

}


/*
 *  Convert a shape expressed in degrees to one in pixels.
 *  Positions are converted using the invert function
 *  Radii are converted by dividing by the CDELT scale value
 *  We support the possibility that positions are expressed in pixels
 *  but radii in angular measure, and vice versa.
 */
void regConvertWorldShape( regShape* Shape, double scale, regInvertFunction invert, int force )
{
 long ii;
 long nradii = 0;
 long npoints = Shape->nPoints;
 double world[2];
 double pixel[2];

 if ( Shape->world_coord || force )
 {
  for ( ii=0; ii<npoints;ii++)
    {
      world[0] = Shape->xpos[ii];
      world[1] = Shape->ypos[ii];
      invert( world, pixel );
      Shape->xpos[ii] = pixel[0];       
      Shape->ypos[ii] = pixel[1];       
    }
 }
 Shape->world_coord = 0;

 if ( Shape->world_size || force )
 {
  nradii = regShapeRadii( Shape );

  for ( ii = 0; ii < nradii; ii++ )
  {
   Shape->radius[ii] = 	Shape->radius[ii] / scale;
  }
 }
 Shape->world_size = 0;
}


regRegion* regCopyRegion( regRegion* inRegion )
{
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };
 regRegion* Region;
 regShape* Shape;
 regShape* inShape;
 regMath glue;

 long lastComponent = 1;

/* Copy shapes */
 
 if ( !inRegion )
  return NULL;

 Region = regCreateRegion( inRegion->xcol[0], inRegion->xcol[1] );

 inShape = inRegion->shape;

 while (inShape != NULL )
    {
      Shape = regCopyShape( inShape );
      if ( inShape->component == lastComponent )
      {
       glue = regAND;
      } else {
       glue = regOR;
      }
      lastComponent = inShape->component;
      regAddShape( Region, glue, Shape );

      inShape = inShape->next;
    }
 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);

 return Region;
}

/* Invert a region

   We come here with a region structure containing a linked list of "shapes"

   Region->shapeA->shapeB->shapeC->shapeD etc

   Each shape has a flag which indicates its relation to the previous shape,
   whether AND or OR.  AND operations take precedence over OR, so we can
   view the shape list as a series of AND terms separated by OR operations:


	P*Q+R*S*T+W*Z = (P*Q)+(R*S*T)+(W*Z)

   This inversion function has two stages.  First it creates a linked
   list of terms.  A term is a set of shapes grouped by AND.  Successive
   entries in the term list are joined by OR.


   From DeMorgan's laws:

   ___________________     _ _ _   _ _ _   _ _ _   _ _ _   _ _ _   _ _ _
   (P*Q)+(R*S*T)+(W*Z)  =  P*R*W + P*R*Z + P*S*W + P*S*Z + P*T*W + P*T*Z +
                           _ _ _   _ _ _   _ _ _   _ _ _   _ _ _   _ _ _
   			   Q*R*W + Q*R*Z + Q*S*W + Q*S*Z + Q*T*W + Q*T*Z 

   Stated more simply, the left hand side is composed of all combinations
   of one shape chosen from each term, with each shape negated.  This is NOT
   the most compact way to write the results of DeMorgan's transformation,
   but is a representation which requires no explicit grouping syntax: AND 
   over OR precedence being adequate.

   The second stage of this function cycles through the term list, pulling
   shapes out, negating them, and adding them to the inverted region shape
   list.
*/
regRegion* regInvert( regRegion* inRegion )
{
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };
 regRegion* Region;
 regShape* Shape;
 regShape* inShape;
 struct regTerm *term;
 struct regTerm *nextTerm;
 struct regTerm *firstTerm;
 struct regTerm *safeTerm;
#ifndef TRUE
#define TRUE 1
#define FALSE 0
 typedef int boolean;
#endif

 boolean done;

/* Copy shapes */
 
 if ( !inRegion )
  return NULL;

 Region = regCreateRegion( inRegion->xcol[0], inRegion->xcol[1] );


 if (inRegion->shape == NULL)
   return Region;

 inShape = inRegion->shape;

 /* Start with field */
 Shape = regCreateNewWorldShape( regFIELD, regInclude, NULL, NULL, 0, NULL, NULL, 1.0, 0 );


 /* parse the region shapes in "terms" connected by OR */

 term = (struct regTerm *) (malloc(sizeof(struct regTerm)));
 term->next = NULL;
 term->prev = NULL;
			  
 firstTerm = term;

 /* initialize first term */
 term->first = inShape;
 term->current = inShape;
 term->last = inShape;

 /* add shapes to the current term */
 while (inShape->next != NULL)
   {
     /* look for OR glue */
     if (inShape->component != inShape->next->component)
       {
	 /* found it - finish off current term, and start new one */
	 term->last = inShape;  /* remember this shape as "last' */

	 /* malloc new term */
	 nextTerm = (struct regTerm *) (malloc(sizeof(struct regTerm)));

	 /* init new term */
	 nextTerm->first = inShape->next;
	 nextTerm->current = inShape->next;
	 term->next = nextTerm;
	 nextTerm->prev = term;
	 nextTerm->next = NULL;

	 /* new term is current term */
	 term = nextTerm;
       }
     inShape = inShape->next;
   }
 term->last = inShape;

 /* now do the inversion math */

 done = FALSE;
 while (!done)
   {
     /* get a shape from each term */
     term = firstTerm;
     do
       {
	 /* remember our current term (we will advance past it) */
	 safeTerm = term;
	 Shape = regCopyShape(term->current );
	 regNegate(Shape);  
	 
	 /* first shape in term is "joined" to predecessor (nothing) by OR */
	 if (term == firstTerm)
	   regAddShape( Region, regOR, Shape );
	 else
	   /* all other terms join shapes by AND */
	   regAddShape( Region, regAND, Shape );

	 term = term->next;
       }
     while (term != NULL);

     /* we have now got a shape from each term */
     term = safeTerm;

     /* advance to next shape in last term */
     if (term->current != term->last)
       term->current = term->current->next;
     else
       {
	 /* already at last shape, reset current term to first */
	 if (term != firstTerm)
	   term->current = term->first;
	 else
	   done = TRUE;

	 while (term != firstTerm)
	   {
	     /* go to previous term and advance or reset as appropriate */
	     term = term->prev;

	     /* any shapes left? */
	     if (term->current != term->last)
	       {
		 /* yes advance and break out of loop */
		 term->current = term->current->next;
		 break;
	       }
	     else
	       /* no - reset to first shape...*/
	       if (term != firstTerm)
		 term->current = term->first;
	       else
		 /* ...unless this is the first term,
		    in which case we are done */
		 done = TRUE;
	   }
       }
   }


 /* dispose of our parsing structures */
 term = firstTerm;
 do
   {
     nextTerm = term->next;
     free(term);
     term = nextTerm;
   }
 while (term != NULL);

 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
 return Region;
}


/* --------------------------------------------------------------------------- */
/* regCompareRegion                                                            */
/*    Return 1 if regions are identical                                        */
/* --------------------------------------------------------------------------- */
int regCompareRegion( regRegion* Region1, regRegion* Region2 )
{
 regShape* Shape1;
 regShape* Shape2;
 int true  = 1;
 int false = 0;
 Shape1 = Region1->shape;
 Shape2 = Region2->shape;

 while (Shape1 != NULL )
 {
  if ( !Shape2 )
   return false;

  if ( Shape1->component != Shape2->component )
   return false;

  if( !regCompareShape( Shape1, Shape2 ) )
   return false;

  Shape1 = Shape1->next;
  Shape2 = Shape2->next;
 }

 if ( Shape2 )
  return false;

 return true;
}



int regCompareShape( regShape* Shape1, regShape* Shape2 )
{
 int false = 0;

 if( Shape1->include != Shape2->include )
  return false;

 return regCompareShapeRaw(Shape1, Shape2);

}

/* Compare shapes ignoring positive or negative */
int regCompareShapeRaw( regShape* Shape1, regShape* Shape2 )
{
 int true  = 1;
 int false = 0;
 long nradii;
 long nangles;
 long i;

 
 if ( Shape1->type != Shape2->type )
  return false;

 
 if(Shape1->user)
   return Shape1->user->isEqual(Shape1,Shape2);

 nradii = regShapeRadii( Shape1 );
 for ( i = 0; i < nradii; i++ )
  if ( Shape1->radius[i] != Shape2->radius[i] )
    return false;

 nangles= regShapeAngles( Shape1 );
 for ( i = 0; i < nangles; i++ )
  if ( Shape1->angle[i] != Shape2->angle[i] )
    return false;

 if ( Shape1->nPoints != Shape2->nPoints )
  return false;

 for ( i = 0; i < Shape1->nPoints; i++ )
 {
  if ( Shape1->xpos[i] != Shape2->xpos[i] )
   return false;
  
  if ( Shape1->ypos[i] != Shape2->ypos[i] )
   return false;
  
 }

 return true;
}


regRegion* regCombineRegion( regRegion* Region1, regRegion* Region2 )
{
 double fx[2] ={ -DBL_MAX, DBL_MAX };
 double fy[2] ={ -DBL_MAX, DBL_MAX };
 regRegion* Region;
 regShape* Shape;
 regShape* inShape;
 regMath glue;

 long lastComponent = 1;

/* Copy shapes */
 if ( !Region1 )
 {
  if ( !Region2 )
   return NULL;
  return regCopyRegion( Region2 );
 }
 Region = regCopyRegion( Region1 );

 inShape = Region2->shape;

 while (inShape != NULL )
    {

      Shape = regCopyShape( inShape );
      if ( inShape->component == lastComponent )
      {
       glue = regAND;
      } else {
       glue = regOR;
      }
      lastComponent = inShape->component;
      regAddShape( Region, glue, Shape );

      inShape = inShape->next;
    }
 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
 return Region;
}

regRegion* regUnionRegion( regRegion* Region1, regRegion* Region2 )
{
  double fx[2] ={ -DBL_MAX, DBL_MAX };
  double fy[2] ={ -DBL_MAX, DBL_MAX };

  regRegion* region;
  regShape*  shape;
  regShape*  inShape;
  regMath    glue;
  
  long lastComponent;
  int  haveMore;
  
  /* Copy shapes */
  if ( !Region1 )
  {
   if ( !Region2 )
    return NULL;
   return regCopyRegion( Region2 );
  }
  
  /* If regions are equal just return copy of one */
  if ( regCompareRegion( Region1, Region2 ))
   return regCopyRegion( Region1 );
  
  /* Make a new region with all the combined components of the input regions */
  /*  - Put USER defined shapes first in the series.                         */

  region = regCreateRegion( Region1->xcol[0], Region1->xcol[1] );

  /* Transfer Region components with USER shapes.                            */
  /*   NOTE: expectation is that USER shapes are first on the component      */
  haveMore = 1;
  inShape = Region1->shape;
  while ( inShape != NULL )
  {
    if ( ( inShape->type == regMASK ) || ( inShape->type == regUSER ) )
    {
      glue = regOR;

      lastComponent = inShape->component;
      while ( inShape && inShape->component == lastComponent )
      {
        shape = regCopyShape( inShape );
        regAddShape( region, glue, shape );
        glue = regAND;
        inShape = inShape->next;
      }
    }
    else
    {
      inShape = regNextComponent( inShape );
    }

    if ( (inShape == NULL)&& haveMore )
    {
      /* scan second region components */
      inShape = Region2->shape;
      haveMore = 0;
    }
  }
  
  /* Transfer Region components Non-USER shapes.                           */
  haveMore = 1;
  inShape = Region1->shape;
  while ( inShape != NULL )
  {
    if ( ( inShape->type == regMASK ) || ( inShape->type == regUSER ) )
    {
      inShape = regNextComponent( inShape );
    }
    else
    {
      glue = regOR;

      lastComponent = inShape->component;
      while ( inShape && inShape->component == lastComponent )
      {
        shape = regCopyShape( inShape );
        regAddShape( region, glue, shape );
        glue = regAND;
        inShape = inShape->next;
      }
    }

    if ( (inShape == NULL)&& haveMore )
    {
      /* scan second region components */
      inShape = Region2->shape;
      haveMore = 0;
    }
  }

  /* re-calculate the region extent */
  regExtent(region, fx, fy, region->xregbounds, region->yregbounds);

  return region;
}


/*
 *  Do two components possibly have nonzero overlap?
 */
int regComponentOverlap( regShape* Shape1, regShape* Shape2 )
{
 long Cpt1, Cpt2;
 int ok = 1;
 regShape* shape;
 if ( Shape1 ) Cpt1 = Shape1->component;
 if ( Shape2 ) Cpt2 = Shape2->component;
 while (Shape1 != NULL && Shape1->component == Cpt1 && ok )
 {
  shape = Shape2;
  while (shape != NULL && shape->component == Cpt2 && ok )
    {
      ok = regShapeOverlap( Shape1, shape );
      shape = shape->next;
    }   
  Shape1 = Shape1->next;
 }      
 return ok;
}

int regIntersectComponent( regRegion* region, regShape* Shape1, regShape* Shape2 )
{
 long Cpt1, Cpt2;
 int ok = 1;
 regShape* shape1;
 regShape* shape2;
 regShape** nshape1;
 regShape** nshape2;
 regMath glue = regOR;
 long n1=0;
 long n2=0;
 long i;
 long j;
 int bIntersectOK = 1;
 long* index1;
 long* index2;
 long* isuser1;
 long* isuser2;

 if ( !Shape1 || !Shape2 )
   return 0;

 Cpt1 = Shape1->component;
 Cpt2 = Shape2->component;

/* Count shapes */
 shape1 = Shape1;
 shape2 = Shape2;

 while (shape1 != NULL && shape1->component == Cpt1 && ok )
 {
  ++n1;
  shape1 = shape1->next;
 }
 while (shape2 != NULL && shape2->component == Cpt2 && ok )
 {
  ++n2;
  shape2 = shape2->next;
 }

 index1  = calloc( n1, sizeof( long ));
 index2  = calloc( n2, sizeof( long ));
 isuser1 = calloc( n1, sizeof( long ));
 isuser2 = calloc( n2, sizeof( long ));
 nshape1 = calloc( n1, sizeof( regShape* ));
 nshape2 = calloc( n2, sizeof( regShape* ));
 shape1  = Shape1;
 shape2  = Shape2;

 /* Make copies of Set 1 shapes.. set flag for all shapes 'ON' */
 for ( i = 0; i < n1 && shape1; i++ )
 {
   nshape1[i] = regCopyShape( shape1 );  
   index1[i] = 1;
   isuser1[i] = 0;
   if (nshape1[i]->type==regMASK || nshape1[i]->type==regUSER)
     isuser1[i] = 1;
   shape1 = shape1->next;
 }

 /* Make copies of Set 2 shapes.. set flag for all shapes 'ON' */
 for ( i = 0; i < n2 && shape2; i++ )
 {
   nshape2[i] = regCopyShape( shape2 );  
   index2[i] = 1;
   isuser2[i] = 0;
   if (nshape2[i]->type==regMASK || nshape2[i]->type==regUSER)
     isuser2[i] = 1;
   shape2 = shape2->next;
 }

 /* Intersect Set 1 with Set 2 - reducing Shapes as possible */
 /* Turns flag 'OFF' for absorbed shapes..                   */
 /*   Code Change  Bug 13500  January 31, 2013               */
 /*   Stop looping if any intersection fails                 */
 for ( i = 0; bIntersectOK && (i < n1); i++ )
 {
   for ( j = 0; bIntersectOK && (j < n2); j++ )
   {
    bIntersectOK = regShapeIntersect( nshape1[i], nshape2[j], &index1[i], &index2[j] );
   }
 }


 /* Now patch things together with positive regions first */
 /*  only if intersection is successful                   */

 if ( bIntersectOK )
 {

  /* User defined shapes from Set 1 which are still 'ON'  */
  /* NOTE: This will do INCLUDE and EXCLUDE, which is     */
  /*       contrary to the other shape types.. OK?        */
  for ( i = 0; i < n1; i++ )
  {
    if ( index1[i] && isuser1[i] )
    {
      regAddShape( region, glue, nshape1[i] );
      glue = regAND;
    }
  }

  /* User defined shapes from Set 2 which are still 'ON'  */
  /* NOTE: This will do INCLUDE and EXCLUDE, which is     */
  /*       contrary to the other shape types.. OK?        */
  for ( i = 0; i < n2; i++ )
  {
    if ( index2[i] && isuser2[i] )
    {
      regAddShape( region, glue, nshape2[i] );
      glue = regAND;
    }
  }

  /* Include shapes from Set 1 which are still 'ON' */
  for ( i = 0; i < n1; i++ )
  {
    if ( index1[i] && nshape1[i]->include==regInclude && !isuser1[i] )
    {
      regAddShape( region, glue, nshape1[i] );
      glue = regAND;
    }
  }

  /* Include shapes from Set 2 which are still 'ON' */
  for ( i = 0; i < n2; i++ )
  {
    if ( index2[i] && nshape2[i]->include == regInclude && !isuser2[i] )
    {
      regAddShape( region, glue, nshape2[i] );
      glue = regAND;
    }
  }

  /* Exclude shapes from Set 1 which are still 'ON' */
  for ( i = 0; i < n1; i++ )
  {
    if ( index1[i] && nshape1[i]->include != regInclude && !isuser1[i] )
    {
      if ( glue == regAND )
	regAddShape( region, glue, nshape1[i] );
    }
  }

  /* Exclude shapes from Set 2 which are still 'ON' */
  for ( i = 0; i < n2; i++ )
  {
    if ( index2[i] && nshape2[i]->include != regInclude && !isuser2[i] )
    {
      if ( glue == regAND )
	regAddShape( region, glue, nshape2[i] );
    }
  }

 } /* end bIntetsectOK */
 

 /* clean up */
 for ( i = 0; i < n1; i++ )
 {
   if ( !index1[i] || (0 == bIntersectOK) )
     regFreeShape( NULL, nshape1[i] );
 }
 
 for ( i = 0; i < n2; i++ )
 {
   if ( !index2[i] || (0 == bIntersectOK) )
     regFreeShape( NULL, nshape2[i] );
 }
 

 free( index1 );
 free( index2 );
 free( nshape1 );
 free( nshape2 );
 free( isuser1 );
 free( isuser2 );

 return bIntersectOK;
}



int regShapeIntersect( regShape* shape1, regShape* shape2, long* index1, long* index2 )
{
 double fx[2] = { -DBL_MAX, DBL_MAX };
 double fy[2] = { -DBL_MAX, DBL_MAX };
 double xpos1[2];
 double ypos1[2];
 double xpos2[2];
 double ypos2[2];
 int ok = 1;

 if ( !shape1 || !shape2 )
   return 0;


 regExtentShapeRaw( shape1, fx, fy, xpos1, ypos1 );
 regExtentShapeRaw( shape2, fx, fy, xpos2, ypos2 );

 ok = regRectangleOverlap( xpos1, ypos1, xpos2, ypos2 );

 if ( ok )
 {

  /* If shapes are identical, just return one of them */
  if ( regCompareShape( shape1, shape2 ))
  {
    *index2 = 0;
    return ok;
  }
  else if ( shape1->include != shape2->include )
  {
   /* The regions are exact opposites, they cancel out
   *  This often happens with bounding areas. */
   if ( regCompareShapeRaw( shape1, shape2 ))
   {
 
     *index1 = 0;
     *index2 = 0;
     ok = 0;
     
     return ok;
   }
  }

  /* Exclusions always overlap */
  if ( shape1->include != regInclude || shape2->include != regInclude ) 
   return ok;

  /* The temporary hack to support only rectangles */
  /* If shape1 is entirely within rectangle 2, ignore rectangle 2. */
  if ( regIsRect( shape2 ))
  {
   if ( regRectangleInside( xpos1, ypos1, xpos2, ypos2 ))
   {
    *index2 = 0;
   }
   else if ( regIsRect( shape1 ))
   {
    /* Replace shape1 with intersected rectangle */
    *index2 = 0;
    if ( shape2->xpos[0] > shape1->xpos[0] ) shape1->xpos[0] = shape2->xpos[0];
    if ( shape2->xpos[1] < shape1->xpos[1] ) shape1->xpos[1] = shape2->xpos[1];
    if ( shape2->ypos[0] > shape1->ypos[0] ) shape1->ypos[0] = shape2->ypos[0];
    if ( shape2->ypos[1] < shape1->ypos[1] ) shape1->ypos[1] = shape2->ypos[1];    
   }
  }
  else if (  regIsRect( shape1 ) && regRectangleInside( xpos2, ypos2, xpos1, ypos1 ))
  {
   *index1 = 0;
  }
  return ok;
 }

 ok = 1;
 if ( shape1->include != regInclude )
 {
  if ( shape2->include != regInclude ) 
   return 1;

   /* Exclude region outside included area of interest - ignore */
   *index1 = 0;
 }
 else if (  shape2->include != regInclude )
 {
   /* Exclude region outside included area of interest - ignore */
   *index2 = 0;
 }
 else
 {
   *index1 = 0;
   *index2 = 0; 
  /*
   * In this case, two included regions have no overlap. 
   * This is the case in which we throw out the whole intersection as null.
   */
   ok = 0;
 }

 return ok;
}

int regIsRect( regShape* shape )
{
 int ok;
 ok = ( shape->type == regRECTANGLE && shape->angle[0] == 0.0 );
 return ok;
}


int regRectangleInside( double* xpos1, double* ypos1, double* xpos2, double* ypos2 )
{
 int ok
  = ( xpos1[0] >= xpos2[0] && xpos1[0] <= xpos2[1] &&
      xpos1[1] >= xpos2[0] && xpos1[1] <= xpos2[1] &&
      ypos1[0] >= ypos2[0] && ypos1[0] <= ypos2[1] &&
      ypos1[1] >= ypos2[0] && ypos1[1] <= ypos2[1] );
 return ok;
}


regRegion* regIntersectRegion( regRegion* Region1, regRegion* Region2 )
{
 double fx[2] = { -DBL_MAX, DBL_MAX };
 double fy[2] = { -DBL_MAX, DBL_MAX };
 long n1 = 0;
 long n2 = 0;
 regRegion *Region = NULL;
 regShape  *Shape1 = NULL;
 regShape  *Shape2 = NULL;


 /* Copy shapes */
 if ( !Region1 )
   return regCopyRegion( Region2 );

 if (!Region2 )
   return regCopyRegion( Region1 );

 /* If regions are equal, just return copy of one */

 if ( regCompareRegion( Region1, Region2 ))
   return regCopyRegion( Region1 );

 Region = regCreateEmptyRegion();
 Shape1 = Region1->shape;

 while (Shape1 != NULL ) /* Loop over components of Region 1 */
 {
  ++n1;
  Shape2 = Region2->shape; 
 
  while( Shape2 != NULL ) /* Loop over components of Region 2 */
  {
    ++n2;
    regIntersectComponent( Region, Shape1, Shape2 );
    Shape2 = regNextComponent( Shape2 );
  }
  Shape1 = regNextComponent( Shape1 );
 }
 
 regExtent(Region, fx, fy, Region->xregbounds, Region->yregbounds);
  
 return Region;
}


int regShapeOverlap( regShape* shape1, regShape* shape2 )
{
 double fx[2] = { -DBL_MAX, DBL_MAX };
 double fy[2] = { -DBL_MAX, DBL_MAX };
 double xpos1[2];
 double ypos1[2];
 double xpos2[2];
 double ypos2[2];
 int ok = 1;

 if ( !shape1 || !shape2 )
   return 0;

/* Exclusions always overlap */
 if ( shape1->include != regInclude || shape2->include != regInclude ) 
   return 1;

 regExtentShape( shape1, fx, fy, xpos1, ypos1 );
 regExtentShape( shape2, fx, fy, xpos2, ypos2 );
 ok = regRectangleOverlap( xpos1, ypos1, xpos2, ypos2 );
 return ok;
}

int regRectangleOverlap( double* xpos1, double* ypos1, double* xpos2, double* ypos2 )
{
 if (   xpos1[1] < xpos2[0] 
     || xpos1[0] > xpos2[1] 
     || ypos1[1] < ypos2[0] 
     || ypos1[0] > ypos2[1] )
     return 0;
 return 1;
}

/*
 *  Old interfaces
 */

regShape *regCreateShape(
			 regRegion *region,
			 regMath glue,
			 regGeometry shape,
			 regFlavor include,
			 double *xpos,
			 double *ypos,
			 long   npoints,
			 double *radius,
			 double *angle
			 )
{
  return regCreateWorldShape( region, glue, shape, include, xpos, ypos, npoints, radius, angle, 0, 0 );
}



regShape *regCreateNewShape(
			 regGeometry shape,
			 regFlavor include,
			 double *xpos,
			 double *ypos,
			 long   npoints,
			 double *radius,
			 double *angle
			 )
{
  return( regCreateWorldShape( NULL, regOR, shape, include, xpos,
			  ypos, npoints, radius, angle, 0, 0 ) );
}


/*
 *   Number of radii for this shape
 */
long regShapeRadii( const regShape* Shape )
{
 long nradii;
  switch( Shape->type )
  {
    case regCIRCLE:
      nradii = 1;
      break;
    case regELLIPSE:
    case regANNULUS:
    case regBOX:
    case regPIE:
      nradii = 2;
      break;
    case regPOINT:
    case regRECTANGLE:
    case regPOLYGON:
    case regSECTOR:
    case regMASK:
    case regUSER:
    default:
      nradii = 0;
      break;
  }
 return nradii;
}


long regShapeAngles( const regShape* Shape )
{
 long nangles;
  switch( Shape->type )
  {
    case regFIELD:
    case regCIRCLE:
    case regANNULUS:
    case regPOINT:
    case regPOLYGON:
      nangles = 0;
      break;
    case regELLIPSE:
    case regBOX:
    case regRECTANGLE:
      nangles = 1;
      break;
    case regPIE:
    case regSECTOR:
      nangles = 2;
      break;
    case regMASK:
    case regUSER:
    default:
      nangles = 0;
      break;
  }
 return nangles;
}


