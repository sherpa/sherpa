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

#include "cxcregion.h"
#include "region_priv.h"
#include <float.h>

char* reg_get_flag_name(int coord);

// *** PUBLIC FUNCTIONS *** //


  

void regAppendShape( regRegion* region, 
		     char* shapeName, 
		     int includeFlag, 
		     int orFlag, 
		     double* xpos, 
		     double* ypos,
		     long nmax, 
		     double* radius, 
		     double* angle, 
		     int flag_coord, 
		     int world_size )
{

  regMath      glue;
  regGeometry  type;
  regFlavor    include;
  long         npoints =0;
  double       fx[2] ={ -DBL_MAX, DBL_MAX };
  double       fy[2] ={ -DBL_MAX, DBL_MAX };
  
  glue    = orFlag ? regOR : regAND;
  include = includeFlag? regInclude : regExclude;
  
  if ( !strcmp( shapeName, "npolygon" ) || !strcmp( shapeName, "NPOLYGON" )) {
    type = regPOLYGON;
    npoints = nmax;
  } else {
    type = reg_shape_name_to_geometry( shapeName ); 
    npoints = reg_shape_find_npoints( type, xpos, ypos, nmax ); 
  }

  if ( type == regMASK ){ return;}

  regCreateShape( region, glue, type, include, xpos, ypos, npoints,
		  radius, angle, flag_coord, world_size );


  regExtent(region, fx, fy, region->xregbounds, region->yregbounds);


}



/* --------------------------------------------------------------------------- */
/* regCopyShape                                                                */
/*   Create a new shape with same definition as the old one.                   */
/*   NOTE:                                                                     */
/*     component, region, next attributes are NOT copied                       */
/* --------------------------------------------------------------------------- */
regShape* regCopyShape( regShape* inShape )
{
    if (!inShape) {
        return NULL;
    }
    
    regShape *shape = NULL;
    shape = inShape->copy(inShape);
    return shape;
}



/* -----------------------------------------------------------------------

 * Returns whether or not the shapes are equal. Inclusion can be ignored
 * in the equality check by including the raw flag. 

   ----------------------------------------------------------------------- */
int regCompareShape( regShape* Shape1, regShape* Shape2, short raw)
{
    // If we care about the inclusion flag
    if (!raw) {
        return Shape1->isEqual(Shape1, Shape2);
    }

    // Otherwise
    short ret;
    regShape *copyShape; 

    // Make a copy of the first shape and switch the flag.
    copyShape = regCopyShape(Shape1);
    copyShape->include = copyShape->include == regInclude ? regExclude : regInclude;

    // If one of them matches then the shapes are the same "shape".
    ret = Shape1->isEqual(Shape1, Shape2) || 
          copyShape->isEqual(copyShape, Shape2); 
    return ret;
}


/*
  +-----------------------------------------------------------------------
  + Simple wrapper to return if a point is inside an individual shape
  +-----------------------------------------------------------------------
*/
int regInsideShape(
		   regShape *shape,
		   double xpos,
		   double ypos
		   )
{
  return ( shape->isInside( shape, xpos, ypos ) );
}


/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */
void regPrintShape(regShape * shape)
{

    if (shape == NULL) {
      printf("ERROR: Input shape is NULL\n");
      return;
    }
    
    char wbuf[80];
    
    // In case of very long polygons

    //printf("regPrintShape() - set buffer size \n"); 
    long size = shape->nPoints > 2 ? (80 + (20*shape->nPoints)) : 120;
    char sbuf[size];

    if(!shape) {
        return;
    }

    printf("%ld\t", shape->component);       
 
    snprintf(wbuf, 80, "(Flag_Coord: %s) (Flag_Radius: %s)", 
            reg_get_flag_name(shape->flag_coord),
            reg_get_flag_name(shape->flag_radius));
    
    shape->toString(shape, sbuf, size);

    printf("%s %s\n", sbuf, wbuf);

}

char* reg_get_flag_name(int coord) {
    
    char *coordunit[] = { "Unknown", "Logical", "Physical", "World" };
    
    if (coord > 3 || coord < 0) {
        return coordunit[0];
    }

    return coordunit[coord];
}



/* Return name of shape. Also return whether included or excluded,
  as function return value
 */

int regShapeGetName(const regShape * shape, char *name, long maxlen)
{

  if ( NULL == shape ) {
    strncpy(name, "Unknown", maxlen);
    return(0);
  }

    strcpy(name, "");
    switch (shape->type) {
    case regCIRCLE:
	strncpy(name, "Circle", maxlen);
	break;
    case regPOINT:
	strncpy(name, "Point", maxlen);
	break;
    case regELLIPSE:
	strncpy(name, "Ellipse", maxlen);
	break;
    case regPIE:
	strncpy(name, "Pie", maxlen);
	break;
    case regSECTOR:
	strncpy(name, "Sector", maxlen);
	break;
    case regANNULUS:
	strncpy(name, "Annulus", maxlen);
	break;
    case regPOLYGON:
	strncpy(name, "Polygon", maxlen);
	break;
    case regBOX:
	strncpy(name, "Box", maxlen);
	if (shape->angle[0] != 0.0)
	    strncpy(name, "RotBox", maxlen);

	break;
    case regRECTANGLE:
	strncpy(name, "Rectangle", maxlen);
	if (shape->angle[0] != 0.0)
	    strncpy(name, "RotRectangle", maxlen);
	break;
    case regFIELD:
	strncpy(name, "Field", maxlen);
	break;
    case regMASK:
       strncpy(name, shape->name, maxlen);
       break;
    default:
	strncpy(name, "Unknown", maxlen);
	break;
    }
    if (shape->include == regExclude)
	return 0;
    else
	return 1;
}


/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */
long regShapeGetPoints(const regShape * shape, double *x, double *y,
		       long dim)
{
    long n;
    long i;
    if (!shape || !x || !y || dim <= 0)
	return 0;

    if (!shape->xpos || !shape->ypos || shape->nPoints <= 0)
	return 0;

    n = shape->nPoints;
/* Only return as many points as are in the array */
    if (n > dim)
	n = dim;

    for (i = 0; i < n; i++) {
	x[i] = shape->xpos[i];
	y[i] = shape->ypos[i];
    }
/* Zero any remaining array elements */
    for (i = n; i < dim; i++) {
	x[i] = 0.0;
	y[i] = 0.0;
    }
    return n;
}


/*
 *  Return the values of radius in the preallocated array.
 *  Return the number of radii as the function value.
 *  nradii can be 0, 1 or 2.
 */
long regShapeGetRadii(const regShape * shape, double *r)
{
    long nradii;
    long i;
    if (!shape || !r)
	return 0;

    nradii = reg_shape_radii(shape);
    for (i = 0; i < nradii; i++)
	r[i] = shape->radius[i];
    return nradii;
}

/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */

long regShapeGetAngles(const regShape * shape, double *angle)
{
    long nangles;
    long i;
    if (!shape || !angle)
	return 0;

    nangles = reg_shape_angles(shape);
    for (i = 0; i < nangles; i++)
	angle[i] = shape->angle[i];
    return nangles;
}

/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */

long regShapeGetNoPoints(const regShape * shape)
{
    if (!shape)
	return 0;

    return shape->nPoints;
}

/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */

long regShapeGetComponent(const regShape * shape)
{
    if (!shape)
	return 0;

    return (long) shape->component;
}


// ***  PRIVATE FUNCTIONS *** //

/* --------------------------------------------------------------------------- */
/* regCreateNewWorldShape                                                      */
/*   Create and populate Shape object with provided inputs.                    */
/*                                                                             */
/*   NOTE: not appropriate for user defined shapes ( regUSER,regMASK)          */
/* --------------------------------------------------------------------------- */
regShape* regCreateNewWorldShape(
			 regGeometry type,
			 regFlavor include,
			 double *xpos,
			 double *ypos,
			 long   npoints,
			 double *radius,
			 double *angle,
                         int flag_coord,    /* xpos ypos in degrees or pixels? */
                         int world_size      /* radius in degrees or pixels? */
			 )
{
  return( regCreateShape( NULL, regOR, type, include, xpos, ypos, npoints, 
			  radius, angle, flag_coord, world_size ) );
}


/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */
regShape *regCreateShape(
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
  double fx[2] ={ -DBL_MAX, DBL_MAX };
  double fy[2] ={ -DBL_MAX, DBL_MAX };

  /* Region library can not create Pixel Mask shapes */
  if ( type == regMASK ) { return newShape; };


  switch ( type )
    {
    case regPOINT:
      newShape = regCreatePoint(include, xpos, ypos, wcoord, wsize); 
      break;
    
    case regCIRCLE:
      newShape = regCreateCircle(include, xpos, ypos, radius, wcoord, wsize); 
      break;
 
    case regANNULUS:
      newShape = regCreateAnnulus(include, xpos, ypos, radius, wcoord, wsize);
      break;
 
    case regELLIPSE:
      newShape = regCreateEllipse(include, xpos, ypos, radius, angle, wcoord, wsize); 
      break;
         
    case regBOX:
      newShape = regCreateBox(include, xpos, ypos, radius, angle, wcoord, wsize); 
      break;

    case regRECTANGLE:
      newShape = regCreateRectangle(include, xpos, ypos, angle, wcoord, wsize); 
      break;

    case regPOLYGON:
      newShape = regCreatePolygon(include, xpos, ypos, npoints, wcoord, wsize);
      break;

    case regPIE:
      newShape = regCreatePie(include, xpos, ypos, radius, angle, wcoord, wsize);
      break;

    case regSECTOR:
      newShape = regCreateSector(include, xpos, ypos, angle, wcoord, wsize);
      break;

    case regFIELD:
      newShape = regCreateField(include, wcoord, wsize);
      break;

   default:
      return(NULL);
    }

  /* attributes get assigned/init when glued to a region */

  if ( newShape == NULL ) return NULL; 

  if ( region != NULL ) {
      regAddShape( region, glue, newShape );
      regExtent(region, fx, fy, region->xregbounds, region->yregbounds);
  }
 
  return( newShape );
}


/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */

long regAddShape( regRegion *region,
		  regMath glue,
		  regShape *shape )
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


/* -----------------------------------------------------------------------

  
   ----------------------------------------------------------------------- */

void regFreeShape( regRegion* region, regShape* atShape )
{
  if( !atShape )
    return;

  if ( atShape->xpos ) free( atShape->xpos );
  if ( atShape->ypos ) free( atShape->ypos );
  if ( atShape->angle ) free( atShape->angle );
  if ( atShape->radius ) free( atShape->radius );
  if ( atShape->sin_theta) free( atShape->sin_theta);
  if ( atShape->cos_theta) free( atShape->cos_theta);

  free( atShape );    
  atShape=NULL;

}
