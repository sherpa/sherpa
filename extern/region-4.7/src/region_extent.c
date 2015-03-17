/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include <float.h>
#include <math.h>
#include "region_priv.h"


#define REGION_MAX_PIXEL_SIZE 0.5

double regAreaPolygon(regShape * shape);
double regBoundsArea(double *fieldx, double *fieldy);
long reg_quadrant(double ang1);
void reg_pie_bounds(double ang1, double ang2, double r2, double r1,
		    double *xoff, double *yoff);
void regSetField(regShape * shape, double *fieldx, double *fieldy);
double regAreaPie(regShape * shape);
int regExtentShapeRaw(regShape * shape, double *fieldx, double *fieldy,
		      double *xpos, double *ypos);
void regSetBounds(double *xpos, double *ypos, double *fieldx,
		  double *fieldy);
double reg_dtheta(double theta, double ang);
double reg_closest_theta(double theta, double ang1, double ang2);
double regComputePixellatedArea(regRegion * region, double *xbounds,
				double *ybounds, double bin);
int regZeroBounds(double *xpos, double *ypos);
int regTrimExtent(double *cxpos, double *cypos, double *sxpos,
		  double *sypos, int cstart);
int regUnionExtent(double *cxpos, double *cypos, double *sxpos,
		   double *sypos, int cstart);
int regBreakdown(double *shape1_x, double *shape1_y,
		 double *shape2_x, double *shape2_y,
		 double *xcoords, double *ycoords);
void setFlags(long num_shapes, regRegion * region, int *intersecting, 
	     double *xpos, double *ypos);
double sumAreas(long num_shapes, regRegion * region, int *intersecting, int *compno,
	int *includes, double *areas, double *xpos, double *ypos, double lbin,
	int *fully_in_field, double *fieldx, double *fieldy, int *union_trim);

double regAreaBox(regShape * shape)
{
    double area;
    area = shape->radius[1] * shape->radius[0];
    return area;
}

double regAreaRectangle(regShape * shape)
{
    double area;
    double xr, yr;
    regRectangleSides(shape, &xr, &yr);
    area = xr * yr;
    return area;
}

double regAreaCircle(regShape * shape)
{
    double area;
    area = PI * shape->radius[0] * shape->radius[0];
    return area;
}

double regAreaEllipse(regShape * shape)
{
    double area;
    area = PI * shape->radius[0] * shape->radius[1];
    return area;
}

double regAreaPie(regShape * shape)
{
    double area;
    double ang1, ang2, theta;
    ang1 = regModAngle(shape->angle[0]);
    ang2 = regModAngle(shape->angle[1]);
    if (ang1 < ang2)
	theta = ang2 - ang1;
    else
	theta = 360.0 - (ang1 - ang2);
    area =
	PI * (theta / 360.0) * (shape->radius[1] * shape->radius[1] -
				shape->radius[0] * shape->radius[0]);
    return area;
}

double regAreaPolygon(regShape * shape)
{

    double area = 0.0;
    long ii;

    /* last point is already there */
    /* eqn is 0.5*sum( x_i*y+i_1 - x_i+1*y_i) */
    for (ii=0;ii<(shape->nPoints -1);ii++) {
      area += ((shape->xpos[ii]*shape->ypos[ii+1]) -
	       (shape->xpos[ii+1]*shape->ypos[ii]) );
    }
    area /= 2.0;

    return fabs(area);
}

double
regShapeAnalyticArea(regShape * shape, double field_area, int *status)
{
    double area = 0.0;

    *status = 1;
    switch (shape->type) {
    case regELLIPSE:
	area = regAreaEllipse(shape);
	break;
    case regPOINT:
	area = 0;
	break;
    case regCIRCLE:
	area = regAreaCircle(shape);
	break;
    case regBOX:
	area = regAreaBox(shape);
	break;
    case regRECTANGLE:
	area = regAreaRectangle(shape);
	break;
    case regFIELD:
	area = field_area;
	break;
    case regANNULUS:
	area =
	    PI * (shape->radius[1] * shape->radius[1] -
		  shape->radius[0] * shape->radius[0]);
	break;
    case regPIE:
	area = regAreaPie(shape);
	break;
    case regPOLYGON:		/* Too complicated */
        area = regAreaPolygon(shape);
        break;
    case regMASK:
    case regUSER:
      area=shape->user->calcArea(shape);
      break;	
    case regSECTOR:		/* Overlaps field always */
    default:
      *status = 0;
      area = 1.0;
      break;
    }
    if (shape->include == regExclude)
	area = -area;
    
    return area;
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regExtentPolygon(regShape * shape, double *xpos, double *ypos)
{
    long ii;
    xpos[0] = shape->xpos[0];
    xpos[1] = shape->xpos[0];
    ypos[0] = shape->ypos[0];
    ypos[1] = shape->ypos[0];

    for (ii = 1; ii < shape->nPoints - 1; ii++) {
	if (shape->xpos[ii] < xpos[0])
	    xpos[0] = shape->xpos[ii];
	if (shape->xpos[ii] > xpos[1])
	    xpos[1] = shape->xpos[ii];
	if (shape->ypos[ii] < ypos[0])
	    ypos[0] = shape->ypos[ii];
	if (shape->ypos[ii] > ypos[1])
	    ypos[1] = shape->ypos[ii];
    }

    return 1;
}

double regBoundsArea(double *fieldx, double *fieldy)
{
    double field_area;
    if (fieldx[1] >= DBL_MAX || fieldx[0] <= -DBL_MAX ||
	fieldy[1] >= DBL_MAX || fieldy[0] <= -DBL_MAX)
	field_area = DBL_MAX;
    else
	field_area = (fieldx[1] - fieldx[0]) * (fieldy[1] - fieldy[0]);
    return field_area;
}

double
regArea(regRegion * region, double *fieldx, double *fieldy, double lbin)
{
  int ii=0;           /* counter for variable initialization */
  double area =0;       /*  area of a single analytic shape */
  double *xpos = NULL;  /* array of shape x positions,
			   lower left & upper right */
  double *ypos = NULL;  /* array of shape y positions,
			   lower left & upper right */
  double *x_segments = NULL;  /* x-position info about segments of
			         overlaping bounding boxes for shapes */
  double *y_segments = NULL;  /* y-position info about segments of
			         overlaping bounding boxes for shapes */
  int status = 0;
  long num_shapes = 0;   /* holds the total number of shapes in region */
  regShape *atShape = NULL;  /* pointer to shape when only one exists */
  int within_field =0;   /* is the region within the field boundaries */
  int *fully_in_field;  /* is each shape within the field boundaries */
  double field_area;  /* total area of the field */
  double total_area = 0.0;  /* total area of the region shapes */
  double *xregbounds = NULL;
  double *yregbounds = NULL;
  double *areas = NULL;  /* this list holds the areas of each shape */
  int *compno = NULL;      /* shape is defined to intersect with another */
  int *intersecting = NULL; /* intersecting indicate if shape really intersects w/ another */
  int *includes = NULL;      /* is the shape regIncluded? */
  long shape_no =0;            /* holds the current shape number */
  regShape *shape = NULL;    /* pointer to current shape */
  long offset =0; /* this is an array item number to begin coords from */
  int idx = 0;   /* shape index , shape_no-1 */
  int* union_trim;  /* 1=unioned 2= trimmed */

/* Calculate field area */

  field_area = regBoundsArea(fieldx, fieldy);
  if (!region)
    return field_area;
  
  xregbounds = calloc(2, sizeof(double));
  yregbounds = calloc(2, sizeof(double));
  
  atShape = region->shape;
  if (atShape == NULL)
    return 0.0;

  /* regExtent returns 1 if all defined shapes fit within field boundaries */
  within_field = regExtent(region, fieldx, fieldy, xregbounds, yregbounds);

  if (atShape->next == NULL && within_field) {
    /* Only one shape; use analytic method for speed? */
    area = regShapeAnalyticArea(atShape, field_area, &status);
    if (area < 0.0)
      area = field_area + area;
    if (status)
      return area;
  }
  
  num_shapes = regGetNoShapes(region);
  areas= calloc(num_shapes,sizeof(double)); /* areas of individual shapes */
  intersecting= calloc(num_shapes,sizeof(int)); /*represents shapes that intersect*/
  compno = calloc(num_shapes, sizeof(int));
  fully_in_field = calloc(num_shapes, sizeof(int));
  includes = calloc(num_shapes, sizeof(int));
  xpos = calloc(2*num_shapes, sizeof(double)); /*shape bound box x-coords*/
  ypos = calloc(2*num_shapes, sizeof(double)); /* shape bound box y-coords */
  x_segments = calloc(14, sizeof(double)); /* x-position info of segments */
  y_segments = calloc(14, sizeof(double)); /* y-position info of segments */
  union_trim = calloc(num_shapes, sizeof(int));

  for (ii=0; ii<num_shapes; ii++)
    {
      union_trim[ii] = 0;
      intersecting[ii] = 0;
      areas[ii] = 0;
      compno[ii] = 0;
      includes[ii] = 0;
      fully_in_field[ii] = 0;
    }

/* Pass 1: determine extent and overlap */

/* We divide up the field into non-overlapping rectangles
   containing the shapes.
   It doesn't matter whether the shapes overlap in intersection or in union:
   any time the shapes overlap we mark off that area as an overlap group
   and compute it by
   pixellating the area and calling regInsideRegion for each pixel.
   If we have a shape that is isolated, we use an analytic expression for
   the shape - but we don't have analytic expressions for sectors (which
   always overlap the field boundaries) and polygons, so we use the pixellated
   approach by making them an overlap group of size 1 (we set the intersecting[])
   field.)

   This algorithm could be improved by adding more special cases of
   intersections we know how to do analytically.
 */

  /* for each shape in the region */
  for (shape_no = 1; shape_no <= num_shapes; shape_no++) {
    idx = shape_no - 1;  /* shape index */
    
    /* get the shape */
    shape = regGetShapeNo(region, shape_no);

    if ( regInclude == shape->include )
      includes[idx] = 1;
    
    /* if there is another shape after this one and its component
       number is the same (ie. it is defined to intersect with the current
       shape) then set this compno flag for current shape
       else
       set the flag if the current component matches the previous shape
       component. (needed for last shape in the region list)
    */
    compno[idx] = shape->component;
    
    /* set a pointer to the shape coordinates */
    offset = 2 * (idx);

    /* get the bounding box of the current shape */
    regExtentShapeRaw(shape, fieldx, fieldy,
		      xpos + offset, ypos + offset);

    /* trim the shape boundaries to the field boundaries */
    if ( regInclude == shape->include )
      {
	fully_in_field[idx] = regTrimExtent(xpos + offset, ypos + offset,
				   fieldx, fieldy, 0);
      }
    
    /* get the analytic area of each individual shape */
    status = 0;
    areas[idx] = regShapeAnalyticArea(shape, field_area, &status);
    if (areas[idx] < 0)
      areas[idx] = 0;
  }

  /* On each iteration, check to see if any bounding boxes overlap.
     If they do, expand a bounding box to include overlapping shapes.
     Keep doing this until there are no more overlapping bounding boxes,
     then add up the areas of the shapes in each bounding box.
  */
 do {
    
    setFlags( num_shapes, region, intersecting, xpos, ypos);

    total_area = sumAreas(num_shapes,region,intersecting,compno,includes,
			  areas,xpos,ypos,lbin,fully_in_field,fieldx,fieldy,union_trim);
    
  } while (total_area == -1);


  free(xregbounds);
  free(yregbounds);
  free(areas);
  free(intersecting);
  free(xpos);
  free(ypos);
  free(x_segments);
  free(y_segments);

  return total_area;
}

/* indicate through the intersecting array which shapes truly intersect with each
   other by setting them to the same number (which is the first shape in
   the region string that other shapes intersect with).  ie. if the idx=0,2,
   and 4 shapes intersect, then intersecting[0] = intersecting[2] = intersecting[4] = 1;
*/
void setFlags(long num_shapes, regRegion * region, int *intersecting,
	      double *xpos, double *ypos)
{
  long shape_no =0;            /* holds the current shape number */
  long prev_no =0;       /* holds the secondary shape number */
  long offset =0; /* this is an array item number to begin coords from */
  long prev_offset =0; /* pointer to coord values for PREVIOUS shape */
                       /* PREVIOUS can be any shape before the current one */
  int idx = 0;   /* shape index , shape_no-1 */
  int pidx = 0;   /* prev. shape index , prev_no-1 */


  /* for each shape in the region
     clear the intersect flags
  */
  for (shape_no = 1; shape_no <= num_shapes; shape_no++)
    {
      idx = shape_no - 1;
      intersecting[idx] = 0;
    }
  
  /* for each shape in the region */
  for (shape_no = 1; shape_no <= num_shapes; shape_no++)
    {
      idx = shape_no - 1;
      
      /* Does this shape overlap with any other (previous shape)? */
      offset = 2 * (idx);  /* set pointer shape bound box coords */
    
      for (prev_no = 1; prev_no < shape_no; prev_no++)
	{
	  pidx = prev_no - 1;
	  
	  /* pointer to PREVIOUS shape coords */
	  prev_offset = 2 * (pidx);

	  /* if this shape overlaps with any others */
	  if (regRectangleOverlap(xpos + offset, ypos + offset,
				  xpos + prev_offset, ypos + prev_offset))
	    {
	      /* if a flag wasn't already set for a PREVIOUS shape,
		 then set it */
	      if (!intersecting[pidx])
		intersecting[pidx] = prev_no;
	      
	      /* set a flag for the CURRENT shape to the first previous
		 shape from the beginning of the region string that it
		 intersects with */
	      intersecting[idx] = prev_no;
	      break;
	    }
	}

    } /* end loop for shape_no */
}

/*
 *  U(CPT1 CPT2 CPT3)...
 *  CPT1 = S11*S12*...S1n
 *  Area = union of cpt areas
 */
double sumAreas(long num_shapes, regRegion * region, int *intersecting, int *compno,
     int *includes, double *areas, double *xpos, double *ypos, double bin,
     int *fully_in_field, double *fieldx, double *fieldy, int *union_trim)
{
  double *c_xpos = NULL;  /* array of component x positions,
			   lower left & upper right */
  double *c_ypos = NULL;  /* array of component y positions,
			   lower left & upper right */
  double total_area = 0.0;  /* total area of the region shapes */
  long shape_no =0;            /* holds the current shape number */
  long other_no =0;       /* holds the secondary shape number */
  long offset =0; /* this is an array item number to begin coords from */
  long coffset =0; /* this is an array item number to begin coords from */
  long other_offset =0; /* pointer to coord values for OTHER shape */
                       /* OTHER can be any shape after the current one */
  int idx = 0;   /* shape index , shape_no-1 */
  int cidx = 0;   /* shape index , shape_no-1 */
  int oidx = 0;   /* other. shape index , other_no-1 */
  int cno =0;               /* component counter */
  int o_cno =0;               /* other component counter */
  int jj =0;               /* counter */
  int num_components =0;   /* total number of components */
  double *comp_areas = NULL;  /* this list holds the areas of each component */
  int dont_compute =0;
  int trim_ret =1;
  int trimmed =0;          /* flag if component shapes were trimmed */
  int unioned =0;          /* flag if components were unioned */
  int* area_done;          /* did area of component x? */
  int* must_pixellate;     /* included shape & excluded shape intersect? */

  total_area = 0.0;

 /* determine the number of components in the list */
  for (shape_no = 1; shape_no <= num_shapes; shape_no++) {

    idx = shape_no - 1;  /* define the index */

    if (compno[idx] > num_components)
      num_components = compno[idx];
  }

  c_xpos = calloc(2*num_components, sizeof(double));
  c_ypos = calloc(2*num_components, sizeof(double)); 
  comp_areas = calloc(num_components+1, sizeof(double));
  area_done = calloc(num_components+1, sizeof(int));
  must_pixellate = calloc(num_components+1, sizeof(int));

  /* initialize */
  for (cno=1; cno <= num_components; cno++)
    {
      comp_areas[ cno ] =0;
      area_done[ cno ] =0;
      must_pixellate[ cno ] =0;
    }
  /* shrink the bounding boxes around each component */
  for (cno=1; cno <= num_components; cno++)
    {
      cidx = cno - 1;  /* define the index */
      coffset = 2 * (cidx); /* pointer to shape bound box coords */

      /* look at each shape on the list.  we'll call each current shape
	 the PRIMARY shape
      */
      for  (shape_no = 1; shape_no <= num_shapes; shape_no++)
	{
	  idx = shape_no - 1;  /* define the index */
	  offset = 2 * (idx); /* pointer to shape bound box coords */

	  /* if this shape isn't part of the component,
	     OR if this shape isn't included, continue to the next shape */
	  if ( compno[idx] != cno || !includes[idx] )
	    continue;

	  /* set the component bounds to be the PRIMARY shape bounds */
	  c_xpos[ coffset ]   = xpos[ offset ];    /* x min */
	  c_ypos[ coffset ]   = ypos[ offset ];    /* y min */
	  c_xpos[ coffset+1 ] = xpos[ offset+1 ];  /* x max */
	  c_ypos[ coffset+1 ] = ypos[ offset+1 ];  /* y max */
	  
	  /* look at all OTHER shapes on the list */
	  for (other_no = 1; other_no <= num_shapes; other_no++)
	    {
	      oidx = other_no - 1;
	      other_offset = 2 * (oidx);

	      /* if the OTHER shape is excluded and it intersects
	         with the included shape, then we must pixellate the
	         included shape/component bounds
	      */
	       if ( !includes[oidx] &&
		    regRectangleOverlap(xpos + offset, ypos + offset,
					xpos + other_offset,
					ypos + other_offset) )
		 {
		   must_pixellate[ cno ] = 1;
		 }
		 
	      /* if this OTHER shape is the same as the PRIMARY shape,
	         OR if this shape isn't part of the component,
	         OR if this shape isn't included
	         continue to the next shape */
	      if ( idx == oidx || compno[oidx] != cno || !includes[ oidx ] )
		continue;

	      /* if this OTHER shape overlaps with the PRIMARY shape */
	      if (regRectangleOverlap(xpos + offset, ypos + offset,
				      xpos + other_offset,
				      ypos + other_offset))
		{
		  /* if the bounds are already the same, call it trimmed */
		  if (xpos[offset]   == xpos[other_offset] &&
		      ypos[offset]   == ypos[other_offset] &&
		      xpos[offset+1] == xpos[other_offset+1] &&
		      xpos[offset+1] == xpos[other_offset+1] )
		    trimmed=1;
		  else
		  /* trim bounds to the intersection */
		    trim_ret =
		      regTrimExtent(xpos + offset, ypos + offset,
				    xpos + other_offset,
				    ypos + other_offset, 0);
		  if( 0 == trim_ret )
		    trimmed = 1;

		  /* reset the component bounds to the intersection bounds */
		  c_xpos[ coffset ]   = xpos[ offset ];    /* x min */
		  c_ypos[ coffset ]   = ypos[ offset ];    /* y min */
		  c_xpos[ coffset+1 ] = xpos[ offset+1 ];  /* x max */
		  c_ypos[ coffset+1 ] = ypos[ offset+1 ];  /* y max */
		  
		}
	      else /* a shape didn't overlap, but it should have */
		{
		  /* implied (default) component area is 0 */
		  area_done[ cno ] = 1;
		  break;
		}
	      
	    } /* end for OTHER shape loop */
	  
	} /* end for PRIMARY shape loop */

      /* if no trimming of shapes in this component
         and area wasn't done (ie. no missing intersect)
      */
      if ( !trimmed && !area_done[ cno ] && !must_pixellate[ cno ])
	for  (idx = 0; idx < num_shapes; idx++)
	  {
	    if ( fully_in_field[ idx ] &&
		 includes[ idx ] && compno[ idx ] == cno )
	      {
		/* this is a stand-alone shape */
		comp_areas[ cno ] += areas[ idx ];
		area_done[ cno ] = 1;
	      }
	  } /* end shape LOOP */
	
    } /* end for COMPONENT loop */


  /* expand bounding boxes around any overlapping components */
  for (cno=1; cno <= num_components; cno++)
    {
      idx = cno - 1;  /* define the index */
      offset = 2 * (idx); /* pointer to shape bound box coords */
      
        for (o_cno=1; o_cno <= num_components; o_cno++)
	  {
	    oidx = o_cno - 1;
	    other_offset = 2 * (oidx);
	    if (idx == oidx) continue;

	    /* if this OTHER component intersects with the PRIMARY component*/
	    if (regRectangleOverlap(c_xpos + offset, c_ypos + offset,
				    c_xpos + other_offset,
				    c_ypos + other_offset))
	      {
		/* expand bounds to encompass each other */
		unioned =
		  regUnionExtent(c_xpos + offset, c_ypos + offset,
				 c_xpos + other_offset,
				 c_ypos + other_offset, 0);
		comp_areas[ cno ] = 0;
		area_done[ cno ] = 0;
		if ( unioned )
		  {
		    unioned = 0;  /* reset the unioned flag */
		    cno = 0;  /* every time there's a union, start over */
		    break;
		  }
	      }
	    
	  } /* end for OTHER component loop */
    } /* end for PRIMARY component loop */

  /* pixellate if necessary, per component */
  for (cno=1; cno <= num_components; cno++)
    {
      dont_compute =0;
      idx = cno - 1;  /* define the index */
      offset = 2 * (idx); /* pointer to shape bound box coords */
      
        for ( o_cno=1; o_cno <= num_components; o_cno++)
	  {
	    oidx = o_cno - 1;
	    other_offset = 2 * (oidx);
	    if (idx == oidx) continue;

	    /* if the OTHER component bounds are the same as
	       the PRIMARY component bounds and the area
	       of the OTHER comoponet was already computed,
	       then dont_compute the area of this component */
	    if ( ((c_xpos[ offset ]   == 0) &&
		  (c_ypos[ offset ]   == 0) &&
		  (c_xpos[ offset+1 ] == 0) &&
		  (c_ypos[ offset+1 ] == 0)) ||
		 ((c_xpos[ offset ]   == c_xpos[ other_offset ] &&
		   c_ypos[ offset ]   == c_ypos[ other_offset ] &&
		   c_xpos[ offset+1 ] == c_xpos[ other_offset+1 ] &&
		   c_ypos[ offset+1 ] == c_ypos[ other_offset+1 ]) &&
		  area_done[ o_cno ] ) )
	      {
		dont_compute =1;
		break;
	      }
   	  } /* end for OTHER component loop */

	if ( area_done[ cno ] || dont_compute )
	  continue;
	/* if the component bounds have a non-zero area
	   and they overlap with the field bounds
	*/
	else if ( ((c_xpos[ offset+1 ] - c_xpos[ offset ] != 0) ||
		   (c_ypos[ offset+1 ] - c_ypos[ offset ] != 0) ) &&
		  regRectangleOverlap(c_xpos + offset, c_ypos + offset,
				      fieldx,fieldy) )
	  {
	    comp_areas[ cno ] =
	      regComputePixellatedArea(region, c_xpos + offset,
				       c_ypos + offset, bin);
	    area_done[ cno ] = 1;
	  }
    } /* end for PRIMARY component loop */

  /* add up the component areas */
  for (jj=1; jj <= num_components; jj++)
    {
      total_area += comp_areas[ jj ];
    }

  if (c_xpos) free(c_xpos);
  if (c_ypos) free(c_ypos);
  if (comp_areas) free(comp_areas);
  if (area_done) free(area_done);

  return total_area;

}

/*
 *  Compute the number of pixels within a complex region bounded by
    xbounds, ybounds 
 */
double
regComputePixellatedArea(regRegion * region, double *xbounds,
			 double *ybounds, double bin)
{
/* OK, do it the hard way, subpixel by subpixel */
    double area;
    double xsize, ysize;
    double xpos, ypos;
    long ii, jj;
    long n;
    xsize = (xbounds[1] - xbounds[0]) / bin + 1;
    ysize = (ybounds[1] - ybounds[0]) / bin + 1;
    n = 0;
    /*printf("Box(%g,%g,%g,%g,0)\n",(xbounds[1]+xbounds[0])/2,
	   (ybounds[1]+ybounds[0])/2,
	   xbounds[1]-xbounds[0],ybounds[1]-ybounds[0]);*/

    for (ii = 0; ii < xsize; ii++) {
	xpos = xbounds[0] + bin * ii;
	for (jj = 0; jj < ysize; jj++) {
	    ypos = ybounds[0] + bin * jj;
	    if (regInsideRegion(region, xpos, ypos))
		n += 1;
	}
    }
    area = n * bin * bin;
    return area;
}

void regSetField(regShape * shape, double *fieldx, double *fieldy)
{
    shape->xpos = calloc(2, sizeof(double));
    shape->ypos = calloc(2, sizeof(double));
    shape->nPoints = 2;
    regSetBounds(shape->xpos, shape->ypos, fieldx, fieldy);
}

void
regSetBounds(double *xpos, double *ypos, double *fieldx, double *fieldy)
{
    xpos[0] = fieldx[0];
    xpos[1] = fieldx[1];
    ypos[0] = fieldy[0];
    ypos[1] = fieldy[1];
}

/* **********************************************************************
 *  Determine extent of a shape
 ********************************************************************* */
int regExtentShape(regShape * shape, double *fieldx, double *fieldy, double *xpos, double *ypos)
{
  int retval = 0;

   if (shape->include != regInclude) 
   {
     /* All exclusions imply the whole field */
     if ( shape->type == regMASK )
     {
       /* except MASK which is defined to be bound by the mask field. */
       retval = regExtentShapeRaw(shape, fieldx, fieldy, xpos, ypos);
     }
     else
       regSetBounds(xpos, ypos, fieldx, fieldy);

     retval = 1;
   }
   else
   {
     retval = regExtentShapeRaw(shape, fieldx, fieldy, xpos, ypos);
   }

   return retval;
}

/* **********************************************************************
 *  Gives shape extent ignoring include/exclude
 ********************************************************************* */
int regExtentShapeRaw(regShape * shape, double *fieldx, double *fieldy, double *xpos, double *ypos)
{
    double xcor[4], ycor[4];
    double ang1, ang2;

    double txpos[2];
    double typos[2];

    /* init to field */
    regSetBounds(xpos, ypos, fieldx, fieldy);

    switch (shape->type) 
    {
    case regCIRCLE:

	xpos[0] = shape->xpos[0] - shape->radius[0];
	xpos[1] = shape->xpos[0] + shape->radius[0];
	ypos[0] = shape->ypos[0] - shape->radius[0];
	ypos[1] = shape->ypos[0] + shape->radius[0];
	break;
    case regELLIPSE:
	regBoxCorners(shape, xcor, ycor);
	regCornerBounds(xcor, ycor, xpos, ypos);
	break;
    case regPOINT:
	xpos[0] = shape->xpos[0];
	xpos[1] = shape->xpos[0];
	ypos[0] = shape->ypos[0];
	ypos[1] = shape->ypos[0];
	break;
    case regBOX:
	regBoxCorners(shape, xcor, ycor);
	regCornerBounds(xcor, ycor, xpos, ypos);
	break;
    case regRECTANGLE:
	regRectangleCorners(shape, xcor, ycor);
	regCornerBounds(xcor, ycor, xpos, ypos);
	break;
    case regSECTOR:
	break;
    case regPOLYGON:
	regExtentPolygon(shape, xpos, ypos);
	break;
    case regANNULUS:
	xpos[0] = shape->xpos[0] - shape->radius[1];
	xpos[1] = shape->xpos[0] + shape->radius[1];
	ypos[0] = shape->ypos[0] - shape->radius[1];
	ypos[1] = shape->ypos[0] + shape->radius[1];
	break;
    case regPIE:
	/* The bounds of an arbitrary pie-annulus are tricky. 
	   The X bounds are given by the closest angles to 0 and 180  included in the pie. */

	ang1 = shape->angle[0];
	ang2 = shape->angle[1];
	reg_pie_bounds(ang1, ang2, shape->radius[0], shape->radius[1],
		       txpos, typos);


	xpos[0] = txpos[0] + shape->xpos[0];
	xpos[1] = txpos[1] + shape->xpos[0];
	ypos[0] = typos[0] + shape->ypos[0];
	ypos[1] = typos[1] + shape->ypos[0];

	break;
    case regFIELD:
	break;

    case regMASK:
    case regUSER:
      shape->user->calcExtent(shape,xpos,ypos);
      break;
    default:
	xpos[0] = 0.0;
	xpos[1] = 0.0;
	ypos[0] = 0.0;
	ypos[1] = 0.0;
	break;
    }
    return 1;

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

void
reg_pie_bounds(double ang1, double ang2, double r2, double r1,
	       double *xoff, double *yoff)
{
    long q1;
    long q2;
    double c1, c2;
    double s1, s2;
    double t1, t2;
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

/*

A1<A2                    A1>A2

Q1Q1  B2  A1   B1  A2    MIN MAX MIN MAX
Q1Q2  A2  A1  F90  MAX   -
Q1Q3  MIN A1  A2   MAX   -
Q1Q4  MIN C0  MIN  MAX   -
Q2Q1  -                  MIN MAX MIN C90
Q2Q2  A2  B1  B2   A1    MIN MAX MIN MAX
Q2Q3  MIN F180 A2   A1    -
Q2Q4  MIN A2  MIN  A1    -
Q3Q1  -                  A1  MAX MIN A2
Q3Q2  -                  C180 MAX MIN MAX
Q3Q3  A1  B2  A2   B1    MIN MAX MIN MAX
Q3Q4  A1  A2  MIN  F270  -
Q4Q1  -                  F0  MAX A1   A2
Q4Q2  -                  A2  MAX A1   MAX
Q4Q3  -                  MIN MAX C270 MAX
Q4Q4  B1  A2  A1   B2    MIN MAX MIN  MAX

 Cn = closest to direction N

 */

double reg_dtheta(double theta, double ang)
{
    double rang = ang - theta;
    if (rang < -180.0)
	rang += 360.0;
    if (rang > 180.0)
	rang -= 360.0;
    return rang;
}

int regZeroBounds(double *xpos, double *ypos)
{
    int retval = 1;
    xpos[0] = 0.0;
    xpos[1] = 0.0;
    ypos[0] = 0.0;
    ypos[1] = 0.0;
    return retval;
}

/* Increase cx to include sx */
int
regUnionExtent(double *cxpos, double *cypos, double *sxpos, double *sypos,
	       int cstart)
{
/* Here we are doing the union, so the ranges get bigger */
    int retval = 0;
    if (cstart || sxpos[0] < cxpos[0]) {
	cxpos[0] = sxpos[0];
	retval = 1;
    }
    if (cstart || sxpos[1] > cxpos[1]) {
	cxpos[1] = sxpos[1];
	retval = 1;
    }
    if (cstart || sypos[0] < cypos[0]) {
	cypos[0] = sypos[0];
	retval = 1;
    }
    if (cstart || sypos[1] > cypos[1]) {
	cypos[1] = sypos[1];
	retval = 1;
    }
    if (retval) {
	if (cxpos[0] > cxpos[1])
	    cxpos[0] = cxpos[1];
	if (cypos[0] > cypos[1])
	    cypos[0] = cypos[1];
    }
    return retval;
}


/* Returns 1 if no trim, 0 if had to trim.
   Trims first pair of args to fit within second */

int
regTrimExtent(double *cxpos, double *cypos, double *sxpos, double *sypos,
	      int cstart)
{
/* Here we are intersecting, so the ranges get smaller */
    int retval = 1;
    int verbose = 0;
    if (cstart < 0) {
	cstart = 0;
	verbose = 1;
    }

    if (cstart || sxpos[0] > cxpos[0]) {
	cxpos[0] = sxpos[0];
	retval = 0;
    }
    if (cstart || sxpos[1] < cxpos[1]) {
	cxpos[1] = sxpos[1];
	retval = 0;
    }
    if (cstart || sypos[0] > cypos[0]) {
	cypos[0] = sypos[0];
	retval = 0;
    }
    if (cstart || sypos[1] < cypos[1]) {
	cypos[1] = sypos[1];
	retval = 0;
    }
    if (!retval) {
	if (cxpos[0] > cxpos[1])
	    cxpos[0] = cxpos[1];
	if (cypos[0] > cypos[1])
	    cypos[0] = cypos[1];
    }
    return retval;
}

/*
  +-----------------------------------------------------------------------
 Now returns false if region extends beyond field
  +-----------------------------------------------------------------------
*/
int
regExtent(regRegion * region,
	  double *fieldx, double *fieldy, double *xpos, double *ypos)
{

    regShape *atShape;
    int retval = 1;
    int start;
    int cstart;
    int state = 0;
    double sxpos[2];
    double sypos[2];
    double cxpos[2];
    double cypos[2];

/* The null region is taken to be the field */
    if (!region) {

	xpos[0] = fieldx[0];
	xpos[1] = fieldx[1];
	ypos[0] = fieldy[0];
	ypos[1] = fieldy[1];
	return retval;
    }

    start = regZeroBounds(xpos, ypos);
    cstart = regZeroBounds(cxpos, cypos);

    atShape = region->shape;
    while (atShape != NULL) {
	do {
	    regExtentShape(atShape, fieldx, fieldy, sxpos, sypos);
	    regTrimExtent(cxpos, cypos, sxpos, sypos, cstart);
	    state = 1;
	    cstart = 0;
	    if (atShape->next == NULL) {
		state = 0;
	    } else if (atShape->next->component != atShape->component) {
		state = 0;
	    }

	    atShape = atShape->next;
	}
	while (state);

/* End of component */
	regUnionExtent(xpos, ypos, cxpos, cypos, start);
	start = 0;
	cstart = regZeroBounds(cxpos, cypos);
    }
/* Trim to field */
     retval = regTrimExtent(xpos, ypos, fieldx, fieldy, 0);
    return (retval);
}


void regCornerBounds(double *xpos, double *ypos, double *xb, double *yb)
{
    long i;
    xb[0] = xpos[0];
    xb[1] = xpos[0];
    yb[0] = ypos[0];
    yb[1] = ypos[0];
    for (i = 1; i <= 3; i++) {
	if (xpos[i] < xb[0])
	    xb[0] = xpos[i];
	if (xpos[i] > xb[1])
	    xb[1] = xpos[i];
	if (ypos[i] < yb[0])
	    yb[0] = ypos[i];
	if (ypos[i] > yb[1])
	    yb[1] = ypos[i];
    }
}

/*
 *  Box corners
 */
void regBoxCorners(regShape * shape, double *xpos, double *ypos)
{
    double dx, dy;
    double xcen, ycen;

    if (shape->type == regELLIPSE) {
	dx = shape->radius[0];
	dy = shape->radius[1];
    } else {
	dx = shape->radius[0] / 2.0;
	dy = shape->radius[1] / 2.0;
    }
    xcen = shape->xpos[0];
    ycen = shape->ypos[0];

    regRotatedCoordsInvert(shape, dx, dy, xcen, ycen, &xpos[0], &ypos[0]);
    regRotatedCoordsInvert(shape, -dx, dy, xcen, ycen, &xpos[1], &ypos[1]);
    regRotatedCoordsInvert(shape, dx, -dy, xcen, ycen, &xpos[2], &ypos[2]);
    regRotatedCoordsInvert(shape, -dx, -dy, xcen, ycen, &xpos[3],
			   &ypos[3]);
}



int regRectangleCorners(regShape * shape, double *xpos, double *ypos)
{
    xpos[0] = shape->xpos[0];
    ypos[0] = shape->ypos[0];
    xpos[1] = shape->xpos[1];
    ypos[1] = shape->ypos[1];


    if (shape->angle[0] == 0.0) {
	xpos[2] = shape->xpos[0];
	ypos[2] = shape->ypos[1];
	xpos[3] = shape->xpos[1];
	ypos[3] = shape->ypos[0];
    } else {
	double xr, yr;
	double xcen = (shape->xpos[1] + shape->xpos[0]) / 2.0;
	double ycen = (shape->ypos[1] + shape->ypos[0]) / 2.0;

	regRotatedCoords(shape, shape->xpos[0], shape->ypos[0], xcen, ycen,
			 &xr, &yr);
	regRotatedCoordsInvert(shape, xr, -yr, xcen, ycen, &xpos[2],
			       &ypos[2]);
	regRotatedCoordsInvert(shape, -xr, yr, xcen, ycen, &xpos[3],
			       &ypos[3]);

    }
    return 1;
}

int regRectangleSides(regShape * shape, double *xr, double *yr)
{
    if (shape->angle[0] == 0.0) {
	*xr = shape->xpos[1] - shape->xpos[0];
	*yr = shape->ypos[1] - shape->ypos[0];
    } else {
	double xcen = (shape->xpos[1] + shape->xpos[0]) / 2.0;
	double ycen = (shape->ypos[1] + shape->ypos[0]) / 2.0;

	regRotatedCoords(shape, shape->xpos[0], shape->ypos[0], xcen, ycen,
			 xr, yr);
	*xr *= 2.0;
	*yr *= 2.0;
    }
    return 1;
}
