/*                                                                
**  Copyright (C) 2007  Smithsonian Astrophysical Observatory 
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

#include <float.h>
#include <math.h>
#include "region_priv.h"

#define REGION_MAX_PIXEL_SIZE 0.5

// Utility functions
double sum_areas(long num_shapes, regRegion * region, int *intersecting, int *compno,
	            int *includes, double *areas, double *xpos, double *ypos, double lbin,
	            int *fully_in_field, double *fieldx, double *fieldy, int *union_trim);
void set_flags(long num_shapes, regRegion * region, int *intersecting, 
	                  double *xpos, double *ypos);

double regArea(regRegion * region, double *fieldx, double *fieldy, double lbin)
{
    int ii = 0;           /* counter for variable initialization */
    double area = 0;      /*  area of a single analytic shape */
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
    int within_field = 0;   /* is the region within the field boundaries */
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
    long offset = 0; /* this is an array item number to begin coords from */
    int idx = 0;   /* shape index , shape_no-1 */
    int* union_trim;  /* 1=unioned 2= trimmed */

    
    /* Calculate field area */
    field_area = reg_bounds_area(fieldx, fieldy);


    // A null region is defined to fill the entire field.
    if (!region) {
        return field_area;
    }
    
    atShape = region->shape;
    
    // An empty region has 0 area.
    if (atShape == NULL) {
        return 0.0;
    }
  
    xregbounds = calloc(2, sizeof(double));
    yregbounds = calloc(2, sizeof(double));
  
    /* regExtent returns 1 if all defined shapes fit within field boundaries */
    within_field = regExtent(region, fieldx, fieldy, xregbounds, yregbounds);

    if (atShape->next == NULL && within_field) {

        /* Only one shape; use analytic method for speed? */
        area = reg_shape_analytic_area(atShape, field_area, &status);
        
        if (area < 0.0 && atShape->type != regPOINT)
        {
            area = field_area + area;
        }

        if (status) {
	  free(xregbounds);
	  free(yregbounds);
	  return area;
        }
    }

  
    num_shapes = regGetNoShapes(region);
    areas = calloc(num_shapes, sizeof(double)); /* areas of individual shapes */
    intersecting = calloc(num_shapes, sizeof(int)); /*represents shapes that intersect*/
    compno = calloc(num_shapes, sizeof(int));
    fully_in_field = calloc(num_shapes, sizeof(int));
    includes = calloc(num_shapes, sizeof(int));
    xpos = calloc(2*num_shapes, sizeof(double)); /*shape bound box x-coords*/
    ypos = calloc(2*num_shapes, sizeof(double)); /* shape bound box y-coords */
    x_segments = calloc(14, sizeof(double)); /* x-position info of segments */
    y_segments = calloc(14, sizeof(double)); /* y-position info of segments */
    union_trim = calloc(num_shapes, sizeof(int));

    for (ii=0; ii<num_shapes; ii++) {
        union_trim[ii] = 0;
        intersecting[ii] = 0;
        areas[ii] = 0;
        compno[ii] = 0;
        includes[ii] = 0;
        fully_in_field[ii] = 0;
    }

    /* Pass 1: determine extent and overlap */

    /*  We divide up the field into non-overlapping rectangles
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

        if ( regInclude == shape->include ) {
            includes[idx] = 1;
        }
    
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
        reg_extent_shape_raw(shape, fieldx, fieldy,
		                  xpos + offset, ypos + offset);

        /* trim the shape boundaries to the field boundaries */
        if ( regInclude == shape->include ) {
	        fully_in_field[idx] = reg_trim_extent(xpos + offset, ypos + offset,
                                                fieldx, fieldy, 0);
        } 
    
        /* get the analytic area of each individual shape */
        status = 0;
        areas[idx] = reg_shape_analytic_area(shape, field_area, &status);
        
        if (areas[idx] < 0) {
            areas[idx] = 0;
        }
    } // End 

    /* On each iteration, check to see if any bounding boxes overlap.
       If they do, expand a bounding box to include overlapping shapes.
       Keep doing this until there are no more overlapping bounding boxes,
       then add up the areas of the shapes in each bounding box.
    */
    do {
        set_flags( num_shapes, region, intersecting, xpos, ypos);
        total_area = sum_areas(num_shapes,region,intersecting,compno,includes,
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
    free(compno);
    free(fully_in_field);
    free(includes);
    free(union_trim);

    return total_area;
}

/*
 *  Compute the number of pixels within a complex region bounded by
    xbounds, ybounds 
 */
double regComputePixellatedArea(regRegion * region, double *xbounds,
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

            if (regInsideRegion(region, xpos, ypos)) {
                n += 1;
            }
        }
    }
    
    area = n * bin * bin;
    return area;
}


/* indicate through the intersecting array which shapes truly intersect with each
   other by setting them to the same number (which is the first shape in
   the region string that other shapes intersect with).  ie. if the idx=0,2,
   and 4 shapes intersect, then intersecting[0] = intersecting[2] = intersecting[4] = 1;
*/
void set_flags(long num_shapes, 
              regRegion * region, 
              int *intersecting,
	          double *xpos, 
              double *ypos)
{
    long shape_no = 0;    /* holds the current shape number */
    long prev_no = 0;     /* holds the secondary shape number */
    long offset = 0;      /* this is an array item number to begin coords from */
    long prev_offset = 0; /* pointer to coord values for PREVIOUS shape */
                          /* PREVIOUS can be any shape before the current one */
    
    int idx = 0;          /* shape index , shape_no-1 */
    int pidx = 0;         /* prev. shape index , prev_no-1 */


    /* for each shape in the region
       clear the intersect flags
    */
    for (shape_no = 1; shape_no <= num_shapes; shape_no++)
    {
        idx = shape_no - 1;
        intersecting[idx] = 0;
    }
  
    /* for each shape in the region */
    for (shape_no = 1; shape_no <= num_shapes; shape_no++) {

        idx = shape_no - 1;
      
        /* Does this shape overlap with any other (previous shape)? */
        offset = 2 * (idx);  /* set pointer shape bound box coords */
    
        for (prev_no = 1; prev_no < shape_no; prev_no++) {
	        pidx = prev_no - 1;
	  
	        /* pointer to PREVIOUS shape coords */
	        prev_offset = 2 * (pidx);

	        /* if this shape overlaps with any others */
	        if (reg_rectangle_overlap(xpos + offset, 
                                    ypos + offset,
				                    xpos + prev_offset, 
                                    ypos + prev_offset))
	        {
	        /* if a flag wasn't already set for a PREVIOUS shape,
		       then set it */
	            if (!intersecting[pidx]) {
		            intersecting[pidx] = prev_no;
                }
	      
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
double sum_areas(long num_shapes, 
                regRegion * region, 
                int *intersecting,
                int *compno,
                int *includes, 
                double *areas, 
                double *xpos, 
                double *ypos, 
                double bin,
                int *fully_in_field, 
                double *fieldx, 
                double *fieldy, 
                int *union_trim)
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

        if (compno[idx] > num_components) {
            num_components = compno[idx];
        }
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
    for (cno=1; cno <= num_components; cno++) {

        cidx = cno - 1;  /* define the index */
        coffset = 2 * (cidx); /* pointer to shape bound box coords */

        /* look at each shape on the list.  we'll call each current shape
	       the PRIMARY shape */
        for (shape_no = 1; shape_no <= num_shapes; shape_no++) {
	        
            idx = shape_no - 1;  /* define the index */
            offset = 2 * (idx); /* pointer to shape bound box coords */

	        /* if this shape isn't part of the component,
	           OR if this shape isn't included, continue to the next shape */
	        
            if ( compno[idx] != cno || !includes[idx] ) {
	            continue;
            }

	        /* set the component bounds to be the PRIMARY shape bounds */
	        c_xpos[ coffset ]   = xpos[ offset ];    /* x min */
	        c_ypos[ coffset ]   = ypos[ offset ];    /* y min */
	        c_xpos[ coffset+1 ] = xpos[ offset+1 ];  /* x max */
	        c_ypos[ coffset+1 ] = ypos[ offset+1 ];  /* y max */
	  
	        /* look at all OTHER shapes on the list */
	        for (other_no = 1; other_no <= num_shapes; other_no++) {
	            oidx = other_no - 1;
	            other_offset = 2 * (oidx);

	            /* if the OTHER shape is excluded and it intersects
	               with the included shape, then we must pixellate the
	               included shape/component bounds
	            */
	            if ( !includes[oidx] &&
		             reg_rectangle_overlap(xpos + offset, ypos + offset,
					                     xpos + other_offset,
					                     ypos + other_offset) )
		        {
		            must_pixellate[ cno ] = 1;
		        }
		 
	            /* if this OTHER shape is the same as the PRIMARY shape,
	               OR if this shape isn't part of the component,
	               OR if this shape isn't included
	               continue to the next shape */
	            if ( idx == oidx || compno[oidx] != cno || !includes[ oidx ] ) {
		            continue;
                }

	            /* if this OTHER shape overlaps with the PRIMARY shape */
	            if (reg_rectangle_overlap(xpos + offset, ypos + offset,
				                        xpos + other_offset,
				                        ypos + other_offset))
		        {
		            /* if the bounds are already the same, call it trimmed */
		            if (xpos[offset] == xpos[other_offset] &&
		                ypos[offset] == ypos[other_offset] &&
		                xpos[offset+1] == xpos[other_offset+1] &&
		                xpos[offset+1] == xpos[other_offset+1] )
                    {
		                trimmed=1;
                    }
		            else
                    {
		                /* trim bounds to the intersection */
		                trim_ret =
		                    reg_trim_extent(xpos + offset, ypos + offset,
				                          xpos + other_offset,
				                          ypos + other_offset, 0);
                    }
		            
                    if( 0 == trim_ret ) {
		                trimmed = 1;
                    }

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
        if ( !trimmed && !area_done[ cno ] && !must_pixellate[ cno ]) {
	        
            for (idx = 0; idx < num_shapes; idx++) {
	    
                if ( fully_in_field[ idx ] &&
		             includes[ idx ] && compno[ idx ] == cno )
	            {
		            /* this is a stand-alone shape */
		            comp_areas[ cno ] += areas[ idx ];
		            area_done[ cno ] = 1;
	            }
	        } /* end shape LOOP */
        }
    } /* end for COMPONENT loop */


    /* expand bounding boxes around any overlapping components */
    for (cno=1; cno <= num_components; cno++) {
        idx = cno - 1;  /* define the index */
        offset = 2 * (idx); /* pointer to shape bound box coords */
      
        for (o_cno=1; o_cno <= num_components; o_cno++) {

	        oidx = o_cno - 1;
	        other_offset = 2 * (oidx);
	        
            if (idx == oidx) {
                continue;
            }
            
	        /* if this OTHER component intersects with the PRIMARY component*/
	        if (reg_rectangle_overlap(c_xpos + offset, 
                                    c_ypos + offset,
				                    c_xpos + other_offset,
				                    c_ypos + other_offset))
	        {
		        /* expand bounds to encompass each other */
		        unioned =
		            reg_union_extent(c_xpos + offset, c_ypos + offset,
				                   c_xpos + other_offset,
				                   c_ypos + other_offset, 0);
		        comp_areas[ cno ] = 0;
		        area_done[ cno ] = 0;

		        if ( unioned ) {
		            unioned = 0;  /* reset the unioned flag */
		            cno = 0;  /* every time there's a union, start over */
		            break;
		        }
	        }
	    
	    } /* end for OTHER component loop */
    } /* end for PRIMARY component loop */

    /* pixellate if necessary, per component */
    for (cno=1; cno <= num_components; cno++) {
        dont_compute =0;
        idx = cno - 1;  /* define the index */
        offset = 2 * (idx); /* pointer to shape bound box coords */
      
        for ( o_cno=1; o_cno <= num_components; o_cno++) {

	        oidx = o_cno - 1;
	        other_offset = 2 * (oidx);
	    
            if (idx == oidx) {
                continue;
            }

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

	    if ( area_done[ cno ] || dont_compute ) {
	        continue;
        }
	
        /* if the component bounds have a non-zero area
	       and they overlap with the field bounds
	    */
	    else if ( ((c_xpos[ offset+1 ] - c_xpos[ offset ] != 0) ||
		           (c_ypos[ offset+1 ] - c_ypos[ offset ] != 0) ) &&
		            reg_rectangle_overlap(c_xpos + offset, c_ypos + offset,
				                        fieldx,fieldy) )
	    {
	        comp_areas[ cno ] =
	            regComputePixellatedArea(region, c_xpos + offset,
				       c_ypos + offset, bin);
	        area_done[ cno ] = 1;
	    }
    } /* end for PRIMARY component loop */

    /* add up the component areas */
    for (jj=1; jj <= num_components; jj++) {
        total_area += comp_areas[ jj ];
    }

    if (c_xpos) free(c_xpos);
    if (c_ypos) free(c_ypos);
    if (comp_areas) free(comp_areas);
    if (area_done) free(area_done);
    if (must_pixellate) free(must_pixellate);
    return total_area;
}


double reg_bounds_area(double *fieldx, double *fieldy)
{
    double field_area;
    if (fieldx[1] >= DBL_MAX || fieldx[0] <= -DBL_MAX ||
	    fieldy[1] >= DBL_MAX || fieldy[0] <= -DBL_MAX) 
    {
	    field_area = DBL_MAX;
    }
    else {
	    field_area = (fieldx[1] - fieldx[0]) * (fieldy[1] - fieldy[0]);
    }
    return field_area;
}

double reg_shape_analytic_area(regShape * shape, double field_area, int *status)
{
    double area = 0.0;
    *status = 1;

    area = shape->calcArea(shape);

    if (shape->type == regSECTOR) {
        *status = 0;
    }
    
    if (shape->include == regExclude) {
	    area = -area;
    }
    
    // This is needed here as fields do not have knowledge of the region boundary.
    if (shape->type == regFIELD) {
        area = field_area;
    }
    
    return area;
}
