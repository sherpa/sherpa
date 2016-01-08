/*
 * Includes operations for regPolygon type shapes.
 */
#include <float.h>
#include "region_priv.h"

regShape* regCopyPolygon( regShape * );
int regIsEqualPolygon( regShape *, regShape * );
double regCalcAreaPolygon( regShape * );
int regCalcExtentPolygon( regShape *, double *, double * );
int regInsidePolygon( regShape *, double, double );
void regToStringPolygon( regShape *, char *, long );

int check_overlap(regShape *);
double reg_calc_area_complex_polygon( regShape* shape );
int reg_poly_winding_num(double* xpos, double* ypos, long nPoints, double x, double y);
double reg_poly_is_left(double x1, double y1, double x2, double y2, double x, double y);

regShape* regCreatePolygon(regFlavor include,
                           double *xpos, 
                           double *ypos,
                           long nPoints,
                           int wcoord,
                           int wsize)
{


    if (!xpos || !ypos) {
        fprintf(stderr, "ERROR: Null input for regCreatePolygon\n");
        return (NULL);
    }


    if (nPoints < 3) {
        fprintf(stderr, "ERROR: Polygons must have at least 3 vertices.\n");
        return (NULL);
    }

    regShape *newShape = NULL;
    long ii;
    int addPoint; 

    // create new shape
    newShape = ( regShape *) calloc ( 1, sizeof( regShape ) );

    // Shape type
    newShape->type = regPOLYGON;
    newShape->name = "Polygon";

    // World coords and inclusion
    newShape->include = include;
    newShape->flag_coord = wcoord;
    newShape->flag_radius = wsize;

    // if the first and last value are not the same, add an additional point
    // at the end. We use the value of addpoint to keep track of whether or not
    // we have to add the additional point.
    addPoint = (xpos[0] != xpos[nPoints-1] || ypos[0] != ypos[nPoints-1]) ? 1 : 0;
    nPoints += addPoint;

    // Allocate memory for points
    newShape->xpos = (double *) calloc(nPoints, sizeof(double));
    newShape->ypos = (double *) calloc(nPoints, sizeof(double));
    newShape->nPoints = nPoints;
    
    // Fill in values (ignoring last point if necessary)
    for (ii=0; ii<nPoints - addPoint; ii++) {
        newShape->xpos[ii] = xpos[ii];
        newShape->ypos[ii] = ypos[ii];
    }
    
    // Add in extra point if we're supposed to.
    if (addPoint) {
        newShape->xpos[nPoints-1] = xpos[0];
        newShape->ypos[nPoints-1] = ypos[0];
    }

    newShape->angle = NULL;
    newShape->radius = NULL;
    
    // Add relevant methods
    newShape->calcArea = regCalcAreaPolygon;
    newShape->calcExtent = regCalcExtentPolygon;
    newShape->copy = regCopyPolygon;
    newShape->isEqual = regIsEqualPolygon;
    newShape->isInside = regInsidePolygon;
    newShape->toString = regToStringPolygon;

    // Verify that the polygon's lines don't just go out and back (i.e. that there is
    // no line segment causing 0 width on the polygon. 
    for (ii=0;ii<nPoints-2;ii++) {
	    if ((newShape->xpos[ii] == newShape->xpos[ii+2]) &&
	        (newShape->ypos[ii] == newShape->ypos[ii+2]) &&
	        (ii+2 != (nPoints-1))) 
        {
	        fprintf(stderr, "WARNING: Polgyon must have finite width; adjacent line segments with"
		            " ends at (%g,%g) overlap completely (index = %lu)\n",
		            newShape->xpos[ii],newShape->ypos[ii], ii );
	        // TODO:
            //  Disallow malformed polygons?
            // return(NULL);
	    }
    }
    
    // Verify that the polygon does not have duplicate points in a row.
    for (ii=0;ii<nPoints-2;ii++) {
	    if ((newShape->xpos[ii] == newShape->xpos[ii+1]) &&
	        (newShape->ypos[ii] == newShape->ypos[ii+1]) ) 
        {
	        fprintf(stderr, "WARNING: Zero length polygon line segment at (%g,%g) (index = %lu).\n",
		            newShape->xpos[ii], newShape->ypos[ii], ii);
	        // TODO:
            //  Disallow malformed polygons?
            // return (NULL);
	    }
    }
    
    return newShape;
} // end regCreatePolygon

// copy
regShape* regCopyPolygon( regShape* shape ) {
    if (shape->type != regPOLYGON) {
	    fprintf( stderr, "ERROR: Attempting to copy %s as a Polygon\n", shape->name);
	    return(NULL);
    }

    return regCreatePolygon(shape->include,
                            shape->xpos,
                            shape->ypos,
                            shape->nPoints,
                            shape->flag_coord,
                            shape->flag_radius);
}

// equals
int regIsEqualPolygon( regShape* thisShape, regShape* otherShape ) { 
    if ( !thisShape && !otherShape ) {
        return 1;
    }

    if ( !thisShape || !otherShape) {
        return 0;
    }
    
    long ii;
    
    if (thisShape->type != regPOLYGON) {
	    fprintf( stderr, "ERROR: not comparing a Polygon\n");
    }

    if (otherShape->type != regPOLYGON) {
        return 0;
    }

    if (thisShape->include != otherShape->include ||
        thisShape->nPoints != otherShape->nPoints)
    {
        return 0;
    }
    
    for (ii=0; ii<thisShape->nPoints; ii++) {
        if (thisShape->xpos[ii] != otherShape->xpos[ii] ||
            thisShape->ypos[ii] != otherShape->ypos[ii])
        {
            return 0;
        }
    }
    return 1;
}

// calcArea
//   As is the case with most theorems for polygons, this formula
//   does not work for complex polygons (i.e., one with overlapping 
//   lines). Currently this is expected behavior. Let's fix this - 
//   either by disallowing complex polygons, forcing a brute force 
//   area computation, or by adding in intersect points to make the
//   calculation work in all cases.
double regCalcAreaPolygon( regShape* shape ) {
    
    // We brute force the area computation for complex polygons
    if (check_overlap(shape)) {
        return reg_calc_area_complex_polygon(shape);   
    }

    double area = 0.0;
    long ii;

    // last point is already there
    // eqn is 0.5*sum( x_i*y+i_1 - x_i+1*y_i)
    for (ii=0;ii<(shape->nPoints -1);ii++) {
        area += ((shape->xpos[ii]*shape->ypos[ii+1]) -
	             (shape->xpos[ii+1]*shape->ypos[ii]) );
    }
    area /= 2.0;

    return fabs(area);
}

double reg_calc_area_complex_polygon( regShape* shape ) {
    
    regRegion *temp;
    regShape *copy;
    double area;

    fprintf(stderr, "WARNING: Calculating area of a complex polygon ");
    fprintf(stderr,"using brute force method. This may take a long time.\n");

    // Create a new region with just the polygon
    temp = regCreateRegion(NULL, NULL);
    copy = shape->copy(shape);

    // Analytic area calculations always computes the area of the interior
    // of the shape.    
    copy->include = regInclude;
    regAddShape(temp, regAND, copy);
    
    // Calc the extent of the polygon then trim the bounds to fit within
    // the original region if available
    regCalcExtentPolygon(shape, temp->xregbounds, temp->yregbounds);
    if (shape->region) {
        reg_trim_extent(temp->xregbounds, 
                        temp->yregbounds,
                        shape->region->xregbounds,
                        shape->region->yregbounds, 0);
    }

    area = regComputePixellatedArea(temp, temp->xregbounds, temp->yregbounds, 1);

    // Free and return
    regFree(temp);
    return area;
}

// calcExtent
int regCalcExtentPolygon( regShape* shape, double* xpos, double* ypos ) {
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

/*
 * Determines if a point (x,y) is inside the polygon using the winding number
 * method.
 *
 *  http://www.engr.colostate.edu/~dga/dga/papers/point_in_polygon.pdf 
 */
int regInsidePolygon( regShape* shape, double xpos, double ypos ) {
    
    int wn = reg_poly_winding_num(shape->xpos, shape->ypos, shape->nPoints, xpos, ypos);
    int ret; 
    
    // wn == 0 iff (x,y) is outside of the polygon.
    ret = wn != 0;
    
    if (shape->include == regExclude) {
        return !ret;
    }
    
    return ret;
}

/*
 * Computes the winding number for the point (x,y) within the polygon described by
 * xpos and ypos.
 */
int reg_poly_winding_num(double* xpos, double* ypos, long nPoints, double x, double y) 
{
    int wn = 0;
    long ii;
    double isLeft;
    
    for (ii=0; ii<nPoints - 1; ii++) {

        isLeft = reg_poly_is_left(xpos[ii], ypos[ii], xpos[ii+1], ypos[ii+1], x, y);
        
        // This needs to be here since we include points on the edges of the polygon as
        // being "inside". The winding number algorithm has some trouble with edge cases,
        // so to make it simple we just include a check to see if the point is on one
        // of the edges.
        if (isLeft == 0) {
            // The points are colinear, need to check if (x,y) is on the edge
            // segment.
            if (((x <= xpos[ii] && x >= xpos[ii+1]) ||
                 (x <= xpos[ii+1] && x >= xpos[ii])) 
                &&
                ((y <= ypos[ii] && y >= ypos[ii+1]) ||
                 (y <= ypos[ii+1] && y >= ypos[ii])))
            {
                return 1;
            }
        }

        // O/w check if the polygon edge has winded around the point.
        if (ypos[ii] <= y) {
            if (ypos[ii+1] > y) {
                // (x,y) is left of edge
                if (isLeft > 0) {
                        wn++;
                }
            }
        }
        else {
            if (ypos[ii+1] <= y) {
                // (x,y) is right of the edge
                if (isLeft < 0) {
                       wn--;
                }
            }
        }
    }

    return wn;
}

/*
 * Determines if the point (x, y) is to the left, right, or on an infinite
 * line defined by (x1,y1), (x2,y2) by evaluating the slope of the line segments
 * between the three points.
 *
 * Returns:
 *  > 0 if point is on the left
 *  = 0 if point is on line
 *  < 0 otherwise
 */
double reg_poly_is_left(double x1, double y1, double x2, double y2, double x, double y) {
    double ret = ((x2 - x1) * (y - y1) - (x - x1) * (y2 - y1));
    return ret;
}

// toString
void regToStringPolygon( regShape* shape, char* ptr, long maxlength ) {

    if (!shape) return;
  
    // Using these to keep track of the length of the word so far, and 
    // the lengths of text that gets added for each point. 
    long length, word;
    long ii;

    length = 0;
    word = 0;

    if (shape->include == regExclude) {
        *ptr++ = '!';
        length++;
    }
    
    long bsize = 80;

    word = snprintf(ptr, maxlength - length, "Polygon(");
    length += word;
    ptr += word;
    
    for (ii=0; ii < shape->nPoints; ii++) {
        char *x = calloc(bsize, sizeof(char));
        char *y = calloc(bsize, sizeof(char));
        reg_print_pos_pair(shape->xpos[ii], shape->ypos[ii], 
                shape->flag_coord, x, y);

        if (ii == 0) {
            word = snprintf(ptr, maxlength - length, "%s,%s", x, y);
        }
        else {
             word = snprintf(ptr, maxlength - length, ",%s,%s", x, y);
        }
        length += word;
        ptr += word;

        free(x);
        free(y);
    }

    word = snprintf(ptr, maxlength - length, ")");
    length += word;
    ptr += word;
}

// This does a check to see if any of the lines of the 
// polgyons over lap (thus making a complex polygon). 
// Make a warning for now
//
int check_overlap(regShape * newShape) {
 
    long nPoints = newShape->nPoints;
    long ii, jj;

    for (ii=0;ii<nPoints-3;ii++) {
	    double x1, x2, y1, y2;

	    x1=newShape->xpos[ii];
	    x2=newShape->xpos[ii+1];
	    y1=newShape->ypos[ii];
	    y2=newShape->ypos[ii+1];

	    for (jj=ii+2;jj<=nPoints-2;jj++) {
	        double x3,x4, y3,y4;
	  
	        // last point is == 1st by definition
	        if ( ( ii == 0 ) && ( jj == nPoints-2 )) continue;
	  
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
	        } 
            else {
	            ua=(x4-x3)*(y1-y3)-(y4-y3)*(x1-x3);
	            ub=(x2-x1)*(y1-y3)-(y2-y1)*(x1-x3);
	            ua/= denom;
	            ub/= denom;

	            if ( (0 <= ua) && ( ua <= 1 ) &&
		             (0 <= ub) && ( ub <= 1 )) 
                {
                    return 1;
	                /* ---- good for debugging
		            fprintf( stderr, "  linesegment #1:  (xlo,ylo) = (%g,%g)\t (xhi,yhi) = (%g,%g)\n",x1,y1,x2,y2);
		            fprintf( stderr, "  linesegment #2:  (xlo,ylo) = (%g,%g)\t (xhi,yhi) = (%g,%g)\n",x3,y3,x4,y4);
		            fprintf( stderr, "  ii=%ld \t jj=%ld\n", ii, jj);
	                */
                
	                // Go ahead and break out; no need to test the remaining data points.
	                // ii=nPoints;
	                // jj=nPoints;
	            }
            } // end else
        } // end jj
    }// end ii
    return 0;
}
