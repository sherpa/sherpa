#include "region_priv.h"
#include <float.h>
#include <ctype.h>


int reg_trim_double(double *, double, double);


void reg_corner_bounds(double *xpos, double *ypos, double *xb, double *yb)
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

void reg_box_corners(regShape * shape, double *xpos, double *ypos)
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

    reg_rotated_coords_invert(shape,  dx,  dy, xcen, ycen, &xpos[0], &ypos[0]);
    reg_rotated_coords_invert(shape, -dx,  dy, xcen, ycen, &xpos[1], &ypos[1]);
    reg_rotated_coords_invert(shape,  dx, -dy, xcen, ycen, &xpos[2], &ypos[2]);
    reg_rotated_coords_invert(shape, -dx, -dy, xcen, ycen, &xpos[3], &ypos[3]);

}


int reg_rectangle_corners(regShape * shape, double *xpos, double *ypos)
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

	reg_rotated_coords(shape, shape->xpos[0], shape->ypos[0], xcen, ycen,
			 &xr, &yr);
	reg_rotated_coords_invert(shape, xr, -yr, xcen, ycen, &xpos[2],
			       &ypos[2]);
	reg_rotated_coords_invert(shape, -xr, yr, xcen, ycen, &xpos[3],
			       &ypos[3]);

    }
    return 1;
}

int reg_rectangle_sides(regShape * shape, double *xr, double *yr)
{
    if (shape->angle[0] == 0.0) {
	*xr = shape->xpos[1] - shape->xpos[0];
	*yr = shape->ypos[1] - shape->ypos[0];
    } else {
	double xcen = (shape->xpos[1] + shape->xpos[0]) / 2.0;
	double ycen = (shape->ypos[1] + shape->ypos[0]) / 2.0;

	reg_rotated_coords(shape, shape->xpos[0], shape->ypos[0], xcen, ycen,
			 xr, yr);
	*xr *= 2.0;
	*yr *= 2.0;
    }
    return 1;
}


double reg_mod_angle( double ang )
{
 double angle;
 angle = fmod( ang, 360.0 );
 if ( angle < 0.0 ) angle += 360.0;
 return angle;
}



void reg_rotated_coords_invert( regShape* shape, double xr, double yr, double xcen, double ycen, double* xp, double* yp )
{
 double ct, st;
 double rotpos[2];
  if ( shape->angle[0] != 0.0 )
    {
      /*double theta = shape->angle[0] * PI / 180.0;*/
      ct = *shape->cos_theta;
      st = *shape->sin_theta;
    } else {
      ct = 1.0;
      st = 0.0;
    }
    rotpos[0] =  ct * xr - st * yr;
    rotpos[1] =  st * xr + ct * yr;
    *xp = rotpos[0] + xcen;
    *yp = rotpos[1] + ycen;

}

void reg_rotated_coords( regShape* shape, double xp, double yp, double xcen, double ycen, double* xr, double* yr )
{
 double ct, st;

  if ( shape->angle[0] != 0.0 )
    {
      /*double theta = -shape->angle[0] * PI / 180.0;*/
      ct = -*shape->cos_theta;
      st = -*shape->sin_theta;
    } else {
      ct = 1.0;
      st = 0.0;
    }
    *xr =  ct * ( xp - xcen ) - st * ( yp - ycen );
    *yr=  st *  ( xp - xcen ) + ct * ( yp - ycen );
}



int reg_shape_overlap( regShape* shape1, regShape* shape2 )
{
    double fx[2] = { -DBL_MAX, DBL_MAX };
    double fy[2] = { -DBL_MAX, DBL_MAX };
    double xpos1[2];
    double ypos1[2];
    double xpos2[2];
    double ypos2[2];
    int ok = 1;

    if ( !shape1 || !shape2 ) {
        return 0;
    }

    /* Exclusions always overlap */
    if ( shape1->include != regInclude || shape2->include != regInclude ) { 
        return 1;
    }

    reg_extent_shape( shape1, fx, fy, xpos1, ypos1 );
    reg_extent_shape( shape2, fx, fy, xpos2, ypos2 );
    ok = reg_rectangle_overlap( xpos1, ypos1, xpos2, ypos2 );
    
    return ok;
}

int reg_rectangle_overlap( double* xpos1, double* ypos1, double* xpos2, double* ypos2 )
{
    // if one is on one side or below the other return false
    if (   xpos1[1] < xpos2[0] 
        || xpos1[0] > xpos2[1] 
        || ypos1[1] < ypos2[0] 
        || ypos1[0] > ypos2[1] )
    {
        return 0;
    }
    return 1;
}


long reg_shape_radii( const regShape* Shape )
{
    long nradii;
    switch( Shape->type ) {
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
        default:
            nradii = 0;
            break;
    }
    return nradii;
}


long reg_shape_angles( const regShape* Shape )
{
    long nangles;
    switch( Shape->type )
    {
        case regELLIPSE:
        case regBOX:
        case regRECTANGLE:
            nangles = 1;
            break;
        case regPIE:
        case regSECTOR:
            nangles = 2;
        break;
        case regFIELD:
        case regCIRCLE:
        case regANNULUS:
        case regPOINT:
        case regPOLYGON:
        case regMASK:
        default:
            nangles = 0;
        break;
    }
    
    return nangles;
}


int reg_case_equal(char *s1, char *s2)
{
    while (1) {
	if (toupper(*s1) != toupper(*s2))
	    return 0;
	if (*s1 == '\0' || *s2 == '\0')
	    return 1;
	s1++;
	s2++;
    }
}


regGeometry reg_shape_name_to_geometry(char *name)
{
    long ntypes = 12;
    char *names[] = { "Circle", "Point", "Ellipse", "Pie", "Sector",
		      "Annulus", "Polygon", "Box", "Rectangle", 
		      "RotBox", "RotRectangle", "Field" };
    
    regGeometry geoms[] = { regCIRCLE, regPOINT, regELLIPSE, regPIE, 
			    regSECTOR, regANNULUS, regPOLYGON, regBOX, 
			    regRECTANGLE, regBOX, regRECTANGLE, regFIELD };
    long ii;

    for (ii = 0; ii < ntypes; ii++) {
	if (reg_case_equal(name, names[ii]))
	  return geoms[ii];
    }

    /* Unknown case */
    return regPOINT;
}

/*
 *  Given an array of positions of size nmax, find the number
 *  of positions actually used. For polygons, this is found because
 *  the last and first points must be identical.
 */
long reg_shape_find_npoints(regGeometry type, double *xpos, double *ypos,
			 long nmax)
{

    long n = 0;
    double x, y;
    int loop;
    switch (type) {
    case regPOLYGON:
	x = xpos[0];
	y = ypos[0];
	n = 1;
	loop = 1;
	while ( n < nmax && loop) {
	  if (( x == xpos[n]) && (y == ypos[n]))
	    loop = 0;
	  else
	    n++;
	}
	if (loop)
	    n = nmax;
	break;
    case regPOINT:
    case regCIRCLE:
    case regANNULUS:
    case regELLIPSE:
    case regBOX:
    case regPIE:
    case regSECTOR:
	n = 1;
	break;
    case regRECTANGLE:
	n = 2;
	break;
    case regFIELD:
    case regMASK:
    default:
	n = 0;
	break;
    }
    return n;
}


int reg_trim_extent(double *cxpos, double *cypos, double *sxpos, double *sypos,
        int cstart)
{
    int outside_field;
    
    if (cstart) {
        cxpos[0] = sxpos[0];
        cxpos[1] = sxpos[1];
        cypos[0] = sypos[0];
        cypos[1] = sypos[1];
        return 0;
    }

    outside_field = reg_trim_double(cxpos, sxpos[0], sxpos[1]);
    outside_field = reg_trim_double(cxpos+1, sxpos[0], sxpos[1]) || outside_field;
    outside_field = reg_trim_double(cypos, sypos[0], sypos[1]) || outside_field;
    outside_field = reg_trim_double(cypos+1, sypos[0], sypos[1]) || outside_field;
    
    return !outside_field;
}

/*
 * Checks if the value at *x is contained by x_1 and x_2, if not it assigns it to 
 * the closer number. Returns true if x is not contained in (x_1, x_2).
 */
int reg_trim_double(double *x, double x_1, double x_2) {
    if (*x < x_1 && *x < x_2) {
        *x = x_1;
        return 1;
    }
    if (*x > x_1 && *x > x_2) {
        *x = x_2;
        return 1;
    }
    return 0;
}

/*
 *  Reverse the polarity
 */
void regNegate( regShape* Shape )
{
  Shape->include = (Shape->include) ? regExclude : regInclude; 
}
