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

#include "cxcregion.h"

/* Private interface to region library */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifndef SZ_CARD
#define SZ_CARD 80
#endif

#ifndef SZ_LARGE
#define SZ_LARGE 32767
#endif


typedef enum
{
  regAND,
  regOR
} regMath;

typedef enum
{
  regExclude,
  regInclude
} regFlavor;

typedef enum
{
  regPOINT,
  regBOX,
  regROTBOX = regBOX,
  regRECTANGLE,
  regROTRECTANGLE = regRECTANGLE,
  regCIRCLE,
  regELLIPSE,
  regPIE,
  regSECTOR,
  regPOLYGON,
  regANNULUS,
  regFIELD,
  regMASK
} regGeometry;


typedef struct regSHAPE
{

  regGeometry  type;
  char        *name;
  regFlavor    include;
  double      *xpos;
  double      *ypos;
  long         nPoints;
  double      *radius;
  double      *angle;
  double      *sin_theta;
  double      *cos_theta;
  long         component;
  void*        spec;        /* Object containing shape specifications */
                            /*    for Pixel Mask == dmBlock           */

  /* Coordinate flags take the values RC_PHYSICAL, RC_LOGICAL, RC_WORLD, RC_UNK  */
  int          flag_coord;  /* Coordinate system for x,y positions (center of circle, etc.) */
  int          flag_radius;   /* Coordinate system for radial sizes (circle radius, etc.) */

  /* Operations */
  /* void is the pointer of regShape */
#define shapeP struct regSHAPE *

  double    (*calcArea)   ( shapeP shape );
  int       (*calcExtent) ( shapeP shape, double* xpos , double* ypos );
  shapeP    (*copy)       ( shapeP shape );
  int       (*isEqual)    ( shapeP thisShape, shapeP otherShape );
  int       (*isInside)   ( shapeP shape, double x, double y );
  void      (*toString)   ( shapeP shape, char * str, long maxlength );

  struct regREGION *region;
  struct regSHAPE  *next;
} regSHAPE;


typedef struct regREGION
{
  struct regSHAPE *shape;  /* Linked-list of Shapes */
  double xregbounds[2];
  double yregbounds[2];
} regREGION;


#define RF_UNK      0
#define RF_SAOIMAGE 1
#define RF_SAOTNG   2
#define RF_PROS     3
#define RF_CXC      4
#define RF_DS9      5
#define RF_DS9_V4   6

#define RC_UNK      0
#define RC_LOGICAL  1
#define RC_PHYSICAL 2
#define RC_WORLD    3


/* Shape Creation Functions  */

regShape* regCreateShape(regRegion *region,
                         regMath glue,
                         regGeometry shape,
                         regFlavor include,
                         double *xpos,
                         double *ypos,
                         long   npoints,
                         double *radius,
                         double *angle,
			 int wcoord,
			 int wsize); 


regShape* regCreatePoint(regFlavor include, 
			 double *xpos, 
			 double *ypos,
			 int wcoord, 
             int wsize);
    
regShape* regCreateEllipse(regFlavor include, 
			   double *xpos,
			   double *ypos,
			   double *radius, 
			   double *angle, 
			   int wcoord, 
               int wsize);

regShape* regCreateCircle(regFlavor include, 
			  double *xpos, 
			  double *ypos,
			  double *radius, 
			  int wcoord, 
              int wsize);

regShape* regCreateBox(regFlavor include, 
		       double *xpos, 
		       double *ypos,
		       double *radius, 
		       double *angle, 
		       int wcoord, 
               int wsize);

regShape* regCreateRectangle(regFlavor include, 
			     double *xpos, 
			     double *ypos, 
			     double *angle, 
			     int wcoord, 
                 int wsize);

regShape* regCreatePie(regFlavor include, 
		       double *xpos, 
		       double *ypos, 
		       double *radius,
		       double *angle, 
		       int wcoord, 
               int wsize);

regShape* regCreateAnnulus(regFlavor include,
                           double *xpos, 
                           double *ypos,
                           double *radius, 
			               int wcoord, 
                           int wsize);

regShape* regCreatePolygon(regFlavor include,
                           double *xpos, 
                           double *ypos,
                           long nPoints, 
			               int wcoord, 
                           int wsize);

regShape* regCreateSector(regFlavor include,
                          double *xpos, 
                          double *ypos,
                          double *angle, 
			              int wcoord, 
                          int wsize);

regShape* regCreateField(regFlavor include, 
			 int wcoord, int wsize);

regShape* regCreateMask(regFlavor include, 
			int wcoord, int wsize);


/* Utility Functions */


// Janine's additions
regShape*   regCreateNewWorldShape( regGeometry type,
				    regFlavor include,
				    double *xpos,
				    double *ypos,
				    long   npoints,
				    double *radius,
				    double *angle,
				    int flag_coord,    /* xpos ypos in degrees or pixels? */
				    int world_size      /* radius in degrees or pixels? */
				    );
void        regNegate(regShape* Shape );


long        regAddShape( regRegion *region,
			 regMath glue,
			 regShape *shape );
  
void        regFreeShape( regRegion* region, regShape* atShape );



regGeometry reg_shape_name_to_geometry(char *name);

void        reg_corner_bounds( double* xpos, double* ypos, double* xb, double* yb );
void        reg_box_corners( regShape* shape, double* xpos, double* ypos );
int         reg_rectangle_corners( regShape *shape, double* xpos, double* ypos );
int         reg_rectangle_sides(regShape * shape, double *xr, double *yr);
void        reg_rotated_coords( regShape* shape, double xp, double yp, double xcen,
			      double ycen, double* xr, double* yr );
void        reg_rotated_coords_invert( regShape* shape, double xr, double yr, double xcen, 
				    double ycen, double* xp, double* yp );
double      reg_mod_angle( double ang );
long        reg_shape_find_npoints(regGeometry type, double *xpos, double *ypos,
				   long nmax);
int         reg_case_equal(char *s1, char *s2);

/* 
 * Erik's additions
 */
int         reg_union_extent(double *, double *, double *, double *, int cstart);
int         reg_extent_shape(regShape * shape, double *fieldx, double *fieldy, 
			     double *xpos, double *ypos);
int         reg_extent_shape_raw(regShape *, double *, double *, double *, double *);
int         reg_trim_extent(double *, double *, double *, double *, int );
int         reg_shape_intersect(regShape *, regShape *, long *, long *);
int         reg_intersect_component(regRegion *, regShape *, regShape *);
double      reg_shape_analytic_area(regShape *, double, int *);
double      reg_bounds_area(double *fieldx, double *fieldy);
int         reg_shape_overlap( regShape *, regShape * );

long        reg_shape_radii( const regShape* Shape );
long        reg_shape_angles( const regShape* Shape );
int         reg_rectangle_overlap(double *, double *, double *, double *);

void        reg_convert_world_shape( regShape* Shape, double scale, 
				     regInvertFunction invert, int force );

regShape*   reg_next_component( regShape * );

void reg_print_pos_pair(double x, double y, int world, char *xbuf, char *ybuf);
void reg_print_pos(double x, int world, char *buf);
void reg_print_radius(double r, int world, char *buf);
char *reg_print_val(double x, char *buf);

