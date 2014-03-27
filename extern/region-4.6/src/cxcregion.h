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

#ifndef REGION_H
#define REGION_H

/* Public interface to region library */

typedef struct regSHAPE regShape;
typedef struct regREGION regRegion;
typedef void (*regInvertFunction)(double*, double*);


/* Public prototypes */
regRegion* regParse( char* buf );
void regFree( regRegion* region );

extern void regPrintRegion( regRegion* region );
extern void regPrintShape( regShape *shape );
extern int regInsideRegion( regRegion *region, double xpos, double ypos );
extern regRegion* regCombineRegion( regRegion* Region1, regRegion* Region2 );

regRegion* regCopyRegion( regRegion* inRegion );

double regComputePixellatedArea(regRegion * region, double *xbounds,
				double *ybounds, double bin);
extern int regRegionToList( regRegion *region, double xmin, double xmax,
		     double ymin, double ymax, double bin, double **xat,
		     double **yat, long  *nat );

extern int regRegionToMask( regRegion *region, double xmin, double xmax, double ymin,
		     double ymax, double bin, short **mask, long   *xlen,
		     long *ylen );
		     
void regConvertWorldShape( regShape* Shape, double scale, regInvertFunction invert, int force );
void regConvertWorldRegion( regRegion* region, double scale, regInvertFunction invert );
void regConvertRegion(
    regRegion* Region,         /* (i) Region which may have shapes in world units */
    double scale,              /* (i) degrees per pixel */
    regInvertFunction invert,      /* (i) Function to convert world to pixel coords */
    int force
 );

void regComposeRegion( const regRegion* region, char* buf, const long maxlen );
int regCompareShape( regShape* Shape1, regShape* Shape2 );
int regCompareShapeRaw( regShape* Shape1, regShape* Shape2 );
int regCompareRegion( regRegion* Region1, regRegion* Region2 );

extern regRegion* regCreateRegion( void*, void* );
extern regRegion* regCreateEmptyRegion( void );
extern void regAppendShape( regRegion* region,
  char* shapeName, int includeFlag, int orFlag, double* xpos, double* ypos,
  long npoints, double* radius, double* angle, int world_coord, int world_size );

/*
 * Routines to get individual shapes in a region
 */
regShape* regGetShapeNo( const regRegion* region, long shapeNo );
long regGetNoShapes( const regRegion* region );

long regShapeRadii( const regShape* Shape );
long regShapeAngles( const regShape* Shape );
extern regShape* regCopyShape( regShape *Shape);

regRegion* regShapeGetRegion( const regShape* shape );
int regShapeGetName( const regShape* shape, char* name, long maxlen );
long regShapeGetPoints( const regShape* shape, double* x, double* y, long dim );
long regShapeGetAngles( const regShape* shape, double* angle );
long regShapeGetRadii( const regShape* shape, double* r );
long regShapeGetNoPoints( const regShape* shape );
long regShapeGetComponent( const regShape* shape );
long regGetMaxPoints( const regRegion* region );
int regInsideShape( regShape *shape, double x, double y);
char* regAllocComposeRegion( const regRegion* region );

double regArea( regRegion* region, double* fieldx, double* fieldy, double bin );
int regExtent( regRegion *region, double* fieldx, double* fieldy, double* xpos,   double* ypos  );
regRegion* regReadAsciiRegion( char* filename, int verbose );

regRegion* regUnionRegion( regRegion* Region1, regRegion* Region2 );
regRegion* regIntersectRegion( regRegion* Region1, regRegion* Region2 );
int regComponentOverlap( regShape* Shape1, regShape* Shape2 );
int regWriteAsciiRegion( char* name, regRegion* region, char** names, long nn );
regRegion* regInvert( regRegion* inRegion );
regRegion* regCreateRectangle( double* x, double* y );
void regResolveField( regRegion* region, double* x, double* y );
#endif
