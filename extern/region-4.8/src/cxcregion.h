/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#ifndef REGION_H
#define REGION_H

/* Public interface to region library */

typedef struct   regSHAPE regShape;
typedef struct   regREGION regRegion;
typedef void     (*regInvertFunction)(double*, double*);


/* Region Create/Free */
regRegion*   regCreateRegion( void*, void* );
regRegion*   regCreateEmptyRegion( void );
regRegion*   regCopyRegion( regRegion* inRegion );
regRegion*   regParse( char* buf );
void         regFree( regRegion* region );


/* Shape Accessors */
int          regShapeGetName( const regShape* shape, char* name, long maxlen );
long         regShapeGetPoints( const regShape* shape, double* x, double* y, long dim );
long         regShapeGetAngles( const regShape* shape, double* angle );
long         regShapeGetRadii( const regShape* shape, double* r );
long         regShapeGetNoPoints( const regShape* shape );
long         regShapeGetComponent( const regShape* shape );


/* Shape Operations */
void        regAppendShape( regRegion* region, char* shapeName, int includeFlag, int orFlag,
			    double* xpos, double* ypos, long npoints, double* radius, 
			    double* angle, int world_coord, int world_size );
int         regCompareShape( regShape* Shape1, regShape* Shape2, short raw );
regShape*   regCopyShape( regShape *Shape);
int         regInsideShape( regShape *shape, double x, double y);
void        regPrintShape( regShape *shape );


/* Region Computation */
double      regArea( regRegion* region, double* fieldx, double* fieldy, double bin );
double      regComputePixellatedArea(regRegion * region, double *xbounds,
				     double *ybounds, double bin);
int         regExtent( regRegion *region, double* fieldx, double* fieldy, 
		       double* xpos, double* ypos  );


/* Region Examination */
int         regCompareRegion( regRegion* Region1, regRegion* Region2 );
int         regInsideRegion( regRegion *region, double xpos, double ypos );
long        regGetMaxPoints( const regRegion* region );
regShape*   regGetShapeNo( const regRegion* region, long shapeNo );
long        regGetNoShapes( const regRegion* region );


/* Region Conversion */
void        regConvertWorldRegion( regRegion* region, double scale, regInvertFunction invert );
void        regConvertRegion(regRegion* Region, double scale, regInvertFunction invert, int force );
regRegion*  regInvert( regRegion* inRegion );
int         regRegionToList( regRegion *region, double xmin, double xmax,
			     double ymin, double ymax, double bin, double **xat,
			     double **yat, long  *nat );
int         regRegionToMask( regRegion *region, double xmin, double xmax, double ymin,
			     double ymax, double bin, short **mask, long   *xlen,
			     long *ylen );
void        regResolveField( regRegion* region, double* x, double* y );
		     

/* Region Logical Combination*/ 
regRegion*  regCombineRegion( regRegion* Region1, regRegion* Region2 );
regRegion*  regUnionRegion( regRegion* Region1, regRegion* Region2 );
regRegion*  regIntersectRegion( regRegion* Region1, regRegion* Region2 );


/* Region Ascii Parsing */
regRegion*  regReadAsciiRegion( char* filename, int verbose );
int         regWriteAsciiRegion( char* name, regRegion* region, char** names, long nn );


/* Region printing */
void        regPrintRegion( regRegion* region );
void        regComposeRegion( const regRegion* region, char* buf, const long maxlen );
char*       regToStringRegion(const regRegion* region);
char*       regAllocComposeRegion(const regRegion* region);


#endif
