/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
/*
 *  Private to region library
 */

/* prototypes */

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
  regMASK,
  regUSER
} regGeometry;


struct regUSERSHAPE
{
  char*       name;        /* Simple name of type (e.g. "Circle" )   */
  void*       spec;        /* Object containing shape specifications */
                           /*  for Mask == dmBlock                   */
  int         flag_coord;  /* (logical, physical, world, l+p,p+w,l+p+w) = {regCoords}  */

  /* Operations */
  /* void is the pointer of regShape */
#define shapeP struct regSHAPE *

  shapeP    (*copy)   ( shapeP shape);
  int       (*isEqual)( shapeP thisShape, shapeP otherShape );
  double    (*calcArea)(shapeP shape);
  int       (*calcExtent)( shapeP shape, double* xpos , double* ypos );
  int       (*isInside)( shapeP shape, double x, double y );
  char*     (*toString)( shapeP shape );
  void      (*free)    ( shapeP shape );
};

struct regSHAPE
{
  regGeometry type;
  regFlavor include;
  double *xpos;
  double *ypos;
  long   nPoints;
  double *radius;
  double *angle;
  double *sin_theta;
  double *cos_theta;
  long   component;
  int world_coord;
  int world_size;
  int   (*inside)( struct regSHAPE *, double , double);
  struct regREGION *region;

  struct regSHAPE *next;
  struct regUSERSHAPE *user;
};



typedef enum
{
  regCHAR,
  regSHORT,
  regUSHORT,
  regLONG,
  regULONG,
  regFLOAT,
  regDOUBLE
} regDataType;



struct regREGION
{
  struct regSHAPE *shape;
  void *xcol[2];  /* column structures, edsColHandles */
  double xregbounds[2];
  double yregbounds[2];
};

#define RF_UNK      0
#define RF_SAOIMAGE 1
#define RF_SAOTNG   2
#define RF_PROS     3
#define RF_CXC      4
#define RF_DS9      5
#define RF_DS9_V4   6

#define RC_PHYSICAL 1
#define RC_LOGICAL  2
#define RC_WORLD    3
#define RC_UNK 0

#include "cxcregion.h"

/* Internal prototypes */

/*
 *  Shapes
 */

extern long regAddShape( regRegion *region, regMath glue, regShape *shape );

extern regShape *regCreateShape( regRegion *region, regMath glue,
				 regGeometry shape, regFlavor include, 
				 double *xpos, double *ypos, long npoints,
				 double *radius, double *angle );

extern regShape *regCreateNewShape(  regGeometry shape, regFlavor include, 
				     double *xpos, double *ypos, long npoints,
				     double *radius, double *angle );

extern regShape *regCreateWorldShape( regRegion *region, regMath glue,
				 regGeometry shape, regFlavor include, 
				 double *xpos, double *ypos, long npoints,
				 double *radius, double *angle, int wcoord, int wsize );

extern regShape *regCreateNewWorldShape(  regGeometry shape, regFlavor include, 
				     double *xpos, double *ypos, long npoints,
				     double *radius, double *angle, int wcoord, int wsize );

void regNegate(regShape* Shape );
int regInsidePoint( regShape *shape, double xpos, double ypos );
int regInsideAnnulus( regShape *shape, double xpos, double ypos );
int regInsideEllipse( regShape *shape, double xpos, double ypos );
int regInsideCircle( regShape *shape, double xpos, double ypos);
int regInsideBox( regShape *shape, double xpos, double ypos);
int regInsideRectangle( regShape *shape, double xpos, double ypos );
int regInsidePie( regShape *shape, double xpos, double ypos);
int regInsideSector( regShape *shape, double xpos, double ypos);
int regInsidePolygon( regShape *shape, double xpos, double ypos );
int regInsideField( regShape *shape, double xpos, double ypos );


char* regPrintVal( double x, char* buf );
void regPrintAngle( double x, char* buf );
void regPrintRadius( double r, int world, char* buf );
void regPrintPos( double x, int world, char* buf );
void regPrintPosPair( double x, double y, int world, char* xbuf, char* ybuf );


int regCaseEqual( char* s1, char* s2 );
regGeometry regShapeNameToGeometry( char* name );

long regShapeFindNPoints( regGeometry type, double* xpos, double* ypos, long nmax );
void regRotatedCoords( regShape* shape, double xp, double yp, double xcen, double ycen, double* xr, double* yr );
void regRotatedCoordsInvert( regShape* shape, double xr, double yr, double xcen, double ycen, double* xp, double* yp );
void regCornerBounds( double* xpos, double* ypos, double* xb, double* yb );
void regBoxCorners( regShape* shape, double* xpos, double* ypos );
int regRectangleCorners( regShape *shape, double* xpos, double* ypos );
int regExtentPolygon( regShape *shape, double* xpos, double* ypos );
int regExtentShape( regShape* shape, double* fieldx, double* fieldy, double* xpos, double* ypos );
int regRectangleSides( regShape *shape, double* xr, double* yr );
double regAreaBox( regShape* shape );
double regAreaRectangle( regShape* shape );
double regAreaCircle( regShape* shape );
double regAreaEllipse( regShape* shape );
double regShapeAnalyticArea( regShape* shape, double f, int* status );
double regModAngle( double ang );
void reg_parse_line( char* buf, long* mode, char** string, long* maxl, long* sys );
char* reg_token_advance( char* ptr, char* colname, char sep );
char* reg_tokens_advance( char* ptr, char* colname, char* sep );
int reg_read_line( FILE* fp, char* buf, long maxlen );
void reg_strcat( char** ptr, long* maxlen, char sep, char* buf );


void regApplyComponent( regRegion* region, regShape* shape, regMath glue );
regShape* regNextComponent( regShape* shape );
int regShapeOverlap( regShape* shape1, regShape* shape2 );
int regRectangleOverlap( double* xpos1, double* ypos1, double* xpos2, double* ypos2 );
void regFreeShape( regRegion* region, regShape* atShape );
