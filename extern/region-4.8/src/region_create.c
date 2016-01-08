/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007,2013)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "cxcregion.h"
#include "region_priv.h"
#include <float.h>


/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
regRegion* regCreateEmptyRegion( void )
{
  return regCreateRegion( NULL,NULL );
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

regRegion* regCreateRegion( void *xcol, void*ycol )
{
  regRegion *region = ( regRegion *) calloc( 1, sizeof( regRegion ) );  
  return( region );  
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

regRegion* regCopyRegion( regRegion* inRegion )
{
  double fx[2] ={ -DBL_MAX, DBL_MAX };
  double fy[2] ={ -DBL_MAX, DBL_MAX };
  regRegion* Region;
  regShape* Shape;
  regShape* inShape;
  regMath glue;
  
  long lastComponent = 1;
  
  if ( !inRegion )
  return NULL;
  
  Region = regCreateRegion(NULL, NULL);
  
  inShape = inRegion->shape;
  
  while (inShape != NULL ) {
    Shape = regCopyShape( inShape );
    if ( inShape->component == lastComponent ) {
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

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/

void regFree( regRegion* region )
{
  regShape *atShape;
  regShape* shape;
  /* Free shapes */
  
  if ( !region )
    return;
  
  atShape = region->shape;

  /* Free shape attributes */
  while (atShape != NULL )
    {
      shape = atShape;
      atShape = atShape->next;
      regFreeShape( region, shape );
      shape = NULL;
      
    }

  free( region );
  region = NULL;
}

