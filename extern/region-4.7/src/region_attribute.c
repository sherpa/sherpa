/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "region_priv.h"

void regFreeShape( regRegion* region, regShape* atShape )
{
  if( !atShape )
    return;
  
  if(atShape->user)
    atShape->user->free( atShape);
  else {
    if ( atShape->xpos ) free( atShape->xpos );
    if ( atShape->ypos ) free( atShape->ypos );
    if ( atShape->angle ) free( atShape->angle );
    if ( atShape->radius ) free( atShape->radius );
    if ( atShape->sin_theta) free( atShape->sin_theta);
    if ( atShape->cos_theta) free( atShape->cos_theta);
  }
  free( atShape );    
  atShape=NULL;

}


/* Attribute stuff */

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
    }
 free( region );
}

