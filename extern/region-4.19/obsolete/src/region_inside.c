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

#include <math.h>
#include <float.h>

#include "region_priv.h"


int regInsideField( regShape *shape, double xpos, double ypos )
{

  int retval = 1;
  
  if ( shape->include == regInclude ) 
      return( retval );
  else
      return( 1 - retval );

}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideCircle( regShape *shape, double xpos, double ypos)
{
  int retval;
  double d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		  ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
 
  retval = (d <= shape->radius[0]) ? 1 : 0;
  if ( shape->include == regInclude ) 
      return( retval );
  else
      return( 1 - retval );
}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideAnnulus( regShape *shape, double xpos, double ypos)
{
  int retval;
  double d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		  ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
 	
 
  retval = (d <= shape->radius[1] && d >= shape->radius[0] ) ? 1 : 0;
  if ( shape->include == regInclude ) 
      return( retval );
  else
      return( 1 - retval );
}





/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideBox( regShape *shape, double xpos, double ypos)
{
  int retval;
  double dx;
  double dy;


  dx = shape->radius[0] / 2.0;
  dy = shape->radius[1] / 2.0;
  
  if ( shape->angle[0] != 0.0 )
    {
      /*double theta = shape->angle[0]*PI/180.0;*/

      double xm = xpos - shape->xpos[0];
      double ym = ypos - shape->ypos[0];
      
      double xp =  *shape->cos_theta * xm + *shape->sin_theta * ym;
      double yp = -*shape->sin_theta * xm + *shape->cos_theta * ym;

     /*      regRotatedCoordsInvert( shape, xpos, ypos, shape->xpos[0], shape->ypos[0], &xp, &yp ); */
     
      
      if ( ( fabs( xp ) <= dx )  &&
	   ( fabs( yp ) <= dy )  )
	     {
	  retval = 1;
	}
      else
	retval = 0;

      
    }
  else
    {

      
      if  ( ( xpos >= shape->xpos[0] - dx ) &&
	    ( xpos <= shape->xpos[0] + dx ) &&
	    ( ypos >= shape->ypos[0] - dy ) &&
	    ( ypos <= shape->ypos[0] + dy ) )
	{
	  retval = 1;
	}
      else
	retval = 0;
    }




  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideRectangle( regShape *shape, double xpos, double ypos )
{
  int retval;

  if ( shape->angle[0] == 0.0 )
    {
      if ( ( xpos >= shape->xpos[0] ) &&
	   ( xpos <= shape->xpos[1] ) &&
	   ( ypos >= shape->ypos[0] ) &&
	   ( ypos <= shape->ypos[1] ) )
	retval = 1;
      else
	retval = 0;

    }
  else
    {

      double xr[2], yr[2];
      double xpr, ypr;
      double xcen = ( shape->xpos[1] + shape->xpos[0] )/2.0;
      double ycen = ( shape->ypos[1] + shape->ypos[0] )/2.0;

      regRotatedCoords( shape, xpos, ypos, xcen, ycen, &xpr, &ypr );
      regRotatedCoords( shape, shape->xpos[0], shape->ypos[0], xcen, ycen, &xr[0], &yr[0] );
      regRotatedCoords( shape, shape->xpos[1], shape->ypos[1], xcen, ycen, &xr[1], &yr[1] );

      if ( ( xpr >= xr[0] ) &&
	   ( xpr <= xr[1] ) &&
	   ( ypr >= yr[0] ) &&
	   ( ypr <= yr[1] ) )
	retval = 1;
      else
	retval = 0;
      


    }

  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideSector( regShape *shape, double xpos, double ypos)
{
  int retval;
  double ang1 = fmod(shape->angle[0], 360.0);
  double ang2 = fmod(shape->angle[1], 360.0);
  double angat;

  if ( ang1 < 0.0 ) ang1 += 360.0;
  if ( ang2 < 0.0 ) ang2 += 360.0;

  angat = atan2( (ypos - shape->ypos[0]) , ( xpos  - shape->xpos[0] ) );
  angat *= 180.0/PI;

  if ( angat < 0.0) angat += 360.0;

  if ( ang1 < ang2 )
    {
      if ( ( angat >= ang1 ) && ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }
  else
    {
      if ( ( angat >= ang1 ) || ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }


  if ( ( xpos == shape->xpos[0] ) && ( ypos == shape->ypos[0] )) retval=1;


  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}

double regModAngle( double ang )
{
 double angle;
 angle = fmod( ang, 360.0 );
 if ( angle < 0.0 ) angle += 360.0;
 return angle;
}
/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsidePie( regShape *shape, double xpos, double ypos)
{
  int retval;
  double d;
  double ang1 = regModAngle(shape->angle[0]);
  double ang2 = regModAngle(shape->angle[1]);
  double angat;

  angat = atan2( (ypos - shape->ypos[0]) , ( xpos  - shape->xpos[0] ) );
  angat = regModAngle( angat * 180.0/PI );

  if ( ang1 < ang2 )
    {
      if ( ( angat >= ang1 ) && ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }
  else
    {
      if ( ( angat >= ang1 ) || ( angat <= ang2 ) )
	retval = 1;
      else
	retval = 0;
    }


  if (( retval ) && ( shape->radius) ) /* sectors have NULL radius */
  {
   d = sqrt( ( xpos - shape->xpos[0] ) * ( xpos - shape->xpos[0] ) +
		  ( ypos - shape->ypos[0] ) * ( ypos - shape->ypos[0] ) );
 	
 
   retval = (d <= shape->radius[1] && d >= shape->radius[0] ) ? 1 : 0;
  } 

  if ( ( xpos == shape->xpos[0] ) && ( ypos == shape->ypos[0] )) 
    {
      if (shape->radius) 
	if ( shape->radius[0] == 0.0 ) retval=1;
    }
  
  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsidePolygon( regShape *shape, double xpos, double ypos )
{

  long counter = 0;
  
  double ytest =0;
  long ii =0;
  int retval =0;
  
  /* x[n-1] = x[0]  and y[n-1] = y[0] */

  if (( ypos == shape->ypos[0] ) && ( xpos == shape->xpos[0] ))
    counter = 1;
  else
    
  for( ii=0; ii< shape->nPoints-1; ii++)
    {
      if (( ypos == shape->ypos[ii+1] ) && ( xpos == shape->xpos[ii+1] ))
	{
	  counter = 1;
	  break;
	}

      if (( shape->ypos[ii] < ypos ) && ( shape->ypos[ii+1] < ypos )) continue;

      if ( shape->xpos[ii] < shape->xpos[ii+1] )
	{
	  if ( ( xpos <= shape->xpos[ii] ) ||
	       ( xpos > shape->xpos[ii+1] ) ) continue;
	}
      else if ( shape->xpos[ii] > shape->xpos[ii+1] )
	{
	  if ( ( xpos <= shape->xpos[ii+1] ) ||
	       ( xpos > shape->xpos[ii] ) ) continue;
	}
      else if ( shape->xpos[ii] == shape->xpos[ii+1])
	{
	  if ( xpos != shape->xpos[ii] ) continue;
	  
	  if ( shape->ypos[ii+1] >= shape->ypos[ii] )
	    {
	      if ( (ypos >= shape->ypos[ii]) &&
		   (ypos <= shape->ypos[ii+1] ))
		{
		  counter=1;
		  break;
		}
	    }
	  else
	    {
	      if ( (ypos <= shape->ypos[ii]) &&
		   (ypos >= shape->ypos[ii+1] ))
		{
		  counter=1;
		  break;
		}

	    }
	  continue;
	}

      ytest =  shape->ypos[ii+1] -  ( shape->xpos[ii+1] - xpos ) *
	( shape->ypos[ii+1] - shape->ypos[ii] ) / 
	( shape->xpos[ii+1] - shape->xpos[ii] ) ;
      
      if ( ytest == ypos )
	{
	  counter = 1;
	  break;
	}
      else if ( ytest < ypos )
	{
	  continue;
	}
      
      counter += 1;


    }


  retval = counter % 2;



  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );

}




/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsidePoint( regShape *shape, double xpos, double ypos )
{
  int retval;
  
  if (  ( xpos == shape->xpos[0] ) && ( ypos == shape->ypos[0] ))
    retval = 1;
  else
    retval = 0;

  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );
}



/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideEllipse( regShape *shape, double xpos, double ypos )
{

  int retval=0;


  double xr, yr;



  if ( shape->angle[0] == 0 )
    {
      
      xr = ( xpos - shape->xpos[0] ) / shape->radius[0] ;
      yr = ( ypos - shape->ypos[0] ) / shape->radius[1] ;
      
    }
  else
    {
      
      double xm, ym;
      
      /*double theta = shape->angle[0] * PI / 180.0;*/
      
      xm = xpos - shape->xpos[0];
      ym = ypos - shape->ypos[0];
      
      xr = ( *shape->cos_theta * xm + *shape->sin_theta * ym)/shape->radius[0];
      yr = (-*shape->sin_theta * xm + *shape->cos_theta * ym)/shape->radius[1];
      
    }
  
  
  if ( ( xr * xr + yr * yr ) <= 1.0 ) retval = 1;
  
  
  if ( shape->include == regInclude )
    return(retval );
  else
    return( 1 - retval );
  

}

/*
  +-----------------------------------------------------------------------
  + Simple wrapper to return if a point is inside an individual shape
  +-----------------------------------------------------------------------
*/


int regInsideShape(
		   regShape *shape,
		   double xpos,
		   double ypos
		   )
{
  return( shape->inside( shape, xpos, ypos ) );
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
int regInsideRegion(
		    regRegion *region, 
		    double xpos,
		    double ypos
		    )
{

  regShape *atShape;
  int retval=0;
  int tmpval;
  int state;
  
  if ( !region )
   return 0;

  if ( ( xpos < region->xregbounds[0] ) || ( xpos > region->xregbounds[1] ) ||
       ( ypos < region->yregbounds[0] ) || ( ypos > region->yregbounds[1] ) )
    return(0);

  atShape = region->shape;
  while( atShape != NULL)
    {

      tmpval = 1;

      do {
	tmpval &= atShape->inside( atShape, xpos, ypos );

	state = 1;
	if ( atShape->next == NULL ) 
	  {
	    state = 0;
	  }
	else if ( atShape->next->component != atShape->component )
	  {
	    state = 0;
	  }

	atShape = atShape->next;
      } while ( state );

      retval |= tmpval;

    }
  
  return( retval );
}



void regRotatedCoordsInvert( regShape* shape, double xr, double yr, double xcen, double ycen, double* xp, double* yp )
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

void regRotatedCoords( regShape* shape, double xp, double yp, double xcen, double ycen, double* xr, double* yr )
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


