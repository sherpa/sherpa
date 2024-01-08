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

#include "regionio.h"

#include <string.h>

/* angle and radius not necessarily always 2 */



regRegion *regReadRegion( dmBlock *inBlock )
{

  regRegion *myRegion;
  regShape *myShape;
  regGeometry myGeo;
  long shapeSize;


  char nameCol[68];

  dmDataType dm_type;
  regDataType reg_type;
  long auxsize;

  dmDescriptor *xposCol = NULL;
  dmDescriptor *yposCol = NULL;
  dmDescriptor *shapeCol = NULL;
  dmDescriptor *compCol = NULL;
  dmDescriptor *radCol = NULL;
  dmDescriptor *angCol = NULL;

  dmDescriptor **auxCol = NULL;

  char **auxdata = NULL;

  long ncols, ii;
  long naux = 0;


  double *xposData;
  double *yposData;
  double radData[2];
  double angData[2];
  long  compData;
  char  shapeData[68];
  char *shapeDataRip;
  long  maxSize;
  long  polySize;

  long  lastComp;

  long  angSize;
  long  radSize;


  myRegion = regCreateRegion(NULL,NULL);
  
  
  ncols = dmTableGetNoCols( inBlock );
  
  
  auxCol = ( dmDescriptor **)calloc( ncols, sizeof( dmDescriptor *) );

  auxdata = ( char **) calloc( ncols, sizeof( char * ));

 
  for (ii=0;ii<ncols; ii++)
    {
      auxCol[naux] = dmTableOpenColumnNo( inBlock, ii+1 );
      if ( auxCol[naux] == NULL )
	return(NULL);

      dmGetName( auxCol[naux], nameCol, 68 );

      if ( strcmp( nameCol, "COMPONENT" ) == 0 )
	compCol = auxCol[naux];
      else if ( strcmp( nameCol, "R" ) == 0 )
	radCol = auxCol[naux];
      else if ( strcmp( nameCol, "ROTANG" ) == 0 )
	angCol = auxCol[naux];
      else if ( strcmp( nameCol, "X" ) == 0 )
	xposCol = auxCol[naux];
      else if ( strcmp( nameCol, "Y" ) == 0 )
	yposCol = auxCol[naux];
      else if ( strcmp( nameCol, "SHAPE" ) == 0 )
	shapeCol = auxCol[naux];
      else
	{
	  long nBytes;

	  dm_type = dmGetDataType( auxCol[naux] );
	  
	  if( dm_type == dmTEXT )
	    auxsize = dmDescriptorGetLength( auxCol[naux] );
	  else
	    auxsize = dmGetArraySize( auxCol[naux] );

	  reg_type = regDM2Type( dm_type);
	  
	  regCreateAttribute( myRegion, nameCol, reg_type, auxsize );


	  switch (reg_type )
	    {
	    case regCHAR:
	      nBytes = sizeof(char); break;
	      
	    case regSHORT:
	    case regUSHORT:
	      nBytes = sizeof(short); break;

	    case regLONG:
	    case regULONG:
	      nBytes = sizeof(long); break;

	    case regFLOAT:
	      nBytes = sizeof(float ); break;
	      
	    case regDOUBLE:
	      nBytes = sizeof(double) ; break;

	    default:
	      break;

	    }

	  auxdata[naux] = (char *) calloc( auxsize, nBytes );
	  
	  naux += 1;
	}
      

    } /* end loop over columns */

  
  /* go thru rows */


  if (  (xposCol == NULL) || ( yposCol == NULL ) ) return( NULL );


  myRegion->xcol[0] = xposCol;
  myRegion->xcol[1] = yposCol;

  maxSize = dmGetArraySize( xposCol );

  /* +1 to handle polygons with last point not same as first */
  xposData = (double *)calloc( maxSize+1, sizeof(double));
  yposData = (double *)calloc( maxSize+1, sizeof(double ));



  lastComp = 0;


  if ( shapeCol == NULL ) strcpy( shapeData, "POINT" );
  if ( compCol == NULL ) compData = 1;
  if ( radCol == NULL ) 
    { 
      radData[0] = 0; radData[1] = 0; 
    }
  else
    {
      radSize = dmGetArraySize( radCol );
      if ( radSize < 2 ) radData[1] = 0;
      if ( radSize > 2 ) radSize = 2;
    }


  if ( angCol == NULL ) 
    { 
      angData[0] = 0; angData[1] = 0; 
    }
  else
    {
      angSize = dmGetArraySize( angCol );
      if ( angSize < 2 ) angData[1] = 0;
      if ( angSize > 2 ) angSize = 2;
    }



  
  do
    {

      regFlavor include;

      dmGetArray_d( xposCol, xposData, maxSize );
      dmGetArray_d( yposCol, yposData, maxSize );
      
      if ( shapeCol != NULL ) dmGetScalar_c( shapeCol, shapeData, 68 );
      if ( compCol != NULL )  compData = dmGetScalar_l( compCol );
      if ( radCol != NULL) dmGetArray_d( radCol, radData, radSize );
      if ( angCol != NULL ) dmGetArray_d( angCol, angData, angSize );


      if ( strchr( shapeData, '!' ) == NULL )
	{
	  include = regInclude;
	  shapeDataRip = shapeData;
	}
      else
	{
	  include = regExclude;
	  shapeDataRip = shapeData+1;
	}

      if ( strcmp( shapeDataRip, "POINT" ) == 0 )
	{
	  myGeo = regPOINT;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "CIRCLE" ) == 0 )
	{
	  myGeo = regCIRCLE;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "BOX" ) == 0 )
	{
	  myGeo = regBOX;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "ROTBOX" ) == 0)
	{
	  myGeo = regROTBOX;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "RECTANGLE" ) == 0 )
	{
	  myGeo = regRECTANGLE;
	  shapeSize = 2;
	}
      else if ( strcmp( shapeDataRip, "ROTRECTANGLE" ) == 0 )
	{
	  myGeo = regROTRECTANGLE;
	  shapeSize = 2;
	}
      else if ( strcmp( shapeDataRip, "ELLIPSE" ) == 0 )
	{
	  myGeo = regELLIPSE;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "PIE" ) == 0 )
	{
	  myGeo = regPIE;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "SECTOR" ) == 0 )
	{
	  myGeo = regSECTOR;
	  shapeSize = 1;
	}
      else if ( strcmp( shapeDataRip, "POLYGON" ) == 0 )
	{
	  long kk;

	  myGeo = regPOLYGON;

	  shapeSize = 0;

	  for( kk=1; kk<maxSize; kk++ )
	    {
	      if ( ( xposData[kk] == xposData[0] )  &&
		   ( yposData[kk] == yposData[0] ) )
		break;
	    }

	  if ( kk == maxSize )
	    {
	      xposData[maxSize] = xposData[0];
	      yposData[maxSize] = yposData[0];

	      shapeSize = maxSize +1;
	    }
	  else
	    {
	      shapeSize = kk ;
	    }


	}
      else
	{
	  return(NULL);
	}


      if ( compData != lastComp )
	{
	  lastComp = compData ;
	  myShape = regCreateShape( myRegion, regOR, 
				    myGeo, include, xposData, yposData, 
				    shapeSize,
				    radData, angData );
	}
      else
	{
	  myShape = regCreateShape( myRegion, regAND, 
				    myGeo, include, xposData, yposData, 
				    shapeSize,
				    radData, angData );
	}

      /* Get / set attributes */


      for (ii=0;ii<auxsize; ii++)
	{
	  switch ( myRegion->attributes[ii].type )
	    {
	    case regCHAR:
	      dmGetScalar_c( auxCol[ii], auxdata[ii], 
			     myRegion->attributes[ii].size  );
	      break;
	    case regSHORT:
	      dmGetArray_s( auxCol[ii], (short *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    case regUSHORT:
	      dmGetArray_us( auxCol[ii], (unsigned short *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    case regLONG:
	      dmGetArray_l( auxCol[ii], (long *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    case regULONG:
	      dmGetArray_ul( auxCol[ii], (unsigned long *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    case regFLOAT:
	      dmGetArray_f( auxCol[ii], (float *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    case regDOUBLE:
	      dmGetArray_d( auxCol[ii], (double *)auxdata[ii], 
			    myRegion->attributes[ii].size  );
	      break;

	    }

	  regSetAttribute( myShape, myRegion->attributes[ii].name, auxdata[ii] );


	  
	}



    } while( dmTableNextRow( inBlock ) != dmNOMOREROWS );


  return( myRegion );

}
