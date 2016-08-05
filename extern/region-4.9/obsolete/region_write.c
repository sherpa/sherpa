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
#define regNULL 0;


/* don't want to TableCreate since may want to append to existing file */

int regWriteRegion(
		   regRegion *region,
		   dmBlock *outBlock
		   )
{

  regShape *atShape;

  dmDescriptor *xposCol;
  dmDescriptor *yposCol;
  dmDescriptor *shapeCol;
  dmDescriptor *radiusCol;
  dmDescriptor *angleCol;
  dmDescriptor *compCol;

  dmDescriptor **auxCols;
  
  long maxDim = 0;

  char shape[20];
  double *xpos;
  double *ypos;
  double  ang[2];
  double  rad[2];

  long ii;


  atShape = region->shape;
  while ( atShape != NULL )
    {
      maxDim = ( maxDim > atShape->nPoints ) ? maxDim : atShape->nPoints;
      atShape = atShape->next;
    }


  xpos = (double *)calloc( maxDim, sizeof(double));
  ypos = (double *)calloc( maxDim, sizeof(double ));



  compCol = dmColumnCreate( outBlock, "COMPONENT", dmLONG, 0, "",
			    "component number" );
  shapeCol = dmColumnCreate( outBlock, "SHAPE", dmTEXT, 16, "", 
			     "shape geometry" );
  xposCol = dmColumnCreateArray( outBlock, "X", dmDOUBLE, 0, "", 
				 "X position", maxDim );
  yposCol = dmColumnCreateArray( outBlock, "Y", dmDOUBLE, 0, "", 
				 "Y position", maxDim );


  region->xcol[0] = xposCol;
  region->xcol[1] = yposCol;



  radiusCol = dmColumnCreateArray( outBlock, "R", dmDOUBLE, 0, "",
				   "radius",2 );
  angleCol = dmColumnCreateArray( outBlock, "ROTANG", dmDOUBLE, 0, "deg",
				  "rotation angle, ccw from X-axis",2 );


  auxCols = ( dmDescriptor**)calloc( region->numAttributes, sizeof( dmDescriptor*));

  for (ii=0;ii<region->numAttributes; ii++)
    {
      dmDataType dt = regType2DM( region->attributes[ii].type );
      if ( dt == dmTEXT )
	auxCols[ii] = dmColumnCreate( outBlock, region->attributes[ii].name,
					   dt,region->attributes[ii].size,
					   "", "");
					   
      else
	auxCols[ii] = dmColumnCreateArray( outBlock, region->attributes[ii].name,
					   dt, 0, "", "", 
					   region->attributes[ii].size );

    }

  dmKeyWrite_c( outBlock, "MFORM1", "X,Y", NULL, NULL);







  atShape = region->shape;
  while( atShape != NULL )
    {

      dmSetScalar_l(  compCol, atShape->component );


      for (ii=0;ii<maxDim;ii++)
	{
	  xpos[ii] = regNULL; ypos[ii] = regNULL;
	}
      ang[0] = regNULL; ang[1] = regNULL;
      rad[0] = regNULL; rad[1] = regNULL;



      switch ( atShape->shape )
	{
	case regPOINT:
	  sprintf( shape, "%sPOINT", atShape->include ? "" : "!" ); 
	  xpos[0] = atShape->xpos[0];
	  ypos[0] = atShape->ypos[0];
	  break;

	case regBOX:
	  if ( atShape->angle[0] == 0.0 )
	    sprintf( shape, "%sBOX", atShape->include ? "" : "!" ); 
	  else
	    sprintf( shape, "%sROTBOX", atShape->include ? "" : "!" ); 

	  xpos[0] = atShape->xpos[0];
	  ypos[0] = atShape->ypos[0];
	  rad[0] = atShape->radius[0];
	  rad[1] = atShape->radius[1];
	  ang[0] = atShape->angle[0];
	  break;
	  
	case regRECTANGLE:
	  if ( atShape->angle[0] == 0.0 )
	    sprintf( shape, "%sRECTANGLE", atShape->include ? "" : "!" ); 
	  else
	    sprintf( shape, "%sROTRECTANGLE", atShape->include ? "" : "!" ); 
	  xpos[0] = atShape->xpos[0]; xpos[1] = atShape->xpos[1];
	  ypos[0] = atShape->ypos[0]; ypos[1] = atShape->ypos[1];
	  ang[0] = atShape->angle[0];
	  break;
	  
	case regCIRCLE:
	  sprintf( shape, "%sCIRCLE", atShape->include ? "" : "!" ); 
	  xpos[0] = atShape->xpos[0];
	  ypos[0] = atShape->ypos[0];
	  rad[0] = atShape->radius[0];
	  break;
	  
	case regELLIPSE:
	  sprintf( shape, "%sELLIPSE", atShape->include ? "" : "!" ); 
	  xpos[0] = atShape->xpos[0];
	  ypos[0] = atShape->ypos[0];
	  rad[0] = atShape->radius[0];
	  rad[1] = atShape->radius[1];
	  ang[0] = atShape->angle[0];
	  break;
	  
	case regPOLYGON:
	  sprintf( shape, "%sPOLYGON", atShape->include ? "" : "!" ); 
	  for (ii=0;ii<atShape->nPoints;ii++)
	    {
	      xpos[ii] = atShape->xpos[ii]; 
	      ypos[ii] = atShape->ypos[ii];
	    }

	  break;
	  
	case regPIE:
	  sprintf( shape, "%sPIE", atShape->include ? "" : "!" ); 
	  xpos[0] = atShape->xpos[0];
	  ypos[0] = atShape->ypos[0];
	  ang[0] = atShape->angle[0];
	  ang[1] = atShape->angle[1];
	  break;
	}

      if(shape->type==regMASK || shape==regUSER)
	continue;
      
      dmSetScalar_c( shapeCol, shape );
      dmSetArray_d( xposCol, xpos, maxDim );
      dmSetArray_d( yposCol, ypos, maxDim );
      dmSetArray_d( radiusCol, rad, 2 );
      dmSetArray_d( angleCol, ang, 2 );



      for (ii=0;ii< region->numAttributes;ii++)
	{
	  if ( atShape->attributes[ii] != NULL )
	    {
	      
	      switch( region->attributes[ii].type )
		{
		case regCHAR:
		  dmSetScalar_c( auxCols[ii], atShape->attributes[ii] );
		  break;
		case regSHORT:
		  dmSetArray_s( auxCols[ii], (short *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;

		case regUSHORT:
		  dmSetArray_us( auxCols[ii], (unsigned short *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;

		case regLONG:
		  dmSetArray_l( auxCols[ii], (long *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;

		case regULONG:
		  dmSetArray_ul( auxCols[ii], (unsigned long *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;

		case regFLOAT:
		  dmSetArray_f( auxCols[ii], (float *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;

		case regDOUBLE:
		  dmSetArray_d( auxCols[ii], (double *)atShape->attributes[ii],
				region->attributes[ii].size );
		  break;



		}

	    }

	}

  
      dmTablePutRow( outBlock, NULL);
    next:
      atShape = atShape->next;
    }


  free(xpos); free(ypos);
      
  return(0);

}









int regWriteMask(
		 short *mask,
		 long xlen,
		 long ylen,
		 char *file,
		 dmBlock **block
		 )
{

  long lAxes[2]; 
  long ll[2] = { 1, 1};
  dmDescriptor *image;

  lAxes[0] = xlen; lAxes[1] = ylen;

  *block = dmImageCreate( file, dmSHORT, lAxes, 2 );

  image = dmImageGetDataDescriptor( *block );

  dmImageDataSetSubArray_s( image, ll, lAxes, mask );


 return 0;

}
