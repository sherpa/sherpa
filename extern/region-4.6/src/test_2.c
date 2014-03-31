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

#include "region.h"


main()
{

  regRegion *myRegion;
  regShape *myShape;
  double xx[2]; 
  double yy[2];
  double angs[2];
  double rads[2];
  long component;
  short status;

  short *mask; long xlen; long ylen; short ii; short jj;

  char message[]="This is my circle";
  short ccd_id;

  /* create a region */
  myRegion = regCreateRegion(NULL, NULL);

  /* Create a simple region with single rectangle element */
  xx[0]=5; yy[0]=5;
  xx[1]=10; yy[1]=10;

  myShape = regCreateShape( myRegion, regOR,  regRECTANGLE, regInclude, xx, 
			    yy, 2, NULL, NULL );  
  regPrintShape( myShape );

  /* test to see if 8,8 is inside region */
  status = regInsideRegion( myRegion, 8, 8);
  printf( "Inside just rectangle = %s\n", status ? "yes": "no");

  /* Remove a circular region from the rectangle */
  xx[0]=7.5; yy[0]=7.5;
  rads[0]=2;
  myShape = regCreateShape( myRegion, regAND,  regCIRCLE, regExclude, xx, yy, 
			    1, rads, NULL );
  regPrintShape( myShape );

  /* test to see if 8,8 is inside combined region */
  status = regInsideRegion( myRegion, 8, 8);
  printf( "Inside combined region = %s\n", status ? "yes": "no");

  /* access shape's inside function directly */
  status = myShape->inside( myShape, 8, 8 );
  printf( "Inside circle region = %s\n", status ? "yes": "no");

  /* create a simple mask */
  regRegionToMask( myRegion, 1, 10, 1, 10, 1, &mask, &xlen, &ylen );
  for (ii=0;ii<xlen;ii++)
    {
      for (jj=0;jj<xlen;jj++)
        {
          printf( "%d ", mask[ii+jj*xlen] );
        }
      printf("\n");
    }

  /* create attributes*/
  regCreateAttribute( myRegion, "chip_id", regSHORT, 1 );
  regCreateAttribute( myRegion, "shape", regCHAR, 20 );

  /* set attributes */
  regSetAttribute( myShape, "shape", message );
  
  myShape=myRegion->shape; /* return to the rectangle */
  
  ccd_id = 10;
  regSetAttribute( myShape, "chip_id", &ccd_id );
  
  /* print attributes */
  regPrintShape( myShape );
  myShape = myShape->next;
  regPrintShape( myShape );


}

