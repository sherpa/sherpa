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



#include "region_priv.h"
#include "regionio.h"
#include <string.h>


regRegion *myRegion;

main()
{
  dmBlock *block;


  regShape *myShape;
  long ii, jj;

  double xx[6], yy[6], angs[2], rad[2];
  short *ox;
  long dx, dy;
  char desc[20];;

  myRegion = regCreateRegion(NULL,NULL);

  regCreateAttribute( myRegion, "description", regCHAR, 20 );

  xx[0] = 10; yy[0]=10;
  rad[0] = 3; rad[1] = 11;
  angs[0] = 120;
  myShape = regCreateShape( myRegion, regOR, regELLIPSE, regInclude, xx, yy, 1, rad, angs );
  strcpy( desc, "atom" );
  regSetAttribute( myShape, "description", desc );

  angs[0] = 60;
  myShape = regCreateShape( myRegion, regAND, regELLIPSE, regExclude, xx, yy, 1, rad, angs );
  strcpy( desc, "example of" );
  regSetAttribute( myShape, "description", desc );

  angs[0] = 120;
  myShape = regCreateShape( myRegion, regOR, regELLIPSE, regExclude, xx, yy, 1, rad, angs );
  strcpy( desc, "XOR = " );
  regSetAttribute( myShape, "description", desc );

  angs[0] = 60;
  myShape = regCreateShape( myRegion, regAND, regELLIPSE, regInclude, xx, yy, 1, rad, angs );
  strcpy( desc, "A^B= A*!B + !A * B" );
  regSetAttribute( myShape, "description", desc );

  xx[0] = 27; yy[0] = 27; rad[0]=7;
  myShape = regCreateShape( myRegion, regOR, regCIRCLE, regInclude, xx, yy, 1, rad, angs );
  strcpy( desc, "Pacman" );
  regSetAttribute( myShape, "description", desc );

  
  angs[0] = 180-45; angs[1] = 180+45;
  myShape = regCreateShape( myRegion, regAND,  regPIE, regExclude, xx, yy, 1, rad, angs );
  strcpy( desc, "take a piece of pie" );
  regSetAttribute( myShape, "description", desc );


  xx[0] = 10; yy[0] = 25;
  rad[0] = 5; rad[1] = 12;
  angs[0] = 45;
  myShape = regCreateShape( myRegion, regOR,  regBOX, regInclude, xx, yy, 1, rad, angs );
  strcpy( desc, "MEG?" );
  regSetAttribute( myShape, "description", desc );



  xx[0] = 25; yy[0] = 5;
  xx[1] = 23; yy[1] = 10;
  xx[2] = 27; yy[2] = 13;
  xx[3] = 31; yy[3] = 10;
  xx[4] = 29; yy[4] = 5;
  xx[5] = 25; yy[5] = 5;
  myShape = regCreateShape( myRegion, regOR,  regPOLYGON, regInclude, xx, yy, 6, rad, angs);
  strcpy( desc, "the pentagon" );
  regSetAttribute( myShape, "description", desc );

 
  regSetAttribute( myShape, "description", desc );

  regCreateAttribute( myRegion, "BLAH", regFLOAT, 2 );
  regSetAttribute( myRegion->shape, "BLAH", xx );
  
  regRegionToMask( myRegion, -2, 35, -2, 35, 1,&ox, &dx, &dy );
  
  block = dmTableCreate( "test.fits[reg]" );
  regWriteRegion(  myRegion, block );
  dmTableClose( block );
  

  block = dmTableOpen( "test.fits[reg]" );
  myRegion = regReadRegion( block );
  dmTableClose( block );



  for (ii=dy-1;ii--;)
    {
      printf("%d\t", ii);
      for (jj=0;jj<dx;jj++)
        {
          printf("%s ", ox[ jj+ii*dx ] ? "*": "." );
        }
      printf("\n");
    }

  block = dmTableOpen( "./evts.fits[events][chipx=1:100]" );

  if ( block != NULL )
    {
      dmDescriptor *xcol = dmTableOpenColumn( block, "x" );
      dmDescriptor *ycol = dmTableOpenColumn( block, "y" );

      double data = dmGetScalar_d( xcol );


    }



  myShape = myRegion->shape;
  while ( myShape != NULL )
    {
      regPrintShape( myShape );
      myShape = myShape->next;
    }





  regWriteMask( ox, dx, dy, "test_mask.fits", &block );
  dmImageClose( block );


}
