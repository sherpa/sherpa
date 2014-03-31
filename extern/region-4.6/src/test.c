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
  long ii, jj;

  double xx[6], yy[6], angs[2], rad[2];
  short *ox;
  long dx, dy;

  xx[0] = 10;
  yy[0] = 10;

  xx[1] = 5;
  yy[1] = 10;

  xx[2] = 15;
  yy[2] = 10;


  xx[3] = 20;
  yy[3] = 5;

  xx[4] = 10;
  yy[4] = 5;


  angs[0] = 120;
  angs[1] =  180;
  rad[0] = 3;
  rad[1]= 11;

  myRegion = regCreateRegion( NULL, NULL);

  myShape = regCreateShape( myRegion, regOR, regELLIPSE, regInclude, xx, yy, 1, rad, angs );

  angs[0] = 60;
  myShape = regCreateShape( myRegion, regAND, regELLIPSE, regExclude, xx, yy, 1, rad, angs );

  angs[0] = 120;
  myShape = regCreateShape( myRegion, regOR, regELLIPSE, regExclude, xx, yy, 1, rad, angs );

  angs[0] = 60;
  myShape = regCreateShape( myRegion, regAND, regELLIPSE, regInclude, xx, yy, 1, rad, angs );

  xx[0] = 27; yy[0] = 27; rad[0]=7;
  myShape = regCreateShape( myRegion, regOR, regCIRCLE, regInclude, xx, yy, 1, rad, angs );
  

  angs[0] = 180-45; angs[1] = 180+45;
  myShape = regCreateShape( myRegion, regAND, regPIE, regExclude, xx, yy, 1, rad, angs );


  xx[0] = 10; yy[0] = 25;
  rad[0] = 5; rad[1] = 12;
  angs[0] = 45;
  myShape = regCreateShape( myRegion, regOR, regBOX, regInclude, xx, yy, 1, rad, angs );

  xx[0] = 25; yy[0] = 5;
  xx[1] = 23; yy[1] = 10;
  xx[2] = 27; yy[2] = 13;
  xx[3] = 31; yy[3] = 10;
  xx[4] = 29; yy[4] = 5;
  xx[5] = 25; yy[5] = 5;
  myShape = regCreateShape( myRegion, regOR, regPOLYGON, regInclude, xx, yy, 6, rad, angs);


  regPrintShape( myShape );

  regRegionToMask( myRegion, 1, 35, 1, 35, 1,&ox, &dx, &dy );
  
  for (ii=dy-1;ii--;)
    {
      printf("%d\t", ii);
      for (jj=0;jj<dx;jj++)
	{
	  printf("%s ", ox[ jj+ii*dx ] ? "*": "." );
	}
      printf("\n");
    }


  


}
