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
#include  <ascdm.h>

extern int regWriteRegion(
                   regRegion *region,
                   dmBlock *outBlock
                   );


extern int regWriteMask(
		 short *mask,
		 long xlen,
		 long ylen,
		 char *file,
		 dmBlock **block
		 );

extern dmDataType regType2DM( regDataType inType );

extern regDataType regDM2Type( dmDataType inType );


extern regRegion *regReadRegion( dmBlock *inBlock );


