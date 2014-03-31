/*                                                                
**  Copyright (C) 1997,2007  Smithsonian Astrophysical Observatory 
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

#define CFAsccsID "%W% %E% %U%"

/*-------------------------------------------------------------------------

   FILE NAME:		%M%

   DEVELOPMENT: 	AXAF Science Center DataModel
 			mnoble@cfa.harvard.edu (Initial version)

   DESCRIPTION:		ASC DataModel filter parser/grammar test program.
   			Use 'quit' token to exit, and 'debug [on|off]' to 
			trace the parsing.

---------------------------------------------------------------------------
---------------------------------------------------------------------------*/

#include "cxcregion.h"
#include <stdio.h>
#include <stdlib.h>
/*-------------------------------------------------------------------------*/


void main(int argc, char **argv)
{
 int v;
 if (argc>2) v = atol( argv[2] );
  /* Allow one command line argument, an input file to be processed */
  regReadAsciiRegion( argv[1], v );
}
