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


dmDataType regType2DM( regDataType inType )
{
  dmDataType retVal;


  switch ( inType )
    {
    case regCHAR:   retVal = dmTEXT;   break;
    case regUSHORT: retVal = dmUSHORT; break;
    case regSHORT:  retVal = dmSHORT;  break;
    case regULONG:  retVal = dmULONG;  break;
    case regLONG:   retVal = dmLONG;   break;
    case regFLOAT:  retVal = dmFLOAT;  break;
    case regDOUBLE: retVal = dmDOUBLE; break;
    default:
      break;
    }

  return( retVal );
}



regDataType regDM2Type( dmDataType inType )
{
  regDataType retVal;

  switch ( inType )
    {
    case dmTEXT:   retVal = regCHAR;   break;
    case dmUSHORT: retVal = regUSHORT; break;
    case dmSHORT:  retVal = regSHORT;  break;
    case dmULONG:  retVal = regULONG;  break;
    case dmLONG:   retVal = regLONG;   break;
    case dmFLOAT:  retVal = regFLOAT;  break;
    case dmDOUBLE: retVal = regDOUBLE; break;
    default:
      break;
    }

  return( retVal );
}
