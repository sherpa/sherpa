/*                                                                
**  Copyright (C) 1998,2007  Smithsonian Astrophysical Observatory 
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


#ifndef DSERROR_STRUCTS_H
#include "dserror_structs.h"
#endif
#define _DSERROR_DATABASE_H

#define dsDATABASEERROFFSET     -8000
#define dsDATABASENUMERRORS     2

#define dsDATABASEEXISTSERR (dsDATABASEERROFFSET - 1)
#define dsDATABASEEXISTSSEV dsERRSEVFATAL
#define dsDATABASEEXISTSSTDMSG "Database %s does not exist\n"

#define dsATTNOTFOUNDERR (dsDATABASEERROFFSET - 2)
#define dsATTNOTFOUNDSEV dsERRSEVFATAL
#define dsATTNOTFOUNDSTDMSG "That attribute was not found in database %s\n"
