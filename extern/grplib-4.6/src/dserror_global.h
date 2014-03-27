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


/************************************************************************/
/*  			GLOBAL DATA REQUIRED				*/
/************************************************************************/
#ifndef _DSERROR_GLOBAL_H
#define _DSERROR_GLOBAL_H

#ifndef _DSERROR_STRUCTS_H
#include "dserror_structs.h"
#endif

#include <stdio.h>

#ifndef _SETJMP_H
#include <setjmp.h>            /* for setjump and longjump */
#endif

#define	ERROR_MAXLINE	4096		/* max line length */


#define ERRLIBALLOCMSG "A memory allocation error has occurred within the ASC Error library.\n"

extern char *err_program_name; 	/* set to argv[0] by init_error_lib in main */
extern jmp_buf dserr_jmpbuf;

/* bit masks for the pipe/tools group */
#define dsPTGRPERR ((dsErrGroup) 0x000F)  /* this is all of the PT group errors together */
#define dsPTDMGRPERR ((dsErrGroup) 0x0001)
#define dsPTINSTRGRPERR ((dsErrGroup) 0x0002)
#define dsPTANALGRPERR ((dsErrGroup) 0x0004)
#define dsPTNONSIGRPERR ((dsErrGroup) 0x0008)
/* bit mask for the database group */
#define dsDBGRPERR ((dsErrGroup) 0x0010)
/* bit mask for the asc fitting engine */
#define dsASCFITERR ((dsErrGroup) 0x0020)

extern FILE *dsErrOutput_p; /* output is directed to this - by default it will
				be set to stderr. */

extern dsErr *dsErrHashMap_a;  /* Array of dsErr's.  This is the hash 
				  map, such that messages and severity's can
				  be found for an error code */
extern dsErrBool dsISLIBINIT; /* flag indicating if the error lib has been 
				 initialized */

extern long dsNUMBEROFERRORS; /* Number of errors in hash map */

extern short dsDEBUGLEV; /* debug level for the error library */

extern dsErrMemList dsErrAllLists_t; /* list of all lists created by user */

extern dsErrMemList dsErrAllInstances_t; /* list of all error instances 
					    created outside of lists */

#endif
