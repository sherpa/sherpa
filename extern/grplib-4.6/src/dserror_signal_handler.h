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


/*************************** signal_handler.h ********************************/
/* header for signal_catching routines */
#ifndef _DSERROR_SIGNAL_HANDLER_H
#define _DSERROR_SIGNAL_HANDLER_H

#ifndef _DSERROR_GLOBAL_H
#include "dserror_global.h"
#endif

/* Include statements */
#include <stdlib.h>
#include <stdio.h>
#include <signal.h>            /* for use by "signal" call in main */
#include <setjmp.h>            /* for setjump and longjump */

#ifndef _DSERROR_TRACE_FCT_H
#include "dserror_trace_fct.h" /* needs exit_upon_error */
#endif

extern void initialize_signal_handler( void );      /* function prototypes */
extern void signal_handler( int );           /* Establishes a signal handler */

#endif /* end _DSERROR_SIGNAL_HANDLER_H */
