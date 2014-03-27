/*                                                                
**  Copyright (C) 1997-2007  Smithsonian Astrophysical Observatory 
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


#ifndef _DSERROR_TRACE_FCT_H
#define _DSERROR_TRACE_FCT_H

#ifndef _DSERROR_GLOBAL_H
#include "dserror_global.h"
#endif

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef DEBUG_FLAGS_H
#define DEBUG_FLAGS_H
typedef struct DEBUG_FLAGS {
  unsigned short int display;     /* how much debug messages to print */
  unsigned short int plot_result; /* how much plotting is to be done */
} DEBUG_FLAGS;
#endif /*DEBUG_FLAGS_H*/

#define NUM_BLANK_SPACES "   "

extern void init_function_stack(char *name, 
				int print_it, 
				int num_fct_to_print);

extern void check_in(const char *name);

extern void check_out(void);

extern int exit_upon_error(long exit_code, 
			   char *format, ...);

extern void print_debug_msg(int print_level,
			    char *format, ...);

extern char *save_string(char *name);

extern void print_function_stack(FILE *fp);
#endif /* end _DSERROR_TRACE_FCT_H */
