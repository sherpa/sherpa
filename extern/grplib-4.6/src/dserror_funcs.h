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


/**********************dserror_funcs.h******************************/
/* Header for ASC error routines, to be included *after* all standard*/
/*  system headers */

#ifndef	__dserror_funcs_h
#define	__dserror_funcs_h

#ifndef _DSERROR_GLOBAL_H
#include "dserror_global.h"
#endif

#ifndef _DSERROR_STRUCTS_H
#include "dserror_structs.h"
#endif

#ifndef _DSERROR_GENERAL_H
#include "dserror_general.h"	/* General errors defined here */
#endif


/* the following are the 4 categories of errors for the 
   pipe/tools team's errors */
#ifndef _DSERROR_PTANALYSIS_H
#include "dserror_ptanalysis.h"
#endif

#ifndef _DSERROR_PTDM_H
#include "dserror_ptdmtools.h"
#endif

#ifndef _DSERROR_PTNONSI_H
#include "dserror_ptnonsi.h"
#endif

#ifndef _DSERROR_PTINSTRUMENT_H
#include "dserror_ptinstrument.h"
#endif

#ifndef _DSERROR_DATABASE_H
#include "dserror_database.h"   /* database team errors defined here */
#endif

#ifndef _SETJMP_H
#include <setjmp.h>            /* for setjump and longjump */
#endif
                            /* 4.3BSD Reno <signal.h> doesn't define SIG_ERR */
#ifndef _DSERROR_SIGNAL_HANDLER_H
#include "dserror_signal_handler.h"
#endif

#ifndef _DSERROR_TRACE_FCT_H
#include "dserror_trace_fct.h"
#endif

#if	defined(SIG_IGN) && !defined(SIG_ERR)
#define	SIG_ERR	((Sigfunc *)-1)
#endif

void	err_dump(const char *, ...);	/* {App misc_source} */
void	err_msg(const char *, ...);
void	err_quit(long, const char *, ...);
void	err_ret(const char *, ...);
void	err_exit(long, const char *, ...);
int     init_error_lib(char *);

/*****************************************************************************/
/* Additions to the error library, post February 1998			     */
/* Author: DLM								     */
/*****************************************************************************/

/*****************************************************************************/
/*  			FUNCTIONS REQUIRED				     */
/*****************************************************************************/

/* Initialization routines */

/* Initialization function for error library, it is a wrapper around the 
   existent one */
extern dsErrCode dsErrInitLib(dsErrGroup error_groups_t,
			      char *tool_name_a);

/* Will return dsErrTrue if lib already initialized */
extern dsErrCode dsErrIsInitialized(void);

/* redirect output to a user specified FILE *.  This could be a logfile,
   or stderr, stdout */
extern void dsErrDirectOutput(FILE *output_p);

/* will free any memory allocated for lists, including any lists created,
   and will free memory in the hash map */
extern void dsErrCloseLib(void);

/* Initialize elements of an error instance.  */
extern dsErrCode dsErrCreateInstance(dsErrInstance **error_instance_p);

/* Reset elements of an error instance */
extern void dsErrResetInstance(dsErrInstance *error_instance_p);
 
/* allocate memory for error lists.  Elements will be initialized. */
extern dsErrCode dsErrCreateList(dsErrList **error_list_p);

/* Addition routines */

/*  Add an error onto the list.  Add it at the end if it is an
    individual error, or the first instance of an accumulation.  Otherwise
    increment the counter in the last instance of the accumulation found
    (it should be the only instance).  Since it might not be able
    to allocate memory for a new node, it needs to return an error
    code(ironic).  Specify whether the error message is to be the generic
    one or a customized one.  Input custom message or element to fill out
    templates in the default generic message. */
extern dsErrCode dsErrAdd(dsErrList *error_list_p,
			  dsErrCode error_code_t,
			  dsErrType error_type_e,		   
			  dsErrMsgType msg_type_e,
			  ...);

/* populate an error instance structure.  This task is similar to dsErrAdd, 
   in that the user specifies the input code, type, message type and message
   elements.  An input error structure will be filled in, following the same
   rules as for dsErrAdd.
*/
extern dsErrCode dsErrSetInstance(dsErrCode error_code_t,
				  dsErrType error_type_e,
				  dsErrMsgType msg_type_e,
				  dsErrInstance *error_instance_p,
				  ...);

/* Add an already populated error instance to the end of the list (the most
   recent errors are added to the end).  The user must have already 
   set up the elements of the error_to_add */
extern dsErrCode dsErrAddToEnd(dsErrList *error_list_p,
			       dsErrInstance *error_to_add_p);


/* Search and retrieval routines */

/* Find and return the data pertaining to the first error in the list */
extern dsErrCode dsErrLook(dsErrList *error_list_p, 
			   dsErrInstance *out_instance_p);

/* Find and return the data pertaining to the Nth error in the list */
extern dsErrCode dsErrLookN(dsErrList *error_list_p,
			    long N, 
			    dsErrInstance *out_instance_p);

/* Find and return the data pertaining to the first instance of an
   error that matches the severity level */
extern dsErrCode dsErrLookSev(dsErrList *error_list_p, 
			      dsErrSeverity err_severity_e,
			      dsErrInstance *out_instance_p);
 
/* Find and return the data pertaining to the first instance of an
   error that matches the error code and type */
extern dsErrCode dsErrLookCode(dsErrList *error_list_p,
			       dsErrCode error_code_t, 
			       dsErrType error_type_e,
			       dsErrInstance *out_instance_p);
 
/* Find and return the data pertaining to the first instance of an
   error that does not matche the error code and type */
extern dsErrCode dsErrLookCodeExcl(dsErrList *error_list_p,
				   dsErrCode error_code_t, 
				   dsErrType error_type_e,
				   dsErrInstance *out_instance_p);

/* Find instances of errors that match the severity level, and put
them into the output list, in the order they are in the input list */
extern dsErrCode dsErrLookAllSev(dsErrList *in_list_p,
				 dsErrSeverity err_severity_e, 
				 dsErrList *out_list_p);
 
/* Find instances of errors that match the error code and type, and put
them into the output list, in the order they are in the input list */
extern dsErrCode dsErrLookAllCode(dsErrList *in_list_p,
				  dsErrCode error_code_t,
				  dsErrType error_type_e,
				  dsErrList *out_list_p);

/* Find instances of errors that do not match the error code and type, and put
them into the output list, in the order they are in the input list */
dsErrCode dsErrLookAllCodeExcl(dsErrList *in_list_p,
			       dsErrCode error_code_t,
			       dsErrType error_type_e,
			       dsErrList *out_list_p);

/* count number of occurences of an error within a list.  */
extern long dsErrGetNumOccur(dsErrList *error_list_p,
			     dsErrCode error_code_t,
			     dsErrType error_type_e);

/* Find and return a pointer to the data pertaining to the first error in the list */
extern dsErrInstance *dsErrPeek(dsErrList *error_list_p);

/* If an instance of an error matches the error code and type, return a pointer to
   the instance, else return NULL. */
extern dsErrInstance *dsErrPeekCode(dsErrList *error_list_p,
				    dsErrCode error_code_t, 
				    dsErrType error_type_e);

/* Find and return a pointer to the data pertaining to the first instance of an
   error that does not match the error code and type */
extern dsErrInstance *dsErrPeekCodeExcl(dsErrList *error_list_p,
					dsErrCode error_code_t, 
					dsErrType error_type_e);

/* The address of the first instance that matches the severity is
   returned, otherwise NULL is returned. */
extern dsErrInstance *dsErrPeekSev(dsErrList *error_list_p, 
				   dsErrSeverity err_severity_e);

/* Return the address of the Nth error in the list */
extern dsErrInstance *dsErrPeekN(dsErrList *error_list_p,
				 long N);

/* Error removal routines */

/* Remove the first error in the error list */
extern dsErrBool dsErrRemove(dsErrList *error_list_p);
 
/* Remove the Nth error in the list */
extern dsErrBool dsErrRemoveN(dsErrList *error_list_p,
			      long N);

/* Remove the first matching instance of an error based on severity */
extern dsErrBool dsErrRemoveSev(dsErrList *error_list_p,
				dsErrSeverity err_severity_e);

/* Remove the first matching instance of an error based on the 
   error code, and it's type */
extern dsErrBool dsErrRemoveCode(dsErrList *error_list_p, 
				 dsErrCode error_code_t, 
				 dsErrType error_type_e);
 
/* Remove all errors in the list */
extern dsErrBool dsErrRemoveAll(dsErrList *error_list_p);

/* Remove all matching errors based on severity */
extern dsErrBool dsErrRemoveAllSev(dsErrList *error_list_p, 
				   dsErrSeverity err_severity_e);
 
/* Remove all matching errors based on error code and type. */
extern dsErrBool dsErrRemoveAllCode(dsErrList *error_list_p, 
				    dsErrCode error_code_t, 
				    dsErrType error_type_e);

/* Error Processing functions */

/* Prints out all the elements of an error list. */
extern void dsErrPrintList(dsErrList *error_list_p,
			   dsErrBool print_accum_e);

/* Prints out the instances error message */
extern void dsErrPrintInstance(dsErrInstance *error_instance_p,
			       dsErrBool print_accum_e);
 
/* Prints out the instances error message, and returns */
extern void dsErrReturnInstance(dsErrInstance *error_instance_p,
				dsErrBool print_accum_e);
 
/* Prints out the instances error message, and exits */
extern void dsErrExitInstance(dsErrInstance *error_instance_p,
			      dsErrBool print_accum_e);
 
/* Prints out the instances error message, and quits */
extern void dsErrQuitInstance(dsErrInstance *error_instance_p,
			      dsErrBool print_accum_e);
 
/* Prints out the instances error message, and dumps core */
extern void dsErrDumpInstance(dsErrInstance *error_instance_p,
			      dsErrBool print_accum_e);

/* Accessor functions */

/* get total number of errors(nodes) in the list */
extern long dsErrGetErrorCt(dsErrList *error_list_p);

/* get the number of fatal errors (nodes) in the list. */
extern long dsErrGetFatalCt(dsErrList *error_list_p);

/* get the number of warnings (nodes) in the list */
extern long dsErrGetWarningCt(dsErrList *error_list_p);

/* get the accumulation count of an error instance.  For individual
   errors this is one, by definition */
extern long dsErrGetInstCt(dsErrInstance *error_instance_p);

/* get the error code of an error instance. */
extern dsErrCode dsErrGetInstCode(dsErrInstance *error_instance_p);

/* get the severity of an error instance. */
extern dsErrSeverity dsErrGetInstSev(dsErrInstance *error_instance_p);

/* get the message of an error instance */
extern dsErrMsg dsErrGetInstMsg(dsErrInstance *error_instance_p);
 
#endif	/* __dserror_funcs_h */

