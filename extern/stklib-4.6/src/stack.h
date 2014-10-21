/*_C_INSERT_SAO_COPYRIGHT_HERE_(1997,2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

#ifndef STKLIB_H     /* This code ensures that this file gets only
		      * included once */
#define STKLIB_H

/* Opaque type definition */

#ifdef STK_SRC
typedef struct {
 long current;
 long size;
 long nmax;
 char** data;
} StackData;
typedef StackData* Stack;
#else
typedef void* Stack;
#endif


/*!

\mainpage CXCDS stack library API

\section intro Introduction
This document is intended as a proposal for IRAF stack
implementation within C programs for Open--IRAF.  The
primary goal is to minimize user interaction with the
library.  IRAF stacks are not stacks in the more general
sense in Computer Science, but rather are lists or queues
commonly used to perform batch processing in IRAF.

The user of these library routines needs to create a
variable to hold the stack of type Stack.  All of the
tools to access this variable are provided in the library,
but the user will need to declare one such variable for
every stack their program uses.  This stack also marks the
user's current position in the stack so that users of the
library may choose to access the stack either via random
access or sequentially.

\section usage Usage

Users wanting to use functions in the stack lib need to

\code
#include "stack.h"
\endcode

Due to unfortunate naming this is very generic and may cause problems
when compiled against other libraries, eg STL.  Care should be taken
on the compile line to get the order of includes correct.

C++ users need to \code extern C \endcode the aforementioned header file.

 */

/*!
\example stktest/stktest.c
*/







/*!
  \struct Stack
  \brief This is opaque structure contains all the information abou the stack
  that has been built and being accessed.
  */


/*!
  @defgroup create Create Stacks
  Functions that can create stacks
  */

/*!
  @defgroup manip Manipulate Stack
  Function that can manipulate stacks
  */

/*!
  @defgroup access Access Stacks
  Functions that access stack elements
  */

/*!
  @defgroup destroy Delete Stacks
  Functions that destroy stacks

  */



/**
  \ingroup create

 * \brief Build a stack.
 */
/**
 * Build a stack from text
 * input.  The text input is a char* string passed in by
 * the user.  This string would most likely have been read in
 * from the task's parameter file, but this is not necessary.
 * This must be the first call made as it initializes and
 * creates the stack on which the other functions operate.  It
 * must be called separately for each stack to be created.
 *
 * \return Returns a Stack if build was successful.  NULL if not.
 *
 */
extern Stack    stk_build(
			  char *list /*!< i: Stack to build */
			  );
/*!
  Same as stk_build_gen without the prefix appending step.

  \return Returns a Stack if build was successful.  NULL if not.
  \ingroup create
 */
extern Stack    stk_build_gen(
			      char *list /*!< i: string to generate stack from */

			      );

/*!
  routine to count the number of
  entries in a stack.  Stacks are counted from 1 to N.


  \return Number of elements in stack
  \ingroup manip
  */
extern int      stk_count(
			  Stack stack /*!< i: Stack to count */
			  );

/*!
  Evaluates current position within the stack.  Stacks indexed from 1 to N.

  \return Current position in stack
  \ingroup manip
  */
extern int      stk_current(Stack stack /*!< i: Stack to get current count */
			    );

/*!
routine to append item(s)
indicated in a text string at the end of the stack.  

\return This
function returns  EXIT_SUCCESS if successful and EXIT_FAILURE if it fails.
\ingroup manip
 */
extern int      stk_append(Stack stack, /*!< i: Input stack */
			   char *descriptor /*!< i: Item to append */
			   );


/*!

  same as stk_append without the path append 

  \return EXIT_SUCCESS or EXIT_FAILURE
  \ingroup manip
  */
extern int      stk_append_gen(Stack stack, char *descriptor);

/*!
  Get the next string off the stack.

  \return char * to next string.
  \ingroup access
  */
extern char    *stk_read_next(Stack stack /*!< i: Stack to read from */
			      );

/*!
  Get the entry-th element off the stack.  If entry <=1 the first element is returned. 
  If entry >= max number of elements in the stack, the last element is returned.
  Stack are counted from 1 to N.

  \return char * to element in stack
  \ingroup access
 */
extern char    *stk_read_num(
			     Stack stack,  /*!< i: Stack to get data from */
			     int entry /*i!<: number to read */
			     );

/*!
  Free memory returned by stk_read_* routines.

  \return char * to element in stack
  \ingroup destroy
*/

extern void stk_read_free( char *name );



/*!
  Rewind the current stack back to the beginning
  \return char * to first element in stack
  \ingroup manip
  */
extern void     stk_rewind(
			   Stack stack /*!< i:  Stack to rewind */
			   );

/*!
  Routine to close a stack.  This is
  necessary to free up the memory used by the stack.  It is
  called when the stack is done being used.  
  \return The function
  returns EXIT_SUCCESS if successful and 
  EXIT_FAILURE if there is an error.
  \ingroup destroy

  */
extern int      stk_close(
			  Stack stack /*!< i: Stack to close */
			  );

/*!
  This routine will create a stack of numerically increasing name.  The input
  char * should have 1 (and only 1) "#" (which is not the first character).
  This routine will replace the "#" with numbers from 1 to int_suf_num

  \return a non-NULL Stack pointer if successful.
  \ingroup create
  */

extern Stack    stk_expand_n(char *list, long int_suf_num);  

/*!
  Deletes the num-th entry from the stack.  If num <= 1 .  If num >
  max number of ....
  Stack are counted from 1 to N.  If the last element of the stack is 
  deleted, the stack
  is NOT freeed.  

  \return error code
  \ingroup manip
 */
extern int      stk_delete_num( 
			       Stack stack,  /*!< i: Stack to remove entry from */
			       int num /*!< i: Element to delete */
			       );

/*!
  Deletes the current entry from the stack.
  If this was the last element on the stack,the stack
  is NOT freeed.  

  \return status
  \ingroup manip
  */  
extern int      stk_delete_current( 
				   Stack stack /*!< i: Stack to remove current element from */
				   );
/*!
  Set the stack to the num-th element.  If num <= 1.... if num > max ...
  \return status
  \ingroup manip
 */
extern int      stk_set_current( 
				Stack stack,  /*!< i: Stack to update */
				int num  /*!< i: Entry to set as current */
				);


extern int      stk_change_num(Stack stack, char *descriptor, int num);
extern int      stk_change_current(Stack stack, char *descriptor);



/*! \internal

  Internal functions not really for public consumption

  */
extern void     stk_disp( Stack stack );
extern void     stk_test( void );


/*  These functions do not exist 
    
    extern int      stk_delete_match(Stack stack, char *descriptor);
    extern int      stk_insert_num(Stack stack, char *descriptor, int num);
    extern int      stk_insert_current(Stack stack, char *descriptor);
    extern int      stk_sort(Stack stack, int az);
    */

#endif



