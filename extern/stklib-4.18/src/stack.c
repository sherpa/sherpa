/*                                                                
**  Copyright (C) 1998-2007,2009,2016-2017,2024  Smithsonian Astrophysical Observatory 
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


/* **************************************************************************
 * FILE NAME: stack.c
 * 
 * DEVELOPEMENT: ASCDS (Data Systems)
 * 
 * DESCRIPTION:
 * 
 * This is a library of functions that can be used to handle IRAF stacks.  By
 * "IRAF stack", I mean a list or queue of text elements specified by the
 * user of a task.  This stack capability will be modelled after the behavior
 * of these stacks in the older, preceding IRAF tasks that run in the cl.
 * 
 * NOTES:
 * 
 * <None>
 * 
 * REVISION HISTORY:
 * 
 * Ref. No.        Date
 * --------        ----
 * 
 * 1.12             05/07/96
 * 2.0              1998 Jul 14 JCM
 *                  03/03/00 
 * ************************************************************************ */

#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/* Include ctype for the isspace function */
#include <ctype.h>
#define STK_SRC
#include "stack.h"
#include "logical.h"
#include <math.h> 

#define STK_INIT_SIZE 100
#define MAX_ITEM_SIZE 1024

/* The old STK_STEP was a comma, but I think semi-colon would be safer */

#define STK_SEP ','
#define STK_SEP_ALT ';'
#define STK_ESC '\\'

/* Declarations of local routines */
static Stack   stk_alloc( long size );
static int     stk_append_entry( Stack stack, char* item );
static int     stk_append_item( Stack stack, char* item,int prepend );
static int     stk_append_prepend(Stack stack, char *entry, int prepend );
static Stack   stk_build_prepend(char *inputList, int prepend );
static char*   stk_cat_string( char* a, char* b );
static char*   stk_copy_string( char* ptr );
static logical stk_next_list_item( char** pptr, char* itembuf, long maxlen );
static void    stk_trim( char* opt );

/* Wrapper methods to support testing of static methods */
Stack          mock_stk_alloc(long size);
int            mock_stk_append_item(Stack stack, char* entry, int prepend);

/*********************************************************************
 * Internal:
 * Allocates memory for a stack.
 *********************************************************************/
static Stack stk_alloc( long size )
{
  Stack stack;
  
  stack = (Stack)malloc( sizeof( StackData ) );
  if ( stack ) 
  {
    
    stack->data = (char**)malloc( size * sizeof( char* ) );
    if ( stack->data ) 
    {
      stack->nmax = size;
      stack->current = 0;
      stack->size = 0;
    }
    else 
    {
      free( stack );
      stack = NULL;
    }
  }
  if ( !stack ) 
  {
    fprintf(stderr, "ERROR: not enough memory to allocate stack\n");
  }
  return( stack );
}

Stack mock_stk_alloc( long size )
{
  return stk_alloc(size);
}

/*********************************************************************
 * Internal:
 * Add a new entry to the stack, reallocating the stack if it needs
 * more memory.
 *********************************************************************/
static int stk_append_entry( Stack stack, char* item )
{
  float xx, xmin, xmax, xdel;
  float yy, ymin, ymax, ydel;
  float xpos, ypos;

  char *region;
  char *ptr;
  char *eqs;

  int rect=0;
  int pol=0;
  int lin=0;

  enum { UNKNOWN, REAL, INTEGER } r_or_i;
  /*   short r_or_i = 0;  real or integer */

  /* 100 is to accomodate the extra region info, may not be 
     long enough in all cases though, could be more robust */
  if ( item && strlen(item) ) 
  {
    region = (char*)calloc(strlen(item)+100,sizeof(char));
  }
  else 
  {
    region = NULL;
  }
    
  if ( (ptr = strstr( item, "rgrid(") ) != NULL )
  {
    eqs = strchr(ptr, ')');
    if ( eqs == NULL )
    {
      xmin=0; xmax=1; xdel=2;
      ymin=0; ymax=1; ydel=2;
      ptr = NULL;
    }
    else
    {
      strncpy( region, ptr, eqs-ptr+1 );
      
      rect = sscanf( region, "rgrid(%f:%f:%f,%f:%f:%f)",
		     &xmin, &xmax,&xdel, &ymin, &ymax, &ydel);
      if ( rect != 6 )
      {
	xmin=0; xmax=1; xdel=2;
	ymin=0; ymax=1; ydel=2;
	ptr = NULL;
      }
      else 
      {
	xmax -= xdel;
	ymax -= ydel;
      }
      r_or_i = REAL;
    }
  }
  else if ( (ptr = strstr(item, "pgrid(")) != NULL )
  {
    eqs = strchr(ptr, ')');
    if ( eqs == NULL )
    {
      xmin=0; xmax=1; xdel=2;
      ymin=0; ymax=1; ydel=2;
      ptr = NULL;
    }
    else
    {
      strncpy( region, ptr, eqs-ptr+1 );
      
      pol = sscanf( region, "pgrid(%f,%f,%f:%f:%f,%f:%f:%f)",
		    &xpos, &ypos, &xmin, &xmax,&xdel, &ymin, &ymax, &ydel);
      if ( pol != 8 )
      {
	xmin=0; xmax=1; xdel=2;
	ymin=0; ymax=1; ydel=2;
	ptr = NULL;
      }
      else 
      {
	xmax -= xdel;
	ymax -= ydel;
      }
      r_or_i = REAL;
    }
  }
  else if ( (ptr = strstr(item, "lgrid(")) != NULL )
  {
    eqs = strchr(ptr, ')');
    if ( eqs == NULL )
    {
      xmin=0; xmax=1; xdel=2;
      ymin=0; ymax=1; ydel=2;
      ptr = NULL;
    }
    else
    {
      strncpy( region, ptr, eqs-ptr+1 );
      
      lin  = sscanf( region, "lgrid(%f:%f:%f)", &xmin, &xmax, &xdel );
      if ( lin != 3 )
      {
	xmin=0; xmax=1; xdel=2;
	ymin=0; ymax=1; ydel=2;
	ptr = NULL;
      }
      ymin = 0; ymax=1; ydel=2;
      /* xmax += xdel; */
      r_or_i = REAL;
    }
  }
  else if ( (ptr = strstr(item, "igrid(")) != NULL )
  {
    eqs = strchr(ptr, ')');
    if ( eqs == NULL )
    {
      xmin=0; xmax=1; xdel=2;
      ymin=0; ymax=1; ydel=2;
      ptr = NULL;
    }
    else
    {
      strncpy( region, ptr, eqs-ptr+1 );
      
      lin  = sscanf( region, "igrid(%f:%f:%f)", &xmin, &xmax, &xdel );
      if ( lin != 3 )
      {
	xmin=0; xmax=1; xdel=2;
	ymin=0; ymax=1; ydel=2;
	ptr = NULL;
      }
      ymin = 0; ymax=1; ydel=2;
      /* xmax += xdel;*/
      r_or_i = INTEGER;
    }
  }
  else
  {
    xmin=0; xmax=1; xdel=2;
    ymin=0; ymax=1; ydel=2;
    ptr = NULL;
  }
  
  for (xx=xmin; xx<=xmax; xx+=xdel )
  {
    for (yy=ymin; yy<=ymax; yy+=ydel )
    {
      if (region)
	strcpy( region, item );	  
	  
      if ( ( ptr != NULL )  && ( *item != '@' ))
      {
	char filter[MAX_ITEM_SIZE];
	if ( rect )
	  sprintf( filter, "rectangle(%g,%g,%g,%g)",
		   xx,yy, xx+xdel, yy+ydel );
	else if ( pol ) 
	  sprintf( filter, "pie(%g,%g,%g,%g,%g,%g)",
		   xpos,ypos,xx,xx+xdel, yy, yy+ydel );
	else { /* assume linear */
	  if ( r_or_i == REAL ) 
	    sprintf( filter, "%g", xx );
	  else
	    sprintf( filter, "%ld", (long)xx );
	}
	
	/* need to combine filters */
	if ( region ) 
	{
	  char *blank;
	  blank = region;
	  while ( *blank ) { *blank = '\0'; blank++;}
	  strncpy( region, item, ptr-item);
	  strcat( region, filter );
	  strcat(region, eqs+1);
	}
	
      }
      
      /* Next entry */
      stack->size++;
      
      /* Allocate more stack if needed */
      if ( stack->size > stack->nmax ) 
      {
	stack->nmax = 2 * stack->size;
	stack->data = (char **) realloc(stack->data,
					stack->nmax * sizeof(char *));
	if (stack->data == NULL)
	{
	  fprintf(stderr, "ERROR: not enough memory\n");
	  return EXIT_FAILURE;
	}
      }
      
      /* Allocate the string copy itself */
      if (region ) 
	stack->data[ stack->size - 1 ] = stk_copy_string( region );
      else
	stack->data[ stack->size - 1 ] = stk_copy_string( "" );
      
    } /* end for yy */
  } /* end for xx */
  
  if (region){ free(region); }
  return EXIT_SUCCESS;
}

/*********************************************************************
 * Internal:
 * Appends the provided 'entry' item to the given stack, optionally 
 * prepending the <path to stack file> on the entry (??).
 *
 *    Stack stack     - U: Stack to be updated with new content
 *    char* entry     - I: Item to process and add to stack
 *    int   prepend   - I: flag control for prepending path on items being appended.
 *
 * entry - supported formats
 *   Stack File:
 *     @{1}\+{0:1}\-{0:1}\s*<path>{0:1}<file>{0:1}<filter>{0:1}
 *         <path>   = absolute or relative path to <file>
 *                      may include pattern matching syntax
 *         <file>   = name of stack file
 *                      may include pattern matching syntax
 *         <filter> = DM filter syntax '[...]'
 *                    ?? filter syntax '<...' ??
 *         '+'      = queues recursive addtion of stack content
 *         '-'      = disables prepend of path to stack records
 *     NOTE: If a DM filter is provided on a '@stack' type entry, that filter
 *           is propogated to the stack contents when appended to the stack.
 *     NOTE: Blank records within the stack file are ignored
 *     NOTE: Records starting with '#' are ignored as comment records.
 *     NOTE: Records starting with '!' will override the prepend setting for
 *            the stack and disable prepend of path for that record only
 *            (or children of record with pattern matching syntax).
 *
 *   Data File or other text entry (eg: region):
 *     <path>{0:1}<file>{0:1}<filter>{0:1}
 *         <path>   = absolute or relative path to <file>
 *                      may include pattern matching syntax
 *         <file>   = name of stack file
 *                      may include pattern matching syntax
 *         <filter> = DM filter syntax '[...]'
 * 
 * A Null argument returns a FAILURE status.
 *********************************************************************/
static int stk_append_item(Stack stack, char *entry, int prepend )
{
  FILE* fileList    = NULL;
  char* prepath     = NULL;
  char* pathedEntry = NULL;
  char* filteredEntry = NULL;
  int   i;
  char  fileEntry[MAX_ITEM_SIZE] = "";
  char* stackfile = NULL;
  char* filter;
  char* ptr;
  int status = EXIT_SUCCESS;
  
  logical recurse = FALSE;
  logical override_prepend = FALSE;
  
  size_t nGlobbed;
  glob_t* globbed;
  
  if ( NULL == stack ) return EXIT_FAILURE;
  if ( NULL == entry ) return EXIT_FAILURE;

  /* allocate glob memory after argument screen. */
  globbed=(glob_t*)calloc(1,sizeof(glob_t));
  
  if (*entry == '@') 
  {
    entry++;

    if ( *entry == '+' ) 
    {
      recurse=1;
      entry++;
    }
    
    if ( *entry == '-' ) 
    {
      prepend=0;
      entry++;
    }

    /* Skip whitespace between @ and filename */
    while( isspace( *entry ) )
      entry++;

    /* Check for and separate any filter expression */
    stackfile = stk_copy_string( entry );
    ptr = strchr( stackfile, '[' );
    
    if ( ptr == NULL )
      ptr = strchr( stackfile, '<' );
    
    if ( ptr == NULL ) 
    {
      filter = stk_copy_string( "" );
    }
    else 
    {
      filter = stk_copy_string( ptr );   /* Copy [.... to filter */
      *ptr = '\0';                       /* Terminate stackfile prior to [ char */
    }
    
    if ((fileList = fopen((stackfile ), "r")) == NULL)
    {
      fprintf(stderr, "# stklib : ERROR: can not open stack '%s'\n", stackfile );

      free(stackfile);
      free(filter);
      free(globbed);
      return EXIT_FAILURE;
    }
    else 
    { /* process stackfile content, add records to stack */
      
      /* setup prepath, removing last entry in / separated path */
      prepath = stk_copy_string( stackfile );
      i = strlen( prepath );
      while( i>= 0 && prepath[i]!='/' ) i--; 
      prepath[i+1] = '\0';

      char* tline = NULL;      /* A new string builder variable to hold an incomplete string */
      
      /* read in each line from the file */
      /*   - fgets reads at most 1 less than size characters. ('-1' is not necessary) */
      while ( fgets(fileEntry, MAX_ITEM_SIZE-1, fileList) != NULL )
      {
        /* build full record from fgets chunks */
        if (!tline)
        {
          /* starting new record */
          tline = (char*)malloc((strlen(fileEntry)+1) * sizeof(char));
          strcpy(tline, fileEntry);
        }
        else
        {
          /* add chunk to tline */
          tline = realloc(tline, (strlen(tline) + strlen(fileEntry) + 1) * sizeof(char));
          tline = strcat(tline, fileEntry);
        }

        /* Check if we need to keep building the record */
        if(tline[strlen(tline) - 1] == '\n')
        {
          /* Line in file ends with newline...*/
          /* - Trim the newline character and any potential whitespace leading to it */
          /* - check for a continuation character too */
          stk_trim(tline);
          if( (strlen(tline) > 0) && (tline[strlen(tline) - 1] == '\\') )
          {
            /* Trim the continuation character from the record. */
            tline[strlen(tline) - 1] = '\0';
            stk_trim(tline); /* Trim any potential whitespace leading up to where the continuation character was */
            continue; /* Go back to the start of the loop to get the next chunk to keep building the record */
          }
        }
        else
        {
          /* Check to make sure we didn't reach the end of the file reading this chunk */
          if(!feof(fileList))
	  {
            /* Go back and get the next chunk if there's still more to be read */
            continue;
	  }
        }
  
        /* tline now holds a complete record from the stackfile.. ready to process */
        ptr = tline;

        /* skip over leading whitespace; trim trailing newline and whitespace */
        while( *ptr == ' ' ) ptr++;
        stk_trim( ptr );

        /* Ignore comments and blank lines */
        if ((*ptr == '#' ) || strlen(ptr) == 0 )   
        {
          free(tline);
          tline = NULL;
          continue;
        }

        /* if record is itself another stack, and recursive is enabled,  */
        /* add its CONTENTS to the stack rather than the record.         */
        if ((*ptr == '@') && ( recurse ))
        {
          char* nextFile = NULL;
          filteredEntry = stk_cat_string( ++ptr, filter );
          if ( *ptr == '/' || !prepend )
            pathedEntry = stk_cat_string( "", filteredEntry );  /* don't append absolute paths */
          else
            pathedEntry = stk_cat_string( prepath, filteredEntry );

          nextFile = stk_cat_string( "@", pathedEntry );
          
          /* free allocated memory here in case the recursive load fails */
          free(filteredEntry);
          free(pathedEntry);
          free(tline);
          tline = NULL;
          
          /* load the substack contents */
          stk_append_item( stack, nextFile, prepend );

          /* can release this now. */
          free(nextFile);
          
          /* go back and build next record */
          continue;
          
        } /* end recursive expand */

        /* check special char to force no path prepend */
        if ( *ptr == '!' )
        {
          override_prepend = TRUE;
          ptr++;
        }

	/* prepend path to stack */
	/*   - do not prepend to absolute paths or if directed not to */
        if ( *ptr == '/' || !prepend || override_prepend )  
        {
          pathedEntry = stk_cat_string( "", ptr );
        }
        else
        {
          pathedEntry = stk_cat_string( prepath, ptr );
        }

        /* use glob to resolve pattern match records.. add each matched item to the stack */
        /* note: if there is no match glob returns the input                              */
        glob( pathedEntry, GLOB_NOCHECK, NULL, globbed );
        for ( nGlobbed = 0; nGlobbed < globbed->gl_pathc; nGlobbed++ )
        {
          /* add any "filter" content pulled from stackfile onto each stack record */
          filteredEntry = stk_cat_string( globbed->gl_pathv[nGlobbed], filter );

          /* add to stack */
          status = stk_append_entry( stack, filteredEntry );

          /* cleanup allocated memory */
          if (filteredEntry)
            free(filteredEntry);
  
          filteredEntry = NULL;
        }
        if (pathedEntry)
          free(pathedEntry);
        pathedEntry = NULL;

        /* restore prepend override flag to default */
        override_prepend = FALSE;
        
        /* we are done with the record now.. clear it for the next iteration */
        free(tline);
        tline = NULL;
        
        globfree( globbed );
        if ( status != EXIT_SUCCESS )
          return EXIT_FAILURE;

      } /* end while( fgets() ) */

      free(prepath);
      fclose(fileList);
    }
    status = EXIT_SUCCESS;
    free( stackfile );
    free( filter );
  }
  else
  { /* not an @ item */
    
    stackfile = stk_copy_string( entry );
    ptr = strchr( stackfile, '[' );
    
    if ( ptr == NULL ) {
      filter = stk_copy_string( "" );
    }
    else 
    {
      filter = stk_copy_string( ptr );   /* Copy [.... to filter */
      *ptr = '\0';                       /* Terminate stackfile prior to [ char */
    }
    
    glob( stackfile, GLOB_NOCHECK, NULL, globbed );
    for ( nGlobbed=0; nGlobbed< globbed->gl_pathc; nGlobbed++)
    {
      if ( ( stackfile[strlen(stackfile)-1] == '/'  ) &&
           ( strncmp( stackfile, globbed->gl_pathv[nGlobbed], strlen(stackfile)-1 ) == 0 ) )
      {
        // #14088 (SL-4): 11/10/2016
        //   if stackfile is a directory which does not exist, retain the "/".
        //   glob() removes it in this case, but does not if directory exists.
        //   So if the result == input except for ending '/'.. set using input
        filteredEntry=stk_cat_string( stackfile, filter );
      }
      else
      {
        filteredEntry=stk_cat_string( globbed->gl_pathv[nGlobbed], filter );
      }
      status = stk_append_entry( stack, filteredEntry);
      free( filteredEntry );
    }
    
    globfree( globbed );
    
    free( stackfile );
    free( filter );
  }
  
  free(globbed);
  return status;
}

int mock_stk_append_item(Stack stack, char *entry, int prepend)
{
  int status;
  status = stk_append_item( stack, entry, prepend );
  return status;
}

/*********************************************************************
 * Internal:
 * This function expands the given entry and appends the results to the 
 * given stack.  Null argument returns FAILURE.
 *********************************************************************/
static int stk_append_prepend(Stack stack, char *entry, int prepend )
{
  int  status = EXIT_SUCCESS;
  char NextItem[MAX_ITEM_SIZE];

  if ( NULL == stack ) return EXIT_FAILURE;
  if ( NULL == entry ) return EXIT_FAILURE;

  if ( strlen(entry) == 0 )
    status = stk_append_entry( stack, "");
  else
  {
    while ( stk_next_list_item( &entry, NextItem, MAX_ITEM_SIZE-1 ) ) 
    {
      status = stk_append_item( stack, NextItem, prepend );
    }
  }

  return( status );
}

/*********************************************************************
 * Internal:
 * Build a file stack by reading a list of files either 
 * in a string or in a file (specified by "@").
 *********************************************************************/
static Stack stk_build_prepend(char *inputList, /* text description of desired stack */ 
			       int prepend      /* Prepend paths by default */ )
{
  Stack stack;
  char  NextItem[MAX_ITEM_SIZE];

  stack = stk_alloc( STK_INIT_SIZE );
  if ( stack == NULL ) {
    return( stack );
  }
  
  if (inputList == NULL) {
    return (stack);
  }
  
  /* Get items from the input list, incrementing the string pointer */
  while ( stk_next_list_item( &inputList, NextItem, MAX_ITEM_SIZE-1 ) ) 
  {
    short retval;
    retval = stk_append_item( stack, NextItem, prepend );
    if (retval!=0) 
    {
      free(stack);
      return(NULL);
    }
  }
  if ( stack->size == 0 ) 
  {
    /*
     * if we have a null string, make it a one item null string stack
     * this is to accomodate situations where a null string may be valid
     * input as a single entry for a filename (i.e. enter null for a
     * program-determined filename
     */
    stk_append_entry( stack, "" );
  }
  
  /* reset to first entry */
  stk_rewind(stack);
  
  return (stack);
}

/*********************************************************************
 * Internal:
 *  Allocate a copy of two concatenated strings
 *********************************************************************/
static char* stk_cat_string( char* a, char* b )
{
  char* value;
  value = (char *) malloc((strlen(a) + strlen(b) + 1) * sizeof(char));
  if (value == NULL) 
  {
    fprintf( stderr, "Out of memory\n" );
  }
  else 
  {
    strcpy( value, a );
    if ( *b != '\0' )
      strcat( value, b );
  }
  return( value );
}

/*********************************************************************
 * Internal:
 * routine to allocate a copy of a string
 *********************************************************************/
static char* stk_copy_string( char* ptr )
{
  char* value;
  
  value = (char *) malloc((strlen(ptr) + 1) * sizeof(char));
  if (value == NULL) 
  {
    fprintf( stderr, "Out of memory\n" );
  }
  else 
  {
    strcpy( value, ptr );
  }
  return( value );
}

/*********************************************************************
 * Internal:
 *  Extract the next item from the string
 *  pptr is the address of the string pointer (i/o)
 *  item is a preallocated string of size maxlen to write the result to (o)
 *
 *  The list is a space or semicolon delimited list, but
 *  delimiters inside [] pairs are ignored.
 *********************************************************************/
static logical stk_next_list_item( char** pptr, char* itembuf, const long maxlen )
{
  char* ptr;
  char* item;
  long numBrackets;
  long numPara, numBrace;
  logical loop;
  logical result;
  logical not_sep;
  char last = '\0';

  ptr = *pptr;  /* Sacrifice efficiency for clarity */
  item = itembuf; /* So we can check itembuf at end */
  while( (isspace( *ptr ) || *ptr == STK_SEP || *ptr == STK_SEP_ALT))
    ptr++;
  
  numBrackets = 0;
  numPara = 0;
  numBrace=0;
  loop = TRUE;
  
  while ( loop ) 
  {
    loop = ( *ptr != '\0' ) && ( *ptr != '\n' );
    
    if ( loop && ( numBrackets == 0 ) && (numPara==0) && (numBrace == 0))
    {
      not_sep = ( *ptr != ' ' ) && ( *ptr != STK_SEP ) && (*ptr != STK_SEP_ALT);
      /* is_esc = (last == STK_ESC); */
      loop = ( not_sep /* || is_esc  */ ) && (item - itembuf < maxlen);
    }
    
    last = *ptr;
    
    if ( loop ) 
    {
      if ( *ptr == '[' ) numBrackets++;
      if ( *ptr == ']' ) numBrackets--;
      
      if ( *ptr == '(' ) numPara++;
      if ( *ptr == ')' ) numPara--;
      
      if ( *ptr == '{' ) numBrace++;
      if ( *ptr == '}' ) numBrace--;
      
      if ( ( !numBrackets ) && (!numPara) && (!numBrace) )
      {
	if ( *ptr == '<' ) numBrackets++;
	if ( *ptr == '>' ) numBrackets--;
      }
      
      if (*ptr == STK_ESC) 
      {
	last = *ptr;
	ptr++;
      }
      
      if ((( *ptr != ' ' ) && ( *ptr != STK_SEP ) && (*ptr != STK_SEP_ALT)) ||
	  ((last == STK_ESC) && (( *ptr == ' ' ) || 
				 ( *ptr == STK_SEP ) ||
				 (*ptr == STK_SEP_ALT))) ) 
      {
	last = *ptr;
	*item++ = *ptr++;
	
      } 
      else
      {
	*item++ = *ptr++;
      }
    }
  }
  
  *item = '\0';  /* Terminate output string */
  
  /* Skip forward to next entry  for tidiness */
  while( isspace( *ptr ) || *ptr == STK_SEP  || *ptr == STK_SEP_ALT )
    ptr++;
  
  *pptr = ptr;   /* Write back the new pointer position */
  result = ( *itembuf != '\0' );  /* True if we found anything */
  
  return( result );
}

/*********************************************************************
 * Internal:
 *  Remove trailing spaces from a string 
 *********************************************************************/
static void stk_trim( char* opt )
{
  int k;
  k = strlen( opt );
  while( k > 0 && ( opt[k-1] == ' ' || opt[k-1] == '\n' )) 
  { 
    k--;
  }
  opt[k] = '\0';
}

/*********************************************************************
 * Internal:
 *  Display contents of stack
 *********************************************************************/
void stk_disp( Stack stack )
{
  long ii;

  if ( stack == NULL ) 
  {
    printf( "Null stack\n" );
    return;
  }
  printf( "------\n" );
  printf( "Stack position: %4d\n", stk_current( stack ) );
  printf( "Stack size:     %4d\n", stk_count( stack ) );
  printf( "Stack allocated:%4ld\n", stack->nmax );
  printf( "Stack entries:\n" );
  for ( ii = 1; ii <= stack->size; ii++ ) 
  {
    printf( "%4ld :%s:\n", ii, stack->data[ii-1] );
  }
  printf( "------\n" );
}

/*
 * ####################################################
 * Here begin the external interface routines
 * ####################################################
 */

/*********************************************************************
 * Function: stk_append()
 * 
 * Description:
 *   Append item(s) to the end of the stack
 * 
 * Returns:
 *   0 = success
 *   1 = failure
 * 
 *********************************************************************/
int stk_append(Stack stack, char *entry )
{
  int status;
  status = stk_append_prepend( stack, entry, 1 );
  return status;
}

int stk_append_gen(Stack stack, char *entry )
{
  int status;
  status = stk_append_prepend( stack, entry, 0 );
  return status;
}

/*********************************************************************
 * Function: stk_build()
 * 
 * Description:
 *   Generate stack from provided file specification.
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack input
 * 
 *********************************************************************/
Stack stk_build(char *inputList /* text description of desired stack */ )
{
 int prepend = 1;
 return stk_build_prepend( inputList, prepend );
}

Stack stk_build_gen(char *inputList /* text description of desired stack */ )
{
 int prepend = 0;
 return stk_build_prepend( inputList, prepend );
}

/*********************************************************************
 * Function: stk_change_current()
 * 
 * Description:
 *   Replace value of current stack element with provided value.
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack input
 *   -1  -  Invalid (NULL) value input
 * 
 *********************************************************************/
int stk_change_current( Stack stack, char *value )
{
  int retval = 1;

  if ( stack )
    retval = stk_change_num( stack, value , stack->current );

  return( retval );
}

/*********************************************************************
 * Function: stk_change_num()
 * 
 * Description:
 *   Replace specified stack element value with provided value.
 *   If specified element < 0, the first element is replaced.
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack input
 *   -1  -  Invalid (NULL) value input
 *   -1  -  Invalid input; index out of range (larger than stack) 
 * 
 *********************************************************************/
int stk_change_num( Stack stack, char *value, int num )
{
  /* check null pointers */
  if ( !stack )
    return( 1 );

  if ( !value )
    return( -1 );

  /* check index against range */
  if ( num > stack->size )
    return( -1 );

  num--;  /* Stacks count 1 to N */

  if (num < 0) num=0;
  
  if ( stack->data[num] ) free( stack->data[num] );
  
  stack->data[num] = (char*)calloc( strlen(value)+1, sizeof(char));

  /* Remove leading white space, consistent with append  */
  while ( (*value == ' ') || (*value == '\t') || (*value == '\n') ) value++;
  
  strcpy( stack->data[num], value );
  
  return(0);
}

/*********************************************************************
 * Function: stk_close()
 * 
 * Description:
 *   Close a stack.
 *   This function closes a stack and frees memory allocated for 
 *   its content. 
 *
 *   Unlike the old version, this one doesn't mind if you close
 *   a NULL stack.
 * 
 * Returns:
 *   EXIT_SUCCESS  - successful.
 * 
 *********************************************************************/
int stk_close( Stack stack )
{
  long ii;

  if (stack != NULL) 
  {
    if ( stack->data != NULL ) {
      for (ii = 0; ii < stack->size; ii++) {
	if (stack->data[ii] != NULL)
	  free(stack->data[ii]);
      }
      free(stack->data);
    }
    free(stack);
  }
  
  return EXIT_SUCCESS;
}

/*********************************************************************
 * Function: stk_count()
 * 
 * Description:
 *   Count and return the number of entries in a stack.
 *   A NULL stack returns 0 (no entries).
 * 
 * Returns:
 *    #  -  integer number of entries.
 * 
 *********************************************************************/
int stk_count( Stack stack )
{
  if ( NULL == stack ) return(0);
  return ((int)stack->size);
}

/*********************************************************************
 * Function: stk_current()
 * 
 * Description:
 *   Give the current position in the stack. 
 *   This function returns an integer indicating the current element
 *   in the stack.  This is the element last read with the
 *   stk_read_num command.  It returns 0 for a fully rewound stack.
 *   If the end of the stack has been reached, it returns the number
 *   of the last entry, no matter how many stk_read_next operations
 *   have been performed since.
 * 
 * Returns:
 *    #  -  current stack element.
 * 
 *********************************************************************/
int stk_current( Stack stack )
{
  int result = 0;
  if ( stack )
  {
    if (stack->current <= stack->size) 
      result = (int)stack->current; 
    else 
      result = (int)stack->size;
  }
  return(result);
}

/*********************************************************************
 * Function: stk_delete_current()
 * 
 * Description:
 *   Remove current stack element.
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack input
 * 
 *********************************************************************/
int stk_delete_current( Stack stack )
{
  if ( !stack )
    return(1);

  return( stk_delete_num( stack, stack->current ));
}

/*********************************************************************
 * Function: stk_delete_num()
 * 
 * Description:
 *   Remove specified stack element.
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack input
 *   -1  -  Invalid input; index out of range (larger than stack) 
 *   -1  -  Invalid input; index out of range (< 0) 
 * 
 *********************************************************************/ 
int stk_delete_num( Stack stack, int num )
{
  if ( !stack )
    return(1);

  int ii;

  if ( num > stack->size ) return(-1);
  if ( num < 0 ) return(-1);
  if ( num == 0 ) num = 1;
  
  for ( ii=num; ii<=stack->size-1; ii++)
  {
    /* make sure there is enough space in target record */
    stack->data[ii-1] = (char*)realloc( stack->data[ii-1], 
					  (strlen(stack->data[ii])+1)*sizeof(char) );

    /* transfer value */
    strcpy( stack->data[ii-1], stack->data[ii] );
  }

  /* release last record */
  free( stack->data[stack->size-1] );
  stack->data[stack->size-1] = NULL;

  stack->size -= 1;

  return(0);
}

/*********************************************************************
 * Function: stk_expand_n()
 * 
 * Description:
 *   Given an integer 'n' and a file string containing a '#', 
 *   this function creates a stack of files from file string 
 *   with the '#' replaced by 1,2,...,n padded by
 *   leading zeroes.
 *   For example, (foo#.fits,11) would be expanded to
 *     foo01.fits, foo02.fits,...,foo11.fits.
 *   Omitting the '#' results in the stack foo.fits being returned.
 *   If n is outside the range 1:STK_INIT_SIZE, n is reset to 1.
 * 
 * Returns:
 *    stack  -  resulting stack
 * 
 *********************************************************************/
Stack stk_expand_n(  char* in_string, long int_suf_num )
{
  Stack stack;                 /* stack to be created & returned */
  char *entry = NULL;          /* names of stack files generated */
  char *pound;                 /* location of '#' in in-string */
  char *prefix;                /* string for chars prior to '#' */ 
  char *suffix;                /* string for chars after '#' */ 
  char *ii_suf_num;            /* chars for each number in file name */
  long i_strlen = 0, ii = 0, jj = 0;
  int tot_num_of_digits = 0;   /* total number of digits in filename */
  int num_of_lead_zeros = 0;   /* number of leading padding zeros */ 
  int num_of_digits_tra = 0;   /* number of trailing digits */ 
  int prepend = 1;  

  stack = stk_alloc( int_suf_num );
  if ( stack == NULL ) {
    fprintf(stderr, "ERROR: not enough memory\n");
    return( stack );
  }
  if ( in_string == NULL ) {
     return ( stack );
  }
  
  /* verify that input integer is within bounds */
  if ( int_suf_num < 1 )  { 
    fprintf( stderr,"Number of input stack items reset to 1 \n");
    int_suf_num = 1;
  }  

  /* locate '#' in input string.  If NULL, return input string  */
  pound = strchr( in_string,'#');
  if ( NULL == pound ) { /* no '#' found , so return in_string as stack */
    fprintf( stderr," No # given, so setting stack to input string\n");
    stack = stk_build( in_string );
    return (stack);
  }

  /* get total number of digits in and length for stack file names */
  tot_num_of_digits = ((log10( (double) int_suf_num )) + 1 ); 
  i_strlen = strlen( in_string) + tot_num_of_digits;

  /* allocate memory for strings */
  prefix = (char*) malloc( (i_strlen + 1) * sizeof( char ));
  suffix = (char*) malloc( (i_strlen + 1) * sizeof( char ));
  ii_suf_num = (char*) malloc( (i_strlen + 1) * sizeof( char ));
  entry = (char*) malloc( (i_strlen + 1) * sizeof( char ));  
  if (( NULL==prefix)||( NULL==suffix)||( NULL == ii_suf_num ||( NULL == entry ))){
    fprintf(stderr, "ERROR: not enough memory\n");
    return (stack);
  }
  
  /* create prefix by copying up to '#' and null-terminating */
  strncpy( prefix, in_string, (size_t)(pound - in_string) );
  prefix[ (size_t)(pound - in_string) ] = '\0';

  /* create suffix by copying past '#' */
  strcpy( suffix, pound + 1 );
  
  /* create each entry and append to stack */
  for ( ii = 1 ; ii <= int_suf_num; ii++ ) {
    num_of_digits_tra = log10( (double) ii ) + 1;
    
    /* pad entry with appropriate number of leading 0's */
    num_of_lead_zeros = tot_num_of_digits - num_of_digits_tra;
    strcpy( entry, prefix );
    for ( jj = 0; jj < num_of_lead_zeros; jj++ ) { 
      strcat( entry, "0");
    }
    sprintf( ii_suf_num,"%ld", ii );
    strcat( entry, ii_suf_num );
    strcat( entry, suffix );
    stk_append_item( stack, entry, prepend ); 
  }
  
  if ( prefix)
    free( prefix );
  if ( suffix)
    free( suffix );
  if( ii_suf_num ) 
    free( ii_suf_num ) ;
  if( entry )
    free( entry );
  
  /* reset to first entry */
  stk_rewind(stack);
  
  return (stack);
}

/*********************************************************************
 * Function: stk_read_free()
 * 
 * Description:
 *   Free memory previously allocated for returned entry value.
 * 
 * Returns:
 *   None
 * 
 *********************************************************************/
void stk_read_free( char *name )
{
  if (name)
    free(name);

  return;
}

/*********************************************************************
 * Function: stk_read_next()
 * 
 * Description:
 *   Read the next stack item. 
 *   This function reads the next stack item.  The first call to 
 *   this function will read the first item in the stack.  Subsequent
 *   calls step through the stack until the end, where it will 
 *   return a NULL.  
 *
 *   Calls to stk_read_num will not affect the sequence of filenames
 *   this function returns.
 * 
 * Returns:
 *   entry = copy of entry.
 * 
 *********************************************************************/
char* stk_read_next(Stack stack /* stack to read next entry from */ )
{
  if ( NULL == stack ) return(NULL);
  
  if (stack->current >= stack->size) {
    stack->current = stack->size;
    return (NULL);
  }
  else 
  {
    stack->current++;
    return( stk_copy_string( stack->data[stack->current-1] ) );
  }
}

/*********************************************************************
 * Function: stk_read_num()
 * 
 * Description:
 *   Read stack entry by number.
 *   This function reads an entry from a stack designated by number.
 *   The number 1 references the first entry in the stack.  An 
 *   attempt to read an entry from beyond the end of the stack
 *   will return a NULL.
 *
 *   Use of this function does NOT affect which item 
 *   stk_read_next will get in any subsequent calls.
 * 
 * Returns:
 *   entry = copy of entry.
 * 
 *********************************************************************/
char* stk_read_num(Stack stack, int entry)
{
  if ( NULL == stack ) return(NULL);
  
  if (entry <= 0) {
    fprintf(stderr, "# stklib : ERROR: stack entries start from 1\n" );
    return NULL;
  } else if (entry > stack->size ) {
    return NULL;
  } else {
    return( stk_copy_string( stack->data[entry-1] ) );
  }
}

/*********************************************************************
 * Function: stk_rewind()
 * 
 * Description:
 *   Reset a stack to the first entry.
 *   This function resets a stack so that subsequent calls to
 *   stk_read_next will start at the first entry again.
 * 
 * Returns:
 *   none
 * 
 *********************************************************************/
void stk_rewind( Stack stack )
{
  if ( NULL == stack ) return;
  stack->current = 0;
}

/*********************************************************************
 * Function: stk_set_current()
 * 
 * Description:
 *   Set stack to specified entry.
 *   Entry > stack size; sets current to last entry
 *   Entry < 1; sets current to first entry
 * 
 * Returns:
 *    0  -  Success
 *    1  -  Invalid (NULL) stack
 *    1  -  Invalid entry; value > stack size.
 *   -1  -  Invalid entry; value < 1.
 * 
 *********************************************************************/
int stk_set_current( Stack stack, int num )
{
  int retval=0;

  if ( !stack )
    return( 1 );

  if ( num > stack->size ) { num = stack->size; retval = 1; }
  if ( ( num < 1 ) && ( stack->size > 1 ) ) { num = 1 ; retval = -1; }

  stack->current = num;

  return( retval );
}
