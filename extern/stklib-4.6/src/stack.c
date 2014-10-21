/*_C_INSERT_SAO_COPYRIGHT_HERE_(1996,1998,2000,2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

/*
 * H**************************************************************************
 * 
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
 * 
 * NOTES:
 * 
 * <None>
 * 
 * REVISION HISTORY:
 * 
 * Ref. No.        Date
 * 
 * --------        ----
 * 
 * 1.12             05/07/96
 * 2.0              1998 Jul 14 JCM
 *                  03/03/00 
 H*************************************************************************
 */


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
static Stack stk_alloc( long size );
static char* stk_copy_string( char* ptr );
static char* stk_cat_string( char* a, char* b );
static logical stk_next_list_item( char** pptr, char* itembuf, long maxlen );
static int stk_append_entry( Stack stack, char* item );
static int stk_append_item( Stack stack, char* item,int prepend );
static int stk_append_prepend(Stack stack, char *entry, int prepend );
static void dmstk_combine_filters( char *, char *);
static Stack
stk_build_prepend(char *inputList, /* text description of desired stack */ 
                  int prepend  /* Prepend paths by default */ );

void dmstk_combine_filters( char *inString, char *filter )
{
  char region[MAX_ITEM_SIZE];

  char *tok;
  int bin, col, ext, filt;

  bin=0; col=0; ext=0; filt=0;


  tok = strtok( inString, "[]");
  while ( tok != NULL)
    {

      bin += ( strstr( tok, "bin " ) || strstr(tok, "BIN ") );
      col += ( strstr( tok, "cols ") || strstr(tok, "COLS ") ||
	       strstr( tok, "columns ") || strstr(tok, "COLUMNS ") );
      
      filt += ( !bin & (strpbrk( tok, "<=>:" ) != NULL));

      
      if ( !bin && !col && !filt ) ext+=1;



      if ( filt && !bin  )
	{
	  sprintf( region, "%s[%s,%s]", region, tok, filter );
	}
      else if ( !filt && bin )
	{
	  sprintf( region, "%s[%s][%s]", region, filter, tok );
	  filt =1;
	}
      else if ( !filt && col )
	{
	  sprintf( region, "%s[%s][%s]", region, filter, tok );
	  filt =1;
	}
      else if ( !filt && ( ext > 2 ) )
	{
	  sprintf( region, "%s[%s][%s]", region, filter, tok );      
	  filt = 1;
	}
      else if ( tok == inString )
	{
	  strcpy(region, tok );
	}
      else if ( ext <= 2 )
	{
	  sprintf( region, "%s[%s]", region, tok );
	}
      else
	{
	  sprintf( region, "%s[%s]", region, filter );
	}

      tok = strtok( NULL, "[]");
    }

  if ( !filt )
	  sprintf( region, "%s[%s]", region, filter );    

  strcpy( inString, region );
      
  return;

}


/* Remove trailing spaces from a string*/
static void stk_trim( char* opt )
{
 int k;
 k = strlen( opt );
 while( k > 0 && ( opt[k-1] == ' ' || opt[k-1] == '\n' )) { 
  k--;
 }
 opt[k] = '\0';
}

/*
 *  Extract the next item from the string
 *  pptr is the address of the string pointer (i/o)
 *  item is a preallocated string of size maxlen to write the result to (o)
 *
 *  The list is a space or semicolon delimited list, but
 *  delimiters inside [] pairs are ignored.
 */
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
 while( (isspace( *ptr ) || *ptr == STK_SEP || *ptr == STK_SEP_ALT)) ptr++;

 numBrackets = 0;
 numPara = 0;
 numBrace=0;
 loop = TRUE;

 while ( loop ) {
    loop = ( *ptr != '\0' ) && ( *ptr != '\n' );

    if ( loop && ( numBrackets == 0 ) && (numPara==0) && (numBrace == 0))
    {
        not_sep = ( *ptr != ' ' ) && ( *ptr != STK_SEP ) && (*ptr != STK_SEP_ALT);
	/* is_esc = (last == STK_ESC); */
	loop = ( not_sep /* || is_esc  */ ) && (item - itembuf < maxlen);
    }

    last = *ptr;

    if ( loop ) {

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

	if (*ptr == STK_ESC) {
	  last = *ptr;
	  ptr++;
	}

	if ((( *ptr != ' ' ) && ( *ptr != STK_SEP ) && 
		     (*ptr != STK_SEP_ALT)) ||
		    ((last == STK_ESC) && (( *ptr == ' ' ) || 
					   ( *ptr == STK_SEP ) ||
					   (*ptr == STK_SEP_ALT))) ) {

	  last = *ptr;
	  *item++ = *ptr++;

	} else {
	  *item++ = *ptr++;
	}

    }
  }

 *item = '\0';  /* Terminate output string */

 /* Skip forward to next entry  for tidiness */
 while( isspace( *ptr ) || *ptr == STK_SEP  || *ptr == STK_SEP_ALT ) ptr++;

 *pptr = ptr;   /* Write back the new pointer position */
 result = ( *itembuf != '\0' );  /* True if we found anything */

 return( result );
}

/*
 *  Internal routine to allocate a copy of a string
 */
static char* stk_copy_string( char* ptr )
{
 char* value;

 value = (char *) malloc((strlen(ptr) + 1) * sizeof(char));
 if (value == NULL) {
   fprintf( stderr, "Out of memory\n" );
 } else {
   strcpy( value, ptr );
 }
 return( value );
}


/*
 *  Internal routine to allocate a copy of two concatenated strings
 */
static char* stk_cat_string( char* a, char* b )
{
 char* value;
 value = (char *) malloc((strlen(a) + strlen(b) + 1) * sizeof(char));
 if (value == NULL) {
   fprintf( stderr, "Out of memory\n" );
 } else {
   strcpy( value, a );
   if ( *b != '\0' ) strcat( value, b );
 }
 return( value );
}

/*
 *  This internal routine allocates memory for a stack.
 */
static Stack stk_alloc( long size )
{
 Stack stack;
 
 stack = (Stack)malloc( sizeof( StackData ) );
 if ( stack ) {

  stack->data = (char**)malloc( size * sizeof( char* ) );
  if ( stack->data ) {
   stack->nmax = size;
   stack->current = 0;
   stack->size = 0;
  } else {
   free( stack );
   stack = NULL;
  }
 }
 if ( !stack ) {
  fprintf(stderr, "ERROR: not enough memory to allocate stack\n");
 }
 return( stack );
}


/*
 * Add a new entry to the stack, reallocating the stack if it needs
 * more memory.
 */

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
  if ( item && strlen(item) ) {
    region = (char*)calloc(strlen(item)+100,sizeof(char));
  } else {
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
	    } else {
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
	    } else {
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

	  if (region) strcpy( region, item );	  
	  
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
	      
	      /* dmstk_combine_filters( region, filter ); */
	      if ( region ) {
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
	  if ( stack->size > stack->nmax ) {
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

    }/* end for xx */
  
  if (region) free(region);
 return EXIT_SUCCESS;
}


static int
stk_append_item(Stack stack, char *entry, int prepend )
{
	FILE           *fileList    = NULL;
	char	       *prepath     = NULL;
	char	       *pathedEntry = NULL;
	char* filteredEntry = NULL;
	int             i;
	char            fileEntry[MAX_ITEM_SIZE] = "";
	char*  stackfile = NULL;
	char*  filter;
	char* ptr;
	int status = EXIT_SUCCESS;

	int recurse=0;

	size_t nGlobbed;
	glob_t *globbed=(glob_t*)calloc(1,sizeof(glob_t));

	if ( NULL == stack ) return EXIT_FAILURE;
	if ( NULL == entry ) return EXIT_FAILURE;

	if (*entry == '@') {
	  entry++;
	  
	  if ( *entry == '+' ) {
	    recurse=1;
	    entry++;
	  }
	  
	  if ( *entry == '-' ) {
	    prepend=0;
	    entry++;
	  }
	  

         while( isspace( *entry ) ) entry++;  /* Skip whitespace between @ and filename */
         stackfile = stk_copy_string( entry );
         ptr = strchr( stackfile, '[' );

	 if ( ptr == NULL )
	   ptr = strchr( stackfile, '<' );

	 if ( ptr == NULL ) {
	  filter = stk_copy_string( "" );
         } else {          
          filter = stk_copy_string( ptr );   /* Copy [.... to filter */
          *ptr = '\0';                       /* Terminate stackfile prior to [ char */
	 }
	 
	 if ((fileList = fopen((stackfile ), "r")) == NULL)
	   {
	     fprintf(stderr, "# stklib : ERROR: can not open stack '%s'\n",
		     stackfile );
	     return EXIT_FAILURE;
	   }
	 
	 
	 /* for each line, add process line read */
	 
	 else {
	   
	   /* setup prepath, removing last entry in / separated path */
	   char *lastline = NULL;

	   prepath = stk_copy_string( stackfile );
	   i = strlen( prepath );
           while( i>= 0 && prepath[i]!='/' ) i--; 
           prepath[i+1] = '\0';
		/* read in each line */

	   while (fgets(fileEntry, MAX_ITEM_SIZE-1, fileList) != NULL) {
	     
	     /* skip over whitespace */
	     ptr = fileEntry;
	     while( *ptr == ' ' ) ptr++;
	     
	     stk_trim( ptr );
	     
	     
	     if ( (*ptr == '@') && ( recurse ))
	       {
		 char nextFile[MAX_ITEM_SIZE];
		 strcpy( nextFile, "@" );
		 if ( prepend )
		   strcat( nextFile, prepath );
		 strcat( nextFile,++ptr);
		 strcat( nextFile, filter);
		 
		 stk_append_item(stack, nextFile, prepend );
		 continue;
	       } /* end recursive expand */
	     
	     if ((*ptr == '#' ) || strlen(ptr) == 0 )   /* Ingore comments and blank line */
	       {
		 continue;
	       }
	     
	     
	     glob( ptr, GLOB_NOCHECK, NULL, globbed );
	     for ( nGlobbed=0; nGlobbed< globbed->gl_pathc; nGlobbed++)
	       {
		 
		 if ( lastline )   /* This allows for a continuation character */
		   {
		     char *ts = stk_cat_string( lastline, globbed->gl_pathv[nGlobbed]);
		     filteredEntry=stk_cat_string( ts, filter );
		     if (ts) free(ts);
		     free(lastline);
		     lastline = NULL;
		     prepend = 0; /* If continuing, don't prepend */
		   }
		 else
		   {
		     filteredEntry=stk_cat_string( globbed->gl_pathv[nGlobbed], filter );
		   }


		 
		 if ( *ptr == '/' || !prepend )  /* don't append absolute paths */
		   {
		     pathedEntry = stk_cat_string( "", filteredEntry );
		   }
		 else if ( *filteredEntry == '!' )  /* special char to force no path pre-pend */
		   { 
		     pathedEntry = stk_cat_string( "", filteredEntry+1);
		   }
		 else
		   {
		     pathedEntry = stk_cat_string( prepath, filteredEntry );
		   }
		 
		 if (( pathedEntry[strlen(pathedEntry)-1] == '\\' ) || 
		     ( strlen(fileEntry)==(MAX_ITEM_SIZE-2)) )/* Allow for a 
								      continuation character */
		   {
		     lastline = (char*)calloc(strlen(pathedEntry)+1, sizeof(char));
		     strcpy( lastline, pathedEntry );
		     if ( lastline[strlen(lastline)-1] == '\\' ) {
		       /* Remove trailing \ */
		       lastline[strlen(lastline)-1] = '\0';
		     }
		   }
		 else
		   {
		     if (lastline) free(lastline);
		     lastline = NULL;
		   }

		 if ( !lastline )
		   status = stk_append_entry( stack, pathedEntry );

		 free(filteredEntry);
		 if (pathedEntry) free(pathedEntry);
		 pathedEntry=NULL;
	       }
	     
	     globfree( globbed );
	     if ( status != EXIT_SUCCESS ) return EXIT_FAILURE;
	     
  	   }
   	   free(prepath);
	   fclose(fileList);
         } 
         status = EXIT_SUCCESS;
	 free( stackfile );
	 free( filter );
	} else { /* not an @ item */
	  
	  stackfile = stk_copy_string( entry );
	  ptr = strchr( stackfile, '[' );
	  
	  if ( ptr == NULL ) {
	    filter = stk_copy_string( "" );
	  } else {          
	    filter = stk_copy_string( ptr );   /* Copy [.... to filter */
	    *ptr = '\0';                       /* Terminate stackfile prior to [ char */
	  }
	  

	  glob( stackfile, GLOB_NOCHECK, NULL, globbed );
	  for ( nGlobbed=0; nGlobbed< globbed->gl_pathc; nGlobbed++)
	    {
	      filteredEntry=stk_cat_string( globbed->gl_pathv[nGlobbed], filter );
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



/*
 *####################################################
 * Here begin the external interface routines
 *####################################################
 */

/*
 * Build a file stack by reading a list of files either 
 * in a string or in a file (specified by "@").
 */

static Stack
stk_build_prepend(char *inputList, /* text description of desired stack */ 
                  int prepend  /* Prepend paths by default */ )
{
	Stack           stack;
	char            NextItem[MAX_ITEM_SIZE];
	stack = stk_alloc( STK_INIT_SIZE );
	
	if ( stack == NULL ) {
	  return( stack );
        }
	
	if (inputList == NULL) {
	  return (stack);
	}

        /* Get items from the input list, incrementing the string pointer */

	while ( stk_next_list_item( &inputList, NextItem, MAX_ITEM_SIZE-1 ) ) {
	  short retval;
	  retval = stk_append_item( stack, NextItem, prepend );
	  if (retval!=0) {
	    free(stack);
	    return(NULL);
	  }
        }
	if ( stack->size == 0 ) {
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



Stack
stk_build(char *inputList /* text description of desired stack */ )
{
 int prepend = 1;
 return stk_build_prepend( inputList, prepend );
}

Stack
stk_build_gen(char *inputList /* text description of desired stack */ )
{
 int prepend = 0;
 return stk_build_prepend( inputList, prepend );
}



/*
 * count the number of entries in a stack. This function counts the total
 * number of entries in an IRAF stack. returns integer number of entries.
 */

int
stk_count(Stack stack /* IRAF stack to be counted */ )
{
	
  if ( NULL == stack ) return(0);
  return ((int)stack->size);
}


/*
 * Give the current position in the stack. This function returns an integer
 * indicating the current element in the stack.  This is the element last
 * read with the stk_read_num command.  It returns 0 for a fully rewound
 * stack.  If the end of the stack has been reached, it returns the number of
 * the last entry, no matter how many stk_read_next operations have been
 * performed since.
 */

int
stk_current(Stack stack /* IRAF stack to check */ )
{

	if ( NULL == stack ) return(0);

	if (stack->current <= stack->size) {
		return ((int)stack->current);
	} else {
		return ((int)stack->size);
	}
}






/*
 * Append descriptor. This function expands the given descriptor and appends
 * the results to the given stack.
 */
static int stk_append_prepend(Stack stack, char *entry, int prepend )
{
 int status =0;
	char            NextItem[MAX_ITEM_SIZE];

	if ( strlen(entry) == 0 )
	  status = stk_append_entry( stack, "");
	else
	  while ( stk_next_list_item( &entry, NextItem, MAX_ITEM_SIZE-1 ) ) {
	    status = stk_append_item( stack, NextItem, prepend );
	  }



 return( status );
}

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




/*
 * read stack entry by number. This function reads an entry from a
 * stack designated by number.  The number 1 references the first entry in
 * the stack.  An attempt to read an entry from beyond the end of the stack
 * will return a NULL.  Use of this function does not affect which item
 * stk_read_next will get in any subsequent calls.
 */

char *
stk_read_num(Stack stack, int entry)
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


/*
 * read the next IRAF stack item. This function reads the next IRAF stack
 * item.  The first call to this function will read the first item in the
 * stack.  Subsequent calls step through the stack until the end, when it
 * will return a NULL.  Calls to stk_read_num will not affect the sequence of
 * filenames this function returns. returns next filename in IRAF stack.
 */

char *
stk_read_next(Stack stack /* stack to read next entry from */ )
{
	if ( NULL == stack ) return(NULL);

	if (stack->current >= stack->size) {
		stack->current = stack->size;
		return (NULL);
	} else {
	        stack->current++;
	        return( stk_copy_string( stack->data[stack->current-1] ) );
	}
}


/*
 * reset a stack to the first entry. This function resets a stack
 * so that subsequent calls to stk_read_next will start at the first entry
 * again.
 */

void
stk_rewind(Stack stack /* stack to be rewound */ )
{
  if ( NULL == stack ) return;
  stack->current = 0;
}


/*
 * close a stack. This function closes a stack, freeing the
 * memory it uses for other tasks. returns EXIT_SUCCESS if successful.
 * Unlike the old version, this one doesn't mind if you close a NULL stack.
 */

int
stk_close(Stack stack /* stack to close */ )
{
	long             ii;


	if (stack != NULL) {
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

/*
 *  Display contents of stack
 */
void stk_disp( Stack stack )
{
 long i;

 if ( stack == NULL ) {
  printf( "Null stack\n" );
  return;
 }
 printf( "------\n" );
 printf( "Stack position: %4d\n", stk_current( stack ) );
 printf( "Stack size:     %4d\n", stk_count( stack ) );
 printf( "Stack allocated:%4ld\n", stack->nmax );
 printf( "Stack entries:\n" );
 for ( i = 1; i <= stack->size; i++ ) {
  printf( "%4ld :%s:\n", i, stack->data[i-1] );
 }
 printf( "------\n" );
}

static void fgets_trim( char* opt, int maxlen, FILE* file )
{
 int k;
 fgets( opt, maxlen, file );
 stk_trim( opt );
 k = strlen( opt );
 while( k > 0 && ( opt[k-1] == ' ' || opt[k-1] == '\n' )) k--;
 opt[k] = '\0';
}


/*
 * Given an integer i and a file string containing a '#', creates a stack
 * of files from file string with the '#' replaced by 1,2,...,i padded by
 * leading zeroes. For example, (foo#.fits,11) would be expanded to
 * foo01.fits, foo02.fits,...,foo11.fits.  Omitting the '#' results in a
 * the stack foo.fits being returned. If i is outside the range
 * (1,STK_INIT_SIZE), i is reset to 1.
 */
Stack 
stk_expand_n(  char* in_string, long int_suf_num )
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


/*
 *  Stack test routine
 */
void stk_test( void )
{
#define CSIZE 80
 char opt[CSIZE];
 char* item;
 int n;
 Stack stack= NULL;
 strcpy( opt, " " );
 while( *opt != 'q' ) {
  printf( "stk: " );
  fscanf( stdin, "%s", opt );
  if ( !strcmp( opt, " " ) ) {
   ;
  } else if ( !strcmp( opt, "build" ) ) {
   fgets_trim( opt, CSIZE, stdin );
   stack = stk_build( opt );
   printf( " \n" );
  } else if ( !strcmp( opt, "close" ) ) {
   stk_close( stack );
   stack = NULL;
   printf( "Close\n" );
  } else if ( !strcmp( opt, "echo" ) ) {
   fgets_trim( opt, CSIZE, stdin );
   printf( "%s\n", opt );
  } else if ( !strcmp( opt, "append" ) ) {
   fgets_trim( opt, CSIZE, stdin );
   stk_append( stack, opt );
   printf( "Rewind\n" );
  } else if ( !strcmp( opt, "disp" ) ) {
   stk_disp( stack );
  } else if ( !strcmp( opt, "rewind" ) ) {
   stk_rewind( stack );
   printf( "Rewind\n" );
  } else if ( !strcmp( opt, "read" ) ) {
   fgets_trim( opt, CSIZE, stdin );
   n = atoi( opt );
   item = stk_read_num( stack, n );   
   printf( "Item %4d = %s\n", n, item );
   free( item );
  } else if ( !strcmp( opt, "next" ) ) {
   item = stk_read_next( stack );
   n = stk_current( stack );
   printf( "Item %4d = %s\n", n, item );
   free( item );
  } else if ( !strcmp( opt, "delete" ) ) {
    fgets_trim( opt, CSIZE, stdin );
    n = atoi( opt );
    if ( n >= 1 )
      stk_delete_num( stack, n );
    else
      stk_delete_current(stack );
  } else if ( !strcmp( opt, "set" ) ) {
    fgets_trim( opt, CSIZE, stdin );
    n = atoi( opt );
    stk_set_current( stack, n );
  } else if ( !strcmp( opt, "change" )) {
    fgets_trim( opt, CSIZE, stdin );
    stk_change_current( stack, opt );
  } else if ( strcmp( opt, "q" ) ) {
    printf( "Unknown opt\n" );
  } 

 }
 printf( "\n" );
}



int stk_delete_current( Stack stack )
{
  return( stk_delete_num( stack, stack->current ));
}

 

int stk_delete_num( Stack stack, int num )
{

  int ii;

  if ( num > stack->size ) return(-1);
  if ( num < 0 ) return(-1);
  if ( num == 0 ) num = 1;
  
  for ( ii=num; ii<=stack->size-1; ii++)
    {
      stack->data[ii-1] = (char*)realloc(  stack->data[ii-1], 
					   (strlen( stack->data[ii])+1)*sizeof(char) );

      strcpy( stack->data[ii-1], stack->data[ii] );
    }

  free( stack->data[stack->size-1] );

  stack->data[stack->size-1] = NULL;


  stack->size -= 1;

  return(0);
}



int stk_change_current( Stack stack, char *value )
{
  return( stk_change_num( stack, value , stack->current ));
  
}



int stk_change_num( Stack stack, char *value, int num )
{
  
  
  if ( num > stack->size )
    return(-1);
  
  if ( !value )
    return(-1);
  
  num--;  /* Stacks count 1 to N */

  if (num < 0) num=0;
  
  if ( stack->data[num] ) free( stack->data[num] );
  
  stack->data[num] = (char*)calloc( strlen(value)+1, sizeof(char));

  /* Remove leading white space, consistent with append  */
  while ( (*value == ' ') || (*value == '\t') || (*value == '\n') ) value++;
  
  strcpy( stack->data[num], value );
  
  return(0);
}



int stk_set_current( Stack stack, int num )
{
  int retval=0;

  if ( num > stack->size ) { num = stack->size; retval = 1; }
  if ( ( num < 1 ) && ( stack->size > 1 ) ) { num = 1 ; retval = -1; }

  stack->current = num;

  return( retval );
}


void stk_read_free( char *name )
{
  if (name) free(name);
  return;
}

