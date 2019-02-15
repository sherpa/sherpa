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

#include "region_priv.h"
#include <string.h>
#include <ctype.h>

#define REGION_TRUNCATED_STR "...[truncated]"

/*
 * Module for string operations with regions. Including:
 * regReadAsciiRegion
 * regWriteAsciiRegion
 * regToStringRegion
 * regComposeRegion
 * regAllocComposeRegion
 * regPrintRegion
 */

int reg_read_line(FILE * fp, char *buf, long maxlen);
void reg_parse_line(char *buf, long *mode, char **stringptr, long *maxlen, 
        long *sys);
void reg_areg_line(FILE * out, regShape * shape, char *shapeName, long nr, long nangle, 
        char *text, int world);
char *reg_token_advance(char *ptr, char *colname, char sep);
void reg_strcat(char **ptr, long *maxlen, char sep, char *buf);
char *reg_tokens_advance(char *ptr, char *colname, char *seps);
void reg_areg_hdr(FILE * out, char *color);
void reg_compose_alloc_shape(regShape * shape, long cpt, char **pbuf,
			  char **pptr, long *size, int alloc);

regRegion* regReadAsciiRegion( char* filename, int verbose )
{
    regRegion *region = NULL;
    FILE *fp = NULL;

    char buf[SZ_LARGE+1];
    long maxlen = SZ_LARGE;
    char *ptr;
    long sys = RC_UNK;
    char *string = NULL;
    long mode = RF_UNK;
    char *smodes[] = { "UNK", "PHYSICAL", "LOGICAL", "WORLD" };
    char *modes[] = { "UNK", "SAOIMAGE", "SAOTNG", "PROS", "CXC", "DS9",
		      "DS9_V4" };
    fp = fopen(filename, "r");

    if (!fp) {
	    return NULL;
    }
    if (verbose >= 1) {
	    fprintf(stderr, "cxcregion: parsing region file %s\n", filename);
    }

    string = calloc(SZ_LARGE, sizeof(char *));
    while (reg_read_line(fp, buf, SZ_LARGE)) {
	    ptr = buf;
	
	    if ((ptr[0] == '#') && (mode != RF_UNK)) {
	        continue;
        }

	    while (*ptr == ' ') {
	        ptr++;
        }
	    if (*ptr) {
	        reg_parse_line(ptr, &mode, &string, &maxlen, &sys);
	    }
    }

    region = regParse(string);
    
    if (verbose >= 4) {
	    fprintf(stderr, "Ascii Region Parse String = %s\n", string);
	    fprintf(stderr, "SYS = %s FORMAT = %s\n", smodes[sys],
		        modes[mode]);
	    regPrintRegion(region);
    }
    free(string);
    return region;
}


/*
 * Writes the ascii region to the file, name. Or to stderr or stdout if those
 * are the specified name.
 *
 * TODO:
 *  I believe that names and nn refer to names for the individual shapes, though
 *  need to confirm that.
 */
int regWriteAsciiRegion( char* name, regRegion* region, char** names, long nn )
{
    int world = 0;
    char color[] = "blue";
    FILE *out;
    char shapeName[21];
    long maxlen = 20;
    long nradii, nangles;
    char text[256];
    regShape *Shape;
    long shapeNo;
    short int close = 0;

    if (!region) {
        return 0;
    }

    if (( strcasecmp( name, "stdout" ) == 0 ) ||
	    ( strcmp( name, "-" ) == 0 )) {
      out = stdout;
    }
    else if ( strcasecmp( name, "stderr") == 0 ) {
      out = stderr;
    }
    else {
        out = fopen(name, "w");
        close = 1;
    }
    if (!out) {
      return 0;
    }

    // File header
    reg_areg_hdr(out, color);
    Shape = region->shape;
    shapeNo = 0;
    while (Shape) {
        // Shapename
        if (Shape->include == regInclude) {
            snprintf(shapeName, maxlen, Shape->name);
        }
        else {
            snprintf(shapeName, maxlen, "-%s", Shape->name);
        }

        // Worth adding a field to the shape to get rid of these?
        nradii = reg_shape_radii(Shape);
        nangles = reg_shape_angles(Shape);

        // If the shape has a name/note add it
        if (names && (shapeNo <= nn)) {
	        sprintf(text, "text=\"%s\"", names[shapeNo]);
        }
        else {
            strcpy(text, " ");
        }

        // Print shape line to file
        reg_areg_line(out, Shape, shapeName, nradii, nangles, text, world);

        Shape = Shape->next;
        shapeNo++;
    }
    
    if (close) {
        fclose(out);
    }
    return 1;
}

/*
 * The same as regToStringRegion, but takes a fixed length array instead of
 * dynamically allocated memory for the string.
 */
void regComposeRegion( const regRegion* region, char* buf, const long maxlen )
{
    if (!region) {
	    snprintf(buf, maxlen, "Null region");
        return;
    }

    regShape *shape;
    char *ptr;
    long cpt = 0;
    strcpy(buf, "");
    long size = maxlen;

    ptr = buf;
    
    shape = region->shape;
    if (!shape) {
        snprintf(buf, maxlen, "Empty region");
        return;
    }
    
	// ptr will be set to NULL if we run out of space
    while (shape != NULL && ptr) {
	    reg_compose_alloc_shape(shape, cpt, &buf, &ptr, &size, 0);
	    cpt = shape->component;
	    shape = shape->next;
    }
} 


/*
 * TODO:
 *  Deprecate regAllocComposeRegion
 */
char* regAllocComposeRegion( const regRegion * region ) {
    char* regstr = regToStringRegion(region);

    // mimic previous behaviour
    if ( strcmp(regstr, "Null region") == 0)
      return NULL;
    else if ( strcmp(regstr, "Empty region") == 0)
      strcpy(regstr, "");

    return regstr;
}


char* regToStringRegion( const regRegion * region ) {

    long init_size = 512;
    long size;
    regShape *Shape;
    char *ptr;
    char *buf;
    long cpt = 0;

    size = init_size;
    buf = calloc(size, sizeof(char));
    ptr = buf;
    
    if (!region) {
        sprintf( buf, "Null region" );
        return buf;
    }

    Shape = region->shape;
    if(!Shape) {
        sprintf( buf, "Empty region" );
        return buf;
    }

    while (Shape) {
	    reg_compose_alloc_shape(Shape, cpt, &buf, &ptr, &size, 1);
	    cpt = Shape->component;
	    Shape = Shape->next;
    }
    return buf;
}

/*
 * Takes the shape, component of the previous shape (0 if there is no
 * previous shape), a pointer to the start the region string, a pointer
 * to the current head of the region string, the maximum size, and if
 * the function should allocate additional space (if needed).
 */
void reg_compose_alloc_shape(regShape * shape, long cpt, char **pbuf,
			  char **pptr, long *size, int alloc)
{
    // There's quite a bit of pointer arithmetic in this method to keep
    // pptr and pbuf pointed at the correct locations. Diagrams help.
    if (!shape) {
      // Do nothing
      return;
    }
  
    if (!*pptr || !*pbuf) {
      fprintf(stderr, "ERROR: No allocated space in buf for reg_compose_alloc_shape\n");
      return;
    }
    
    long tsize = *pptr - *pbuf;  // Length so far
    long maxlen = *size - tsize; // Length remaining
    char *shapestring; // holds the current shape

    // These are guesses as to the maximum length of an individual shape
    // given what we are shooting for in terms of accuracy.
    int safe = 80;
    if (shape->nPoints > 2) {
      safe = 40 + 20 * shape->nPoints;
    }
//    printf("MCD TEMP: reg_compose_alloc_shape()\n");
//    printf("MCD TEMP: + size = %ld\n",*size);
//    printf("MCD TEMP: + tsize = %ld\n",tsize);
//    printf("MCD TEMP: + maxlen = %ld\n",maxlen);
//    printf("MCD TEMP: + safe   = %d\n",safe);
    
    
    // Handle cases if we have gone over the max length.
    // (Length so far) + (estimated next size) > (REMAINING length) s/b (TOTAL length)
    while (tsize > *size - safe) {
      // If we can allocate resources
      if (alloc) {
	// Double array size
	*size = *size * 2; 
	*pbuf = realloc(*pbuf, *size * sizeof(char));
	
	maxlen = *size - tsize; // new max length
	*pptr = *pbuf + tsize;
	    }
      // Otherwise we add REGION_TRUNCATED_STR and return what we have
      else {
	safe = strlen( REGION_TRUNCATED_STR ) + 1;
	
	// Move the head back enough to print the truncated string
	while (*size - (*pptr - *pbuf) < safe) {
	  *pptr--;
	}
	
	//        printf("MCD TEMP: + before pbuf length = %ld\n", strlen(*pbuf) );
	if (*pptr - *pbuf >= 0) {
	  sprintf( *pptr, REGION_TRUNCATED_STR );
	  //          printf("MCD TEMP: + after pbuf length = %ld\n", strlen(*pbuf) );
	}
	else {
	  // Not enough space was allocated to even print ...truncated, this
	  // is obnoxious because it will show up for every shape (if you don't 
	  // stop when *pptr is null), but we believe it to be acceptable 
	  // behavior for the time being.
	  fprintf(stderr, 
		  "WARNING: Not enough space allocated to print region (%lu chars)",
		  *size);
	}
	
	// Set pptr to NULL to let caller know we're out of space
	      *pptr = NULL;
	      return;
      }
    }
    
    if (cpt > 0) {
      // Draw this out and I promise it makes sense
      *(*pptr)++ = (cpt == shape->component ? '&' : '|');
      tsize++;
    }
    
    shapestring = malloc(maxlen * sizeof(char));
    shape->toString(shape, shapestring, maxlen);  
//    printf("MCD TEMP: + maxlen (allowed) = %ld\n", maxlen );
//    printf("MCD TEMP: + new shape string = %s\n", shapestring );
//    printf("MCD TEMP: + new shape length (actual) = %ld\n", strlen(shapestring) );
    
    // Reusing tsize here
    tsize = snprintf(*pptr, maxlen, "%s", shapestring);
    free(shapestring);
    
    // Move pointer to end of the string 
    *pptr += maxlen < tsize ? maxlen : tsize;

    //    printf("MCD TEMP: + end pbuf = %s\n", *pbuf );
}


void regPrintRegion( regRegion* region )
{
    regShape *Shape;
    
    if (!region) {
        return;
    }

    if (region->shape == NULL) {
        return;
    }
    int ii = 1;
    Shape = region->shape;

    while (Shape != NULL) {
        regPrintShape(Shape);
        Shape = Shape->next;
    }
}


/*
 * Read a line from an ASCII file
 */
int reg_read_line(FILE * fp, char *buf, long maxlen)
{
    char *ptr;
    char *eptr;
    ptr = fgets(buf, maxlen, fp);
    if (!ptr) {
        fclose(fp);
        return 0;
    }

    if (!strncmp(buf, "SIMPLE  =", 9)) {
        /* FITS file; read another way */
        fclose(fp);
        return 0;
    }
    eptr = buf + strlen(buf) - 1;
    if (*eptr == '\n') {
        *eptr = '\0';
    }
    return 1;
}


/*
 *  Parse a line from the ASCII region file.
 *  We return the type (in variable &mode) attempting to recognize
 *  the region type.
 *  If the region type is already identified we use it to help parse.
 *  The string buffer is appended to.
 *  Eventually it contains the entire region which can be very long.
 */
void reg_parse_line(char *buf, long *mode, char **stringptr, long *maxlen, long *sys)
{
    char *ptr;
    char *optr;
    char sep = 0;
    int include = 1;

    char shape[SZ_LARGE];
    int state = 1;
    int parse_more;
    char *tptr;
    char *iptr;


    if (!*stringptr) {
	    *stringptr = calloc(SZ_LARGE, sizeof(char));
    }

    iptr = *stringptr;
    if (**stringptr) {
	    state = 0;
    }
    ptr = buf;

    /* If already decided it's ds9 mode, then skip remaining comments */
    if ((*mode == RF_DS9 || *mode == RF_DS9_V4) && ('#' == *ptr)) {
	    return;
    }

    if (!strncmp(ptr, "# filename", 10)) {	/* Possible saotng comment */
	    *mode = RF_SAOTNG;
	    return;
    } 
    else if (!strncmp(ptr, "# Region file format: DS9 version 4", 35)) {
        *mode = RF_DS9_V4;
    } 
    else if (!strncmp(ptr, "# Region file format: DS9", 25)) {
	    *mode = RF_DS9;
    } 
    else if (!strncmp(ptr, "XPA", 3)) {
	    *mode = RF_SAOTNG;
	    return;
    } else if (*ptr == '#') {
	    *mode = RF_SAOIMAGE;
	    return;
    } 
    else if (*ptr == '+' || *ptr == '-') {
	    if (state && (*mode != RF_DS9 && *mode != RF_DS9_V4)) { /* First entry */
	        *mode = RF_SAOTNG;
        }
	
        /* Parse this shape */
	    include = (*ptr == '+');
	    if (state) {
	        if (include) {
		        ptr++;
            }
	        else {
		        strcpy(*stringptr, "FIELD");
            }
	    }

	    ptr = reg_token_advance(ptr, shape, '#');

	    /* Strip trailing blanks */
	    tptr = shape + strlen(shape) - 1;
	    while (*tptr == ' ' || *tptr == '\t') {
	        *tptr-- = '\0';
        }
	    reg_strcat(stringptr, maxlen, 0, shape);
        
        /*
            if ( ptr && *ptr ) printf( "TNG Format: Ignoring %s\n", ptr );
        */
    } 
    else if (*mode == RF_DS9 || *mode == RF_DS9_V4) {
        /* In DS9 format, the "global" line gives us some crucial info. */
	    parse_more = 1;
	    while ( parse_more )
	    {
	        if (!strncmp(ptr, "global", 6)) {
		        if (!strncmp(ptr + 7, "coordsys=", 9)) {
		            ptr = ptr + 16;
		            if (!strncmp(ptr, "physical", 8)) {
			            *sys = RC_PHYSICAL;
                    }
		            else if (!strncmp(ptr, "image", 6)) {
			            *sys = RC_LOGICAL;
                    }
		            else if (!strncmp(ptr, "fk", 2)) {
			            *sys = RC_WORLD;
                    }
		        }
		        return;
	        }

	        if (*mode == RF_DS9_V4)
	        {
		        if (!strncmp(ptr, "physical", 8)) {
		            *sys = RC_PHYSICAL;
		            return;
		        }
		        else if (!strncmp(ptr, "image", 6)) {
		            *sys = RC_LOGICAL;
		            return;
		        }
		        else if (!strncmp(ptr, "fk", 2)) {
		            *sys = RC_WORLD;
		            return;
		        }
	        }
	    
            /* In DS9 format, the "local" line is ignored */
	        if (!strncmp(ptr, "local", 6)) {
		        return;
            }

            /* Otherwise we have a standard DS9 shape line. 
               First we strip off any coord system indicator */
	        if (*mode == RF_DS9) {
		        tptr = reg_token_advance(ptr, shape, ';');
		        if (reg_case_equal(shape, "PHYSICAL")) {
		            *sys = RC_PHYSICAL;
		            ptr = tptr;
		        } 
                else if (reg_case_equal(shape, "IMAGE")) {
		            *sys = RC_LOGICAL;
		            ptr = tptr;
		        } 
                else if (!strncmp(ptr, "fk", 2)) {
		            *sys = RC_WORLD;
		            ptr = tptr;
		        }
		        while (*ptr == ' ' || *ptr == ';') {
		            ptr++;
                }
            }

            /* Next we strip off the actual shape */
	        ptr = reg_tokens_advance(ptr, shape, "#;");
	        sep = state ? ' ' : '+';
	        state = 1;
	    
            /* Strip trailing blanks */
	        tptr = shape + strlen(shape) - 1;
	        while (*tptr == ' ' || *tptr == '\t') {
		        *tptr-- = '\0';
            }
	        if (*shape == '-') {
		        sep = ' ';
            }
            /* Append the shape to the region */
            reg_strcat(stringptr, maxlen, sep, shape);

            /* Need to check if there are more regions in the line after a ';' char */
	        if ( *ptr == '#' ) {
		        while ( *ptr != ';' && *ptr != '\0' )
		        ptr++;
	        }
	        if ( *ptr != ';' ) {
		        parse_more = 0;
            }
	        else {
		        ptr++;
		        while( *ptr == ' ' || *ptr == '\t' ) {
                    ptr++;
                }
		        if ( *ptr == '#' ) {
		            parse_more = 0;
                }
	        }
	    } // end while parsemore
    } // end else *mode == RF_DS9 || *mode == RF_DS9_V4
    else {
	    optr = ptr;
	    ptr = reg_token_advance(ptr, shape, ' ');
	    if (shape[strlen(shape) - 1] == ';') {
	        *mode = RF_PROS;	/* Pros regions have  coordsys; shape */
	        if (reg_case_equal(shape, "PHYSICAL;")) {
		        *sys = RC_PHYSICAL;
            }
	        else if (reg_case_equal(shape, "LOGICAL;")) {
		        *sys = RC_LOGICAL;
            }
	        else {
		        *sys = RC_WORLD;
            }
	        while (*ptr == ' ') {
		        ptr++;
            }
	    } else {
	        ptr = optr;		/* Reset pointer */
	        if (!strchr(shape, '(')) {
		        *mode = RF_PROS;
            }
	        *sys = RC_LOGICAL;
	    }
	    while (*ptr != 0) {
	        optr = ptr;
	        ptr = reg_token_advance(ptr, shape, ' ');
	        while (*ptr == ' ') {
		        ptr++;
            }
	        if (*shape == '^') {
		        *shape = '|';	/* We don't handle exclusive-or correctly */
            }
	        if (*shape == '&' || *shape == '|' || *shape == '!') {
		        sep = state ? 0 : ')';
		        reg_strcat(stringptr, maxlen, sep, shape);
		        state = 1;
	        } 
            else if (*shape == ';') {
		        if (!state) {
		            reg_strcat(stringptr, maxlen, ')', NULL);
                }

		        ptr = optr + 1;	/* Reparse since there may not be a space after the semi */
	        } 
            else if (isalpha(*shape)) {
		        sep = state ? ' ' : '+';
		        reg_strcat(stringptr, maxlen, sep, shape);
		        state = 1;
	        } 
            else if (isdigit(*shape) || *shape == '+' ||
		            (*shape == '-' && isdigit(shape[2]))) {
		        sep = state ? '(' : ',';
		        reg_strcat(stringptr, maxlen, sep, shape);
		        state = 0;
		        /* Check for celestial format values */
		        if (strchr(shape, ':') || strchr(shape, '\'')) {
		            *sys = RC_WORLD;
                }
		        if (shape[strlen(shape) - 1] == 'd') {
		            *sys = RC_WORLD;
                }
	        }
	    } // End while 
	    if (!state) {
	        reg_strcat(stringptr, maxlen, ')', NULL);
        }    
    }
}


void reg_areg_line(FILE * out, regShape * shape, char *shapeName, long nr, long nangles,
        char *text, int world)
{
    long maxlen = 2048;
    char *buf;
    char *ptr;
    long i;

    double *xpos = shape->xpos;
    double *ypos = shape->ypos;
    long np = shape->nPoints;
    double *radii = shape->radius;
    double *angles = shape->angle;

    buf = calloc(maxlen + (20 * np), sizeof(char));
    ptr = buf;

    if (!strcmp(shapeName, "Rectangle") || !strcmp(shapeName, "-Rectangle")) {
        double xcen, ycen, xlen, ylen;
        double ang = 0.0;
        char *shp = (*shapeName == '-' ) ? "-Box" : "Box";
        xcen = ( xpos[0]+xpos[1] ) / 2.0;
        ycen = ( ypos[0]+ypos[1] ) / 2.0;
        xlen = fabs(xpos[0]-xpos[1]);
        ylen = fabs(ypos[0]-ypos[1]);
      
        ptr += sprintf(buf, "%s(%g,%g,%g,%g,%g) # %s", shp, xcen,
		       ycen, xlen, ylen, ang, text);
      
    }
    else {
        ptr += sprintf(ptr, "%s(", shapeName);

        // Add Points
	for (i = 0; i < np; i++) {
	  ptr += sprintf(ptr, "%g,%g,", xpos[i], ypos[i]);
	}
	
        // Add Radii
        for (i = 0; i < nr; i++) {
            ptr += sprintf(ptr, "%g,", radii[i]);
        }

        // Add Angles
        for (i = 0; i < nangles; i++) {
            ptr += sprintf(ptr, "%g,", angles[i]);
        }

	ptr--;
	sprintf(ptr, ") # %s", text);
    }

    fprintf(out, "%s\n", buf);
    free(buf);
}


/*
 *  Append a separator and a string to the string buffer.
 *  Reallocate the string buffer as needed.
 */
void reg_strcat(char **ptr, long *maxlen, char sep, char *buf)
{
    char *string;
    long next;
    char *pos;
    long need = 2;
    string = *ptr;
    
    if (buf) {
	    need += strlen(buf);
    }
    next = strlen(string);
    
    if (next + need > *maxlen) {
	    *maxlen += SZ_LARGE;
	    *ptr = realloc(string, *maxlen * sizeof(char));
	    string = *ptr;
    }
    pos = string + next;
    
    if (sep != '\0') {
        *pos++ = sep;
    }
    
    if (buf) {
        strcpy(pos, buf); 
    }
    else {
        *pos = '\0';
    }
}


char *reg_token_advance(char *ptr, char *colname, char sep)
{
    char *optr = colname;
    while (*ptr != sep && *ptr != '\0') {
        *optr++ = *ptr++;
    }
    
    *optr = '\0';
    return ptr;
}


char *reg_tokens_advance(char *ptr, char *colname, char *seps)
{
    char *optr = colname;
    while (!strchr(seps, *ptr)) {
        *optr++ = *ptr++;
    }
    *optr = '\0';
    return ptr;
}


void reg_areg_hdr(FILE * out, char *color)
{
    char buf[512];
    fprintf(out, "%s\n", "# Region file format: DS9 version 4.1");
    sprintf(buf,
	    "global color=%s dashlist=8 3 width=1 font=\"helvetica 10 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1", color);
    fprintf(out, "%s\n", buf);
    fprintf(out, "physical\n");   // now a global property in ds9 v4
}

void reg_print_pos_pair(double x, double y, int world, char *xbuf, char *ybuf)
{
    if (world == RC_WORLD)
	reg_print_pos(x / 15.0, world, xbuf);
    else
	reg_print_pos(x, world, xbuf);

    reg_print_pos(y, world, ybuf);
}

void reg_print_pos(double x, int world, char *buf)
{
    char tmp[80];
    char *ptr;
    long h, m;
    double s;
    double fs;
    long is;
    int south;
    if (world == RC_WORLD) {
	south = x < 0;
	s = x * 3600.0;
	if (south)
	    s = -s;
	is = (long) s;
	if (s - is > 0.9999) {	/* What's important is the precision we'll print to */
	    is++;
	    fs = 0;
	} else {
	    fs = s - is;
	}
	m = is / 60;
	is = is % 60;
	h = m / 60;
	m = m % 60;
	reg_print_val(fs, tmp);
	ptr = tmp;
	if (*ptr == '0')
	    ptr++;
	if (south)
	    sprintf(buf, "-%02ld:%02ld:%02ld%s", h, m, is, ptr);
	else
	    sprintf(buf, "%02ld:%02ld:%02ld%s", h, m, is, ptr);
    } else {
	reg_print_val(x, buf);
    }
}

void reg_print_radius(double r, int world, char *buf)
{
    double x;
    if (world == RC_WORLD && r < 1.0) {
	    x = r * 60.0;
	    reg_print_val(x, buf);
	    strcat(buf, "'");
    } else
	    reg_print_val(r, buf);
}

/*
 *  Pretty print a double value, stripping trailing zeros from
 *  mantissa.
 */
char *reg_print_val(double x, char *buf)
{
    char tmp[80];
    char *ptr;
    char *decimal;
    long j = 0;

    sprintf(tmp, "%g", x);
    ptr = strpbrk(tmp, "eE");
    if (ptr)
	j = ptr - tmp;
    else
	j = strlen(tmp) - 1;

/* Strip trailing mantissa zeros */
    decimal = strchr(tmp, '.');
    if (decimal && j > decimal - tmp) {

	while (j > 0 && tmp[j] == '0')
	    tmp[j--] = '\0';

	if (j > 0 && tmp[j] == '.')
	    tmp[j] = '\0';
    }

    strcpy(buf, tmp);
    if (ptr && j < ptr - tmp)
	strcat(buf, ptr);

    return buf;
}
