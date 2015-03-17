/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/
#include "region_priv.h"
#include <string.h>
#include <ctype.h>

#define REGION_TRUNCATED_STR "...[truncated]"

void regComposeAllocShape(regShape * shape, long cpt, char **pbuf,
			  char **pptr, long *size, int alloc);
void reg_areg_line(FILE * out, char *shapeName, double *pos, long nmax,
		   long np, double *r, long nr, double *angle, long nangle,
		   char *text, int world);
void reg_areg_hdr(FILE * out, char *color);
void regComposeShape(regShape * Shape, long cpt, char *buf, char **pptr,
		     long maxlen);

/*
 *  return the maximum number of points of any
 *  shape in the region. Usually 2, but if polygons
 *  are present may be arb.
 */
long regGetMaxPoints(const regRegion * region)
{
    regShape *Shape;
    long n = 0;
    if (!region)
	return 0;

    Shape = region->shape;
    while (Shape) {
	if (Shape->nPoints > n)
	    n = Shape->nPoints;

	Shape = Shape->next;
    }
    return n;
}

long regGetNoShapes(const regRegion * region)
{
    regShape *Shape;
    long n = 0;
    if (!region)
	return 0;

    Shape = region->shape;
    while (Shape) {
	n++;
	Shape = Shape->next;
    }
    return n;
}

regShape *regGetShapeNo(const regRegion * region, long shapeNo)
{
    regShape *Shape;
    long no = 1;

    if (!region)
	return NULL;
    Shape = region->shape;
    while (no < shapeNo) {
	if (!Shape)
	    return NULL;
	Shape = Shape->next;
	no++;
    }
    return Shape;
}

regRegion *regShapeGetRegion(const regShape * shape)
{
    if (!shape)
	return NULL;

    return shape->region;
}

long regShapeGetNoPoints(const regShape * shape)
{
    if (!shape)
	return 0;

    return shape->nPoints;
}

long regShapeGetComponent(const regShape * shape)
{
    if (!shape)
	return 0;

    return (long) shape->component;
}

/*
 *  Return the values of radius in the preallocated array.
 *  Return the number of radii as the function value.
 *  nradii can be 0, 1 or 2.
 */
long regShapeGetRadii(const regShape * shape, double *r)
{
    long nradii;
    long i;
    if (!shape || !r)
	return 0;

    nradii = regShapeRadii(shape);
    for (i = 0; i < nradii; i++)
	r[i] = shape->radius[i];
    return nradii;
}


long regShapeGetAngles(const regShape * shape, double *angle)
{
    long nangles;
    long i;
    if (!shape || !angle)
	return 0;

    nangles = regShapeAngles(shape);
    for (i = 0; i < nangles; i++)
	angle[i] = shape->angle[i];
    return nangles;
}

/*
 *  Given an array of positions of size nmax, find the number
 *  of positions actually used. For polygons, this is found because
 *  the last and first points must be identical.
 */
long regShapeFindNPoints(regGeometry type, double *xpos, double *ypos,
			 long nmax)
{

    long n = 0;
    double x, y;
    int loop;
    switch (type) {
    case regPOLYGON:
	x = xpos[0];
	y = ypos[0];
	n = 1;
	loop = 1;
	while (n < nmax && loop) {
	  if (( x == xpos[n]) && (y == ypos[n]))
	    loop = 0;
	  else
	    n++;
	}
	if (loop)
	    n = nmax;
	break;
    case regPOINT:
    case regCIRCLE:
    case regANNULUS:
    case regELLIPSE:
    case regBOX:
    case regPIE:
    case regSECTOR:
	n = 1;
	break;
    case regRECTANGLE:
	n = 2;
	break;
    case regFIELD:
    case regMASK:
    case regUSER:
    default:
	n = 0;
	break;
    }
    return n;
}



long regShapeGetPoints(const regShape * shape, double *x, double *y,
		       long dim)
{
    long n;
    long i;
    if (!shape || !x || !y || dim <= 0)
	return 0;

    if (!shape->xpos || !shape->ypos || shape->nPoints <= 0)
	return 0;

    n = shape->nPoints;
/* Only return as many points as are in the array */
    if (n > dim)
	n = dim;

    for (i = 0; i < n; i++) {
	x[i] = shape->xpos[i];
	y[i] = shape->ypos[i];
    }
/* Zero any remaining array elements */
    for (i = n; i < dim; i++) {
	x[i] = 0.0;
	y[i] = 0.0;
    }
    return n;
}


int regCaseEqual(char *s1, char *s2)
{
    while (1) {
	if (toupper(*s1) != toupper(*s2))
	    return 0;
	if (*s1 == '\0' || *s2 == '\0')
	    return 1;
	s1++;
	s2++;
    }
}


/* Return name of shape. Also return whether included or excluded,
  as function return value
 */

int regShapeGetName(const regShape * shape, char *name, long maxlen)
{

  if ( NULL == shape ) {
    strncpy(name, "Unknown", maxlen);
    return(0);
  }

    strcpy(name, "");
    switch (shape->type) {
    case regCIRCLE:
	strncpy(name, "Circle", maxlen);
	break;
    case regPOINT:
	strncpy(name, "Point", maxlen);
	break;
    case regELLIPSE:
	strncpy(name, "Ellipse", maxlen);
	break;
    case regPIE:
	strncpy(name, "Pie", maxlen);
	break;
    case regSECTOR:
	strncpy(name, "Sector", maxlen);
	break;
    case regANNULUS:
	strncpy(name, "Annulus", maxlen);
	break;
    case regPOLYGON:
	strncpy(name, "Polygon", maxlen);
	break;
    case regBOX:
	strncpy(name, "Box", maxlen);
	if (shape->angle[0] != 0.0)
	    strncpy(name, "RotBox", maxlen);

	break;
    case regRECTANGLE:
	strncpy(name, "Rectangle", maxlen);
	if (shape->angle[0] != 0.0)
	    strncpy(name, "RotRectangle", maxlen);
	break;
    case regFIELD:
	strncpy(name, "Field", maxlen);
	break;
    case regMASK:
    case regUSER:   
       strncpy(name, shape->user->name, maxlen);
       break;
    default:
	strncpy(name, "Unknown", maxlen);
	break;
    }
    if (shape->include == regExclude)
	return 0;
    else
	return 1;
}

regGeometry regShapeNameToGeometry(char *name)
{
    long ntypes = 12;
    char *names[] = { "Circle", "Point", "Ellipse", "Pie", "Sector",
	"Annulus", "Polygon",
        "Box", "Rectangle", "RotBox", "RotRectangle", "Field"
    };

    regGeometry geoms[] =
	{ regCIRCLE, regPOINT, regELLIPSE, regPIE, regSECTOR,
	regANNULUS, regPOLYGON, regBOX, regRECTANGLE,
	regBOX, regRECTANGLE, regFIELD
    };
    long i;

    for (i = 0; i < ntypes; i++) {
	if (regCaseEqual(name, names[i]))
	    return geoms[i];
    }
    /* Unknown case */
    return regPOINT;
}


void regComposeRegion(const regRegion * region, char *buf, long maxlen)
{
    regShape *Shape;
    char *ptr;
    long cpt = 0;
    strcpy(buf, "");
    if (!region)
	return;

    ptr = buf;
    Shape = region->shape;
    while (Shape != NULL && ptr) {
	/* ptr will be set to NULL if we run out of space */
	regComposeShape(Shape, cpt, buf, &ptr, maxlen);
	cpt = Shape->component;
	Shape = Shape->next;
    }
}

char *regAllocComposeRegion(const regRegion * region)
{
    long init_size = 512;
    long size;
    regShape *Shape;
    char *ptr;
    char *buf;
    long cpt = 0;
    if (!region)
	return NULL;

    size = init_size;
    buf = calloc(size, sizeof(char));
    ptr = buf;

    Shape = region->shape;
    while (Shape != NULL) {
	regComposeAllocShape(Shape, cpt, &buf, &ptr, &size, 1);
	cpt = Shape->component;
	Shape = Shape->next;
    }
    return buf;
}

void regComposeShape(regShape * shape, long cpt, char *buf, char **pptr,
		     long maxlen)
{
    regComposeAllocShape(shape, cpt, &buf, pptr, &maxlen, 0);
}

void regComposeAllocShape(regShape * shape, long cpt, char **pbuf,
			  char **pptr, long *size, int alloc)
{
    char *buf = *pbuf;
    char *ptr = *pptr;
    long ii;
    long i;
    long safe = 80;		/* Space for writing text; allow for longer precision now used */
    long bsize = 80;		/* Space for writing text */
    long nradii;
    long nangles;
    char *radius[2];
    char *angle[2];
    char **xpoint = NULL;
    char **ypoint = NULL;
    long maxlen = *size;
    long tsize;
    char *tmpstr;

    if (shape && shape->nPoints > 2)
	safe = 40 + 20 * shape->nPoints;
    if (ptr - buf > maxlen - safe) {
	if (alloc) {
	    tsize = ptr - buf;	/* Length used */

	    if (maxlen < safe)
		maxlen = safe * 2;
	    else
		maxlen = maxlen * 2;	/* New max length */
	    buf = realloc(buf, maxlen * sizeof(char));
	    *size = maxlen;	/* Return new max length */
	    *pbuf = buf;	/* Return new buffer start */
	    ptr = buf + tsize;
	    *pptr = ptr;	/* Update new buffer position */
	} else {
	    /* String is fixed size, and we are out of space - sorry! */

	    safe = strlen( REGION_TRUNCATED_STR ) + 1;
	    while (ptr - buf > maxlen - safe)
		ptr--;

	    if (ptr - buf >= 0)
	        sprintf( ptr, REGION_TRUNCATED_STR );

	    /* Set pptr to NULL to let caller know we're out of space */
	    *pptr = NULL;

	    return;
	}
    }
 
    if (cpt > 0) {
	if (cpt == shape->component)
	    *ptr++ = '&';
	else
	    *ptr++ = '|';
    }

    if (shape->include == regExclude)
	*ptr++ = '!';

    nradii = regShapeRadii(shape);
    nangles = regShapeAngles(shape);
    for (i = 0; i < nradii; i++) {
	radius[i] = calloc(bsize, sizeof(char));
	regPrintRadius(shape->radius[i], shape->world_size, radius[i]);
    }

    for (i = 0; i < nangles; i++) {
	angle[i] = calloc(bsize, sizeof(char));
	regPrintAngle(shape->angle[i], angle[i]);
    }

    if (shape->nPoints > 0) {
	xpoint = (char **) calloc(shape->nPoints, sizeof(char *));
	ypoint = (char **) calloc(shape->nPoints, sizeof(char *));
    }
    for (i = 0; i < shape->nPoints; i++) {
	xpoint[i] = calloc(bsize, sizeof(char));
	ypoint[i] = calloc(bsize, sizeof(char));
	regPrintPosPair(shape->xpos[i], shape->ypos[i],
			shape->world_coord, xpoint[i], ypoint[i]);

    }

    switch (shape->type) {
    case regPOINT:
	ptr += sprintf(ptr, "Point(%s,%s)", xpoint[0], ypoint[0]);
	break;

    case regCIRCLE:
	ptr += sprintf(ptr, "Circle(%s,%s,%s)",
		       xpoint[0], ypoint[0], radius[0]);
	break;

    case regANNULUS:
	ptr += sprintf(ptr, "Annulus(%s,%s,%s,%s)",
		       xpoint[0], ypoint[0], radius[0], radius[1]);
	break;

    case regELLIPSE:
	ptr += sprintf(ptr, "Ellipse(%s,%s,%s,%s,%s)",
		       xpoint[0], ypoint[0], radius[0], radius[1],
		       angle[0]);
	break;
    case regPOLYGON:
	ptr += sprintf(ptr, "Polygon(");
	for (ii = 0; ii < shape->nPoints; ii++) {
	    if (ii == 0)
		ptr += sprintf(ptr, "%s,%s", xpoint[ii], ypoint[ii]);
	    else
		ptr += sprintf(ptr, ",%s,%s", xpoint[ii], ypoint[ii]);
	}
	ptr += sprintf(ptr, ")");
	break;

    case regPIE:
	ptr += sprintf(ptr, "Pie(%s,%s,%s,%s,%s,%s)",
		       xpoint[0], ypoint[0], radius[0], radius[1],
		       angle[0], angle[1]);
	break;

    case regSECTOR:
	ptr += sprintf(ptr, "Sector(%s,%s,%s,%s)",
		       xpoint[0], ypoint[0], angle[0], angle[1]);
	break;

    case regBOX:

	if (shape->angle[0] == 0.0)
	    ptr += sprintf(ptr, "Box(%s,%s,%s,%s)",
			   xpoint[0], ypoint[0], radius[0], radius[1]);
	else
	    ptr += sprintf(ptr, "RotBox(%s,%s,%s,%s,%s)",
			   xpoint[0], ypoint[0], radius[0], radius[1],
			   angle[0]);

	break;

    case regRECTANGLE:
	if (shape->angle[0] == 0.0)
	    ptr += sprintf(ptr, "Rectangle(%s,%s,%s,%s)",
			   xpoint[0], ypoint[0], xpoint[1], ypoint[1]);
	else
	    ptr += sprintf(ptr, "RotRectangle(%s,%s,%s,%s,%s)",
			   xpoint[0], ypoint[0], xpoint[1], ypoint[1],
			   angle[0]);
	break;


    case regFIELD:
	ptr += sprintf(ptr, "Field()");
	break;

    case regMASK:
    case regUSER:
      tmpstr = shape->user->toString( shape );      
      ptr += sprintf(ptr, tmpstr);
      free(tmpstr);
      break;

    }

/* Free the strings */
    for (i = 0; i < shape->nPoints; i++) {
	free(xpoint[i]);
	free(ypoint[i]);
    }
    free(xpoint);
    free(ypoint);
    for (i = 0; i < nangles; i++)
	free(angle[i]);
    for (i = 0; i < nradii; i++)
	free(radius[i]);


    *pptr = ptr;
}


void regPrintRegion(regRegion * region)
{
    regShape *Shape;

    if (!region) {
      printf("Null Region\n");
    } else {
      if (region->xcol[0] && region->xcol[1])
	printf("Region: Pointers x= %p y=%p\n", region->xcol[0],
	       region->xcol[1]);



      printf(" \n");
      Shape = region->shape;
      while (Shape != NULL) {
	regPrintShape(Shape);
	Shape = Shape->next;
      }
    }
}

/*
  +-----------------------------------------------------------------------
  +-----------------------------------------------------------------------
*/
void regPrintShape(regShape * shape)
{
    int ii;
    char buf[120];
    char *coordunit[] = { "pixel", "deg" };
    char wbuf[40];
    char *tmpstr;

    if(!shape)
      return;
    printf("%ld\t", shape->component);       
 
    if (shape->include == regExclude)
	printf("!");
    sprintf(wbuf, "(Pos: %s, Size: %s)", coordunit[shape->world_coord],
	    coordunit[shape->world_size]);

    switch (shape->type) {
    case regPOINT:
	sprintf(buf, "Point(%f, %f) %s", shape->xpos[0], shape->ypos[0],
		wbuf);
	break;

    case regCIRCLE:
	sprintf(buf, "Circle(%f, %f, %f) %s",
		shape->xpos[0], shape->ypos[0], shape->radius[0], wbuf);
	break;

    case regANNULUS:
	sprintf(buf, "Annulus(%f, %f, %f, %f) %s",
		shape->xpos[0], shape->ypos[0], shape->radius[0],
		shape->radius[1], wbuf);
	break;

    case regELLIPSE:
	sprintf(buf, "Ellipse(%f, %f, %f, %f, %f) %s",
		shape->xpos[0], shape->ypos[0], shape->radius[0],
		shape->radius[1], shape->angle[0], wbuf);
	break;
    case regPOLYGON:
	printf("POLYGON     %s\t", wbuf);
	for (ii = 0; ii < shape->nPoints; ii++)
	    printf("(%f,%f)\t", shape->xpos[ii], shape->ypos[ii]);
	sprintf(buf, " ");
	break;
    case regPIE:
	sprintf(buf, "Pie(%f, %f, %f, %f, %f, %f) %s",
		shape->xpos[0], shape->ypos[0], shape->radius[0],
		shape->radius[1], shape->angle[0], shape->angle[1], wbuf);
	break;
    case regSECTOR:
	sprintf(buf, "Sector(%f, %f, %f, %f) %s",
		shape->xpos[0], shape->ypos[0], shape->angle[0],
		shape->angle[1], wbuf);
	break;

    case regBOX:

	if (shape->angle[0] == 0.0)
	    sprintf(buf, "Box(%f, %f, %f, %f) %s",
		    shape->xpos[0], shape->ypos[0], shape->radius[0],
		    shape->radius[1], wbuf);
	else
	    sprintf(buf, "RotBox(%f, %f, %f, %f, %f) %s",
		    shape->xpos[0], shape->ypos[0], shape->radius[0],
		    shape->radius[1], shape->angle[0], wbuf);

	break;

    case regRECTANGLE:
	if (shape->angle[0] == 0.0)
	    sprintf(buf, "Rectangle(%f, %f, %f, %f) %s",
		    shape->xpos[0], shape->ypos[0], shape->xpos[1],
		    shape->ypos[1], wbuf);
	else
	    sprintf(buf, "RotRectangle(%f, %f, %f, %f, %f) %s",
		    shape->xpos[0], shape->ypos[0], shape->xpos[1],
		    shape->ypos[1], shape->angle[0], wbuf);
	break;

    case regFIELD:
	sprintf(buf, "Field() %s", wbuf);
	break;


    case regMASK:
    case regUSER:
      tmpstr = shape->user->toString( shape );
      sprintf(buf,  tmpstr);
      free(tmpstr);
      break;
    }

    printf("%s\n", buf);

    if (shape->region == NULL)
	return;
}


/*
 *  Returns allocated double arrays *xat and *yat
 */
int regRegionToList(regRegion * region,
		    double xmin,
		    double xmax,
		    double ymin,
		    double ymax,
		    double bin, double **xat, double **yat, long *nat)
{
    long ii, jj;
    double xpos, ypos;
    double xlen = (xmax - xmin) / bin + 1;
    double ylen = (ymax - ymin) / bin + 1;
    int stat;
    long Buffer = 200;

    *nat = 0;

    *xat = (double *) calloc(Buffer, sizeof(double));
    *yat = (double *) calloc(Buffer, sizeof(double));


    for (ii = 0; ii < xlen; ii++) {
	xpos = xmin + bin * ii;
	for (jj = 0; jj < ylen; jj++) {
	    ypos = ymin + bin * jj;

	    stat = regInsideRegion(region, xpos, ypos);

	    if (stat) {

		*nat += 1;
		if ((*nat >= Buffer)) {
		    Buffer = 2 * (*nat / Buffer + 1) * Buffer;
		    *xat =
			(double *) realloc(*xat, Buffer * sizeof(double));
		    *yat =
			(double *) realloc(*yat, Buffer * sizeof(double));
		}

		*(*xat + (*nat - 1)) = xpos;
		*(*yat + (*nat - 1)) = ypos;

	    }

	}
    }


    return 0;
}


/*
 *  Returns allocated short array *mask
 */

int regRegionToMask(regRegion * region,
		    double xmin,
		    double xmax,
		    double ymin,
		    double ymax,
		    double bin, short **mask, long *xlen, long *ylen)
{
    long ii, jj;

    *xlen = (xmax - xmin) / bin + 1;
    *ylen = (ymax - ymin) / bin + 1;

    *mask = (short *) calloc((*xlen) * (*ylen), sizeof(short));
    if (*mask == NULL)
	return (-1);

    for (ii = 0; ii < *xlen; ii++)
	for (jj = 0; jj < *ylen; jj++) {
	    *(*mask + (ii + jj * (*xlen))) =
		regInsideRegion(region, xmin + bin * ii, ymin + bin * jj);

	}

    return (0);
}


void regPrintPosPair(double x, double y, int world, char *xbuf, char *ybuf)
{
    if (world)
	regPrintPos(x / 15.0, world, xbuf);
    else
	regPrintPos(x, world, xbuf);

    regPrintPos(y, world, ybuf);
}

void regPrintPos(double x, int world, char *buf)
{
    char tmp[80];
    char *ptr;
    long h, m;
    double s;
    double fs;
    long is;
    int south;
    if (world) {
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
	regPrintVal(fs, tmp);
	ptr = tmp;
	if (*ptr == '0')
	    ptr++;
	if (south)
	    sprintf(buf, "-%02ld:%02ld:%02ld%s", h, m, is, ptr);
	else
	    sprintf(buf, "%02ld:%02ld:%02ld%s", h, m, is, ptr);
    } else {
	regPrintVal(x, buf);
    }
}

void regPrintRadius(double r, int world, char *buf)
{
    double x;
    if (world && r < 1.0) {
	x = r * 60.0;
	regPrintVal(x, buf);
	strcat(buf, "'");
    } else
	regPrintVal(r, buf);
}

void regPrintAngle(double x, char *buf)
{
    regPrintVal(x, buf);
}

/*
 *  Pretty print a double value, stripping trailing zeros from
 *  mantissa.
 */
char *regPrintVal(double x, char *buf)
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


regRegion *regReadAsciiRegion(char *filename, int verbose)
{
    regRegion *region = NULL;
    FILE *fp = NULL;

    char buf[SZ_LARGE];
    long maxlen = SZ_LARGE;
    char *ptr;
    long sys = RC_UNK;
    char *string = NULL;
    long mode = RF_UNK;
    char *smodes[] = { "UNK", "PHYSICAL", "LOGICAL", "WORLD" };
    char *modes[] = { "UNK", "SAOIMAGE", "SAOTNG", "PROS", "CXC", "DS9",
		      "DS9_V4" };
    fp = fopen(filename, "r");
    if (!fp)
	return NULL;

    if (verbose >= 1)
	fprintf(stderr, "cxcregion: parsing region file %s\n", filename);


    string = calloc(SZ_LARGE, sizeof(char *));
    while (reg_read_line(fp, buf, SZ_LARGE)) {
	ptr = buf;
	
	if ((ptr[0] == '#') && (mode != RF_UNK))
	  continue;

	while (*ptr == ' ')
	    ptr++;
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
 *  Parse a line from the ASCII region file.
 *  We return the type (in variable &mode) attempting to recognize
 *  the region type.
 *  If the region type is already identified we use it to help parse.
 *  The string buffer is appended to.
 *  Eventually it contains the entire region which can be very long.
 */
void reg_parse_line(char *buf, long *mode, char **stringptr, long *maxlen,
		    long *sys)
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


    if (!*stringptr)
	*stringptr = calloc(SZ_LARGE, sizeof(char));

    iptr = *stringptr;
    if (**stringptr)
	state = 0;
    ptr = buf;

    /* If already decided it's ds9 mode, then skip remaining comments */
    if ((*mode == RF_DS9 || *mode == RF_DS9_V4) && ('#' == *ptr))
	return;

    if (!strncmp(ptr, "# filename", 10)) {	/* Possible saotng comment */
	*mode = RF_SAOTNG;
	return;
    } else if (!strncmp(ptr, "# Region file format: DS9 version 4", 35)) {
        *mode = RF_DS9_V4;
    } else if (!strncmp(ptr, "# Region file format: DS9", 25)) {
	*mode = RF_DS9;
    } else if (!strncmp(ptr, "XPA", 3)) {
	*mode = RF_SAOTNG;
	return;
    } else if (*ptr == '#') {
	*mode = RF_SAOIMAGE;
	return;
    } else if (*ptr == '+' || *ptr == '-') {
	if (state && (*mode != RF_DS9 && *mode != RF_DS9_V4)) /* First entry */
	    *mode = RF_SAOTNG;
	/* Parse this shape */
	include = (*ptr == '+');
	if (state) {
	    if (include)
		ptr++;
	    else
		strcpy(*stringptr, "FIELD");
	}

	ptr = reg_token_advance(ptr, shape, '#');

	/* Strip trailing blanks */
	tptr = shape + strlen(shape) - 1;
	while (*tptr == ' ' || *tptr == '\t')
	    *tptr-- = '\0';
	reg_strcat(stringptr, maxlen, 0, shape);
/*
    if ( ptr && *ptr ) printf( "TNG Format: Ignoring %s\n", ptr );
 */
    } else if (*mode == RF_DS9 || *mode == RF_DS9_V4) {
/* In DS9 format, the "global" line gives us some crucial info. */
	parse_more = 1;
	while ( parse_more )
	{
	    if (!strncmp(ptr, "global", 6)) {
		if (!strncmp(ptr + 7, "coordsys=", 9)) {
		    ptr = ptr + 16;
		    if (!strncmp(ptr, "physical", 8))
			*sys = RC_PHYSICAL;
		    else if (!strncmp(ptr, "image", 6))
			*sys = RC_LOGICAL;
		    else if (!strncmp(ptr, "fk", 2))
			*sys = RC_WORLD;
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
	    if (!strncmp(ptr, "local", 6))
		return;

/* Otherwise we have a standard DS9 shape line. 
   First we strip off any coord system indicator */
	    
	    if (*mode == RF_DS9)
	      {
		tptr = reg_token_advance(ptr, shape, ';');
		if (regCaseEqual(shape, "PHYSICAL")) {
		  *sys = RC_PHYSICAL;
		  ptr = tptr;
		} else if (regCaseEqual(shape, "IMAGE")) {
		  *sys = RC_LOGICAL;
		  ptr = tptr;
		} else if (!strncmp(ptr, "fk", 2)) {
		  *sys = RC_WORLD;
		  ptr = tptr;
		}
		while (*ptr == ' ' || *ptr == ';')
		  ptr++;
	      }
/* Next we strip off the actual shape */
	    
	    ptr = reg_tokens_advance(ptr, shape, "#;");
	    sep = state ? ' ' : '+';
	    state = 1;
	    /* Strip trailing blanks */
	    tptr = shape + strlen(shape) - 1;
	    while (*tptr == ' ' || *tptr == '\t')
		*tptr-- = '\0';
	    if (*shape == '-')
		sep = ' ';
/* Append the shape to the region */
	    
	    reg_strcat(stringptr, maxlen, sep, shape);

/* Need to check if there are more regions in the line after a ';' char */
	    if ( *ptr == '#' )
	    {
		while ( *ptr != ';' && *ptr != '\0' )
		    ptr++;
	    }
	    if ( *ptr != ';' )
		parse_more = 0;
	    else
	    {
		ptr++;
		while( *ptr == ' ' || *ptr == '\t' ) ptr++;
		if ( *ptr == '#' )
		    parse_more = 0;
	    }
	}

    } else {
	optr = ptr;
	ptr = reg_token_advance(ptr, shape, ' ');
	if (shape[strlen(shape) - 1] == ';') {
	    *mode = RF_PROS;	/* Pros regions have  coordsys; shape */
	    if (regCaseEqual(shape, "PHYSICAL;"))
		*sys = RC_PHYSICAL;
	    else if (regCaseEqual(shape, "LOGICAL;"))
		*sys = RC_LOGICAL;
	    else
		*sys = RC_WORLD;
	    while (*ptr == ' ')
		ptr++;
	} else {
	    ptr = optr;		/* Reset pointer */
	    if (!strchr(shape, '('))
		*mode = RF_PROS;

	    *sys = RC_LOGICAL;
	}
	while (*ptr != 0) {
	    optr = ptr;
	    ptr = reg_token_advance(ptr, shape, ' ');
	    while (*ptr == ' ')
		ptr++;
	    if (*shape == '^')
		*shape = '|';	/* We don't handle exclusive-or correctly */
	    if (*shape == '&' || *shape == '|' || *shape == '!') {
		sep = state ? 0 : ')';
		reg_strcat(stringptr, maxlen, sep, shape);
		state = 1;
	    } else if (*shape == ';') {
		if (!state)
		    reg_strcat(stringptr, maxlen, ')', NULL);

		ptr = optr + 1;	/* Reparse since there may not be a space after the semi */
	    } else if (isalpha(*shape)) {
		sep = state ? ' ' : '+';
		reg_strcat(stringptr, maxlen, sep, shape);
		state = 1;
	    } else if (isdigit(*shape) || *shape == '+' ||
		       (*shape == '-' && isdigit(shape[2]))) {
		sep = state ? '(' : ',';
		reg_strcat(stringptr, maxlen, sep, shape);
		state = 0;
		/* Check for celestial format values */
		if (strchr(shape, ':') || strchr(shape, '\''))
		    *sys = RC_WORLD;
		if (shape[strlen(shape) - 1] == 'd')
		    *sys = RC_WORLD;
	    }
	}
	if (!state)
	    reg_strcat(stringptr, maxlen, ')', NULL);
    }
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
    if (buf)
	need += strlen(buf);
    next = strlen(string);
    if (next + need > *maxlen) {
	*maxlen += SZ_LARGE;
	*ptr = realloc(string, *maxlen * sizeof(char));
	string = *ptr;
    }
    pos = string + next;
    if (sep != '\0')
	*pos++ = sep;
    if (buf)
	strcpy(pos, buf);
    else
	*pos = '\0';
}

char *reg_token_advance(char *ptr, char *colname, char sep)
{
    char *optr = colname;
    while (*ptr != sep && *ptr != '\0')
	*optr++ = *ptr++;
    *optr = '\0';
    return ptr;
}

char *reg_tokens_advance(char *ptr, char *colname, char *seps)
{
    char *optr = colname;
    while (!strchr(seps, *ptr))
	*optr++ = *ptr++;
    *optr = '\0';
    return ptr;
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
    if (*eptr == '\n')
	*eptr = '\0';
    return 1;
}


int regWriteAsciiRegion(char *name, regRegion * region, char **names,
			long nn)
{

    int world = 0;
    char color[] = "blue";
    long arrayDim = 0;
    long n;
    long axlen[1] = { 1 };
    double *pos = NULL;
    FILE *out;
    double *x;
    double *y;
    long no;
    double r[2];
    double angle[2];
    char shapeName[20];
    char tmpname[20];
    long maxlen = 20;
    long nradii, nangles;
    int include;
    long npoints;
    short cpt = 1;
    char text[256];
    regShape *Shape;

    arrayDim = regGetMaxPoints(region);
    if (arrayDim < 1)
	arrayDim = 1;

    if (( strcasecmp( name, "stdout" ) == 0 ) ||
	( strcmp( name, "-" ) == 0 )) 
      out = stdout;
    else if ( strcasecmp( name, "stderr") == 0 )
      out = stderr;
    else  out = fopen(name, "w");
    if (!out)
      return 0;

    n = regGetNoShapes(region);
    axlen[0] = arrayDim;

    pos = (double *) calloc(2 * arrayDim, sizeof(double));
    x = &pos[0];
    y = &pos[arrayDim];

    reg_areg_hdr(out, color);

    for (no = 1; no <= n; no++) {
	r[0] = 0.0;
	r[1] = 0.0;
	angle[0] = 0.0;
	angle[1] = 0.0;
	Shape = regGetShapeNo(region, no);
	if (Shape) {
	    include = regShapeGetName(Shape, tmpname, maxlen);
	    if (!include)
		sprintf(shapeName, "-%s", tmpname);
	    else
		strcpy(shapeName, tmpname);

	    npoints = regShapeGetPoints(Shape, x, y, arrayDim);
	    nradii = regShapeGetRadii(Shape, r);
	    nangles = regShapeGetAngles(Shape, angle);
	    cpt = regShapeGetComponent(Shape);
	    if (names && no <= nn)
		sprintf(text, "text=\"%s\"", names[no - 1]);
	    else
		strcpy(text, " ");

	    /* wow ... is this wrong */
	    reg_areg_line(out, shapeName, pos, arrayDim, npoints, r, nradii,
			  angle, nangles, text, world);
	}
    }

    if (pos)
	free(pos);
    fclose(out);
    return 1;
}




void reg_areg_line(FILE * out, char *shapeName, double *pos, long nmax,
		   long np, double *r, long nr, double *angle, long nangle,
		   char *text, int world)
{

    long maxlen = 2048;
    long safe = 100;
    char *buf;
    char *ptr;
    long dp;
    buf = calloc(maxlen, sizeof(char));

    /* 
       circle(x,y,r)
       box   (x,y,xl,yl)
       rotbox(x,y,xl,yl,a)
       ellipse(x,y,xl,yl,a)
       field()
     */

    if (!strcmp(shapeName, "Polygon") || !strcmp(shapeName, "-Polygon")) {
	int i;
	ptr = buf;
	ptr += sprintf(ptr, "physical;%s(", shapeName);
	for (i = 0; i < np; i++) {
	    ptr += sprintf(ptr, "%g,%g,", pos[i], pos[nmax + i]);
	    dp = ptr - buf;
	    if (dp > maxlen - safe) {
		maxlen = 2 * maxlen;
		buf = realloc(buf, maxlen * sizeof(char));
		ptr = buf + dp;
	    }
	}
	ptr--;
	sprintf(ptr, ") # %s", text);
    } else if (!strcmp(shapeName, "Rectangle") || !strcmp(shapeName, "-Rectangle")) {
      double xcen, ycen, xlen, ylen;
      double ang = 0.0;
      char *shp = (*shapeName == '-' ) ? "-Box" : "Box";
      xcen = ( pos[0]+pos[1] ) /2.0;
      ycen = ( pos[nmax]+pos[nmax+1] ) /2.0;
      xlen = fabs(pos[0]-pos[1]);
      ylen = fabs(pos[nmax]-pos[nmax+1]);
      
      sprintf(buf, "physical;%s(%g,%g,%g,%g,%g) # %s", shp, xcen,
	      ycen, xlen, ylen, ang, text);
      
    } else if (nr == 2 && nangle == 1) {
	sprintf(buf, "physical;%s(%g,%g,%g,%g,%g) # %s", shapeName, pos[0],
		pos[nmax], r[0], r[1], angle[0], text);
    } else if (nr == 1 && nangle == 0) {
	sprintf(buf, "physical;%s(%g,%g,%g) # %s", shapeName, pos[0],
		pos[nmax], r[0], text);
    } else if (nr == 0 && nangle == 1 && np == 2 ) {
	sprintf(buf, "physical;%s(%g,%g,%g,%g) # %s", shapeName, pos[0],
		pos[2], pos[nmax], pos[3], text);
    } else if (nr == 0 && nangle == 0 && np == 0 ) {
        sprintf(buf, "physical;%s() # %s", shapeName, text);
    } else if ( np == 1 && nr == 0 && nangle == 0 ) {
        sprintf(buf, "physical;%s(%g,%g) # %s", 
		shapeName, pos[0], pos[nmax], text );
    } else if (nr == 0 && nangle == 2 && np == 1 ) { /* sector */
	sprintf(buf, "physical;%s(%g,%g,%g,%g) # %s", shapeName, pos[0],
	        pos[nmax], angle[0], angle[1], text);
    } else if (nr == 2 && nangle == 2 && np == 1 ) { /* pie */
	sprintf(buf, "physical;%s(%g,%g,%g,%g,%g,%g) # %s", shapeName, pos[0],
	        pos[nmax], r[0], r[1], angle[0], angle[1], text);
    } else {
	sprintf(buf, "Unknown");
	printf("Shape not supported\tnp=%ld\tnangle=%ld\tnr=%ld\n",np,nangle,nr);
    }
    fprintf(out, "%s\n", buf);
    free(buf);
}

void reg_areg_hdr(FILE * out, char *color)
{
    char buf[512];
    fprintf(out, "%s\n", "# Region file format: DS9 version 3.0");
    sprintf(buf,
	    "global color=%s font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0",
	    color);
    fprintf(out, "%s\n", buf);
}
