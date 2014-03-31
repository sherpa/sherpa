/*** File wcslib/zpxpos.c
 *** October 31, 2012
 *** By Frank Valdes, valdes@noao.edu
 *** Modified from tnxpos.c by Jessica Mink, jmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** After IRAF mwcs/wfzpx.x
 *** Copyright (C) 1998-2012
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: jmink@cfa.harvard.edu
           Postal address: Jessica Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "wcs.h"

#define	TOL 1e-13
#define SPHTOL 0.00001
#define BADCVAL 0.0
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* wfzpx -- wcs function driver for the zenithal / azimuthal polynomial.
 *    zpxinit (header, wcs)
 *    zpxclose (wcs)
 *    zpxfwd (xpix, ypix, wcs, xpos, ypos)	Pixels to WCS
 *    zpxrev (xpos, ypos, wcs, xpix, ypix)	WCS to pixels
 */

#define	max_niter	500
#define	SZ_ATSTRING	2000
static void wf_gsclose();

/* zpxinit -- initialize the zenithal/azimuthal polynomial forward or
 * inverse transform. initialization for this transformation consists of,
 * determining which axis is ra / lon and which is dec / lat, computing the
 * celestial longitude and colatitude of the native pole, reading in the the
 * native longitude of the pole of the celestial coordinate system longpole
 * from the attribute list, precomputing the euler angles and various
 * intermediary functions of the reference coordinates, reading in the
 * projection parameter ro from the attribute list, reading in up to ten
 * polynomial coefficients, and, for polynomial orders greater than 2 computing
 * the colatitude and radius of the  first point of inflection. if longpole is
 * undefined then a value of 180.0 degrees is assumed. if ro is undefined a
 * value of 180.0 / pi is assumed. if the polynomial coefficients are all zero
 * then an error condition is posted. if the order of the polynomial is 2 or
 * greater and there is no point of inflection an error condition is posted.
 * the zpx projection with an order of 1 and 0th and 1st coefficients of 0.0
 * and 1.0 respectively is equivalent to the arc projtection. in order to
 * determine the axis order, the parameter "axtype={ra|dec} {xlon|xlat}" must
 * have been set in the attribute list for the function. the longpole and ro
 * parameters may be set in either or both of the axes attribute lists, but the
 * value in the ra axis attribute list takes precedence.
 */

int
zpxinit (header, wcs)

const char *header;	/* FITS header */
struct WorldCoor *wcs;	/* pointer to WCS structure */
{
    int i, j;
    struct IRAFsurface *wf_gsopen();
    char key[8], *str1, *str2, *lngstr, *latstr, *header1;
    double zd1, d1, zd2,d2, zd, d, r;
    extern void wcsrotset();

    /* allocate space for the attribute strings */
    str1 = malloc (SZ_ATSTRING);
    str2 = malloc (SZ_ATSTRING);
    if (!hgetm (header, "WAT1", SZ_ATSTRING, str1)) {
	/*  this is a kludge to handle NOAO archived data where the first
	 *  WAT cards are in the primary header and this code does not
	 *  implement the inheritance convention.  since zpx is largely an
	 *  NOAO system and it doesn't make sense for WAT1 to be missing if
	 *  ctype is ZPX, this block is only triggered with this kludge.
	 *  there had to be a few changes to defeat the caching of the
	 *  index of the header string so that the added cards are also
	 *  found.
	 */
	
	header1 = malloc (strlen(header)+200);
        strcpy (header1, "WAT1_001= 'wtype=zpx axtype=ra projp0=0. projp1=1. projp2=0. projp3=337.74 proj'WAT2_001= 'wtype=zpx axtype=dec projp0=0. projp1=1. projp2=0. projp3=337.74 pro'");
	strcat (header1, header);
	hgetm (header1, "WAT1", SZ_ATSTRING, str1);
	hgetm (header1, "WAT2", SZ_ATSTRING, str2);
	free (header1);
    }
    hgetm (header, "WAT2", SZ_ATSTRING, str2);

    lngstr = malloc (SZ_ATSTRING);
    latstr = malloc (SZ_ATSTRING);

    /* determine the native longitude of the pole of the celestial
	coordinate system corresponding to the FITS keyword longpole.
	this number has no default and should normally be set to 180
	degrees. search both axes for this quantity. */

    if (wcs->longpole > 360.0) {
	if (!igetr8 (str1, "longpole", &wcs->longpole)) {
	    if (!igetr8 (str2, "longpole", &wcs->longpole))
		wcs->longpole = 180.0;
	    }
	}

    /*  Fetch the ro projection parameter which is the radius of the
	generating sphere for the projection. if ro is absent which
	is the usual case set it to 180 / pi. search both axes for
	this quantity. */

    if (!igetr8 (str1, "ro", &wcs->rodeg)) {
	if (!igetr8 (str2, "ro", &wcs->rodeg))
	    wcs->rodeg = 180.0 / PI;
	}

    /* Fetch the zenithal polynomial coefficients. */
    for (i = 0; i < 10; i++) {
	sprintf (key,"projp%d",i);
	if (!igetr8 (str1, key, &wcs->prj.p[i]))
	    wcs->prj.p[i] = 0.0;
    }

    /*  Fetch the longitude correction surface. note that the attribute
	string may be of any length so the length of atvalue may have
	to be adjusted. */

    if (!igets (str1, "lngcor", SZ_ATSTRING, lngstr)) {
	if (!igets (str2, "lngcor", SZ_ATSTRING, lngstr))
	    wcs->lngcor = NULL;
	else
	    wcs->lngcor = wf_gsopen (lngstr);
	}
    else
	wcs->lngcor = wf_gsopen (lngstr);

    /*  Fetch the latitude correction surface. note that the attribute
	string may be of any length so the length of atvalue may have
	to be adjusted. */

    if (!igets (str2, "latcor", SZ_ATSTRING, latstr)) {
	if (!igets (str1, "latcor", SZ_ATSTRING, latstr))
	    wcs->latcor = NULL;
	else
	    wcs->latcor = wf_gsopen (latstr);
	}
    else
	wcs->latcor = wf_gsopen (latstr);

    /* Determine the number of ZP coefficients */
    for (i = 9; i >= 0 && wcs->prj.p[i] == 0.; i--);
    wcs->zpnp = i;

    if (i >= 3) {
        /* Find the point of inflection closest to the pole. */
	zd1 = 0.;
	d1 = wcs->prj.p[1];

	/* Find the point where the derivative first goes negative. */
	for (i = 1; i<= 180; i++) {
	    zd2 = PI * i / 180.0;
	    d2 = 0.;
	    for (j = wcs->zpnp; j >= 1; j--) {
		d2 = d2 * zd2 + j * wcs->prj.p[j];
		}
	    if (d2 <= 0.)
		break;
	    zd1 = zd2;
	    d1 = d2;
	    }

	/* Find where the derivative is 0. */
	if (d2 <= 0.0) {
	    for (i = 1; i <= 10; i++) {
		zd = zd1 - d1 * (zd2 - zd1) / (d2 - d1);
		d = 0.;
		for (j = wcs->zpnp; j >= 1; j--) {
		    d = d * zd + j * wcs->prj.p[j];
		    }
		if (fabs(d) < TOL)
		    break;
		if (d < 0.) {
		    zd2 = zd;
		    d2 = d;
		    }
		else {
		    zd1 = zd;
		    d1 = d;
		    }
		}
	    }

	/* No negative derivative. */
	else 
	    zd = PI;

	r = 0.;
	for (j = wcs->zpnp; j >= 0; j--)
	    r = r * zd + wcs->prj.p[j];
	wcs->zpzd = zd;
	wcs->zpr = r;
	}

    /* Compute image rotation */
    wcsrotset (wcs);

    /* free working space. */
    free (str1);
    free (str2);
    free (lngstr);
    free (latstr);

    /* Return 1 if there are no correction coefficients */
    if (wcs->latcor == NULL && wcs->lngcor == NULL)
	return (1);
    else
	return (0);
}


/* zpxpos -- forward transform (physical to world) gnomonic projection. */

int
zpxpos (xpix, ypix, wcs, xpos, ypos)

double	xpix, ypix;	/*i physical coordinates (x, y) */
struct WorldCoor *wcs;	/*i pointer to WCS descriptor */
double	*xpos, *ypos;	/*o world coordinates (ra, dec) */
{
    int	i, j, k, ira, idec;
    double x, y, r, phi, theta, costhe, sinthe, dphi, cosphi, sinphi, dlng, z;
    double colatp, coslatp, sinlatp, longp;
    double xs, ys, ra, dec, xp, yp;
    double a, b, c, d, zd, zd1, zd2, r1, r2, rt, lambda;
    double wf_gseval();

    /* Convert from pixels to image coordinates */
    xpix = xpix - wcs->crpix[0];
    ypix = ypix - wcs->crpix[1];

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	x = xpix * wcs->cd[0] + ypix * wcs->cd[1];
	y = xpix * wcs->cd[2] + ypix * wcs->cd[3];
	}

    else {

	/* Check axis increments - bail out if either 0 */
	if (wcs->cdelt[0] == 0.0 || wcs->cdelt[1] == 0.0) {
	    *xpos = 0.0;
	    *ypos = 0.0;
	    return 2;
	    }

	/* Scale using CDELT */
	xs = xpix * wcs->cdelt[0];
	ys = ypix * wcs->cdelt[1];

	/* Take out rotation from CROTA */
	if (wcs->rot != 0.0) {
	    double cosr = cos (degrad (wcs->rot));
	    double sinr = sin (degrad (wcs->rot));
	    x = xs * cosr - ys * sinr;
	    y = xs * sinr + ys * cosr;
    	    }
	else {
	    x = xs;
	    y = ys;
	    }
	}

    /* Get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}
    colatp = degrad (90.0 - wcs->crval[idec]);
    coslatp = cos(colatp);
    sinlatp = sin(colatp);
    longp = degrad(wcs->longpole);

    /*  Compute native spherical coordinates phi and theta in degrees from the
	projected coordinates. this is the projection part of the computation */
    k = wcs->zpnp;
    if (wcs->lngcor != NULL)
	xp = x + wf_gseval (wcs->lngcor, x, y);
    else
	xp = x;
    if (wcs->latcor != NULL)
	yp = y + wf_gseval (wcs->latcor, x, y);
    else
	yp = y;
    x = xp;
    y = yp;
    r = sqrt (x * x + y * y) / wcs->rodeg;

    /* Solve */

    /* Constant no solution */
    if (k < 1) {
        *xpos = BADCVAL;
        *ypos = BADCVAL;
	return (1);
	}

    /* Linear */
    else if (k == 1) {
        zd = (r - wcs->prj.p[0]) / wcs->prj.p[1];
	}

    /* Quadratic */
    else if (k == 2) {

        a = wcs->prj.p[2];
        b = wcs->prj.p[1];
        c = wcs->prj.p[0] - r;
	d = b * b - 4. * a * c;
	if (d < 0.) {
	    *xpos = BADCVAL;
	    *ypos = BADCVAL;
	    return (1);
	    }
	d = sqrt (d);

	/* Choose solution closest to the pole */
	zd1 = (-b + d) / (2. * a);
	zd2 = (-b - d) / (2. * a);
	if (zd1 < zd2)
	    zd = zd1;
	else
	    zd = zd2;
	if (zd < -TOL) {
	    if (zd1 > zd2)
		zd = zd1;
	    else
		zd = zd2;
	    }
	if (zd < 0.) {
	    if (zd < -TOL) {
		*xpos = BADCVAL;
		*ypos = BADCVAL;
		return (1);
		}
	    zd = 0.;
	    }
	else if (zd > PI) {
	    if (zd > (PI + TOL)) {
		*xpos = BADCVAL;
		*ypos = BADCVAL;
		return (1);
		}
	    zd = PI;
	    }
	}

    /* Higher order solve iteratively */
    else {

        zd1 = 0.;
	r1 = wcs->prj.p[0];
	zd2 = wcs->zpzd;
	r2 = wcs->zpr;

	if (r < r1) {
	    if (r < (r1 - TOL)) {
		*xpos = BADCVAL;
		*ypos = BADCVAL;
		return (1);
		}
	    zd = zd1;
	    }
	else if (r > r2) {
	    if (r > (r2 + TOL)) {
		*xpos = BADCVAL;
		*ypos = BADCVAL;
		return (1);
		}
	    zd = zd2;
	    }
	else {
	    for (j=0; j<100; j++) {
	        lambda = (r2 - r) / (r2 - r1);
		if (lambda < 0.1)
		    lambda = 0.1;
		else if (lambda > 0.9)
		    lambda = 0.9;
		zd = zd2 - lambda * (zd2 - zd1);
		rt = 0.;
		for (i=k; i>=0; i--)
		    rt = (rt * zd) + wcs->prj.p[i];
		if (rt < r) {
		    if ((r - rt) < TOL)
		        break;
		    r1 = rt;
		    zd1 = zd;
		    }
		else {
		    if ((rt - r) < TOL)
		        break;
		    r2 = rt;
		    zd2 = zd;
		    }
		lambda = zd2 - zd1;
		lambda = fabs (zd2 - zd1);
		if (fabs (zd2 - zd1) < TOL)
		    break;
		}
	    }
	}

    /* Compute phi */
    if (r == 0.0)
	phi = 0.0;
    else
	phi = atan2 (x, -y);

    /* Compute theta */
    theta = PI / 2 - zd;

    /*  Compute the celestial coordinates ra and dec from the native
	coordinates phi and theta. this is the spherical geometry part
	of the computation */

    costhe = cos (theta);
    sinthe = sin (theta);
    dphi = phi - longp;
    cosphi = cos (dphi);
    sinphi = sin (dphi);

    /* Compute the ra */
    x = sinthe * sinlatp - costhe * coslatp * cosphi;
    if (fabs (x) < SPHTOL)
	x = -cos (theta + colatp) + costhe * coslatp * (1.0 - cosphi);
    y = -costhe * sinphi;
    if (x != 0.0 || y != 0.0)
	dlng = atan2 (y, x);
    else
	dlng = dphi + PI ;
    ra =  wcs->crval[ira] + raddeg(dlng);

    /* normalize ra */
    if (wcs->crval[ira] >= 0.0) {
	if (ra < 0.0)
	    ra = ra + 360.0;
	}
    else {
	if (ra > 0.0)
	    ra = ra - 360.0;
	}
    if (ra > 360.0)
	ra = ra - 360.0;
    else if (ra < -360.0)
	ra = ra + 360.0;

    /* compute the dec */
    if (fmod (dphi, PI) == 0.0) {
	dec = raddeg(theta + cosphi * colatp);
	if (dec > 90.0)
	    dec = 180.0 - dec;
	if (dec < -90.0)
	    dec = -180.0 - dec;
	}
    else {
	z = sinthe * coslatp + costhe * sinlatp * cosphi;
	if (fabs(z) > 0.99) {
	    if (z >= 0.0)
		dec = raddeg(acos (sqrt(x * x + y * y)));
	    else
		dec = raddeg(-acos (sqrt(x * x + y * y)));
	    }
	else
		dec = raddeg(asin (z));
	}

    /* store the results */
    *xpos  = ra;
    *ypos = dec;
    return (0);
}


/* zpxpix -- inverse transform (world to physical) for the zenithal
 * azimuthal polynomial projection.
 */

int
zpxpix (xpos, ypos, wcs, xpix, ypix)

double	xpos, ypos;	/*i world coordinates (ra, dec) */
struct WorldCoor *wcs;	/*i pointer to WCS descriptor */
double	*xpix, *ypix;	/*o physical coordinates (x, y) */
{
    int	i, ira, idec, niter;
    double ra, dec, cosdec, sindec, cosra, sinra, x, y, phi, theta;
    double s, r, dphi, z, dpi, dhalfpi, twopi, tx;
    double xm, ym, f, fx, fy, g, gx, gy, denom, dx, dy;
    double colatp, coslatp, sinlatp, longp, sphtol;
    double wf_gseval(), wf_gsder();

    /* get the axis numbers */
    if (wcs->coorflip) {
	ira = 1;
	idec = 0;
	}
    else {
	ira = 0;
	idec = 1;
	}

    /*  Compute the transformation from celestial coordinates ra and
	dec to native coordinates phi and theta. this is the spherical
	geometry part of the transformation */

    ra  = degrad (xpos - wcs->crval[ira]);
    dec = degrad (ypos);
    cosra = cos (ra);
    sinra = sin (ra);
    cosdec = cos (dec);
    sindec = sin (dec);
    colatp = degrad (90.0 - wcs->crval[idec]);
    coslatp = cos (colatp);
    sinlatp = sin (colatp);
    if (wcs->longpole == 999.0)
	longp = degrad (180.0);
    else
	longp = degrad(wcs->longpole);
    dpi = PI;
    dhalfpi = dpi * 0.5;
    twopi = PI + PI;
    sphtol = SPHTOL;

    /* Compute phi */
    x = sindec * sinlatp - cosdec * coslatp * cosra;
    if (fabs(x) < sphtol)
	x = -cos (dec + colatp) + cosdec * coslatp * (1.0 - cosra);
    y = -cosdec * sinra;
    if (x != 0.0 || y != 0.0)
	dphi = atan2 (y, x);
    else
	dphi = ra - dpi;
    phi = longp + dphi;
    if (phi > dpi)
	phi = phi - twopi;
    else if (phi < -dpi)
	phi = phi + twopi;

    /* Compute theta */
    if (fmod (ra, dpi) == 0.0) {
	theta = dec + cosra * colatp;
	if (theta > dhalfpi)
	    theta = dpi - theta;
	if (theta < -dhalfpi)
	    theta = -dpi - theta;
	}
    else {
	z = sindec * coslatp + cosdec * sinlatp * cosra;
	if (fabs (z) > 0.99) {
	    if (z >= 0.0)
		theta = acos (sqrt(x * x + y * y));
	    else
		theta = -acos (sqrt(x * x + y * y));
	    }
	else
	    theta = asin (z);
	}

    /*  Compute the transformation from native coordinates phi and theta
	to projected coordinates x and y */

    s = dhalfpi - theta;
    r = 0.;
    for (i=9; i>=0; i--)
        r = r * s + wcs->prj.p[i];
    r = wcs->rodeg * r;

    if (wcs->lngcor == NULL && wcs->latcor == NULL) {
	if (wcs->coorflip) {
	    y  = r * sin (phi);
	    x = -r * cos (phi);
	} else {
	    x  = r * sin (phi);
	    y = -r * cos (phi);
	}
    } else {
	xm  = r * sin (phi);
	ym = -r * cos (phi);
	x = xm;
	y = ym;
	niter = 0;
	while (niter < max_niter) {
	    if (wcs->lngcor != NULL) {
		f = x + wf_gseval (wcs->lngcor, x, y) - xm;
		fx = wf_gsder (wcs->lngcor, x, y, 1, 0);
		fx = 1.0 + fx;
		fy = wf_gsder (wcs->lngcor, x, y, 0, 1);
		}
	    else {
		f = x - xm;
		fx = 1.0 ;
		fy = 0.0;
		}
	    if (wcs->latcor != NULL) {
		g = y + wf_gseval (wcs->latcor, x, y) - ym;
		gx = wf_gsder (wcs->latcor, x, y, 1, 0);
		gy = wf_gsder (wcs->latcor, x, y, 0, 1);
		gy = 1.0 + gy;
		}
	    else {
		g = y - ym;
		gx = 0.0 ;
		gy = 1.0;
		}

	    denom = fx * gy - fy * gx;
	    if (denom == 0.0)
		break;
	    dx = (-f * gy + g * fy) / denom;
	    dy = (-g * fx + f * gx) / denom;
	    x = x + dx;
	    y = y + dy;
	    if (MAX(MAX(fabs(dx),fabs(dy)),MAX(fabs(f),fabs(g))) < 2.80e-8)
		break;

	    niter = niter + 1;
	}

	/* Reverse x and y if axes flipped */
	if (wcs->coorflip) {
	    tx = x;
	    x = y;
	    y = tx;
	}
    }

    /* Scale and rotate using CD matrix */
    if (wcs->rotmat) {
	*xpix = x * wcs->dc[0] + y * wcs->dc[1];
	*ypix = x * wcs->dc[2] + y * wcs->dc[3];
	}

    else {

	/* Correct for rotation */
	if (wcs->rot!=0.0) {
	    double cosr = cos (degrad (wcs->rot));
	    double sinr = sin (degrad (wcs->rot));
	    *xpix = x * cosr + y * sinr;
	    *ypix = y * cosr - x * sinr;
	    }
	else {
	    *xpix = x;
	    *ypix = y;
	    }

	/* Scale using CDELT */
	if (wcs->xinc != 0.)
	    *xpix = *xpix / wcs->xinc;
	if (wcs->yinc != 0.)
	    *ypix = *ypix / wcs->yinc;
	}

    /* Convert to pixels  */
    *xpix = *xpix + wcs->xrefpix;
    *ypix = *ypix + wcs->yrefpix;

    return (0);
}


/* ZPXCLOSE -- free up the distortion surface pointers */

void
zpxclose (wcs)

struct WorldCoor *wcs;		/* pointer to the WCS descriptor */

{
    if (wcs->lngcor != NULL)
	wf_gsclose (wcs->lngcor);
    if (wcs->latcor != NULL)
	wf_gsclose (wcs->latcor);
    return;
}


/* wf_gsclose -- procedure to free the surface descriptor */

static void
wf_gsclose (sf)

struct IRAFsurface *sf;	/* the surface descriptor */

{
    if (sf != NULL) {
	if (sf->xbasis != NULL)
	    free (sf->xbasis);
	if (sf->ybasis != NULL)
	    free (sf->ybasis);
	if (sf->coeff != NULL)
	    free (sf->coeff);
	free (sf);
	}
    return;
}

/*
 * Mar  8 2011  Created from tnxpos.c and wfzpx.x
 *
 * Oct 31 2012	End comment on line 346 after pole; fix code thereafter
 */
