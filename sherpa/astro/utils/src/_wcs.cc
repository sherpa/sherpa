//
//  Copyright (C) 2009, 2025-2026
//  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//


#include "sherpa/extension.hh"

extern "C" {
#include "wcs.h"
}

static PyObject* pix2world( PyObject* self, PyObject* args )
{
  DoubleArray x0;
  DoubleArray x1;
  DoubleArray crpix;
  DoubleArray crval;
  DoubleArray cdelt;
  double crota = 0.0, epoch = 0.0, equinox = 0.0;
  const char *ctype = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"sO&O&O&O&O&ddd",
			  &ctype,
			  CONVERTME(DoubleArray), &x0,
			  CONVERTME(DoubleArray), &x1,
			  CONVERTME(DoubleArray), &crpix,
			  CONVERTME(DoubleArray), &crval,
			  CONVERTME(DoubleArray), &cdelt,
			  &crota,
			  &equinox,
			  &epoch ) )
    return NULL;

  if (crpix.get_size() != 2 ||
      crval.get_size() != 2 ||
      cdelt.get_size() != 2) {
    return PyErr_Format( PyExc_TypeError,
			 (char*)"crpix, crval, cdelt must be 2-element arrays" );
  }

  if ( x0.get_size() != x1.get_size() ) {
    return PyErr_Format( PyExc_TypeError,
			 (char*)"input pixel arrays x0,x1 not of equal length" );
  }

  long nsets = long(x0.get_size());

  char    ct1[80] = "";
  char    ct2[80] = "";

  if ( strcmp(ctype,"LINEAR") == 0 )
  {
    strcpy(ct1, "LINEAR");
    strcpy(ct2, "LINEAR");
  }
  else if ( strcmp(ctype,"TAN-P") == 0 )
  {
    strcpy(ct1, "LONG-TAN");
    strcpy(ct2, "NPOL-TAN");
  }
  else
  {
    strcpy(ct1, "RA---TAN");
    strcpy(ct2, "DEC--TAN");
  }


  // Initialize WCS
  WorldCoor* wcs = wcskinit( 2,           // number of values per set
			     nsets,       // number of value sets.
			     ct1,         // Transform flavor
			     ct2,         // Transform flavor
			     crpix[0],
			     crpix[1],
			     crval[0],
			     crval[1],
			     NULL,
			     cdelt[0],
			     cdelt[1],
			     crota,
			     int(equinox),
			     epoch
			     );
  if (!wcs) {
    return PyErr_Format( PyExc_RuntimeError,
			 (char*) "Unable to initialze the WCS structure" );
  }

  // ctype is expected to be one of: "WCS", "LINEAR", or "TAN-P".
  //
  if ( strcmp(ctype,"LINEAR") == 0 )
  {
    // The above will initialize wcs->sysout to "FK5" which will then cause
    // pix2wcs() to bound check the xpos output value and adjust it by +/-
    // 360.0 if the value is either positive or negative.  To avoid this,
    // the sysout value must be 'WCS_LINEAR'
    // This behaviour change occurs for WCSSUBS-3.7.0.. probably since 3.5.7
    // I do not see an init routine that will not do this..
    wcs->sysout = WCS_LINEAR;
  }
  else if ( strcmp(ctype,"TAN-P") == 0 )
  {
    // LONGPOLE and LATPOLE should be populated from the header keywords
    // LONP# and LATP#, where # is associated with the transform.
    // Hardcode for now, but should modify to do better (add parameters?)
    wcs->longpole = 270.;
    wcs->latpole  =   0.;
    if ( wcs->longpole != 0.0 )
      wcs->cel.ref[2] = wcs->longpole;
    if ( wcs->latpole != 0.0 )
      wcs->cel.ref[3] = wcs->latpole;
  }

  for (int ii = 0; ii < nsets; ii++)
  {
    // Convert to world coords
    pix2wcs(wcs, x0[ii], x1[ii], &x0[ii], &x1[ii]);
  }

  wcsfree(wcs);

  return Py_BuildValue( (char*)"(NN)",
			x0.return_new_ref(),
			x1.return_new_ref() );
}


static PyObject* world2pix( PyObject* self, PyObject* args )
{
  DoubleArray x0;
  DoubleArray x1;
  DoubleArray crpix;
  DoubleArray crval;
  DoubleArray cdelt;
  double crota = 0.0, epoch = 0.0, equinox = 0.0;
  const char *ctype = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"sO&O&O&O&O&ddd",
			  &ctype,
			  CONVERTME(DoubleArray), &x0,
			  CONVERTME(DoubleArray), &x1,
			  CONVERTME(DoubleArray), &crpix,
			  CONVERTME(DoubleArray), &crval,
			  CONVERTME(DoubleArray), &cdelt,
			  &crota,
			  &equinox,
			  &epoch ) )
    return NULL;

  if (crpix.get_size() != 2 ||
      crval.get_size() != 2 ||
      cdelt.get_size() != 2) {
    return PyErr_Format( PyExc_TypeError,
			 (char*)"crpix, crval, cdelt must be 2-element arrays" );
  }

  if ( x0.get_size() != x1.get_size() ) {
    return PyErr_Format( PyExc_TypeError,
			 (char*)"input pixel arrays x0,x1 not of equal length" );
  }

  long nsets = long(x0.get_size());

  char    ct1[80] = "";
  char    ct2[80] = "";

  // ctype is expected to be one of: "WCS", "LINEAR", or "TAN-P".
  //
  if ( strcmp(ctype,"LINEAR") == 0 )
  {
    strcpy(ct1, "LINEAR");
    strcpy(ct2, "LINEAR");
  }
  else
  {
    strcpy(ct1, "RA---TAN");
    strcpy(ct2, "DEC--TAN");
  }


  // Initialize WCS
  WorldCoor* wcs = wcskinit( 2,           // number of values per set
			     nsets,       // number of value sets.
			     ct1,         // Transform flavor
			     ct2,         // Transform flavor
			     crpix[0],
			     crpix[1],
			     crval[0],
			     crval[1],
			     NULL,
			     cdelt[0],
			     cdelt[1],
			     crota,
			     int(equinox),
			     epoch
			     );
  if (!wcs) {
    return PyErr_Format( PyExc_RuntimeError,
			 (char*) "Unable to initialze the WCS structure" );
  }

  int scale;
  for (int ii = 0; ii < nsets; ii++)
  {
    // Convert to world coords; do not report if scale is not 0 (i.e.
    // the conversion is "out of bounds").
    wcs2pix(wcs, x0[ii], x1[ii], &x0[ii], &x1[ii], &scale);
  }

  wcsfree(wcs);

  return Py_BuildValue( (char*)"(NN)",
			x0.return_new_ref(),
			x1.return_new_ref() );
}


static PyMethodDef WcsFcts[] = {
  FCTSPECDOC(pix2world, pix2world,
	     "Pixel to World Coordinates\n\n"
	     "Parameters\n----------\n"
	     "ctype : {'LINEAR', 'WCS', 'WCS-P'}\n"
	     "    The coordinate transformation ('WCS' is a tangent-plane\n"
	     "    projection).\n"
	     "x0, x1 : array\n"
	     "    The coordinates to convert (each a 1D array of n points).\n"
	     "crpix, crval, cdelt : array\n"
	     "    These are 2-element arrays.\n"
	     "crota : float\n"
	     "equinox : int\n"
	     "epoch : float\n\n"
	     "Returns\n-------\n"
	     "w0, w1 : ndarray\n\n"
	     "Notes\n-----\n"
	     "The x0 and x1 arrays will be updated to match w0, w1 if they\n"
	     "are ndarrays of real numbers. This behaviour should not be relied on.\n\n"
	     "Examples\n--------\n"
	     ">>> x0 = [20.5, 5.2, 3.4]\n"
	     ">>> x1 = [10.5, 7.0, 15.8]\n"
	     ">>> ctype = 'LINEAR'\n"
	     ">>> crpix = [20.5, 10.5]\n"
	     ">>> crval = [400, -200]\n"
	     ">>> cdelt = [-4.0, 2.0]\n"
	     ">>> crota = 0    # unused with LINEAR\n"
	     ">>> equinox = 0  # unused with LINEAR\n"
	     ">>> epoch = 0    # unused with LINEAR\n"
	     ">>> w0, w1 = pix2world(ctype, x0, x1, crpix, crval, cdelt, crota, equinox, epoch)\n"
	     ">>> print(w0)\n"
	     "[400.  461.2 468.4]\n"
	     ">>> print(w1)\n"
	     "[-200.  -207.  -189.4]\n"

	     ),
  FCTSPECDOC(world2pix, world2pix,
	     "World to Pixel Coordinates\n\n"
	     "Parameters\n----------\n"
	     "ctype : {'LINEAR', 'WCS', 'WCS-P'}\n"
	     "    The coordinate transformation ('WCS' is a tangent-plane\n"
	     "    projection).\n"
	     "w0, w1 : array\n"
	     "    The coordinates to convert (each a 1D array of n points).\n"
	     "crpix, crval, cdelt : array\n"
	     "    These are 2-element arrays.\n"
	     "crota : float\n"
	     "equinox : int\n"
	     "epoch : float\n\n"
	     "Returns\n-------\n"
	     "x0, x1 : ndarray\n\n"
	     "Notes\n-----\n"
	     "The w0 and w1 arrays will be updated to match w0, w1 if they\n"
	     "are ndarrays of real numbers. This behaviour should not be relied on.\n\n"
	     "Examples\n--------\n"
	     ">>> w0 = [400.0, 461.2, 468.4]\n"
	     ">>> w1 = [-200.0, -207, -189.4]\n"
	     ">>> ctype = 'LINEAR'\n"
	     ">>> crpix = [20.5, 10.5]\n"
	     ">>> crval = [400, -200]\n"
	     ">>> cdelt = [-4.0, 2.0]\n"
	     ">>> crota = 0    # unused with LINEAR\n"
	     ">>> equinox = 0  # unused with LINEAR\n"
	     ">>> epoch = 0    # unused with LINEAR\n"
	     ">>> x0, x1 = world2pix(ctype, w0, w1, crpix, crval, cdelt, crota, equinox, epoch)\n"
	     ">>> print(x0)\n"
	     "[20.5  5.2  3.4]\n"
	     ">>> print(x1)\n"
	     "[10.5  7.  15.8]\n"
	     ),

  { NULL, NULL, 0, NULL }

};

SHERPAMODDOC(_wcs, WcsFcts,
	     "Low-level conversion between 2D Astronomical coordinate systems.\n\n"
	     "This module is designed to support Chandra images and has limited\n"
	     "support for other coordinate systems. These routines are not intended\n"
	     "to be called directly, but indirectly via the sherpa.astro.io.wcs.WCS\n"
	     "class.\n"

	     )
