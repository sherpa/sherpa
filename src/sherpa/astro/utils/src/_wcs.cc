// 
//  Copyright (C) 2009  Smithsonian Astrophysical Observatory
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
#include <sstream>
#include <iostream>

extern "C" {
#include "wcs.h"

  void init_wcs();
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
  
  if ( ( cdelt.get_size() != crpix.get_size() ) ||
       ( crpix.get_size() != crval.get_size() ) ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "cdelt: " << cdelt.get_size()
	<< " vs crpix: " << crpix.get_size()
	<< " vs crval: " << crval.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( x0.get_size() != x1.get_size() ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"input pixel arrays x0,x1 not of equal length" );
    return NULL;
  }
  
  long nsets = long(x0.get_size());

  char    ct1[80] = "";
  char    ct2[80] = "";

  if ( !crpix || !crval || !cdelt ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"WCS params failed" );
    return NULL;
  }
  
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
			     ct1,         // Tranform flavor
			     ct2,         // Tranform flavor
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
  
  if ( ( cdelt.get_size() != crpix.get_size() ) ||
       ( crpix.get_size() != crval.get_size() ) ) {
    std::ostringstream err;
    err << "input array sizes do not match, "
	<< "cdelt: " << cdelt.get_size()
	<< " vs crpix: " << crpix.get_size()
	<< " vs crval: " << crval.get_size();
    PyErr_SetString( PyExc_TypeError, err.str().c_str() );
    return NULL;
  }

  if ( x0.get_size() != x1.get_size() ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"input pixel arrays x0,x1 not of equal length" );
    return NULL;
  }
  
  long nsets = long(x0.get_size());

  char    ct1[80] = "";
  char    ct2[80] = "";

  if ( !crpix || !crval || !cdelt ) {
    PyErr_SetString( PyExc_TypeError,
		     (char*)"WCS params failed" );
    return NULL;
  }

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
			     ct1,         // Tranform flavor
			     ct2,         // Tranform flavor
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
  
  int scale;
  for (int ii = 0; ii < nsets; ii++)
  {
    // Convert to world coords 
    wcs2pix(wcs, x0[ii], x1[ii], &x0[ii], &x1[ii], &scale);
  }    

  wcsfree(wcs);
  
  return Py_BuildValue( (char*)"(NN)",
			x0.return_new_ref(),
			x1.return_new_ref() );
}


static PyMethodDef WcsFcts[] = {
  { (char*)"pix2world", (PyCFunction)pix2world,
    METH_VARARGS, (char*)"Pixel to World Coordinates" },
  { (char*)"world2pix", (PyCFunction)world2pix,
    METH_VARARGS, (char*)"World to Pixel Coordinates" },
  
  { NULL, NULL, 0, NULL }

};

SHERPAMOD(_wcs, WcsFcts)
