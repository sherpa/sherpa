//  Copyright (C) 2007, 2015-2018, 2019, 2020
//      Smithsonian Astrophysical Observatory
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

// INIT_XSPEC is used by xspecmodelfct() in xspec_extension.hh, so it needs to
// be defined before that file is included
int _sherpa_init_xspec_library();
#define INIT_XSPEC _sherpa_init_xspec_library

// Have sherpa include first so that Python.h is first, to avoid warning
// messages about redefining _XOPEN_SOURCE
#include "sherpa/astro/xspec_extension.hh"
#include <iostream>
#include <fstream>

// The symbols listed in XSPEC version 12.9.1
// at https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixExternal.html
// are given below. Note that this is the C/FORTRAN interface, not the
// more-featureful FunctionUtility module.
//
// Functions which are used below:
// FNINIT	Initializes data directory locations needed by the models. See below for a fuller description.
// FGABND	Get an element abundance.
// FGCHAT	Get current chatter level setting for model functions' output verbosity.
// FPCHAT	Set the chatter level. Default is 10, higher chatter levels produce more output.
// FGMSTR	Get a model string value (see XSPEC xset command).
// FPMSTR	Set a model string value.
// FGDATD	Get the model .dat files path.
// FPDATD	Set the model .dat files path.
// FGMODF	Get the model ion data path.
// FPSLFL	Load values of a file solar abundance table (see abund command).
// FGSOLR	Get the solar abundance table setting.
// FPSOLR	Set the solar abundance table.
// FGXSCT	Get the cross section table setting.
// FPXSCT	Set the cross section table.
// csmgh0	Get the cosmology H$_0$ setting (see the cosmo command).
// csmph0	Set H$_0$.
// csmgl0	Get $\Lambda_0$.
// csmpl0	Set $\Lambda_0$.
// csmgq0	Get q$_0$.
// csmpq0	Put q$_0$.
// xs_getVersion (or xgvers)	Retrieve XSPEC's version string.
//
// Functions which are not wrapped as their functionality is available:
// RFLABD	Read abundance data from a file, then load and set this to be the current abundance table. (Essentially this combines a file read with the FPSLFL and FPSOLR functions.)
//
// Functions not wrapped as not felt to be that useful:
// fzsq	Computes the luminosity distance, (c/H$_0$)*fzsq. The function is valid for small values of q$_0$*z for the case of no cosmological constant and uses the approximation of Pen (1999 ApJS 120, 49) for the case of a cosmological constant and a flat Universe. The function is not valid for non-zero cosmological constant if the Universe is not flat.
//
// Functions not wrapped since they are not useful as is (they need
// functionality from 12.9.1 to set the XFLT keywords):
// DGFILT	Get a particular XFLT keyword value from a data file.
// DGNFLT	Get the number of XFLT keywords in a data file.
//
// Other symbols in xsFortran.h are:
// DGQFLT       Does a XFLT keyword exist?
// PDBVAL       Set a database value
//
// Symbols in 12.9.1/HEASOFT 6.22 but not in 12.9.0/HEASOFT 6.19
// FGABNZ
// FGTABN
// FGTABZ
// FGELTI
// FGNELT
// FGABFL
// FPABFL
// FGAPTH
// FPAPTH
// csmpall
// DPFILT
// DCLFLT
// GDBVAL
// CDBASE
// FGATDV
// FPATDV
//
// These seem unlikely to be useful for Sherpa
// xs_getChat
// xs_write
// xs_read
//
// These are numeric functions which we should have available elsewhere
// xs_erf
// xs_erfc
// gammap
// gammq
//
#include "xsFortran.h"

// TODO: is this defined in an XSPEC header file?
#define ABUND_SIZE (30) // number of elements in Solar Abundance table

// C_<model> are declared here; the other models are defined in
// functionMap.h but that requires using the XSPEC build location
// rather than install location.
//
#include "funcWrappers.h"

extern "C" {

#ifdef XSPEC_12_10_1
void agnsed_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void qsosed_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifdef XSPEC_12_11_0
void agnslim_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifndef XSPEC_12_9_1
void xsaped_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xsblbd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbbrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbmc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbrms_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_9_1
void xsbvpe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifndef XSPEC_12_10_0
void c6mekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6vmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void cemekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compls_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xstitg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void disk_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskir_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdskb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#ifndef XSPEC_12_10_1
void xsdili_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif
void diskm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void disko_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskpbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdiskpn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsxpdec_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void ezdiskbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_9_1
void xsgaul_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void grad_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifdef XSPEC_12_10_0
// Note: have dropped the leading 'c_' for this model
void xsgrbcomp(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

void jet_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xsgrbm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void spin_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifdef XSPEC_12_10_1
void kyrline_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifndef XSPEC_12_9_1
void xslorz_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmeka_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void nsa_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsagrav_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsatmos_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_10_0
void nsmax_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xspegp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsp1tr_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsposm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_9_1
void xsrays_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xredge_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsrefsch_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void srcut_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void sresc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#ifdef XSPEC_12_10_0
void ssa_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif
void xsstep_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_9_1
void xsvape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xsbrmv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_9_1
void xsvmek_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifndef XSPEC_12_9_1
void xsvrys_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xszbod_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbrm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_10_0
void acisabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xscnst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscycl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdust_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsedge_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsexp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xshecu_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xshrfl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsntch_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsphab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsplab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xscred_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xssmdg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsspln_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xssssi_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_10_0
void swind1_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xsred_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvphb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsabsw_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xswnab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsxirf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void mszdst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszedg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszhcu_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszabp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszphb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_10_0
void zxipcf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void xszcrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void msldst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvfe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvph_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszwnb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);


#ifndef XSPEC_12_9_1
void xsbvvp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvvap_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

void zigm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifndef XSPEC_12_10_1
void logpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif
void eplogpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagnf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void pexmon_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// additive
void xscompmag(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xscomptb(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

#ifndef XSPEC_12_10_0
void nsmaxg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsx_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

//multiplicative
void xsphei_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xslyman_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

#ifdef XSPEC_12_9_1
void ismabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
// Note: have dropped the leading 'c_' for this model
void slimbbmodel(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
#endif

#ifdef XSPEC_12_11_0
void ismdust_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void olivineabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// Note: have dropped the leading 'c_' for this model
// TODO: should look at whether we want to fix this as now have a
//       number of c_xxx functions in model.dat
void beckerwolff(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
#endif


// XSPEC table models; in XSPEC 12.10.1 these have been consolidated
// into the tabint routine, but that is only available to C++, and
// so is defined in sherpa/include/sherpa/astro/xspec_extension.hh
// In XSPEC 12.11.0 tabint is available in C scope (but is again
// defined in xspec_extension.hh as it doesn't need to be visible
// here).
//
#ifndef XSPEC_12_10_1
void xsatbl(float* ear, int ne, float* param, const char* filenm, int ifl,
	    float* photar, float* photer);
void xsmtbl(float* ear, int ne, float* param, const char* filenm, int ifl,
	    float* photar, float* photer);
#endif

// XSPEC convolution models
//

void rgsxsrc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

#ifdef XSPEC_12_10_1
void kyconv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

#ifdef XSPEC_12_11_0
void thcompf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
#endif

}

// This routine could be called when the module is being initialized,
// but this would cause XSPEC module initialization even if XSPEC
// functionality is not used. So each function/model has to ensure
// that they call _sherpa_init_xspec_library before calling any
// XSPEC routine.
//
// Sun's C++ compiler complains if this is declared static
int _sherpa_init_xspec_library()
{

  static bool init = false;

  if ( init )
    return EXIT_SUCCESS;


  if ( !getenv("HEADAS") ) {
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "HEADAS environment variable is not set" );
    return EXIT_FAILURE;
  }

  // Stream buffer for redirected stdout
  std::streambuf* cout_sbuf = NULL;
  std::ofstream fout;

  try {

    // Redirect std::cout during call to FNINIT; this causes
    // us to swallow annoying messages about having set the
    // XSPEC abundance and cross section by default.

    // We have to do it this way because XSPEC code sends message
    // to XSPEC stream buffers, that in turn talk to std::cin,
    // std::cout, and std::cerr.  When we are not given a toggle
    // to turn off output to output streams, we can instead redirect
    // the stream itself.

    // In this case, the only stream we want to redirect is stdout.

    // If creating any stream buffer fails for any reason,
    // make sure we can still use original std::cout stream buffer.

    // Note:  since we are *not* redirecting std::cerr, any error
    // messages are still seen; if, for example, HEADAS is not set,
    // user will still see error message printed to stderr (i.e., to
    // screen by default).  We still want such error messages to be
    // seen!  So do *not* redirect std::cerr.
    //
    // Unfortunately it appears that XSPEC does not use stderr for its
    // error messages, but stdout, so the following code will hide
    // any error messages from FNINIT (XSPEC versions 12.8.2 and 12.9.0).
    // Perhaps the screen output could be saved, rather than redirected
    // to /dev/null, and checked for the string '***Error: ', which
    // appears to be used in the error messages from XSPEC.

    cout_sbuf = std::cout.rdbuf();
    fout.open("/dev/null");
    if (cout_sbuf != NULL && fout.is_open())
      std::cout.rdbuf(fout.rdbuf()); // temporary redirect stdout to /dev/null

    // Initialize XSPEC model library
    FNINIT();

    // Get back original std::cout
    if (cout_sbuf != NULL) {
      std::cout.clear();
      std::cout.rdbuf(cout_sbuf);
    }
    fout.clear();
    fout.close();

    // Try to minimize model chatter for normal operation.
    FPCHAT( 0 );

    // Set cosmology initial values to XSPEC initial values
    csmph0( 70.0 );
    csmpq0( 0.0 );
    csmpl0( 0.73 );

    // Work around a XSPEC 12.10.0 issue where the atomdb version is
    // hardcoded to 3.0.7 for the models-only build but it should be
    // 3.0.9. This is fixed in the yet-to-be-released 12.10.1 (hence
    // the attempt to only restrict to 12.10.0).
    //
#if defined (XSPEC_12_10_0) && !defined (XSPEC_12_10_1)
    char atomdbVersion1210[6] = "3.0.9";
    atomdbVersion1210[5] = '\0';
    FPATDV(atomdbVersion1210);
#endif

  } catch(...) {

    // Get back original std::cout
    if (cout_sbuf != NULL) {
      std::cout.clear();
      std::cout.rdbuf(cout_sbuf);
    }
    fout.clear();
    fout.close();

    // Raise appropriate error message that XSPEC initialization failed.
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "check HEADAS environment variable" );
    return EXIT_FAILURE;

  }

  init = true;

  return EXIT_SUCCESS;

}

static PyObject* get_version( PyObject *self )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char version[256];
  int retval;

  try {

    retval = xs_getVersion(version, 256);

  } catch(...) {

    PyErr_SetString( PyExc_LookupError,
		     (char*)"could not get XSPEC version string" );
    return NULL;

  }

  if( retval < 0 ) {
    PyErr_SetString( PyExc_LookupError,
		     (char*)"XSPEC version string was truncated" );
    return NULL;
  }

  return Py_BuildValue( (char*)"s", version );

}

static PyObject* get_chatter( PyObject *self )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  int chatter = 0;

  try {

    chatter = FGCHAT();

  } catch(...) {

    PyErr_SetString( PyExc_LookupError,
		     (char*)"could not get XSPEC chatter level" );
    return NULL;

  }

  return Py_BuildValue( (char*)"i", chatter );

}


static PyObject* get_abund( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* abund = NULL;
  char* element = NULL;
  PyObject *retval = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"|s", &element ) )
    return NULL;

  try {

    abund = FGSOLR();

  } catch(...) {

    PyErr_SetString( PyExc_LookupError,
		     (char*)"could not get XSPEC solar abundance" );
    return NULL;

  }

  if( !element ) {

    retval = (PyObject*) Py_BuildValue( (char*)"s", abund );

  } else {

    float abundVal = 0.0;
    std::streambuf *cerr_sbuf = NULL;
    std::ostringstream fout;

    try {

      cerr_sbuf = std::cerr.rdbuf();

      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(fout.rdbuf());

      abundVal = FGABND(element);

      // Get back original std::cerr
      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(cerr_sbuf);


    } catch(...) {

      // Get back original std::cerr
      if (cerr_sbuf != NULL)
	std::cerr.rdbuf(cerr_sbuf);

      PyErr_Format( PyExc_ValueError,
		    (char*)"could not get XSPEC abundance for '%s'",
		  element);
      return NULL;

    }

    if( fout.str().size() > 0 ) {
      PyErr_Format( PyExc_TypeError,
		    (char*)"could not find element '%s'", element);
      return NULL;
    }

    retval = (PyObject*) Py_BuildValue( (char*)"f", abundVal );
  }

  return retval;
}


static PyObject* get_cosmo( PyObject *self )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  float h0;
  float l0;
  float q0;

  try {

    h0 = csmgh0();
    l0 = csmgl0();
    q0 = csmgq0();

  } catch(...) {

    PyErr_SetString( PyExc_LookupError,
		     (char*)"could not get XSPEC cosmology settings" );
    return NULL;

  }

  return Py_BuildValue( (char*)"fff", h0, q0, l0 );

}


static PyObject* get_cross( PyObject *self )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* cross = NULL;

  try {

    cross = FGXSCT();

  } catch(...) {

    PyErr_SetString( PyExc_LookupError,
		     (char*)"could not get XSPEC cross-section" );
    return NULL;

  }

  return Py_BuildValue( (char*)"s", cross );

}


static PyObject* set_chatter( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  int chatter = 0;

  if ( !PyArg_ParseTuple( args, (char*)"i", &chatter ) )
    return NULL;

  try {

    FPCHAT( chatter );

  } catch(...) {

    PyErr_Format( PyExc_ValueError,
                  (char*)"could not set XSPEC chatter level to %d",
                  chatter);
    return NULL;

  }

  Py_RETURN_NONE;

}


static PyObject* set_abund( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* table = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"s", &table ) )
    return NULL;

  try {

    FPSOLR( table, &status );

  } catch(...) {

    status = 1;

  }

  // if abundance table name fails, try it as a filename
  if( status ) {

    std::ifstream fileStream(table);
    std::vector<float> vals(ABUND_SIZE, 0);
    size_t count(0);

    try {

      float element;
      fileStream.exceptions(std::ios_base::failbit);

      while (count < ABUND_SIZE && fileStream >> element) {
	vals[count] = element;
	++count;
      }

      status = 0;
    }
    catch ( std::exception& ) {

      if( !fileStream.eof() ) {
      	PyErr_Format( PyExc_ValueError,
                      (char*)"Cannot read file '%s'.  It may not exist or contains invalid data",
                      table);
	return NULL;
      }

      status = 1;
    }

    try {

      FPSOLR((char*)"file", &status);
      FPSLFL( &vals[0], ABUND_SIZE, &status );

    } catch(...) {

      status = 1;

    }
  }

  if ( 0 != status ) {
    PyErr_Format( PyExc_ValueError,
		  (char*)"could not set XSPEC abundance to %s",
                  table );
    return NULL;
  }

  Py_RETURN_NONE;

}


static PyObject* set_cosmo( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  float h0;
  float l0;
  float q0;

  if ( !PyArg_ParseTuple( args, (char*)"fff", &h0, &q0, &l0 ) )
    return NULL;

  try {

    csmph0( h0 );
    csmpl0( l0 );
    csmpq0( q0 );

  } catch(...) {

    PyErr_Format( PyExc_ValueError,
                  (char*)"could not set XSPEC cosmology settings to "
                  "H0=%g q0=%g Lambda0=%g",
                  h0, l0, q0);
    return NULL;

  }

  Py_RETURN_NONE;

}


static PyObject* set_cross( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* csection = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"s", &csection ) )
    return NULL;

  try {

    FPXSCT( csection, &status );

  } catch(...) {

    status = 1;

  }

  if ( 0 != status ) {
    PyErr_Format( PyExc_ValueError,
                  (char*)"could not set XSPEC photoelectric "
                  "cross-section to '%s'",
                  csection);
    return NULL;
  }

  Py_RETURN_NONE;

}


static PyObject* set_xset( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* str_name = NULL;
  char* str_value = NULL;
  int status = 0;

  if ( !PyArg_ParseTuple( args, (char*)"ss", &str_name, &str_value ) )
    return NULL;

  try {

    FPMSTR( str_name, str_value );

  } catch(...) {

    status = 1;

  }

  if ( 0 != status ) {
    PyErr_Format( PyExc_ValueError,
		  (char*)"could not set XSPEC model strings '%s: %s'",
		  str_name, str_value);
    return NULL;
  }

  Py_RETURN_NONE;

}

static PyObject* get_xset( PyObject *self, PyObject *args  )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* str_name = NULL;
  char* str_value = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s", &str_name ) )
    return NULL;

  try {

    str_value = FGMSTR( str_name );

  } catch(...) {

    PyErr_Format( PyExc_KeyError,
		  (char*)"could not get XSPEC model string '%s'",
		  str_name);
    return NULL;

  }

  return Py_BuildValue( (char*)"s", str_value );

}

// Perhaps this should be expanded to some of the other routines
// that return a string, rather than just the paths?
//
static PyObject* get_xspec_path( const char *label, char *getfunc() )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* str_value = NULL;
  try {
    str_value = getfunc();
  } catch(...) {

    std::ostringstream emsg;
    emsg << "could not get XSPEC " << label << " path";
    PyErr_SetString( PyExc_LookupError,
                     emsg.str().c_str() );
    return NULL;

  }

  return Py_BuildValue( (char*)"s", str_value );

}

static PyObject* get_manager_data_path( PyObject *self )
{
  return get_xspec_path("manager", FGDATD);
}

static PyObject* get_model_data_path( PyObject *self )
{
  return get_xspec_path("model", FGMODF);
}

static PyObject* set_manager_data_path( PyObject *self, PyObject *args )
{

  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  char* path = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s", &path ) )
    return NULL;

  try {

    FPDATD( path );

  } catch(...) {

    std::ostringstream emsg;
    emsg << "could not set XSPEC manager path to '" << path << "'";
    PyErr_SetString( PyExc_ValueError,
                     emsg.str().c_str() );
    return NULL;
  }

  Py_RETURN_NONE;

}

static PyMethodDef XSpecMethods[] = {
  { (char*)"get_xsversion", (PyCFunction)get_version, METH_NOARGS, NULL },
  { (char*)"get_xschatter", (PyCFunction)get_chatter, METH_NOARGS, NULL },
  FCTSPEC(set_xschatter, set_chatter),
  FCTSPEC(get_xsabund, get_abund),
  FCTSPEC(set_xsabund, set_abund),
  FCTSPEC(set_xscosmo, set_cosmo),
  { (char*)"get_xscosmo", (PyCFunction)get_cosmo, METH_NOARGS, NULL },
  { (char*)"get_xsxsect", (PyCFunction)get_cross, METH_NOARGS, NULL },
  FCTSPEC(set_xsxsect, set_cross),
  FCTSPEC(set_xsxset, set_xset),
  FCTSPEC(get_xsxset, get_xset),
  { (char*)"get_xspath_manager",
    (PyCFunction)get_manager_data_path, METH_NOARGS, NULL },
  { (char*)"get_xspath_model",
    (PyCFunction)get_model_data_path, METH_NOARGS, NULL },
  FCTSPEC(set_xspath_manager, set_manager_data_path),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_NORM( agnsed, 16 ),
  XSPECMODELFCT_NORM( qsosed, 7 ),
#endif

#ifdef XSPEC_12_11_0
  XSPECMODELFCT_NORM( agnslim, 15 ),
#endif

#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_apec, 4 ),
  XSPECMODELFCT_C_NORM( C_bapec, 5 ),
#else
  XSPECMODELFCT_NORM( xsaped, 4 ),
  XSPECMODELFCT_NORM( xsbape, 5 ),
#endif
  XSPECMODELFCT_NORM( xsblbd, 2 ),
  XSPECMODELFCT_NORM( xsbbrd, 2 ),
  XSPECMODELFCT_C_NORM( C_xsbexrav, 10 ),
  XSPECMODELFCT_C_NORM( C_xsbexriv, 12 ),
  XSPECMODELFCT_C_NORM( C_brokenPowerLaw, 4 ),
  XSPECMODELFCT_C_NORM( C_broken2PowerLaw, 6 ),
  XSPECMODELFCT_C_NORM( C_sirf, 10 ),
  XSPECMODELFCT_NORM( xsbmc, 4 ),
  XSPECMODELFCT_NORM( xsbrms, 2 ),
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_brnei, 7 ),
#endif
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_bvapec, 17 ),
#else
  XSPECMODELFCT_NORM( xsbvpe, 17 ),
#endif
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_bvrnei, 19 ),
  XSPECMODELFCT_C_NORM( C_c6mekl, 11 ),
  XSPECMODELFCT_C_NORM( C_c6pmekl, 11 ),
  XSPECMODELFCT_C_NORM( C_c6pvmkl, 24 ),
  XSPECMODELFCT_C_NORM( C_c6vmekl, 24 ),
#else
  XSPECMODELFCT_NORM( c6mekl, 11 ),
  XSPECMODELFCT_NORM( c6pmekl, 11 ),
  XSPECMODELFCT_NORM( c6pvmkl, 24 ),
  XSPECMODELFCT_NORM( c6vmekl, 24 ),
#endif
  XSPECMODELFCT_NORM( cemekl, 7 ),
  XSPECMODELFCT_C_NORM( C_cemVMekal, 20 ),
  XSPECMODELFCT_C_NORM( C_xscflw, 6 ),
  XSPECMODELFCT_NORM( compbb, 4 ),
  XSPECMODELFCT_NORM( compls, 3 ),
  XSPECMODELFCT_C_NORM( C_xscompps, 20 ),
  XSPECMODELFCT_NORM( compst, 3 ),
  XSPECMODELFCT_NORM( xstitg, 6 ),
  XSPECMODELFCT_C_NORM( C_cutoffPowerLaw, 3 ),
  XSPECMODELFCT_NORM( disk, 4 ),
  XSPECMODELFCT_NORM( diskir, 9 ),
  XSPECMODELFCT_NORM( xsdskb, 2 ),
#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_diskline, 6 ),
#else
  XSPECMODELFCT_NORM( xsdili, 6 ),
#endif
  XSPECMODELFCT_NORM( diskm, 5 ),
  XSPECMODELFCT_NORM( disko, 5 ),
  XSPECMODELFCT_NORM( diskpbb, 3 ),
  XSPECMODELFCT_NORM( xsdiskpn, 3 ),
  XSPECMODELFCT_C_NORM( C_equil, 4 ),
  XSPECMODELFCT_NORM( xsxpdec, 2 ),
  XSPECMODELFCT_NORM( ezdiskbb, 2 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_gaussianLine, 3 ),
#else
  XSPECMODELFCT_NORM( xsgaul, 3 ),
#endif
  XSPECMODELFCT_C_NORM( C_gnei, 6 ),
  XSPECMODELFCT_NORM( grad, 7 ),
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( xsgrbcomp, 10 ),
#endif
  XSPECMODELFCT_NORM( xsgrbm, 4 ),
  XSPECMODELFCT_C_NORM( C_kerrbb, 10 ),

#ifdef XSPEC_12_11_0
  XSPECMODELFCT_C_NORM( C_zkerrbb, 10 ),
#endif

#ifdef XSPEC_12_10_0
  /* From an email from Craig Gordon at HEASARC:
     12.10.0 and later: for kerrd call C_kerrd. For earlier versions kerrd should call C_kerrdisk.
     (the model.dat file gives C_kerrdisk up to 12.10.1b)
  */
  XSPECMODELFCT_C_NORM( C_kerrd, 8 ),
#else
  XSPECMODELFCT_C_NORM( C_kerrdisk, 8 ),
#endif
  XSPECMODELFCT_NORM( spin, 10 ),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_NORM( kyrline, 12 ),
#endif

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_laor, 6 ),
#else
  XSPECMODELFCT_C_NORM( C_xslaor, 6 ),
#endif
  XSPECMODELFCT_C_NORM( C_laor2, 8 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_lorentzianLine, 3 ),
  XSPECMODELFCT_C_NORM( C_meka, 5 ),
  XSPECMODELFCT_C_NORM( C_mekal, 6 ),
#else
  XSPECMODELFCT_NORM( xslorz, 3 ),
  XSPECMODELFCT_NORM( xsmeka, 5 ),
  XSPECMODELFCT_NORM( xsmekl, 6 ),
#endif
  XSPECMODELFCT_C_NORM( C_xsmkcf, 6 ),
  XSPECMODELFCT_C_NORM( C_nei, 5 ),
  XSPECMODELFCT_C_NORM( C_nlapec, 4 ),
  XSPECMODELFCT_C_NORM( C_npshock, 7 ),
  XSPECMODELFCT_NORM( nsa, 5 ),
  XSPECMODELFCT_NORM( nsagrav, 4 ),
  XSPECMODELFCT_NORM( nsatmos, 5 ),
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_nsmax, 4 ),
#else
  XSPECMODELFCT_NORM( nsmax, 4 ),
#endif
  XSPECMODELFCT_C_NORM( C_xsnteea, 16 ),
  XSPECMODELFCT_C_NORM( C_nthcomp, 6 ),
  XSPECMODELFCT_NORM( xspegp, 4 ),
  XSPECMODELFCT_C_NORM( C_xspexrav, 8 ),
  XSPECMODELFCT_C_NORM( C_xspexriv, 10 ),
  XSPECMODELFCT_NORM( xsp1tr, 11 ),
  XSPECMODELFCT_C_NORM( C_powerLaw, 2 ),
  XSPECMODELFCT_NORM( xsposm, 1 ),
  XSPECMODELFCT_C_NORM( C_pshock, 6 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM(C_raysmith, 4 ),
#else
  XSPECMODELFCT_NORM( xsrays, 4 ),
#endif
  XSPECMODELFCT_NORM( xredge, 3 ),
  XSPECMODELFCT_NORM( xsrefsch, 14 ),

  XSPECMODELFCT_C_NORM( C_sedov, 6 ),
  XSPECMODELFCT_NORM( srcut, 3 ),
  XSPECMODELFCT_NORM( sresc, 3 ),
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_NORM( ssa, 3 ),
#endif
  XSPECMODELFCT_NORM( xsstep, 3 ),

#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_vapec, 16 ),
#else
  XSPECMODELFCT_NORM( xsvape, 16 ),
#endif
  XSPECMODELFCT_NORM( xsbrmv, 3 ),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_vcph, 18 ),
#endif

  XSPECMODELFCT_C_NORM( C_vequil, 15 ),
  XSPECMODELFCT_C_NORM( C_vgnei, 18 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_vmeka, 18 ),
  XSPECMODELFCT_C_NORM( C_vmekal, 19 ),
#else
  XSPECMODELFCT_NORM( xsvmek, 18 ),
  XSPECMODELFCT_NORM( xsvmkl, 19 ),
#endif
  XSPECMODELFCT_C_NORM( C_xsvmcf, 19 ),
  XSPECMODELFCT_C_NORM( C_vnei, 17 ),
  XSPECMODELFCT_C_NORM( C_vnpshock, 19 ),
  XSPECMODELFCT_C_NORM( C_vpshock, 18 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_vraysmith, 15 ),
#else
  XSPECMODELFCT_NORM( xsvrys, 15 ),
#endif
  XSPECMODELFCT_C_NORM( C_vsedov, 18 ),
  XSPECMODELFCT_NORM( xszbod, 3 ),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_zBrokenPowerLaw, 5 ),
#endif

  XSPECMODELFCT_NORM( xszbrm, 3 ),
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_zcutoffPowerLaw, 4),
#endif
  XSPECMODELFCT_C_NORM( C_xszgau, 4 ),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_zLogpar, 5 ),
#endif

  XSPECMODELFCT_C_NORM( C_zpowerLaw, 3 ),
  XSPECMODELFCT_C( C_xsabsori, 6 ),

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C( C_acisabs, 8 ),
#else
  XSPECMODELFCT( acisabs, 8 ),
#endif

  XSPECMODELFCT( xscnst, 1 ),
  XSPECMODELFCT( xscabs, 1 ),
  XSPECMODELFCT( xscycl, 5 ),
  XSPECMODELFCT( xsdust, 2 ),
  XSPECMODELFCT( xsedge, 2 ),
  XSPECMODELFCT( xsabsc, 1 ),
  XSPECMODELFCT( xsexp, 3 ),
  XSPECMODELFCT_C( C_gaussianAbsorptionLine, 3 ),
  XSPECMODELFCT( xshecu, 2 ),
  XSPECMODELFCT( xshrfl, 8 ),
  XSPECMODELFCT( xsntch, 3 ),
  XSPECMODELFCT( xsabsp, 2 ),
  XSPECMODELFCT( xsphab, 1 ),
  XSPECMODELFCT( xsplab, 2 ),
  XSPECMODELFCT_C( C_xspwab, 3 ),
  XSPECMODELFCT( xscred, 1 ),
  XSPECMODELFCT( xssmdg, 4 ),
  XSPECMODELFCT_C( C_superExpCutoff, 2 ),
  XSPECMODELFCT( xsspln, 6 ),
  XSPECMODELFCT( xssssi, 1 ),

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C( C_swind1, 4 ),
#else
  XSPECMODELFCT( swind1, 4 ),
#endif

  XSPECMODELFCT_C( C_tbabs, 1 ),
  XSPECMODELFCT_C( C_tbgrain, 6 ),
  XSPECMODELFCT_C( C_tbvabs, 42 ),
  XSPECMODELFCT( xsred, 1 ),
  XSPECMODELFCT( xsabsv, 18 ),
  XSPECMODELFCT( xsvphb, 18 ),
  XSPECMODELFCT( xsabsw, 1 ),
  XSPECMODELFCT( xswnab, 2 ),
  XSPECMODELFCT( xsxirf, 13 ),
  XSPECMODELFCT( mszdst, 4 ),
  XSPECMODELFCT( xszedg, 3 ),
  XSPECMODELFCT( xszhcu, 3 ),
  XSPECMODELFCT( xszabp, 3 ),
  XSPECMODELFCT( xszphb, 2 ),

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C( C_zxipcf, 4 ),
#else
  XSPECMODELFCT( zxipcf, 4 ),
#endif

  XSPECMODELFCT( xszcrd, 2 ),
  XSPECMODELFCT( msldst, 4 ),
  XSPECMODELFCT_C( C_ztbabs, 2 ),
  XSPECMODELFCT( xszvab, 19 ),
  XSPECMODELFCT( xszvfe, 5 ),
  XSPECMODELFCT( xszvph, 19 ),
  XSPECMODELFCT( xszabs, 2 ),
  XSPECMODELFCT( xszwnb, 3 ),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_cph, 5 ),
#endif

  XSPECMODELFCT_C_NORM( C_cplinear, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqpair, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqth, 21 ),
  XSPECMODELFCT_C_NORM( C_xscompth, 21 ),
#ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM( C_bvvapec, 34 ),
  XSPECMODELFCT_C_NORM( C_vvapec, 33 ),
#else
  XSPECMODELFCT_NORM( xsbvvp, 34 ),
  XSPECMODELFCT_NORM( xsvvap, 33 ),
#endif
#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_bvvrnei, 36 ),
#endif
  XSPECMODELFCT( zigm, 3 ),
  // New XSPEC 12.7.1 models
  XSPECMODELFCT_C_NORM( C_gaussDem, 7 ),
  XSPECMODELFCT_C_NORM( C_vgaussDem, 20 ),
  XSPECMODELFCT_NORM( eplogpar, 3 ),
#ifdef XSPEC_12_10_1
  XSPECMODELFCT_C_NORM( C_logpar, 4 ),
#else
  XSPECMODELFCT_NORM( logpar, 4 ),
#endif
  XSPECMODELFCT_NORM( optxagn, 14 ),
  XSPECMODELFCT_NORM( optxagnf, 12 ),
  XSPECMODELFCT_NORM( pexmon, 8 ),

  // Models from 12.8.0, 12.8.1, and 12.8.2
  //
  // additive
  XSPECMODELFCT_C_NORM( C_agauss, 3 ),
  XSPECMODELFCT_C_NORM( C_zagauss, 4 ),
  XSPECMODELFCT_C_NORM( xscompmag, 9 ), // DJB thinks it's okay to use the C++ wrapper for C
  XSPECMODELFCT_C_NORM( xscomptb, 7 ), // DJB thinks it's okay to use the C++ wrapper for C

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_C_NORM( C_nsmaxg, 6 ),
  XSPECMODELFCT_C_NORM( C_nsx, 6 ),
#else
  XSPECMODELFCT_NORM( nsmaxg, 6 ),
  XSPECMODELFCT_NORM( nsx, 6 ),
#endif

  XSPECMODELFCT_C_NORM( C_rnei, 6 ),
  XSPECMODELFCT_C_NORM( C_vrnei, 18 ),
  XSPECMODELFCT_C_NORM( C_vvrnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvgnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvnei, 34 ),
  XSPECMODELFCT_C_NORM( C_vvnpshock, 36 ),
  XSPECMODELFCT_C_NORM( C_vvpshock, 35 ),
  XSPECMODELFCT_C_NORM( C_vvsedov, 35 ),

#ifdef XSPEC_12_11_0
  // We do not have a direct interface to the c_func routines, so
  // take advantage of the fact XSPEC provides multiple APIs
  // and use the C one.
  // XSPECMODELFCT_C_NORM( c_beckerwolff, 13 ),
  XSPECMODELFCT_C_NORM( beckerwolff, 13 ),
#endif

  //multiplicative
  XSPECMODELFCT( xsphei, 3 ),
  XSPECMODELFCT( xslyman, 4 ),
  XSPECMODELFCT_C( xszbabs, 4 ), // DJB thinks it's okay to use the C++ wrapper for C

  // XSPEC table models
#ifdef XSPEC_12_10_1
  XSPECTABLEMODEL,
#else
  XSPECTABLEMODEL_NORM( xsatbl ),
  XSPECTABLEMODEL( xsmtbl ),
#endif

  // XSPEC convolution models
  XSPECMODELFCT_CON(C_cflux, 3),
  XSPECMODELFCT_CON(C_cpflux, 3),

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_CON(C_gsmooth, 2),
#else
  XSPECMODELFCT_CON(C_xsgsmt, 2),
#endif

  XSPECMODELFCT_CON(C_ireflct, 7),
  XSPECMODELFCT_CON(C_kdblur, 4),
  XSPECMODELFCT_CON(C_kdblur2, 6),
  XSPECMODELFCT_CON(C_spinconv, 7),

#ifdef XSPEC_12_10_1
  XSPECMODELFCT_CON_F77(kyconv, 12),
#endif

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_CON(C_lsmooth, 2),
#else
  XSPECMODELFCT_CON(C_xslsmt, 2),
#endif

  XSPECMODELFCT_CON(C_PartialCovering, 1),
  XSPECMODELFCT_CON(C_rdblur, 4),
  XSPECMODELFCT_CON(C_reflct, 5),

  XSPECMODELFCT_CON_F77(rgsxsrc, 1),
  XSPECMODELFCT_CON(C_simpl, 3),
  XSPECMODELFCT_CON(C_zashift, 1),
  XSPECMODELFCT_CON(C_zmshift, 1),

#ifdef XSPEC_12_11_0
  XSPECMODELFCT_CON_F77(thcompf, 4),
#endif

  // Models from 12.9.1
  //
  //
  #ifdef XSPEC_12_9_1
  XSPECMODELFCT_C_NORM(C_btapec, 6),
  XSPECMODELFCT_C_NORM(C_bvtapec, 18),
  XSPECMODELFCT_C_NORM(C_bvvtapec, 35),
  XSPECMODELFCT_C_NORM(C_tapec, 5),
  XSPECMODELFCT_C_NORM(C_vtapec, 17),
  XSPECMODELFCT_C_NORM(C_vvtapec, 34),

  XSPECMODELFCT_C_NORM(C_carbatm, 4),
  XSPECMODELFCT_C_NORM(C_hatm, 4),
  XSPECMODELFCT(ismabs, 31),

  XSPECMODELFCT_C_NORM(slimbbmodel, 10),
  XSPECMODELFCT_C_NORM(C_snapec, 7),
  XSPECMODELFCT_C(C_tbfeo, 4),
  XSPECMODELFCT_C(C_tbgas, 2),
  XSPECMODELFCT_C(C_tbpcf, 3),
  XSPECMODELFCT_C(C_tbrel, 42),
  XSPECMODELFCT_C_NORM(C_voigtLine, 4),
  XSPECMODELFCT_C(C_xscatmodel, 4),

  XSPECMODELFCT_CON(C_clumin, 4),
  XSPECMODELFCT_CON(C_rfxconv, 5),
  XSPECMODELFCT_CON(C_vashift, 1),
  XSPECMODELFCT_CON(C_vmshift, 1),
  XSPECMODELFCT_CON(C_xilconv, 6),
  #endif

#ifdef XSPEC_12_10_0
  XSPECMODELFCT_NORM(jet, 16),
#endif

#ifdef XSPEC_12_11_0
  XSPECMODELFCT(ismdust, 3),
  XSPECMODELFCT(olivineabs, 2),
  XSPECMODELFCT_C(C_logconst, 1),
  XSPECMODELFCT_C(C_log10con, 1),
#endif

  { NULL, NULL, 0, NULL }

};

static struct PyModuleDef xspec_module = {
        PyModuleDef_HEAD_INIT,
        "_xspec",
        NULL,
        -1,
        XSpecMethods,
};

PyMODINIT_FUNC PyInit__xspec(void) {
  import_array();
  return PyModule_Create(&xspec_module);
}
