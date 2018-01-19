//  Copyright (C) 2007, 2015-2018  Smithsonian Astrophysical Observatory
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

// C_<model> are declared here
#include "funcWrappers.h"

extern "C" {

void xsaped_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsblbd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbbrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbmc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbrms_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbvpe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6mekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6pvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void c6vmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void cemekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compls_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void compst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xstitg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void disk_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskir_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdskb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdili_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void disko_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void diskpbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsdiskpn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsxpdec_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void ezdiskbb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsgaul_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void grad_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsgrbm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void spin_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xslorz_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmeka_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsmekl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsa_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsagrav_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsatmos_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsmax_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xspegp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsp1tr_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsposm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsrays_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xredge_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsrefsch_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void srcut_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void sresc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsstep_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvape_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsbrmv_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvmek_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvmkl_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvrys_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbod_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbrm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void acisabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
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
void swind1_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
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
void zxipcf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszcrd_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void msldst_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvab_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvfe_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszvph_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszabs_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszwnb_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// New XSPEC 12.7 models

void xsbvvp_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xsvvap_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void zigm_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// New XSPEC 12.7.1 models

void logpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void eplogpar_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagn_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void optxagnf_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void pexmon_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

// Models from 12.8.0, 12.8.1, and 12.8.2

// additive
void xscompmag(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void xscomptb(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void nsmaxg_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void nsx_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);

//multiplicative
void xsphei_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xslyman_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void xszbabs(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

// XSPEC table models
void xsatbl(float* ear, int ne, float* param, const char* filenm, int ifl, 
	    float* photar, float* photer);
void xsmtbl(float* ear, int ne, float* param, const char* filenm, int ifl, 
	    float* photar, float* photer);

// XSPEC convolution models: ordering below matches that of model.dat,
// apart from those added in 12.9.1 have been moved into the
// ifdef section below.
//
void C_cflux(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_cpflux(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xsgsmt(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_ireflct(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_kdblur(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_kdblur2(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_spinconv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xslsmt(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_PartialCovering(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_rdblur(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_reflct(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
// rgsxsrc is the only convolution-style model that uses the Fortran interface
void rgsxsrc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);
void C_simpl(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zashift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_zmshift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

// Models from 12.9.1
//
//
#ifdef XSPEC_12_9_1
void C_btapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_bvtapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_bvvtapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_tapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vtapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vvtapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

void C_carbatm(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_hatm(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
  // void ismabs(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);  // is this the Fortran interface?
  // void c_slimbh(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr); // note the lower-case c in the prefix
void C_snapec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_tbfeo(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_tbgas(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_tbpcf(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_tbrel(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

void C_clumin(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_rfxconv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vashift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_vmshift(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);
void C_xilconv(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr);

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

// for documentation
#define SEEALSODOC "\nSee also\n--------\n"
#define NOTESDOC "\nNotes\n-----\n"
#define REFERENCESDOC "\nReferences\n----------\n"
#define EXAMPLESDOC "\nExamples\n--------\n"
#define PARAMETERSDOC "\nParameters\n----------\n"
#define RETURNSDOC "\nReturns\n-------\n"

static PyMethodDef XSpecMethods[] = {
  { (char*)"get_xsversion", (PyCFunction)get_version, METH_NOARGS,
    (char*) "get_xsversion()\n\n"
            "Return the version of the X-Spec model library in use.\n"
            RETURNSDOC
            "version : str\n"
            "   The version of the X-Spec model library used\n"
            "   by Sherpa [1]_.\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/docs/xanadu/xspec/\n"
            EXAMPLESDOC "\n"
            ">>> get_xsversion()\n'12.9.1p'\n\n"},

  { (char*)"get_xschatter", (PyCFunction)get_chatter, METH_NOARGS,
    (char*) "get_xschatter()\n\n"
            "Return the chatter level used by X-Spec.\n"
            RETURNSDOC
            "chatter : int\n"
            "   The chatter setting used by the X-Spec routines.\n"
            SEEALSODOC
            "set_xschatter : Set the X-Spec chatter level.\n"
            EXAMPLESDOC "\n"
            ">>> get_xschatter()\n0\n\n"},
  { (char*)"set_xschatter", (PyCFunction)set_chatter, METH_VARARGS,
    (char*) "set_xschatter(level)\n\n"
            "Set the chatter level used by X-Spec.\n\n"
            "Set the chatter setting used by the X-Spec routines\n"
            "for determining what information gets printed to the\n"
            "screen. It is equivalent to the X-Spec ``chatter``\n"
            "command [1]_.\n"
            PARAMETERSDOC
            "level : int\n"
            "   The higher the value of ``level``, the more screen output will\n"
            "   be created by X-Spec routines. A value of ``0`` hides\n"
            "   most information while ``25`` will generate a lot of\n"
            "   debug output.\n"
            SEEALSODOC
            "get_xschatter : Return the X-Spec chatter setting.\n"
            NOTESDOC
            "The default chatter setting used by Sherpa is ``0``, which\n"
            "is lower than - so, creates less screen output - the\n"
            "default value used by X-Spec (``10``).\n\n"
            "There is no way to change the X-Spec \"log chatter\"\n"
            "setting.\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSchatter.html\n"
            "       Note that this may refer to a newer version than the\n"
            "       compiled version used by Sherpa; use `get_xsversion` to\n"
            "       check.\n\n"
            EXAMPLESDOC "\n"
            "Set the chatter level to the default used by X-Spec:\n\n"
            ">>> set_xschatter(10)\n\n"},

  { (char*)"get_xsabund", (PyCFunction)get_abund, METH_VARARGS,
    (char*) "get_xsabund(element=None)\n\n"
            "Return the X-Spec abundance setting or elemental abundance.\n"
            PARAMETERSDOC
            "element : str, optional\n"
            "   When not given, the abundance table name is returned.\n"
            "   If a string, then it must be an element name from:\n"
            "   'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',\n"
            "   'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',\n"
            "   'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',\n"
            "   'Cu', 'Zn'. Case is important.\n"
            RETURNSDOC
            "val : str or float\n"
            "   When ``element`` is ``None``, the abundance table name is\n"
            "   returned (see `set_xsabund`); the string 'file' is\n"
            "   used when the abundances were read from a file. A\n"
            "   numeric value is returned when an element name is\n"
            "   given. This value is the elemental abundance relative\n"
            "   to H.\n"
            SEEALSODOC
            "set_xsabund : Set the X-Spec abundance table.\n"
            EXAMPLESDOC "\n"
            "Return the current abundance setting, which in this case\n"
            "is 'angr', the default value for X-Spec:\n\n"
            ">>> get_xsabund()\n'angr'\n\n"
            "The `set_xsabund` function has been used to read in the\n"
            "abundances from a file, so the routine now returns the\n"
            "string 'file':\n\n"
            ">>> set_xsabund('abund.dat')\n"
            ">>> get_xsabund()\n'file'\n\n"
            ">>> get_xsabund('He')\n0.09769999980926514\n\n"},
  { (char*)"set_xsabund", (PyCFunction)set_abund, METH_VARARGS,
    (char*) "set_xsabund(abundance)\n\n"
            "Set the elemental abundances used by X-Spec models.\n\n"
            "Set the abundance table used in the X-Spec plasma emission and\n"
            "photoelectric absorption models. It is equivalent to the X-Spec\n"
            "``abund`` command [1]_.\n"
            PARAMETERSDOC
            "abundance : str\n"
            "   A file name, format described below, or one of the\n"
            "   pre-defined names listed below.\n"
            SEEALSODOC
            "get_xsabund : Return the X-Spec abundance setting or elemental abundance.\n"
            "get_xsversion : Return the version of X-Spec used by this module.\n"
            "set_xschatter : Control the screen output of X-Spec functions and models.\n"
            NOTESDOC
            "The pre-defined abundance tables are:\n\n"
            " - 'angr', from [2]_\n"
            " - 'aspl', from [3]_\n"
            " - 'feld', from [4]_, except for elements not listed which\n"
            "   are given 'grsa' abundances\n"
            " - 'aneb', from [5]_\n"
            " - 'grsa', from [6]_\n"
            " - 'wilm', from [7]_, except for elements not listed which\n"
            "   are given zero abundance\n"
            " - 'lodd', from [8]_\n\n"
            "The values for these tables are given at [1]_.\n\n"
            "Data files should be in ASCII format, containing a single\n"
            "numeric (floating-point) column of the abundance values,\n"
            "relative to Hydrogen.\n\n"
            "The screen output of this function is controlled by the\n"
            "X-Spec chatter setting (`set_xschatter`).\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html\n"
            "       Note that this may refer to a newer version than the\n"
            "       compiled version used by Sherpa; use `get_xsversion` to\n"
            "       check.\n\n"
            ".. [2] Anders E. & Grevesse N. (1989, Geochimica et\n"
            "       Cosmochimica Acta 53, 197)\n"
            "       http://adsabs.harvard.edu/abs/1989GeCoA..53..197A\n\n"
            ".. [3] Asplund M., Grevesse N., Sauval A.J. & Scott P.\n"
            "       (2009, ARAA, 47, 481)\n"
            "       http://adsabs.harvard.edu/abs/2009ARA%26A..47..481A\n\n"
            ".. [4] Feldman U.(1992, Physica Scripta 46, 202)\n"
            "       http://adsabs.harvard.edu/abs/1992PhyS...46..202F\n\n"
            ".. [5] Anders E. & Ebihara (1982, Geochimica et Cosmochimica\n"
            "       Acta 46, 2363)\n"
            "       http://adsabs.harvard.edu/abs/1982GeCoA..46.2363A\n\n"
            ".. [6] Grevesse, N. & Sauval, A.J. (1998, Space Science\n"
            "       Reviews 85, 161)\n"
            "       http://adsabs.harvard.edu/abs/1998SSRv...85..161G\n\n"
            ".. [7] Wilms, Allen & McCray (2000, ApJ 542, 914)\n"
            "       http://adsabs.harvard.edu/abs/2000ApJ...542..914W\n\n"
            ".. [8] Lodders, K (2003, ApJ 591, 1220)\n"
            "       http://adsabs.harvard.edu/abs/2003ApJ...591.1220L\n\n"
            EXAMPLESDOC "\n"
            ">>> set_xsabund('lodd')\n"
            " Solar Abundance Vector set to lodd:  Lodders, K. ApJ 591, 1220 (2003)\n\n"
            ">>> set_xsabund('abund.dat')\n"
            " Solar Abundance Vector set to file:  User defined abundance vector / no description specified\n\n"},

  { (char*)"set_xscosmo", (PyCFunction)set_cosmo, METH_VARARGS,
    (char*) "set_xscosmo(h0, q0, l0)\n\n"
            "Set the cosmological parameters used by X-Spec models.\n\n"
            "Set the cosmological parameters (H_0, q_0, lambda_0) used\n"
            "by X-Spec. It is equivalent to the X-Spec ``cosmo``\n"
            "command [1]_. The default values are h0=70, q0=0, and l0=0.73\n"
            PARAMETERSDOC
            "h0 : number\n"
            "   The Hubble constant in km/s/Mpc.\n"
            "q0 : number\n"
            "   The deceleration parameter.\n"
            "l0 : number\n"
            "   The cosmological constant. If this is non-zero\n"
            "   then the q0 parameter is ignored and the Universe\n"
            "   is assumed to be flat.\n"
            SEEALSODOC
            "get_xscosmo : Return the X-Spec cosmology settings.\n"
            "get_xsversion : Return the version of X-Spec used by this module.\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XScosmo.html\n"
            "       Note that this may refer to a newer version than the\n"
            "       compiled version used by Sherpa; use `get_xsversion` to\n"
            "       check.\n\n"
            EXAMPLESDOC "\n"
            ">>> set_xscosmo(73, 0, 0.73)\n"},
  { (char*)"get_xscosmo", (PyCFunction)get_cosmo, METH_NOARGS,
    (char*) "get_xscosmo()\n\n"
            "Return the X-Spec cosmology settings.\n"
            RETURNSDOC
            "(h0, q0, l0)\n"
            "   The Hubble constant, in km/s/Mpc, the deceleration\n"
            "   parameter, and the cosmological constant.\n"
            SEEALSODOC
            "set_xscosmo : Set the X-Spec cosmology settings.\n"
            EXAMPLESDOC "\n"
            ">>> get_xscosmo()\n(70.0, 0.0, 0.7300000190734863)\n\n"},

  { (char*)"get_xsxsect", (PyCFunction)get_cross, METH_NOARGS,
    (char*) "get_xsxsect()\n\n"
            "Return the cross sections used by X-Spec models.\n\n"
            "Return the name of the X-Spec photoelectric absorption "
            "cross-sections setting.\n"
            RETURNSDOC
            "val : str\n"
            "   The value of the photoelectric absorption setting:\n"
            "   one of 'bcmc', 'obcm', and 'vern'.\n"
            SEEALSODOC
            "set_xsxsect : Set the X-Spec photoelectric absorption cross-sections.\n"
            EXAMPLESDOC "\n"
            ">>> get_xsxsect()\n'bcmc'\n\n"},
  { (char*)"set_xsxsect", (PyCFunction)set_cross, METH_VARARGS,
    (char*) "set_xsxsect(name)\n\n"
            "Set the cross sections used by X-Spec models.\n\n"
            "Set the X-Spec photoelectric absorption cross-sections\n"
            "setting, which changes the cross-sections used by all\n"
            "X-Spec absorption models *except* for `XSwabs`. It is\n"
            "equivalent to the X-Spec ``xsect`` command [1]_.\n"
            PARAMETERSDOC
            "name : { 'bcmc', 'obcm', 'vern' }\n"
            "   The options are: 'bcmc' from [2]_ with a new\n"
            "   He cross-section based on [3]_; 'obcm' which is,\n"
            "   the same as 'bcmc', but with the He cross-section\n"
            "   from [2]_, or 'vern' [4]_.\n"
            SEEALSODOC
            "get_xsxsect : Return the name of the X-Spec photoelectric absorption cross-sections.\n"
            "get_xsversion : Return the version of X-Spec used by this module.\n"
            "set_xschatter : Control the screen output of X-Spec functions and models.\n"
            NOTESDOC
            "The screen output of this function is controlled by the\n"
            "X-Spec chatter setting (`set_xschatter`).\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSxsect.html\n"
            "       Note that this may refer to a newer version than the\n"
            "       compiled version used by Sherpa; use `get_xsversion` to\n"
            "       check.\n\n"
            ".. [2] Balucinska-Church & McCammon (1992; Ap.J.400, 699).\n"
            "       http://adsabs.harvard.edu/abs/1992ApJ...400..699B\n\n"
            ".. [3] Yan, M., Sadeghpour, H. R., Dalgarno, A. 1998,\n"
            "       Ap.J. 496, 1044\n"
            "       http://adsabs.harvard.edu/abs/1998ApJ...496.1044Y\n\n"
            ".. [4] Verner et. al., 1996, Ap.J., 465, 487.\n"
            "       http://adsabs.harvard.edu/abs/1996ApJ...465..487V\n\n"
            EXAMPLESDOC "\n"
            ">>> set_xsxsect('vern')\n"
            " Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996\n\n"},

  { (char*)"set_xsxset", (PyCFunction)set_xset, METH_VARARGS,
    (char*) "set_xsxset(name, value)\n\n"
            "Set an X-Spec XSET variable to a value.\n\n"
            "Set variables used by X-Spec models. It is equivalent to the\n"
            "X-Spec ``xset`` command [1]_, but only for setting the model\n"
            "database settings. See `set_xsabund`, `set_xscosmo`, and\n"
            "`set_xsxsect` for the other settings.\n"
            PARAMETERSDOC
            "name : str\n"
            "   The name of the setting. It is case insensitive, and\n"
            "   is converted to upper case. There is no check that the\n"
            "   name is valid.\n"
            "value : str\n"
            "   The new value of the setting. It must be given as a\n"
            "   string.\n"
            SEEALSODOC
            "get_xsxset : Return the value of X-Spec model settings.\n"
            "get_xsversion : Return the version of X-Spec used by this module.\n"
            "set_xsabund : Set the X-Spec abundance table.\n"
            "set_xschatter : Set the X-Spec chatter level.\n"
            "set_xscosmo : Set the X-Spec cosmology settings.\n"
            "set_xsxsect : Set the X-Spec photoelectric absorption cross-sections.\n"
            NOTESDOC
            "The available settings are listed at [1]_. Not all the\n"
            "X-Spec model types are supported by Sherpa - for instance\n"
            "X-Spec \"mixing models\" - so changing some of these settings\n"
            "will make no difference. The X-Spec chatter setting can be\n"
            "increased with `set_xschatter` if it is not clear if\n"
            "a setting is being used.\n"
            REFERENCESDOC "\n"
            ".. [1] http://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html\n"
            "       Note that this may refer to a newer version than the\n"
            "       compiled version used by Sherpa; use `get_xsversion` to\n"
            "       check.\n\n"
            EXAMPLESDOC "\n"
	    ">>> set_xsxset('NEIVERS', '2.0')\n\n"
	    ">>> set_xsxset('NEIAPECROOT', '/data/spectral/modelData/APEC_nei_v11')\n\n"
            ">>> set_xsxset('POW_EMIN', '0.5')\n"
            ">>> set_xsxset('POW_EMAX', '2.0')\n\n"},
  { (char*)"get_xsxset", (PyCFunction)get_xset, METH_VARARGS,
    (char*) "get_xsxset(name)\n\n"
            "Return the value of an X-Spec XSET variable.\n\n"
            PARAMETERSDOC
            "name : str\n"
            "   The name of the setting. There is no check that\n"
            "   the name is valid.\n"
            RETURNSDOC
            "val : str\n"
            "   Returns the value set by a previous call to `set_xsxset`\n"
            "   or the empty string, if not set.\n"
            SEEALSODOC
            "set_xsxset : Set a X-Spec model setting.\n"
            NOTESDOC
            "Due to the way X-Spec model settings work, `get_xsxset`\n"
            "will only return a value if it has previously been set\n"
            "with a call to `set_xsxset`. There is no way to retrive\n"
            "the default value of a setting.\n\n"},

  { (char*)"get_xspath_manager",
    (PyCFunction)get_manager_data_path, METH_NOARGS,
    (char*) "get_xspath_manager()\n\n"
            "Return the path to the files describing the XSPEC models.\n"
            RETURNSDOC
            "path : str\n"
            "   The path to the manager directory containing the various\n"
            "   *.dat files used by XSPEC.\n"
            SEEALSODOC
            "get_xspath_model : Return the path to the model data files.\n"
            "set_xspath_manager : Set the path to the files describing the XSPEC models.\n"
            EXAMPLESDOC "\n"
            ">>> get_xspath_manager()\n"
            "'/usr/local/heasoft-6.22/x86_64-unknown-linux-gnu-libc2.24/../spectral/manager'\n\n"},

  { (char*)"get_xspath_model",
    (PyCFunction)get_model_data_path, METH_NOARGS,
    (char*) "get_xspath_model()\n\n"
            "Return the path to the model data files.\n"
            RETURNSDOC
            "path : str\n"
            "   The path to the directory containing the files used by\n"
            "   the XSPEC models.\n"
            SEEALSODOC
            "get_xspath_manager : Return the path to the files describing the XSPEC models.\n"
            EXAMPLESDOC "\n"
            ">>> get_xspath_model()\n"
            "'/usr/local/heasoft-6.22/x86_64-unknown-linux-gnu-libc2.24/../spectral/modelData'\n\n"},

  { (char*)"set_xspath_manager",
    (PyCFunction)set_manager_data_path, METH_VARARGS,
    (char*) "set_xspath_manager(path)\n\n"
            "Set the path to the files describing the XSPEC models.\n"
            PARAMETERSDOC
            "path : str\n"
            "   The new path.\n"
            SEEALSODOC
            "get_xspath_manager : Return the path to the files describing the XSPEC models.\n"
            EXAMPLESDOC "\n"
            ">>> set_xspath_manager('/data/xspec/spectral/manager')\n\n"},

  XSPECMODELFCT_NORM( xsaped, 4 ),
  XSPECMODELFCT_NORM( xsbape, 5 ),
  XSPECMODELFCT_NORM( xsblbd, 2 ),
  XSPECMODELFCT_NORM( xsbbrd, 2 ),
  XSPECMODELFCT_C_NORM( C_xsbexrav, 10 ),
  XSPECMODELFCT_C_NORM( C_xsbexriv, 12 ),
  XSPECMODELFCT_C_NORM( C_brokenPowerLaw, 4 ),
  XSPECMODELFCT_C_NORM( C_broken2PowerLaw, 6 ),
  XSPECMODELFCT_C_NORM( C_sirf, 10 ),
  XSPECMODELFCT_NORM( xsbmc, 4 ),
  XSPECMODELFCT_NORM( xsbrms, 2 ),
  XSPECMODELFCT_NORM( xsbvpe, 17 ),
  XSPECMODELFCT_NORM( c6mekl, 11 ),
  XSPECMODELFCT_NORM( c6pmekl, 11 ),
  XSPECMODELFCT_NORM( c6pvmkl, 24 ),
  XSPECMODELFCT_NORM( c6vmekl, 24 ),
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
  XSPECMODELFCT_NORM( xsdili, 6 ),
  XSPECMODELFCT_NORM( diskm, 5 ),
  XSPECMODELFCT_NORM( disko, 5 ),
  XSPECMODELFCT_NORM( diskpbb, 3 ),
  XSPECMODELFCT_NORM( xsdiskpn, 3 ),
  XSPECMODELFCT_C_NORM( C_equil, 4 ),
  XSPECMODELFCT_NORM( xsxpdec, 2 ),
  XSPECMODELFCT_NORM( ezdiskbb, 2 ),
  XSPECMODELFCT_NORM( xsgaul, 3 ),
  XSPECMODELFCT_C_NORM( C_gnei, 6 ),
  XSPECMODELFCT_NORM( grad, 7 ),
  XSPECMODELFCT_NORM( xsgrbm, 4 ),
  XSPECMODELFCT_C_NORM( C_kerrbb, 10 ),
  XSPECMODELFCT_C_NORM( C_kerrdisk, 8 ),
  XSPECMODELFCT_NORM( spin, 10 ),
  XSPECMODELFCT_C_NORM( C_xslaor, 6 ),
  XSPECMODELFCT_C_NORM( C_laor2, 8 ),
  XSPECMODELFCT_NORM( xslorz, 3 ),
  XSPECMODELFCT_NORM( xsmeka, 5 ),
  XSPECMODELFCT_NORM( xsmekl, 6 ),
  XSPECMODELFCT_C_NORM( C_xsmkcf, 6 ),
  XSPECMODELFCT_C_NORM( C_nei, 5 ),
  XSPECMODELFCT_C_NORM( C_nlapec, 4 ),
  XSPECMODELFCT_C_NORM( C_npshock, 7 ),
  XSPECMODELFCT_NORM( nsa, 5 ),
  XSPECMODELFCT_NORM( nsagrav, 4 ),
  XSPECMODELFCT_NORM( nsatmos, 5 ),
  XSPECMODELFCT_NORM( nsmax, 4 ),
  XSPECMODELFCT_C_NORM( C_xsnteea, 16 ),
  XSPECMODELFCT_C_NORM( C_nthcomp, 6 ),
  XSPECMODELFCT_NORM( xspegp, 4 ),
  XSPECMODELFCT_C_NORM( C_xspexrav, 8 ),
  XSPECMODELFCT_C_NORM( C_xspexriv, 10 ),
  XSPECMODELFCT_NORM( xsp1tr, 11 ),
  XSPECMODELFCT_C_NORM( C_powerLaw, 2 ),
  XSPECMODELFCT_NORM( xsposm, 1 ),
  XSPECMODELFCT_C_NORM( C_pshock, 6 ),
  XSPECMODELFCT_NORM( xsrays, 4 ),
  XSPECMODELFCT_NORM( xredge, 3 ),
  XSPECMODELFCT_NORM( xsrefsch, 14 ),
  XSPECMODELFCT_C_NORM( C_sedov, 6 ),
  XSPECMODELFCT_NORM( srcut, 3 ),
  XSPECMODELFCT_NORM( sresc, 3 ),
  XSPECMODELFCT_NORM( xsstep, 3 ),
  XSPECMODELFCT_NORM( xsvape, 16 ),
  XSPECMODELFCT_NORM( xsbrmv, 3 ),
  XSPECMODELFCT_C_NORM( C_vequil, 15 ),
  XSPECMODELFCT_C_NORM( C_vgnei, 18 ),
  XSPECMODELFCT_NORM( xsvmek, 18 ),
  XSPECMODELFCT_NORM( xsvmkl, 19 ),
  XSPECMODELFCT_C_NORM( C_xsvmcf, 19 ),
  XSPECMODELFCT_C_NORM( C_vnei, 17 ),
  XSPECMODELFCT_C_NORM( C_vnpshock, 19 ),
  XSPECMODELFCT_C_NORM( C_vpshock, 18 ),
  XSPECMODELFCT_NORM( xsvrys, 15 ),
  XSPECMODELFCT_C_NORM( C_vsedov, 18 ),
  XSPECMODELFCT_NORM( xszbod, 3 ),
  XSPECMODELFCT_NORM( xszbrm, 3 ),
  XSPECMODELFCT_C_NORM( C_xszgau, 4 ),
  XSPECMODELFCT_C_NORM( C_zpowerLaw, 3 ),
  XSPECMODELFCT_C( C_xsabsori, 6 ),
  XSPECMODELFCT( acisabs, 8 ),
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
  XSPECMODELFCT( swind1, 4 ),
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
  XSPECMODELFCT( zxipcf, 4 ),
  XSPECMODELFCT( xszcrd, 2 ),
  XSPECMODELFCT( msldst, 4 ),
  XSPECMODELFCT_C( C_ztbabs, 2 ),
  XSPECMODELFCT( xszvab, 19 ),
  XSPECMODELFCT( xszvfe, 5 ),
  XSPECMODELFCT( xszvph, 19 ),
  XSPECMODELFCT( xszabs, 2 ),
  XSPECMODELFCT( xszwnb, 3 ),
  // New XSPEC 12.7 models
  XSPECMODELFCT_C_NORM( C_cplinear, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqpair, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqth, 21 ),
  XSPECMODELFCT_C_NORM( C_xscompth, 21 ),
  XSPECMODELFCT_NORM( xsbvvp, 34 ),
  XSPECMODELFCT_NORM( xsvvap, 33 ),
  XSPECMODELFCT( zigm, 3 ),
  // New XSPEC 12.7.1 models
  XSPECMODELFCT_C_NORM( C_gaussDem, 7 ),
  XSPECMODELFCT_C_NORM( C_vgaussDem, 20 ),
  XSPECMODELFCT_NORM( eplogpar, 3 ),
  XSPECMODELFCT_NORM( logpar, 4 ),
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
  XSPECMODELFCT_NORM( nsmaxg, 6 ),
  XSPECMODELFCT_NORM( nsx, 6 ),
  XSPECMODELFCT_C_NORM( C_rnei, 6 ),
  XSPECMODELFCT_C_NORM( C_vrnei, 18 ),
  XSPECMODELFCT_C_NORM( C_vvrnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvgnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvnei, 34 ),
  XSPECMODELFCT_C_NORM( C_vvnpshock, 36 ),
  XSPECMODELFCT_C_NORM( C_vvpshock, 35 ),
  XSPECMODELFCT_C_NORM( C_vvsedov, 35 ),

  //multiplicative
  XSPECMODELFCT( xsphei, 3 ),
  XSPECMODELFCT( xslyman, 4 ),
  XSPECMODELFCT_C( xszbabs, 4 ), // DJB thinks it's okay to use the C++ wrapper for C

  // XSPEC table models
  XSPECTABLEMODEL_NORM( xsatbl ),
  XSPECTABLEMODEL_NORM( xsmtbl ),

  // XSPEC convolution models
  XSPECMODELFCT_CON(C_cflux, 3),
  XSPECMODELFCT_CON(C_cpflux, 3),
  XSPECMODELFCT_CON(C_xsgsmt, 2),
  XSPECMODELFCT_CON(C_ireflct, 7),
  XSPECMODELFCT_CON(C_kdblur, 4),
  XSPECMODELFCT_CON(C_kdblur2, 6),
  XSPECMODELFCT_CON(C_spinconv, 7),
  XSPECMODELFCT_CON(C_xslsmt, 2),
  XSPECMODELFCT_CON(C_PartialCovering, 1),
  XSPECMODELFCT_CON(C_rdblur, 4),
  XSPECMODELFCT_CON(C_reflct, 5),
  XSPECMODELFCT_CON_F77(rgsxsrc, 1),
  XSPECMODELFCT_CON(C_simpl, 3),
  XSPECMODELFCT_CON(C_zashift, 1),
  XSPECMODELFCT_CON(C_zmshift, 1),

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

  XSPECMODELFCT_C_NORM(C_snapec, 7),
  XSPECMODELFCT_C_NORM(C_tbfeo, 4),
  XSPECMODELFCT_C_NORM(C_tbgas, 2),
  XSPECMODELFCT_C_NORM(C_tbpcf, 3),
  XSPECMODELFCT_C_NORM(C_tbrel, 42),

  XSPECMODELFCT_CON(C_clumin, 4),
  XSPECMODELFCT_CON(C_rfxconv, 5),
  XSPECMODELFCT_CON(C_vashift, 1),
  XSPECMODELFCT_CON(C_vmshift, 1),
  XSPECMODELFCT_CON(C_xilconv, 6),
  #endif

  { NULL, NULL, 0, NULL }

};


#ifdef PY3

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

#else

PyMODINIT_FUNC
init_xspec(void) {
  import_array();
  Py_InitModule( (char*)"_xspec", XSpecMethods );
}
#endif
