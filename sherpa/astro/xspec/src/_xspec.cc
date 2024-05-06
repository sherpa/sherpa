//  Copyright (C) 2007, 2015 - 2024
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
// TODO:: switch to C++ FuntionUtility interface rather than use xsFortran.h
//
// funcWrappers: C_<model> are declared here; the other models are defined in
// functionMap.h
//
#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/funcWrappers.h>
#include <XSFunctions/functionMap.h>

// TODO: is this defined in an XSPEC header file?
#define ABUND_SIZE (30) // number of elements in Solar Abundance table


// The XSPEC initialization used to be done lazily - that is, only
// when the first routine from XSPEC was about to be called - but
// the module is now set up so that we need to know the version
// of XSPEC being used when the sherpa.astro.xspec module is
// being created. As the version requires a run-time check (that
// is, the get_version call is made) then we know that we need to
// initialize the XSPEC code when the Python module is installed.
// So we no-longer need to support the lazy loading.
//

static int _sherpa_init_xspec_library()
{

  // This routine is only meant to be called once
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

    // We used to set the chatter to 0 but this misses useful info
    // (like XSPEC can't find the data files) that would be reported
    // by XSPEC, so use the default XSPEC setting. It does not appear
    // that this value is read in from ~/.xspec/Xspec.init so set
    // it here so we have repeatable behavior.
    //
    FPCHAT( 10 );

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

  // Start model definitions

  XSPECMODELFCT_C_NORM(C_agauss, 3),               // XSagauss
  XSPECMODELFCT_NORM(agnsed, 16),                  // XSagnsed
  XSPECMODELFCT_NORM(agnslim, 15),                 // XSagnslim
  XSPECMODELFCT_C_NORM(C_apec, 4),                 // XSapec
  XSPECMODELFCT_C_NORM(C_bapec, 5),                // XSbapec
  XSPECMODELFCT_C_NORM(C_btapec, 6),               // XSbtapec
  XSPECMODELFCT_NORM(xsblbd, 2),                   // XSbbody
  XSPECMODELFCT_NORM(xsbbrd, 2),                   // XSbbodyrad
  XSPECMODELFCT_C_NORM(C_xsbexrav, 10),            // XSbexrav
  XSPECMODELFCT_C_NORM(C_xsbexriv, 12),            // XSbexriv
  XSPECMODELFCT_C_NORM(C_brokenPowerLaw, 4),       // XSbknpower
  XSPECMODELFCT_C_NORM(C_broken2PowerLaw, 6),      // XSbkn2pow
  XSPECMODELFCT_NORM(xsbmc, 4),                    // XSbmc
  XSPECMODELFCT_NORM(xsbrms, 2),                   // XSbremss
  XSPECMODELFCT_C_NORM(C_brnei, 7),                // XSbrnei
  XSPECMODELFCT_C_NORM(C_bvapec, 17),              // XSbvapec
  XSPECMODELFCT_C_NORM(C_bvrnei, 19),              // XSbvrnei
  XSPECMODELFCT_C_NORM(C_bvtapec, 18),             // XSbvtapec
  XSPECMODELFCT_C_NORM(C_bvvapec, 34),             // XSbvvapec
  XSPECMODELFCT_C_NORM(C_bvvrnei, 36),             // XSbvvrnei
  XSPECMODELFCT_C_NORM(C_bvvtapec, 35),            // XSbvvtapec
  XSPECMODELFCT_C_NORM(C_c6mekl, 11),              // XSc6mekl
  XSPECMODELFCT_C_NORM(C_c6pmekl, 11),             // XSc6pmekl
  XSPECMODELFCT_C_NORM(C_c6pvmkl, 24),             // XSc6pvmkl
  XSPECMODELFCT_C_NORM(C_c6vmekl, 24),             // XSc6vmekl
  XSPECMODELFCT_C_NORM(C_carbatm, 4),              // XScarbatm
  XSPECMODELFCT_NORM(cemekl, 7),                   // XScemekl
  XSPECMODELFCT_C_NORM(C_cemVMekal, 20),           // XScevmkl
  XSPECMODELFCT_C_NORM(C_xscflw, 6),               // XScflow
  XSPECMODELFCT_NORM(compbb, 4),                   // XScompbb
  XSPECMODELFCT_C_NORM(xscompmag, 9),              // XScompmag
  XSPECMODELFCT_NORM(compls, 3),                   // XScompLS
  XSPECMODELFCT_C_NORM(C_xscompps, 20),            // XScompPS
  XSPECMODELFCT_NORM(compst, 3),                   // XScompST
  XSPECMODELFCT_C_NORM(xscomptb, 7),               // XScomptb
  XSPECMODELFCT_C_NORM(C_xscompth, 21),            // XScompth
  XSPECMODELFCT_NORM(xstitg, 6),                   // XScompTT
  XSPECMODELFCT_C_NORM(C_cph, 5),                  // XScph
  XSPECMODELFCT_C_NORM(C_cplinear, 21),            // XScplinear
  XSPECMODELFCT_C_NORM(C_cutoffPowerLaw, 3),       // XScutoffpl
  XSPECMODELFCT_NORM(disk, 4),                     // XSdisk
  XSPECMODELFCT_NORM(diskir, 9),                   // XSdiskir
  XSPECMODELFCT_NORM(xsdskb, 2),                   // XSdiskbb
  XSPECMODELFCT_C_NORM(C_diskline, 6),             // XSdiskline
  XSPECMODELFCT_NORM(diskm, 5),                    // XSdiskm
  XSPECMODELFCT_NORM(disko, 5),                    // XSdisko
  XSPECMODELFCT_NORM(diskpbb, 3),                  // XSdiskpbb
  XSPECMODELFCT_NORM(xsdiskpn, 3),                 // XSdiskpn
  XSPECMODELFCT_NORM(eplogpar, 3),                 // XSeplogpar
  XSPECMODELFCT_C_NORM(C_xseqpair, 21),            // XSeqpair
  XSPECMODELFCT_C_NORM(C_xseqth, 21),              // XSeqtherm
  XSPECMODELFCT_C_NORM(C_equil, 4),                // XSequil
  XSPECMODELFCT_NORM(xsxpdec, 2),                  // XSexpdec
  XSPECMODELFCT_NORM(ezdiskbb, 2),                 // XSezdiskbb
  XSPECMODELFCT_C_NORM(C_gaussianLine, 3),         // XSgaussian
  XSPECMODELFCT_C_NORM(C_gaussDem, 7),             // XSgadem
  XSPECMODELFCT_C_NORM(C_gnei, 6),                 // XSgnei
  XSPECMODELFCT_NORM(grad, 7),                     // XSgrad
  XSPECMODELFCT_C_NORM(xsgrbcomp, 10),             // XSgrbcomp
  XSPECMODELFCT_C_NORM(xsgrbjet, 14),              // XSgrbjet
  XSPECMODELFCT_NORM(xsgrbm, 4),                   // XSgrbm
  XSPECMODELFCT_C_NORM(C_hatm, 4),                 // XShatm
  XSPECMODELFCT_NORM(jet, 16),                     // XSjet
  XSPECMODELFCT_C_NORM(C_kerrbb, 10),              // XSkerrbb
  XSPECMODELFCT_C_NORM(C_kerrd, 8),                // XSkerrd
  XSPECMODELFCT_C_NORM(C_spin, 10),                // XSkerrdisk

  XSPECMODELFCT_CON_F77(kyconv, 12),               // XSkyconv

  XSPECMODELFCT_NORM(kyrline, 12),                 // XSkyrline
  XSPECMODELFCT_C_NORM(C_laor, 6),                 // XSlaor
  XSPECMODELFCT_C_NORM(C_laor2, 8),                // XSlaor2
  XSPECMODELFCT_C_NORM(C_logpar, 4),               // XSlogpar
  XSPECMODELFCT_C_NORM(C_lorentzianLine, 3),       // XSlorentz
  XSPECMODELFCT_C_NORM(C_meka, 5),                 // XSmeka
  XSPECMODELFCT_C_NORM(C_mekal, 6),                // XSmekal
  XSPECMODELFCT_C_NORM(C_xsmkcf, 6),               // XSmkcflow
  XSPECMODELFCT_C_NORM(C_nei, 5),                  // XSnei
  XSPECMODELFCT_C_NORM(C_nlapec, 4),               // XSnlapec
  XSPECMODELFCT_C_NORM(C_npshock, 7),              // XSnpshock
  XSPECMODELFCT_NORM(nsa, 5),                      // XSnsa
  XSPECMODELFCT_NORM(nsagrav, 4),                  // XSnsagrav
  XSPECMODELFCT_NORM(nsatmos, 5),                  // XSnsatmos
  XSPECMODELFCT_C_NORM(C_nsmax, 4),                // XSnsmax
  XSPECMODELFCT_C_NORM(C_nsmaxg, 6),               // XSnsmaxg
  XSPECMODELFCT_C_NORM(C_nsx, 6),                  // XSnsx
  XSPECMODELFCT_C_NORM(C_xsnteea, 16),             // XSnteea
  XSPECMODELFCT_C_NORM(C_nthcomp, 6),              // XSnthComp
  XSPECMODELFCT_NORM(optxagn, 14),                 // XSoptxagn
  XSPECMODELFCT_NORM(optxagnf, 12),                // XSoptxagnf
  XSPECMODELFCT_NORM(xspegp, 4),                   // XSpegpwrlw
  XSPECMODELFCT_NORM(pexmon, 8),                   // XSpexmon
  XSPECMODELFCT_C_NORM(C_xspexrav, 8),             // XSpexrav
  XSPECMODELFCT_C_NORM(C_xspexriv, 10),            // XSpexriv
  XSPECMODELFCT_NORM(xsp1tr, 11),                  // XSplcabs
  XSPECMODELFCT_C_NORM(C_powerLaw, 2),             // XSpowerlaw
  XSPECMODELFCT_NORM(xsposm, 1),                   // XSposm
  XSPECMODELFCT_C_NORM(C_pshock, 6),               // XSpshock
  XSPECMODELFCT_NORM(qsosed, 7),                   // XSqsosed
  XSPECMODELFCT_C_NORM(C_raysmith, 4),             // XSraymond
  XSPECMODELFCT_NORM(xredge, 3),                   // XSredge
  XSPECMODELFCT_NORM(xsrefsch, 14),                // XSrefsch
  XSPECMODELFCT_C_NORM(C_rnei, 6),                 // XSrnei
  XSPECMODELFCT_C_NORM(C_sedov, 6),                // XSsedov
  XSPECMODELFCT_C_NORM(C_sirf, 10),                // XSsirf
  XSPECMODELFCT_C_NORM(slimbbmodel, 10),           // XSslimbh
  XSPECMODELFCT_C_NORM(C_snapec, 7),               // XSsnapec
  XSPECMODELFCT_NORM(srcut, 3),                    // XSsrcut
  XSPECMODELFCT_NORM(sresc, 3),                    // XSsresc
  XSPECMODELFCT_NORM(ssa, 3),                      // XSssa
  XSPECMODELFCT_NORM(xsstep, 3),                   // XSstep
  XSPECMODELFCT_C_NORM(C_tapec, 5),                // XStapec
  XSPECMODELFCT_C_NORM(C_vapec, 16),               // XSvapec
  XSPECMODELFCT_NORM(xsbrmv, 3),                   // XSvbremss
  XSPECMODELFCT_C_NORM(C_vcph, 18),                // XSvcph
  XSPECMODELFCT_C_NORM(C_vequil, 15),              // XSvequil
  XSPECMODELFCT_C_NORM(C_vgaussDem, 20),           // XSvgadem
  XSPECMODELFCT_C_NORM(C_vgnei, 18),               // XSvgnei
  XSPECMODELFCT_C_NORM(C_vmeka, 18),               // XSvmeka
  XSPECMODELFCT_C_NORM(C_vmekal, 19),              // XSvmekal
  XSPECMODELFCT_C_NORM(C_xsvmcf, 19),              // XSvmcflow
  XSPECMODELFCT_C_NORM(C_vnei, 17),                // XSvnei
  XSPECMODELFCT_C_NORM(C_vnpshock, 19),            // XSvnpshock
  XSPECMODELFCT_C_NORM(C_voigtLine, 4),            // XSvoigt
  XSPECMODELFCT_C_NORM(C_vpshock, 18),             // XSvpshock
  XSPECMODELFCT_C_NORM(C_vraysmith, 15),           // XSvraymond
  XSPECMODELFCT_C_NORM(C_vrnei, 18),               // XSvrnei
  XSPECMODELFCT_C_NORM(C_vsedov, 18),              // XSvsedov
  XSPECMODELFCT_C_NORM(C_vtapec, 17),              // XSvtapec
  XSPECMODELFCT_C_NORM(C_vvapec, 33),              // XSvvapec
  XSPECMODELFCT_C_NORM(C_vvgnei, 35),              // XSvvgnei
  XSPECMODELFCT_C_NORM(C_vvnei, 34),               // XSvvnei
  XSPECMODELFCT_C_NORM(C_vvnpshock, 36),           // XSvvnpshock
  XSPECMODELFCT_C_NORM(C_vvpshock, 35),            // XSvvpshock
  XSPECMODELFCT_C_NORM(C_vvrnei, 35),              // XSvvrnei
  XSPECMODELFCT_C_NORM(C_vvsedov, 35),             // XSvvsedov
  XSPECMODELFCT_C_NORM(C_vvtapec, 34),             // XSvvtapec
  XSPECMODELFCT_C_NORM(C_vvwDem, 37),              // XSvvwdem
  XSPECMODELFCT_C_NORM(C_vwDem, 21),               // XSvwdem
  XSPECMODELFCT_C_NORM(C_wDem, 8),                 // XSwdem
  XSPECMODELFCT_C_NORM(C_zagauss, 4),              // XSzagauss
  XSPECMODELFCT_NORM(xszbod, 3),                   // XSzbbody
  XSPECMODELFCT_C_NORM(C_zBrokenPowerLaw, 5),      // XSzbknpower
  XSPECMODELFCT_NORM(xszbrm, 3),                   // XSzbremss
  XSPECMODELFCT_C_NORM(C_zcutoffPowerLaw, 4),      // XSzcutoffpl
  XSPECMODELFCT_C_NORM(C_xszgau, 4),               // XSzgauss
  XSPECMODELFCT_C_NORM(C_zkerrbb, 10),             // XSzkerrbb
  XSPECMODELFCT_C_NORM(C_zLogpar, 5),              // XSzlogpar
  XSPECMODELFCT_C_NORM(C_zpowerLaw, 3),            // XSzpowerlw

  XSPECMODELFCT_C(C_xsabsori, 6),                  // XSabsori
  XSPECMODELFCT_C(C_acisabs, 8),                   // XSacisabs
  XSPECMODELFCT(xscnst, 1),                        // XSconstant
  XSPECMODELFCT(xscabs, 1),                        // XScabs
  XSPECMODELFCT(xscycl, 5),                        // XScyclabs
  XSPECMODELFCT(xsdust, 2),                        // XSdust
  XSPECMODELFCT(xsedge, 2),                        // XSedge
  XSPECMODELFCT(xsabsc, 1),                        // XSexpabs
  XSPECMODELFCT(xsexp, 3),                         // XSexpfac
  XSPECMODELFCT_C(C_gaussianAbsorptionLine, 3),    // XSgabs
  XSPECMODELFCT(xsphei, 3),                        // XSheilin
  XSPECMODELFCT(xshecu, 2),                        // XShighecut
  XSPECMODELFCT(xshrfl, 8),                        // XShrefl
#ifdef XSPEC_12_12_1
  XSPECMODELFCT_DBL(ismabs, 31),                   // XSismabs
  XSPECMODELFCT_DBL(ismdust, 3),                   // XSismdust
#else
  XSPECMODELFCT(ismabs, 31),                       // XSismabs
  XSPECMODELFCT(ismdust, 3),                       // XSismdust
#endif
  XSPECMODELFCT_C(C_logconst, 1),                  // XSlogconst
  XSPECMODELFCT_C(C_log10con, 1),                  // XSlog10con
  XSPECMODELFCT(xslyman, 4),                       // XSlyman
  XSPECMODELFCT(xsntch, 3),                        // XSnotch
#ifdef XSPEC_12_12_1
  XSPECMODELFCT_DBL(olivineabs, 2),                // XSolivineabs
#else
  XSPECMODELFCT(olivineabs, 2),                    // XSolivineabs
#endif
  XSPECMODELFCT(xsabsp, 2),                        // XSpcfabs
  XSPECMODELFCT(xsphab, 1),                        // XSphabs
  XSPECMODELFCT(xsplab, 2),                        // XSplabs
  XSPECMODELFCT_C(C_xspwab, 3),                    // XSpwab
  XSPECMODELFCT(xscred, 1),                        // XSredden
  XSPECMODELFCT(xssmdg, 4),                        // XSsmedge
  XSPECMODELFCT_C(C_superExpCutoff, 2),            // XSspexpcut
  XSPECMODELFCT(xsspln, 6),                        // XSspline
  XSPECMODELFCT(xssssi, 1),                        // XSSSS_ice
  XSPECMODELFCT_C(C_swind1, 4),                    // XSswind1
  XSPECMODELFCT_C(C_tbabs, 1),                     // XSTBabs
  XSPECMODELFCT_C(C_tbfeo, 4),                     // XSTBfeo
  XSPECMODELFCT_C(C_tbgas, 2),                     // XSTBgas
  XSPECMODELFCT_C(C_tbgrain, 6),                   // XSTBgrain
  XSPECMODELFCT_C(C_tbvabs, 42),                   // XSTBvarabs
  XSPECMODELFCT_C(C_tbpcf, 3),                     // XSTBpcf
  XSPECMODELFCT_C(C_tbrel, 42),                    // XSTBrel
  XSPECMODELFCT(xsred, 1),                         // XSuvred
  XSPECMODELFCT(xsabsv, 18),                       // XSvarabs
  XSPECMODELFCT(xsvphb, 18),                       // XSvphabs
  XSPECMODELFCT(xsabsw, 1),                        // XSwabs
  XSPECMODELFCT(xswnab, 2),                        // XSwndabs
  XSPECMODELFCT(xsxirf, 13),                       // XSxion
  XSPECMODELFCT_C(C_xscatmodel, 4),                // XSxscat
  XSPECMODELFCT_C(xszbabs, 4),                     // XSzbabs
  XSPECMODELFCT(mszdst, 4),                        // XSzdust
  XSPECMODELFCT(xszedg, 3),                        // XSzedge
  XSPECMODELFCT(xszhcu, 3),                        // XSzhighect
  XSPECMODELFCT(zigm, 3),                          // XSzigm
  XSPECMODELFCT(xszabp, 3),                        // XSzpcfabs
  XSPECMODELFCT(xszphb, 2),                        // XSzphabs
  XSPECMODELFCT(zxipab, 5),                        // XSzxipab
  XSPECMODELFCT_C(C_zxipcf, 4),                    // XSzxipcf
  XSPECMODELFCT(xszcrd, 2),                        // XSzredden
  XSPECMODELFCT(msldst, 4),                        // XSzsmdust
  XSPECMODELFCT_C(C_ztbabs, 2),                    // XSzTBabs
  XSPECMODELFCT(xszvab, 19),                       // XSzvarabs
  XSPECMODELFCT(xszvfe, 5),                        // XSzvfeabs
  XSPECMODELFCT(xszvph, 19),                       // XSzvphabs
  XSPECMODELFCT(xszabs, 2),                        // XSzwabs
  XSPECMODELFCT(xszwnb, 3),                        // XSzwndabs

  XSPECMODELFCT_CON(C_cflux, 3),                   // XScflux
  XSPECMODELFCT_CON(C_clumin, 4),                  // XSclumin
#ifdef XSPEC_12_13_0
  XSPECMODELFCT_CON(C_cglumin, 4),                 // XScglumin
#endif
  XSPECMODELFCT_CON(C_cpflux, 3),                  // XScpflux
  XSPECMODELFCT_CON(C_gsmooth, 2),                 // XSgsmooth
  XSPECMODELFCT_CON(C_ireflct, 7),                 // XSireflect
  XSPECMODELFCT_CON(C_kdblur, 4),                  // XSkdblur
  XSPECMODELFCT_CON(C_kdblur2, 6),                 // XSkdblur2
  XSPECMODELFCT_CON(C_spinconv, 7),                // XSkerrconv
  XSPECMODELFCT_CON(C_lsmooth, 2),                 // XSlsmooth
  XSPECMODELFCT_CON(C_PartialCovering, 1),         // XSpartcov
  XSPECMODELFCT_CON(C_rdblur, 4),                  // XSrdblur
  XSPECMODELFCT_CON(C_reflct, 5),                  // XSreflect
  XSPECMODELFCT_CON(C_rfxconv, 5),                 // XSrfxconv
  XSPECMODELFCT_CON_F77(rgsxsrc, 1),               // XSrgsxsrc
  XSPECMODELFCT_CON(C_simpl, 3),                   // XSsimpl
  XSPECMODELFCT_CON_F77(thcompf, 4),               // XSthcomp
  XSPECMODELFCT_CON(C_vashift, 1),                 // XSvashift
  XSPECMODELFCT_CON(C_vmshift, 1),                 // XSvmshift
  XSPECMODELFCT_CON(C_xilconv, 6),                 // XSxilconv
  XSPECMODELFCT_CON(C_zashift, 1),                 // XSzashift
  XSPECMODELFCT_CON(C_zmshift, 1),                 // XSzmshift

  XSPECMODELFCT_C_NORM(beckerwolff, 13),           // XSbwcycl

  // End model definitions

  // XSPEC table models
  KWSPEC(tabint, sherpa::astro::xspec::xspectablemodel),

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
  // Ensure the XSPEC library is initialized.
  if ( EXIT_SUCCESS != _sherpa_init_xspec_library() )
    return NULL;

  import_array();
  return PyModule_Create(&xspec_module);
}
