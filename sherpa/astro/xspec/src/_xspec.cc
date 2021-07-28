//  Copyright (C) 2007, 2015-2018, 2019, 2020, 2021, 2022, 2023
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

#include <XSFunctions/Utilities/xsFortran.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSutility.h>

// funcWrappers: C_<model> are declared here; the other models are defined in
// functionMap.h
//
#include <XSFunctions/funcWrappers.h>
#include <XSFunctions/functionMap.h>

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

  // Redirect the stdout channel for the duration of the FNINIT call.
  // The replacement is set to ostringstream so that we can find out
  // what was sent to it in case of emergency (but would need code
  // changes to do this). This used to be done by manually replacing
  // the std::cout buffer.
  //
  // Note that only FNINIT needs this, as the FunctionUtility calls do
  // not create string output (at least for 12.12.x/12.13.0).
  //
  std::ostream* outStream = IosHolder::outHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			&tmpStream,
			IosHolder::errHolder());

  try {
    // Initialize XSPEC model library
    FNINIT();

  } catch(...) {
    IosHolder::setStreams(IosHolder::inHolder(),
			  outStream,
			  IosHolder::errHolder());

    // Raise appropriate error message that XSPEC initialization failed.
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "check HEADAS environment variable" );
    return EXIT_FAILURE;
  }

  IosHolder::setStreams(IosHolder::inHolder(),
			outStream,
			IosHolder::errHolder());

  // We used to set the chatter to 0 but this misses useful info
  // (like XSPEC can't find the data files) that would be reported
  // by XSPEC, so use the default XSPEC setting. It does not appear
  // that this value is read in from ~/.xspec/Xspec.init so set
  // it here so we have repeatable behavior.
  //
  FunctionUtility::xwriteChatter( 10 );

  // Set cosmology initial values to XSPEC initial values
  FunctionUtility::setH0( 70.0 );
  FunctionUtility::setq0( 0.0 );
  FunctionUtility::setlambda0( 0.73 );

  init = true;
  return EXIT_SUCCESS;

}

static PyObject* get_chatter( PyObject *self )
{
  return Py_BuildValue( (char*)"i", FunctionUtility::xwriteChatter() );
}


// TODO: we could send in an integer for the Z number (ie either name or number)
static PyObject* get_abund( PyObject *self, PyObject *args )
{

  char* element = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"|s", &element ) )
    return NULL;

  // Not asked for an element so return the table name
  if ( !element ) {
    return (PyObject*) Py_BuildValue( (char*)"s", FunctionUtility::ABUND().c_str() );
  }

  // Get the specific abundance. Unfortunately getAbundance reports an
  // error to stderr when an invalid element is used, so we need to
  // hide this.
  //
  std::ostream* errStream = IosHolder::errHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			&tmpStream);

  float abundVal = FunctionUtility::getAbundance(string(element));

  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			errStream);

  // Was there an error?
  //
  if( tmpStream.str().size() > 0 ) {
    return PyErr_Format( PyExc_TypeError, // TODO: change from TypeError to ValueError?
			 (char*)"could not find element '%s'", element);
  }

  return (PyObject*) Py_BuildValue( (char*)"f", abundVal );

}


static PyObject* get_cosmo( PyObject *self )
{
  // Assume these can not throw errors
  float h0 = FunctionUtility::getH0();
  float l0 = FunctionUtility::getlambda0();
  float q0 = FunctionUtility::getq0();

  return Py_BuildValue( (char*)"fff", h0, q0, l0 );

}


static PyObject* set_chatter( PyObject *self, PyObject *args )
{

  int chatter = 0;

  if ( !PyArg_ParseTuple( args, (char*)"i", &chatter ) )
    return NULL;

  FunctionUtility::xwriteChatter(chatter);
  Py_RETURN_NONE;

}


// Based on xsFortran::FPSOLR
//
// TODO: add a version where we can send in an array of numbers
//
static PyObject* set_abund( PyObject *self, PyObject *args )
{

  char* table = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"s", &table ) )
    return NULL;

  string tableName = string(table);
  tableName = XSutility::lowerCase(tableName);

  if (tableName == "file") {
    FunctionUtility::ABUND(tableName);
    Py_RETURN_NONE;
  }

  if (FunctionUtility::checkAbund(tableName)) {
    FunctionUtility::ABUND(tableName);
    Py_RETURN_NONE;
  }

  // If we've got here then try to read the data from a file
  //
  const size_t nelems = FunctionUtility::NELEMS();
  std::vector<float> vals(nelems, 0);
  size_t count(0);

  std::ifstream fileStream(table);

  try {
    float element;
    fileStream.exceptions(std::ios_base::failbit);

    while (count < nelems && fileStream >> element) {
      vals[count] = element;
      ++count;
    }
  }
  catch ( std::exception& ) {

    // We do not error out if it was an eofbit failure, as we assume
    // this just means that the file contains < NELEMS elements, and
    // the missing elements will be set to 0 abundance. We could try
    // and explicitly handle this case from the exception, but let's
    // follow XSPEC and check the eof status.
    //
    if (!fileStream.eof()) {
      return PyErr_Format( PyExc_ValueError,
			   (char*)"Cannot read file '%s'.  It may not exist or contains invalid data",
			   table);
    }
  }

  FunctionUtility::ABUND("file");
  FunctionUtility::abundanceVectors("file", vals);

  Py_RETURN_NONE;

}


static PyObject* set_cosmo( PyObject *self, PyObject *args )
{

  float h0;
  float l0;
  float q0;

  if ( !PyArg_ParseTuple( args, (char*)"fff", &h0, &q0, &l0 ) )
    return NULL;

  FunctionUtility::setH0(h0);
  FunctionUtility::setq0(q0);
  FunctionUtility::setlambda0(l0);
  Py_RETURN_NONE;

}


static PyObject* set_cross( PyObject *self, PyObject *args )
{

  char* csection = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s", &csection ) )
    return NULL;

  string tableName = string(csection);
  tableName = XSutility::lowerCase(tableName);

  // On failure, checkXsect catches a YellowAlert which creates output
  // to the error channel, so we over-ride it for this call.
  //
  std::ostream* errStream = IosHolder::errHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			&tmpStream);
  const bool known = FunctionUtility::checkXsect(tableName);
  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			errStream);

  if (!known) {
    return PyErr_Format( PyExc_ValueError,
			 (char*)"could not set XSPEC photoelectric "
			 "cross-section to '%s'",
			 csection);
  }

  FunctionUtility::XSECT(tableName);
  Py_RETURN_NONE;

}


// TODO: We could have a seperate "reset" command
//
static PyObject* set_xset( PyObject *self, PyObject *args )
{

  char* str_name = NULL;
  char* str_value = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"ss", &str_name, &str_value ) )
    return NULL;

  string name = XSutility::upperCase(string(str_name));
  if (name == "INITIALIZE") {
    FunctionUtility::eraseModelStringDataBase();
  } else {
    FunctionUtility::setModelString(name, string(str_value));
  }
  Py_RETURN_NONE;

}

static PyObject* get_xset( PyObject *self, PyObject *args  )
{

  char* str_name = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"s", &str_name ) )
    return NULL;

  static string value;
  value = FunctionUtility::getModelString(string(str_name));
  if (value == FunctionUtility::NOT_A_KEY()) {
    value.erase();
  }

  return Py_BuildValue( (char*)"s", value.c_str() );

}

template <const std::string& get()>
static PyObject* get_xspec_string( PyObject *self ) {
  return Py_BuildValue( (char*)"s", get().c_str() );
}

static PyObject* set_manager_data_path( PyObject *self, PyObject *args )
{

  char* path = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"s", &path ) )
    return NULL;

  FunctionUtility::managerPath(string(path));
  Py_RETURN_NONE;

}

#define NOARGSPEC(name, func) \
  { (char *)#name, (PyCFunction)func, METH_NOARGS, NULL }


static PyMethodDef XSpecMethods[] = {
  // Utility routines
  //
  NOARGSPEC(get_xsversion, get_xspec_string<XSutility::xs_version>),
  NOARGSPEC(get_xschatter, get_chatter),
  FCTSPEC(set_xschatter, set_chatter),
  FCTSPEC(get_xsabund, get_abund),
  FCTSPEC(set_xsabund, set_abund),
  FCTSPEC(set_xscosmo, set_cosmo),
  NOARGSPEC(get_xscosmo, get_cosmo),
  NOARGSPEC(get_xsxsect, get_xspec_string<FunctionUtility::XSECT>),

  FCTSPEC(set_xsxsect, set_cross),
  FCTSPEC(set_xsxset, set_xset),
  FCTSPEC(get_xsxset, get_xset),
  NOARGSPEC(get_xspath_manager, get_xspec_string<FunctionUtility::managerPath>),
  NOARGSPEC(get_xspath_model, get_xspec_string<FunctionUtility::modelDataPath>),
  FCTSPEC(set_xspath_manager, set_manager_data_path),

  // Models
  //
  XSPECMODELFCT_NORM( agnsed, 16 ),
  XSPECMODELFCT_NORM( qsosed, 7 ),

  XSPECMODELFCT_NORM( agnslim, 15 ),

  XSPECMODELFCT_C_NORM( C_apec, 4 ),
  XSPECMODELFCT_C_NORM( C_bapec, 5 ),
  XSPECMODELFCT_NORM( xsblbd, 2 ),
  XSPECMODELFCT_NORM( xsbbrd, 2 ),
  XSPECMODELFCT_C_NORM( C_xsbexrav, 10 ),
  XSPECMODELFCT_C_NORM( C_xsbexriv, 12 ),
  XSPECMODELFCT_C_NORM( C_brokenPowerLaw, 4 ),
  XSPECMODELFCT_C_NORM( C_broken2PowerLaw, 6 ),
  XSPECMODELFCT_C_NORM( C_sirf, 10 ),
  XSPECMODELFCT_NORM( xsbmc, 4 ),
  XSPECMODELFCT_NORM( xsbrms, 2 ),
  XSPECMODELFCT_C_NORM( C_brnei, 7 ),
  XSPECMODELFCT_C_NORM( C_bvapec, 17 ),
  XSPECMODELFCT_C_NORM( C_bvrnei, 19 ),
  XSPECMODELFCT_C_NORM( C_c6mekl, 11 ),
  XSPECMODELFCT_C_NORM( C_c6pmekl, 11 ),
  XSPECMODELFCT_C_NORM( C_c6pvmkl, 24 ),
  XSPECMODELFCT_C_NORM( C_c6vmekl, 24 ),
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
  XSPECMODELFCT_C_NORM( C_diskline, 6 ),
  XSPECMODELFCT_NORM( diskm, 5 ),
  XSPECMODELFCT_NORM( disko, 5 ),
  XSPECMODELFCT_NORM( diskpbb, 3 ),
  XSPECMODELFCT_NORM( xsdiskpn, 3 ),
  XSPECMODELFCT_C_NORM( C_equil, 4 ),
  XSPECMODELFCT_NORM( xsxpdec, 2 ),
  XSPECMODELFCT_NORM( ezdiskbb, 2 ),
  XSPECMODELFCT_C_NORM( C_gaussianLine, 3 ),
  XSPECMODELFCT_C_NORM( C_gnei, 6 ),
  XSPECMODELFCT_NORM( grad, 7 ),
  XSPECMODELFCT_C_NORM( xsgrbcomp, 10 ),
  XSPECMODELFCT_NORM( xsgrbm, 4 ),
  XSPECMODELFCT_C_NORM( C_kerrbb, 10 ),

  XSPECMODELFCT_C_NORM( C_zkerrbb, 10 ),

  XSPECMODELFCT_C_NORM( C_kerrd, 8 ),
  XSPECMODELFCT_C_NORM( C_spin, 10 ),

  XSPECMODELFCT_NORM( kyrline, 12 ),

  XSPECMODELFCT_C_NORM( C_laor, 6 ),
  XSPECMODELFCT_C_NORM( C_laor2, 8 ),
  XSPECMODELFCT_C_NORM( C_lorentzianLine, 3 ),
  XSPECMODELFCT_C_NORM( C_meka, 5 ),
  XSPECMODELFCT_C_NORM( C_mekal, 6 ),
  XSPECMODELFCT_C_NORM( C_xsmkcf, 6 ),
  XSPECMODELFCT_C_NORM( C_nei, 5 ),
  XSPECMODELFCT_C_NORM( C_nlapec, 4 ),
  XSPECMODELFCT_C_NORM( C_npshock, 7 ),
  XSPECMODELFCT_NORM( nsa, 5 ),
  XSPECMODELFCT_NORM( nsagrav, 4 ),
  XSPECMODELFCT_NORM( nsatmos, 5 ),
  XSPECMODELFCT_C_NORM( C_nsmax, 4 ),
  XSPECMODELFCT_C_NORM( C_xsnteea, 16 ),
  XSPECMODELFCT_C_NORM( C_nthcomp, 6 ),
  XSPECMODELFCT_NORM( xspegp, 4 ),
  XSPECMODELFCT_C_NORM( C_xspexrav, 8 ),
  XSPECMODELFCT_C_NORM( C_xspexriv, 10 ),
  XSPECMODELFCT_NORM( xsp1tr, 11 ),
  XSPECMODELFCT_C_NORM( C_powerLaw, 2 ),
  XSPECMODELFCT_NORM( xsposm, 1 ),
  XSPECMODELFCT_C_NORM( C_pshock, 6 ),
  XSPECMODELFCT_C_NORM(C_raysmith, 4 ),
  XSPECMODELFCT_NORM( xredge, 3 ),
  XSPECMODELFCT_NORM( xsrefsch, 14 ),

  XSPECMODELFCT_C_NORM( C_sedov, 6 ),
  XSPECMODELFCT_NORM( srcut, 3 ),
  XSPECMODELFCT_NORM( sresc, 3 ),
  XSPECMODELFCT_NORM( ssa, 3 ),
  XSPECMODELFCT_NORM( xsstep, 3 ),

  XSPECMODELFCT_C_NORM( C_vapec, 16 ),
  XSPECMODELFCT_NORM( xsbrmv, 3 ),

  XSPECMODELFCT_C_NORM( C_vcph, 18 ),

  XSPECMODELFCT_C_NORM( C_vequil, 15 ),
  XSPECMODELFCT_C_NORM( C_vgnei, 18 ),
  XSPECMODELFCT_C_NORM( C_vmeka, 18 ),
  XSPECMODELFCT_C_NORM( C_vmekal, 19 ),
  XSPECMODELFCT_C_NORM( C_xsvmcf, 19 ),
  XSPECMODELFCT_C_NORM( C_vnei, 17 ),
  XSPECMODELFCT_C_NORM( C_vnpshock, 19 ),
  XSPECMODELFCT_C_NORM( C_vpshock, 18 ),
  XSPECMODELFCT_C_NORM( C_vraysmith, 15 ),
  XSPECMODELFCT_C_NORM( C_vsedov, 18 ),
  XSPECMODELFCT_NORM( xszbod, 3 ),

  XSPECMODELFCT_C_NORM( C_zBrokenPowerLaw, 5 ),

  XSPECMODELFCT_NORM( xszbrm, 3 ),
  XSPECMODELFCT_C_NORM( C_zcutoffPowerLaw, 4),
  XSPECMODELFCT_C_NORM( C_xszgau, 4 ),

  XSPECMODELFCT_C_NORM( C_zLogpar, 5 ),

  XSPECMODELFCT_C_NORM( C_zpowerLaw, 3 ),
  XSPECMODELFCT_C( C_xsabsori, 6 ),

  XSPECMODELFCT_C( C_acisabs, 8 ),

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

  XSPECMODELFCT_C( C_swind1, 4 ),

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

  XSPECMODELFCT_C( C_zxipcf, 4 ),

  XSPECMODELFCT( xszcrd, 2 ),
  XSPECMODELFCT( msldst, 4 ),
  XSPECMODELFCT_C( C_ztbabs, 2 ),
  XSPECMODELFCT( xszvab, 19 ),
  XSPECMODELFCT( xszvfe, 5 ),
  XSPECMODELFCT( xszvph, 19 ),
  XSPECMODELFCT( xszabs, 2 ),
  XSPECMODELFCT( xszwnb, 3 ),

  XSPECMODELFCT_C_NORM( C_cph, 5 ),

  XSPECMODELFCT_C_NORM( C_cplinear, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqpair, 21 ),
  XSPECMODELFCT_C_NORM( C_xseqth, 21 ),
  XSPECMODELFCT_C_NORM( C_xscompth, 21 ),
  XSPECMODELFCT_C_NORM( C_bvvapec, 34 ),
  XSPECMODELFCT_C_NORM( C_vvapec, 33 ),
  XSPECMODELFCT_C_NORM( C_bvvrnei, 36 ),
  XSPECMODELFCT( zigm, 3 ),
  // New XSPEC 12.7.1 models
  XSPECMODELFCT_C_NORM( C_gaussDem, 7 ),
  XSPECMODELFCT_C_NORM( C_vgaussDem, 20 ),
  XSPECMODELFCT_NORM( eplogpar, 3 ),
  XSPECMODELFCT_C_NORM( C_logpar, 4 ),
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

  XSPECMODELFCT_C_NORM( C_nsmaxg, 6 ),
  XSPECMODELFCT_C_NORM( C_nsx, 6 ),

  XSPECMODELFCT_C_NORM( C_rnei, 6 ),
  XSPECMODELFCT_C_NORM( C_vrnei, 18 ),
  XSPECMODELFCT_C_NORM( C_vvrnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvgnei, 35 ),
  XSPECMODELFCT_C_NORM( C_vvnei, 34 ),
  XSPECMODELFCT_C_NORM( C_vvnpshock, 36 ),
  XSPECMODELFCT_C_NORM( C_vvpshock, 35 ),
  XSPECMODELFCT_C_NORM( C_vvsedov, 35 ),

  // We do not have a direct interface to the c_func routines, so
  // take advantage of the fact XSPEC provides multiple APIs
  // and use the C one.
  // XSPECMODELFCT_C_NORM( c_beckerwolff, 13 ),
  XSPECMODELFCT_C_NORM( beckerwolff, 13 ),

  //multiplicative
  XSPECMODELFCT( xsphei, 3 ),
  XSPECMODELFCT( xslyman, 4 ),
  XSPECMODELFCT_C( xszbabs, 4 ), // DJB thinks it's okay to use the C++ wrapper for C

  // XSPEC table models
  KWSPEC(tabint, sherpa::astro::xspec::xspectablemodel),

  // XSPEC convolution models
  XSPECMODELFCT_CON(C_cflux, 3),
  XSPECMODELFCT_CON(C_cpflux, 3),

  XSPECMODELFCT_CON(C_gsmooth, 2),

  XSPECMODELFCT_CON(C_ireflct, 7),
  XSPECMODELFCT_CON(C_kdblur, 4),
  XSPECMODELFCT_CON(C_kdblur2, 6),
  XSPECMODELFCT_CON(C_spinconv, 7),

  XSPECMODELFCT_CON_F77(kyconv, 12),

  XSPECMODELFCT_CON(C_lsmooth, 2),

  XSPECMODELFCT_CON(C_PartialCovering, 1),
  XSPECMODELFCT_CON(C_rdblur, 4),
  XSPECMODELFCT_CON(C_reflct, 5),

  XSPECMODELFCT_CON_F77(rgsxsrc, 1),
  XSPECMODELFCT_CON(C_simpl, 3),
  XSPECMODELFCT_CON(C_zashift, 1),
  XSPECMODELFCT_CON(C_zmshift, 1),

  XSPECMODELFCT_CON_F77(thcompf, 4),

  // Models from 12.9.1
  //
  //
  XSPECMODELFCT_C_NORM(C_btapec, 6),
  XSPECMODELFCT_C_NORM(C_bvtapec, 18),
  XSPECMODELFCT_C_NORM(C_bvvtapec, 35),
  XSPECMODELFCT_C_NORM(C_tapec, 5),
  XSPECMODELFCT_C_NORM(C_vtapec, 17),
  XSPECMODELFCT_C_NORM(C_vvtapec, 34),

  XSPECMODELFCT_C_NORM(C_carbatm, 4),
  XSPECMODELFCT_C_NORM(C_hatm, 4),

#ifdef XSPEC_12_12_1
  XSPECMODELFCT_DBL(ismabs, 31),
#else
  XSPECMODELFCT(ismabs, 31),
#endif

  XSPECMODELFCT_C_NORM(slimbbmodel, 10),
  XSPECMODELFCT_C_NORM(C_snapec, 7),
  XSPECMODELFCT_C(C_tbfeo, 4),
  XSPECMODELFCT_C(C_tbgas, 2),
  XSPECMODELFCT_C(C_tbpcf, 3),
  XSPECMODELFCT_C(C_tbrel, 42),
  XSPECMODELFCT_C_NORM(C_voigtLine, 4),
  XSPECMODELFCT_C(C_xscatmodel, 4),

  XSPECMODELFCT_CON(C_clumin, 4),
#ifdef XSPEC_12_13_0
  XSPECMODELFCT_CON(C_cglumin, 4),
#endif
  XSPECMODELFCT_CON(C_rfxconv, 5),
  XSPECMODELFCT_CON(C_vashift, 1),
  XSPECMODELFCT_CON(C_vmshift, 1),
  XSPECMODELFCT_CON(C_xilconv, 6),

  XSPECMODELFCT_NORM(jet, 16),

#ifdef XSPEC_12_12_1
  XSPECMODELFCT_DBL(ismdust, 3),
  XSPECMODELFCT_DBL(olivineabs, 2),
#else
  XSPECMODELFCT(ismdust, 3),
  XSPECMODELFCT(olivineabs, 2),
#endif // XSPEC_12_12_1

  XSPECMODELFCT_C(C_logconst, 1),
  XSPECMODELFCT_C(C_log10con, 1),

  XSPECMODELFCT_C_NORM( xsgrbjet, 14 ),  // follow xsgrbcomp and drop the leading c_
  XSPECMODELFCT_C_NORM( C_vvwDem, 37 ),
  XSPECMODELFCT_C_NORM( C_vwDem, 21 ),
  XSPECMODELFCT_C_NORM( C_wDem, 8 ),

  XSPECMODELFCT(zxipab, 5),

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
