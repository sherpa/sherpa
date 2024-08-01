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

// Documentation for the "XSPEC internal functions" is at
// https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/internal/XspecInternalFunctionsGuide.html
//
// xsfortran is only needed to support FNINIT as there's (as of XSPEC 12.13.0)
// no version of this functionality in FunctionUtility.
//
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
  //
  std::ostream* outStream = IosHolder::outHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			&tmpStream,
			IosHolder::errHolder());

  try {
    // Initialize XSPEC model library.
    FNINIT();

  } catch(...) {
    IosHolder::setStreams(IosHolder::inHolder(),
			  outStream,
			  IosHolder::errHolder());

    // The contents of tmpStream could be inspected to see if it
    // contains useful information for the user, but at this point of
    // the initialization it is not obvious that it would provide any
    // extra information.
    //
    PyErr_SetString( PyExc_ImportError,
		     (char*)"XSPEC initialization failed; "
		     "check HEADAS environment variable" );
    return EXIT_FAILURE;
  }

  IosHolder::setStreams(IosHolder::inHolder(),
			outStream,
			IosHolder::errHolder());

  // Set a number of values to their XSPEC defaults (as of XSPEC
  // 12.14.1, but they have not changed for a long time). It appears
  // these are not read in from the user's ~/.xspec/Xspec.init file.
  //
  FunctionUtility::xwriteChatter( 10 );
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


// Return the "nice" name for an abundance table - e.g.
// "Anders E. & Grevesse N. Geochimica et Cosmochimica Acta 53, 197 (1989)"
// for "angr".
//
static PyObject* get_abund_doc( PyObject *self, PyObject *args )
{
  char *name = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"s", &name ) )
    return NULL;

  std::string doc = FunctionUtility::abundDoc(std::string(name));

  return Py_BuildValue( (char*)"s", doc.c_str() );
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
  FCTSPEC(get_xsabund_doc, get_abund_doc),
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

  // Start model definitions

  XSPECMODELFCT_C_NORM(C_agauss, 3),               // XSagauss
  XSPECMODELFCT_NORM(agnsed, 16),                  // XSagnsed
  XSPECMODELFCT_NORM(agnslim, 15),                 // XSagnslim
  XSPECMODELFCT_C_NORM(C_apec, 4),                 // XSapec
  XSPECMODELFCT_C_NORM(C_bapec, 5),                // XSbapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bcempow, 8),              // XSbcempow
  XSPECMODELFCT_C_NORM(C_bcheb6, 12),              // XSbcheb6
  XSPECMODELFCT_C_NORM(C_bcie, 6),                 // XSbcie
  XSPECMODELFCT_C_NORM(C_bcoolflow, 7),            // XSbcoolflow
  XSPECMODELFCT_C_NORM(C_bcph, 6),                 // XSbcph
  XSPECMODELFCT_C_NORM(C_bequil, 5),               // XSbequil
  XSPECMODELFCT_C_NORM(C_bexpcheb6, 12),           // XSbexpcheb6
  XSPECMODELFCT_C_NORM(C_bgaussDem, 8),            // XSbgadem
  XSPECMODELFCT_C_NORM(C_bgnei, 7),                // XSbgnei
  XSPECMODELFCT_C_NORM(C_bnei, 6),                 // XSbnei
  XSPECMODELFCT_C_NORM(C_bsnapec, 8),              // XSbsnapec
#endif
  XSPECMODELFCT_C_NORM(C_btapec, 6),               // XSbtapec
  XSPECMODELFCT_NORM(xsblbd, 2),                   // XSbbody
  XSPECMODELFCT_NORM(xsbbrd, 2),                   // XSbbodyrad
  XSPECMODELFCT_C_NORM(C_xsbexrav, 10),            // XSbexrav
  XSPECMODELFCT_C_NORM(C_xsbexriv, 12),            // XSbexriv
  XSPECMODELFCT_C_NORM(C_brokenPowerLaw, 4),       // XSbknpower
  XSPECMODELFCT_C_NORM(C_broken2PowerLaw, 6),      // XSbkn2pow
  XSPECMODELFCT_NORM(xsbmc, 4),                    // XSbmc
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bnpshock, 8),             // XSbnpshock
  XSPECMODELFCT_C_NORM(C_bpshock, 7),              // XSbpshock
#endif
  XSPECMODELFCT_NORM(xsbrms, 2),                   // XSbremss
  XSPECMODELFCT_C_NORM(C_brnei, 7),                // XSbrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bsedov, 7),               // XSbsedov
#endif
  XSPECMODELFCT_C_NORM(C_bvapec, 17),              // XSbvapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bvcempow, 21),            // XSbvcempow
  XSPECMODELFCT_C_NORM(C_bvcheb6, 25),             // XSbvcheb6
  XSPECMODELFCT_C_NORM(C_bvcie, 18),               // XSbvcie
  XSPECMODELFCT_C_NORM(C_bvcoolflow, 20),          // XSbvcoolflow
  XSPECMODELFCT_C_NORM(C_bvcph, 19),               // XSbvcph
  XSPECMODELFCT_C_NORM(C_bvequil, 16),             // XSbvequil
  XSPECMODELFCT_C_NORM(C_bvexpcheb6, 25),          // XSbvexpcheb6
  XSPECMODELFCT_C_NORM(C_bvgaussDem, 21),          // XSbvgadem
  XSPECMODELFCT_C_NORM(C_bvgnei, 19),              // XSbvgnei
  XSPECMODELFCT_C_NORM(C_bvnei, 18),               // XSbvnei
  XSPECMODELFCT_C_NORM(C_bvnpshock, 20),           // XSbvnpshock
  XSPECMODELFCT_C_NORM(C_bvpshock, 19),            // XSbvpshock
#endif
  XSPECMODELFCT_C_NORM(C_bvrnei, 19),              // XSbvrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bvsedov, 19),             // XSbvsedov
#endif
  XSPECMODELFCT_C_NORM(C_bvtapec, 18),             // XSbvtapec
  XSPECMODELFCT_C_NORM(C_bvvapec, 34),             // XSbvvapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bvvcie, 35),              // XSbvvcie
  XSPECMODELFCT_C_NORM(C_bvvgaussDem, 37),         // XSbvvgadem
  XSPECMODELFCT_C_NORM(C_bvvgnei, 36),             // XSbvvgnei
  XSPECMODELFCT_C_NORM(C_bvvnei, 35),              // XSbvvnei
  XSPECMODELFCT_C_NORM(C_bvvnpshock, 37),          // XSbvvnpshock
  XSPECMODELFCT_C_NORM(C_bvvpshock, 36),           // XSbvvpshock
#endif
  XSPECMODELFCT_C_NORM(C_bvvrnei, 36),             // XSbvvrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bvvsedov, 36),            // XSbvvsedov
#endif
  XSPECMODELFCT_C_NORM(C_bvvtapec, 35),            // XSbvvtapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_bvvwDem, 38),             // XSbvvwdem
  XSPECMODELFCT_C_NORM(C_bvwDem, 22),              // XSbvwdem
  XSPECMODELFCT_C_NORM(C_bwDem, 9),                // XSbwdem
#endif
  XSPECMODELFCT_C_NORM(C_c6mekl, 11),              // XSc6mekl
  XSPECMODELFCT_C_NORM(C_c6pmekl, 11),             // XSc6pmekl
  XSPECMODELFCT_C_NORM(C_c6pvmkl, 24),             // XSc6pvmkl
  XSPECMODELFCT_C_NORM(C_c6vmekl, 24),             // XSc6vmekl
  XSPECMODELFCT_C_NORM(C_carbatm, 4),              // XScarbatm
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_cemMekal, 7),             // XScemekl
  XSPECMODELFCT_C_NORM(C_cempow, 7),               // XScempow
#else
  XSPECMODELFCT_NORM(cemekl, 7),                   // XScemekl
#endif
  XSPECMODELFCT_C_NORM(C_cemVMekal, 20),           // XScevmkl
  XSPECMODELFCT_C_NORM(C_xscflw, 6),               // XScflow
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_cheb6, 11),               // XScheb6
  XSPECMODELFCT_C_NORM(C_cie, 5),                  // XScie
#endif
  XSPECMODELFCT_NORM(compbb, 4),                   // XScompbb
  XSPECMODELFCT_C_NORM(xscompmag, 9),              // XScompmag
  XSPECMODELFCT_NORM(compls, 3),                   // XScompLS
  XSPECMODELFCT_C_NORM(C_xscompps, 20),            // XScompPS
  XSPECMODELFCT_NORM(compst, 3),                   // XScompST
  XSPECMODELFCT_C_NORM(xscomptb, 7),               // XScomptb
  XSPECMODELFCT_C_NORM(C_xscompth, 21),            // XScompth
  XSPECMODELFCT_NORM(xstitg, 6),                   // XScompTT
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_coolflow, 6),             // XScoolflow
#endif
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
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_eebremss, 4),             // XSeebremss
#endif
  XSPECMODELFCT_NORM(eplogpar, 3),                 // XSeplogpar
  XSPECMODELFCT_C_NORM(C_xseqpair, 21),            // XSeqpair
  XSPECMODELFCT_C_NORM(C_xseqth, 21),              // XSeqtherm
  XSPECMODELFCT_C_NORM(C_equil, 4),                // XSequil
  XSPECMODELFCT_NORM(xsxpdec, 2),                  // XSexpdec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_expcheb6, 11),            // XSexpcheb6
#endif
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
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_vcempow, 20),             // XSvcempow
  XSPECMODELFCT_C_NORM(C_vcheb6, 24),              // XSvcheb6
  XSPECMODELFCT_C_NORM(C_vcie, 17),                // XSvcie
  XSPECMODELFCT_C_NORM(C_vcoolflow, 19),           // XSvcoolflow
#endif
  XSPECMODELFCT_C_NORM(C_vcph, 18),                // XSvcph
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_vexpcheb6, 24),           // XSvexpcheb6
#endif
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
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C_NORM(C_vvcie, 34),               // XSvvcie
  XSPECMODELFCT_C_NORM(C_vvgaussDem, 36),          // XSvvgadem
#endif
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
