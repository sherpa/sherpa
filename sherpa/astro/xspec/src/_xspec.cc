//  Copyright (C) 2007, 2015 - 2025
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

  XSPECMODELFCT_C(C_agauss, 2),                    // XSagauss
  XSPECMODELFCT(agnsed, 15),                       // XSagnsed
  XSPECMODELFCT(agnslim, 14),                      // XSagnslim
  XSPECMODELFCT_C(C_apec, 3),                      // XSapec
  XSPECMODELFCT_C(C_bapec, 4),                     // XSbapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bcempow, 7),                   // XSbcempow
  XSPECMODELFCT_C(C_bcheb6, 11),                   // XSbcheb6
  XSPECMODELFCT_C(C_bcie, 5),                      // XSbcie
  XSPECMODELFCT_C(C_bcoolflow, 6),                 // XSbcoolflow
  XSPECMODELFCT_C(C_bcph, 5),                      // XSbcph
  XSPECMODELFCT_C(C_bequil, 4),                    // XSbequil
  XSPECMODELFCT_C(C_bexpcheb6, 11),                // XSbexpcheb6
  XSPECMODELFCT_C(C_bgaussDem, 7),                 // XSbgadem
  XSPECMODELFCT_C(C_bgnei, 6),                     // XSbgnei
  XSPECMODELFCT_C(C_bnei, 5),                      // XSbnei
  XSPECMODELFCT_C(C_bsnapec, 7),                   // XSbsnapec
#endif
  XSPECMODELFCT_C(C_btapec, 5),                    // XSbtapec
  XSPECMODELFCT(xsblbd, 1),                        // XSbbody
  XSPECMODELFCT(xsbbrd, 1),                        // XSbbodyrad
  XSPECMODELFCT_C(C_xsbexrav, 9),                  // XSbexrav
  XSPECMODELFCT_C(C_xsbexriv, 11),                 // XSbexriv
  XSPECMODELFCT_C(C_brokenPowerLaw, 3),            // XSbknpower
  XSPECMODELFCT_C(C_broken2PowerLaw, 5),           // XSbkn2pow
  XSPECMODELFCT(xsbmc, 3),                         // XSbmc
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bnpshock, 7),                  // XSbnpshock
  XSPECMODELFCT_C(C_bpshock, 6),                   // XSbpshock
#endif
  XSPECMODELFCT(xsbrms, 1),                        // XSbremss
  XSPECMODELFCT_C(C_brnei, 6),                     // XSbrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bsedov, 6),                    // XSbsedov
#endif
  XSPECMODELFCT_C(C_bvapec, 16),                   // XSbvapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bvcempow, 20),                 // XSbvcempow
  XSPECMODELFCT_C(C_bvcheb6, 24),                  // XSbvcheb6
  XSPECMODELFCT_C(C_bvcie, 17),                    // XSbvcie
  XSPECMODELFCT_C(C_bvcoolflow, 19),               // XSbvcoolflow
  XSPECMODELFCT_C(C_bvcph, 18),                    // XSbvcph
  XSPECMODELFCT_C(C_bvequil, 15),                  // XSbvequil
  XSPECMODELFCT_C(C_bvexpcheb6, 24),               // XSbvexpcheb6
  XSPECMODELFCT_C(C_bvgaussDem, 20),               // XSbvgadem
  XSPECMODELFCT_C(C_bvgnei, 18),                   // XSbvgnei
  XSPECMODELFCT_C(C_bvnei, 17),                    // XSbvnei
  XSPECMODELFCT_C(C_bvnpshock, 19),                // XSbvnpshock
  XSPECMODELFCT_C(C_bvpshock, 18),                 // XSbvpshock
#endif
  XSPECMODELFCT_C(C_bvrnei, 18),                   // XSbvrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bvsedov, 18),                  // XSbvsedov
#endif
  XSPECMODELFCT_C(C_bvtapec, 17),                  // XSbvtapec
  XSPECMODELFCT_C(C_bvvapec, 33),                  // XSbvvapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bvvcie, 34),                   // XSbvvcie
  XSPECMODELFCT_C(C_bvvgaussDem, 36),              // XSbvvgadem
  XSPECMODELFCT_C(C_bvvgnei, 35),                  // XSbvvgnei
  XSPECMODELFCT_C(C_bvvnei, 34),                   // XSbvvnei
  XSPECMODELFCT_C(C_bvvnpshock, 36),               // XSbvvnpshock
  XSPECMODELFCT_C(C_bvvpshock, 35),                // XSbvvpshock
#endif
  XSPECMODELFCT_C(C_bvvrnei, 35),                  // XSbvvrnei
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bvvsedov, 35),                 // XSbvvsedov
#endif
  XSPECMODELFCT_C(C_bvvtapec, 34),                 // XSbvvtapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_bvvwDem, 37),                  // XSbvvwdem
  XSPECMODELFCT_C(C_bvwDem, 21),                   // XSbvwdem
  XSPECMODELFCT_C(C_bwDem, 8),                     // XSbwdem
#endif
  XSPECMODELFCT_C(C_c6mekl, 10),                   // XSc6mekl
  XSPECMODELFCT_C(C_c6pmekl, 10),                  // XSc6pmekl
  XSPECMODELFCT_C(C_c6pvmkl, 23),                  // XSc6pvmkl
  XSPECMODELFCT_C(C_c6vmekl, 23),                  // XSc6vmekl
  XSPECMODELFCT_C(C_carbatm, 3),                   // XScarbatm
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_cemMekal, 6),                  // XScemekl
#else
  XSPECMODELFCT(cemekl, 6),                        // XScemekl
#endif
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_cempow, 6),                    // XScempow
#endif
  XSPECMODELFCT_C(C_cemVMekal, 19),                // XScevmkl
  XSPECMODELFCT_C(C_xscflw, 5),                    // XScflow
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_cheb6, 10),                    // XScheb6
  XSPECMODELFCT_C(C_cie, 4),                       // XScie
#endif
  XSPECMODELFCT(compbb, 3),                        // XScompbb
  XSPECMODELFCT_C(xscompmag, 8),                   // XScompmag
  XSPECMODELFCT(compls, 2),                        // XScompLS
  XSPECMODELFCT_C(C_xscompps, 19),                 // XScompPS
  XSPECMODELFCT(compst, 2),                        // XScompST
  XSPECMODELFCT_C(xscomptb, 6),                    // XScomptb
  XSPECMODELFCT_C(C_xscompth, 20),                 // XScompth
  XSPECMODELFCT(xstitg, 5),                        // XScompTT
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_coolflow, 5),                  // XScoolflow
#endif
  XSPECMODELFCT_C(C_cph, 4),                       // XScph
  XSPECMODELFCT_C(C_cplinear, 20),                 // XScplinear
  XSPECMODELFCT_C(C_cutoffPowerLaw, 2),            // XScutoffpl
  XSPECMODELFCT(disk, 3),                          // XSdisk
  XSPECMODELFCT(diskir, 8),                        // XSdiskir
  XSPECMODELFCT(xsdskb, 1),                        // XSdiskbb
  XSPECMODELFCT_C(C_diskline, 5),                  // XSdiskline
  XSPECMODELFCT(diskm, 4),                         // XSdiskm
  XSPECMODELFCT(disko, 4),                         // XSdisko
  XSPECMODELFCT(diskpbb, 2),                       // XSdiskpbb
  XSPECMODELFCT(xsdiskpn, 2),                      // XSdiskpn
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_eebremss, 3),                  // XSeebremss
#endif
  XSPECMODELFCT(eplogpar, 2),                      // XSeplogpar
  XSPECMODELFCT_C(C_xseqpair, 20),                 // XSeqpair
  XSPECMODELFCT_C(C_xseqth, 20),                   // XSeqtherm
  XSPECMODELFCT_C(C_equil, 3),                     // XSequil
  XSPECMODELFCT(xsxpdec, 1),                       // XSexpdec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_expcheb6, 10),                 // XSexpcheb6
#endif
  XSPECMODELFCT(ezdiskbb, 1),                      // XSezdiskbb
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_FeKfromSevenLorentzians, 0),   // XSfeklor
#endif
  XSPECMODELFCT_C(C_gaussianLine, 2),              // XSgaussian
  XSPECMODELFCT_C(C_gaussDem, 6),                  // XSgadem
  XSPECMODELFCT_C(C_gnei, 5),                      // XSgnei
  XSPECMODELFCT(grad, 6),                          // XSgrad
  XSPECMODELFCT_C(xsgrbcomp, 9),                   // XSgrbcomp
  XSPECMODELFCT_C(xsgrbjet, 13),                   // XSgrbjet
  XSPECMODELFCT(xsgrbm, 3),                        // XSgrbm
  XSPECMODELFCT_C(C_hatm, 3),                      // XShatm
  XSPECMODELFCT(jet, 15),                          // XSjet
  XSPECMODELFCT_C(C_kerrbb, 9),                    // XSkerrbb
  XSPECMODELFCT_C(C_kerrd, 7),                     // XSkerrd
  XSPECMODELFCT_C(C_spin, 9),                      // XSkerrdisk

  XSPECMODELFCT_CON_F77(kyconv, 12),               // XSkyconv

  XSPECMODELFCT(kyrline, 11),                      // XSkyrline
  XSPECMODELFCT_C(C_laor, 5),                      // XSlaor
  XSPECMODELFCT_C(C_laor2, 7),                     // XSlaor2
  XSPECMODELFCT_C(C_logpar, 3),                    // XSlogpar
  XSPECMODELFCT_C(C_lorentzianLine, 2),            // XSlorentz
  XSPECMODELFCT_C(C_meka, 4),                      // XSmeka
  XSPECMODELFCT_C(C_mekal, 5),                     // XSmekal
  XSPECMODELFCT_C(C_xsmkcf, 5),                    // XSmkcflow
  XSPECMODELFCT_C(C_nei, 4),                       // XSnei
  XSPECMODELFCT_C(C_nlapec, 3),                    // XSnlapec
  XSPECMODELFCT_C(C_npshock, 6),                   // XSnpshock
  XSPECMODELFCT(nsa, 4),                           // XSnsa
  XSPECMODELFCT(nsagrav, 3),                       // XSnsagrav
  XSPECMODELFCT(nsatmos, 4),                       // XSnsatmos
  XSPECMODELFCT_C(C_nsmax, 3),                     // XSnsmax
  XSPECMODELFCT_C(C_nsmaxg, 5),                    // XSnsmaxg
  XSPECMODELFCT_C(C_nsx, 5),                       // XSnsx
  XSPECMODELFCT_C(C_xsnteea, 15),                  // XSnteea
  XSPECMODELFCT_C(C_nthcomp, 5),                   // XSnthComp
  XSPECMODELFCT(optxagn, 13),                      // XSoptxagn
  XSPECMODELFCT(optxagnf, 11),                     // XSoptxagnf
  XSPECMODELFCT(xspegp, 3),                        // XSpegpwrlw
  XSPECMODELFCT(pexmon, 7),                        // XSpexmon
  XSPECMODELFCT_C(C_xspexrav, 7),                  // XSpexrav
  XSPECMODELFCT_C(C_xspexriv, 9),                  // XSpexriv
  XSPECMODELFCT(xsp1tr, 10),                       // XSplcabs
  XSPECMODELFCT_C(C_powerLaw, 1),                  // XSpowerlaw
  XSPECMODELFCT(xsposm, 0),                        // XSposm
  XSPECMODELFCT_C(C_pshock, 5),                    // XSpshock
  XSPECMODELFCT(qsosed, 6),                        // XSqsosed
  XSPECMODELFCT_C(C_raysmith, 3),                  // XSraymond
  XSPECMODELFCT(xredge, 2),                        // XSredge
  XSPECMODELFCT(xsrefsch, 13),                     // XSrefsch

#ifdef XSPEC_12_15_0
  XSPECMODELFCT_CON(C_rgsExtendedSource, 2),       // XSrgsext
#endif
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_CON(C_rgsxsrc, 1),                 // XSrgsxsrc
#else
  XSPECMODELFCT_CON_F77(rgsxsrc, 1),               // XSrgsxsrc
#endif

  XSPECMODELFCT_C(C_rnei, 5),                      // XSrnei
  XSPECMODELFCT_C(C_sedov, 5),                     // XSsedov
  XSPECMODELFCT_C(C_sirf, 9),                      // XSsirf
  XSPECMODELFCT_C(slimbbmodel, 9),                 // XSslimbh
  XSPECMODELFCT_C(C_snapec, 6),                    // XSsnapec
  XSPECMODELFCT(srcut, 2),                         // XSsrcut
  XSPECMODELFCT(sresc, 2),                         // XSsresc
  XSPECMODELFCT(ssa, 2),                           // XSssa
#ifdef XSPEC_12_14_1
  XSPECMODELFCT(sssed, 14),                        // XSsssed
#endif
  XSPECMODELFCT(xsstep, 2),                        // XSstep
  XSPECMODELFCT_C(C_tapec, 4),                     // XStapec
#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_vagauss, 2),                   // XSvagauss
#endif
  XSPECMODELFCT_C(C_vapec, 15),                    // XSvapec
  XSPECMODELFCT(xsbrmv, 2),                        // XSvbremss
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_vcempow, 19),                  // XSvcempow
  XSPECMODELFCT_C(C_vcheb6, 23),                   // XSvcheb6
  XSPECMODELFCT_C(C_vcie, 16),                     // XSvcie
  XSPECMODELFCT_C(C_vcoolflow, 18),                // XSvcoolflow
#endif
  XSPECMODELFCT_C(C_vcph, 17),                     // XSvcph
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_vexpcheb6, 23),                // XSvexpcheb6
#endif
  XSPECMODELFCT_C(C_vequil, 14),                   // XSvequil
  XSPECMODELFCT_C(C_vgaussDem, 19),                // XSvgadem
#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_vgaussianLine, 2),             // XSvgaussian
#endif
  XSPECMODELFCT_C(C_vgnei, 17),                    // XSvgnei
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_vlorentzianLine, 2),           // XSvlorentz
#endif
  XSPECMODELFCT_C(C_vmeka, 17),                    // XSvmeka
  XSPECMODELFCT_C(C_vmekal, 18),                   // XSvmekal
  XSPECMODELFCT_C(C_xsvmcf, 18),                   // XSvmcflow
  XSPECMODELFCT_C(C_vnei, 16),                     // XSvnei
  XSPECMODELFCT_C(C_vnpshock, 18),                 // XSvnpshock
  XSPECMODELFCT_C(C_voigtLine, 3),                 // XSvoigt
  XSPECMODELFCT_C(C_vpshock, 17),                  // XSvpshock
  XSPECMODELFCT_C(C_vraysmith, 14),                // XSvraymond
  XSPECMODELFCT_C(C_vrnei, 17),                    // XSvrnei
  XSPECMODELFCT_C(C_vsedov, 17),                   // XSvsedov
  XSPECMODELFCT_C(C_vtapec, 16),                   // XSvtapec
  XSPECMODELFCT_C(C_vvapec, 32),                   // XSvvapec
#ifdef XSPEC_12_14_0
  XSPECMODELFCT_C(C_vvcie, 33),                    // XSvvcie
  XSPECMODELFCT_C(C_vvgaussDem, 35),               // XSvvgadem
#endif
  XSPECMODELFCT_C(C_vvgnei, 34),                   // XSvvgnei
  XSPECMODELFCT_C(C_vvnei, 33),                    // XSvvnei
  XSPECMODELFCT_C(C_vvnpshock, 35),                // XSvvnpshock
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_vvoigtLine, 3),                // XSvvoigt
#endif
  XSPECMODELFCT_C(C_vvpshock, 34),                 // XSvvpshock
  XSPECMODELFCT_C(C_vvrnei, 34),                   // XSvvrnei
  XSPECMODELFCT_C(C_vvsedov, 34),                  // XSvvsedov
  XSPECMODELFCT_C(C_vvtapec, 33),                  // XSvvtapec
  XSPECMODELFCT_C(C_vvwDem, 36),                   // XSvvwdem
  XSPECMODELFCT_C(C_vwDem, 20),                    // XSvwdem
  XSPECMODELFCT_C(C_wDem, 7),                      // XSwdem
  XSPECMODELFCT_C(C_zagauss, 3),                   // XSzagauss
  XSPECMODELFCT(xszbod, 2),                        // XSzbbody
  XSPECMODELFCT_C(C_zBrokenPowerLaw, 4),           // XSzbknpower
  XSPECMODELFCT(xszbrm, 2),                        // XSzbremss
  XSPECMODELFCT_C(C_zcutoffPowerLaw, 3),           // XSzcutoffpl
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_zFeKfromSevenLorentzians, 1),  // XSzfeklor
#endif
  XSPECMODELFCT_C(C_xszgau, 3),                    // XSzgauss
  XSPECMODELFCT_C(C_zkerrbb, 9),                   // XSzkerrbb
  XSPECMODELFCT_C(C_zLogpar, 4),                   // XSzlogpar
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_zlorentzianLine, 3),           // XSzlorentz
#endif
  XSPECMODELFCT_C(C_zpowerLaw, 2),                 // XSzpowerlw
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_zvlorentzianLine, 3),          // XSzvlorentz
  XSPECMODELFCT_C(C_zvoigtLine, 4),                // XSzvoigt
  XSPECMODELFCT_C(C_zvvoigtLine, 4),               // XSzvvoigt
#endif

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
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_lorentzianAbsorptionLine, 3),  // XSlorabs
#endif
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
#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_vgaussianAbsorptionLine, 3),   // XSvgabs
#endif
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_vlorentzianAbsorptionLine, 3), // XSvlorabs
  XSPECMODELFCT_C(C_voigtAbsorptionLine, 4),       // XSvoigtabs
#endif
  XSPECMODELFCT(xsvphb, 18),                       // XSvphabs
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_vvoigtAbsorptionLine, 4),      // XSvvoigtabs
#endif
  XSPECMODELFCT(xsabsw, 1),                        // XSwabs
  XSPECMODELFCT(xswnab, 2),                        // XSwndabs
  XSPECMODELFCT(xsxirf, 13),                       // XSxion
  XSPECMODELFCT_C(C_xscatmodel, 4),                // XSxscat
  XSPECMODELFCT_C(xszbabs, 4),                     // XSzbabs
  XSPECMODELFCT(mszdst, 4),                        // XSzdust
  XSPECMODELFCT(xszedg, 3),                        // XSzedge
#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_zgaussianAbsorptionLine, 4),   // XSzgabs
#endif
  XSPECMODELFCT(xszhcu, 3),                        // XSzhighect
  XSPECMODELFCT(zigm, 3),                          // XSzigm
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_zlorentzianAbsorptionLine, 4), // XSzlorabs
#endif
  XSPECMODELFCT(xszabp, 3),                        // XSzpcfabs
  XSPECMODELFCT(xszphb, 2),                        // XSzphabs
#ifdef XSPEC_12_15_0
  XSPECMODELFCT_C(C_zvlorentzianAbsorptionLine, 4), // XSzvlorabs
  XSPECMODELFCT_C(C_zvoigtAbsorptionLine, 5),      // XSzvoigtabs
  XSPECMODELFCT_C(C_zvvoigtAbsorptionLine, 5),     // XSzvvoigtabs
#endif
  XSPECMODELFCT(zxipab, 5),                        // XSzxipab
  XSPECMODELFCT_C(C_zxipcf, 4),                    // XSzxipcf
  XSPECMODELFCT(xszcrd, 2),                        // XSzredden
  XSPECMODELFCT(msldst, 4),                        // XSzsmdust
  XSPECMODELFCT_C(C_ztbabs, 2),                    // XSzTBabs

#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_zvagauss, 3),                  // XSzvagauss
#endif

  XSPECMODELFCT(xszvab, 19),                       // XSzvarabs
  XSPECMODELFCT(xszvfe, 5),                        // XSzvfeabs
#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_zvgaussianAbsorptionLine, 4),  // XSzvgabs
#endif

#ifdef XSPEC_12_14_1
  XSPECMODELFCT_C(C_zvgaussianLine, 3),            // XSzvgaussian
#endif

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
  XSPECMODELFCT_CON(C_simpl, 3),                   // XSsimpl
  XSPECMODELFCT_CON_F77(thcompf, 4),               // XSthcomp
  XSPECMODELFCT_CON(C_vashift, 1),                 // XSvashift
  XSPECMODELFCT_CON(C_vmshift, 1),                 // XSvmshift
  XSPECMODELFCT_CON(C_xilconv, 6),                 // XSxilconv
  XSPECMODELFCT_CON(C_zashift, 1),                 // XSzashift
  XSPECMODELFCT_CON(C_zmshift, 1),                 // XSzmshift

  XSPECMODELFCT_C(beckerwolff, 12),                // XSbwcycl

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
