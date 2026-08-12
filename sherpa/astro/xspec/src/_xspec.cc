//  Copyright (C) 2007, 2015-2026
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


// Case-insensitive search.
//
static bool is_latest(const std::string &value)
{
  std::string check = value;
  std::transform(check.begin(), check.end(), check.begin(),
		 [](unsigned char c) { return std::tolower(c); });
  return (check == "latest");
}


template <const std::string &getfunc(),
          void setfunc(const std::string&)>
static void set_if_latest(const std::map<std::string, std::string> versionMap,
                          const std::string &key,
                          const std::string &def)
{

  if (!is_latest(getfunc())) {
    return;
  }

  std::map<std::string, std::string>::const_iterator it =
    versionMap.find(key);
  if (it == versionMap.end()) {
    setfunc(def);
  } else {
    setfunc(it->second);
  }
}


// As of XSPEC 12.15.1 the APEC/ATOMDB, NEI, and SPEX versions can be
// set to "latest" and then converted to a value by the XSPEC
// initialization code.
//
// The "latest" string is converted to a version by
//
// - checking for $HEADAS/spectral/modelData/latest.txt file
//
// - if not present, use a hard-coded value (for XSPEC 12.15.1 these
//   values are different depending on if the XSPEC program or
//   model libraries are in use).
//
// The "latest.txt" file is a text key-value file (with no apparent
// error checking) with lines like
//
//     ATOMDB_VERSION: 3.1.3
//     SPEX_VERSION: 3.0.8
//     NEI_VERSION: 3.1.3
//
// Note that the file may not be present because, for instance, the
// "XSPEC data" conda package has not been installed from the HEASOFT
// conda channel.
//
// Note that these "versions" map to
//
//      FunctionUtility::atomdbVersion()
//      FunctionUtility::neiVersion()
//      FunctionUtility::spexVersion()
//
// and to the following "string model" keywords
//
//      APECROOT
//      NEIAPECROOT
//      SPEXROOT
//
// It appears that setting one of these will set the other, at least
// for XSPEC 12.15.1.
//
// There does not appear to be any methods in FunctionUtility that
// provide easy access to handling this (and, even if there were, they
// would not be available to XSPEC libraries prior to 12.15.1), which
// is why there's a lot of "DIY" code here that hopefully matches
// XSPEC.
//
static void validateVersions()
{

  // Are any keywords set to "latest"?
  //
  bool foundLatest = false;
  foundLatest |= is_latest(FunctionUtility::atomdbVersion());
  foundLatest |= is_latest(FunctionUtility::neiVersion());
#ifdef XSPEC_12_15_0
  foundLatest |= is_latest(FunctionUtility::spexVersion());
#endif
  if (!foundLatest) {
    return;
  }

  // Read in the versions from the latest file.
  // There appears to be no validity check, and
  // the lines are expected to be
  //    <key>_VERSION: <version_str>
  //
  const std::string versionPath(FunctionUtility::modelDataPath() +
				"latest.txt");
  std::map<std::string, std::string> versionMap;

  std::ifstream versionFile(versionPath.c_str());
  if (versionFile) {
    const std::string suffix = "_VERSION:";
    std::string lKey, lVersion;
    while (versionFile >> lKey >> lVersion) {
      if (lKey.ends_with(suffix)) {
	const std::string lName = lKey.substr(0, lKey.size() - suffix.size());
	versionMap[lName] = lVersion;
      }
    }
  }

  // It is not obvious what to use as the default value (if the above
  // has no match) since this code supports multiple XSPEC versions.
  // Use the defaults from XSPEC 12.15.1/HEADAS 6.36 which added
  // support for the "latest.txt" file.
  //
  set_if_latest<FunctionUtility::atomdbVersion,
                FunctionUtility::atomdbVersion>
    (versionMap, "ATOMDB", "3.1.3");
  set_if_latest<FunctionUtility::neiVersion,
                FunctionUtility::neiVersion>
    (versionMap, "NEI", "3.1.3");
#ifdef XSPEC_12_15_0
  set_if_latest<FunctionUtility::spexVersion,
                FunctionUtility::spexVersion>
    (versionMap, "SPEX", "3.08");
#endif

}


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

  // Convert "latest" version numbers. Ideally this would only be done
  // if XSPEC 12.15.1 or later were in use, but users can have a XSPEC
  // 12.15.1 ~/.xspec/Xspec.init file - e.g. with lines like
  //
  //      ATOMDB_VERSION:  latest
  //
  // but still be building against XSPEC 12.15.0 (or earlier).
  //
  validateVersions();

  init = true;
  return EXIT_SUCCESS;

}

static PyObject* get_chatter( PyObject *self )
{
  return Py_BuildValue( (char*)"i", FunctionUtility::xwriteChatter() );
}


// TODO:
//   we could send in an integer for the Z number (ie either name
//   or number) but that seems a bit excessive, as the user can
//   get a dict of abundances keyed by the element name.
//
// See also: get_abund_from_table
//
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
  if( !tmpStream.str().empty() ) {
    return PyErr_Format( PyExc_TypeError, // TODO: change from TypeError to ValueError?
			 (char*)"could not find element '%s'", element);
  }

  return (PyObject*) Py_BuildValue( (char*)"f", abundVal );

}


// See also: get_abund
//
// It is simpler to have separate routines rather than to try to deal
// with the multiple options in one routine.
//
static PyObject* get_abund_from_table( PyObject *self, PyObject *args )
{

  // This requires both the table and element name.
  //
  char* table = NULL;
  char* element = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"ss", &table, &element ) )
    return NULL;

  // Get the specific abundance. Unfortunately getAbundance reports an
  // error to stderr when an invalid element is used, so we need to
  // hide this. However it does throw an error if the table is unknown.
  //
  std::ostream* errStream = IosHolder::errHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			&tmpStream);

  float abundVal = 0.0;
  try {
    abundVal = FunctionUtility::getAbundance(string(table),
					     string(element));
  } catch (FunctionUtility::NoInitializer&) {

    IosHolder::setStreams(IosHolder::inHolder(),
			  IosHolder::outHolder(),
			  errStream);

    return PyErr_Format( PyExc_ValueError,
			 "Unknown abundance table '%s'",
			 table );
  }

  IosHolder::setStreams(IosHolder::inHolder(),
			IosHolder::outHolder(),
			errStream);

  // Was there an error?
  //
  if( !tmpStream.str().empty() ) {
    // No backwards compatibility to worry about, so use the sensible
    // error type (ValueError rather than TypeError as used by
    // get_abund).
    //
    return PyErr_Format( PyExc_ValueError,
			 (char*)"could not find element '%s' in table '%s'",
			 element, table );
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
static PyObject* set_abund( PyObject *self, PyObject *args )
{

  char* table = NULL;
  if ( !PyArg_ParseTuple( args, (char*)"s", &table ) )
    return NULL;

  string tableName = string(table);
  tableName = XSutility::lowerCase(tableName);

  if (tableName == "file") {
    // Can not use this if no abundances have been loaded (otherwise
    // XSPEC has been known to crash).
    //
    if (!FunctionUtility::abundChanged()) {
      PyErr_SetString( PyExc_ValueError,
		       (char*)"Abundances have not been read in from a file or array" );
      return NULL;
    }

    FunctionUtility::ABUND(tableName);
    Py_RETURN_NONE;
  }

  if (FunctionUtility::checkAbund(tableName)) {
    FunctionUtility::ABUND(tableName);
    Py_RETURN_NONE;
  }

  // If we've got here then try to read the data from a file. This
  // could be done with a call to FunctionUtility::readNewAbundances()
  // but
  // - it doesn't seem to support reading a file with less then
  //   NELEMS elements,
  // - and if it did it's not clear how to handle the screen output
  //   that (may) be created in that case.
  //
  // So we essentially repeat the readNewAbundaces code here, which
  // has the advantage of not having to throw an error which we then
  // have to catch.
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
  FunctionUtility::abundChanged(true);

  Py_RETURN_NONE;

}


// Handle a vector of abundances. It must be the right size.
// To match set_abund when given a file name we set the
// abundances to "file". This means that a user can not
// load up a set of abundances and *NOT* use them; they
// would have to reset the abundance table after loading.
//
// It looks like we could label these vectors with any value,
// such as "tbl1" or "aneb", rather than "file", which would
// allow multiple tables to be loaded. However, that is for
// later work to see if it is worthwhile (the XSPEC code doesn't
// make it clear how "open" the namespace is here)
//
static PyObject* set_abund_vector( PyObject *self, PyObject *args )
{
  sherpa::astro::xspec::FloatArray vector;
  if ( !PyArg_ParseTuple( args, (char*)"O&",
			  (converter)sherpa::convert_to_contig_array< sherpa::astro::xspec::FloatArray >,
			  &vector ) )
    return NULL;

  size_t nelem = FunctionUtility::NELEMS();
  size_t nvector = static_cast<size_t>(vector.get_size());

  // Rather than worry about what to do with either too many or too
  // few values, just error out.
  //
  if ( nvector != nelem ) {
    return PyErr_Format( PyExc_ValueError,
			 (char*)"Array must contain %d elements, not %d",
			 nelem, nvector );
  }

  std::vector<float> vals(nelem);
  std::copy(&vector[0], &vector[0] + nelem, &vals[0]);

  // Hide the screen output from this call.
  //
  std::ostream* outStream = IosHolder::outHolder();
  std::ostringstream tmpStream;
  IosHolder::setStreams(IosHolder::inHolder(),
			&tmpStream,
			IosHolder::errHolder());

  FunctionUtility::ABUND("file");

  IosHolder::setStreams(IosHolder::inHolder(),
			outStream,
			IosHolder::errHolder());

  FunctionUtility::abundanceVectors("file", vals);
  FunctionUtility::abundChanged(true);

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


static PyObject* clear_xset( PyObject *self )
{
  FunctionUtility::eraseModelStringDataBase();
  Py_RETURN_NONE;
}

static PyObject* set_xset( PyObject *self, PyObject *args )
{

  char* str_name = NULL;
  char* str_value = NULL;

  if ( !PyArg_ParseTuple( args, (char*)"ss", &str_name, &str_value ) )
    return NULL;

  // Sending in INITIALIZE will reset the database but
  // - users can now use the clear_xsxset() routine
  // - using INITIALIZE for this has been marked as deprecated in
  //   4.17.1
  //
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

  if ( !PyArg_ParseTuple( args, (char*)"|s", &str_name ) )
    return NULL;

  // If no argument is given then we return a dictionary
  // of all items.
  //
  if ( str_name == NULL ) {

    PyObject *d = PyDict_New();
    for (const auto& item : FunctionUtility::modelStringDataBase()) {
      PyObject *value = PyUnicode_FromString(item.second.c_str());
      PyDict_SetItemString(d, item.first.c_str(), value);
      Py_DECREF(value);
    }

    return d;
  }

  // Treat an unknown key as an error.
  //
  static string value;
  value = FunctionUtility::getModelString(string(str_name));
  if (value == FunctionUtility::NOT_A_KEY()) {
    PyErr_SetString( PyExc_KeyError, str_name );
    return NULL;
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
  FCTSPEC(get_xsabund_table, get_abund_from_table),
  FCTSPEC(get_xsabund_doc, get_abund_doc),
  FCTSPEC(set_xsabund, set_abund),
  FCTSPEC(set_xsabund_vector, set_abund_vector),
  FCTSPEC(set_xscosmo, set_cosmo),
  NOARGSPEC(get_xscosmo, get_cosmo),
  NOARGSPEC(get_xsxsect, get_xspec_string<FunctionUtility::XSECT>),

  FCTSPEC(set_xsxsect, set_cross),
  NOARGSPEC(clear_xsxset, clear_xset),
  FCTSPEC(set_xsxset, set_xset),
  FCTSPEC(get_xsxset, get_xset),

  // The set commands are not wrapped yet as it's not clear how well
  // the system handles these changes (e.g. it doesn't seem to update
  // the stored abundances if you change one or both of the abundance
  // settings). The cross-section file should also be accessible in a
  // similar manner, but the XSPEC API does not provide access to this
  // (at least for XSPEC 12.12.1).
  //
  // Also, abundPath is essentially managerPath, but we provide access
  // to it as it could be changed (but not by any routine we currently
  // provide access to).
  //
  NOARGSPEC(get_abundance_file, get_xspec_string<FunctionUtility::abundanceFile>),
  NOARGSPEC(get_xspath_abundance, get_xspec_string<FunctionUtility::abundPath>),
  // FCTSPEC(set_abundance_file, set_xspec_string<FunctionUtility::abundanceFile>),
  // FCTSPEC(set_xspath_abundance, set_xspec_string<FunctionUtility::abundPath>),

  NOARGSPEC(get_xspath_manager, get_xspec_string<FunctionUtility::managerPath>),
  NOARGSPEC(get_xspath_model, get_xspec_string<FunctionUtility::modelDataPath>),
  FCTSPEC(set_xspath_manager, set_manager_data_path),

  // Start model definitions
#include "_xspec.hh"
  // End model definitions

  // XSPEC table models
  KWSPECDOC(tabint, sherpa::astro::xspec::xspectablemodel, \
  "The XSPEC table model interface.\n\n" \
  "Parameters\n" \
  "----------\n" \
  "pars : array-like\n" \
  "   The model parameters (size depends on the table model).\n" \
  "xlo, xhi : array-like\n" \
  "   The energy bins (keV) if in ascending order, wavelength (Angstrom)\n" \
  "   if in descending order. ehi can be None, otherwise it must have\n" \
  "   the same size as xlo (which is 1D) and xhi[i] > xlo[i] even when\n" \
  "   the units are in Angstroms.\n" \
  "filename : str\n" \
  "   The name of the XSPEC model file.\n" \
  "tabtype : {'add', 'mul', 'exp'}\n" \
  "   The table type (additive, multiplicative, exponential).\n" \
  "\n" \
  "Returns\n" \
  "-------\n" \
  "modelvals : np.ndarray\n" \
  "   The model values. Matches the size of xlo but if xhi is None then\n" \
  "   the last element is 0.0 (since the xlo values then define the bin\n" \
  "   edges in this case).\n" \
  "\n"
),

  { NULL, NULL, 0, NULL }

};

static struct PyModuleDef xspec_module = {
        PyModuleDef_HEAD_INIT,
        "_xspec",
        PyDoc_STR("Direct access to the XSPEC models and ancillary routines."
                  "\n\n"
                  ".. versionchanged:: 4.19.0\n"
                  "   Models are now accessed using the user name (e.g. 'apec') rather\n"
                  "   than the actual function name (e.g. 'C_apec') from the\n"
                  "   $HEADAS/spectral/manager/model.dat file.\n"
                  "\n"
                  ),
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
