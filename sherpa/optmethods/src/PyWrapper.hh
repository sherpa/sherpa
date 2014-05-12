#ifndef PyWrapper_hh
#define PyWrapper_hh

#include <cstdlib>

namespace sherpa {
  //
  // The file Python.h defines Py_PYTHON_H so if the user is not
  // going to hook up into python, then the following iddef will
  // define the class PyObject so the user can use the class NelderMead
  // without Python.
  //
#ifndef Py_PYTHON_H

  class PyObject {

  };

#define Py_DECREF(arg) (arg=NULL)
#define Py_INCREF(arg) (arg=NULL)

#endif

  class PyWrapper {
  public:

    ~PyWrapper( ) { Py_DECREF( py_function ); }

    PyWrapper( PyObject* ptr ) : py_function( ptr ) { 
      Py_INCREF( py_function );
    }

    PyObject* get_callback_function_ptr( ) { return py_function; }

  private:
    PyObject *py_function;

  };

typedef void (*usrfuncproto)( int n, double* p, double& f, int& e, PyWrapper* x );

}


#endif
