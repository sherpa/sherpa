// 
//  Copyright (C) 2008,2009  Smithsonian Astrophysical Observatory
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
