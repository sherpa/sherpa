//
//  Copyright (C) 2007, 2020, 2021  Smithsonian Astrophysical Observatory
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

// This define is needed for "#i" argument to PyArg_ParseTuple in init_optfcn
// and must be made before including Python.h
#define PY_SSIZE_T_CLEAN

#include <Python.h>

#include <sherpa/extension.hh>
#include "tstoptfct.hh"

static PyObject *Ackley( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Ackley<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Ackley Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Booth( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Booth<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Booth Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Bohachevsky1( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Bohachevsky1<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Bohachevsky1 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Bohachevsky2( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Bohachevsky2<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Bohachevsky2 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Bohachevsky3( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Bohachevsky3<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Bohachevsky3 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *BoxBetts( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BoxBetts<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for BoxBetts Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Branin( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Branin<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Branin Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Branin2( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Branin2<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Branin2 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Chichinadze( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Chichinadze<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Chichinadze Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Cola( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Cola<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Cola Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Colville( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Colville<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Colville Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *dcs( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::dcs<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for dcs Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *decanom( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::decanom<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for decanom Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *dodecal( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::dodecal<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for dodecal Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *DixonPrice( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::DixonPrice<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for DixonPrice Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Easom( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Easom<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Easom Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *factor( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::factor<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for factor Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Func1( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Func1<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Func1 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *GoldsteinPrice( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::GoldsteinPrice<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for GoldsteinPrice Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Griewank( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Griewank<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Griewank Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Hansen( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Hansen<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Hansen Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Hartman6( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Hartman6<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Hartman6 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Himmelblau( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Himmelblau<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Himmelblau Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Holzman1( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Holzman1<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Holzman1 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Holzman2( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Holzman2<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Holzman2 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Judge( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Judge<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Judge Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Levy( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Levy<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Levy Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *McCormick( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::McCormick<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for McCormick Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *McKinnon( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::McKinnon<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for McKinnon Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Michalewicz( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Michalewicz<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Michalewicz Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Paviani( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Paviani<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Paviani Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Rastrigin( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Rastrigin<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Rastrigin Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *seqp( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::seqp<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for seqp Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Shekel5( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Shekel5<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Shekel5 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Shekel7( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Shekel7<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Shekel7 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Shekel10( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Shekel10<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Shekel10 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *ShekelModified( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::ShekelModified<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for ShekelModified Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Shubert( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Shubert<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Shubert Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *SixHumpCamel( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::SixHumpCamel<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for SixHumpCamel Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Trecanni( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Trecanni<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Trecanni Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

static PyObject *Trefethen4( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  if ( EXIT_SUCCESS != fvec.create( 1, &npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval = 0;
  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Trefethen4<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Trefethen4 Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );
}

////////////////////////////////////////////////////////////////////////////////
static PyObject *bard( PyObject *self, PyObject *args ) {
  DoubleArray xpar, fvec;
  if ( !PyArg_ParseTuple( args, "O&", CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 15 * npar / 3;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Bard<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				   NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for bard function" );
      return NULL;
    }
  }
  {
    tstoptfct::Bard<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Bard Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *beale( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 3 * npar / 2;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Beale<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				    NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for beale function" );
      return NULL;
    }
  }
  {
    tstoptfct::Beale<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Beale Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *biggs( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 6;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Biggs<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				    NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for biggs function" );
      return NULL;
    }
  }
  {
    tstoptfct::Biggs<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Biggs Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *box3d( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 6;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Box3d<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				    NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for box3d function" );
      return NULL;
    }
  }
  {
    tstoptfct::Box3d<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Box3d Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *broyden_banded( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BroydenBanded<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					    ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for broyden_banded function" );
      return NULL;
    }
  }
  {
    tstoptfct::BroydenBanded<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for BroydenBanded Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *broyden_tridiagonal( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BroydenTridiagonal<double,void*>( mfct, npar, &xpar[0],
						 &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for broyden_tridiagonal function" );
      return NULL;
    }
  }
  {
    tstoptfct::BroydenTridiagonal<double,void*>( npar, &xpar[0], fval, ierr,
						 NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for BroydenTridiagonal Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *brown_almost_linear( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BrownAlmostLinear<double,void*>( mfct, npar, &xpar[0], &fvec[0],
						ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for brown_almost_linear function" );
      return NULL;
    }
  }
  {
    tstoptfct::BrownAlmostLinear<double,void*>( npar, &xpar[0], fval, ierr,
						NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for BrownAlmostLinear Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *brown_badly_scaled( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar + npar / 2;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BrownBadlyScaled<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					       ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for brown_badly_scaled function" );
      return NULL;
    }
  }
  {
    tstoptfct::BrownBadlyScaled<double,void*>( npar, &xpar[0], fval, ierr,
					       NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Brownbadlyscaled Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *brown_dennis( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 20;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::BrownDennis<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
					  NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for brown_dennis function" );
      return NULL;
    }
  }
  {
    tstoptfct::BrownDennis<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for brown_dennis Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *chebyquad( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Chebyquad<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
					NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for chebyquad function" );
      return NULL;
    }
  }
  {
    tstoptfct::Chebyquad<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for chebyquad Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *discrete_boundary( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::DiscreteBoundary<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for discrete_boundary function" );
      return NULL;
    }
  }
  {
    tstoptfct::DiscreteBoundary<double,void*>( npar, &xpar[0], fval, ierr,
					       NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for discrete_boundary_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *discrete_integral( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::DiscreteIntegral<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					       ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for discrete_integral function" );
      return NULL;
    }
  }
  {
    tstoptfct::DiscreteIntegral<double,void*>( npar, &xpar[0], fval, ierr,
					       NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for discrete_integral_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *freudenstein_roth( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::FreudensteinRoth<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for freudenstein_roth function" );
      return NULL;
    }
  }
  {
    tstoptfct::FreudensteinRoth<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for freudenstein_roth_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *gaussian( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 15;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Gaussian<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				       NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for gaussian function" );
      return NULL;
    }
  }
  {
    tstoptfct::Gaussian<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Gaussian Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *gulf_research_development( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::GulfResearchDevelopment<double,void*>( mfct, npar, &xpar[0],
						      &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for gulf_research_development function" );
      return NULL;
    }
  }
  {
    tstoptfct::GulfResearchDevelopment<double,void*>( npar, &xpar[0], fval,
						      ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Gulf_Research_Development Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *helical_valley( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::HelicalValley<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					    ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for helical_valley function" );
      return NULL;
    }
  }
  {
    tstoptfct::HelicalValley<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for helical_valley_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *jennrich_sampson( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 10 * npar / 2;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::JennrichSampson<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					      ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for jennrich_sampson function" );
      return NULL;
    }
  }
  {
    tstoptfct::JennrichSampson<double,void*>( npar, &xpar[0], fval, ierr,
					      NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for jennrich_sampson_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *kowalik_osborne( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 11;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::KowalikOsborne<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					     ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for kowalikosborne function" );
      return NULL;
    }
  }
  {
    tstoptfct::KowalikOsborne<double,void*>( npar, &xpar[0], fval, ierr,
					     NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for kowalikosborne_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *linear_fullrank( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::LinearFullRank<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					     ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank function" );
      return NULL;
    }
  }
  {
    tstoptfct::LinearFullRank<double,void*>( npar, &xpar[0], fval, ierr,
					     NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank1_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *linear_fullrank1( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::LinearFullRank1<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					      ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank1 function" );
      return NULL;
    }
  }
  {
    tstoptfct::LinearFullRank1<double,void*>( npar, &xpar[0], fval, ierr,
					      NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *linear_fullrank0col0rows( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::LinearFullRank0cols0rows<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank0col0rows function" );
      return NULL;
    }
  }
  {
    tstoptfct::LinearFullRank0cols0rows<double,void*>( npar, &xpar[0], fval,
						       ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for linear_fullrank0col0rows_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *meyer( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 16 * npar / 3;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Meyer<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for meyer function" );
      return NULL;
    }
  }
  {
    tstoptfct::Meyer<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for meyer_fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *osborne1( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 33;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Osborne1<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for osborne1 function" );
      return NULL;
    }
  }
  {
    tstoptfct::Osborne1<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for osborne1_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *osborne2( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 65;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Osborne2<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				       NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for osborne2 function" );
      return NULL;
    }
  }
  {
    tstoptfct::Osborne2<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for osborne2_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

 }

 static PyObject *penaltyI( PyObject *self, PyObject *args ) {

   DoubleArray xpar, fvec;

   if ( !PyArg_ParseTuple( args,
			   "O&",
			   CONVERTME(DoubleArray), &xpar ) )
     return NULL;
   npy_intp npar = xpar.get_size( );
   npy_intp mfct = npar + 1;
   if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
     PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
     return NULL;
   }
   double fval;

   int ierr = EXIT_SUCCESS;
   {
     tstoptfct::PenaltyI<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
					NULL );
     if ( EXIT_SUCCESS != ierr ) {
       PyErr_SetString( PyExc_ValueError, "error returned for penaltyI function" );
       return NULL;
     }
   }
   {
     tstoptfct::PenaltyI<double,void*>( npar, &xpar[0], fval, ierr, NULL );
     if ( EXIT_SUCCESS != ierr ) {
       PyErr_SetString( PyExc_ValueError, "error returned for penaltyI_fct function" );
       return NULL;
     }
   }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *penaltyII( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 2 * npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::PenaltyII<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
					NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for penaltyII function" );
      return NULL;
    }
  }
  {
    tstoptfct::PenaltyII<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for penaltyII_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *powell_badly_scaled( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::PowellBadlyScaled<double,void*>( mfct, npar, &xpar[0], &fvec[0],
						ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for powell_badly_scaled function" );
      return NULL;
    }
  }
  {
    tstoptfct::PowellBadlyScaled<double,void*>( npar, &xpar[0], fval, ierr,
						NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for powell_badly_scaled_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *powell_singular( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::PowellSingular<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					     ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for powell_singular function" );
      return NULL;
    }
  }
  {
    tstoptfct::PowellSingular<double,void*>( npar, &xpar[0], fval, ierr,
					     NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for powell_singular_fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *rosenbrock( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Rosenbrock<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
					 NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for rosenbrock function" );
      return NULL;
    }
  }
  {
    tstoptfct::Rosenbrock<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Rosenbrock Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *trigonometric( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Trigonometric<double,void*>( mfct, npar, &xpar[0], &fvec[0],
					    ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for trigonmetric function" );
      return NULL;
    }
  }
  {
    tstoptfct::Trigonometric<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Trigonmetric Fct function" );
      return NULL;
    }
  }
  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *variably_dimensioned( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = npar + 2;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::VariablyDimensioned<double,void*>( mfct, npar, &xpar[0],
						  &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for varibly_dimensioned function" );
      return NULL;
    }
  }
  {
    tstoptfct::VariablyDimensioned<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for variably dimensioned Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *watson( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 31;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Watson<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for watson function" );
      return NULL;
    }
  }
  {
    tstoptfct::Watson<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for watson Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *wood( PyObject *self, PyObject *args ) {

  DoubleArray xpar, fvec;

  if ( !PyArg_ParseTuple( args,
			  "O&",
			  CONVERTME(DoubleArray), &xpar ) )
    return NULL;
  npy_intp npar = xpar.get_size( );
  npy_intp mfct = 6;
  if ( EXIT_SUCCESS != fvec.create( 1, &mfct ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'fvec'" );
    return NULL;
  }
  double fval;

  int ierr = EXIT_SUCCESS;
  {
    tstoptfct::Wood<double,void*>( mfct, npar, &xpar[0], &fvec[0], ierr,
				   NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for wood function" );
      return NULL;
    }
  }
  {
    tstoptfct::Wood<double,void*>( npar, &xpar[0], fval, ierr, NULL );
    if ( EXIT_SUCCESS != ierr ) {
      PyErr_SetString( PyExc_ValueError, "error returned for Wood Fct function" );
      return NULL;
    }
  }

  return Py_BuildValue( "dN", fval, fvec.return_new_ref() );

}

static PyObject *init_optfcn( PyObject *self, PyObject *args ) {

  int npar;
  Py_ssize_t name_length;
  char* name;

  if ( !PyArg_ParseTuple( args,
			  "s#i",
			  &name,
			  &name_length,
			  &npar ) )
    return NULL;

  DoubleArray xpar, lo, hi;
  npy_intp my_npar = npar;
  if ( EXIT_SUCCESS != xpar.create( 1, &my_npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'xpar'" );
    return NULL;
  }
  if ( EXIT_SUCCESS != lo.create( 1, &my_npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'lo'" );
    return NULL;
  }
  if ( EXIT_SUCCESS != hi.create( 1, &my_npar ) ) {
    PyErr_SetString( PyExc_ValueError, "Unable to create 'hi'" );
    return NULL;
  }

  int mfct;
  double answer;
  if ( 0 == strncmp( &name[0], "Ackley", name_length ) )
    tstoptfct::AckleyInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Booth", name_length ) )
    tstoptfct::BoothInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Bohachevsky1", name_length ) )
    tstoptfct::Bohachevsky1Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Bohachevsky2", name_length ) )
    tstoptfct::Bohachevsky2Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Bohachevsky3", name_length ) )
    tstoptfct::Bohachevsky3Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "BoxBetts", name_length ) )
    tstoptfct::BoxBettsInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Branin", name_length ) )
    tstoptfct::BraninInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Branin2", name_length ) )
    tstoptfct::Branin2Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Chichinadze", name_length ) )
    tstoptfct::ChichinadzeInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Cola", name_length ) )
    tstoptfct::ColaInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Colville", name_length ) )
    tstoptfct::ColvilleInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "dcs", name_length ) )
    tstoptfct::dcsInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "decanom", name_length ) )
    tstoptfct::decanomInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "dodecal", name_length ) )
    tstoptfct::dodecalInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "DixonPrice", name_length ) )
    tstoptfct::DixonPriceInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Easom", name_length ) )
    tstoptfct::EasomInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "factor", name_length ) )
    tstoptfct::factorInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Func1", name_length ) )
    tstoptfct::Func1Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "GoldsteinPrice", name_length ) )
    tstoptfct::GoldsteinPriceInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Griewank", name_length ) )
    tstoptfct::GriewankInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Hansen", name_length ) )
    tstoptfct::HansenInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Hartman6", name_length ) )
    tstoptfct::Hartman6Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Himmelblau", name_length ) )
    tstoptfct::HimmelblauInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Holzman1", name_length ) )
    tstoptfct::Holzman1Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Holzman2", name_length ) )
    tstoptfct::Holzman2Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Judge", name_length ) )
    tstoptfct::JudgeInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Levy", name_length ) )
    tstoptfct::LevyInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "McCormick", name_length ) )
    tstoptfct::McCormickInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "McKinnon", name_length ) )
    tstoptfct::McKinnonInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Michalewicz", name_length ) )
    tstoptfct::MichalewiczInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Paviani", name_length ) )
    tstoptfct::PavianiInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Rastrigin", name_length ) )
    tstoptfct::RastriginInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "seqp", name_length ) )
    tstoptfct::seqpInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Shekel5", name_length ) )
    tstoptfct::Shekel5Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Shekel7", name_length ) )
    tstoptfct::Shekel7Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Shekel10", name_length ) )
    tstoptfct::Shekel10Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "ShekelModified", name_length ) )
    tstoptfct::ShekelModifiedInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Shubert", name_length ) )
    tstoptfct::ShubertInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "SixHumpCamel", name_length ) )
    tstoptfct::SixHumpCamelInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Trecanni", name_length ) )
    tstoptfct::TrecanniInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "Trefethen4", name_length ) )
    tstoptfct::Trefethen4Init( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
    //////////////////////////////////////////////////////////////////////
  else if ( 0 == strncmp( &name[0], "rosenbrock", name_length ) )
    tstoptfct::RosenbrockInit( npar, mfct, answer, &xpar[0],
			       &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "freudenstein_roth", name_length ) )
    tstoptfct::FreudensteinRothInit( npar, mfct, answer, &xpar[0],
				     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "powell_badly_scaled", name_length ) )
    tstoptfct::PowellBadlyScaledInit( npar, mfct, answer, &xpar[0],
				      &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "brown_badly_scaled", name_length ) )
    tstoptfct::BrownBadlyScaledInit( npar, mfct, answer, &xpar[0],
				     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "beale", name_length ) )
    tstoptfct::BealeInit( npar, mfct, answer, &xpar[0],
				     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "jennrich_sampson", name_length ) )
    tstoptfct::JennrichSampsonInit( npar, mfct, answer, &xpar[0],
				    &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "helical_valley", name_length ) )
    tstoptfct::HelicalValleyInit( npar, mfct, answer, &xpar[0],
				    &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "bard", name_length ) )
    tstoptfct::BardInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "gaussian", name_length ) )
    tstoptfct::GaussianInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "meyer", name_length ) )
    tstoptfct::MeyerInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "gulf_research_development",
			  name_length ) )
    tstoptfct::GulfResearchDevelopmentInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "box3d", name_length ) )
    tstoptfct::Box3dInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "powell_singular", name_length ) )
    tstoptfct::PowellSingularInit( npar, mfct, answer, &xpar[0],
				   &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "wood", name_length ) )
    tstoptfct::WoodInit( npar, mfct, answer, &xpar[0],
			 &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "kowalik_osborne", name_length ) )
    tstoptfct::KowalikOsborneInit( npar, mfct, answer, &xpar[0],
				   &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "brown_dennis", name_length ) )
    tstoptfct::BrownDennisInit( npar, mfct, answer, &xpar[0],
				&lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "osborne1", name_length ) )
    tstoptfct::Osborne1Init( npar, mfct, answer, &xpar[0],
			     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "biggs", name_length ) )
    tstoptfct::BiggsInit( npar, mfct, answer, &xpar[0],
			  &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "osborne2", name_length ) )
    tstoptfct::Osborne2Init( npar, mfct, answer, &xpar[0],
			     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "watson", name_length ) )
    tstoptfct::WatsonInit( npar, mfct, answer, &xpar[0],
			   &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "penaltyI", name_length ) )
    tstoptfct::PenaltyIInit( npar, mfct, answer, &xpar[0],
			     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "penaltyII", name_length ) )
    tstoptfct::PenaltyIIInit( npar, mfct, answer, &xpar[0],
			      &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "variably_dimensioned", name_length ) )
    tstoptfct::VariablyDimensionedInit( npar, mfct, answer, &xpar[0],
					&lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "trigonometric", name_length ) )
    tstoptfct::TrigonometricInit( npar, mfct, answer, &xpar[0],
				 &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "brown_almost_linear", name_length ) )
    tstoptfct::BrownAlmostLinearInit( npar, mfct, answer, &xpar[0],
				      &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "discrete_boundary", name_length ) )
    tstoptfct::DiscreteBoundaryInit( npar, mfct, answer, &xpar[0],
				     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "discrete_integral", name_length ) )
    tstoptfct::DiscreteIntegralInit( npar, mfct, answer, &xpar[0],
				     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "broyden_tridiagonal", name_length ) )
    tstoptfct::BroydenTridiagonalInit( npar, mfct, answer, &xpar[0],
				       &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "broyden_banded", name_length ) )
    tstoptfct::BroydenBandedInit( npar, mfct, answer, &xpar[0],
				       &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "linear_fullrank", name_length ) )
    tstoptfct::LinearFullRankInit( npar, mfct, answer, &xpar[0],
				  &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "linear_fullrank1", name_length ) )
    tstoptfct::LinearFullRank1Init( npar, mfct, answer, &xpar[0],
				  &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "linear_fullrank0cols0rows", name_length ) )
    tstoptfct::LinearFullRank0cols0rowsInit( npar, mfct, answer, &xpar[0],
					     &lo[0], &hi[0] );
  else if ( 0 == strncmp( &name[0], "chebyquad", name_length ) )
    tstoptfct::ChebyquadInit( npar, mfct, answer, &xpar[0], &lo[0], &hi[0] );
  else
    fprintf( stderr, "init_optfcn: Unrecognized test function '%s'\n",
	     &name[0] );

  if ( xpar.get_size() != lo.get_size() ||
       xpar.get_size() != hi.get_size() ) {
    PyErr_Format( PyExc_ValueError,
		  "init_optfcn: Incompatible array sizes "
	          "xpar=%d, lo=%d, hi=%d\n",
	          (int) xpar.get_size(), (int) lo.get_size(),
                  (int) hi.get_size() );
    return NULL;
  }

  return Py_BuildValue( "NNNd", xpar.return_new_ref(),
			lo.return_new_ref(), hi.return_new_ref(), answer );

}

// A listing of our methods/functions:
static PyMethodDef _tstoptfct_methods[] = {
  // name, function, argument type, docstring

  { "init", init_optfcn, METH_VARARGS, "init starting params and bounds" },
  { "Ackley", Ackley, METH_VARARGS, "ackley function vector" },
  { "Booth", Booth, METH_VARARGS, "booth function vector" },
  { "Bohachevsky1", Bohachevsky1, METH_VARARGS, "bohachevsky1 function vector" },
  { "Bohachevsky2", Bohachevsky2, METH_VARARGS, "bohachevsky2 function vector" },
  { "Bohachevsky3", Bohachevsky3, METH_VARARGS, "bohachevsky3 function vector" },
  { "BoxBetts", BoxBetts, METH_VARARGS, "boot function vector" },
  { "Branin", Branin, METH_VARARGS, "branin function vector" },
  { "Branin2", Branin2, METH_VARARGS, "branin2 function vector" },
  { "Chichinadze", Chichinadze, METH_VARARGS, "boot function vector" },
  { "Cola", Cola, METH_VARARGS, "cola function vector" },
  { "Colville", Colville, METH_VARARGS, "colville function vector" },
  { "dcs", dcs, METH_VARARGS, "dcs function vector" },
  { "decanom", decanom, METH_VARARGS, "decanom function vector" },
  { "dodecal", dodecal, METH_VARARGS, "dodecal function vector" },
  { "DixonPrice", DixonPrice, METH_VARARGS, "DixonPrice function vector" },
  { "Easom", Easom, METH_VARARGS, "Easom function vector" },
  { "factor", factor, METH_VARARGS, "factor function vector" },
  { "Func1", Func1, METH_VARARGS, "Func1 function vector" },
  { "GoldsteinPrice", GoldsteinPrice, METH_VARARGS, "GoldsteinPrice function vector" },
  { "Griewank", Griewank, METH_VARARGS, "Griewank function vector" },
  { "Hansen", Hansen, METH_VARARGS, "Hansen function vector" },
  { "Hartman6", Hartman6, METH_VARARGS, "Hartman6 function vector" },
  { "Himmelblau", Himmelblau, METH_VARARGS, "Himmelblau function vector" },
  { "Holzman1", Holzman1, METH_VARARGS, "Holzman1 function vector" },
  { "Holzman2", Holzman2, METH_VARARGS, "Holzman2 function vector" },
  { "Judge", Judge, METH_VARARGS, "Judge function vector" },
  { "Levy", Levy, METH_VARARGS, "Levy function vector" },
  { "McCormick", McCormick, METH_VARARGS, "McCormick function vector" },
  { "McKinnon", McKinnon, METH_VARARGS, "McKinnon function vector" },
  { "Michalewicz", Michalewicz, METH_VARARGS, "Michalewicz function vector" },
  { "Paviani", Paviani, METH_VARARGS, "Paviani function vector" },
  { "Rastrigin", Rastrigin, METH_VARARGS, "Rastrigin function vector" },
  { "seqp", seqp, METH_VARARGS, "seqp function vector" },
  { "Shekel5", Shekel5, METH_VARARGS, "Shekel5 function vector" },
  { "Shekel7", Shekel7, METH_VARARGS, "Shekel7 function vector" },
  { "Shekel10", Shekel10, METH_VARARGS, "Shekel10 function vector" },
  { "ShekelModified", ShekelModified, METH_VARARGS, "ShekelModified function vector" },
  { "Shubert", Shubert, METH_VARARGS, "Shubert function vector" },
  { "SixHumpCamel", SixHumpCamel, METH_VARARGS, "SixHumpCamel function vector" },
  { "Trecanni", Trecanni, METH_VARARGS, "Trecanni function vector" },
    { "Trefethen4", Trefethen4, METH_VARARGS, "Trefethen4 function vector" },
  ////////////////////////////////////////////////////////////////////////
  { "bard", bard, METH_VARARGS, "bard function vector" },
  { "beale", beale, METH_VARARGS, "beale function vector" },
  { "biggs", biggs, METH_VARARGS, "biggs function vector" },
  { "box3d", box3d, METH_VARARGS, "box3d function vector" },
  { "broyden_banded", broyden_banded, METH_VARARGS, "broyden_banded function vector" },
  { "broyden_tridiagonal", broyden_tridiagonal, METH_VARARGS, "broyden_tridiagonal function vector" },
  { "brown_almost_linear", brown_almost_linear, METH_VARARGS, "brown_almost_linear function vector" },
  { "brown_badly_scaled", brown_badly_scaled, METH_VARARGS, "brown_badly_scaled function vector" },

  { "brown_dennis", brown_dennis, METH_VARARGS, "brown_dennis function vector" },
  { "chebyquad", chebyquad, METH_VARARGS, "chebyquad function vector" },

  { "discrete_boundary", discrete_boundary, METH_VARARGS, "discrete_boundary function vector" },

  { "discrete_integral", discrete_integral, METH_VARARGS, "discrete_integral function vector" },

  { "freudenstein_roth", freudenstein_roth, METH_VARARGS, "freudenstein_roth function vector" },

  { "gaussian", gaussian, METH_VARARGS, "gaussian function vector" },

  { "gulf_research_development", gulf_research_development, METH_VARARGS, "gulf_research_development function vector" },

  { "helical_valley", helical_valley, METH_VARARGS, "helical_valley function vector" },

  { "jennrich_sampson", jennrich_sampson, METH_VARARGS, "jennrich_sampson function vector" },

  { "kowalik_osborne", kowalik_osborne, METH_VARARGS, "kowalik_osborne function vector" },

  { "linear_fullrank", linear_fullrank, METH_VARARGS, "linear_fullrank function vector" },

  { "linear_fullrank1", linear_fullrank1, METH_VARARGS, "linear_fullrank1 function vector" },

  { "linear_fullrank0cols0rows", linear_fullrank0col0rows, METH_VARARGS, "linear_fullrank0cols0rows function vector" },

  { "meyer", meyer, METH_VARARGS, "meyer function vector" },

  { "osborne1", osborne1, METH_VARARGS, "osborne1 function vector" },

  { "osborne2", osborne2, METH_VARARGS, "osborne2 function vector" },

  { "penaltyI", penaltyI, METH_VARARGS, "penaltyI function vector" },

  { "penaltyII", penaltyII, METH_VARARGS, "penaltyII function vector" },

  { "powell_badly_scaled", powell_badly_scaled, METH_VARARGS, "powell_badly_scaled function vector" },

  { "powell_singular", powell_singular, METH_VARARGS, "powellsingular function vector" },

  { "rosenbrock", rosenbrock, METH_VARARGS, "rosenbrock function vector" },

  { "trigonometric", trigonometric, METH_VARARGS, "trigonometric function vector" },

  { "variably_dimensioned", variably_dimensioned, METH_VARARGS, "Variably_Dimensioned function vector" },


  { "watson", watson, METH_VARARGS, "Watson function vector" },


  { "wood", wood, METH_VARARGS, "Wood function vector" },


  { NULL, NULL, 0, NULL }                 /* sentinel */
};

SHERPAMOD( _tstoptfct, _tstoptfct_methods );
