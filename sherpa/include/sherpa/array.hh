// 
//  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_array_hh__
#define __sherpa_array_hh__

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>


namespace sherpa {


  template <typename CType, int ArrayType>
  class Array {

  public:

    ~Array() { Py_XDECREF( array ); }

    Array()
      : array( NULL ),
	data( NULL ),
	stride( 0 ),
	size( 0 )
    { }

    // Allows using an instance as a truth value
    operator void*() const { return array; }

    int create( int ndim, const npy_intp* dims, CType* cdata=NULL )
    {
      return init( PyArray_New( &PyArray_Type, ndim,
				const_cast< npy_intp* >( dims ), ArrayType,
				NULL, cdata, 0, NPY_ARRAY_CARRAY, NULL ) );
    }

    int from_obj( PyObject* obj, bool contiguous=false );

    int zeros( int ndim, const npy_intp* dims )
    {
      return init( PyArray_Zeros( ndim, const_cast< npy_intp* >( dims ),
				  PyArray_DescrFromType( ArrayType ), 0 ) ); 
    }

    PyObject* borrowed_ref()
    {
      return array;
    }

    PyObject* new_ref()
    {
      Py_XINCREF( array );
      return array;
    }

    PyObject* return_new_ref()
    {
      return PyArray_Return( static_cast< PyArrayObject* >
			     ( static_cast< void* >( new_ref() ) ) );
    }

    npy_intp get_size() const
    {
      return size;
    }

    int get_ndim() const
    {
      return PyArray_NDIM( (PyArrayObject*) array );
    }

    const npy_intp* get_dims() const
    {
      return PyArray_DIMS( (PyArrayObject*) array );
    }

    // What should the ref field be set to: it seems like 1 will check that
    // there are no references to the array, which seems to make sense
    // here.
    //
    void resize1d( npy_intp newsize ) {
      npy_intp dims[1] = { newsize };
      PyArray_Dims newdims;
      newdims.ptr = dims;
      newdims.len = 1;

      (void) PyArray_Resize( (PyArrayObject*) array, &newdims,
                             1, NPY_ANYORDER );

      // We know what the size value shoule be (newsize), but
      // check with Python.
      data = PyArray_BYTES( (PyArrayObject*) array );
      size = PyArray_SIZE( (PyArrayObject*) array );
      stride = PyArray_STRIDE( (PyArrayObject*) array, 0 );
    }

    const CType& operator[]( npy_intp index ) const
    {
      return *( static_cast< CType* >
		( static_cast< void* >( data + index*stride ) ) );
    }

    CType& operator[]( npy_intp index )
    {
      return ( const_cast< CType& >
	       ( static_cast< const Array& >( *this )[ index ] ) );
    }

  private:

    int init( PyObject* new_array );

    PyObject* array;
    char* data;
    npy_intp stride;
    npy_intp size;

  };


  template <typename CType, int ArrayType>
  int Array< CType, ArrayType >::from_obj( PyObject* obj,
					   bool contiguous )
  {

    bool decref = false;

    if ( PyArray_Check( obj ) &&
	 !PyArray_CanCastSafely( PyArray_TYPE(  (PyArrayObject*) obj ), ArrayType ) ) {

      // Force downcasting (e.g double->float)
      obj = PyArray_CastToType( (PyArrayObject*)obj,
				PyArray_DescrFromType( ArrayType ), 0 );
      decref = true;

    }

    int rv = init
      ( PyArray_FROMANY( obj, ArrayType, 0, 0,
			 ( contiguous ? NPY_ARRAY_CARRAY : NPY_ARRAY_BEHAVED ) ) );

    if ( decref ) {
      Py_XDECREF( obj );
    }
    return rv;

  }


  template <typename CType, int ArrayType>
  int Array< CType, ArrayType >::init( PyObject* new_array )
  {

    if ( NULL == new_array )
      return EXIT_FAILURE;

    if ( PyArray_NDIM(  (PyArrayObject*) new_array ) > 1 ) {

      PyErr_SetString( PyExc_TypeError,
		       (char*)"array must have 0 or 1 dimensions" );
      Py_DECREF( new_array );
      return EXIT_FAILURE;

    }

    Py_XDECREF( array );
    array = new_array;

    data = PyArray_BYTES(  (PyArrayObject*) array );

    if ( 0 == PyArray_NDIM(  (PyArrayObject*) array ) )
      stride = 0;
    else
      stride = PyArray_STRIDE(  (PyArrayObject*) array, 0 );

    size = PyArray_SIZE( (PyArrayObject*) array );

    return EXIT_SUCCESS;

  }


  template <typename ArrayType>
  int convert_to_array( PyObject* obj, void* dest_addr )
  {

    ArrayType& array = *(static_cast< ArrayType* >( dest_addr ));

    if ( EXIT_SUCCESS != array.from_obj( obj ) )
      return 0;

    return 1;

  }


  template <typename ArrayType>
  int convert_to_contig_array( PyObject* obj, void* dest_addr )
  {

    ArrayType& array = *(static_cast< ArrayType* >( dest_addr ));

    if ( EXIT_SUCCESS != array.from_obj( obj, true ) )
      return 0;

    return 1;

  }


  template <typename ArrayType>
  int array_or_none( PyObject* obj, void* dest_addr )
  {

    if ( Py_None == obj )
      return 1;

    return convert_to_array< ArrayType >( obj, dest_addr );

  }


} /* namespace sherpa */


#endif /* __sherpa_array_hh__ */
