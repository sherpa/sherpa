#ifndef myArray_hh
#define myArray_hh

// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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


#include <iostream>
#include <stdexcept>
#include <vector>

namespace sherpa {

#ifdef NDEBUG

  // production code should run silently!
#define my_cmp_sizes( arg )
#define PRINTARRAY4( os, arg, suffix )
#define PRINTARRAY( arg )
  // production code should run silently!

#else

#define my_cmp_sizes( arg ) this->check_size( arg )
#define STRINGIFY(x) #x
#define PRINTARRAY4( os, arg, num, suffix ) arg.myprint( os, num, STRINGIFY( arg ), suffix );
#define PRINTARRAY( arg ) PRINTARRAY4( std::cout, arg, 6, "\n" )
 
#endif

  ////////////////////////////////// Array1D //////////////////////////////////

  template < typename Real >
  class Array1D : public std::vector< Real > {
    //
    // base class std::vector has a non-virtual destructor 
    // therefore: Thou shall not inherit from std::vector!
    //
  public:

    Array1D( int arg=0 ) : std::vector< Real >( arg ) { }

    ~Array1D( ) { }

    Array1D( const Array1D< Real >& arg ) {
      this->operator=( arg );
    }

    Array1D< Real >& operator = ( const Array1D< Real >& rhs ) {
      if ( this == &rhs )
	return *this;
      std::vector< Real >::operator=( rhs );
      return *this;
    }

    Array1D< Real >& operator = ( Real rhs ) {
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) = rhs;
      return *this;
    }

    const Real& operator [] ( int arg ) const {
#ifdef NDEBUG
      return std::vector< Real >::operator[]( arg );
#else
      return this->at( arg );
#endif
    }

    Real& operator [] ( int arg ) {
#ifdef NDEBUG
      return std::vector< Real >::operator[]( arg );
#else
      return this->at( arg );
#endif
    }

    inline const Real& operator () ( int arg ) const { 
      return this->operator[]( arg );
    }

    inline Real& operator () ( int arg ) { return this->operator[]( arg ); }
    
    Array1D< Real >& operator += ( const Array1D< Real >& rhs ) {
      my_cmp_sizes( rhs );
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) += rhs[ ii ];
      return *this;
    }

    Array1D< Real >& operator -= ( const Array1D< Real >& rhs ) {
      my_cmp_sizes( rhs );
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) -= rhs[ ii ];
      return *this;
    }

    Array1D< Real >& operator *= ( const Array1D< Real >& rhs ) {
      my_cmp_sizes( rhs );
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) *= rhs[ ii ];
      return *this;
    }

    Array1D< Real >& operator /= ( const Array1D< Real >& rhs ) {
      my_cmp_sizes( rhs );
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) /= rhs[ ii ];
      return *this;
    }

    Array1D< Real >& operator += ( Real rhs ) {
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) += rhs;
      return *this;
    }

    Array1D< Real >& operator -= ( Real rhs ) {
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) -= rhs;
      return *this;
    }

    Array1D< Real >& operator *= ( Real rhs ) {
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) *= rhs;
      return *this;
    }

    Array1D< Real >& operator /= ( Real rhs ) {
      for ( size_t ii = 0; ii < rhs.size( ); ++ii )
	this->operator[]( ii ) /= rhs;
      return *this;
    }

    std::ostream& myprint( std::ostream& os, int num=4, const char* prefix=NULL,
			   const char* suffix=NULL ) const {

      if ( this->size( ) > 0 ) {
	if ( NULL != prefix )
	  os << prefix << "[0]=";
	os << this->operator[]( 0 );
      }

      for ( size_t ii = 1; ii < this->size( ); ++ii ) {
	if ( 0 == ii % num )
	  os << '\n';
	else
	  os << ", ";
	if ( NULL != prefix )
	  os << prefix << '[' << ii << "]=";
	os << this->operator[]( ii );
      }

      if ( NULL != suffix )
	os << suffix;

      return os;
    }

  private:

    void check_size( const Array1D< Real >& arg ) const {
      if ( this->size() != arg.size( ) )
	throw std::runtime_error( "Operation attempted on different size arrays" );
    }

  };                                                            // class Array1D

  template< typename Real >
  std::ostream& operator << ( std::ostream& os, const Array1D< Real >& rhs ) {
    return rhs.myprint( os );
  }

  template< typename Real >
  bool operator == ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return std::equal( lhs.begin(), lhs.end(), rhs.begin() );
  }
  template< typename Real >
  bool operator < ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return std::lexicographical_compare( lhs.begin(), lhs.end(),
					 rhs.begin(), rhs.end() );
  }
  template< typename Real >
  bool operator != ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return !( lhs == rhs );
  }
  template< typename Real >
  bool operator > ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return rhs < lhs;
  }
  template< typename Real >
  bool operator <= ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return !( rhs < lhs);
  }
  template< typename Real >
  bool operator >= ( const Array1D< Real >& lhs, const Array1D< Real >& rhs ) {
    return !( lhs < rhs );
  }
  ////////////////////////////////// Array1D //////////////////////////////////

  ////////////////////////////////// Array2D //////////////////////////////////
  //
  // A simple 2d array class written for NelderMead/MultiDirSearch.
  // Note Array2d as written, using vector of vector, does not guaranteed
  // contiguous data. Check out one of the followings if contiguous 
  // in memory is necessary:
  //
  //    http://math.nist.gov/tnt/index.html
  //    http://www.oonumerics.org/blitz/
  //    http://boost.org/libs/multi_array/doc/user.html
  //
  //
  template < typename T >
  class Array2d {

    friend std::ostream& operator << ( std::ostream& os,
				       const Array2d< T >& a ) {
      return a.print( os );
    }

  public:

    virtual ~Array2d( ) { }

    Array2d( int r=0, int c=0 ) : nrow( r ), ncol( c ),
				  array( r, std::vector< T >( c ) ) { }

    std::vector< T >& operator[] ( int arg ) { 
      return array[ arg ];
    }
    const std::vector< T >& operator[] ( int arg ) const { 
      return array[ arg ];
    }

    int ncols( ) const { return ncol; }

    int nrows( ) const { return nrow; }

    static void copy_vector( int num, const std::vector< T >& from,
			     std::vector< T >& to ) {

      for ( int ii = 0; ii < num; ++ii )
	to[ ii ] = from[ ii ];

    }

    std::ostream& print( std::ostream& os ) const {

      for ( int ii = 0; ii < nrow; ++ii ) {
	os << array[ ii ][ 0 ];
	for ( int jj = 1; jj < ncol; ++jj )
	  os << ' ' << array[ ii ][ jj ];
	if ( nrow - 1 != ii )
	  os << '\n';
      }
      return os;
    }

    // resize the two dimentional array
    virtual void resize( int r, int c ) {
      array.resize( r );
      for( int ii = 0; ii < r; ++ii ) 
	array[ ii ].resize( c );
      nrow = r;
      ncol = c;
    }

  protected:

    T& get( int r, int c ) { return array[ r ][ c ]; }
    const T& get( int r, int c ) const { return array[ r ][ c ]; }

    void set( int r, int c, const T& val ) { array[ r ][ c ] = val; }

  private:

    int nrow, ncol;
    std::vector< std::vector< T > > array;

    Array2d& operator = (Array2d const&);  // declare but, purposely, not define
    Array2d( Array2d const& );             // declare but, purposely, not define

  };                                                            // class Array2D
  ////////////////////////////////// Array2D //////////////////////////////////

}                                                            // namespace sherpa

#endif

#ifdef testArray
#include "StopWatch.hh"

int compare( const void* a, const void* b ) {
  int ai = *((int*)a);
  int bi = *((int*)b);
  if( ai < bi )
    return -1;
  else if( ai > bi )
    return 1;
  else 
    return 0;
}

template < typename Type >
void tst_qsort( int num, Type* ptr ) {
  qsort( ptr, num, sizeof( Type ), &compare );
}
template < typename Type >
void tst_qsort( int niter, int alreadysorted, int num,
		sherpa::Array1D< Type >& array ) {

  std::clock_t start = std::clock( );
  for ( int ii = 0; ii < niter; ++ii ) {
    for ( int jj = 0; jj < num; ++jj )
      array[ jj ] = rand();
    tst_qsort( num, &array[0] );
    if ( alreadysorted )
      tst_qsort( num, &array[0] ); // check already sorted
  }
  std::clock_t total = clock( ) - start; 
  std::cout << "qsort took " << double( total ) / CLOCKS_PER_SEC << " secs\n";
}

template < typename Type >
void tst_sort( int num, sherpa::Array1D< Type >& ptr ) {

  sort( ptr.begin(), ptr.end() );

}
template < typename Type >
void tst_sort( int niter, int alreadysorted, int num,
	       sherpa::Array1D< Type >& array ) {

  std::clock_t start = std::clock( );
  for ( int ii = 0; ii < niter; ++ii ) {
    for ( int jj = 0; jj < num; ++jj )
      array[ jj ] = rand();
    tst_sort( num, array );
    if ( alreadysorted )
      tst_sort( num, array );
  }
  std::clock_t total = clock( ) - start; 
  std::cout << "sort took " << double( total ) / CLOCKS_PER_SEC << " secs\n";

}

int main( int argc, char* argv[] ) {

  int seed = 1357;
  int niter = 10;
  int num=1000000;

  int c, alreadysorted = 1;
  while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
    while ( c = *++argv[ 0 ] )
      switch( c ) {
      case 's':
	alreadysorted = 0;
	break;
      default:
	fprintf( stderr, "%s: illegal option '%c'\n", argv[ 0 ], c );
	fprintf( stderr, "Usage %s [ -g ] [ -u ] [ npar ]\n", argv[ 0 ] );
	return EXIT_FAILURE;
      }

  if ( argc == 1 )
    num = atoi( *argv );

  sherpa::Array1D< int > tmp1( num ), tmp2( num );

  std::cout << "num = " << num << "\tniter = " << niter << "\talreadysorted = " << alreadysorted << "\n";
  srand( seed );        // make sure we are always using the same set of numbers
  tst_qsort( niter, alreadysorted, num, tmp1 );
  srand( seed );        // make sure we are always using the same set of numbers
  tst_sort( niter, alreadysorted, num, tmp2 );

  if ( tmp1 != tmp2 )
    std::cerr << "tmp1 != tmp2\n";

  return 0;

}

#endif

#ifdef testArray2d

#include <iostream>

#ifdef testArrayNd
#include "ArrayNd.hh"
#include "Array.h"
#endif
#include "StopWatch.hh"

template < typename Type >
void initLaplacian( Type& uu, int r, int c ) {

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

}

template < typename Type >
void timeme( Type& uu, Type& laplacian, int r, int c, const char* header ) {

  StopWatch stopwatch( header );

  initLaplacian( uu, r, c );

  for( int kk=0; kk < 10; kk++ ) {
    for( int ii = 1; ii < r - 1; ii++ ) {
      for( int jj = 1; jj < c - 1; jj++ ) {
        laplacian[ ii ][ jj ] = - uu[ ii - 1 ][ jj ] - uu[ ii ][ jj - 1 ] +
          4.0 * uu[ ii ][ jj ] - uu[ ii ][ jj + 1 ] - uu[ ii + 1 ][ jj ];
      }
    }
  }

  double sum = 0.0;
  for ( int ii = 1; ii < r - 1; ++ii )
    for ( int jj = 1; jj < c - 1; jj++ )
      sum += laplacian[ ii ][ jj ];

  fprintf( stdout, "%.12f\t", sum );

}

template < typename Type >
Type*** alloc3d( int nx, int ny, int nz ) {

  Type ***ptr = new int** [ nx ];

  for( int x = 0; x < nx; ++x ) {
    ptr[ x ] = new Type* [ ny ];

    for( int y = 0; y < ny; ++ y )
      ptr[ x ][ y ] = new Type[ nz ];
  }

  return ptr;

}

template < typename Type >
void del3d( Type*** ptr, int nx, int ny ) {
  for( int x = 0; x < nx; ++x ) {
    for( int y = 0; y < ny; ++y )
      delete [] ptr[ x ][ y ], ptr[ x ][ y ] = 0;
    delete [] ptr[ x ], ptr[ x ] = 0;
  }
  delete [] ptr;
  ptr = 0;
}

template < typename Type >
Type** alloc2d( int r, int c ) {
  Type** ptr = new Type*[ r ];
  for ( int ii = 0; ii < r; ++ii )
    ptr[ ii ] = new Type[ c ];
  return ptr;
}

template < typename Type >
void del2d( Type** ptr, int r ) {

  for ( int ii = 0; ii < r; ++ii )
    delete [] ptr[ ii ];
  delete [] ptr;

}

void timeclassic( int r, int c ) {

  double** uu = alloc2d< double >( r, c );
  double** laplacian = alloc2d< double >( r, c );

  timeme( uu, laplacian, r, c, "classic" );

  del2d( laplacian, r );
  del2d( uu, r );

}

void timebracket( int r, int c ) {

  sherpa::Array2d< double > uu( r, c ), laplacian( r, c );

  for ( int ii = 0; ii < r; ++ii ) {
    for ( int jj = 0; jj < c; ++jj )
      uu[ ii ][ jj ] = ii * ii + jj * jj;
  }

  timeme( uu, laplacian, r, c, "sherpa::Array2d[][]" );

}

#ifdef testArrayNd
void timemulti_array( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  MultiArray::multi_array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "multi_array[][]" );

}

void timeArray( int r, int c ) {

  std::vector<size_t> dims( 2, r );
  MyArray::Array< double, 2 > laplacian( dims ), uu( dims );

  timeme( uu, laplacian, r, c, "Array[][]" );

}

void timeBavestrelliArray( int r, int c ) {

  bavestrelli::Array<double, 2> laplacian( bavestrelli::ArraySizes(r)(c) ),
    uu( bavestrelli::ArraySizes(r)(c) );

  timeme( uu, laplacian, r, c, "bavestrelli::Array[][]" );

}
void timeJBArray( int r, int c ) {

  Array::array2<double> uu(r,c), laplacian(r,c);

  timeme( uu, laplacian, r, c, "jb::Array[][]" );

}
#endif

int main( int argc, char* argv[] ) {

  const int num = 4000;

  int row = num;
  int col = num;
  if ( argc == 3 ) {
    row = atoi( argv[1] );
    col = atoi( argv[2] );
  }


  timeclassic( row, col );
  timebracket( row, col );
#ifdef testArrayNd
  timemulti_array( row, col );
  timeArray( row, col );
  timeBavestrelliArray( row, col );
  timeJBArray( row, col );
#endif
  return 0;

}

#endif

// To test Array1D:
//
// cp myArray.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArray tmp.cc; rm -f tmp.cc; a.out
// cp myArray.hh tmp.cc; g++ -Wall -ansi -O3 -DtestArray -DNDEBUG tmp.cc; rm -f tmp.cc; a.out
//
// num = 1000000   niter = 10      alreadysorted = 1
// qsort took 2.68 secs
// sort took 1 secs
//
// C++ stl sort is faster then C lib qsort!

// To test Array2D:
//
// cp myArray.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArray2d -DNDEBUG tmp.cc; rm -f tmp.cc; a.out
// cp myArray.hh tmp.cc; g++ -Wall -ansi -pedantic -O3 -DtestArrayNd -DtestArray2d -DNDEBUG tmp.cc; rm -f tmp.cc; a.out
