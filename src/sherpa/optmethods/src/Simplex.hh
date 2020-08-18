#ifndef Simplex_hh
#define Simplex_hh

#include <cmath>
#include <vector>

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


#include "sherpa/myArray.hh"

namespace sherpa {

  class Simplex : public sherpa::Array2d< double > {

  public:

    Simplex( int r=0, int c=0 ) : sherpa::Array2d< double >( r, c ),
				  key( r, 0.0 ) { }

    static double
    calc_standard_deviation_square( int num, const std::vector<double>& ptr );

    // copy the c-th col of the Array2d to the pointer 'cpto'
    void copy_col( int c, std::vector< double >& cpto ) {
      const int myncol = ncols( );
      if ( c < 0 || c >= static_cast< int >( myncol ) )
	throw std::runtime_error( "index out of bounds" );
      const int mynrow = nrows( );
      for ( int ii = 0; ii < mynrow; ++ii )
	cpto[ ii ] = this->operator[]( ii )[ c ];
    }

    // copy from the pointer 'cpfrom' to the c-th col of the array
    void copy_col( const std::vector< double >& cpfrom, int c ) {
      const int myncol = ncols( );
      if ( c < 0 || c >= static_cast< int >( myncol ) )
	throw std::runtime_error( "index out of bounds" );
      const int mynrow = nrows( );
      for ( int ii = 0; ii < mynrow; ++ii )
	this->operator[]( ii )[ c ] = cpfrom[ ii ];
    }

    // copy the r-th row of the Array2d to the pointer 'cpto'
    void copy_row( int r, std::vector< double >& cpto ) {
      const int mynrow = nrows( );
      if ( r < 0 || r >= static_cast< int >( mynrow ) )
	throw std::runtime_error( "index out of bounds" );
      sherpa::Array2d<double>::copy_vector( ncols( ), this->operator[]( r ),
					    cpto );
      /*
      const int myncol = ncols( );
      for ( int ii = 0; ii < myncol; ++ii )
	cpto[ ii ] = this->operator[]( r )[ ii ];
      */
    }

    // Copy from 'cpfrom' to the r-th row of the array
    void copy_row( const std::vector< double >& cpfrom, int r ) {
      const int mynrow = nrows( );
      if ( r < 0 || r >= static_cast< int >( mynrow ) )
	throw std::runtime_error( "index out of bounds" );
      /*
      const int myncol = ncols( );
      for ( int ii = 0; ii < myncol; ++ii )
	this->operator[]( r )[ ii ] = cpfrom[ ii ];
      */
      sherpa::Array2d<double>::copy_vector( ncols( ), cpfrom,
					    this->operator[]( r ) );
    }

    bool check_convergence( double tolerance, double tol_sqr,
			    int finalsimplex=0 );

    void init_simplex( int initsimplex, const std::vector<double>& par,
		       const std::vector<double>& step ) {

      const int npar = npars( );
      std::vector<double> mystep( npar + 1 );
      check_step( npar, step, mystep );

      for ( int ii = 0; ii < npar; ++ii )
	set( 0, ii, par[ ii ] );

      switch( initsimplex ) {
      case 1:
	SpendleyHextHimsworth_simplex( mystep, par );
	break;
      default:
	dtn_simplex( mystep, par );
      }

    }

    int npars( ) const { return ncols( ) - 1; }

    void print_simplex( ) const;

    void resize( int r, int c ) {
      key.resize( r );
      sherpa::Array2d<double>::resize( r, c );
    }

    void sort( );

  protected:

    void check_step( int npar, const std::vector<double>& step,
		     std::vector<double>& mystep ) {

      int allzero = 0;
      for ( int ii = 0; ii < npar; ++ii ) {
	mystep[ ii ] = step[ ii ];
	if ( 0.0 == step[ ii ] )
	  ++allzero;
      }
      if ( npar == allzero )
	for ( int ii = 0; ii < npar; ++ii )
	  mystep[ ii ] = 1.0;

    }                                                             // check_step

  private:

    std::vector<double> key;

    bool are_fct_vals_close_enough( double tolerance ) const;

    bool is_max_length_small_enough( double tolerance ) const;

    bool is_stddev_small_enough( double tolerance, double tol_sqr );

    Simplex& operator = (Simplex const&); // declare but, purposely, not define
    Simplex( Simplex const& );            // declare but, purposely, not define

    void dtn_simplex( const std::vector<double>& step,
		      const std::vector<double>& par ) {

      const int npar = npars( );
      for ( int ii = 0; ii < npar; ++ii ) {
	for ( int jj = 0; jj < npar; ++jj )
	  set( ii + 1, jj, par[ jj ] );
	set( ii + 1, ii, par[ ii ] + step[ ii ] );
      }

    }                                                            // dtn_simplex
    
    void print_vertex( std::ostream& os, size_t npar,
		       const std::vector<double>& vertex ) const;
    
    //
    // construct a regular simplex, see: Spendley, Hext and Himsworth,
    // "Sequential Application of Simplex Designs in Optimization and
    // "Evolutionary Operation", Technometics, Vol 4, No 4, Nov 1962
    // pages 441-461.
    //
    void SpendleyHextHimsworth_simplex( const std::vector<double>& step,
					const std::vector<double>& par ) {

      const int npar = npars( );
      double nparsqrt2 = npar * sqrt( 2.0 );
      double sqrtnpar1 = std::sqrt( double( npar + 1 ) );
      double pn = ( sqrtnpar1 - 1 + npar ) / nparsqrt2;
      double qn = ( sqrtnpar1 - 1 ) / nparsqrt2;
	
      for ( int ii = 1; ii <= npar; ++ii )
	for (int jj = 0; jj < npar; ++jj )
	  if ( ii - 1 == jj )
	    set( ii, jj, pn + par[ jj ] );
	  else
	    set( ii, jj, qn + par[ jj ] );

    }                                          // SpendleyHextHimsworth_simplex

  };                                                           // class Simplex

  class mySimplex : public sherpa::Simplex {

  public:

    mySimplex( int r, int c, const std::vector<double>& step,
	       const std::vector<double>& par ) : sherpa::Simplex( r, c ) {
      const int npar = c - 1;
      std::vector<double> mystep( npar, 1.0 );
      check_step( npar, step, mystep );

      for ( int ii = 0; ii < npar; ++ii )
	set( 0, ii, par[ ii ] );
      
      for ( int ii = 0; ii < npar; ++ii ) {
	for ( int jj = 0; jj < npar; ++jj )
	  set( ii + 1, jj, par[ jj ] );
	set( ii + 1, ii, par[ ii ] + mystep[ ii ] );
      }
      
    }                                                  // mySimplex constructor

  private:

    // declare but, purposely, not define
    mySimplex& operator = (mySimplex const&);
    mySimplex( mySimplex const& );        // declare but, purposely, not define

  };                                                         // class mySimplex

  class SpendleyHextHimsworthSimplex: public sherpa::Simplex {

  public:

    SpendleyHextHimsworthSimplex( int r, int c,
				  const std::vector<double>& par )
      : sherpa::Simplex( r, c ) {
	
      const int npar = c - 1;
      double nparsqrt2 = npar * sqrt( 2.0 );
      double sqrtnpar1 = std::sqrt( double( npar + 1 ) );
      double pn = ( sqrtnpar1 - 1 + npar ) / nparsqrt2;
      double qn = ( sqrtnpar1 - 1 ) / nparsqrt2;
	
      for ( int ii = 0; ii < npar; ++ii )
	set( 0, ii, par[ ii ] );

      for ( int ii = 1; ii <= npar; ++ii )
	for (int jj = 0; jj < npar; ++jj )
	  if ( ii - 1 == jj )
	    set( ii, jj, pn + par[ jj ] );
	  else
	    set( ii, jj, qn + par[ jj ] );
	
    }

  private:

    // declare but, purposely, not define
    SpendleyHextHimsworthSimplex& operator = 
      (SpendleyHextHimsworthSimplex const&);
    // declare but, purposely, not define
    SpendleyHextHimsworthSimplex( SpendleyHextHimsworthSimplex const& );
      
  };                                            // SpendleyHextHimsworthSimplex

}                                                                  // namespace

#endif

#ifdef testSimplex

class Junk {
public:

  Junk( int r, int c, const std::vector<double>& step,
	const std::vector<double>& par ) {

    simplex.reset( new sherpa::mySimplex( r, c, step, par ) );

    std::cout << simplex->get( r - 1, c - 1 ) << '\n';

  }

private:

  std::auto_ptr< sherpa::Simplex > simplex;

};

int main( ) {

  int nrow = 14;
  int ncol = 5;
  std::vector<double> step( ncol, 0.0 ), par( ncol, 0.0 );
  Junk( nrow, ncol, step, par );

}

#endif
//
// cp Simplex.hh tmp.cc; g++ -g -Wall -ansi -pedantic -O3 -DtestSimplex tmp.cc; rm tmp.cc; myvalgrind a.out
//
