#include <cmath>
#include <cstdio>

#include "sherpa/fcmp.hh"

#include "Simplex.hh"

using namespace sherpa;

bool Simplex::are_fct_vals_close_enough( double tolerance ) const {

  
  const int npar = npars( );
  if ( 0 == sao_fcmp( get( 0, npar ), get( nrows( ) - 1, npar ), tolerance ) )
    return true;

  return false;

}

double Simplex::calc_standard_deviation_square( int num,
						const std::vector<double>& ptr ) {

  //
  // The standard deviation algorithm is due to Donald E. Knuth (1998).
  // The art of Computer Programming, volume2: Seminumerical Algorithms,
  // 3rd edn., p 232
  //
  double mean = 0.0, stddev = 0.0;
  for ( int ii = 0; ii < num; ++ii ) {
    double delta = ptr[ ii ] - mean;
    mean += delta / double( ii + 1 );
    stddev += delta * ( ptr[ ii ] - mean );
  }

  if ( 1 != num )
    stddev /= double( num - 1);

  return stddev;

}

bool Simplex::check_convergence( double tolerance, double tol_sqr,
				 int finalsimplex ) {

  switch( finalsimplex ) {
  case 0:
    if ( false == is_max_length_small_enough( tolerance ) )
      return false;
    return true;
  case 2:
    {
      if ( false == is_max_length_small_enough( tolerance ) )
	return false;
      bool stddev = is_stddev_small_enough( tolerance, tol_sqr );
      bool fctval = are_fct_vals_close_enough( tolerance );
      return stddev && fctval;
    }
    break;
  default:
    {
      if ( false == is_max_length_small_enough( tolerance ) )
	return false;
      bool stddev = is_stddev_small_enough( tolerance, tol_sqr );
      bool fctval = are_fct_vals_close_enough( tolerance );
      return stddev || fctval;
    }
  }

  return false;

}                                                          // check_convergence

bool Simplex::is_max_length_small_enough( double tol )  const {

  const int index_smallest = 0;
  const int npar = npars( );
  double maxof_x_i_minus_x_min = -1.0; // norm is always a positive number.
  for ( int ii = 0; ii <= npar; ++ii ) {
    double tmp = 0.0;
    if ( ii != index_smallest )
      for ( int jj = 0; jj < npar; ++jj )
	tmp += ( get( ii, jj ) - get( index_smallest, jj ) ) *
	  ( get( ii, jj ) - get( index_smallest, jj ) );

    maxof_x_i_minus_x_min = std::max( maxof_x_i_minus_x_min, tmp );
  }
  double norm_min = 0.0;
  for ( int ii = 0; ii < npar; ++ii )
    norm_min += get( index_smallest, ii ) * get( index_smallest, ii );
  norm_min = norm_min > 1.0 ? norm_min : 1.0;
  if ( maxof_x_i_minus_x_min <= tol * norm_min )
    return true;

  return false;

}                                                 // is_max_length_small_enough

bool Simplex::is_stddev_small_enough( double tolerance, double tol_sqr ) {

  const int npar = npars( );
  this->copy_col( npar, key );
  double std_dev_sqr = calc_standard_deviation_square( ncols( ), key );
  if ( sao_fcmp( std_dev_sqr, tol_sqr, tolerance ) <= 0 )
    return true;

  return false;

}                                                     // is_stddev_small_enough

void Simplex::print_simplex( ) const {

  const int npar = npars( );
  for ( int ii = 0; ii <= npar; ++ii ) {
    const std::vector<double>& vertex = this->operator[ ]( ii );
    print_vertex( std::cout, npar, vertex );
  }

}                                                              // print_simplex

//
// This method is redundant, call Opt::print_par instead
//
void Simplex::print_vertex( std::ostream& os, size_t npar, 
			    const std::vector<double>& vertex ) const {
  os.precision( 6 );
  os << "f( " << std::scientific << vertex[ 0 ];
  for ( size_t ii = 1; ii < npar; ++ii )
    os << ", " << std::scientific << vertex[ ii ];
  os << " ) = " << vertex[ npar ] << '\n';
  return;
}                                                               // print_vertex

//
// This method should be deprecated now that Array2d has been
// switched to use vector of vector, one should use one of stl
// sort routines instead.
//
void Simplex::sort( ) {

  const int mynrow = nrows( );
  const int myncol = ncols( );
  const int fvalindex = npars( );

  for ( int jj = 1; jj < mynrow; ++jj ) {

    for ( int ii = 0; ii < myncol; ++ii )
      key[ ii ] = get( jj, ii );

    int ii = jj;
    for ( ; ii > 0 && get( ii - 1, fvalindex ) > key[ fvalindex ]; --ii ) {

      for ( int kk = 0; kk < myncol; ++kk )
	set( ii, kk, get( ii - 1, kk ) );

    }

    for ( int kk = 0; kk < myncol; ++kk )
      set( ii, kk, key[ kk ] );

  }

}                                                                       // sort
