#include <cmath>
#include <stdexcept>

#include "MulDirSearch.hh"

using namespace sherpa;

template< typename Fct, typename Data >
int MulDirSearch<Fct,Data>::operator( )( double* model_par, int verbose,
					 int initsimplex,
					 std::vector< int >& finalsimplex,
					 double tolerance,
					 const double* step, int maxnfev,
					 int& nfev,
					 double& fmin, double* g, double* h ) {

  try {

    nfev = 0;

    int err_status = EXIT_SUCCESS;

    sherpa::Simplex contraction( sherpa::Opt<Fct,Data>::npar,
				 sherpa::Opt<Fct,Data>::npar + 1 ),
      expansion( sherpa::Opt<Fct,Data>::npar,
		 sherpa::Opt<Fct,Data>::npar + 1 ),
      reflection( sherpa::Opt<Fct,Data>::npar,
		  sherpa::Opt<Fct,Data>::npar + 1 );

    DirectSearch<Fct,Data>::init_simplex( initsimplex, model_par, step );
    err_status =
      DirectSearch<Fct,Data>::eval_init_simplex( maxnfev, model_par, nfev );
    if ( EXIT_SUCCESS != err_status )
      return err_status;

    double tol_sqr = tolerance * tolerance;
    std::vector< double > fctvals( sherpa::Opt<Fct,Data>::npar + 1 );

    for ( ; nfev < maxnfev; ) {

      // 1. Order.
      // The best vertex is labeled as x , so that f(x ) = arg min{ f(x ) }.
      //                                1             1                i
      //
      DirectSearch<Fct,Data>::simplex.order_simplex( );
      int index_smallest =
	DirectSearch<Fct,Data>::simplex.get_smallest_index( );

      // print_vertex( index_smallest, index_smallest );
      if ( DirectSearch<Fct,Data>::simplex.check_convergence( tolerance,
							      tol_sqr,
							      fctvals ) )
	break;

      for ( bool replaced = false; false == replaced; ) {

	// 2. Reflect.
	err_status = 
	  reflect( index_smallest, maxnfev, nfev, reflection );
	if ( EXIT_SUCCESS != err_status )
	  break;

	reflection.order_simplex( );
	int index_reflect = reflection.get_smallest_index( );

	//           (i)
	// If min { f   } < f   go to step 3; otherwise, go to step 4.
	//           f       1
	if ( reflection[ index_reflect ][ sherpa::Opt<Fct,Data>::npar ] < 
	     DirectSearch<Fct,Data>::simplex[ index_smallest ][ sherpa::Opt<Fct,Data>::npar ] ) {

	  // 3. Expand.
	  err_status = expand( index_smallest, maxnfev, nfev, expansion );

	  expansion.order_simplex( );
	  int index_expand = expansion.get_smallest_index( );

	  //            (i)           (i)
	  // If min  { f   } < min { f   }, then accept the expanded,
	  //            e         i   r
	  //
	  //                          (i)
	  // i.e., x  is replaced by x    for i = 2, ..., n + 1; 
	  //        i                 e
	  if ( expansion[ index_expand ][ sherpa::Opt<Fct,Data>::npar ] < 
	       reflection[ index_reflect ][ sherpa::Opt<Fct,Data>::npar ] ) {

	    for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
	      if ( ii != index_smallest ) {
		sherpa::Array2d<double>::Row rowindex = expansion[ index ];
		DirectSearch<Fct,Data>::simplex.copy_row( &rowindex[0], ii );
		++index;
	      }

	  } else {

	    // otherwise, accept the reflected simplex, i.e., x  is
	    //                                                 i
	    //              (i)
	    // replaced by x   .
	    //              r
	    for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
	      if ( ii != index_smallest ) {
		sherpa::Array2d<double>::Row rowindex = reflection[ index ];
		DirectSearch<Fct,Data>::simplex.copy_row( &rowindex[0], ii );
		++index;
	      }

	  }

	  // In either case, terminate the iteration.
	  replaced = true;
	} else {

	  // 4. Contract.
	  err_status = contract( index_smallest, maxnfev, nfev, contraction );
	  if ( EXIT_SUCCESS != err_status )
	    break;

	  //                                   (i)
	  // Terminate the iteration if min { f   } < f ; otherwise, 
	  //                               i   c       1
	  // return to step 2.
	  contraction.order_simplex( );
	  int index_contraction = contraction.get_smallest_index( );
	  replaced = contraction[ index_contraction ][ sherpa::Opt<Fct,Data>::npar ] < DirectSearch<Fct,Data>::simplex[ index_smallest ][ sherpa::Opt<Fct,Data>::npar ] ? true : false;

	}


      }                  // for ( bool replaced = false; false != replaced; ) {

    }                                            // for ( ; nfev < maxnfev; ) {

    DirectSearch<Fct,Data>::simplex.order_simplex( );
    int index_smallest = DirectSearch<Fct,Data>::simplex.get_smallest_index( );
    for ( int ii = 0; ii < sherpa::Opt<Fct,Data>::npar; ++ii )
      model_par[ ii ] =
	DirectSearch<Fct,Data>::simplex[ index_smallest ][ ii ];
    fmin = DirectSearch<Fct,Data>::simplex[ index_smallest ][
					    sherpa::Opt<Fct,Data>::npar ];

    return err_status;

  } catch( std::runtime_error& re ) {
    throw re;
  } catch( std::exception& e ) {
    throw e;
  }

}                                                                // operator( )

template< typename Fct, typename Data >
int MulDirSearch<Fct,Data>::contract( int index_smallest, int maxnfev, int& nfev,
			    sherpa::Array2d< double >& contraction ) {

  int err_status = EXIT_SUCCESS;
  
  //                                     (i)
  // Calculate the contracted vertices, x    = x  + gamma ( x  - x ),
  //                                     c      1            i    1
  //
  // for i = 2, ..., n + 1, and
  for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
    if ( ii != index_smallest ) {
      for ( int jj = 0; jj < sherpa::Opt<Fct,Data>::npar; ++jj )
	contraction[ index ][ jj ] = sherpa::DirectSearch<Fct,Data>::simplex[ index_smallest ][ jj ] +
	  sherpa::DirectSearch<Fct,Data>::contraction_coef * ( sherpa::DirectSearch<Fct,Data>::simplex[ ii ][ jj ] -
			       sherpa::DirectSearch<Fct,Data>::simplex[ index_smallest ][ jj ] );
      ++index;
    }

  //           (i)      (i)
  // evaluate f   = f( x   ).
  //
  for ( int ii = 0; ii < sherpa::Opt<Fct,Data>::npar; ++ii ) {
    sherpa::Array2d<double>::Row rowii = contraction[ ii ];
    eval_user_func( maxnfev, &rowii[0], contraction[ ii ][ sherpa::Opt<Fct,Data>::npar ], nfev,
		    err_status );
    if ( EXIT_SUCCESS != err_status )
      return err_status;
  }

  //                                 (i)     (i)
  // For i = 2, ..., n + 1, replace x    by x   .
  //                                 i       c
  for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
    if ( ii != index_smallest ) {
      sherpa::Array2d<double>::Row rowindex = contraction[ index ];
      DirectSearch<Fct,Data>::simplex.copy_row( &rowindex[0], ii );
      ++index;
    }

  return err_status;

}

template< typename Fct, typename Data >
int MulDirSearch<Fct,Data>::expand( int index_smallest, int maxnfev, int& nfev,
			  sherpa::Array2d< double >& expansion ) {

  int err_status = EXIT_SUCCESS;
  
  //                                 (i)
  // Compute the expanded vertices, x   = x  + chi ( x  - x ),
  //                                 e     1          1    i
  // for i = 2, ..., n + 1, and
  for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
    if ( ii != index_smallest ) {
      for ( int jj = 0; jj < sherpa::Opt<Fct,Data>::npar; ++jj )
	expansion[ index ][ jj ] = sherpa::DirectSearch<Fct,Data>::simplex[ index_smallest ][ jj ] +
	  sherpa::DirectSearch<Fct,Data>::expansion_coef * ( sherpa::DirectSearch<Fct,Data>::simplex[ index_smallest ][ jj ] - 
			     sherpa::DirectSearch<Fct,Data>::simplex[ ii ][ jj ] );
      ++index;
    }

  //           (i)       (i)
  // evaluate f    = f( x   ).
  //           e         e
  for ( int ii = 0; ii < sherpa::Opt<Fct,Data>::npar; ++ii ) {
    sherpa::Array2d<double>::Row rowii = expansion[ ii ];
    eval_user_func( maxnfev, &rowii[0], expansion[ ii ][ sherpa::Opt<Fct,Data>::npar ], nfev,
		    err_status );
    if ( EXIT_SUCCESS != err_status )
      return err_status;
  }

  return err_status;

}

template< typename Fct, typename Data >
int MulDirSearch<Fct,Data>::reflect( int index_smallest, int maxnfev, int& nfev,
			   sherpa::Array2d< double >& reflection ) {

  int err_status = EXIT_SUCCESS;

  //                                  (i)
  // Define the n reflected vertices x   = 2 x - x , for i = 2, ..., n+1,
  //                                  r       1   i
  for ( int ii = 0, index = 0; ii <= sherpa::Opt<Fct,Data>::npar; ++ii )
    if ( ii != index_smallest ) {
      for ( int jj = 0; jj < sherpa::Opt<Fct,Data>::npar; ++jj )
	reflection[ index ][ jj ] = 2.0 * sherpa::DirectSearch<Fct,Data>::simplex[ index_smallest ][ jj ] -
	  sherpa::DirectSearch<Fct,Data>::simplex[ ii ][ jj ];
      ++index;
    }

  //            (i)       (i)
  // Evaluate f     = f( x   )
  //           r          r 
  for ( int ii = 0; ii < sherpa::Opt<Fct,Data>::npar; ++ii ) {
    sherpa::Array2d<double>::Row rowii = reflection[ ii ];
    eval_user_func( maxnfev, &rowii[0], reflection[ ii ][ sherpa::Opt<Fct,Data>::npar ], nfev,
		    err_status );
    if ( EXIT_SUCCESS != err_status )
      return err_status;
  }

  return err_status;

}

#ifdef testMulDirSearch

#include <stdexcept>

#include "fcmp.h"
#include "tstoptfct.hh"
#include "sherpa/functor.hh"

template <typename Real>
void print_pars( const char* name, int nfev, Real stat, Real answer,
		 int n, const std::vector< Real >& x,
		 Real tol=
		 1.0e4*std::sqrt( std::numeric_limits< Real >::epsilon() ) ) {

 
  std::cout << "MulDirSearch_" << name << '\t';
  if ( 0 == _sao_fcmp( stat, answer, std::sqrt(tol) ) )
    std::cout << nfev << '\t';
  else
    std::cout << -nfev << '\t';
  std::cout << answer << '\t';
  std::cout << stat << '\t';
  std::cout << x[0];
  for ( int ii = 1; ii < n; ++ii )
    std::cout << ',' << x[ii];
  std::cout << '\n';

}

template< typename Init, typename Fct >
void justdoit( Init init, Fct fct, int npars, 
	       std::vector< double >& pars, std::vector< double >& lo,
	       std::vector< double >& hi, double tol, const char* header ) {

  int mfcts;
  double answer;
  double xprob = 0.9, scale=1.0;
  int begin_strategy = 0, end_strategy = 10;

  init( npars, mfcts, answer, &pars[0], &lo[0], &hi[0] );

  sherpa::MulDirSearch< Fct, void* > mds( npars, &pars[0], &lo[0], &hi[0],
					  fct, NULL );

  int nfev, verbose=0, initsimplex=0, maxnfev=npars*256;
  double fmin;
  std::vector< int > finalsimplex( 2, 1 );
  std::vector< double > step( npars, 1.2 );
  mds( &pars[0], verbose, initsimplex, finalsimplex, tol, &step[0],
       maxnfev, nfev, fmin, NULL, NULL );
    
  print_pars( header, nfev, fmin, answer, npars, pars );

}

void tstuncopt( int npars, double tol ) {

  const int verbose=0, size=npars*32;
  std::vector< double > par( size, 0 ), lo( size, -1.0e2 ), hi( size, 1.0e2 );

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Rosenbrock<double,void*> );

    justdoit( fct_ptr( tstoptfct::RosenbrockInit<double> ), fct, npars, par,
	      lo, hi, tol, "Rosenbrock" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::FreudensteinRoth<double,void*> );

    justdoit( fct_ptr( tstoptfct::FreudensteinRothInit<double> ),
	      fct, npars, par, lo, hi, tol, "FreudensteinRoth" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PowellBadlyScaled<double,void*> );

    justdoit( fct_ptr( tstoptfct::PowellBadlyScaledInit<double> ),
	      fct, npars, par, lo, hi, tol, "PowellBadlyScaled" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownBadlyScaled<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownBadlyScaledInit<double> ),
	      fct, 2, par, lo, hi, tol, "BrownBadlyScaled" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Beale<double,void*> );

    justdoit( fct_ptr( tstoptfct::BealeInit<double> ),
	      fct, npars, par, lo, hi, tol, "Beale" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::JennrichSampson<double,void*> );

    justdoit( fct_ptr( tstoptfct::JennrichSampsonInit<double> ),
	      fct, npars, par, lo, hi, tol,
	      "JennrichSampson" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::HelicalValley<double,void*> );

    justdoit( fct_ptr( tstoptfct::HelicalValleyInit<double> ),
	      fct, 3, par, lo, hi, tol, "HelicalValley" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Bard<double,void*> );

    justdoit( fct_ptr( tstoptfct::BardInit<double> ),
	      fct, 3, par, lo, hi, tol, "Bard" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Gaussian<double,void*> );

    justdoit( fct_ptr( tstoptfct::GaussianInit<double> ),
	      fct, 3, par, lo, hi, tol, "Gaussian" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Meyer<double,void*> );

    justdoit( fct_ptr( tstoptfct::MeyerInit<double> ),
	      fct, 3, par, lo, hi, tol, "Meyer" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::GulfResearchDevelopment<double,void*> );

    justdoit( fct_ptr( tstoptfct::GulfResearchDevelopmentInit<double> ),
	      fct, 3, par, lo, hi, tol, "GulfResearchDevelopment" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Box3d<double,void*> );

    justdoit( fct_ptr( tstoptfct::Box3dInit<double> ),
	      fct, 3, par, lo, hi, tol, "Box3d" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PowellSingular<double,void*> );

    justdoit( fct_ptr( tstoptfct::PowellSingularInit<double> ),
	      fct, 4, par, lo, hi, tol, "PowellSingular" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Wood<double,void*> );

    justdoit( fct_ptr( tstoptfct::WoodInit<double> ),
	      fct, 4, par, lo, hi, tol, "Wood" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::KowalikOsborne<double,void*> );

    justdoit( fct_ptr( tstoptfct::KowalikOsborneInit<double> ),
	      fct, 4, par, lo, hi, tol, "KowalikOsborne" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownDennis<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownDennisInit<double> ),
	      fct, 4, par, lo, hi, tol, "BrownDennis" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Osborne1<double,void*> );

    justdoit( fct_ptr( tstoptfct::Osborne1Init<double> ),
	      fct, 5, par, lo, hi, tol, "Osborne1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Biggs<double,void*> );

    justdoit( fct_ptr( tstoptfct::BiggsInit<double> ),
	      fct, 6, par, lo, hi, tol, "Biggs" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Osborne2<double,void*> );

    justdoit( fct_ptr( tstoptfct::Osborne2Init<double> ),
	      fct, 11, par, lo, hi, tol, "Osborne2" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Watson<double,void*> );

    justdoit( fct_ptr( tstoptfct::WatsonInit<double> ),
	      fct, 6, par, lo, hi, tol, "Watson" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PenaltyI<double,void*> );

    justdoit( fct_ptr( tstoptfct::PenaltyIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyI" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::PenaltyII<double,void*> );

    justdoit( fct_ptr( tstoptfct::PenaltyIIInit<double> ),
	      fct, 4, par, lo, hi, tol, "PenaltyII" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::VariablyDimensioned<double,void*> );

    justdoit( fct_ptr( tstoptfct::VariablyDimensionedInit<double> ),
	      fct, npars, par, lo, hi, tol, "VariablyDimensioned" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Trigonometric<double,void*> );

    justdoit( fct_ptr( tstoptfct::TrigonometricInit<double> ),
	      fct, npars, par, lo, hi, tol, "Trigonometric" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BrownAlmostLinear<double,void*> );

    justdoit( fct_ptr( tstoptfct::BrownAlmostLinearInit<double> ),
	      fct, npars, par, lo, hi, tol, "BrownAlmostLinear" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::DiscreteBoundary<double,void*> );

    justdoit( fct_ptr( tstoptfct::DiscreteBoundaryInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteBoundary" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::DiscreteIntegral<double,void*> );

    justdoit( fct_ptr( tstoptfct::DiscreteIntegralInit<double> ),
	      fct, npars, par, lo, hi, tol, "DiscreteIntegral" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BroydenTridiagonal<double,void*> );

    justdoit( fct_ptr( tstoptfct::BroydenTridiagonalInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenTridiagonal" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::BroydenBanded<double,void*> );

    justdoit( fct_ptr( tstoptfct::BroydenBandedInit<double> ),
	      fct, npars, par, lo, hi, tol, "BroydenBanded" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRankInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank1<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRank1Init<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank1" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::LinearFullRank0cols0rows<double,void*> );

    justdoit( fct_ptr( tstoptfct::LinearFullRank0cols0rowsInit<double> ),
	      fct, npars, par, lo, hi, tol, "LinearFullRank0cols0rows" );
  }

  {
    FctPtr< void, int, double*, double&, int&, void* >
      fct( tstoptfct::Chebyquad<double,void*> );

    justdoit( fct_ptr( tstoptfct::ChebyquadInit<double> ),
	      fct, 9, par, lo, hi, tol, "Chebyquad" );
  }

  return;

}

int main( int argc, char* argv[] ) {

  try {

    int c, passit = 0;
    while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
      while ( c = *++argv[ 0 ] )
	switch( c ) {
	case 'p':
	  passit = 1;
	  break;
	default:
	  fprintf( stderr, "%s: illegal option '%c'\n", argv[ 0 ], c );
	  fprintf( stderr, "Usage %s [ -p ] [ npars ]\n", argv[ 0 ] );
	  return EXIT_FAILURE;
      }


    int npars=6;
    if ( argc == 1 )
      npars = atoi( *argv );
    
    if ( npars % 2 || npars < 2 ) {
      printf( "The minimum value for the free parameter must be an even "
	      "and it is greater then 2\n" );
      return EXIT_FAILURE;
    }


    std::vector< int > finalsimplex( 1 );
    finalsimplex[0] = 1;
    int nfev=0, mfcts=2, initsimplex=0;
    int maxnfev = npars * npars * 4096, verbose=0;
    double fmin = 0.0, tol = 1.0e-8, answer;

    int size = npars * npars * 4;
    std::vector< double > pars( size ), step( size ), lo( size ),
      hi( size );

    for ( int ii = 0; ii < npars; ++ii )
      step[ ii ] = 0.4;

    std::cout << "#:tol=" << tol << '\n';
    std::cout << "# A negative value for the nfev signifies that the "
      "optimization method did not converge\n#\n";
    std::cout << "name\tnfev\tanswer\tstat\tpars\nS\tN\tN\tN\tN\n";

    tstuncopt( npars, tol );

    return 0;

  } catch( std::exception& e ) {

    std::cerr << e.what( ) << '\n';
    return 1;

  }

}

/*
gcc -g  -Wall -pedantic -ansi -c ../../utils/src/gsl/fcmp.c
g++ -g  -Wall -pedantic -ansi -c -I../../include/ -I../../utils/src/gsl DirectSearch.cc
g++ -g -Wall -pedantic -ansi -c -O3 -I../../include/ -I../../utils/src/gsl Simplex.cc
g++ -g -Wall -pedantic -ansi -I../../include/ -I../tests -I../../utils/src/gsl -DtestMulDirSearch MulDirSearch.cc DirectSearch.o Simplex.o fcmp.o
valgrind --tool=memcheck --leak-check=yes --show-reachable=yes a.out
==22409== Memcheck, a memory error detector for x86-linux.
==22409== Copyright (C) 2002-2004, and GNU GPL'd, by Julian Seward et al.
==22409== Using valgrind-2.2.0, a program supervision framework for x86-linux.
==22409== Copyright (C) 2000-2004, and GNU GPL'd, by Julian Seward et al.
==22409== For more details, rerun with: -v
==22409==
#:tol=1e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name    nfev    answer  stat    pars
S       N       N       N       N
MulDirSearch_Rosenbrock 147457  0       0.00613462      0.955176,0.91214,0.967174,0.935322,0.944886,0.892601
MulDirSearch_FreudensteinRoth   -7507   0       146.955 11.3717,-0.899234,11.3636,-0.899769,11.364,-0.899722
MulDirSearch_PowellBadlyScaled  1555    0       9.73129e-09     1.05757e-05,9.45562,1.01709e-05,9.83196,1.16402e-05,8.5909
MulDirSearch_BrownBadlyScaled   938     0       3.33455e-11     1e+06,1.99999e-06
MulDirSearch_Beale      13411   0       1.40412e-05     2.99514,0.49875,2.9931,0.498293,2.99607,0.498986
MulDirSearch_JennrichSampson    3091    373.086 373.087 0.257818,0.257831,0.257974,0.257672,0.257723,0.257924
MulDirSearch_HelicalValley      19234   0       1.06026e-05     1,0.00205256,0.00325314
MulDirSearch_Bard       7546    0.00821487      0.00822208      0.0833934,1.16526,2.31277
MulDirSearch_Gaussian   -94     1.12793e-08     1.06248e-07     0.399114,1.00113,0
MulDirSearch_Meyer      -147457 87.9458 848562  0.02,5006.16,298.058
MulDirSearch_GulfResearchDevelopment    16636   0       2.35182e-07     2.59537,40.2849,0.79618
MulDirSearch_Box3d      -1144   0       0.0906338       386.726,0.77728,-1.2392
MulDirSearch_PowellSingular     5189    0       1.83411e-05     0.0481079,-0.00481118,0.0297988,0.0298687
MulDirSearch_Wood       138181  0       0.000315773     0.990713,0.981475,1.00938,1.01889
MulDirSearch_KowalikOsborne     -4205   0.000307505     0.000353182     0.187826,0.407499,0.217093,0.22412
MulDirSearch_BrownDennis        853     85822.2 85822.2 -11.5923,13.2029,-0.403557,0.236705
MulDirSearch_Osborne1   -10776  5.46489e-05     0.000181403     0.358441,1.25255,-0.766925,0.0105316,0.0278215
MulDirSearch_Biggs      147457  0       2.99441e-33     1,10,1,5,4,3
MulDirSearch_Osborne2   -6612   0.0401377       0.446397        1.08132,0.0715816,0.44248,0.589331,0.286894,9.91299,5,7,2,4.5,5.5
MulDirSearch_Watson     -137659 0.00228767      0.00252001      -0.0223851,1.00553,-0.162291,1.01996,-1.22388,0.865028
MulDirSearch_PenaltyI   -237    2.24997e-05     3.36137e-05     0.119206,-0.217437,0.127521,0.415026
MulDirSearch_PenaltyII  -197    9.37629e-06     9.65518e-06     0.200086,-0.0710923,0.512439,0.547212
MulDirSearch_VariablyDimensioned        3295    0       4.88261e-07     0.999814,0.99991,0.999394,1.00019,1.00019,1.00007
MulDirSearch_Trigonometric      6463    0       1.49431e-06     1.29047,0.429007,0.327551,0.277009,0.245061,0.222426
MulDirSearch_BrownAlmostLinear  -12403  1       6.55632e-06     0.938981,0.939062,0.939395,0.939319,0.939024,1.36545
MulDirSearch_DiscreteBoundary   4711    0       4.67695e-06     -0.0919674,-0.171894,-0.157053,-0.135171,-0.103965,-0.06027
MulDirSearch_DiscreteIntegral   499     0       2.54001e-08     -0.0654367,-0.118258,-0.154585,-0.170013,-0.157339,-0.106096
MulDirSearch_DiscreteBoundary   5059    0       6.97475e-06     -0.0922681,-0.172597,-0.157948,-0.135862,-0.104621,-0.0606138
MulDirSearch_BroydenTridiagonal 2815    0       3.47055e-07     -0.576107,-0.695996,-0.680344,-0.643073,-0.556548,-0.366102
MulDirSearch_BroydenBanded      787     0       1.75071e-08     -0.428313,-0.476598,-0.519666,-0.558073,-0.593447,-0.593447
MulDirSearch_LinearFullRank     859     0       1.76755e-08     -1.00006,-0.999943,-1.00004,-0.999928,-0.999964,-1.00005
MulDirSearch_LinearFullRank1    283     1.15385 1.15385 1,1,1,-0.423626,0.161128,-0.813394
MulDirSearch_LinearFullRank0cols0rows   295     2.66667 2.66667 1,1.1206,1,0.192902,-1.13589,1
MulDirSearch_Chebyquad  3739    0       5.32477e-06     0.0672183,0.365938,0.290552,0.635054,0.709738,0.932999
==22409==
==22409== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 15 from 1)
==22409== malloc/free: in use at exit: 0 bytes in 0 blocks.
==22409== malloc/free: 836697 allocs, 836697 frees, 88000272 bytes allocated.
==22409== For counts of detected errors, rerun with: -v
==22409== No malloc'd blocks -- no leaks are possible.
*/
#endif
