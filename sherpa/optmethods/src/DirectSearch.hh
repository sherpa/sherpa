#ifndef DirectSearch_hh
#define DirectSearch_hh

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


//
// The abstract base class for the Nelder-Mead and Multi-Directional-Search
// algorithms.  The algorithms are described in the following two papers:
//
// Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright, Paul E. Wright
// "Convergence Properties of the Nelder-Mead Simplex Algorithm in Low
// Dimensions", SIAM Journal on Optimization,Vol. 9, No. 1 (1998),
// pages 112-147.
// http://citeseer.ist.psu.edu/3996.html
//
// Wright, M. H. (1996) "Direct Search Methods: Once Scorned, Now Respectable"
// in Numerical Analysis 1995 (Proceedings of the 1995 Dundee Biennial
// Conference in Numerical Analysis) (D.F. Griffiths and G.A. Watson, eds.),
// 191-208, Addison Wesley Longman, Harlow, United Kingdom.
// http://citeseer.ist.psu.edu/155516.html
//
// Jan 2008 D. T. Nguyen
//

#include <cmath>

#include <fcmp.h>

#include "Opt.hh"
#include "Simplex.hh"

namespace sherpa {

  template< typename Func, typename Data >
  class DirectSearch : public sherpa::OptFunc< Func, Data > {

  public:

    DirectSearch( int numpar, double* par, const double* lo, const double* hi,
		  Func func, Data xtra, double contraction=0.5,
		  double expansion=2.0, double reflection=1.0 )
      : sherpa::OptFunc< Func, Data >( numpar, par, lo, hi, func, xtra ),
	simplex( numpar + 1, numpar + 1 ), contraction_coef( contraction ),
	expansion_coef( expansion ), reflection_coef( reflection ) {
    
	check_coefficients( );

    }
    //
    // The nearly universal choices used in the standard Nelder-Mead
    // algorithm are
    //
    // reflection(rho)=1.0, expansion(chi)=2.0, 
    // contraction(gamma)=0.5, shrink(sigma)=0.5
    // 
    //

    virtual int operator( )( double* model_par, int verbose, int initsimplex,
			     std::vector< int >& finalsimplex, 
			     double tolerance, const double* step, int maxnfev,
			     int& nfev, double& fmin ) = 0;

  protected:

    //
    // simplex[ 0...npar ][ 0...npar-1 ] contains the npar + 1 vertices
    // of the simplex, the trial parameters.
    // simplex[ 0...npar ][ npar ] contains the function values of the simplex
    //
    sherpa::Simplex simplex;

    const double contraction_coef, expansion_coef, reflection_coef;

    int eval_init_simplex( int maxnfev, double* model_par, int& nfev ) {

      //
      // make sure the initial simplex vertices are within bounds.
      //

      /*
      for ( int ii = 1; ii <= npar; ++ii )
	for ( int jj = 0; jj < npar; ++jj ) {
	  if ( simplex[ ii ][ jj ] < Opt<Func,Data>::lo_bound[ jj ] )
	    simplex[ ii ][ jj ] = Opt<Func,Data>::lo_bound[ jj ] *
	    (1.0 + myrandom( 0.0, 0.5 ) );
	  if ( simplex[ ii ][ jj ] > Opt<Func,Data>::hi_bound[ jj ] )
	    simplex[ ii ][ jj ] = Opt<Func,Data>::hi_bound[ jj ] *
	    (1.0 - myrandom( 0.0, 0.5 ) );
	}
      */

      for ( int ii = 1; ii < NPAR; ++ii )

	for ( int jj = 0; jj < NPAR; ++jj ) {
      
	  if ( simplex[ ii ][ jj ] < Opt<Func,Data>::lo_bound[ jj ] )
	    if ( Opt<Func,Data>::hi_bound[ jj ] -
		 Opt<Func,Data>::lo_bound[ jj ] < 10.0 )
	      simplex[ ii ][ jj ] =
		Opt<Func,Data>::lo_bound[ jj ] +
		( Opt<Func,Data>::hi_bound[ jj ] -
		  Opt<Func,Data>::lo_bound[ jj ] ) / 4.0;
	    else
	      simplex[ ii ][ jj ] =
		std::min( model_par[ ii ] + 0.01 * fabs( model_par[ ii ] ),
			  Opt<Func,Data>::hi_bound[ jj ] );
	  
	  if ( simplex[ ii ][ jj ] > Opt<Func,Data>::hi_bound[ jj ] )
	    if ( Opt<Func,Data>::hi_bound[ jj ] -
		 Opt<Func,Data>::lo_bound[ jj ] < 10.0 )
	      simplex[ ii ][ jj ] = Opt<Func,Data>::lo_bound[ jj ] +
		( Opt<Func,Data>::hi_bound[ jj ] -
		  Opt<Func,Data>::lo_bound[ jj ] ) / 4.0;
	    else
	      simplex[ ii ][ jj ] = 
		std::max( Opt<Func,Data>::lo_bound[ jj ], model_par[ ii ] -
			  0.01 * fabs( model_par[ ii ] ) );
	  
	}

      int err_status = EXIT_SUCCESS;
      for ( int ii = 0; ii <= NPAR; ++ii ) {
	sherpa::Array2d<double>::Row rowii = simplex[ ii ];
	eval_user_func( maxnfev, &rowii[0],
			simplex[ ii ][ NPAR ], nfev,
			err_status );
	if ( EXIT_SUCCESS != err_status )
	  break;
      }

      return err_status;
      
    }                                                      // eval_init_simplex


    void init_simplex( int initsimplex, const double* model_par,
		       const double* step ) {

      int allzero = 0;
      for ( int ii = 0; ii < NPAR; ++ii )
	if ( 0.0 == step[ ii ] )
	  ++allzero;
      if ( NPAR == allzero )
	throw std::runtime_error( "Step cannot be all zeros (this would "
				  "yield an initial degenerate simplex)" );
      
      for ( int ii = 0; ii < NPAR; ++ii )
	simplex( 0, ii ) = model_par[ ii ];
      
      switch( initsimplex ) {
      case 1:
	SpendleyHextHimsworth_simplex( step, model_par );
	break;
      default:
	dtn_simplex( step, model_par );
      }

    }                                                           // init_simplex

  private:

    //
    // reflection_coef > 0, expansion_coef > 1,
    // expansion_coef > reflection_coef, 0 < contraction_coef < 1,
    // and 0 < shrinkage_coef < 1.
    //
    void check_coefficients( ) const {
      
      if ( reflection_coef <= 0.0 )
	throw std::runtime_error( "The reflection coefficient must be > 0" );
      if ( expansion_coef <= 1.0 )
	throw std::runtime_error( "The expansion coefficient must be > 1" );
      if ( contraction_coef <= 0.0 || contraction_coef >= 1.0 )
	throw std::runtime_error( "The contraction coefficient must be "
				  "within (0,1)" );

    }                                                     // check_coefficients


    void dtn_simplex( const double* step, const double* model_par ) {

      for ( int ii = 0; ii < NPAR; ++ii ) {
	for ( int jj = 0; jj < NPAR; ++jj )
	  simplex[ ii + 1 ][ jj ] = model_par[ jj ];
	simplex[ ii + 1 ][ ii ] = model_par[ ii ] + step[ ii ];
      }

    }                                                            // dtn_simplex


    //
    // construct a regular simplex, see: Spendley, Hext and Himsworth,
    // "Sequential Application of Simplex Designs in Optimization and
    // "Evolutionary Operation", Technometics, Vol 4, No 4, Nov 1962
    // pages 441-461.
    //
    void SpendleyHextHimsworth_simplex( const double* step,
					const double* model_par )  {

      double nparsqrt2 = NPAR * sqrt( 2.0 );
      double sqrtnpar1 = std::sqrt( double( NPAR + 1 ) );
      double pn = ( sqrtnpar1 - 1 + NPAR ) / nparsqrt2;
      double qn = ( sqrtnpar1 - 1 ) / nparsqrt2;
	
      for ( int ii = 1; ii <= NPAR; ++ii )
	for (int jj = 0; jj < NPAR; ++jj )
	  if ( ii - 1 == jj )
	    simplex[ ii ][ jj ] = pn + model_par[ jj ];
	  else
	    simplex[ ii ][ jj ] = qn + model_par[ jj ];

    }                                          // SpendleyHextHimsworth_simplex

  };               // class DirectSearch : public sherpa::OptFunc< Func, Data >

#define SIMPLEX DirectSearch<Func,Data>::simplex

}                                                          // namespace sherpa

#endif
