#ifndef Opt_hh
#define Opt_hh

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


#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
#include <stdexcept>

namespace sherpa {

  class OptErr {

    friend std::ostream& operator << ( std::ostream& os, const OptErr& opte ) {
      return opte.print( os ); }

  public:

    enum Err { Success, Input, OutOfBound, MaxFev, UsrFunc, Unknown };

    OptErr( OptErr::Err e ) : err( e ) { }

    Err err;

  private:

    std::ostream& print( std::ostream& os ) const {

      const char* msg[] = {
	"No error",
	"Input error",
	"Parameter is out of bound",
	"Max number of function evaluation",
	"User Function error",
	"Unknown error"
      };

      os << msg[ err ];

      return os;

    }

  };

  template< typename real >
  class Bounds {

  public:

    Bounds( const std::vector<real>& l, const std::vector<real>& u ) :
      lb(l), ub(u) { }

    // bool are_pars_outside_limits( int npar, const std::vector<real>& par )
    //   const {
    //   for ( int ii = 0; ii < npar; ++ii )
    //     if ( par[ ii ] < lb[ ii ] || par[ ii ] > ub[ ii ] )
    //       return true;
    //   return false;
    // }

    const std::vector<real>& get_lb( ) const { return lb; }
    const std::vector<real>& get_ub( ) const { return ub; }

  private:

    const std::vector<real>& lb;
    const std::vector<real>& ub;

  };

  template< typename Data, typename real >
  class Opt {

  public:

    virtual ~Opt( ) { }

    Opt( Data data ) : usr_data( data ) { }

    //
    // If lo <= xpar <= hi return false, else return true
    //
    static bool are_pars_outside_limits( int npar,
                                         const std::vector<real>& par,
                                         const sherpa::Bounds<real>& limits ) {
      const std::vector<real>& low = limits.get_lb();
      const std::vector<real>& high = limits.get_ub();
      for ( int ii = 0; ii < npar; ++ii )
	if ( par[ ii ] < low[ ii ] || par[ ii ] > high[ ii ] )
	  return true;
      return false;
    }

    virtual real eval_func( int maxnfev, const sherpa::Bounds<real>& limits,
                            int npar, std::vector<real>& par, int& nfev ) {
      std::cerr << "Opt::eval_func define me!\n";
      return 0.0; }

    virtual int minimize( int maxnfev, const sherpa::Bounds<real>& limits,
			  double tol, int npar, std::vector<real>& par,
			  double& fmin, int& nfev ) {
      std::cerr << "Opt::minimize define me!\n";
      return EXIT_FAILURE;
    }


    static std::ostream& print_par( std::ostream& os,
                                    const std::vector<real>& par ) {
      const size_t npar = par.size( ) - 1;
      os.precision( 6 );
      os << "f( " << std::scientific << par[ 0 ];
      for ( size_t ii = 1; ii < npar; ++ii )
	os << ", " << std::scientific << par[ ii ];
      os << " ) = " << par[ npar ] << '\n';
      return os;
    }

  protected:

    Data get_usr_data( ) { return usr_data; }


  private:

    Data usr_data;

  };                                                               // class Opt

  template< typename Func, typename Data, typename real >
  class OptFunc : public Opt<Data, real> {

  public:

    virtual ~OptFunc( ) { }

    OptFunc( Func func, Data data, int mfct=0 ) : Opt< Data, real>( data ),
                                                  usr_func( func ),
                                                  mfcts( mfct ) { }


    virtual real eval_func( int maxnfev, const sherpa::Bounds<real>& limits,
                            int npar, std::vector<real>& par, int& nfev ) {

      if ( sherpa::Opt<Data, real>::are_pars_outside_limits( npar, par,
                                                             limits ) ) {
	par[ npar ] = std::numeric_limits< real >::max( );
	return par[ npar ];
      }

      ++nfev;

      int ierr = EXIT_SUCCESS;
      usr_func( npar, &par[0], par[npar], ierr, sherpa::Opt<Data, real>::get_usr_data() );
      if ( EXIT_SUCCESS != ierr )
	throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
      if ( nfev >= maxnfev )
	throw sherpa::OptErr( sherpa::OptErr::MaxFev );

      if ( nfev >= maxnfev )
	ierr = sherpa::OptErr::MaxFev;

      return par[ npar ];

    }                                                         // eval_user_func

    int minimize( int maxnfev, const sherpa::Bounds<real>& limits,
		  real tol, int npar, std::vector<real>& par, real& fmin,
		  int& nfev ) {
      int ierr = EXIT_SUCCESS;
      fmin = eval_func( maxnfev, limits, npar, par, nfev );
      return ierr;
    }

  protected:

    Func get_func( ) { return usr_func; }

  private:

    Func usr_func;
    const int mfcts;

    OptFunc& operator = (OptFunc const&); // declare but, purposely, not define
    OptFunc( OptFunc const& );            // declare but, purposely, not define

  };                                                           // class OptFunc

}                                                           // namespace sherpa

#endif                                                        // #ifndef Opt_hh

#ifdef testOpt
#include "sherpa/functor.hh"
#include "tests/tstoptfct.hh"

template< typename Init, typename Func >
void justdoit( Init init, Func fct, int npar, std::vector<double>& par ) {

  int mfcts, maxnfev = 2, nfev = 0;
  double answer, fmin;
  std::vector< double > lo( npar ), hi( npar );
  init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );

  sherpa::Bounds<double> bounds(lo, hi);
  sherpa::OptFunc< Func, const sherpa::Bounds<double>&, double>
    optfunc( fct, bounds );

  try {

    for ( int ii = 0; ii < maxnfev + 1; ++ii )
      fmin = optfunc.eval_func( maxnfev, npar, par, nfev );

  } catch( sherpa::OptErr& oe ) {

    std::cerr << oe << '\n';
    std::cerr << "nfev = " << nfev << '\n';

  }

}

int main( int argc, char* argv[] ) {

  int npar = 2;
  if ( argc == 2 )
    npar = atoi( argv[1] );

  std::vector<double> par( npar + 1, 0.0 );
  sherpa::FctPtr< void, int, double*, double&, int&,
    const sherpa::Bounds<double>& >
    fct( tstoptfct::Rosenbrock<double,const sherpa::Bounds<double>&> );

  justdoit( sherpa::fct_ptr( tstoptfct::RosenbrockInit<double> ),
	    fct, npar, par );

  return 0;

}

#endif                                                        // #ifdef testOpt
//
// cp Opt.hh tmp.cc; g++ -I.. -I../../include/ -g -Wall -ansi -pedantic -O3 -DtestOpt tmp.cc; rm tmp.cc; myvalgrind a.out
//
