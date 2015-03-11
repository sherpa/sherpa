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

  class Opt {

  public:

    typedef std::vector<double> myvec;
    typedef std::pair< const myvec, const myvec > mypair;

    virtual ~Opt( ) { }
    
    Opt( ) { }

    virtual double eval_func( int maxnfev, const Opt::mypair& limits, int npar,
			      Opt::myvec& par, int& nfev ) {
      std::cerr << "Opt::eval_func define me!\n";
      return 0.0; }

    virtual int minimize( int maxnfev, const sherpa::Opt::mypair& limits,
			  double tol, int npar, sherpa::Opt::myvec& par,
			  double& fmin, int& nfev ) {
      std::cerr << "Opt::minimize define me!\n";
      return EXIT_FAILURE;
    }

    //
    // If lo <= xpar <= hi return false, else return true
    //
    static bool are_pars_outside_limits( int npar, const mypair& limits,
					 const myvec& par ) {
      const myvec& low = limits.first;
      const myvec& high = limits.second;
      for ( int ii = 0; ii < npar; ++ii )
	if ( par[ ii ] < low[ ii ] || par[ ii ] > high[ ii ] )
	  return true;
      return false;

    }

    static std::ostream& print_par( std::ostream& os, const myvec& par ) {
      const size_t npar = par.size( ) - 1;
      os.precision( 6 );
      os << "f( " << std::scientific << par[ 0 ];
      for ( size_t ii = 1; ii < npar; ++ii )
	os << ", " << std::scientific << par[ ii ];
      os << " ) = " << par[ npar ] << '\n';
      return os;
    }

  private:

  };                                                               // class Opt

  template< typename Func, typename Data >
  class OptFunc : public Opt {

  public:

    virtual ~OptFunc( ) { }

    OptFunc( Func func, Data data, int mfct=0 ) :
      usr_func( func ), usr_data( data ), mfcts( mfct ) { }


    virtual double eval_func( int maxnfev, const Opt::mypair& limits, int npar,
			      Opt::myvec& par, int& nfev ) {
      
      if ( sherpa::Opt::are_pars_outside_limits( npar, limits, par ) ) {
	par[ npar ] = std::numeric_limits< double >::max( );
	return par[ npar ];
      }

      ++nfev;
      
      int ierr = EXIT_SUCCESS;
      usr_func( npar, &par[0], par[npar], ierr, usr_data );
      if ( EXIT_SUCCESS != ierr )
	throw sherpa::OptErr( sherpa::OptErr::UsrFunc );
      if ( nfev >= maxnfev )
	throw sherpa::OptErr( sherpa::OptErr::MaxFev );

      if ( nfev >= maxnfev )
	ierr = sherpa::OptErr::MaxFev;

      return par[ npar ];

    }                                                         // eval_user_func

    int minimize( int maxnfev, const sherpa::Opt::mypair& limits,
		  double tol, int npar, sherpa::Opt::myvec& par, double& fmin,
		  int& nfev ) {
      int ierr = EXIT_SUCCESS;
      fmin = eval_func( maxnfev, limits, npar, par, nfev );
      return ierr;
    }

  protected:

    Func get_func( ) { return usr_func; }
    Data get_data( ) { return usr_data; }

  private:

    Func usr_func;
    Data usr_data;
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

  const sherpa::Opt::mypair limits( lo, hi );

  sherpa::OptFunc< Func, void* > optfunc( fct, NULL );

  try {

    for ( int ii = 0; ii < maxnfev + 1; ++ii )
      fmin = optfunc.eval_func( maxnfev, limits, npar, par, nfev );

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
  sherpa::FctPtr< void, int, double*, double&, int&, void* >
    fct( tstoptfct::Rosenbrock<double,void*> );

  justdoit( sherpa::fct_ptr( tstoptfct::RosenbrockInit<double> ),
	    fct, npar, par );

  return 0;

}

#endif                                                        // #ifdef testOpt
//
// cp Opt.hh tmp.cc; g++ -I.. -I../../include/ -g -Wall -ansi -pedantic -O3 -DtestOpt tmp.cc; rm tmp.cc; myvalgrind a.out
//
