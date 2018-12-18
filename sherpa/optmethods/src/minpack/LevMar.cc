#ifdef testLevMar

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

#include "LevMar.hh"

#include "tests/tstopt.hh"

#ifdef USE_ADEPT
typedef void (*FctVecAdept)(int m, int n, adept::adouble *x, adept::adouble *fvec, int& ierr, mybounds);
typedef void (*tstLMder)( Init init, FctVecAdept fct, int npar, std::vector<double>& par, std::vector<double>& lo, std::vector<double>& hi, double tol, const char* fct_name, int npop, int maxfev, double xprob, double sfactor );
#define ARG1 tstLMder
#define ARG2 adept::adouble
#define ARG3 FctVecAdept
#else
#define ARG1 tstFctVec
#define ARG2 double
#define ARG3 FctVec
#endif

template<typename aFunc>
void tstlm( Init init, aFunc fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
	    double tol, const char* fct_name, int npop, int maxfev,
	    double xprob, double sfactor ) {

  try {

    int mfcts;
    double answer;
    int maxnfev=128*npar, nfev=0;
    double fmin=0.0;

    init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
    int nprint = 0;
    double factor=100.0;
    std::vector<double> covarerr( 8 * npar, 0.0 );
    std::vector<double> fjac( mfcts * npar );
    sherpa::Bounds<double> bounds(lo, hi);

#ifdef USE_ADEPT
    minpack::LevMarDerAdept< aFunc, mybounds, double > lm( fct, bounds, mfcts, npar );
    int njev = 0;
    lm( npar, tol, tol, tol, maxnfev, factor, nprint, par, nfev, njev,
        fmin, bounds, fjac, covarerr );
    double mytol = 1.0e4*sqrt( std::numeric_limits<double>::epsilon());
    print_pars( "lmder_", fct_name, nfev, fmin, answer, npar, par,
         	&covarerr[0], mytol, njev );
#else
    minpack::LevMarDif< aFunc, mybounds, double > lm( fct, bounds, mfcts );
    double epsfcn = std::sqrt( std::numeric_limits< double >::epsilon( ) );
    lm( npar, tol, tol, tol, maxnfev, epsfcn, factor, nprint, par, nfev, fmin,
        bounds, fjac, covarerr );
    print_pars( "lmdif_", fct_name, nfev, fmin, answer, npar, par,
		&covarerr[0] );

#endif


  } catch( const sherpa::OptErr& oe ) {

    std::cerr << oe << '\n';

  }

}

int main( int argc, char* argv[] ) {

  int npar=16;
  if ( argc == 2 )
    npar = atoi( argv[1] );

  if ( npar % 2 || npar < 2 ) {
    std::cerr << "The minimum value for the free parameter must be an even "
      "and it is greater then 2\n";
    return EXIT_FAILURE;
  }

  double tol = std::sqrt( std::numeric_limits< double >::epsilon() );
  int maxfev=128, npop=0;
  double xprob=0.0, sfactor=0.0;
  std::cout << "#\n#:npar = " << npar << "\n";
  std::cout << "#:tol=" << tol << '\n';
  std::cout << "# A negative value for the nfev signifies that the "
    "optimization method did not converge\n#\n";
  std::cout << "name\tnfev\tnjev\tanswer\tfval\tpar\terr\nS\tN\tN\tN\tN\tN\n";
  tst_unc_opt<ARG1, ARG2>( npar, tol, tstlm<ARG3>, npop, maxfev, xprob,
                           sfactor);

  return 0;

}
#endif
/*

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/data/scialg/staff/dtn/soft/autodiff/adept-1.1/lib

g++ -I. -I../../../include -I../.. -I/data/scialg/staff/dtn/soft/autodiff/adept-1.1/include -Wall -ansi -pedantic -O3 -g -fopenmp -DtestLevMar -DUSE_ADEPT LevMar.cc -L/data/scialg/staff/dtn/soft/autodiff/adept-1.1/lib -lm -ladept -o tst_lmder

g++ -I. -I../../../include -I../.. -Wall -ansi -pedantic -O3 -g -DtestLevMar LevMar.cc -o tst_lmdif


 */
