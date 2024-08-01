#ifdef testMinim

//
//  Copyright (C) 2020, 2021  Smithsonian Astrophysical Observatory
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

#include <cstring>
#include <memory>

#include "minim.hh"
#include "tests/tstopt.hh"


void tstminim( Init init, Fct fct, int npar, sherpa::Array1D<double>& par,
	    sherpa::Array1D<double>& lo, sherpa::Array1D<double>& hi,
	    double tol, const char* fct_name, int reflect, int maxfev,
	    double c1, double c2 ) {

  try {

    char header[64];

    //
    // you may think you are clever by eliminating the following overhead
    // and simply use the vector par, but believe me it this is necessary
    //
    std::vector<double> mypar( npar, 0.0 );

    std::vector<double> step( npar * npar * 4 );
    for ( int ii = 0; ii < npar; ++ii )
      step[ ii ] = 1.2;

    std::vector< int > finalsimplex(npar, 0);

    int mfcts;
    double answer;

    init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
    std::copy( &par[0], &par[npar], &mypar[0] );
    sherpa::Bounds<double> bounds(lo, hi);

    sherpa::Minim<Fct, const sherpa::Bounds<double>&, double>* nm = NULL;
    
    if (reflect)
      nm = new sherpa::Minim<Fct, const sherpa::Bounds<double>&, double>( fct, bounds);
    else
      nm = new sherpa::MinimNoReflect<Fct, const sherpa::Bounds<double>&, double>( fct, bounds );

#if (__cplusplus < 201103L)
    std::auto_ptr< sherpa::Minim<Fct, const sherpa::Bounds<double>&, double> > mynm(nm);
#else
    std::unique_ptr< sherpa::Minim<Fct, const sherpa::Bounds<double>&, double> > mynm(nm);
#endif

    int verbose=0, maxnfev=npar*npar*maxfev, nfev=0;
    double fmin;
    int initsimplex=1;
    mynm->operator()( verbose, maxnfev, tol, npar, initsimplex, finalsimplex, bounds,
      step, mypar, nfev, fmin );

    strcpy( header, "Minim_" );
    std::copy( &mypar[0], &mypar[npar], &par[0] );
    print_pars( header, fct_name, nfev, fmin, answer, npar, par );

  } catch( const sherpa::OptErr& oe ) {

    std::cerr << oe << '\n';

  }

  return;

}

int main( int argc, char* argv[] ) {

  try {

    int c, uncopt = 1, globalopt = 1, reflect = 1;
    while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
      while ( (c = *++argv[ 0 ]) )
	switch( c ) {
	case 'u':
	  uncopt = 0;
	  break;
	case 'g':
	  globalopt = 0;
	  break;
        case 'r':
          reflect = 0;
          break;
	default:
	  fprintf( stderr, "%s: illegal option '%c'\n", argv[ 0 ], c );
	  fprintf( stderr, "Usage %s [ -g ] [ -u ] [ npar ]\n", argv[ 0 ] );
	  return EXIT_FAILURE;
      }


    int npar=6;
    if ( argc == 1 )
      npar = atoi( *argv );

    if ( npar % 2 || npar < 2 ) {
      printf( "The minimum value for the free parameter must be an even "
	      "and it is greater then 2\n" );
      return EXIT_FAILURE;
    }

    double tol = 1.0e-8;
    std::cout << "#\n#:npar = " << npar << "\n";
    std::cout << "#:tol=" << tol << '\n';
    std::cout << "# A negative value for the nfev signifies that the "
      "optimization method did not converge\n#\n";
    std::cout << "name\tnfev\tanswer\tstat\tpar\nS\tN\tN\tN\tN\n";

    int maxfev=1024;
    double c1=0.0, c2=0.0;
    if ( uncopt )
      tst_unc_opt<tstFct, double>( npar, tol, tstminim, reflect, maxfev, c1, c2 );

    if ( globalopt )
      tst_global( npar, tol, tstminim, reflect, maxfev, c1, c2 );

    return EXIT_SUCCESS;

  } catch( std::exception& e ) {

    std::cerr << e.what( ) << '\n';
    return EXIT_FAILURE;

  }

  return 0;
}
#endif

/*
g++ -o minim -DtestMinim -DNDEBUG -Wall -ansi -pedantic -O3 -I../../include -I.. minim.cc

(base) [dtn@devel12 src]$ time minim > minim.reflect

real0m0.081s
user0m0.079s
sys0m0.002s
(base) [dtn@devel12 src]$ time minim -r > minim.noreflect

real0m0.162s
user0m0.160s
sys0m0.002s
(base) [dtn@devel12 src]$ diff minim.reflect minim.noreflect 
11c11
< Minim_BrownBadlyScaled215107.60681e-201e+06,1.99991e-06,1e+06,2.0001e-06,1e+06,2e-06
---
> Minim_BrownBadlyScaled454701.10634e-171e+06,1.99977e-06,1e+06,1.99819e-06,1e+06,2e-06
41,44c41,44
< Minim_McCormick93-1.9132-1.91322-0.547198,-1.5472
< Minim_BoxBetts9408.5606e-161,10,1
< Minim_Paviani1497-45.778-45.62429.32565,9.39821,9.34252,9.46265,9.35731,9.36868,9.28165,9.299,9.40574,9.21271
< Minim_GoldsteinPrice92337.848e-14,-1
---
> Minim_McCormick4098-1.9132-nan-nan,-nan
> Minim_BoxBetts9901.02735e-141,10,1
> Minim_Paviani-102402-45.778-2.755035.78989,7.49795,5.06157,7.23616,4.91172,5.56942,7.55769,6.10409,6.90818,7.98112
> Minim_GoldsteinPrice8633-2.85576e-11,-1
53,54c53,54
< Minim_SixHumpCamel91-1.03-1.031630.089842,-0.712656
< Minim_Branin970.3978890.3978879.42478,2.475
---
> Minim_SixHumpCamel89-1.03-1.031630.089842,-0.712656
> Minim_Branin810.3978890.3978879.42478,2.475
67c67
< Minim_McKinnon-4098-0.25-2.376e+11-100,-100
---
> Minim_McKinnon-4099-0.25-inf-1.46426e+76,-1.06339e+77

(base) [dtn@devel12 src]$ valgrind minim -r > minim.noreflect
==10907== Memcheck, a memory error detector
==10907== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==10907== Using Valgrind-3.14.0 and LibVEX; rerun with -h for copyright info
==10907== Command: minim -r
==10907== 
==10907== 
==10907== HEAP SUMMARY:
==10907==     in use at exit: 0 bytes in 0 blocks
==10907==   total heap usage: 23,929 allocs, 23,929 frees, 3,684,216 bytes allocated
==10907== 
==10907== All heap blocks were freed -- no leaks are possible
==10907== 
==10907== For counts of detected and suppressed errors, rerun with: -v
==10907== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
(base) [dtn@devel12 src]$ valgrind  minim > minim.reflect
==10935== Memcheck, a memory error detector
==10935== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==10935== Using Valgrind-3.14.0 and LibVEX; rerun with -h for copyright info
==10935== Command: minim
==10935== 
==10935== 
==10935== HEAP SUMMARY:
==10935==     in use at exit: 0 bytes in 0 blocks
==10935==   total heap usage: 93,345 allocs, 93,345 frees, 113,814,936 bytes allocated
==10935== 
==10935== All heap blocks were freed -- no leaks are possible
==10935== 
==10935== For counts of detected and suppressed errors, rerun with: -v
==10935== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

*/
