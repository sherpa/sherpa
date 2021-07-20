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


void tstminim( Init init, Fct fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
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
    for ( int jj = 0; jj < npar; ++jj )
      mypar[ jj ] = par[ jj ];
    sherpa::Bounds<double> bounds(lo, hi);

    sherpa::Minim<Fct, const sherpa::Bounds<double>&, double>* nm = NULL;
    
    if (reflect)
      nm = new sherpa::Minim<Fct, const sherpa::Bounds<double>&, double>( fct, bounds);
    else
      nm = new sherpa::MinimNoReflect<Fct, const sherpa::Bounds<double>&, double>( fct, bounds );

#if (__STDC_VERSION__ < 201112L)
    std::auto_ptr< sherpa::Minim<Fct, const sherpa::Bounds<double>&, double> > mynm(nm);
#else
    std::unique_ptr< sherpa::Minim<Fct, const sherpa::Bounds<double>&, double> > mynm(nm);
#endif

    int verbose=0, maxnfev=npar*npar*maxfev, nfev=0;
    double fmin;
    int initsimplex=1;
    mynm->operator()( verbose, maxnfev, tol, npar, initsimplex, finalsimplex, lo, hi,
      step, mypar, nfev, fmin );

    strcpy( header, "Minim_" );
    print_pars( header, fct_name, nfev, fmin, answer, npar, mypar );

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

(sherpa-4.14-dev) [dtn@devel12 src]$ time minim -r > noreflect

real0m0.160s
user0m0.158s
sys0m0.002s
(sherpa-4.14-dev) [dtn@devel12 src]$ time valgrind minim > reflect
==6304== Memcheck, a memory error detector
==6304== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==6304== Using Valgrind-3.14.0 and LibVEX; rerun with -h for copyright info
==6304== Command: minim
==6304== 
==6304== 
==6304== HEAP SUMMARY:
==6304==     in use at exit: 0 bytes in 0 blocks
==6304==   total heap usage: 93,345 allocs, 93,345 frees, 113,812,184 bytes allocated
==6304== 
==6304== All heap blocks were freed -- no leaks are possible
==6304== 
==6304== For counts of detected and suppressed errors, rerun with: -v
==6304== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

real0m1.876s
user0m1.820s
sys0m0.049s

(sherpa-4.14-dev) [dtn@devel12 src]$ time valgrind minim -r > noreflect
==6316== Memcheck, a memory error detector
==6316== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==6316== Using Valgrind-3.14.0 and LibVEX; rerun with -h for copyright info
==6316== Command: minim -r
==6316== 
==6316== 
==6316== HEAP SUMMARY:
==6316==     in use at exit: 0 bytes in 0 blocks
==6316==   total heap usage: 23,929 allocs, 23,929 frees, 3,681,464 bytes allocated
==6316== 
==6316== All heap blocks were freed -- no leaks are possible
==6316== 
==6316== For counts of detected and suppressed errors, rerun with: -v
==6316== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

real0m3.821s
user0m3.781s
sys0m0.040s

(sherpa-4.14-dev) [dtn@devel12 src]$ diff reflect noreflect 
11c11
< Minim_BrownBadlyScaled	2151	0	7.60681e-20	1e+06,1.99991e-06,1e+06,2.0001e-06,1e+06,2e-06
---
> Minim_BrownBadlyScaled	4547	0	1.10634e-17	1e+06,1.99977e-06,1e+06,1.99819e-06,1e+06,2e-06
41,44c41,44
< Minim_McCormick	93	-1.9132	-1.91322	-0.547198,-1.5472
< Minim_BoxBetts	94	0	8.5606e-16	1,10,1
< Minim_Paviani	1497	-45.778	-45.6242	9.32565,9.39821,9.34252,9.46265,9.35731,9.36868,9.28165,9.299,9.40574,9.21271
< Minim_GoldsteinPrice	92	3	3	7.848e-14,-1
---
> Minim_McCormick	4098	-1.9132	-nan	-nan,-nan
> Minim_BoxBetts	114	0	1.43046e-16	1,10,1
> Minim_Paviani	-102402	-45.778	-2.75503	5.78989,7.49795,5.06157,7.23616,4.91172,5.56942,7.55769,6.10409,6.90818,7.98112
> Minim_GoldsteinPrice	-132	3	84	1.8,0.2
49,51c49,51
< Minim_Levy5	-317	-11.504	38.8617	3.96386,3.99615,3.99623,3.99624,3.99945
< Minim_Levy6	-306	-11.504	47.8504	3.96386,3.99615,3.99623,3.99623,3.99624,3.99945
< Minim_Levy7	-383	-11.504	56.8392	3.96386,3.99615,3.99623,3.99623,3.99623,3.99624,3.99945
---
> Minim_Levy5	-263	-11.504	38.8617	3.96386,3.99615,3.99623,3.99624,3.99945
> Minim_Levy6	-345	-11.504	47.8504	3.96386,3.99615,3.99623,3.99623,3.99624,3.99945
> Minim_Levy7	-391	-11.504	56.8392	3.96386,3.99615,3.99623,3.99623,3.99623,3.99624,3.99945
53,54c53,54
< Minim_SixHumpCamel	91	-1.03	-1.03163	0.089842,-0.712656
< Minim_Branin	97	0.397889	0.397887	9.42478,2.475
---
> Minim_SixHumpCamel	93	-1.03	-1.03163	0.089842,-0.712656
> Minim_Branin	81	0.397889	0.397887	9.42478,2.475
67c67
< Minim_McKinnon	-4098	-0.25	-2.376e+11	-100,-100
---
> Minim_McKinnon	-4099	-0.25	-inf	-1.46426e+76,-1.06339e+77

*/
