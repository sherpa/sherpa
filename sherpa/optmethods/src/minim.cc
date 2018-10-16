#ifdef testMinim

#include <cstring>

#include "minim.hh"
#include "tests/tstopt.hh"


void tst_minim( ) {

  const int npar=2;
  std::vector<int> finalsimplex(npar);
  std::vector<double> x(npar), vc(npar*(npar+1)/2, 0.0);
  std::vector<double> lb(npar, -10.0), ub(npar, 10.0), step(npar, 1.0);
  x[0] = -1.0;
  x[1] = 1.2;

  int maxfev=1024, iprint=0, neval=0, initsimplex=1;
  double func=0.0, tol=1.0e-6;

  const sherpa::Bounds<double> bounds(lb, ub);
  sherpa::Minim<Fct, const sherpa::Bounds<double>&, double>
    minim( tstoptfct::Rosenbrock, bounds );
  minim( iprint, maxfev, tol, npar, initsimplex, finalsimplex, lb, ub,
         step, x, neval, func );

  // std::cout << "ifault = " << ifault << ", neval = " << neval << '\n';
  // std::cout << "f(" << x[0];
  // for ( int ii = 1; ii < npar; ++ii )
  //   std::cout << ", " << x[ii];
  // std::cout << ") = " << func << "\nvc = ";;
  // for ( int ii = 0; ii < npar*(npar+1)/2; ++ii )
  //   std::cout << vc[ii] << '\t';
  // std::cout << '\n';

}

void tstminim( Init init, Fct fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
	    double tol, const char* fct_name, int npop, int maxfev,
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
    sherpa::Minim< Fct, const sherpa::Bounds<double>&, double > nm( fct,
                                                                    bounds );

    int verbose=0, maxnfev=npar*npar*maxfev, nfev;
    double fmin;
    int initsimplex=1;
    nm( verbose, maxnfev, tol, npar, initsimplex, finalsimplex, lo, hi,
        step, mypar, nfev, fmin );

    strcpy( header, "Minim_" );
    print_pars( header, fct_name, nfev, fmin, answer, npar, mypar );

  } catch( const sherpa::OptErr& oe ) {

    std::cerr << oe << '\n';

  }

  return;

}

int main( int argc, char* argv[] ) {

  tst_minim();

  try {

    int c, uncopt = 1, globalopt = 1;
    while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
      while ( (c = *++argv[ 0 ]) )
	switch( c ) {
	case 'u':
	  uncopt = 0;
	  break;
	case 'g':
	  globalopt = 0;
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

    int npop=0, maxfev=1024;
    double c1=0.0, c2=0.0;
    if ( uncopt )
      tst_unc_opt( npar, tol, tstminim, npop, maxfev, c1, c2 );

    if ( globalopt )
      tst_global( npar, tol, tstminim, npop, maxfev, c1, c2 );

    return EXIT_SUCCESS;

  } catch( std::exception& e ) {

    std::cerr << e.what( ) << '\n';
    return EXIT_FAILURE;

  }

  return 0;
}
#endif
//
// g++ -o minim -DtestMinim -Wall -ansi -pedantic -O3 -I../../include -I.. minim.cc
//
/*
==9066== Memcheck, a memory error detector
==9066== Copyright (C) 2002-2012, and GNU GPL'd, by Julian Seward et al.
==9066== Using Valgrind-3.8.1 and LibVEX; rerun with -h for copyright info
==9066== Command: minim
==9066==
#
#:npar = 6
#:tol=1e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	stat	par
S	N	N	N	N
Minim_Rosenbrock	-961	0	1.29567	0.184248,0.0284826,1.63673,2.68041,0.529682,0.278636
Minim_FreudensteinRoth	-574	0	146.953	11.4128,-0.896805,11.4128,-0.896805,11.4128,-0.896805
Minim_PowellBadlyScaled	1003	0	1.63198e-06	1.50587e-05,6.64168,8.25696e-06,12.1085,1.30173e-05,7.68102
Minim_BrownBadlyScaled	4547	0	1.10634e-17	1e+06,1.99977e-06,1e+06,1.99819e-06,1e+06,2e-06
Minim_Beale	399	0	2.70425e-16	3,0.5,3,0.5,3,0.5
Minim_JennrichSampson	861	373.086	373.087	0.257825,0.257825,0.257825,0.257825,0.257825,0.257825
Minim_HelicalValley	145	0	1.12431e-16	1,4.47078e-09,6.95162e-09
Minim_Bard	154	0.00821487	0.00821488	0.0824106,1.13304,2.34369
Minim_Gaussian	112	1.12793e-08	1.12793e-08	0.398956,1.00002,-4.65898e-10
Minim_Meyer	3081	87.9458	87.9459	0.00560964,6181.35,345.224
Minim_GulfResearchDevelopment	308	0	7.68697e-05	19.4537,-6.56898,1.03413
Minim_Box3d	214	0	8.7639e-16	1,10,1
Minim_PowellSingular	243	0	1.16082e-08	-0.00181298,0.000181175,0.00351953,0.00353512
Minim_Wood	316	0	2.6502e-15	1,1,1,1
Minim_KowalikOsborne	255	0.000307505	0.000307506	0.192807,0.19129,0.123059,0.136066
Minim_BrownDennis	332	85822.2	85822.2	-11.5944,13.2036,-0.403439,0.236779
Minim_Osborne1	519	5.46489e-05	5.49625e-05	0.376183,2.03372,-1.56315,0.0130563,0.0217551
Minim_Biggs	162	0	4.61328e-10	1.60027,10.0005,1.00001,4.99958,3.99978,2.99958
Minim_Osborne2	-3415	0.0401377	0.0772305	1.22656,0.322674,0.568139,0.537683,0.532828,1.67286,1.32382,4.26204,2.426,4.58298,5.66337
Minim_Watson	576	0.00228767	0.00228767	-0.0157251,1.01243,-0.232992,1.26043,-1.51373,0.992996
Minim_PenaltyI	-274	9.37629e-06	4.50241e-05	-0.198834,-0.221917,0.390035,-0.0954568
Minim_PenaltyII	-153	9.37629e-06	9.98143e-06	0.200004,0.400644,0.422337,-0.0415098
Minim_VariablyDimensioned	376	0	1.61504e-24	1,1,1,1,1,1
Minim_Trigonometric	600	0	5.22623e-14	0.0121122,0.0117014,0.0113404,-0.122016,-0.0940155,0.0104681
Minim_BrownAlmostLinear	-362	1	1.19581e-14	1,1,1,1,1,1
Minim_DiscreteBoundary	281	0	1.75795e-18	-0.0898882,-0.167863,-0.15361,-0.132462,-0.102058,-0.0592968
Minim_DiscreteIntegral	283	0	7.94022e-19	-0.0654635,-0.118166,-0.154627,-0.169992,-0.15727,-0.106031
Minim_BroydenTridiagonal	346	0	3.40409e-17	-0.576058,-0.69593,-0.680249,-0.642988,-0.556421,-0.366025
Minim_BroydenBanded	449	0	3.44158e-16	-0.428303,-0.476596,-0.519653,-0.558073,-0.593437,-0.593437
Minim_LinearFullRank	326	0	1.77494e-30	-1,-1,-1,-1,-1,-1
Minim_LinearFullRank1	247	1.15385	1.15385	2.97377,2.2046,2.13584,0.15931,0.59993,-2.8661
Minim_LinearFullRank0cols0rows	244	2.66667	2.66667	2.14801,1.95638,0.261166,0.8645,-1.56419,2.61241
Minim_Chebyquad	641	0	8.64238e-06	0.04472,0.205762,0.230925,0.424449,0.491337,0.591077,0.761322,0.804917,0.956254
Minim_McCormick	4098	-1.91	-nan	-nan,-nan
Minim_BoxBetts	99	0	1.02735e-14	1,10,1
Minim_Paviani	-102402	-45.7	-2.75503	5.78989,7.49795,5.06157,7.23616,4.91172,5.56942,7.55769,6.10409,6.90818,7.98112
Minim_GoldsteinPrice	86	3	3	-2.85576e-11,-1
Minim_Shekel5	-171	-10.1532	-5.10077	7.99958,7.99964,7.99958,7.99964
Minim_Shekel7	-179	-10.4029	-5.12882	7.99951,7.99962,7.9995,7.99961
Minim_Shekel10	-189	-10.5364	-5.17565	7.99948,7.99945,7.99946,7.99944
Minim_Levy4	-803	-21.502	9.85125	2.64769,1.99586,0.670413,6.99797
Minim_Levy5	-317	-11.504	38.8617	3.96386,3.99615,3.99623,3.99624,3.99945
Minim_Levy6	-306	-11.504	47.8504	3.96386,3.99615,3.99623,3.99623,3.99624,3.99945
Minim_Levy7	-383	-11.504	56.8392	3.96386,3.99615,3.99623,3.99623,3.99623,3.99624,3.99945
Minim_Griewank	-79	0	4.91141	100.481,-97.6456
Minim_SixHumpCamel	89	-1.03	-1.03163	0.089842,-0.712656
Minim_Branin	81	0.397889	0.397887	9.42478,2.475
Minim_Shubert	-96	-24.06	-14.6909	8.82716,5.79179
Minim_Hansen	113	-176.54	-176.542	4.97648,4.85806
Minim_Cola	-866	12.8154	256.203	1.94861,-0.57903,0.373911,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.673911,-0.0382,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108,0.0739108
Minim_Ackley	-85	0	19.3325	16.9988,16.9988
Minim_Bohachevsky1	122	0	0	-1.20688e-13,-4.84744e-14
Minim_Bohachevsky2	122	0	0	4.72187e-13,-1.38867e-13
Minim_Bohachevsky3	119	0	0	-7.20017e-14,5.35598e-14
Minim_Easom	-4096	-1	-0	25.55,25.55
Minim_Rastrigin	-97	0	7.95966	1.98991,1.98991
Minim_Michalewicz2	77	-1.8013	-1.8013	2.20291,1.5708
Minim_Michalewicz5	-323	-4.68766	-4.3749	2.20291,1.5708,2.21933,1.92306,0.996677
Minim_Michalewicz10	-938	-9.66015	-7.54148	2.20303,1.57081,1.28501,1.39237,1.72048,1.57081,2.22106,1.5867,1.65572,1.5708

==9066==
==9066== HEAP SUMMARY:
==9066==     in use at exit: 0 bytes in 0 blocks
==9066==   total heap usage: 24,021 allocs, 24,021 frees, 3,682,352 bytes allocated
==9066==
==9066== All heap blocks were freed -- no leaks are possible
==9066==
==9066== For counts of detected and suppressed errors, rerun with: -v
==9066== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 6 from 6
 */
