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

void tstlm( Init init, FctVec fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
	    double tol, const char* fct_name, int npop, int maxfev,
	    double xprob, double sfactor ) {

  try {

    int mfcts;
    double answer;

    init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
    sherpa::Bounds<double> bounds(lo, hi);
    minpack::LevMarDif< FctVec, mybounds, double > lm( fct, bounds, mfcts );

    int maxnfev=128*npar, nfev;
    double fmin=0.0;

    int nprint = 0;
    double epsfcn = std::sqrt( std::numeric_limits< double >::epsilon( ) );
    double factor=100.0;
    std::vector<double> covarerr( 8 * npar );
    std::vector<double> fjac( mfcts * npar );

    lm( npar, tol, tol, tol, maxnfev, epsfcn, factor, nprint, par, nfev, fmin,
        bounds, fjac, covarerr );

    // lm.minimize( &par[0], tol, maxnfev, nfev, fmin );

    print_pars( "lmdif_", fct_name, nfev, fmin, answer, npar, par,
		&covarerr[0] );

  } catch( const sherpa::OptErr& oe ) {

    std::cerr << oe << '\n';

  }

}

typedef struct  {
    int m;
    double *y;
    double *xmin;
    double *xmax;
} fcndata_t;

typedef int (*fcnLMDER)(int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag, fcndata_t*);


int fcn(int m, int n, const double *x, double *fvec, double *fjac,
        int ldfjac, int iflag, fcndata_t* p) {

  /*      subroutine fcn for lmder example. */

  int i;
  double tmp1, tmp2, tmp3, tmp4;
  const double *y = ((fcndata_t*)p)->y;
#ifdef BOX_CONSTRAINTS
  const double *xmin = ((fcndata_t*)p)->xmin;
  const double *xmax = ((fcndata_t*)p)->xmax;
  int j;
  double xb[3];
  double jacfac[3];
  double xmiddle, xwidth, th;

  for (j = 0; j < 3; ++j) {
    xmiddle = (xmin[j]+xmax[j])/2.;
    xwidth = (xmax[j]-xmin[j])/2.;
    th =  tanh((x[j]-xmiddle)/xwidth);
    xb[j] = (xmin[j]+xmax[j])/2. + th * xwidth;
    jacfac[j] = 1. - th * th;
  }
  x = xb;
#endif

  // assert(m == 15 && n == 3);

  if (iflag == 0) {
    /*      insert print statements here when nprint is positive. */
    /* if the nprint parameter to lmder is positive, the function is
       called every nprint iterations with iflag=0, so that the
       function may perform special operations, such as printing
       residuals. */
    return 0;
  }

  if (iflag != 2) {
    /* compute residuals */
    for (i=0; i < 15; ++i) {
      tmp1 = i + 1;
      tmp2 = 15 - i;
      tmp3 = (i > 7) ? tmp2 : tmp1;
      fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
    }
  } else {
    /* compute Jacobian */
    for (i=0; i<15; ++i) {
      tmp1 = i + 1;
      tmp2 = 15 - i;
      tmp3 = (i > 7) ? tmp2 : tmp1;
      tmp4 = (x[1]*tmp2 + x[2]*tmp3); tmp4 = tmp4*tmp4;
      fjac[i + ldfjac*0] = -1.;
      fjac[i + ldfjac*1] = tmp1*tmp2/tmp4;
      fjac[i + ldfjac*2] = tmp1*tmp3/tmp4;
    }
#    ifdef BOX_CONSTRAINTS
    for (j = 0; j < 3; ++j) {
      for (i=0; i < 15; ++i) {
        fjac[i + ldfjac*j] *= jacfac[j];
      }
    }
#    endif
  }
  return 0;
}


void tstlmder( ) {
  const int m = 15;
  const int n = 3;

  double xx[n] = {1.0, 1.0, 1.0};
  double y[m] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
                  3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

  double xmin[n] = {0., 0.1, 0.5};
  double xmax[n] = {2., 1.5, 2.3};
  std::vector<double> low(n);
  std::vector<double> high(n);
  std::vector<double> x(n);
  std::vector<double> fjac(m*n);
  for ( int ii = 0; ii < n; ++ii ) {
    low[ii] = xmin[ii];
    high[ii] = xmax[ii];
    x[ii] = xx[ii];
  }

  fcndata_t data;
  data.m = m;
  data.y = y;

  data.xmin = xmin;
  data.xmax = xmax;

  double ftol = std::sqrt( std::numeric_limits<double>::epsilon( ) );
  double xtol = std::sqrt( std::numeric_limits<double>::epsilon( ) );
  double epsfcn = ftol;
  double factor = 100.0;
  double gtol = 0.;
  double fmin;
  int nprint = 0;
  int nfev=0, njev=0, maxfev=400, rank;

  const sherpa::Bounds<double> bounds( low, high );

  minpack::LevMarDer< fcnLMDER, fcndata_t*, double > lm( fcn, &data, m );
  lm.fitme( n, ftol, xtol, gtol, maxfev, epsfcn, factor, nprint, x, nfev, njev,
            fmin, fjac, bounds, rank );

  double answer[n] = {0.0836609, 1.17808, 2.3 };

  for ( int ii = 0; ii < n; ++ii )
    if ( 0 != sao_fcmp( x[ii], answer[ii], 1.0e-5 ) )
      std::cerr << "answer (" << answer[ii ] << ") != x (" << x[ii] << ")\n";

  double myfjac[n][n] =
    {
      { 0.000153669, 0.00299417, -0.00278012 },
      { 0.00299417, 0.102588, -0.0986002 },
      { -0.00278012, -0.0986002, 0.0952316 }
    };

  for (int i=0; i<n; ++i)
    for (int j=0; j<n; ++j)
      if ( 0 != sao_fcmp( fjac[ i  * m + j ], myfjac[i][j], 1.0e-5 ) )
        std::cerr << "fjac " << fjac[ i  * m + j ] << " !=  myfjac " << myfjac[i][j] << '\n';

}

// void tstlmder_bounds( ) {
//   const int m = 15;
//   const int n = 3;

//   std::vector<double> fjac( m * n );
//   double xx[n] = {1.0, 1.0, 1.0};
//   double y[m] = {1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1,
//                   3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34, 2.1, 4.39};

//   double xmin[n] = {0., 0.1, 0.5};
//   double xmax[n] = {2., 1.5, 2.3};
//   std::vector<double> low(n);
//   std::vector<double> high(n);
//   std::vector<double> x(n);
//   for ( int ii = 0; ii < n; ++ii ) {
//     low[ii] = xmin[ii];
//     high[ii] = xmax[ii];
//     x[ii] = xx[ii];
//   }

//   fcndata_t data;
//   data.m = m;
//   data.y = y;

//   data.xmin = xmin;
//   data.xmax = xmax;

//   double ftol = std::sqrt( std::numeric_limits<double>::epsilon( ) );
//   double xtol = std::sqrt( std::numeric_limits<double>::epsilon( ) );
//   double epsfcn = ftol;
//   double factor = 100.0;
//   double gtol = 0.;
//   double fmin;
//   int nprint = 0;
//   int nfev=0, njev=0, maxfev=400, rank;

//   const sherpa::Bounds<double> bounds( low, high );

//   minpack::LevMarDer< fcnLMDER, fcndata_t*, double > lm( fcn, &data, m, n );
//   lm.fitme( n, ftol, xtol, gtol, maxfev, epsfcn, factor, nprint, x, nfev, njev,
//             fmin, bounds, rank );

//   double answer[n] = {0.08241058, 1.133037, 2.343695};
//   for ( int ii = 0; ii < n; ++ii )
//     if ( 0 != sao_fcmp( x[ii], answer[ii], 1.0e-6 ) )
//       std::cerr << "answer (" << answer[ii ] << ") != x (" << x[ii] << ")\n";

//   double myfjac[n][n] =
//     {
//       { 0.0001531202, 0.002869941, -0.002656662 },
//       {0.002869941, 0.09480935, -0.09098995},
//       {-0.002656662, -0.09098995, 0.08778727}
//     };

//   for (int i=0; i<n; ++i)
//     for (int j=0; j<n; ++j)
//       if ( 0 != sao_fcmp( fjac[ i  * m + j ], myfjac[i][j], 1.0e-6 ) )
//         std::cerr << "fjac " << fjac[ i  * m + j ] << " !=  myfjac " << myfjac[i][j] << '\n';

// }



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
  std::cout << "name\tnfev\tanswer\tfval\tpar\terr\nS\tN\tN\tN\tN\tN\n";
  tst_unc_opt( npar, tol, tstlm, npop, maxfev, xprob, sfactor );

  tstlmder();

  return 0;

}
#endif

/*
g++ -ansi -pedantic -Wall -O3 -I. -I../../../include -I../.. -DtestLevMar.cc LevMar.cc -o levmar

==23858== Memcheck, a memory error detector
==23858== Copyright (C) 2002-2012, and GNU GPL'd, by Julian Seward et al.
==23858== Using Valgrind-3.8.1 and LibVEX; rerun with -h for copyright info
==23858== Command: a.out
==23858==
#
#:npar = 16
#:tol=1.49012e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	fval	par	err
S	N	N	N	N	N
lmdif_Rosenbrock	278	0	0	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	1,2.00262,1,2.00262,1,2.00262,1,2.00262,1,2.00262,1,2.00262,1,2.00262,1,2.00262
lmdif_FreudensteinRoth	-142	0	391.874	11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849,11.4122,-0.896849	22488.5,1680.69,22488.5,1680.68,22488.5,1680.68,22488.5,1680.69,22488.5,1680.69,22488.5,1680.69,22488.5,1680.69,22488.5,1680.69
lmdif_PowellBadlyScaled	291	0	6.16298e-32	1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615	1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0
lmdif_BrownBadlyScaled	256	0	0	1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1e-06
lmdif_Beale	121	0	1.46521e-25	3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5	2.49987,0.653901,2.49987,0.653901,2.49987,0.653901,2.49987,0.653901,2.49987,0.653901,2.49987,0.653901,2.49987,0.653901,2.49987,0.653901
lmdif_JennrichSampson	213	994.896	994.897	0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829	475.937,475.894,475.937,475.894,475.937,475.894,475.937,475.894,475.937,475.894,475.937,475.894,475.937,475.894,475.937,475.894
lmdif_HelicalValley	35	0	9.69231e-33	1,-6.18577e-18,0	0.1,0.631452,1
lmdif_Bard	21	0.00821487	0.00821488	0.0824106,1.13304,2.34369	0.472941,11.7696,11.3259
lmdif_Gaussian	13	1.12793e-08	1.12793e-08	0.398956,1.00002,-3.34829e-13	0.65051,3.76595,0
lmdif_Meyer	-387	87.9458	868.804	0.00742627,5949.59,337.349	1.91166e-06,0,0.00676072
lmdif_GulfResearchDevelopment	384	0	3.50541e-08	5.37938,36.6053,0.984953	39819,41145.7,1837.58
lmdif_Box3d	29	0	2.27799e-30	1,10,1	2.16805,37.3451,1.82013
lmdif_PowellSingular	287	0	4.83914e-58	3.98842e-15,-3.98842e-16,1.78612e-15,1.78612e-15	0,0.1,0.447214,0
lmdif_Wood	326	0	6.01112e-28	1,1,1,1	0.534245,1.06654,0.534359,1.06654
lmdif_KowalikOsborne	82	0.000307505	0.000307506	0.192806,0.191308,0.123062,0.136074	1.72542,29.6232,12.1971,13.5846
lmdif_BrownDennis	-516	85822.2	107689	-7.17716,11.686,-0.975378,-0.29509	0.0903093,0.0237209,0.339559,0.569962
lmdif_Osborne1	93	5.46489e-05	5.46489e-05	0.37541,1.93582,-1.46466,0.0128675,0.0221228	1.48308,157.674,158.705,0.321153,0.64027
lmdif_Biggs	7	0	0	1,10,1,5,4,3	0,204.524,64.0215,274.359,302.262,211.652
lmdif_Osborne2	148	0.0401377	0.0401683	1.30997,0.431458,0.633631,0.599303,0.753912,0.905584,1.36503,4.8248,2.39882,4.56887,5.67537	0.675508,0.859748,0.481427,1.30486,1.88631,5.0202,6.71742,17.2964,1.14824,1.05994,0.599554
lmdif_Watson	-36	0.00228767	0.00260576	1.88024e-22,1.01364,-0.244321,1.37383,-1.68571,1.09812	1,0.759043,5.37417,14.8652,16.8648,6.7074
lmdif_PenaltyI	-126	9.37629e-06	2.24998e-05	0.250021,0.250012,0.249996,0.250002	273.871,273.843,273.868,273.864
lmdif_PenaltyII	514	9.37629e-06	9.37903e-06	0.199999,0.215019,0.467088,0.514755	1,1757.2,1771.19,2236.35
lmdif_VariablyDimensioned	154	0	1.21393e-24	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	0.999665,0.99866,0.996984,0.994633,0.991605,0.987892,0.983487,0.978381,0.972563,0.966018,0.958733,0.950688,0.941863,0.932236,0.921779,0.910462
lmdif_Trigonometric	527	0	0	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0	1.00006,0.999085,0.998112,0.99714,0.99617,0.995201,0.994235,0.993271,0.992308,0.991348,0.990389,0.989432,0.988477,0.987523,0.986572,0.985622
lmdif_BrownAlmostLinear	54	1	1	-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,17.8906	0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0
lmdif_DiscreteBoundary	69	0	3.73713e-34	-0.0462988,-0.0908017,-0.0889284,-0.086703,-0.0840814,-0.0810162,-0.0774561,-0.0733461,-0.0686267,-0.0632337,-0.0570977,-0.0501438,-0.0422906,-0.0334499,-0.0235258,-0.0124139	1.86639,3.60664,3.48577,3.36158,3.2335,3.10082,2.96264,2.81783,2.66494,2.50202,2.3264,2.13427,1.9198,1.67333,1.3762,0.981057
lmdif_DiscreteIntegral	69	0	1.59314e-32	-0.0284861,-0.0550797,-0.0795978,-0.101833,-0.121548,-0.138475,-0.152302,-0.162673,-0.169172,-0.171318,-0.168542,-0.160174,-0.145417,-0.123314,-0.0927079,-0.0521848	0.995075,0.990626,0.9866,0.982949,0.979634,0.976625,0.973909,0.971491,0.969403,0.967718,0.966569,0.966179,0.966917,0.969379,0.974523,0.983895
lmdif_BroydenTridiagonal	86	0	7.34587e-26	-0.580265,-0.707105,-0.707103,-0.707097,-0.707083,-0.70705,-0.70697,-0.706777,-0.70631,-0.705184,-0.702468,-0.69593,-0.680249,-0.642988,-0.556421,-0.366025	0.206486,0.227555,0.227556,0.227558,0.227562,0.227573,0.227599,0.227661,0.227813,0.228181,0.22908,0.23129,0.236692,0.249269,0.273268,0.288683
lmdif_BroydenBanded	103	0	9.97564e-24	-0.428303,-0.476596,-0.519652,-0.558099,-0.592506,-0.624504,-0.623239,-0.62142,-0.619616,-0.618226,-0.617518,-0.617732,-0.617901,-0.617982,-0.618896,-0.586311	0.210531,0.185079,0.165401,0.150094,0.137961,0.127799,0.128205,0.128803,0.1294,0.129858,0.130092,0.130023,0.129968,0.129943,0.129615,0.140127
lmdif_LinearFullRank	35	0	7.88861e-31	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
lmdif_LinearFullRank1	34	3.63636	3.63636	-176.185,-87.5926,-161.918,-43.2963,-320.277,-80.4592,210.955,-21.1481,-58.5169,-159.638,86.4791,-39.7296,125.641,105.978,52.753,-29.5472	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0016159
lmdif_LinearFullRank0cols0rows	34	5.13793	5.13793	1,102.112,305.423,51.5561,-65.5527,153.211,-281.776,26.2781,110.795,-32.2763,-45.1477,77.1057,0.221514,-140.388,46.907,1	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00209255,0
lmdif_Chebyquad	84	0	3.50827e-25	0.0442053,0.199491,0.235619,0.416047,0.5,0.583953,0.764381,0.800509,0.955795	0.370106,2.06127,1.55376,2.76913,2.92582,2.76943,1.55377,2.06202,0.369586
==23858==
==23858== HEAP SUMMARY:
==23858==     in use at exit: 0 bytes in 0 blocks
==23858==   total heap usage: 345 allocs, 345 frees, 115,024 bytes allocated
==23858==
==23858== All heap blocks were freed -- no leaks are possible
==23858==
==23858== For counts of detected and suppressed errors, rerun with: -v
==23858== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 6 from 6)
*/
