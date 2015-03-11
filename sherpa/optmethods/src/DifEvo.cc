#ifdef testDifEvo

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


#include "DifEvo.hh"

#include "minpack/LevMar.hh"
#include "NelderMead.hh"

#include "tests/tstopt.hh"

void tstde( Init init, Fct fct, int npar, std::vector<double>& par,
	    std::vector<double>& lo, std::vector<double>& hi,
	    double tol, const char* fct_name, int npop, int maxfev,
	    double xprob, double sfactor ) {

  int nfev, mfcts=0, seed=1357, verbose = 0, size = npop * npar,
    maxnfev = maxfev * npar * size;
  double fmin, answer=0.0;
  //
  // you may think you are clever by eliminating the following overhead
  // and simply use the vector par, but believe me it is necessary to call
  // with de_nm and de with mypar!
  //
  std::vector<double> mypar( npar, 0.0 );

  init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
  sherpa::DifEvo< Fct, void*,
    sherpa::NelderMead< Fct, void* > > de_nm( fct, NULL );
  for ( int ii = 0; ii < npar; ++ii )
    mypar[ ii ] = par[ ii ];
  de_nm( verbose, maxnfev, tol, size, seed, xprob, sfactor, npar, lo, hi,
	 mypar, nfev, fmin );
  print_pars( "DifEvo_nm_", fct_name, nfev, fmin, answer, npar, mypar );

  init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
  sherpa::DifEvo< Fct, void*,
    sherpa::OptFunc< Fct, void* > > de( fct, NULL );
  for ( int ii = 0; ii < npar; ++ii )
    mypar[ ii ] = par[ ii ];
  de( verbose, maxnfev, tol, size, seed, xprob, sfactor, npar, lo, hi,
      mypar, nfev, fmin );
  print_pars( "DifEvo_", fct_name, nfev, fmin, answer, npar, mypar );

}

void tstde_lm( Init init, FctVec fct, int npar, std::vector<double>& par,
	       std::vector<double>& lo, std::vector<double>& hi,
	       double tol, const char* fct_name, int npop, int maxfev,
	       double xprob, double sfactor ) {

  int nfev, mfcts=0, seed=1357, verbose = 0, size = npop * npar,
    maxnfev = maxfev * npar * size;
  double fmin, answer=0.0;
  //
  // you may think you are clever by eliminating the following overhead
  // and simply use the vector par, but believe me it is necessary to call
  // with de_nm and de with mypar!
  //
  std::vector<double> mypar( npar, 0.0 );

  init( npar, mfcts, answer, &par[0], &lo[0], &hi[0] );
  sherpa::DifEvo< FctVec, void*,
    minpack::LevMar< FctVec, void* > > de_lm( fct, NULL, mfcts );
  for ( int ii = 0; ii < npar; ++ii )
    mypar[ ii ] = par[ ii ];
  de_lm( verbose, maxnfev, tol, size, seed, xprob, sfactor, npar, lo, hi,
	 mypar, nfev, fmin );
  print_pars( "DifEvo_lm_", fct_name, nfev, fmin, answer, npar, mypar );
}

int main( int argc, char* argv[] ) {

  int c, uncopt = 1, globalopt = 1;
  while ( --argc > 0 && (*++argv)[ 0 ] == '-' )
    while ( c = *++argv[ 0 ] )
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

  int npar=2;
  if ( argc == 1 )
    npar = atoi( *argv );

  double tol=1.0e-6;
  std::cout << "#\n#:npar = " << npar << "\n";
  std::cout << "#:tol=" << tol << '\n';
  std::cout << "#\n# A negative value for the nfev signifies that the "
    "optimization method did not converge\n#\n";
  std::cout << "name\tnfev\tanswer\tstat\tpar\nS\tN\tN\tN\tN\n";


  int npop=16, maxfev=64;
  double xprob=0.9, sfactor=1.0;

  if ( uncopt ) {
    tst_unc_opt( npar, tol, tstde, npop, maxfev, xprob, sfactor );
    tst_unc_opt( npar, tol, tstde_lm, npop, maxfev, xprob, sfactor );
  }
  if ( globalopt )
    tst_global( npar, tol, tstde, npop, maxfev, xprob, sfactor );
    

  return 0;

}

/*
==1609== Memcheck, a memory error detector
==1609== Copyright (C) 2002-2009, and GNU GPL'd, by Julian Seward et al.
==1609== Using Valgrind-3.5.0 and LibVEX; rerun with -h for copyright info
==1609== Command: tstde
==1609==
#
#:npar = 2
#:tol=1e-06
#
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	stat	par
S	N	N	N	N
DifEvo_nm_Rosenbrock	132	0	9.9559e-08	0.999758,0.999536
DifEvo_Rosenbrock	-1269	0	0.479983	1.69281,2.86563
DifEvo_nm_FreudensteinRoth	1620	0	3.84554e-07	5.00002,4.00001
DifEvo_FreudensteinRoth	-507	0	105.64	18.2546,-0.775174
DifEvo_nm_PowellBadlyScaled	191	0	1.24641e-09	1.13276e-05,8.82805
DifEvo_PowellBadlyScaled	-874	0	1	0,305017
DifEvo_nm_BrownBadlyScaled	391	0	2.96214e-07	1e+06,2.00048e-06
DifEvo_BrownBadlyScaled	-209	0	9.99998e+11	1,1
DifEvo_nm_Beale	877	0	3.15355e-07	2.99976,0.500055
DifEvo_Beale	-1978	0	4.76487	1,-0.733235
DifEvo_nm_JennrichSampson	504	124.362	124.362	0.257808,0.257827
DifEvo_JennrichSampson	-4096	124.362	2020	-71514,-67403.1
DifEvo_nm_HelicalValley	311	0	9.44041e-08	1.00002,6.81226e-05,0.000131188
DifEvo_HelicalValley	-1630	0	2500	-1,0,0
DifEvo_nm_Bard	678	0.00821487	0.00821488	0.0823833,1.13289,2.34378
DifEvo_Bard	-467	0.00821487	17.8095	1,1,-7.76037e+11
DifEvo_nm_Gaussian	-748	1.12793e-08	2.31995e-08	0.398944,0.999662,5.61966e-05
DifEvo_Gaussian	-1486	1.12793e-08	3.88811e-06	0.4,1,0
DifEvo_nm_Meyer	237	87.9458	87.9459	0.00560995,6181.3,345.222
DifEvo_Meyer	-437	87.9458	1.39812e+09	138.542,75704.7,16626.3
DifEvo_nm_GulfResearchDevelopment	804	0	7.16939e-09	14.1092,31.4397,1.21371
DifEvo_GulfResearchDevelopment	2562	0	0.0013055	522.104,-99.5005,1.6251
DifEvo_nm_Box3d	758	0	1.42751e-08	1.00002,10.0026,1.00001
DifEvo_Box3d	3141	0	3.97157e-12	98.2308,98.4796,1.02445e-06
DifEvo_nm_PowellSingular	3344	0	3.76824e-12	-0.00014145,1.40511e-05,-0.000254254,-0.000253504
DifEvo_PowellSingular	-4520	0	95	-3,1,0,-3
DifEvo_nm_Wood	534	0	5.18912e-07	0.999954,0.999933,0.999929,0.9999
DifEvo_Wood	-1237	0	19192	-3,-1,-3,-1
DifEvo_nm_KowalikOsborne	3113	0.000307505	0.000307506	0.192802,0.191566,0.123083,0.136232
DifEvo_KowalikOsborne	-8435	0.000307505	0.0017944	0.25,-7658.21,-8244.23,-3350.5
DifEvo_nm_BrownDennis	203	85822.2	85822.2	-11.5931,13.2026,-0.40607,0.23835
DifEvo_BrownDennis	-2468	85822.2	7.92669e+06	25,5,-5,-1
DifEvo_nm_Osborne1	-2350	5.46489e-05	7.18413e-05	0.369941,1.52707,-1.05201,0.0118033,0.0245497
DifEvo_Osborne1	-4846	5.46489e-05	0.879026	0.5,1.5,-1,0.01,0.02
DifEvo_nm_Biggs	5787	0	0	1,10,1,5,4,3
DifEvo_Biggs	4544	0	0	1,10,1,5,4,3
DifEvo_nm_Osborne2	1769	0.0401377	0.0401787	1.3089,0.429306,0.632834,0.595968,0.748637,0.919646,1.35078,4.8479,2.39916,4.57012,5.67618
DifEvo_Osborne2	-42394	0.0401377	0.382866	1.0763,968.461,0.4195,0.577258,0.281223,581.856,5.01397,7,-0.118868,4.5,5.5
DifEvo_nm_Watson	431	0.00228767	0.00228769	-0.015653,1.01242,-0.233011,1.26064,-1.51402,0.993226
DifEvo_Watson	-3617	0.00228767	30	0,0,0,0,0,0
DifEvo_nm_PenaltyI	-508	9.37629e-06	2.24998e-05	0.249964,0.249672,0.250394,0.250001
DifEvo_PenaltyI	-3794	9.37629e-06	885.063	1,2,3,4
DifEvo_nm_PenaltyII	790	9.37629e-06	9.38448e-06	0.200008,0.148759,0.50296,0.517361
DifEvo_PenaltyII	-3609	9.37629e-06	2.34001	0.5,0.5,0.5,0.5
DifEvo_nm_VariablyDimensioned	340	0	1.6999e-07	0.999875,1.00022
DifEvo_VariablyDimensioned	-208	0	21	1,0
DifEvo_nm_Trigonometric	108	0	4.31325e-08	1.35221,0.631523
DifEvo_Trigonometric	469	0	0.00622617	925639,404486
DifEvo_nm_BrownAlmostLinear	-978	1	7.20419e-08	0.499704,2.00067
DifEvo_BrownAlmostLinear	-872	1	2.8125	0.5,0.5
DifEvo_nm_DiscreteBoundary	179	0	2.48442e-07	-0.138079,-0.181779
DifEvo_DiscreteBoundary	-1294	0	0.128254	-0,-0.222222
DifEvo_nm_DiscreteIntegral	92	0	7.65449e-08	-0.128008,-0.159333
DifEvo_DiscreteIntegral	-1566	0	0.0250711	0,-0.222222
DifEvo_nm_BroydenTridiagonal	1030	0	4.37316e-07	-0.445255,-0.365894
DifEvo_BroydenTridiagonal	-1130	0	13	-1,-1
DifEvo_nm_BroydenBanded	1625	0	4.21199e-07	-0.427387,-0.4272
DifEvo_BroydenBanded	-1313	0	72	-1,-1
DifEvo_nm_LinearFullRank	72	0	1.40066e-07	-1.00026,-0.999727
DifEvo_LinearFullRank	-1132	0	8	1,1
DifEvo_nm_LinearFullRank1	475	0.2	0.2	3.17795,-1.28898
DifEvo_LinearFullRank1	-296	0.2	29	1,1
DifEvo_nm_LinearFullRank0cols0rows	4096	2	2	1,1
DifEvo_LinearFullRank0cols0rows	4096	2	2	1,1
DifEvo_nm_Chebyquad	3072	0	7.26574e-07	0.0443548,0.199942,0.235651,0.41659,0.500148,0.583794,0.764891,0.799883,0.955769
DifEvo_Chebyquad	-5305	0	0.028883	0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9
DifEvo_lm_Rosenbrock	54	0	0	1,1
DifEvo_lm_FreudensteinRoth	-23	0	48.9843	11.4174,-0.896464
DifEvo_lm_PowellBadlyScaled	50	0	6.75462e-29	1.09816e-05,9.10615
DifEvo_lm_BrownBadlyScaled	46	0	1.97215e-31	1e+06,2e-06
DifEvo_lm_Beale	23	0	1.83079e-26	3,0.5
DifEvo_lm_JennrichSampson	41	124.362	124.362	0.257837,0.257818
DifEvo_lm_HelicalValley	35	0	9.69231e-33	1,-6.18577e-18,0
DifEvo_lm_Bard	21	0.00821487	0.00821488	0.0824106,1.13304,2.34369
DifEvo_lm_Gaussian	9	1.12793e-08	1.12793e-08	0.398956,1.00002,-3.34829e-13
DifEvo_lm_Meyer	537	87.9458	87.9459	0.00560966,6181.34,345.224
DifEvo_lm_GulfResearchDevelopment	1729	0	1.7747e-21	50,25,1.5
DifEvo_lm_Box3d	25	0	2.49066e-22	1,10,1
DifEvo_lm_PowellSingular	287	0	4.83914e-58	3.98842e-15,-3.98842e-16,1.78612e-15,1.78612e-15
DifEvo_lm_Wood	-31	0	7.87697	-0.969487,0.950066,-0.968,0.948317
DifEvo_lm_KowalikOsborne	62	0.000307505	0.000307506	0.192813,0.191149,0.123032,0.136001
DifEvo_lm_BrownDennis	1055	85822.2	85823	-11.5697,13.1946,-0.404388,0.236564
DifEvo_lm_Osborne1	93	5.46489e-05	5.46489e-05	0.37541,1.93582,-1.46466,0.0128675,0.0221228
DifEvo_lm_Biggs	7	0	0	1,10,1,5,4,3
DifEvo_lm_Osborne2	148	0.0401377	0.0401683	1.30997,0.431458,0.633631,0.599303,0.753912,0.905584,1.36503,4.8248,2.39882,4.56887,5.67537
DifEvo_lm_Watson	-36	0.00228767	0.00260576	1.88024e-22,1.01364,-0.244321,1.37383,-1.68571,1.09812
DifEvo_lm_PenaltyI	-118	9.37629e-06	2.24998e-05	0.249976,0.249907,0.249945,0.250203
DifEvo_lm_PenaltyII	585	9.37629e-06	9.37629e-06	0.199999,0.191382,0.480075,0.518822
DifEvo_lm_VariablyDimensioned	19	0	0	1,1
DifEvo_lm_Trigonometric	25	0	0	0,0
DifEvo_lm_BrownAlmostLinear	-7	1	1.97215e-31	0.5,2
DifEvo_lm_DiscreteBoundary	13	0	7.89631e-33	-0.138333,-0.18186
DifEvo_lm_DiscreteIntegral	13	0	9.62965e-34	-0.128247,-0.159268
DifEvo_lm_BroydenTridiagonal	16	0	2.56873e-29	-0.445209,-0.366025
DifEvo_lm_BroydenBanded	19	0	1.26733e-26	-0.427305,-0.427305
DifEvo_lm_LinearFullRank	6	0	0	-1,-1
DifEvo_lm_LinearFullRank1	6	0.2	0.2	1,-0.2
DifEvo_lm_LinearFullRank0cols0rows	3	2	2	1,1
DifEvo_lm_Chebyquad	84	0	3.50827e-25	0.0442053,0.199491,0.235619,0.416047,0.5,0.583953,0.764381,0.800509,0.955795
DifEvo_nm_McCormick	948	-1.91	-1.91322	-0.546776,-1.54733
DifEvo_McCormick	913	-1.91	-1.91321	-0.547323,-1.54402
DifEvo_nm_BoxBetts	1978	0	1.99794e-08	1.00001,10.0021,0.99997
DifEvo_BoxBetts	1174	0	9.28277e-06	0.998957,10.0246,1.00266
DifEvo_nm_Paviani	503	-45.7	-45.7784	9.35094,9.35057,9.35122,9.35054,9.35029,9.35182,9.351,9.34902,9.35211,9.35062
DifEvo_Paviani	53542	-45.7	-45.776	9.35841,9.3382,9.34925,9.34385,9.34023,9.34306,9.3535,9.35683,9.34874,9.36638
DifEvo_nm_GoldsteinPrice	940	3	3	-3.46632e-05,-1.00005
DifEvo_GoldsteinPrice	991	3	3.01109	-0.00384429,-1.00518
DifEvo_nm_Shekel5	-786	-10.1532	-5.10077	7.99967,7.99985,7.99985,7.99987
DifEvo_Shekel5	-690	-10.1532	-5.10076	8,8,8,8
DifEvo_nm_Shekel7	3552	-10.4029	-10.4029	4.00089,4.00086,3.99937,3.99956
DifEvo_Shekel7	-868	-10.4029	-5.1288	8,8,8,8
DifEvo_nm_Shekel10	-2449	-10.5364	-5.17564	7.99963,7.99923,7.99955,7.99932
DifEvo_Shekel10	6884	-10.5364	-10.5326	4.0042,3.99791,4.00391,3.99837
DifEvo_nm_Levy4	1457	-21.502	-21.3924	0.67054,1.00058,1.00162,-9.75198
DifEvo_Levy4	4384	-21.502	-21.4998	0.997884,1.02367,1.02269,-9.75371
DifEvo_nm_Levy5	-2147	-11.504	-10.5048	1.00002,1.00001,0.994246,1.00009,-4.25488
DifEvo_Levy5	12857	-11.504	-11.5042	1.00068,0.996768,1.00957,0.993801,-4.75458
DifEvo_nm_Levy6	5178	-11.504	-11.5044	0.999965,1.00132,0.999997,0.998609,0.999548,-4.75449
DifEvo_Levy6	18778	-11.504	-11.5038	0.998553,0.983666,1.00953,0.999233,0.998124,-4.75395
DifEvo_nm_Levy7	-2275	-11.504	-5.39837	1.32966,0.990614,1.00139,0.997245,1.00045,0.999861,-1.75922
DifEvo_Levy7	25584	-11.504	-11.5014	1.00281,0.954325,0.998442,1.00281,0.990281,1.00492,-4.75403
DifEvo_nm_Griewank	86	0	0.00986565	6.28101,0.00143405
DifEvo_Griewank	1112	0	0.00954858	-3.1916,4.49578
DifEvo_nm_SixHumpCamel	1118	-1.03	-1.03163	0.0901221,-0.712666
DifEvo_SixHumpCamel	1301	-1.03	-1.03162	-0.090866,0.71167
DifEvo_nm_Branin	490	0.397889	0.397888	9.42459,2.4742
DifEvo_Branin	758	0.397889	0.397916	3.14261,2.27908
DifEvo_nm_Shubert	570	-24.06	-24.0625	-0.491212,-0.491489
DifEvo_Shubert	-319	-24.06	-20.789	-0.520742,-0.633668
DifEvo_nm_Hansen	685	-176.54	-176.541	4.9764,4.85846
DifEvo_Hansen	644	-176.54	-176.15	-1.29284,-1.43547
DifEvo_nm_Cola	295936	12.8154	12.8668	3.37086,-0.968193,0.0776433,1.34515,0.666296,0.0593934,3.33401,0.676154,3.30015,-2.02212,1.23676,-0.597915,3.38159,-1.12001,3.13852,1.42521,1.48432
DifEvo_Cola	-149662	12.8154	13.1554	1.71311,0.650242,-0.769167,-1.45086,0.415661,-2.7686,-1.82584,-1.83686,-2.91437,0.726565,-2.29956,-2.24795,-2.44269,-1.31676,-3.14506,-1.93621,-0.344594
DifEvo_nm_Ackley	602	0	8.55324e-07	-3.01211e-07,2.68131e-08
DifEvo_Ackley	-1260	0	0.0309546	0.00442412,0.00897125
DifEvo_nm_Bohachevsky1	1226	0	3.12923e-07	-2.3973e-05,-9.52514e-05
DifEvo_Bohachevsky1	1236	0	0.00840666	0.021335,-0.00753833
DifEvo_nm_Bohachevsky2	692	0	2.03681e-06	-7.80486e-05,-0.000275493
DifEvo_Bohachevsky2	-1302	0	0.0159607	-0.00842493,-0.0242418
DifEvo_nm_Bohachevsky3	1471	0	3.04e-08	-0.000120621,7.79945e-05
DifEvo_Bohachevsky3	1327	0	0.00373264	0.00331815,-0.014323
DifEvo_nm_Easom	1296	-1	-0.999999	3.14085,3.14218
DifEvo_Easom	-758	-1	-2.52976e-09	6.86302,0.821265
DifEvo_nm_Rastrigin	398	0	2.45568e-07	-1.45021e-05,-3.20543e-05
DifEvo_Rastrigin	-373	0	0.73034	0.00521053,0.060817
DifEvo_nm_Michalewicz2	508	-1.8013	-1.8013	2.20273,1.57087
DifEvo_Michalewicz2	768	-1.8013	-1.80122	2.20083,1.57136
DifEvo_nm_Michalewicz5	1083	-4.68766	-4.64589	2.20292,1.57053,1.2851,1.11375,1.72044
DifEvo_Michalewicz5	13344	-4.68766	-4.68754	2.205,1.57137,1.28544,1.92321,1.72067
DifEvo_nm_Michalewicz10	12511	-9.66015	-9.61838	2.20326,1.57071,1.28499,1.11359,1.72045,1.57073,1.45438,1.75605,1.65573,1.57082
DifEvo_Michalewicz10	52873	-9.66015	-9.55588	2.20803,1.57377,1.28519,1.92256,1.72089,1.57051,1.4554,1.75435,1.28298,1.21835
==1609==
==1609== HEAP SUMMARY:
==1609==     in use at exit: 0 bytes in 0 blocks
==1609==   total heap usage: 919,570 allocs, 919,570 frees, 157,350,820 bytes allocated
==1609==
==1609== All heap blocks were freed -- no leaks are possible
==1609==
==1609== For counts of detected and suppressed errors, rerun with: -v
==1609== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)
*/

#endif
