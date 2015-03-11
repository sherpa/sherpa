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
    
    minpack::LevMar< FctVec, void* > lm( fct, NULL, mfcts );

    int maxnfev=128*npar, nfev;
    double fmin=0.0;

    int nprint = 0;
    double epsfcn = 1.0e-8, factor=100.0;
    std::vector<double> covarerr( 8 * npar );

    lm( npar, tol, tol, tol, maxnfev, epsfcn, factor, nprint, lo, hi, par,
	nfev, fmin, covarerr );

    // lm.minimize( &par[0], tol, maxnfev, nfev, fmin );

    print_pars( "lmdif_", fct_name, nfev, fmin, answer, npar, par,
		&covarerr[0] );

  } catch( const sherpa::OptErr& oe ) {
    
    std::cerr << oe << '\n';
    
  }

}

int main( int argc, char* argv[] ) {

  int npar=16;
  if ( argc == 2 )
    npar = atoi( argv[1] );

  if ( npar % 2 || npar < 2 ) {
    printf( "The minimum value for the free parameter must be an even "
            "and it is greater then 2\n" );
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

  return 0;
  
}

/*
==2299== Memcheck, a memory error detector
==2299== Copyright (C) 2002-2009, and GNU GPL'd, by Julian Seward et al.
==2299== Using Valgrind-3.5.0 and LibVEX; rerun with -h for copyright info
==2299== Command: tstlm
==2299==
#
#:npar = 16
#:tol=1.49012e-08
# A negative value for the nfev signifies that the optimization method did not converge
#
name	nfev	answer	fval	par	err
S	N	N	N	N	N
lmdif_Rosenbrock	278	0	0	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026,1,2.0026
lmdif_FreudensteinRoth	-142	0	391.874	11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842,11.4123,-0.896842	26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23,26456.6,1977.23
lmdif_PowellBadlyScaled	274	0	5.93125e-29	1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615,1.09816e-05,9.10615	1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0,1.09816e-05,0
lmdif_BrownBadlyScaled	256	0	0	1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06,1e+06,2e-06	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1e-06
lmdif_Beale	121	0	7.16896e-26	3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5,3,0.5	2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891,2.4998,0.653891
lmdif_JennrichSampson	213	994.896	994.897	0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829,0.25782,0.257829	474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133,474.176,474.133
lmdif_HelicalValley	35	0	9.71834e-33	1,-6.19407e-18,8.90209e-67	0.1,0.0628319,1
lmdif_Bard	21	0.00821487	0.00821488	0.0824106,1.13304,2.34369	0.47294,11.7694,11.3256
lmdif_Gaussian	13	1.12793e-08	1.12793e-08	0.398956,1.00002,-2.24701e-13	0.650508,3.76589,0
lmdif_Meyer	-387	87.9458	802.357	0.00732639,5960.66,337.728	1.88816e-06,0,0.00676482
lmdif_GulfResearchDevelopment	386	0	3.85848e-08	5.05075,36.9509,0.969567	35901.4,39786,1775.51
lmdif_Box3d	29	0	1.99834e-30	1,10,1	2.16802,37.3447,1.82013
lmdif_PowellSingular	283	0	1.69653e-57	5.43907e-15,-5.43907e-16,2.51559e-15,2.51559e-15	0,0.1,0.447214,0
lmdif_Wood	326	0	3.6586e-28	1,1,1,1	0.534248,1.06653,0.534362,1.06653
lmdif_KowalikOsborne	82	0.000307505	0.000307506	0.192808,0.191257,0.123052,0.136051	1.72541,29.6237,12.1973,13.5845
lmdif_BrownDennis	-513	85822.2	101279	-7.85984,11.9153,-0.538758,0.228339	0.0708236,0.022317,0.403789,0.595972
lmdif_Osborne1	93	5.46489e-05	5.46489e-05	0.37541,1.93582,-1.46466,0.0128675,0.0221228	1.48313,157.678,158.709,0.321146,0.6403
lmdif_Biggs	7	0	0	1,10,1,5,4,3	0,204.531,64.0246,274.365,302.265,211.655
lmdif_Osborne2	148	0.0401377	0.0401683	1.30997,0.431458,0.633631,0.599303,0.753912,0.905583,1.36503,4.8248,2.39882,4.56887,5.67537	0.675507,0.859784,0.481429,1.30504,1.88633,5.02027,6.71767,17.2939,1.14824,1.05995,0.599553
lmdif_Watson	-36	0.00228767	0.00260576	-1.31773e-22,1.01364,-0.244317,1.37382,-1.68569,1.09812	1,0.759045,5.37423,14.8655,16.8653,6.70762
lmdif_PenaltyI	-130	9.37629e-06	2.24998e-05	0.250001,0.250017,0.25001,0.250001	273.85,273.867,273.864,273.866
lmdif_PenaltyII	514	9.37629e-06	9.37902e-06	0.199999,0.215021,0.467087,0.514754	1,1757.21,1771.18,2236.35
lmdif_VariablyDimensioned	154	0	1.14323e-25	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1	0.999665,0.998661,0.996986,0.994636,0.991608,0.987895,0.983491,0.978384,0.972565,0.96602,0.958733,0.950688,0.941863,0.932237,0.921783,0.910471
lmdif_Trigonometric	595	0	0	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0	1.00005,0.999251,0.998452,0.997656,0.99686,0.996066,0.995273,0.994481,0.99369,0.992901,0.992113,0.991326,0.99054,0.989756,0.988973,0.988191
lmdif_BrownAlmostLinear	54	1	1	-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,-0.0556601,17.8906	0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0.966227,0
lmdif_DiscreteBoundary	69	0	3.37602e-34	-0.0462988,-0.0908017,-0.0889284,-0.086703,-0.0840814,-0.0810162,-0.0774561,-0.0733461,-0.0686267,-0.0632337,-0.0570977,-0.0501438,-0.0422906,-0.0334499,-0.0235258,-0.0124139	1.86639,3.60664,3.48577,3.36158,3.2335,3.10082,2.96264,2.81783,2.66494,2.50202,2.3264,2.13427,1.9198,1.67333,1.3762,0.981057
lmdif_DiscreteIntegral	69	0	1.06592e-32	-0.0284861,-0.0550797,-0.0795978,-0.101833,-0.121548,-0.138475,-0.152302,-0.162673,-0.169172,-0.171318,-0.168542,-0.160174,-0.145417,-0.123314,-0.0927079,-0.0521848	0.995075,0.990626,0.986601,0.982949,0.979634,0.976625,0.973909,0.971491,0.969403,0.967718,0.966569,0.966179,0.966918,0.969379,0.974523,0.983895
lmdif_BroydenTridiagonal	86	0	2.09653e-26	-0.580265,-0.707105,-0.707103,-0.707097,-0.707083,-0.70705,-0.70697,-0.706777,-0.70631,-0.705184,-0.702468,-0.69593,-0.680249,-0.642988,-0.556421,-0.366025	0.206484,0.227554,0.227554,0.227556,0.227561,0.227571,0.227597,0.227659,0.227811,0.228179,0.229078,0.231288,0.23669,0.249267,0.273266,0.288681
lmdif_BroydenBanded	103	0	3.54099e-24	-0.428303,-0.476596,-0.519652,-0.558099,-0.592506,-0.624504,-0.623239,-0.62142,-0.619616,-0.618226,-0.617518,-0.617732,-0.617901,-0.617982,-0.618896,-0.586311	0.210528,0.185077,0.165399,0.150092,0.137958,0.127797,0.128203,0.128801,0.129398,0.129856,0.13009,0.130021,0.129966,0.129941,0.129613,0.140125
lmdif_LinearFullRank	35	0	1.97215e-31	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
lmdif_LinearFullRank1	35	3.63636	3.63636	-4883.01,-822.335,516.085,-518.777,53.0791,-21.7553,435.962,157.453,-22.3739,11.9037,267.798,-109.554,-80.2355,230.792,-143.096,63.6255	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0016159
lmdif_LinearFullRank0cols0rows	35	5.13793	5.13793	1,-2412.19,-251.123,262.525,490.5,-327.764,375.448,374.867,-15.7958,-172.354,160.149,200.502,51.161,-173.107,-141.919,1	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00209255,0
lmdif_Chebyquad	93	0	3.70419e-23	0.0442053,0.199491,0.235619,0.416047,0.5,0.583953,0.764381,0.800509,0.955795	0.370101,2.06122,1.55372,2.76901,2.92569,2.76925,1.55373,2.06184,0.369676
==2299==
==2299== HEAP SUMMARY:
==2299==     in use at exit: 0 bytes in 0 blocks
==2299==   total heap usage: 399 allocs, 399 frees, 384,556 bytes allocated
==2299==
==2299== All heap blocks were freed -- no leaks are possible
==2299==
==2299== For counts of detected and suppressed errors, rerun with: -v
==2299== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 4 from 4)
 */

#endif
