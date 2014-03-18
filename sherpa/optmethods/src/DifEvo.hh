#ifndef DifEvo_hh
#define DifEvo_hh

//
// An implementation of the Differential Evolution for continous
// function optimization, an algorithm by Kenneth Price and Rainer Storn.
// See: http://www.icsi.berkeley.edu/~storn/code.html
//
//  /***************************************************************
//  **                                                            **
//  **        D I F F E R E N T I A L     E V O L U T I O N       **
//  **                                                            **
//  ** Program: de.c                                              **
//  ** Version: 3.6                                               **
//  **                                                            **
//  ** Authors: Dr. Rainer Storn                                  **
//  **          c/o ICSI, 1947 Center Street, Suite 600           **
//  **          Berkeley, CA 94707                                **
//  **          Tel.:   510-642-4274 (extension 192)              **
//  **          Fax.:   510-643-7684                              **
//  **          E-mail: storn@icsi.berkeley.edu                   **
//  **          WWW: http://http.icsi.berkeley.edu/~storn/        **
//  **          on leave from                                     **
//  **          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6          **
//  **          D-81739 Muenchen, Germany                         **
//  **          Tel:    636-40502                                 **
//  **          Fax:    636-44577                                 **
//  **          E-mail: rainer.storn@zfe.siemens.de               **
//  **                                                            **
//  **          Kenneth Price                                     **
//  **          836 Owl Circle                                    **
//  **          Vacaville, CA 95687                               **
//  **          E-mail: kprice@solano.community.net               ** 
//  **                                                            **
//  ** This program implements some variants of Differential      **
//  ** Evolution (DE) as described in part in the techreport      **
//  ** tr-95-012.ps of ICSI. You can get this report either via   **
//  ** ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z  **
//  ** or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
//  ** A more extended version of tr-95-012.ps is submitted for   **
//  ** publication in the Journal Evolutionary Computation.       ** 
//  **                                                            **
//  ** You may use this program for any purpose, give it to any   **
//  ** person or change it according to your needs as long as you **
//  ** are referring to Rainer Storn and Ken Price as the origi-  **
//  ** nators of the the DE idea.                                 **
//  ** If you have questions concerning DE feel free to contact   **
//  ** us. We also will be happy to know about your experiences   **
//  ** with DE and your suggestions of improvement.               **
//  **                                                            **
//  ***************************************************************/
//
// Differential Evolution Solver Class
// Based on algorithms developed by Dr. Rainer Storn & Kenneth Price
// Written By: Lester E. Godwin
//             PushCorp, Inc.
//             Dallas, Texas
//             972-840-0208 x102
//             godwin@pushcorp.com
// Created: 6/8/98
// Last Modified: 6/8/98
// Revision: 1.0
//

#include "sherpa/MersenneTwister.h"

#include "Opt.hh"
#include "Simplex.hh"

namespace sherpa {

  template < typename Func, typename Data, typename Algo >
  class DifEvo : public sherpa::Opt {
    
  public:

    typedef DifEvo<Func,Data,Algo> MyDifEvo;
    typedef void (MyDifEvo::*StrategyFuncPtr)( int, double, double, int,
					       const sherpa::Simplex&,
					       const std::vector<double>&,
					       MTRand&, std::vector<double>& );

    enum Strategy { Best1Exp, Rand1Exp, RandToBest1Exp, Best2Exp, Rand2Exp,
		    Best1Bin, Rand1Bin, RandToBest1Bin, Best2Bin, Rand2Bin };

    DifEvo( Func func, Data xdata, int mfct=0 )
      : Opt( ), usr_func( func ), usr_data( xdata ),
	local_opt( func, xdata, mfct ), strategy_func_ptr( 0 ) { }
	
    int operator( )( int verbose, int maxnfev, double tol, int population_size,
		     int seed, double cross_over_probability, 
		     double scale_factor, int npar,
		     const std::vector<double>& low,
		     const std::vector<double>& high,
		     std::vector<double>& par, int& nfev, double& fmin ) {

      int ierr = EXIT_SUCCESS;

      nfev = 0;
      fmin = std::numeric_limits< double >::max( );
      std::vector<double> mypar( npar + 1, 0.0 );
      for ( int ii = 0; ii < npar; ++ii )
	mypar[ ii ] = par[ ii ];

      try {

	const sherpa::Opt::mypair limits( low, high );
	if ( sherpa::Opt::are_pars_outside_limits( npar, limits, par ) )
	  throw sherpa::OptErr( sherpa::OptErr::OutOfBound );

	ierr = difevo( verbose, maxnfev, tol, population_size, seed,
		       cross_over_probability, scale_factor, npar, limits,
		       mypar, nfev );

      } catch( OptErr& oe ) {

	if ( verbose )
	  std::cerr << oe << '\n';
	ierr = oe.err;

      } catch( std::runtime_error& re ) {

	if ( verbose )
	  std::cerr << re.what( ) << '\n';
	ierr = OptErr::Unknown;

      } catch( std::exception& e ) {

	if ( verbose )
	  std::cerr << e.what( ) << '\n';
	ierr = OptErr::Unknown;

      }

      for ( int ii = 0; ii < npar; ++ii )
	par[ ii ] = mypar[ ii ];
      fmin = mypar[ npar ];
      return ierr;

    }


  private:
    Func usr_func;
    Data usr_data;
    Algo local_opt;
    StrategyFuncPtr strategy_func_ptr;

    void choose_strategy( int strategy ) {

      switch ( strategy  ) {
      case Best1Exp:
	// strategy DE0 not in the paper
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1exp;
	break;
      case Rand1Exp:
	// strategy DE1 in the techreport
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand1exp;
	break;
      case RandToBest1Exp:
	// similiar to DE2 but generally better
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::randtobest1exp;
	break;
      case Best2Exp:
	// is another powerful strategy worth trying
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best2exp;
	break;
      case Rand2Exp:
	// seems to be a robust optimizer for many functions
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand2exp;
	break;
      case Best1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1bin;
	break;
      case Rand1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand1bin;
	break;
      case RandToBest1Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::randtobest1bin;
	break;
      case Best2Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best2bin;
	break;
      case Rand2Bin:
	// Essentially same strategies but BINOMIAL CROSSOVER
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::rand2bin;
	break;
      default:
	strategy_func_ptr = &DifEvo<Func,Data,Algo>::best1exp;
	break;

      }

      //
      // Choice of strategy
      // We have tried to come up with a sensible naming-convention: DE/x/y/z
      // DE :  stands for Differential Evolution
      // x  :  a string which denotes the vector to be perturbed
      // y  :  number of difference vectors taken for perturbation of x
      // z  :  crossover method (exp = exponential, bin = binomial)
      //
      // There are some simple rules which are worth following:
      // 1)  F is usually between 0.5 and 1 (in rare cases > 1)
      // 2)  CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to be 
      //     tried first
      // 3)  To start off NP = 10*D is a reasonable choice. Increase NP if 
      //     misconvergence happens.                                       
      // 4)  If you increase NP, F usually has to be decreased
      // 5)  When the DE/best... schemes fail DE/rand... usually works and
      //     vice versa
      //
   
    }

    int difevo( int verbose, int maxnfev, double tol, int population_size,
		int seed, double cross_over_probability, double scale_factor,
		int npar, const sherpa::Opt::mypair& limits,
		std::vector<double>& par, int& nfev ) {
      

      int ierr = EXIT_SUCCESS;
      par[ npar ]  = std::numeric_limits< double >::max( );
      population_size = std::abs( population_size );

      MTRand mt_rand( seed );

      //
      // For each row of the 2d-array population and children:
      // the columns [ 0, npar - 1 ] contain the parameters, and
      // the column npar contains the function values:
      // (*usrfunc)( population(ii,0), ... population(ii,npar-1) ) =
      //                                                   population(ii,npar);
      //
      // The array shall have dimension: population(population_size, npar + 1)
      //
      // Will use the class Simplex since it has all the the infrastructure
      // that is needed to check for convergence although it is not a simplex
      // in the classic sense.
      //
      const std::vector<double>& low = limits.first;
      const std::vector<double>& high = limits.second;
      sherpa::Simplex population( population_size, npar + 1 );
      for ( int ii = 0; ii < population_size; ++ii ) {
	for ( int jj = 0; jj < npar; ++jj )
	  population[ ii ][ jj ] =
	    low[ jj ] + ( high[ jj ] - low[ jj ] ) * mt_rand.randDblExc( );
	population[ ii ][ npar ] = std::numeric_limits< double >::max( );
      }

      //
      // allocate an extra element to store the function value
      //
      std::vector<double> trial_solution( npar + 1 );
      const double tol_sqr = tol * tol;
      const int simplex_tst = 0;

      ierr = local_opt.minimize( maxnfev - nfev, limits, tol, npar, par, 
				 par[ npar ], nfev );
      if ( EXIT_SUCCESS != ierr )
	return ierr;

      for ( ; nfev < maxnfev; ) {

	for ( int candidate=0; candidate < population_size && nfev < maxnfev;
	      ++candidate ) {

	  population.copy_row( candidate, trial_solution );

	  for ( int strategy = 0; strategy < 10; ++strategy ) {

	    choose_strategy( strategy );

	    (this->*strategy_func_ptr)( candidate, cross_over_probability,
					scale_factor, npar, population, par,
					mt_rand, trial_solution );

	    trial_solution[ npar ] = 
	      local_opt.eval_func( maxnfev, limits, npar, trial_solution,
				   nfev );

	    if ( trial_solution[ npar ] <
		 population[ candidate ][ npar ] ) {
	      population.copy_row( trial_solution, candidate );

	      if ( trial_solution[ npar ] < par[ npar ] ) {

		ierr = local_opt.minimize( maxnfev - nfev, limits, tol, npar,
					   trial_solution,
					   trial_solution[ npar ], nfev );
		if ( EXIT_SUCCESS != ierr )
		  return ierr;

		sherpa::Array2d<double>::copy_vector( npar + 1, trial_solution,
						      par );
		if ( verbose > 1 )
		  sherpa::Opt::print_par( std::cout, par );

	      }  // if ( trial_solution[ npar ] < par[ npar ] ) {

	      population.sort( );
	      if ( population.check_convergence( tol, tol_sqr, simplex_tst ) )
		return EXIT_SUCCESS;

	    }                  // if ( trial_solution[ npar ] < population( ...

	  }              // for ( int strategy = 0; strategy < 10; ++strategy )

	}              // for ( int candidate=0; candidate < population_size &&

      }                                            // for ( ; nfev < maxnfev; )

      return ierr;

    }                                                                // difevo


    //
    // EXPONENTIAL CROSSOVER
    //
    // DE/best/1/exp
    // Our oldest strategy but still not bad. However, we have found several
    // optimization problems where misconvergence occurs.
    //    
    void best1exp( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution )  {

      int r1, r2;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; mt_rand.rand( ) < xprob && ii < npar; ++ii ) {
	trial_solution[ n ] = par[ n ] +
	  sfactor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = ( n + 1 ) % npar;
      }

      return;
    }

    //
    // DE/rand/1/exp
    // This is one of my favourite strategies. It works especially well when
    // the "bestit[]"-schemes experience misconvergence. Try e.g.
    // F=0.7 and CR=0.5 as a first guess.
    //
    void rand1exp( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2, &r3 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; mt_rand.rand( ) < xprob && ii < npar; ++ii ) {
	trial_solution[ n ] = population[ r1 ][ n ] +
	  + sfactor * ( population[ r2 ][ n ] - population[ r3 ][ n ] );
	n = (n + 1) % npar;
      }
      
      return;
      
    }

    //
    // DE/rand-to-best/1/exp
    // This strategy seems to be one of the best strategies. Try F=0.85 and
    // CR=1. If you get misconvergence try to increase NP. If this doesn't
    // help you should play around with all three control variables.
    //
    void randtobest1exp( int candidate, double xprob, double sfactor, int npar,
			 const sherpa::Simplex& population,
			 const std::vector<double>& par, MTRand& mt_rand,
			 std::vector<double>& trial_solution ) {

      int r1, r2;  
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; mt_rand.rand( ) < xprob && ii < npar; ++ii ) {
	trial_solution[n] += sfactor * ( par[ n ] - trial_solution[ n ] ) +
	  sfactor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void best2exp( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3, r4;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2, &r3,
		      &r4 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; mt_rand.rand( ) < xprob && ii < npar; ++ii ) {
	trial_solution[n] = par[ n ] + 
	  sfactor * ( population[ r1 ][ n ] + population[ r2 ][ n ] -
			   - population[ r3 ][ n ] - population[ r4 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void rand2exp( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3, r4, r5;
      select_samples( candidate, population.nrows( ), mt_rand,
		      &r1, &r2, &r3, &r4, &r5 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; mt_rand.rand( ) < xprob && ii < npar; ++ii ) {
	trial_solution[n] = population[ r1 ][ n ] +
	  sfactor * ( population[ r2 ][ n ] + population[ r3 ][ n ] -
			   population[ r4 ][ n ] - population[ r5 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void best1bin( int candidate, double xprob, double sfactor, int npar,
			  const sherpa::Simplex& population,
			  const std::vector<double>& par, MTRand& mt_rand,
			  std::vector<double>& trial_solution ) {

      int r1, r2;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2 );
      int n = mt_rand.randInt( npar - 1 );  
      for ( int ii = 0; ii < npar; ++ii ) {
	if ( mt_rand.rand( ) < xprob || npar - 1 == ii )
	  trial_solution[n] = par[ n ] +
	    sfactor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void rand1bin( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2, &r3 );
      int n = mt_rand.randInt( npar - 1 );  
      for ( int ii = 0; ii < npar; ++ii ) {
	if ( mt_rand.rand( ) < xprob || npar - 1 == ii )
	  trial_solution[n] = population[ r1 ][ n ] +
	    sfactor * ( population[ r2 ][ n ] - population[ r3 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void randtobest1bin( int candidate, double xprob, double sfactor, int npar,
			 const sherpa::Simplex& population,
			 const std::vector<double>& par, MTRand& mt_rand,
			 std::vector<double>& trial_solution ) {

      int r1, r2;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2 );
      int n = mt_rand.randInt( npar - 1 );  
      for ( int ii = 0; ii < npar; ++ii ) {
	if ( mt_rand.rand( ) < xprob || npar - 1 == ii )
	  trial_solution[n] +=
	    sfactor * (par[ n ] - trial_solution[ n ] ) +
	    sfactor * ( population[ r1 ][ n ] - population[ r2 ][ n ] );
	n = (n + 1) % npar;
      }
  
      return;

    }

    void best2bin( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3, r4;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2, &r3,
		      &r4 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; ii < npar; ++ii ) {
	if ( mt_rand.rand( ) < xprob || npar - 1 == ii )
	  trial_solution[n] = par[ n ] +
	    sfactor * ( population[ r1 ][ n ] + population[ r2 ][ n ] -
			     population[ r3 ][ n ] - population[ r4 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }

    void rand2bin( int candidate, double xprob, double sfactor, int npar,
		   const sherpa::Simplex& population,
		   const std::vector<double>& par, MTRand& mt_rand,
		   std::vector<double>& trial_solution ) {

      int r1, r2, r3, r4, r5;
      select_samples( candidate, population.nrows( ), mt_rand, &r1, &r2, &r3,
		      &r4, &r5 );
      int n = mt_rand.randInt( npar - 1 );
      for ( int ii = 0; ii < npar; ++ii ) {
	// perform npar binomial trials
	if ( mt_rand.rand( ) < xprob || npar - 1 == ii )
	  trial_solution[n] = population[ r1 ][ n ] + 
	    sfactor * ( population[ r2 ][ n ] + population[ r3 ][ n ] -
			population[ r4 ][ n ] - population[ r5 ][ n ] );
	n = (n + 1) % npar;
      }

      return;

    }


    static void select_samples( int candidate, int npop, MTRand& mt_rand,
				int* r1, int* r2=0, int* r3=0, int* r4=0,
				int* r5=0 ) {
      if ( r1 ) {
	do {
	  *r1 = mt_rand.randInt( npop - 1 );
	} while (*r1 == candidate);
      }
  
      if ( r2 ) {
	do {
	  *r2 = mt_rand.randInt( npop - 1 );
	} while ( (*r2 == candidate) || (*r2 == *r1) );
      }
  
      if ( r3 ) {
	do {
	  *r3 = mt_rand.randInt( npop - 1 );
	}
	while ( (*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1) );
      }
  
      if ( r4 ) {
	do {
	  *r4 = mt_rand.randInt( npop - 1 );
	} while ( (*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) ||
		  (*r4 == *r1) );
      }
  
      if ( r5 ) {
	do {
	  *r5 = mt_rand.randInt( npop - 1 );
	} while ( (*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3) ||
		  (*r5 == *r2) || (*r5 == *r1) );
      }  
  
      return;

    }

  };                                                            // class DifEvo

}                                                                  // namespace

#endif
