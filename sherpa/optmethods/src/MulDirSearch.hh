#ifndef MulDirSearch_hh
#define MulDirSearch_hh

//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_

//
// The implementation is based on the paper:
//
// Wright, M. H. (1996) "Direct Search Methods: Once Scorned, Now Respectable"
// in Numerical Analysis 1995 (Proceedings of the 1995 Dundee Biennial
// Conference in Numerical Analysis) (D.F. Griffiths and G.A. Watson, eds.),
// 191-208, Addison Wesley Longman, Harlow, United Kingdom.
// http://citeseer.ist.psu.edu/155516.html
//
// Feb 2008 Original version written by D. T. Nguyen
//
// There is still a lot of work to be done for this class, it is not
// quite ready for public consumption :(
//

#include "DirectSearch.hh"
#include "Simplex.hh"

namespace sherpa {

  
  template< typename Fct, typename Data >
  class MulDirSearch : public sherpa::DirectSearch< Fct, Data > {

  public:
    
    MulDirSearch( int numpar, double* par, const double* lo, const double* hi,
		  Fct func, Data xtra )
      : DirectSearch< Fct, Data >( numpar, par, lo, hi, func, xtra ) { }
    
    int operator( )( double* model_par, int verbose, int initsimplex,
		     std::vector< int >& finalsimplex,
		     double tolerance, const double* step,
		     int maxnfev, int& nfev, double& fmin,
		     double* g=NULL, double* h=NULL );
    
  private:

    int contract( int index_smallest, int maxnfev, int& nfev,
		  sherpa::Array2d< double >& contraction );

    int expand( int index_smallest, int maxnfev, int& nfev,
		sherpa::Array2d< double >& expansion );

    int reflect( int index_smallest, int maxnfev, int& nfev,
		 sherpa::Array2d< double >& reflection );

  };

}                                                          // namespace sherpa

#endif
