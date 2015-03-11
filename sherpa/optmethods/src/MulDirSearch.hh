#ifndef MulDirSearch_hh
#define MulDirSearch_hh

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
