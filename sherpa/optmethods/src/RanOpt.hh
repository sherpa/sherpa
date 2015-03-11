#ifndef RanOpt_hh
#define RanOpt_hh

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


#include "sherpa/MersenneTwister.h"

namespace sherpa {

  /*
  class MyMTRand : public MTRand {

  public:

    MyMTRand( const uint32 oneSeed ) : MTRand( oneSeed ) { }

    // return a pseudo random number in the [ 0, 1 ] real-interval
    double random_number( ) { return rand( ); }

    // return a pseudo random number in the ( 0, 1 ) real-interval
    double random_number_ee( ) { return randDblExc( ); }

    // return a pseudo random number in the [ 0, high ] integer-interval
    MTRand::uint32 random_number( const MTRand::uint32 n ) { 
      return randInt( n ); }

    // return a pseudo random number in the ( lo, high ) real-interval
    double random_number_ee( double low, double high ) { 
      return low + ( high - low ) * randDblExc( );
    }

    double gaussian_deviate( const double mean, const double stddev ) {
      return randNorm( mean, stddev ); }

  };
  */

}                                                           // namespace sherpa

#endif
