#ifndef RanOpt_hh
#define RanOpt_hh

//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_

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
