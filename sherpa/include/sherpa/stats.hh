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

#ifndef __sherpa_stats_hh__
#define __sherpa_stats_hh__

#include <cstdlib>
#include <cmath>

#include <sherpa/utils.hh>


namespace sherpa { namespace stats {


  //#define MODELFUDGEVAL 1.0e-25

  //
  // fvec is needed for chi square and lmdif, it is not needed for cStat
  //
  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_cstat_stat( IndexType num, const ConstArrayType& yraw,
		       const ConstArrayType& model,
		       const ConstArrayType& error,
		       const ConstArrayType& syserror,
		       const ConstArrayType& weight,
		       ArrayType& fvec, DataType& stat,
		       DataType& trunc_value ) {

    DataType mymodel;
    for ( IndexType ii = num - 1; ii >= 0; --ii ) {

      if ( model[ ii ] > 0.0) 
	mymodel = model[ ii ];
      else {
	if (trunc_value > 0)
	  mymodel = trunc_value;
	else
	  return EXIT_FAILURE;
      }

      if( yraw[ ii ] > 0.0 )
	fvec[ ii ] = mymodel - yraw[ ii ] +
	  yraw[ ii ] * ( std::log( yraw[ ii ] / mymodel ) );
      else if ( yraw[ ii ] == 0.0 )
	fvec[ ii ] = mymodel;
      else
	return EXIT_FAILURE;
      
      if ( weight )
	fvec[ ii ] *= weight[ ii ];
     
    }
 
    stat = 2.0 *
      sherpa::utils::kahan_sum< ArrayType, DataType, IndexType >( num, fvec );

    DataType sqrt2 = std::sqrt( 2.0 );
    for ( IndexType ii = num - 1; ii >= 0; --ii )
      fvec[ ii ] = sqrt2 * std::sqrt( std::fabs(fvec[ii]) );

    return EXIT_SUCCESS;

  }


  //
  // fvec is needed for chi square and lmdif, it is not needed for Cash
  //
  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_cash_stat( IndexType num, const ConstArrayType& yraw,
		      const ConstArrayType& model, const ConstArrayType& error,
		      const ConstArrayType& syserror,
		      const ConstArrayType& weight, ArrayType& fvec,
		      DataType& stat, DataType& trunc_value ) {

    DataType mymodel, d;
    for ( IndexType ii = num - 1; ii >= 0; --ii ) {
      
      if ( model[ ii ] > 0.0) 
	mymodel = model[ ii ];
      else {
	if (trunc_value > 0)
	  mymodel = trunc_value;
	else
	  return EXIT_FAILURE;
      }
      
      if ( 0.0 == yraw[ ii ] )
      	d = mymodel;
      else
	d = mymodel - ( yraw[ ii ] * std::log( mymodel ) );
      
      if ( weight )
	d *= weight[ ii ];
 
      fvec[ ii ] = d;

    }
 
    stat = 2.0 *
      sherpa::utils::kahan_sum< ArrayType, DataType, IndexType >( num, fvec );

    {
      DataType junkcstat;
      return calc_cstat_stat( num, yraw, model, error, syserror, weight, fvec,
			      junkcstat, trunc_value );
    }

    return EXIT_SUCCESS;
 
  }


  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2gehrels_errors( IndexType num, const ConstArrayType& yraw,
			       ArrayType& err ) {

    const DataType fudge = 0.75;
    const DataType sqrt_1_75 = 1.0 + std::sqrt( fudge );
  
    for ( IndexType ii = num - 1; ii >= 0; --ii ) {
   
      DataType tmp = yraw[ ii ] + fudge;

      if ( tmp < 0.0 )
	err[ ii ] = sqrt_1_75;
      else
	err[ ii ] = 1.0 + std::sqrt( tmp );

    }

    return EXIT_SUCCESS;

  }

  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_stat( IndexType num, const ConstArrayType& weight,
		 ArrayType& fvec, DataType& stat ) {

    if ( weight )
      for ( IndexType ii = num - 1; ii >= 0; --ii )
	if ( weight[ ii ] >= 0.0 )
	  fvec[ ii ] *= std::sqrt( weight[ ii ] );
	else
	  return EXIT_FAILURE;

    stat = sherpa::utils::enorm2< ArrayType, DataType, IndexType >( num, fvec );    
    /*
    stat = 0.0;
    for ( IndexType ii = num - 1; ii >= 0; --ii )
      stat += fvec[ ii ] * fvec[ ii ];
    */

    return EXIT_SUCCESS;

  }

  //
  // upon return from the function calculate_chi2_stat, fvec 
  // (an array of length num) shall contain the following values:
  //
  //               sqrt( weight[ ii ] )  *  ( data[ ii ] - model[ ii ] )
  // fvec[ ii ] =  -----------------------------------------------------
  //                   sqrt( syserror[ ii ]^2  +  error[ ii ]^2 )
  //
  // the return value from the function is the statistic
  //
  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2_stat( IndexType num, const ConstArrayType& yraw,
		      const ConstArrayType& model,
		      const ConstArrayType& error, 
		      const ConstArrayType& syserror,
		      const ConstArrayType& weight, ArrayType& fvec,
		      DataType& stat, DataType& trunc_value ) {

    for ( IndexType ii = num - 1; ii >= 0; --ii ) {

      fvec[ ii ] = model[ ii ] - yraw[ ii ];

      DataType err = error[ ii ];
    
      if ( syserror ) {
	err *= error[ ii ];
	err += syserror[ ii ] * syserror[ ii ];
	err = std::sqrt( err );
      }

      // the error for Chi2 is 0, so a sanity check is required.
      if ( 0.0 != err )
	fvec[ ii ] /= err;
    
    }

    return sherpa::stats::calc_stat< ArrayType, ConstArrayType, DataType, IndexType >( num, weight, fvec, stat );
  
  }


  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2constvar_errors( IndexType num, const ConstArrayType& yraw,
				ArrayType& err ) {

    DataType mu = sherpa::utils::kahan_sum< ArrayType, DataType, IndexType >
      ( num, yraw );

    if ( mu < 0.0 )
      mu = DataType( num );
    mu = std::sqrt( mu / DataType( num ) );

    for ( IndexType ii = num - 1; ii >= 0; --ii )
      err[ ii ] = mu;

    return EXIT_SUCCESS;

  }


  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2datavar_errors( IndexType num, const ConstArrayType& yraw,
				ArrayType& err ) {

    for ( IndexType ii = num - 1; ii >= 0; --ii )
      if ( yraw[ ii ] <= 0.0 )
	return EXIT_FAILURE;
      else
	err[ ii ] = std::sqrt( yraw[ ii ] );
    
    return EXIT_SUCCESS;

  }


  //
  // upon return from the function calculate_ChiModVar_stat, fvec 
  // (an array of length num) shall contain the following values:
  //
  //               sqrt( weight[ ii ] )  *  ( data[ ii ] - model[ ii ] )
  // fvec[ ii ] =  -----------------------------------------------------
  //                   sqrt( syserror[ ii ]^2  +  error[ ii ]^2 )
  //
  // the return value from the function is the statistics
  //
  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2modvar_stat( IndexType num, const ConstArrayType& yraw,
			    const ConstArrayType& model,
			    const ConstArrayType& error,
			    const ConstArrayType& syserror,
			    const ConstArrayType& weight,
			    ArrayType& fvec, DataType& stat,
			    DataType& trunc_value ) {

    for ( IndexType ii = num - 1; ii >= 0; --ii ) {

      fvec[ ii ] = yraw[ ii ] - model[ ii ];
    
      DataType err_sqr = model[ ii ];
      if ( err_sqr < 1.0 )
	err_sqr = 1.0;
      // err_sqr is guarranteed to be >= 1.0

      if ( syserror )
	err_sqr += syserror[ ii ] * syserror[ ii ];

      // syserror^2 is always >= 0.0, so need to check if err_sqr < 0.0
      DataType err = std::sqrt( err_sqr );

      // the error is guaranteed to be > 0.0 so err could never be equal to 0.0
      fvec[ ii ] /= err;
    
    }

    return sherpa::stats::calc_stat< ArrayType, ConstArrayType, DataType, IndexType >( num, weight, fvec, stat );

  }


  // Implement XSPEC variation of data variance as a separate
  // option here.  Instead of returning failure when yraw[ii]
  // is zero or less, just set the error to 1.0 instead.  This
  // is added simply to accommodate those used to this XSPEC
  // behavior.
  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_chi2xspecvar_errors( IndexType num, const ConstArrayType& yraw,
				ArrayType& err ) {

    for ( IndexType ii = num - 1; ii >= 0; --ii )
      if ( yraw[ ii ] <= 0.0 )
	err[ ii ] = 1.0;
      else
	err[ ii ] = std::sqrt( yraw[ ii ] );
    
    return EXIT_SUCCESS;

  }


  template <typename ArrayType, typename ConstArrayType, typename DataType,
	    typename IndexType>
  inline int calc_lsq_stat( IndexType num, const ConstArrayType& yraw,
		      const ConstArrayType& model,
		      const ConstArrayType& error, 
		      const ConstArrayType& syserror,
		      const ConstArrayType& weight, ArrayType& fvec,
		      DataType& stat, DataType& trunc_value ) {

    for ( IndexType ii = num - 1; ii >= 0; --ii )
      fvec[ ii ] = model[ ii ] - yraw[ ii ];

    stat = sherpa::utils::enorm2< ArrayType, DataType, IndexType >( num, fvec );
    return EXIT_SUCCESS;

  }


}  }  /* namespace stats, namespace sherpa */


#endif /* __sherpa_stats_hh__ */
