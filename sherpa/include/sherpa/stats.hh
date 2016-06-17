// 
//  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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

      // weight is not given in the following url
      // https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html
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
                             const ConstArrayType& model,
                             const ConstArrayType& error,
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
  inline int calc_chi2xspecvar_errors( IndexType num,
                                       const ConstArrayType& yraw,
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


  template <typename DataType, typename IndexType>
  inline int my_calc_w_stat( IndexType num, const DataType* src_raw,
                             const DataType* src_model,
                             const DataType* bkg_raw,
                             const DataType* backscale_ratio,
                             DataType* fvec, 
                             const DataType src_exp_time,
                             const DataType my_bkg_exp_time,
                             const DataType trunc_value ) {

    // 
    // heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSappendixStatistics.html
    //
    //
    // S     = src_raw[ ii ]
    //  i
    //
    // B     = bkg_raw[ ii ]
    //  i
    //
    // t     = src_exp_time
    //  s
    //
    // t     = bkg_exp_time
    //  b
    //
    // t  m  = src_model[ ii ]
    //  s  i
    //

    for ( IndexType ii = num - 1; ii >= 0; --ii ) {

      //
      // Must scale the background area to the source area
      //
      DataType bkg_exp_time = my_bkg_exp_time * backscale_ratio[ ii ];
      DataType src_bkg_time = src_exp_time + bkg_exp_time;
      DataType ln_ts_div_src_bkg_time = 
        std::log( src_exp_time / src_bkg_time );
      DataType ln_tb_div_src_bkg_time = 
        std::log( bkg_exp_time / src_bkg_time );

      DataType msubi = src_model[ ii ] / src_exp_time;
      DataType tb_msubi = bkg_exp_time * msubi;
      DataType ts_msubi = src_exp_time * msubi;

      //
      // If any bin has S and/or B  zero then its contribution to W (W )
      //                 i        i                                   i 
      //
      // is calculated as a special case.
      //
      if ( 0 == src_raw[ ii ] ) {

        fvec[ ii ] = ts_msubi - bkg_raw[ ii ] * ln_tb_div_src_bkg_time;
        //
        // So, if S is zero then:
        //         i
        //
        //              W  =  t  m  - B  ln( t  / ( t  + t ) )
        //               i     s  i    i      b      s    b
        //
        //
        // W_i = t_sm_i-B_i\ln{(t_b/(t_s+t_b))}
        //


      } else if ( 0 == bkg_raw[ ii ] ) {

        if ( msubi < src_raw[ ii ] / src_bkg_time )
          fvec[ ii ] = - tb_msubi - src_raw[ ii ] * ln_ts_div_src_bkg_time;
        else {
          DataType log_src_model = 
            src_model[ ii ] <= 0 ? trunc_value : std::log( src_model[ ii ] );
          fvec[ ii ] = src_model[ ii ] + 
            src_raw[ ii ] * ( std::log( src_raw[ ii ] ) - log_src_model - 1 );
        }

        //
        // If B is zero then there are two special cases. 
        //     i
        //
        // If m  < S  / ( t  + t ) then:
        //     i    i      s    b
        //
        //              W  = - t  m  - S  ln( t  / ( t  + t ) )
        //               i      b  i    i      s      s    b
        //
        // otherwise:
        //
        //              W  = t  m  + S  ( ln( S ) - ln( t  m  ) - 1 )
        //               i    s  i    i        i         s  i 
        //
        // If $m_i < S_i/(t_s+t_b)$ then:
        //
        //   W_i = -t_bm_i-S_i\ln{(t_s/(t_s+t_b))}
        //
        // otherwise:
        //
        //   W_i = t_sm_i+S_i(\ln{S_i}-\ln{(t_sm_i)}-1)
        //

      } else {

        DataType raw_data = src_raw[ ii ] + bkg_raw[ ii ];
        DataType src_bkg_time_msubi = src_bkg_time * msubi;
        DataType tmp1 = src_bkg_time_msubi - raw_data;
        DataType tmp2 = src_bkg_time * bkg_raw[ ii ] * msubi;
        DataType dsubi = std::sqrt( tmp1 * tmp1 + 4.0 * tmp2 );
        DataType fsubi = ( raw_data - src_bkg_time_msubi + dsubi ) / 
          ( 2.0 * src_bkg_time );
        DataType log_model_srcexptime_fsubi =
          src_model[ ii ] + src_exp_time * fsubi <= 0 ? trunc_value :
          std::log( src_model[ ii ] + src_exp_time * fsubi );
        DataType log_bkgexptime_fsubi =
          bkg_exp_time * fsubi <= 0 ? trunc_value : 
          std::log( bkg_exp_time * fsubi );
        fvec[ ii ] = src_model[ ii ] + src_bkg_time * fsubi -
          src_raw[ ii ] * log_model_srcexptime_fsubi -
          bkg_raw[ ii ] * log_bkgexptime_fsubi -
          src_raw[ ii ] * ( 1.0 - std::log( src_raw[ ii ] ) ) -
          bkg_raw[ ii ] * ( 1.0 - std::log( bkg_raw[ ii ] ) );
        
        //
        //  W   =  t  m  + ( t  + t ) f - S ln( t  m + t  f ) - B ln( t  f ) -
        //   i      s  i      s    b   i   i     s  i   s  i     i     b  i
        //
        //                S  ( 1 - ln( S ) - B ( 1 - ln( B ) )
        //                 i            i     i           i 
        //
        // where
        //
        //                S  + B - ( t  + t ) m  + d
        //                 i    i     s    b   i    i
        //         f  =  -----------------------------
        //          i              2 ( t  + t  )
        //                             s     b
        //
        // and
        //                                            2
        //        d  = sqrt( [ ( t  + t ) m - S  - B ]  + 4 ( t  + t  ) B  m )
        //         i              s    b   i   i    i          s    b    i  i
        //
        // Solving for the $f_i$ and substituting gives the profile likelihood:
        //
        // W = 2\sum_{i=1}^N t_sm_i+(t_s+t_b)f_i-S_i\ln{(t_sm_i+t_sf_i)} -B_i\ln{(t_bf_i)}-S_i(1-\ln{S_i})-B_i(1-\ln{B_i})
        //
        // where
        //
        // f_i = {{S_i+B_i-(t_s+t_b)m_i + d_i}\over{2(t_s+t_b)}}
        //
        // and
        //
        // d_i = \sqrt{[(t_s+t_b)m_i-S_i-B_i]^2+4(t_s+t_b)B_im_i}
        //

      }

    }

    return EXIT_SUCCESS;

  }

template <typename ArrayType, typename ConstArrayType, typename DataType,
          typename IndexType, typename iArrayType>
  inline int calc_wstat_stat( IndexType num, const ConstArrayType& yraw,
                              const ConstArrayType& model,
                              const iArrayType& data_size,
                              const ConstArrayType& exposure_time,
                              const ConstArrayType& bkg,
                              const ConstArrayType& backscale_ratio,
                              ArrayType& fvec, DataType& stat,
                              const DataType trunc_value ) {

  const IndexType num_data_sets = data_size.get_size( );

  int offset = 0;

  for ( IndexType ii = 0; ii < num_data_sets; ++ii ) { 

    const double* src_raw = &yraw[ offset ];
    const double* src_model = &model[ offset ];
    const double* bkg_raw = &bkg[ offset ];
    const double* ratio_backscale = &backscale_ratio[ offset ];


    double resp_exp_time = exposure_time[ 2 * ii ];
    double bkg_exp_time  = exposure_time[ 2 * ii + 1 ];

    my_calc_w_stat( data_size[ ii ], src_raw,
                    src_model, bkg_raw, ratio_backscale, &fvec[ offset ],
                    resp_exp_time, bkg_exp_time, trunc_value );

    offset += data_size[ ii ];

  }

  stat = 2.0 *
    sherpa::utils::kahan_sum< ArrayType, DataType, IndexType >( num, fvec );

  DataType sqrt2 = std::sqrt( 2.0 );
  for ( IndexType ii = num - 1; ii >= 0; --ii )
    fvec[ ii ] = sqrt2 * std::sqrt( std::fabs(fvec[ii]) );
  
  return EXIT_SUCCESS;
}

}  }  /* namespace stats, namespace sherpa */


#endif /* __sherpa_stats_hh__ */
