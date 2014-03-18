//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#ifndef __sherpa_astro_utils_hh__
#define __sherpa_astro_utils_hh__

#include "sherpa/extension.hh"
#include "fcmp.h"
#include <cstdlib>
#include <cfloat>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;

#ifndef MID
#define MID( a, b ) (( a + b ) / 2.0 )
#endif

namespace sherpa { namespace astro { namespace utils {


  template <typename ArrayType, typename ConstArrayType, typename IndexType>
  void arf_fold( IndexType num,
		 const ConstArrayType& in1, const ConstArrayType& in2,
		 ArrayType& out )
  {

    for ( IndexType ii = 0; ii < num; ii++ )
      out[ ii ] = in1[ ii ] * in2[ ii ];

  }


  //
  // Function to perform XSPEC-style convolution, using an RMF
  // Extracted from the old Sherpa code for doing so
  //
  // Arguments:
  //
  // len_source is the number of energy bins used to calculate the source
  // source is an array that has calculated source model values
  //    (photons/s/cm**2/keV -- usual units of XSPEC models)
  // num_chans is the number of channels into which photons are
  //    redistributed
  // counts is an array that has the number of counts you
  //    get after "folding" source through the response; this array
  //    is what you compare to the data when calculating chi-squared,
  //    so num_chans better be equal to the number of data points or
  //    bins in the data.
  //
  // The caller must ensure that, on input, counts is filled with zeros
  // and the various integer arrays contain no negative values.
  //
  // Relationships between argument sizes (the function will verify
  // these):
  //
  //   len_num_groups	== len_source
  //   len_first_chan	== sum( num_groups )
  //   len_num_chans	== sum( num_groups )
  //   len_response	== sum( num_chans )
  //

  template <typename ConstIntType, typename ConstFloatType,
	    typename FloatType, typename IndexType, typename UIndexType>
  int rmf_fold( IndexType len_source, const ConstFloatType *source, 
		IndexType len_num_groups, const ConstIntType *num_groups,
		IndexType len_first_chan, const ConstIntType *first_chan,
		IndexType len_num_chans, const ConstIntType *num_chans,
		IndexType len_response, const ConstFloatType *resp,
		IndexType len_counts, FloatType *counts,
		UIndexType offset)
  {

    //int flag = 0;
    if ( ( len_num_groups != len_source ) || 
	 ( len_first_chan != len_num_chans ))
      return EXIT_FAILURE;

    // In the loop below, ii is the counter for moving through the
    // num_groups vector; the position in the num_groups vector is
    // always the same as the position in the energy bounds vectors
    // xlo and xhi.

    // The value at num_groups[ii] is the number of channel groups
    // associated with the energy bin ii.  Thus we use group_counter
    // to indicate the current position in the first_chan and
    // num_chans vectors, as we can have multiple groups associated
    // with energy bin ii.

    register IndexType ii;

    // The position in the response array is yet different; if there
    // are two channel groups associated with energy bin ii, and the
    // number of channels in the first group is five, and the number
    // of channels in the second group is fifteen: then to advance to
    // the next energy bin, we must:
    //
    //  - Increment ii by one;
    //  - Increment group_counter by two;
    //  - And increment response_counter by twenty (five + fifteen)

    //register IndexType resp_counter = 0;

    // How many groups are in the current energy bin?
    register IndexType current_num_groups = 0;

    // How many channels are in the current group?
    register IndexType current_num_chans = 0;

    // What is the current channel of the output (counts) array?
    //register IndexType current_chan = 0;

    register FloatType source_bin_ii;
    register const ConstFloatType *resp_tmp = resp;
    register const ConstIntType *first_chan_tmp = first_chan;
    register const ConstIntType *num_chans_tmp = num_chans;
    register FloatType *counts_tmp = counts;

    for ( ii = 0; ii < len_source; ii++ ) {
      
      // ii is the current energy bin
      source_bin_ii = source[ ii ];

      current_num_groups = num_groups[ ii ];

      while( current_num_groups ) {
	
	if ( ( IndexType(first_chan_tmp - first_chan) >= len_num_chans ) ||
	     ( UIndexType(*first_chan_tmp) < offset ) )
	  return EXIT_FAILURE;
	counts_tmp = counts + *first_chan_tmp - offset;
	current_num_chans = *num_chans_tmp;
	first_chan_tmp++;
	num_chans_tmp++;
	
	if ( ( (IndexType(counts_tmp-counts) + current_num_chans) > len_counts )
	     ||
	     ( (IndexType(resp_tmp-resp) + current_num_chans) > len_response ) )
	  return EXIT_FAILURE;
	
	while ( current_num_chans ) {

	  *counts_tmp += *resp_tmp * source_bin_ii;
	  counts_tmp++;
	  resp_tmp++;
	  current_num_chans--;

	}
	current_num_groups--;

      }

    } // end for ii
    
    return EXIT_SUCCESS;

  }

  template <typename ConstIntType, typename IndexType, typename IntType>
  bool is_in( const ConstIntType *noticed_chans, IndexType& size,
	      IntType& lo, IntType& hi ) {

    IntType lochan = noticed_chans[0],
            hichan = noticed_chans[size-1];
    
    // assumes that noticed_chans is sorted
    if( lo < lochan && hi > hichan )
      return true;
    
    if( binary_search( noticed_chans, noticed_chans + size, lo ) )
      return true;
    
    if( binary_search( noticed_chans, noticed_chans + size, hi ) )
      return true;
    
    // consider the case of fragmented channels
    if( lochan < lo && lo < hichan && hi > hichan )
      return true;
    
    if( lochan < hi && hi < hichan && lo < lochan )
      return true;

    // consider the 'hidden' fragmented case
    // noticed_chans = 1-249, 400-500, 601-1024
    // lo=250 hi=600
    if( lochan < lo && lo < hichan && lochan < hi && hi < hichan ) {
      const ConstIntType *idx = upper_bound( noticed_chans,
					  noticed_chans + size, lo );
      if( idx == (noticed_chans + size) )
	return false;
      
      // determine if any 'middle' channels are noticed in hidden interval
      if( *idx < hi )
	return true;
    }
    
    return false;
    
//     for( IndexType ii = 0; ii < size; ++ii ) {
//       IntType chan = noticed_chans[ ii ];
//       if( lo <= chan && chan <= hi )
// 	return true;
//     }
//     return false;
  }

  template <typename ConstIntType, typename IntType,
	    typename ConstFloatType, typename FloatType,
	    typename BoolType, typename IndexType>
  int _filter_resp(const ConstIntType *noticed_chans,
		   IndexType len_not_chans,
		   const ConstIntType *n_grp, IndexType len_num_groups,
		   const ConstIntType *f_chan, IndexType len_num_chans,
		   const ConstIntType *n_chan,
		   const ConstFloatType *matrix, IndexType len_response,
		   unsigned int offset,
		   vector<IntType>& grp,
		   vector<IntType>& fch,
		   vector<IntType>& nch,
		   vector<FloatType>& rsp,
		   BoolType *mask) {

    register IndexType response_counter = 0, group_counter = 0;
    register IntType current_num_chans = 0, current_chan = 0,
      current_num_groups = 0, lo = 0, hi = 0;
    
    for( IndexType ii = 0; ii < len_num_groups; ++ii ) {
      
      current_num_groups = n_grp[ ii ];
     
      register IntType ngrp = 0;
 
      while( current_num_groups > 0 ) {
	
	if ( ( group_counter >= len_num_chans ) ||
	     ( f_chan[ group_counter ] < (IntType)offset ) )
	  return EXIT_FAILURE;
	
	current_num_chans = n_chan[ group_counter ];
	current_chan = f_chan[ group_counter ] ;
	
	if ( response_counter + (IndexType)current_num_chans > len_response )
	  return EXIT_FAILURE;

	lo = current_chan;
	// if f_chan values are indices, convert to channels
	if(!offset)
	  lo++;
	hi = lo + current_num_chans;
	
	if( is_in( noticed_chans, len_not_chans, lo, hi ) ) {
	  fch.push_back( current_chan );
	  nch.push_back( current_num_chans );
	  mask[ ii ] = true;
	  ngrp++;

	  while ( current_num_chans > 0 ) {
	    rsp.push_back(matrix[ response_counter ]);
	    response_counter++;
	    current_num_chans--;
	  }
	}
	else {
	  response_counter += current_num_chans;
	}

	group_counter++;
	current_num_groups--;
	
      } // end while
      
      if( ngrp > 0 )
	grp.push_back( ngrp );
      
    } // end for ii
    
    return EXIT_SUCCESS;
  }

  template <typename ConstFloatArrayType, typename IndexType>
  void _sum(const ConstFloatArrayType& data, IndexType start,
	    IndexType stop, SherpaFloat& val) {
    
    val = 0.0;
    for( IndexType ii = start; ii < stop; ii++ )
      val += data[ii];
  }

  template <typename ConstFloatArrayType, typename IndexType>
  void _sum_sq(const ConstFloatArrayType& data, IndexType start,
	       IndexType stop, SherpaFloat& val) {
    
    val = 0.0;
    for( IndexType ii = start; ii < stop; ii++ )
      val += ( data[ii] * data[ii] );
    
    val = sqrt( val );
  }

  template <typename ConstFloatArrayType, typename IndexType>  
  void _max(const ConstFloatArrayType& data, IndexType start,
		 IndexType stop, SherpaFloat& val) {

    SherpaFloat max = data[start];
    for( IndexType ii = start; ii < stop - 1; ii++ )
      max = std::max( max, data[ ii + 1 ] );
    
    val = max;
  }

  template <typename ConstFloatArrayType, typename IndexType>
  void _min(const ConstFloatArrayType& data, IndexType start,
		 IndexType stop, SherpaFloat& val) {

    SherpaFloat min = data[start];
    for( IndexType ii = start; ii < stop - 1; ii++ )
      min = std::min( min, data[ ii + 1 ] );
    
    val = min;
  }

  template <typename ConstFloatArrayType, typename IndexType>  
  void _middle(const ConstFloatArrayType& data, IndexType start,
	       IndexType stop, SherpaFloat& val) {
    
    SherpaFloat min = data[start];
    SherpaFloat max = data[start];
    for( IndexType ii = start; ii < stop - 1; ii++ ) {
      min = std::min( min, data[ ii + 1 ] );
      max = std::max( max, data[ ii + 1 ] );
    }
    val = MID( min, max );
  }
  
  template <typename ConstFloatArrayType, typename FloatArrayType,
	    typename ConstIntArrayType, typename IndexType>
  int _do_group( IndexType len_data, const ConstFloatArrayType& data,
		 IndexType len_group, const ConstIntArrayType& group,
		 FloatArrayType& grouped, const char *type )
  {

    typedef void (*fptr)( const ConstFloatArrayType&, IndexType, IndexType,
			  SherpaFloat&);
    string func(type);
    map<string, fptr> funcs;
    SherpaFloat val;
    
    funcs["sum"] = _sum;
    funcs["_sum_sq"] = _sum_sq;
    funcs["_max"] = _max;
    funcs["_min"] = _min;
    funcs["_middle"] = _middle;

    vector< IndexType > pick_pts;
    
    for( IndexType ii = 0; ii < len_group; ii++ )
      //if( group[ ii ] == 1 )
      // include channels where grouping == 0 so the filter will catch large
      // energy bins
      if( group[ ii ] >= 0 )
	pick_pts.push_back( ii );
    pick_pts.push_back( len_group );

    npy_intp dim = npy_intp( pick_pts.size( ) - 1 );
    if ( EXIT_SUCCESS != grouped.create( 1, &dim ) )
      return EXIT_FAILURE;
    
    for( size_t ii = 0; ii < pick_pts.size( ) - 1; ii++ ) {
      IndexType start = pick_pts[ ii ];
      IndexType stop = pick_pts[ ii + 1 ];
      
      if ( stop > len_data )
	return EXIT_FAILURE;

      if ( func == "_make_groups" ) {
	grouped[ ii ] = data[0] + (SherpaFloat) ii;
	continue;
      }
      
      funcs[func]( data, start, stop, val );
      grouped[ ii ] = val;
      
    } // end ii
    
    return EXIT_SUCCESS;
    
  }


  double interpolate( double x, double x0, double x1,
		      double y0, double y1, double tol ) {
    double m;
    
    if( 0 == _sao_fcmp( x1, x0, tol ) )
      m = 0.0;
    else
      m = ( y1 - y0 ) / ( x1 - x0 );

    if( 0 == _sao_fcmp( x, x0, tol ) &&
	0 == _sao_fcmp( x, x1, tol ) )
      return ((y0 + y1)/2.0);
    
    if( 0 == _sao_fcmp( x, x0, tol ) )
      return y0;
    
    if( 0 == _sao_fcmp( x, x1, tol ) )
      return y1;
    
    else 
      return (y0 + (x-x0)*m);
  }
    
  template <typename FloatArrayType, typename IndexType>
  int _shrink_specresp( FloatArrayType& specresp, FloatArrayType& arf_lo,
			IndexType len_arf,
			FloatArrayType& rmf_lo, FloatArrayType& result,
			IndexType len_rmf ) {
    
    int ii = 0, jj;
    double tol = DBL_EPSILON;

    for( jj = 0; jj < len_rmf; ++jj ) {

      switch ( _sao_fcmp( rmf_lo[ jj ], arf_lo[ ii ], tol ) ) {
	
      case 0:
	// case where rmf_lo[ jj ] == arf_lo[ ii ]
	result[jj] = specresp[ii];
	ii++;

	if( ii == len_arf )
	  return EXIT_SUCCESS;
	
	break;

      case 1:
	// case where rmf_lo[ jj ] > arf_lo[ ii ]
	int res;
	while( (res = _sao_fcmp( rmf_lo[ jj ], arf_lo[ ii ], tol )) == 1 ) {
	  ii++;	  
	  if( ii == len_arf )
	    return EXIT_SUCCESS;
	}

	// case where rmf_lo[ jj ] == arf_lo[ ii ]
	if( res == 0 ) {
	  result[jj] = specresp[ii];
	  ii++;
	  if( ii == len_arf )
	    return EXIT_SUCCESS;
	}

	// case where arf_lo[ ii - 1 ] < rmf_lo[ jj ] < arf_lo[ ii ] and ii > 0
	else if( ii > 0 ) {
	  result[jj] = interpolate( rmf_lo[ jj ], arf_lo[ ii - 1 ],
				    arf_lo[ ii ], specresp[ ii - 1 ],
				    specresp[ ii ], tol );

	  //result[jj] = specresp[ ii - 1 ];
	  //result[jj] = specresp[ ii ];
	  ii++;
	  if( ii == len_arf )
	    return EXIT_SUCCESS;
	}
	else {
	  return EXIT_FAILURE;
	}
	break;	
	
      default:
      	// rmf_lo[ jj ] < arf_lo[ ii ] or _sao_fcmp failed
	return EXIT_FAILURE;
	
      }
    }
    
    return EXIT_SUCCESS;
  } 
  

}  }  } /* namespace utils, namespace astro, namespace sherpa */


#endif /* __sherpa_astro_utils_hh__ */
