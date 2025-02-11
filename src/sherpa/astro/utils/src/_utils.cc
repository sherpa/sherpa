//
//  Copyright (C) 2009, 2017, 2021, 2022, 2024
//  Smithsonian Astrophysical Observatory
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

#include "sherpa/extension.hh"
#include "sherpa/astro/utils.hh"
#include <sstream>
#include <iostream>
#include <stdexcept>

extern "C" {
  void init_utils();
}

typedef sherpa::Array< npy_bool, NPY_BOOL > BoolArray;

namespace sherpa { namespace astro { namespace utils {

  template <typename ArrayType>
  PyObject* arf_fold( PyObject* self, PyObject* args )
  {

    ArrayType source;
    ArrayType effarea;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&",
			    (converter)convert_to_array< ArrayType >,
			    &source,
			    (converter)convert_to_array< ArrayType >,
			    &effarea ) )
      return NULL;

    npy_intp nelem = source.get_size();

    if ( effarea.get_size() != nelem ) {
      ostringstream err;
      err << "input array sizes do not match, "
	  << "source: " << nelem << " vs effarea: " << effarea.get_size();
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    ArrayType result;
    if ( EXIT_SUCCESS != result.create( source.get_ndim(),
					source.get_dims() ) )
      return NULL;

    arf_fold( nelem, effarea, source, result );

    return result.return_new_ref();

  }


  template <typename FloatArrayType, typename IntArrayType>
  PyObject* rmf_fold( PyObject* self, PyObject* args )
  {

    FloatArrayType source;
    IntArrayType num_groups;
    IntArrayType first_chan;
    IntArrayType num_chans;
    FloatArrayType response;
    long len_counts;
    unsigned int offset;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&lI",
			    (converter)convert_to_array< FloatArrayType >,
			    &source,
			    (converter)convert_to_array< IntArrayType >,
			    &num_groups,
			    (converter)convert_to_array< IntArrayType >,
			    &first_chan,
			    (converter)convert_to_array< IntArrayType >,
			    &num_chans,
			    (converter)convert_to_array< FloatArrayType >,
			    &response,
			    &len_counts,
			    &offset) )
      return NULL;

    //
    // The rmf_fold function will validate argument sizes
    //

    npy_intp dim = npy_intp( len_counts );
    FloatArrayType counts;
    if ( EXIT_SUCCESS != counts.zeros( 1, &dim ) )
      return NULL;

    if ( EXIT_SUCCESS != rmf_fold( source.get_size(), &source[0],
				   num_groups.get_size(), &num_groups[0],
				   first_chan.get_size(), &first_chan[0],
				   num_chans.get_size(), &num_chans[0],
				   response.get_size(), &response[0],
				   counts.get_size(), &counts[0],
				   npy_uintp(offset))) {

      PyErr_SetString( PyExc_ValueError,
		       (char*)"RMF data is invalid or inconsistent" );
      return NULL;

    }

    return counts.return_new_ref();

  }

  template <typename FloatArrayType, typename IntArrayType>
  PyObject* do_group( PyObject* self, PyObject* args )
  {

    FloatArrayType data;
    FloatArrayType grouped;
    IntArrayType group;
    const char* type;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&s",
			    (converter)convert_to_array< FloatArrayType >,
			    &data,
			    (converter)convert_to_array< IntArrayType >,
			    &group,
			    &type))

    return NULL;

    if ( data.get_size() != group.get_size() ) {
      ostringstream err;
      err << "input array sizes do not match, "
	  << "data: " << data.get_size() << " vs group: " << group.get_size();
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    try {
      if( EXIT_SUCCESS != _do_group(data.get_size(), data,
				    group.get_size(), group,
				    grouped, type) ) {
	PyErr_SetString( PyExc_ValueError,
			 (char*)"group data is invalid or inconsistent" );
	return NULL;
      }
    } catch ( std::out_of_range& ) {
      ostringstream err;
      err << "unsupported group function: " << type;
      PyErr_SetString( PyExc_ValueError, err.str().c_str() );
      return NULL;
    }

  return grouped.return_new_ref();

  }


  template <typename FloatArrayType>
  PyObject* shrink_specresp( PyObject* self, PyObject* args )
  {

    FloatArrayType specresp;
    FloatArrayType arf_lo;
    FloatArrayType rmf_lo;
    FloatArrayType result;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&O&",
			    (converter)convert_to_array< FloatArrayType >,
			    &specresp,
			    (converter)convert_to_array< FloatArrayType >,
			    &arf_lo,
			    (converter)convert_to_array< FloatArrayType >,
			    &rmf_lo))

    return NULL;

    if ( specresp.get_size() != arf_lo.get_size() ) {
      ostringstream err;
      err << "input array sizes do not match, "
	  << "specresp: " << specresp.get_size()
	  << " vs arf_lo: " << arf_lo.get_size();
      PyErr_SetString( PyExc_TypeError, err.str().c_str() );
      return NULL;
    }

    if ( rmf_lo.get_size() > specresp.get_size() ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"RMF is higher resolution than ARF.  Need to expand, not shrink effective area" );
      return NULL;
    }

    if ( EXIT_SUCCESS != result.create( rmf_lo.get_ndim(),
					rmf_lo.get_dims() ) )
      return NULL;

    if( EXIT_SUCCESS != _shrink_specresp( specresp,
					  arf_lo,
					  arf_lo.get_size(),
					  rmf_lo,
					  result,
					  rmf_lo.get_size() ) ) {
      PyErr_SetString( PyExc_ValueError,
		       (char*)"shrinking effective area failed" );
      return NULL;
    }

  return result.return_new_ref();

  }

  template <typename FloatArrayType, typename IntArrayType,
	    typename IntType, typename FloatType>
  PyObject* filter_resp(PyObject* self, PyObject* args) {

    IntArrayType noticed_chans;
    IntArrayType n_grp;
    IntArrayType f_chan;
    IntArrayType n_chan;
    FloatArrayType matrix;

    BoolArray mask;
    vector<FloatType> resp_buf;
    vector<IntType> f_buf;
    vector<IntType> n_buf;
    vector<IntType> grp_buf;

    FloatArrayType resp;
    IntArrayType grp;
    IntArrayType fchan;
    IntArrayType nchan;

    npy_intp gdims[1];
    npy_intp rdims[1];
    npy_intp cdims[1];

    unsigned int offset;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&O&O&O&I",
			    (converter)convert_to_contig_array< IntArrayType >,
			    &noticed_chans,
			    (converter)convert_to_array< IntArrayType >,
			    &n_grp,
			    (converter)convert_to_array< IntArrayType >,
			    &f_chan,
			    (converter)convert_to_array< IntArrayType >,
			    &n_chan,
			    (converter)convert_to_array< FloatArrayType >,
			    &matrix,
			    &offset) )
      return NULL;

    if (0 == noticed_chans.get_size()) {
      PyErr_SetString( PyExc_ValueError,
		       (char*)"There are no noticed channels" );
      return NULL;
    }

    if( EXIT_SUCCESS != mask.zeros( 1, n_grp.get_dims() ))
      return NULL;

    grp_buf.reserve(size_t(n_grp.get_size()));
    f_buf.reserve(size_t(f_chan.get_size()));
    n_buf.reserve(size_t(n_chan.get_size()));
    resp_buf.reserve(size_t(matrix.get_size()));

    if( EXIT_SUCCESS != _filter_resp( &noticed_chans[0],
				      noticed_chans.get_size(),
				      &n_grp[0], n_grp.get_size(), &f_chan[0],
				      f_chan.get_size(), &n_chan[0], &matrix[0],
				      matrix.get_size(), offset, grp_buf, f_buf,
				      n_buf, resp_buf, &mask[0]) ) {
      PyErr_SetString( PyExc_TypeError,
		       (char*)"response filter failed" );
      return NULL;
    }

    gdims[0] = npy_intp(grp_buf.size());
    rdims[0] = npy_intp(resp_buf.size());
    cdims[0] = npy_intp(f_buf.size());

    if( EXIT_SUCCESS != resp.create( 1, rdims ))
      return NULL;

    if( EXIT_SUCCESS != fchan.create( 1, cdims ))
      return NULL;

    if( EXIT_SUCCESS != nchan.create( 1, cdims ))
      return NULL;

    if( EXIT_SUCCESS != grp.create( 1, gdims ))
      return NULL;

    for(int ii = 0; ii < gdims[0]; ++ii) {
      grp[ii] = grp_buf[ii];
    }

    for(int ii = 0; ii < cdims[0]; ++ii) {
      fchan[ii] = f_buf[ii];
      nchan[ii] = n_buf[ii];
    }

    for(int ii = 0; ii < rdims[0]; ++ii)
      resp[ii] = resp_buf[ii];

    return Py_BuildValue( (char*)"NNNNN",
			 grp.return_new_ref(),
			 fchan.return_new_ref(),
			 nchan.return_new_ref(),
			 resp.return_new_ref(),
			 mask.return_new_ref());

  }


  static PyObject* _expand_grouped_mask(PyObject* self, PyObject* args) {

    BoolArray mask;
    BoolArray res;
    IntArray group;

    if ( !PyArg_ParseTuple( args, (char*)"O&O&",
			    (converter)convert_to_array< BoolArray >,
			    &mask,
			    (converter)convert_to_array< IntArray >,
			    &group ) )
      return NULL;

    if (mask.get_ndim() != 1) {
      PyErr_SetString( PyExc_ValueError, (char*)"mask array must be 1D" );
      return NULL;
    }

    if (group.get_ndim() != 1) {
      PyErr_SetString( PyExc_ValueError, (char*)"group array must be 1D" );
      return NULL;
    }

    // Use TypeError for the size=0 checks as this is what has historically been
    // used for this check (when only mask_size was checked).
    //
    npy_intp mask_size = mask.get_size();
    npy_intp group_size = group.get_size();
    if( mask_size == 0 ) {
      PyErr_SetString( PyExc_TypeError, (char*)"mask array has no elements" );
      return NULL;
    }

    if( group_size == 0 ) {
      PyErr_SetString( PyExc_TypeError, (char*)"group array has no elements" );
      return NULL;
    }

    if (group[0] < 0) {
      PyErr_SetString( PyExc_ValueError, (char*)"The first element of group is negative" );
      return NULL;
    }

    if( EXIT_SUCCESS != res.zeros( group.get_ndim(), group.get_dims() ))
      return NULL;

    npy_intp jj = 0;

    // The group array starts a group with a value >= 0 and continues the
    // group with a value < 0. The first element therefore starts a new
    // group.
    //
    // The assumption is that the difference, in time, to just setting
    //     res[ii] = mask[jj]
    // is insignificant.
    //
    if( mask[0] )
      res[0] = true;

    for( npy_intp ii = 1; ii < group_size; ++ii ) {
      // Are we starting a new group?
      //
      if( group[ ii ] >= 0 ) {
	++jj;
	if (jj >= mask_size) {
	  PyErr_SetString( PyExc_ValueError, (char*) "More groups than mask elements" );
	  return NULL;
	}
      }
      if( mask[jj] )
	res[ii] = true;
    }

    if (jj != (mask_size - 1)) {
      PyErr_SetString( PyExc_ValueError, (char*) "More mask elements than groups" );
      return NULL;
    }

    return Py_BuildValue( (char*)"N", res.return_new_ref() );
  }

  template <typename ArrayType, typename DataType>
  PyObject* is_in(PyObject* self, PyObject* args) {

    ArrayType chans;
    DataType lo;
    DataType hi;
    int size;
    bool val;

    if ( !PyArg_ParseTuple( args, (char*)"O&II",
			    (converter)convert_to_contig_array< ArrayType >,
			    &chans, &lo, &hi) )
      return NULL;

    size = int(chans.get_size());
    val = is_in(&chans[0], size, lo, hi);

    return Py_BuildValue( (char*)"O", PyBool_FromLong(long(val)) );
  }

}  }  } /* namespace utils, namespace astro, namespace sherpa */

static PyMethodDef UtilsFcts[] = {

  FCTSPEC( arf_fold, sherpa::astro::utils::arf_fold< SherpaFloatArray > ),
  FCTSPEC( rmf_fold, (sherpa::astro::utils::rmf_fold< SherpaFloatArray,
		      SherpaUIntArray >) ),

  FCTSPECDOC( do_group, (sherpa::astro::utils::do_group<SherpaFloatArray, IntArray>),
	      "Group the array using OGIP standards.\n\n"
	      "Parameters\n"
	      "----------\n"
	      "data : array_like\n"
	      "    The data to group.\n"
	      "group : array_like\n"
	      "    The OGIP grouping data: 1 indicates the start of a group and\n"
	      "    -1 continues the group.\n"
	      "name : {'sum', '_sum_sq', '_max', '_min', '_middle', '_make_groups'}\n"
	      "    The grouping scheme to combine values within a group.\n\n"
	      "Returns\n"
	      "-------\n"
	      "grouped : array\n"
	      "    The grouped data. It will be smaller than data unless group only\n"
	      "    contains 1's.\n\n"
	      "Examples\n"
	      "--------\n\n"
	      "Group the array [1, 2, 3, 4, 5, 6] into groups of length 2, 1, and 3,\n"
	      "using different grouping schemes:\n\n"
	      ">>> data = [1, 2, 3, 4, 5, 6]\n"
	      ">>> group = [1, -1, 1, 1, -1, -1]\n"
	      ">>> do_group(data, group, '_make_groups')\n"
	      "[1, 2, 3]\n"
	      ">>> do_group(data, group, 'sum')\n"
	      "[3, 3, 15]\n"
	      ">>> do_group(data, group, '_min')\n"
	      "[1, 3, 4]\n"
	      ">>> do_group(data, group, '_max')\n"
	      "[2, 3, 6]\n"
	      ">>> do_group(data, group, '_middle')\n"
	      "[1.5, 3. , 5. ]\n" ),

  FCTSPEC( shrink_effarea, (sherpa::astro::utils::shrink_specresp<
			    SherpaFloatArray >) ),

  FCTSPEC( filter_resp, (sherpa::astro::utils::filter_resp< SherpaFloatArray,
			 SherpaUIntArray, SherpaUInt, SherpaFloat >)),

  FCTSPECDOC( expand_grouped_mask, sherpa::astro::utils::_expand_grouped_mask,
	      "Expand a mask array to match the ungrouped data.\n\n"
	      "The mask array size must match the number of groups in the\n"
	      "group array.\n\n"
	      "Parameters\n"
	      "----------\n"
	      "mask : array_like\n"
	      "    The mask array to expand (treated as a boolean array).\n"
	      "group : array_like\n"
	      "    The OGIP grouping data: 1 indicates the start of a group and\n"
	      "    -1 continues the group.\n\n"
	      "Returns\n"
	      "-------\n"
	      "full_mask : array\n"
	      "    The mask data (booleans) expanded to match the size of group.\n\n"
	      "Examples\n"
	      "--------\n\n"
	      ">>> mask = [True, False, False, True]\n"
	      ">>> group = [1, -1, 1, 1, -1, 1, -1, -1]\n"
	      ">>> expand_grouped_mask(mask, group)\n"
	      "[True, True, False, False, False, True, True, True]\n\n"
              ">>> expand_grouped_mask([True, False], [1, 1])\n"
              "[True, False]\n"
              ">>> expand_grouped_mask([1, 0], [1, 1])\n"
              "[True, False]\n" ),

  FCTSPEC( is_in, (sherpa::astro::utils::is_in<SherpaUIntArray, SherpaUInt>)),

  { NULL, NULL, 0, NULL }

};


SHERPAMODDOC(_utils, UtilsFcts,
	     "Routines for handling Astronomy data.\n\n"
	     "OGIP documentation:\n"
	     " - \"The OGIP Spectral File Format\", https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007/ogip_92_007.html\n\n"
	     " - \"The OGIP Spectral File Format Addendum: Changes log\", https://heasarc.gsfc.nasa.gov/docs/heasarc/ofwg/docs/spectra/ogip_92_007a/ogip_92_007a.html\n\n"
	     )
