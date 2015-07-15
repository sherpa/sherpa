// 
//  Copyright (C) 2009, 2015  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_astro_xspec_extension_hh__
#define __sherpa_astro_xspec_extension_hh__

// Have sherpa includes first so that Python.h is first, to avoid warning
// messages about redefining _XOPEN_SOURCE
#include <sherpa/extension.hh>
#include <sherpa/constants.hh>
#include <cfloat>
#include <vector>
#include <sstream>
#include <iostream>
#include "sherpa/fcmp.hh"

namespace sherpa { namespace astro { namespace xspec {


typedef sherpa::Array< float, NPY_FLOAT > FloatArray;
typedef float FloatArrayType;

// XSpec models can be called from Sherpa using either
//   - a single array for the grid
//   - two arrays for the grid
//
// The first form is straightforward, since this is what the XSpec
// models expect. In this case a 0 is added to the end of the array
// returned by XSpec since Sherpa expects that the input and output
// arrays have the same length.
//
// When given two arrays, there is the possibility of a non-contiguous
// grid - e.g. if the low values are [0.1,0.2,0.6,0.7] and the high
// values are [0.2,0.3,0.7,0.8] then the range 0.3 to 0.6 is not
// required. This requires passing through the arrays to find any
// gaps, and then dealing with them. The original approach was to
// just create a grid as if there were no gaps; that is, the model would
// be evaluated on the grid [0.1,0.2,0.6,0.7,0.8] but record the position
// and correct widths for the "gaps" and then re-run the model to
// "fill in" the output array. This is being changed so that the grid
// is created with extra edges, so that a single call is made, and then
// the un-wanted bins are removed. This is in part because XSpec models
// may require a large set-up time (so it's possibly an optimisation),
// but also because some models involve interpolating the answer onto
// the output grid, so that unexpected results may be returned for
// bins at the end of a contiguous section (the new scheme is not
// guaranteed to fix this, but is likely to be better). Another
// alternative would be to run the model on each contiguous section
// (this has not been tried).
//
// The convention is that if the grids are in ascending order they
// are in keV - the units required by XSpec - otherwise they are
// in Angstrom, so must be converted to keV. Note that the constraint
// xhi > xlo is expected to hold even when the units are Angstrom -
// e.g. xlo = 112.7, 103.3, 95.4, ...
//      xhi = 124.0, 112.7, 103.3, ...
//

// When creating the flux array to be sent to XSpec, the array is filled
// in with 0's (that is, the zeros array method is used, rather than
// create). This is to ensure that a "sensible" answer is returned on
// error (all 0's), as the XSpec model API does not provide a way to
// return an error status. There is (likely; not checked) a run-time cost
// to using zeros rather than create, but safety is better than
// performance here. As a micro-optimization, the error array is
// now just created as a C++ vector, rather than a NumPy array,
// and so is no-longer set to 0.

// The spectrum number is set to 1. It used to be 0, but some code
// (a user model, so not included in the XSpec model library being
// built against) has been seen to behave strangely with a value of 0,
// so a value of 1 is being used "just in case". A keyword argument
// could be added so that the user can override this, but it is
// only really woth doing once write access to the XFLT keywords
// is added (and then working out how to take advantage of it; perhaps
// a wrapper model that provides a parameter-like interface to the value
// so that a model instance can be associated with a particular dataset).
// This might be enough to then support the smaug model.

template <npy_intp NumPars, bool HasNorm,
void (*XSpecFunc)( float* ear, int* ne, float* param, int* ifl,
		float* photar, float* photer )>
PyObject* xspecmodelfct( PyObject* self, PyObject* args )
{

#ifdef INIT_XSPEC
	if ( EXIT_SUCCESS != INIT_XSPEC() )
		return NULL;
#endif

        FloatArray pars;
        FloatArray xlo;
        FloatArray xhi;

        // The code used to read in xlo and xhi as DoubleArray, which
        // would then be downcast to FloatArray. It seems to make sense
        // to do the downcasting here.
        if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
			(converter)convert_to_contig_array< FloatArray >,
			&pars,
			(converter)convert_to_contig_array< FloatArray >,
			&xlo,
			(converter)convert_to_contig_array< FloatArray >,
			&xhi ) )
		return NULL;

	npy_intp npars = pars.get_size();

	if ( NumPars != npars ) {
		std::ostringstream err;
		err << "expected " << NumPars << " parameters, got " << npars;
		PyErr_SetString( PyExc_TypeError, err.str().c_str() );
		return NULL;
	}

	int nelem = int( xlo.get_size() );

        if ( nelem < 2 ) {
          std::ostringstream err;
          err << "input array must have at least 2 elements, found " << nelem;
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        if( xhi && (nelem != int(xhi.get_size())) ) {
          std::ostringstream err;
          err << "input arrays are not the same size: " << nelem
              << " and " << int( xhi.get_size() );
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

	int ifl = 1;

	bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;
        FloatArray *x1 = &xlo;
        FloatArray *x2 = &xhi;
        if (is_wave && xhi) {
            x1 = &xhi;
            x2 = &xlo;
        }

        // Are there any non-contiguous bins?
        std::vector<int> gaps_index;
        std::vector<FloatArrayType> gaps_edges;
        if (xhi) {
          const int gap_found = is_wave ? 1 : -1;
          for (int i = 0; i < nelem-1; i++) {
            int cmp = sao_fcmp((*x2)[i], (*x1)[i+1], FLT_EPSILON);
            if (cmp == gap_found) {
              gaps_index.push_back(i);
              gaps_edges.push_back((*x2)[i]);
            } else if (cmp != 0) {
              std::ostringstream err;
              // not convinced this is understandable to users, particularly
              // if the grid is in Angstrom. It is also possible that the
              // format used isn't sufficient to show the problem, but I
              // do not want to tweak the format here just yet.
              err << "Grid cells overlap: cell " << i
                  << " (" << (*x1)[i] << " to " << (*x2)[i] << ")"
                  << " and cell " << (i+1)
                  << " (" << (*x1)[i+1] << " to " << (*x2)[i+1] << ")";
              PyErr_SetString( PyExc_ValueError, err.str().c_str() );
              return NULL;
            }
          }
        }

        int ngaps = (int) gaps_edges.size();

        // The size of the energy array sent to XSpec
        int ngrid = nelem;
        if (xhi) {
          ngrid += 1 + ngaps;
        }

        // XSpec traditionally refers to the input energy grid as ear.
        std::vector<FloatArrayType> ear(ngrid);

        // The grid is created, converted from Angstrom to Energy
        // (if required), and then checked for being monotonic.
        // The multiple loops are not necessarily as efficient
        // as a single loop, but simpler to write.
        //
        {
          // Process the contiguous sections by looping through
          // the gaps_index/edges arrays.
          int start = 0;
          for (int j = 0 ; j < ngaps; j++) {
            int end = gaps_index[j] + 1;
            for(int i = start; i < end; i++) {
              ear[i + j] = (*x1)[i];
            }
            ear[end + j] = gaps_edges[j];
            start = end;
          }

          // need to do the last contiguous grid
          for(int i = start; i < nelem; i++) {
            ear[i + ngaps] = (*x1)[i];
          }

          // Add on the last bin value if needed
          if (xhi) {
            ear[ngrid - 1] = (*x2)[nelem - 1];
          }
        }

        // Originally the conversion was done using the Double input
        // array. It's now done using the float array sent to the
        // model. Is this likely to be a source of error or confusion
        // (the numeric result may be slightly different than previously
        // die to the different order of the casting)?
        //
        if (is_wave) {
          float hc = (float) (sherpa::constants::c_ang<SherpaFloat>() *
                              sherpa::constants::h_kev<SherpaFloat>());
          for (int i = 0; i < ngrid; i++) {
            if (ear[i] <= 0.0) {
              std::ostringstream err;
              err << "Wavelength must be > 0, sent " << ear[i];
              PyErr_SetString( PyExc_ValueError, err.str().c_str() );
              return NULL;
            }
            ear[i] = hc / ear[i];
          }
        }

        // Check for monotonic (could be included in the above, but
        // this is much simpler to write here).
        // Should this be done, or just let the user get invalid
        // results?
        //
        // The earlier check with sao_fcmp catches some of these,
        // but not all of them (if xhi is not given, or if xlo==xhi
        // for any bin).
        //
        for (int i = 0; i < ngrid - 1; i++) {
          if (ear[i] >= ear[i+1]) {
            std::ostringstream err;
            err << "Grid is not monotonic: " << ear[i] << " to " <<
              ear[i+1];
            PyErr_SetString( PyExc_ValueError, err.str().c_str() );
            return NULL;
          }
        }

        // Although the XSpec model expects the flux/fluxerror arrays
        // to have size ngrid-1, the return array has to match the
        // input size.
        npy_intp dims[1] = { ngrid };
        if (xhi)
          dims[0]--;

        FloatArray result;
        if ( EXIT_SUCCESS != result.zeros( 1, dims ) )
          return NULL;

        // Since the flux error is discarded, it does not need to be a
        // NumPy array. Ideally it would be set to zeros for safety.
        std::vector<FloatArrayType> error(dims[0]);

	// Even though the XSPEC model function is Fortran, it could call
	// C++ functions, so swallow exceptions here

	try {

          int npts = ngrid - 1;
          XSpecFunc( &ear[0], &npts, &pars[0], &ifl, &result[0], &error[0] );

	} catch(...) {

          PyErr_SetString( PyExc_ValueError,
                           (char*)"XSPEC model evaluation failed" );
          return NULL;

	}

        // Remove gaps
        if (ngaps > 0) {
          // We can skip copying the first contiguous block
          // as it is the identity transform.
          //
          int start = gaps_index[0] + 1;
          for (int j = 1 ; j < ngaps; j++) {
            int end = gaps_index[j] + 1;
            for(int i = start; i < end; i++) {
              result[i] = result[i + j];
            }
            start = end;
          }

          // need to do the last contiguous grid
          for(int i = start; i < nelem; i++) {
            result[i] = result[i + ngaps];
          }

          // Resize the data.
          result.resize1d(nelem);
        }

	// Apply normalization if required
	if ( HasNorm )
          for (int i = 0; i < nelem; i++)
            result[i] *= pars[NumPars - 1];

	return result.return_new_ref();

}


template <npy_intp NumPars, bool HasNorm,
void (*XSpecFunc)( const double* energy, int nFlux,
		const double* params, int spectrumNumber,
		double* flux, double* fluxError,
		const char* initStr )>
PyObject* xspecmodelfct_C( PyObject* self, PyObject* args )
{

#ifdef INIT_XSPEC
	if ( EXIT_SUCCESS != INIT_XSPEC() )
		return NULL;
#endif

	DoubleArray pars;
	DoubleArray xlo;
	DoubleArray xhi;

	if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
			(converter)convert_to_contig_array< DoubleArray >,
			&pars,
			(converter)convert_to_contig_array< DoubleArray >,
			&xlo,
			(converter)convert_to_contig_array< DoubleArray >,
			&xhi ) )
		return NULL;

	npy_intp npars = pars.get_size();

	if ( NumPars != npars ) {
		std::ostringstream err;
		err << "expected " << NumPars << " parameters, got " << npars;
		PyErr_SetString( PyExc_TypeError, err.str().c_str() );
		return NULL;
	}

	int nelem = int( xlo.get_size() );

        if ( nelem < 2 ) {
          std::ostringstream err;
          err << "input array must have at least 2 elements, found " << nelem;
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        if( xhi && (nelem != int(xhi.get_size())) ) {
          std::ostringstream err;
          err << "input arrays are not the same size: " << nelem
              << " and " << int( xhi.get_size() );
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        int ifl = 1;

        bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;
        DoubleArray *x1 = &xlo;
        DoubleArray *x2 = &xhi;
        if (is_wave && xhi) {
            x1 = &xhi;
            x2 = &xlo;
        }

        // Are there any non-contiguous bins?
        std::vector<int> gaps_index;
        std::vector<SherpaFloat> gaps_edges;
        if (xhi) {
          const int gap_found = is_wave ? 1 : -1;
          for (int i = 0; i < nelem-1; i++) {
            int cmp = sao_fcmp((*x2)[i], (*x1)[i+1], DBL_EPSILON);
            if (cmp == gap_found) {
              gaps_index.push_back(i);
              gaps_edges.push_back((*x2)[i]);
            } else if (cmp != 0) {
              std::ostringstream err;
              // not convinced this is understandable to users, particularly
              // if the grid is in Angstrom. It is also possible that the
              // format used isn't sufficient to show the problem, but I
              // do not want to tweak the format here just yet.
              err << "Grid cells overlap: cell " << i
                  << " (" << (*x1)[i] << " to " << (*x2)[i] << ")"
                  << " and cell " << (i+1)
                  << " (" << (*x1)[i+1] << " to " << (*x2)[i+1] << ")";
              PyErr_SetString( PyExc_ValueError, err.str().c_str() );
              return NULL;
            }
          }
        }

        int ngaps = (int) gaps_edges.size();

        // The size of the energy array sent to XSpec
        int ngrid = nelem;
        if (xhi) {
          ngrid += 1 + ngaps;
        }

        // XSpec traditionally refers to the input energy grid as ear.
        std::vector<SherpaFloat> ear(ngrid);

        // The grid is created, converted from Angstrom to Energy
        // (if required), and then checked for being monotonic.
        // The multiple loops are not necessarily as efficient
        // as a single loop, but simpler to write.
        //
        {
          // Process the contiguous sections by looping through
          // the gaps_index/edges arrays.
          int start = 0;
          for (int j = 0 ; j < ngaps; j++) {
            int end = gaps_index[j] + 1;
            for(int i = start; i < end; i++) {
              ear[i + j] = (*x1)[i];
            }
            ear[end + j] = gaps_edges[j];
            start = end;
          }

          // need to do the last contiguous grid
          for(int i = start; i < nelem; i++) {
            ear[i + ngaps] = (*x1)[i];
          }

          // Add on the last bin value if needed
          if (xhi) {
            ear[ngrid - 1] = (*x2)[nelem - 1];
          }
        }

        if (is_wave) {
          double hc = (sherpa::constants::c_ang<SherpaFloat>() *
                       sherpa::constants::h_kev<SherpaFloat>());
          for (int i = 0; i < ngrid; i++) {
            if (ear[i] <= 0.0) {
              std::ostringstream err;
              err << "Wavelength must be > 0, sent " << ear[i];
              PyErr_SetString( PyExc_ValueError, err.str().c_str() );
              return NULL;
            }
            ear[i] = hc / ear[i];
          }
        }

        // Check for monotonic (could be included in the above, but
        // this is much simpler to write here).
        // Should this be done, or just let the user get invalid
        // results?
        //
        // The earlier check with sao_fcmp catches some of these,
        // but not all of them (if xhi is not given, or if xlo==xhi
        // for any bin).
        //
        for (int i = 0; i < ngrid - 1; i++) {
          if (ear[i] >= ear[i+1]) {
            std::ostringstream err;
            err << "Grid is not monotonic: " << ear[i] << " to " <<
              ear[i+1];
            PyErr_SetString( PyExc_ValueError, err.str().c_str() );
            return NULL;
          }
        }

        // Although the XSpec model expects the flux/fluxerror arrays
        // to have size ngrid-1, the return array has to match the
        // input size.
        npy_intp dims[1] = { ngrid };
        if (xhi)
          dims[0]--;

	DoubleArray result;
	if ( EXIT_SUCCESS != result.zeros( 1, dims ) )
		return NULL;

        // Since the flux error is discarded, it does not need to be a
        // NumPy array. Should be set to zeros for safety.
        std::vector<SherpaFloat> error(dims[0]);

	try {

          int npts = ngrid - 1;
          XSpecFunc( &ear[0], npts, &pars[0], ifl, &result[0], &error[0], NULL );

	} catch(...) {

          PyErr_SetString( PyExc_ValueError,
                           (char*)"XSPEC model evaluation failed" );
          return NULL;

	}

        // Remove gaps
        if (ngaps > 0) {
          // We can skip copying the first contiguous block
          // as it is the identity transform.
          //
          int start = gaps_index[0] + 1;
          for (int j = 1 ; j < ngaps; j++) {
            int end = gaps_index[j] + 1;
            for(int i = start; i < end; i++) {
              result[i] = result[i + j];
            }
            start = end;
          }

          // need to do the last contiguous grid
          for(int i = start; i < nelem; i++) {
            result[i] = result[i + ngaps];
          }

          // Resize the data.
          result.resize1d(nelem);
        }

	// Apply normalization if required
	if ( HasNorm )
          for (int i = 0; i < nelem; i++)
            result[i] *= pars[NumPars - 1];

	return result.return_new_ref();

}

// For convolution models; assumed to have C++ linkage
//
// Should this just error out if gaps are found, since this is
// unlikely to work well with a convolution-style model,
// whatever scheme is used.
//
template <npy_intp NumPars,
void (*XSpecFunc)( const double* energy, int nFlux,
		const double* params, int spectrumNumber,
		double* flux, double* fluxError,
		const char* initStr )>
PyObject* xspecmodelfct_con( PyObject* self, PyObject* args )
{

#ifdef INIT_XSPEC
	if ( EXIT_SUCCESS != INIT_XSPEC() )
		return NULL;
#endif

	DoubleArray pars;
	DoubleArray xlo;
	DoubleArray xhi;
	DoubleArray fluxes;
	DoubleArray *x;

        // The arguments are parsed as
        //   pars, xlo, fluxes
        //   pars, xlo, xhi, fluxes
        //
	if ( !PyArg_ParseTuple( args, (char*)"O&O&O&|O&",
			(converter)convert_to_contig_array< DoubleArray >,
			&pars,
			(converter)convert_to_contig_array< DoubleArray >,
			&xlo,
			(converter)convert_to_contig_array< DoubleArray >,
			&xhi,
			(converter)convert_to_contig_array< DoubleArray >,
			&fluxes ) )
		return NULL;

	npy_intp npars = pars.get_size();

	if ( NumPars != npars ) {
		std::ostringstream err;
		err << "expected " << NumPars << " parameters, got " << npars;
		PyErr_SetString( PyExc_TypeError, err.str().c_str() );
		return NULL;
	}

	int nelem = int( xlo.get_size() );

        if ( nelem < 2 ) {
          std::ostringstream err;
          err << "input array must have at least 2 elements, found " << nelem;
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        // swap arguments if needed
        if (xhi && !fluxes) {
          DoubleArray *tmp = &fluxes;
          fluxes = xhi;
          xhi = *tmp;
        }

        if( xhi && (nelem != int(xhi.get_size())) ) {
          std::ostringstream err;
          err << "input arrays are not the same size: " << nelem
              << " and " << int( xhi.get_size() );
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        // For now require the fluxes array to have the same
        // size as the input grid. If xhi is not given then
        // technically fluxes should be one less, but this is
        // likely to cause problems (as it doesn't match how
        // the rest of the interface works).
        if( nelem != int(fluxes.get_size()) ) {
          std::ostringstream err;
          err << "flux array does not match the input grid: " << nelem
              << " and " << int( fluxes.get_size() );
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

	double hc = (sherpa::constants::c_ang<SherpaFloat>() *
			sherpa::constants::h_kev<SherpaFloat>());
	bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

	std::vector<int> gaps;
	std::vector<double> gap_widths;

	// The XSPEC functions expect the input array to be of length nFlux+1
	int near = nelem;
	if( xhi ) {
		near++;

                // See comments in xspecmodelfct template.

		for (int i = 0; i < nelem-1; i++) {
			int cmp;
			if ( is_wave ) {
				cmp = sao_fcmp(xlo[i], xhi[i+1], DBL_EPSILON);
			} else {
				cmp = sao_fcmp(xhi[i], xlo[i+1], DBL_EPSILON);
			}
			if (0 != cmp) {
				gaps.push_back(i);
				double width = fabs(xhi[i] - xlo[i]);
				if( is_wave ) {
					width = hc / width;
				}
				gap_widths.push_back(width);
			}
		}
	}

	std::vector<SherpaFloat> ear(near);

	for( int ii = 0; ii < nelem; ii++ ) {
		if( is_wave ) {

			// wave analysis swaps edges, e.g. wave_hi <--> energy_lo
			// if xhi is available use it
			x = (xhi) ? &xhi : &xlo;

			if ( 0.0 == (*x)[ii] ) {
				PyErr_SetString( PyExc_ValueError,
						(char*)"XSPEC model evaluation failed, division by zero" );
				return NULL;
			}
			ear[ ii ] = ( SherpaFloat ) (hc / (*x)[ ii ]);
		}
		else
			ear[ ii ] = ( SherpaFloat ) xlo[ ii ];
	}

	if( xhi ) {

		if( is_wave ) {

			// wave analysis swaps edges, e.g. wave_lo <--> energy_hi
			// use xlo

			if ( 0.0 == xlo[ xlo.get_size() - 1 ] ) {
				PyErr_SetString( PyExc_ValueError,
						(char*)"XSPEC model evaluation failed, division by zero" );
				return NULL;
			}
			ear[ near - 1 ] = ( SherpaFloat ) (hc / xlo[ xlo.get_size() - 1 ]);
		}
		else
			ear[ near - 1 ] = ( SherpaFloat ) xhi[ xhi.get_size() - 1 ];

	}
	else
		nelem--;

	DoubleArray result;
	if ( EXIT_SUCCESS != result.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

        // Copy over the flux array:
        for (int i = 0; i < nelem; i++)
          result[i] = fluxes[i];

        // Since the flux error is discarded, it does not need to be a
        // NumPy array. Should be set to zeros for safety.
        std::vector<SherpaFloat> error(nelem);

	// Swallow C++ exceptions

	try {

                int ifl = 1;
		XSpecFunc( &ear[0], nelem, &pars[0], ifl, &result[0], &error[0], NULL );

                // See comments in xspecmodelfct template.

		while(!gaps.empty()) {
			std::vector<SherpaFloat> ear2(3);
			std::vector<SherpaFloat> result2(2);
			std::vector<SherpaFloat> error2(2);
                        double bin_width = gap_widths.back();
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + bin_width;
                        ear2[2] = ear2[1] + bin_width;
                        result2[0] = fluxes[bin_number]; // as this is a convolution
                        result2[1] = fluxes[bin_number+1];
			int ear2_nelem = 2;
			XSpecFunc( &ear2[0], ear2_nelem, &pars[0], ifl,
                                   &result2[0], &error2[0], NULL );

                        // Since error array is thrown away could just use
                        // error2 in the call above and not bother copying
                        // it over
                        result[bin_number] = result2[0];
                        error[bin_number] = error2[0];

			gaps.pop_back();
			gap_widths.pop_back();
		}

	} catch(...) {

		PyErr_SetString( PyExc_ValueError,
				(char*)"XSPEC convolution model evaluation failed" );
		return NULL;

	}

	// The XSPEC functions expect the output array to be of length nFlux
	// (one less than the input array), so set the last element to
	// zero to avoid having random garbage in it
	if( !xhi )
		result[ result.get_size() - 1 ] = 0.0;

	return result.return_new_ref();

}

template <bool HasNorm,
void (*XSpecFunc)( float* ear, int ne, float* param,
		const char* filenm, int ifl, float* photar,
		float* photer )>
PyObject* xspectablemodel( PyObject* self, PyObject* args, PyObject *kwds )
{

#ifdef INIT_XSPEC
	if ( EXIT_SUCCESS != INIT_XSPEC() )
		return NULL;
#endif

	FloatArray pars;
	DoubleArray xlo;
	DoubleArray xhi;
	DoubleArray *x;
	char *filename;
	static char *kwlist[] = {(char*)"pars", (char*)"xlo", (char*)"xhi",
			(char*)"filename", NULL};
	if ( !PyArg_ParseTupleAndKeywords( args, kwds, (char*)"O&O&|O&s", kwlist,
			(converter)convert_to_contig_array< FloatArray >,
			&pars,
			(converter)convert_to_contig_array< DoubleArray >,
			&xlo,
			(converter)convert_to_contig_array< DoubleArray >,
			&xhi,
			&filename) )
		return NULL;

	npy_intp npars = pars.get_size();

	int nelem = int( xlo.get_size() );

        if ( nelem < 2 ) {
          std::ostringstream err;
          err << "input array must have at least 2 elements, found " << nelem;
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

        if( xhi && (nelem != int(xhi.get_size())) ) {
          std::ostringstream err;
          err << "input arrays are not the same size: " << nelem
              << " and " << int( xhi.get_size() );
          PyErr_SetString( PyExc_TypeError, err.str().c_str() );
          return NULL;
        }

	// FIXME how to handle the spectrum number??
	int ifl = 1;

	double hc = (sherpa::constants::c_ang<SherpaFloat>() *
			sherpa::constants::h_kev<SherpaFloat>());
	bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

	std::vector<int> gaps;
	std::vector<double> gap_widths;

	// The XSPEC functions expect the input array to be of length ne+1
	int near = nelem;
	if( xhi ) {
		near++;

                // See comments in xspecmodelfct template.

		for (int i = 0; i < nelem-1; i++) {
			if (0 != sao_fcmp(xhi[i], xlo[i+1], DBL_EPSILON)) {
				gaps.push_back(i);
				double width = fabs(xhi[i] - xlo[i]);
				if( is_wave ) {
					width = hc / width;
				}
				gap_widths.push_back(width);
			}
		}
	}

        // TODO: should be documented why this can not use C++ vector
	//std::vector<FloatArrayType> ear(near);
	float *ear = NULL;
	ear = (float*)malloc(near*sizeof(float));

	for( int ii = 0; ii < nelem; ii++ ) {
		if( is_wave ) {

			// wave analysis swaps edges, e.g. wave_hi <--> energy_lo
			// if xhi is available use it
			x = (xhi) ? &xhi : &xlo;

			if ( 0.0 == (*x)[ii] ) {
				PyErr_SetString( PyExc_ValueError,
						(char*)"XSPEC model evaluation failed, division by zero" );
				return NULL;
			}
			ear[ ii ] = ( FloatArrayType ) (hc / (*x)[ ii ]);
		}
		else
			ear[ ii ] = ( FloatArrayType ) xlo[ ii ];
	}

	if( xhi ) {

		if( is_wave ) {

			// wave analysis swaps edges, e.g. wave_lo <--> energy_hi
			// use xlo

			if ( 0.0 == xlo[ xlo.get_size() - 1 ] ) {
				PyErr_SetString( PyExc_ValueError,
						(char*)"XSPEC model evaluation failed, division by zero" );
				return NULL;
			}
			ear[ near - 1 ] = ( FloatArrayType ) (hc / xlo[ xlo.get_size() - 1 ]);
		}
		else
			ear[ near - 1 ] = ( FloatArrayType ) xhi[ xhi.get_size() - 1 ];

	}
	else
		nelem--;

	FloatArray result;
	if ( EXIT_SUCCESS != result.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

        // Since the flux error is discarded, it does not need to be a
        // NumPy array. Should be set to zeros for safety.
        std::vector<FloatArrayType> error(nelem);

	// Even though the XSPEC model function is Fortran, it could call
	// C++ functions, so swallow exceptions here

	try {

		XSpecFunc( ear, nelem, &pars[0], filename, ifl, &result[0],
				&error[0] );

                // See comments in xspecmodelfct template.

		while(!gaps.empty()) {
                        // Arrays are small, so why use malloc?
                        float ear2[3] = { 0.0, 0.0, 0.0 };
                        float result2[2] = { 0.0, 0.0 };
                        float error2[2] = { 0.0, 0.0 };
                        double bin_width = gap_widths.back();
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + bin_width;
			ear2[2] = ear2[1] + bin_width;
			int ear2_nelem = 2;
                        XSpecFunc( ear2, ear2_nelem, &pars[0], filename, ifl,
                                   result2, error2);

                        // Since error array is thrown away could just use
                        // error2 in the call above and not bother copying
                        // it over
                        result[bin_number] = result2[0];
                        error[bin_number] = error2[0];

                        gaps.pop_back();
			gap_widths.pop_back();
		}

	} catch(...) {

		PyErr_SetString( PyExc_ValueError,
				(char*)"XSPEC table model evaluation failed" );
		return NULL;

	}

	// Apply normalization if required
	if ( HasNorm )
		for ( int ii = 0; ii < nelem; ii++ )
			result[ii] *= pars[npars - 1];

	// The XSPEC functions expect the output array to be of length ne
	// (one less than the input array), so set the last element to
	// zero to avoid having random garbage in it
	if( !xhi )
		result[ result.get_size() - 1 ] = 0.0;

	if( ear ) free(ear);

	return result.return_new_ref();

}

} } } /* namespace xspec, namespace astro, namespace sherpa */


#define _XSPECFCTSPEC(name, npars, has_norm) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct< npars, has_norm, \
				name##_ >))

#define XSPECMODELFCT(name, npars)  _XSPECFCTSPEC(name, npars, false)
#define XSPECMODELFCT_NORM(name, npars)  _XSPECFCTSPEC(name, npars, true)

#define XSPECMODELFCT_C(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, false, name >))

#define XSPECMODELFCT_C_NORM(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, true, name >))

#define XSPECMODELFCT_CON(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_con< npars, name >))

#define _XSPECTABLEMODELSPEC(name, has_norm) \
		{ (char*)#name, \
	(PyCFunction)((PyCFunctionWithKeywords)sherpa::astro::xspec::xspectablemodel< has_norm, name >), \
	METH_VARARGS|METH_KEYWORDS, \
	NULL }

#define XSPECTABLEMODEL(name) \
		_XSPECTABLEMODELSPEC(name, false)

#define XSPECTABLEMODEL_NORM(name) \
		_XSPECTABLEMODELSPEC(name, true)


#endif /* __sherpa_astro_xspec_extension_hh__ */
