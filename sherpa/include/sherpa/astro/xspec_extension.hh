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
// gaps, and then dealing with them. The approach is to
// just create a grid as if there were no gaps; that is, the model would
// be evaluated on the grid [0.1,0.2,0.6,0.7,0.8] but record the position
// and correct widths for the "gaps" and then re-run the model to
// "fill in" the output array, on a bin-by-bin basis. In this case
// the model would be called a second time to evaluate the model
// for the grid [0.2,0.3] and this value inserted into the second
// bin of the flux array created by the first call. This can lead
// to errors (in XSPEC 12.8.2 the apec-style models will crash when
// called with only two bins; this has been fixed in 12.9.0).

// When creating the flux and flux error arrays to be sent to
// the XSpec routines, the arrays are filled with 0's - that is
// the zeros array method is used, rather than create. This is
// to ensure that a "sensible" answer is returned on error (all
// 0's), as the XSpec model API does not provide a way to return
// an error status. There is (likely; not checked) a run-time cost
// to using zeros rather than create, but safety is better than
// performance here.

// The spectrum number is set to 1. It used to be 0, but some code
// (a user model, so not included in the XSpec model library being
// built against) has been seen to behave strangely with a value of 0,
// so a value of 1 is being used "just in case". This is also the
// approach that XSPEC (the application) uses. A keyword argument
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
	DoubleArray xlo;
	DoubleArray xhi;
	DoubleArray *x;

	if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
			(converter)convert_to_contig_array< FloatArray >,
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

	double hc = (sherpa::constants::c_ang<SherpaFloat>() *
			sherpa::constants::h_kev<SherpaFloat>());
	bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

	std::vector<int> gaps;
	std::vector<double> gap_widths;

	// The XSPEC functions expect the input array to be of length ne+1
	int near = nelem;
	if( xhi ) {
		near++;

		// Suppose the data were filtered, such that there is a gap
		// in the middle of the energy array.  In that case *only*,
		// we will find that xlo[i+1] != xhi[i].  However, XSPEC models
		// expect that xlo[i+1] == xhi[i].
		//
		// So, if we pass in filtered data and xlo[i+1] != xhi[i],
		// then at energy bin i we will end up calculating an energy
		// flux that is far too great.  We will correct that by gathering
		// information to allow us to recalculate individual bins, with
		// boundaries xlo[i], xhi[i], to correct for cases where
		// boundaries xlo[i], xlo[i+1] results in a bin that is too big.
		//
		// We will gather the locations of the gaps here, and calculate
		// actual widths based on xhi[i] - xlo[i] downstream.
		//
		// If we are working in wavelength space we will also correct for that.
		// SMD 11/21/12.

		for (int i = 0; i < nelem-1; i++) {
			double cmp;
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

	std::vector<FloatArrayType> ear(near);

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

	// The XSPEC functions require fluxError to be non-NULL, so we create
	// it but discard it after the computation is done
	FloatArray error;
	if ( EXIT_SUCCESS != error.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

	// Even though the XSPEC model function is Fortran, it could call
	// C++ functions, so swallow exceptions here

	try {

		XSpecFunc( &ear[0], &nelem, &pars[0], &ifl, &result[0], &error[0] );

		// If there were gaps in the energy array, because of locations
		// where xlo[i+1] != xhi[i], then this is place where we recalculate
		// energy fluxes for those bins *only*.
		//
		// For each such location in the energy grid, construct a new
		// 2-bin energy array, such that the 2-bin array is [xlo[i],
		// xhi[i]].  This is accomplished by:
		//
		// ear2[0] = ear[location of gap]
		// ear2[1] = ear[location of gap] + abs(xhi[location of gap] -
		//                                      xlo[location of gap])
		// The locations of the gaps, and the actual widths of the energy
		// bins at those locations, were calculated above.  So use the
		// gaps and gap_widths vectors here to recalculate energy fluxes
		// at affected bins *only*. SMD 11/21/12

		while(!gaps.empty()) {
			std::vector<FloatArrayType> ear2(2);
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + gap_widths.back();
			int ear2_nelem = 1;
			XSpecFunc( &ear2[0], &ear2_nelem, &pars[0], &ifl, &result[bin_number],
					&error[bin_number]);

			gaps.pop_back();
			gap_widths.pop_back();
		}

	} catch(...) {

		PyErr_SetString( PyExc_ValueError,
				(char*)"XSPEC model evaluation failed" );
		return NULL;

	}

	// Apply normalization if required
	if ( HasNorm )
		for ( int ii = 0; ii < nelem; ii++ )
			result[ii] *= pars[NumPars - 1];

	// The XSPEC functions expect the output array to be of length ne
	// (one less than the input array), so set the last element to
	// zero to avoid having random garbage in it
	if( !xhi )
		result[ result.get_size() - 1 ] = 0.0;

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
	DoubleArray *x;

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

	double hc = (sherpa::constants::c_ang<SherpaFloat>() *
			sherpa::constants::h_kev<SherpaFloat>());
	bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;

	std::vector<int> gaps;
	std::vector<double> gap_widths;

	// The XSPEC functions expect the input array to be of length nFlux+1
	int near = nelem;
	if( xhi ) {
		near++;

		// Suppose the data were filtered, such that there is a gap
		// in the middle of the energy array.  In that case *only*,
		// we will find that xlo[i+1] != xhi[i].  However, XSPEC models
		// expect that xlo[i+1] == xhi[i].
		//
		// So, if we pass in filtered data and xlo[i+1] != xhi[i],
		// then at energy bin i we will end up calculating an energy
		// flux that is far too great.  We will correct that by gathering
		// information to allow us to recalculate individual bins, with
		// boundaries xlo[i], xhi[i], to correct for cases where
		// boundaries xlo[i], xlo[i+1] results in a bin that is too big.
		//
		// We will gather the locations of the gaps here, and calculate
		// actual widths based on xhi[i] - xlo[i] downstream.
		//
		// If we are working in wavelength space we will also correct for that.
		// SMD 11/21/12.

		for (int i = 0; i < nelem-1; i++) {
			double cmp;
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

	// The XSPEC functions require fluxError to be non-NULL, so we create
	// it but discard it after the computation is done
	DoubleArray error;
	if ( EXIT_SUCCESS != error.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

	// Swallow C++ exceptions

	try {

                int ifl = 1;
		XSpecFunc( &ear[0], nelem, &pars[0], ifl, &result[0], &error[0], NULL );

		// If there were gaps in the energy array, because of locations
		// where xlo[i+1] != xhi[i], then this is place where we recalculate
		// energy fluxes for those bins *only*.
		//
		// For each such location in the energy grid, construct a new
		// 2-bin energy array, such that the 2-bin array is [xlo[i],
		// xhi[i]].  This is accomplished by:
		//
		// ear2[0] = ear[location of gap]
		// ear2[1] = ear[location of gap] + abs(xhi[location of gap] -
		//                                      xlo[location of gap])
		// The locations of the gaps, and the actual widths of the energy
		// bins at those locations, were calculated above.  So use the
		// gaps and gap_widths vectors here to recalculate energy fluxes
		// at affected bins *only*. SMD 11/21/12

		while(!gaps.empty()) {
			std::vector<SherpaFloat> ear2(2);
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + gap_widths.back();
			int ear2_nelem = 1;
			XSpecFunc( &ear2[0], ear2_nelem, &pars[0], ifl, &result[bin_number],
					&error[bin_number], NULL );

			gaps.pop_back();
			gap_widths.pop_back();
		}

	} catch(...) {

		PyErr_SetString( PyExc_ValueError,
				(char*)"XSPEC model evaluation failed" );
		return NULL;

	}

	// Apply normalization if required
	if ( HasNorm )
		for ( int ii = 0; ii < nelem; ii++ )
			result[ii] *= pars[NumPars - 1];

	// The XSPEC functions expect the output array to be of length nFlux
	// (one less than the input array), so set the last element to
	// zero to avoid having random garbage in it
	if( !xhi )
		result[ result.get_size() - 1 ] = 0.0;

	return result.return_new_ref();

}

// Handle convolution models, which are assumed to have C-style
// linkage.
//
// This template does not support non-contiguous grids.
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
        //   pars, fluxes, xlo
        //   pars, fluxes, xlo, xhi
        // where fluxes is the spectrum that is to be convolved
        // by the model.
        //
	if ( !PyArg_ParseTuple( args, (char*)"O&O&O&|O&",
			(converter)convert_to_contig_array< DoubleArray >,
			&pars,
			(converter)convert_to_contig_array< DoubleArray >,
			&fluxes,
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

	// The XSPEC functions expect the input array to be of length nFlux+1
	int near = nelem;
	if( xhi ) {
		near++;

                // raise an error if the grid is not contiguous
		for (int i = 0; i < nelem-1; i++) {
			double cmp;
			if ( is_wave ) {
				cmp = sao_fcmp(xlo[i], xhi[i+1], DBL_EPSILON);
			} else {
				cmp = sao_fcmp(xhi[i], xlo[i+1], DBL_EPSILON);
			}
			if (0 != cmp) {
                          PyErr_SetString( PyExc_ValueError,
                                           (char*)"XSPEC convolution model requires a contiguous grid" );
                          return NULL;
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

        for ( int ii = 0; ii < nelem; ii++ )
          result[ii] = fluxes[ii];

	// The XSPEC functions require fluxError to be non-NULL, so we create
	// it but discard it after the computation is done
	DoubleArray error;
	if ( EXIT_SUCCESS != error.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

	// Swallow C++ exceptions

	try {

                int ifl = 1;
		XSpecFunc( &ear[0], nelem, &pars[0], ifl, &result[0], &error[0], NULL );

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

		// Suppose the data were filtered, such that there is a gap
		// in the middle of the energy array.  In that case *only*,
		// we will find that xlo[i+1] != xhi[i].  However, XSPEC models
		// expect that xlo[i+1] == xhi[i].
		//
		// So, if we pass in filtered data and xlo[i+1] != xhi[i],
		// then at energy bin i we will end up calculating an energy
		// flux that is far too great.  We will correct that by gathering
		// information to allow us to recalculate individual bins, with
		// boundaries xlo[i], xhi[i], to correct for cases where
		// boundaries xlo[i], xlo[i+1] results in a bin that is too big.
		//
		// We will gather the locations of the gaps here, and calculate
		// actual widths based on xhi[i] - xlo[i] downstream.
		//
		// If we are working in wavelength space we will also correct for that.
		// SMD 11/21/12.

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

	// The XSPEC functions require fluxError to be non-NULL, so we create
	// it but discard it after the computation is done
	FloatArray error;
	if ( EXIT_SUCCESS != error.zeros( xlo.get_ndim(), xlo.get_dims() ) )
		return NULL;

	// Even though the XSPEC model function is Fortran, it could call
	// C++ functions, so swallow exceptions here

	try {

		XSpecFunc( ear, nelem, &pars[0], filename, ifl, &result[0],
				&error[0] );

		// If there were gaps in the energy array, because of locations
		// where xlo[i+1] != xhi[i], then this is place where we recalculate
		// energy fluxes for those bins *only*.
		//
		// For each such location in the energy grid, construct a new
		// 2-bin energy array, such that the 2-bin array is [xlo[i],
		// xhi[i]].  This is accomplished by:
		//
		// ear2[0] = ear[location of gap]
		// ear2[1] = ear[location of gap] + abs(xhi[location of gap] -
		//                                      xlo[location of gap])
		// The locations of the gaps, and the actual widths of the energy
		// bins at those locations, were calculated above.  So use the
		// gaps and gap_widths vectors here to recalculate energy fluxes
		// at affected bins *only*. SMD 11/21/12

		while(!gaps.empty()) {
			float *ear2 = NULL;
			ear = (float*)malloc(2*sizeof(float));
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + gap_widths.back();
			int ear2_nelem = 1;
			try {
				XSpecFunc( ear2, ear2_nelem, &pars[0], filename, ifl,
						&result[bin_number], &error[bin_number]);
			} catch(...) {
				if (ear2) free(ear2);
				throw;
			}
			gaps.pop_back();
			gap_widths.pop_back();
			if (ear2) free(ear2);
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
