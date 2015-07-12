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

        // Since the flux error is discarded, it does not need to be a
        // NumPy array. Should be set to zeros for safety.
        std::vector<FloatArrayType> error(nelem);

	// Even though the XSPEC model function is Fortran, it could call
	// C++ functions, so swallow exceptions here

	try {

		XSpecFunc( &ear[0], &nelem, &pars[0], &ifl, &result[0], &error[0] );

		// If there were gaps in the energy array, because of locations
		// where xlo[i+1] != xhi[i], then this is place where we recalculate
		// energy fluxes for those bins *only*.
		//
		// For each such location in the energy grid, construct a new
		// 3-bin energy array:
		//
		// ear2[0] = ear[location of gap]
		// ear2[1] = ear[0] + abs(xhi[location of gap] -
		//                        xlo[location of gap])
		// ear2[2] = ear[1] + abs(xhi[location of gap] -
		//                        xlo[location of gap])
		// The locations of the gaps, and the actual widths of the energy
		// bins at those locations, were calculated above.  So use the
		// gaps and gap_widths vectors here to recalculate energy fluxes
		// at affected bins *only*. SMD 11/21/12
                //
                // Changed from using 2 to 3 bins - because XSpec models
                // that call the SUMDEM routins, such as apec, crash if
                // called with only 2 bins in 12.8.2q. The choice for
                // the value of third bin is technically arbitrary, but it
                // does actually make a difference for some XSpec models
                // (presumably because they interpolate answers from an
                // internally-based grid onto the user-requested grid).
                // July 9 2015
		while(!gaps.empty()) {
			std::vector<FloatArrayType> ear2(3);
                        std::vector<FloatArrayType> result2(2);
                        std::vector<FloatArrayType> error2(2);
                        double bin_width = gap_widths.back();
			int bin_number = gaps.back();
			ear2[0] = ear[bin_number];
			ear2[1] = ear2[0] + bin_width;
			ear2[2] = ear2[1] + bin_width;
                        result2[0] = 0.0; // in case of error
                        result2[1] = 0.0;
			int ear2_nelem = 2;
			XSpecFunc( &ear2[0], &ear2_nelem, &pars[0], &ifl,
                                   &result2[0], &error2[0] );

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
	DoubleArray fluxes;
	DoubleArray *x;

	if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&O&",
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
        if( fluxes && (nelem != int(fluxes.get_size())) ) {
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

	if (fluxes) {
		for ( int ii = 0; ii < nelem; ii++ )
			result[ii] = fluxes[ii];
	}

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
                        result2[0] = 0.0; // in case of error
                        result2[1] = 0.0;
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
