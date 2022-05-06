//
//  Copyright (C) 2009, 2015, 2017, 2020, 2021, 2022
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
#include <algorithm>
#include <stdexcept>

#include "sherpa/fcmp.hh"

// We should be able to just include funcType.h but our XSPEC conda
// builds, at least for testing/12.11.1, do not include this file,
// so we just include what we need.
//
// #ifdef XSPEC_12_12_0
// #include "XSFunctions/Utilities/funcType.h"
// #else
// #include "funcType.h"
// #endif

#include "xsTypes.h"   // get Real typedef

extern "C" {

        typedef void (xsf77Call) (const float* energyArray,
                                  const int& Nenergy,
                                  const float* parameterValues,
                                  const int& spectrumNumber,
                                  float* flux,
                                  float* fluxError);

        typedef void (xsF77Call) (const double* energyArray,
                                  const int& Nenergy,
                                  const double* parameterValues,
                                  const int& spectrumNumber,
                                  double* flux,
                                  double* fluxError);

        typedef void (xsccCall)   (const Real* energyArray,
                                   int Nenergy,
                                   const Real* parameterValues,
                                   int spectrumNumber,
                                   Real* flux,
                                   Real* fluxError,
                                   const char* initString);
}

// Prior to XSPEC 12.10.1, the table models were split into different
// functions. These functions are defined in _xspec.cc.
//
// In 12.10.1 they were consolidated into a single function, tabint,
// and so the declaration was moved here. The function was only
// available in C++ scope.
//
// In XSPEC 12.11.0 (the next one after 12.10.1), the tabint function
// was moved into C scope. In XSPEC 12.12.1 the signature was changed
// to mark more arguments as const.
//
#ifdef XSPEC_12_10_1
#ifdef XSPEC_12_11_0
extern "C" {
#endif

#ifndef XSPEC_12_12_1
  void tabint(float* ear, int ne, float* param,
	      int npar, const char* filenm, int ifl,
	      const char* tabtyp, float* photar, float* photer);
#else
  void tabint(const float* ear, const int ne, const float* param,
	      const int npar, const char* filenm, int ifl,
	      const char* tabtyp, float* photar, float* photer);
#endif

#ifdef XSPEC_12_11_0
}
#endif

#endif


namespace sherpa { namespace astro { namespace xspec {

      // Track the type of Python error we want to throw:
      //  - NoError when the python error is already set
      //  - TypeError when a python TypeError is to be raised
      //  - ValueError when a python ValueError is to be raised
      //
      struct NoError : std::runtime_error
      {
	NoError(const char* msg) : std::runtime_error(msg) {}
	NoError(const std::string& msg) : std::runtime_error(msg) {}
      };

      struct ValueError : std::runtime_error
      {
	ValueError(const char* msg) : std::runtime_error(msg) {}
	ValueError(const std::string& msg) : std::runtime_error(msg) {}
      };

      struct TypeError : std::runtime_error
      {
	TypeError(const char* msg) : std::runtime_error(msg) {}
	TypeError(const std::string& msg) : std::runtime_error(msg) {}
      };

typedef sherpa::Array< float, NPY_FLOAT > FloatArray;
typedef float FloatArrayType;

// Try and support the use of std::transform while still building
// against C++-98 compilers.
//
#if __cplusplus > 199711L
#define CONVERTARRAY(orig, out, npts)					\
        std::transform(std::begin(orig), std::end(orig), std::begin(out), \
                       [](const double val) -> FloatArrayType { return static_cast<FloatArrayType>(val); });
#else
#define CONVERTARRAY(orig, out, npts)					\
	for (int i = 0; i < npts; i++) { \
          out[i] = static_cast<FloatArrayType>(orig[i]); \
        }
#endif


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
// gaps, and then dealing with them. Any gaps are identified, and
// bins added to the arrays sent to the XSPEC model. After calling
// the models, the excess bins are removed. This has been discussed
// with Keith Arnaud as a sensible approach. An alternative would be
// to call the model on each contiguous section, but the issue here
// is that there may be a non-negligible set-up cost within the models.

// When creating the flux and flux error arrays to be sent to
// the XSpec routines, the arrays are filled with 0's - that is
// the zeros array method is used, rather than create. This is
// to ensure that a "sensible" answer is returned on error (all
// 0's), as the XSpec model API does not provide a way to return
// an error status. There is (likely; not checked) a run-time cost
// to using zeros rather than create, but safety is better than
// performance here. There are also problems with some models in
// XSPEC 12.8.2 (not in 12.9.0) where they would crash if called with
// only a single bin.

// The convention is that if the input grids are in ascending order they
// are in keV - the units required by XSpec - otherwise they are
// in Angstrom, so must be converted to keV. Note that the constraint
// xhi > xlo is expected to hold even when the units are Angstrom -
// e.g. xlo = 112.7, 103.3, 95.4, ...
//      xhi = 124.0, 112.7, 103.3, ...
// where the values are in Angstroms.

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


// Perform sanity checks and then create the grid of points that
// is sent to XSPEC. Complications include
//    - convert from Angstrom to KeV (this is based purely on the
//      grid being in descending order, and not due to any sort of
//      unit checking)
//    - if sent low and high bin edges, this needs to be converted
//      into a single array AND any gaps need a bin inserted
//      (these extra bins are removed by finalize_grid)
//
// This is only ever used with a SherpaFloat (aka double), and its
// internals rely on a double, due to sao_fcmp, but it's
// left as a template for now.
//
template <typename CType, int ArrayType>
static void create_grid(const sherpa::Array<CType, ArrayType> &xlo,
			const sherpa::Array<CType, ArrayType> &xhi,
			std::vector<CType> &ear,
			std::vector<int> &gaps_index ) {
  std::vector<CType> gaps_edges;

  int nelem = int( xlo.get_size() );
  if ( nelem < 2 ) {
    std::ostringstream err;
    err << "input array must have at least 2 elements, found " << nelem;
    throw TypeError(err.str());
  }

  if( xhi && (nelem != int(xhi.get_size())) ) {
    std::ostringstream err;
    err << "input arrays are not the same size: " << nelem
        << " and " << int( xhi.get_size() );
    throw TypeError(err.str());
  }

  bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;
  const sherpa::Array<CType, ArrayType> *x1 = (is_wave && xhi) ? &xhi : &xlo;
  const sherpa::Array<CType, ArrayType> *x2 = (is_wave && xhi) ? &xlo : &xhi;

  // assume the gaps array is empty on input

  if (xhi) {
    const int gap_found = is_wave ? 1 : -1;
    for (int i = 0; i < nelem-1; i++) {
      int cmp = sao_fcmp((*x2)[i], (*x1)[i+1], DBL_EPSILON);
      if (cmp == gap_found) {
        gaps_index.push_back(i);
        gaps_edges.push_back((*x2)[i]);
        /*** DO NOT INCLUDE THIS CHECK YET, AS UNSURE IF
             IT IS GOING TO CAUSE PROBLEMS, GIVEN THAT
             ARF/RMF GRIDS CAN BE POORLY DEFINED
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
        return false;
        ***/
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
  // Do a two-step conversion; create the array and then a type
  // conversion (which is excessive if no grid points are added
  // in, and the input is in keV).
  //
  ear.assign(ngrid, 0);

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

  // Convert from wavelength (Angstrom) to energy (keV)
  //
  if (is_wave) {
    CType hc = (sherpa::constants::c_ang<CType>() *
                 sherpa::constants::h_kev<CType>());
    for (int i = 0; i < ngrid; i++) {
      if (ear[i] <= 0.0) {
        std::ostringstream err;
        err << "Wavelength must be > 0, sent " << ear[i];
	throw ValueError(err.str());
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
  /*** DO NOT INCLUDE FOR THE SAME REASON AS ABOVE, AS
       UNSURE ABOUT ARF/RMF GRIDS
  for (int i = 0; i < ngrid - 1; i++) {
    if (ear[i] >= ear[i+1]) {
      std::ostringstream err;
      err << "Grid is not monotonic: " << ear[i] << " to " <<
        ear[i+1];
      PyErr_SetString( PyExc_ValueError, err.str().c_str() );
      return false;
    }
  }
  ***/

} /* create_grid */


// Similar to create_grid except that it does not support grids that are
// not contiguous (an error is raised).
//
template <typename CType, int ArrayType>
static void create_contiguous_grid(const sherpa::Array<CType, ArrayType> &xlo,
				   const sherpa::Array<CType, ArrayType> &xhi,
				   std::vector<CType> &ear) {

  int nelem = int( xlo.get_size() );
  if ( nelem < 2 ) {
    std::ostringstream err;
    err << "input array must have at least 2 elements, found " << nelem;
    throw TypeError(err.str());
  }

  if( xhi && (nelem != int(xhi.get_size())) ) {
    std::ostringstream err;
    err << "input arrays are not the same size: " << nelem
        << " and " << int( xhi.get_size() );
    throw TypeError(err.str());
  }

  bool is_wave = (xlo[0] > xlo[nelem-1]) ? true : false;
  const sherpa::Array<CType, ArrayType> *x1 = (is_wave && xhi) ? &xhi : &xlo;
  const sherpa::Array<CType, ArrayType> *x2 = (is_wave && xhi) ? &xlo : &xhi;

  // Are there any non-contiguous bins? The check lets through
  // overlapping bins.
  if (xhi) {
    const int gap_found = is_wave ? 1 : -1;
    for (int i = 0; i < nelem-1; i++) {
      int cmp = sao_fcmp((*x2)[i], (*x1)[i+1], DBL_EPSILON);
      if (cmp == gap_found) {
	/*** Maybe confusing to users
	     std::ostringstream err;
	     err << "Grid cells are not contiguous: cell " << i
	     << " (" << (*x1)[i] << " to " << (*x2)[i] << ")"
	     << " and cell " << (i+1)
	     << " (" << (*x1)[i+1] << " to " << (*x2)[i+1] << ")";
	     PyErr_SetString( PyExc_ValueError, err.str().c_str() );
	***/
	throw ValueError("XSPEC convolution model requires a contiguous grid");
      }
    }
  }

  // The size of the energy array sent to XSpec
  int ngrid = nelem;
  if (xhi) {
    ngrid += 1;
  }

  // XSpec traditionally refers to the input energy grid as ear.
  // Do a two-step conversion; create the array and then a type
  // conversion (which is excessive if no grid points are added
  // in, and the input is in keV).
  //
  ear.assign(ngrid, 0);

  // The grid is created, converted from Angstrom to Energy
  // (if required), and then checked for being monotonic.
  // The code has been kept similar to create_grid.
  //
  for(int i = 0; i < nelem; i++) {
    ear[i] = (*x1)[i];
  }

  // Add on the last bin value if needed
  if (xhi) {
    ear[ngrid - 1] = (*x2)[nelem - 1];
  }

  // Convert from Angstrom to keV
  if (is_wave) {
    CType hc = (sherpa::constants::c_ang<CType>() *
                 sherpa::constants::h_kev<CType>());
    for (int i = 0; i < ngrid; i++) {
      if (ear[i] <= 0.0) {
        std::ostringstream err;
        err << "Wavelength must be > 0, sent " << ear[i];
	throw ValueError(err.str());
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
  /*** DO NOT INCLUDE FOR THE SAME REASON AS ABOVE, AS
       UNSURE ABOUT ARF/RMF GRIDS
  for (int i = 0; i < ngrid - 1; i++) {
    if (ear[i] >= ear[i+1]) {
      std::ostringstream err;
      err << "Grid is not monotonic: " << ear[i] << " to " <<
        ear[i+1];
      PyErr_SetString( PyExc_ValueError, err.str().c_str() );
      return false;
    }
  }
  ***/

} /* create_contiguous_grid */


// Remove any elements that were inserted to deal with gaps in the grid.
//
template <typename CType, int ArrayType>
static void finalize_grid(int nelem,
			  sherpa::Array< CType, ArrayType > &result,
			  std::vector<int> &gaps_index) {

  // Remove gaps. It is assumed there are some.
  int ngaps = (int) gaps_index.size();
  // if (ngaps <= 0) { return; }  this check assumed to be made by the caller

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

} /* finalize_grid */


template <typename T>
static void create_output(int nbins, T &a, T &b) {

  npy_intp dims[1] = { nbins };

  // It's never clear from the Python documentation whether we should throw
  // some form of memory-allocation error here, but let's throw a ValueError
  // (if it fails then we are likely in serious trouble anyway so it's not
  // worth worrying about the details).
  //
  if ( (EXIT_SUCCESS != a.zeros( 1, dims )) ||
       (EXIT_SUCCESS != b.zeros( 1, dims )) )
    throw ValueError("Unable to create output");

} /* create_output */

      class PyArgTupleBase {
      public:
        PyArgTupleBase( ) {
#ifdef INIT_XSPEC
	if ( EXIT_SUCCESS != INIT_XSPEC() )
          throw std::runtime_error("Unable to initialize XSpec");
#endif
        }

      };

      template <typename RealArray>
      class PyArgTuple3 : public PyArgTupleBase {
      public:
        PyArgTuple3( PyObject* args, npy_intp NumPars, RealArray& pars,
                     DoubleArray& xlo, DoubleArray& xhi ) :
          PyArgTupleBase( ) {
          //
          // The grid arrays could be cast to FloatArray here, saving
          // conversion later on in this routine. However, that can then
          // lead to differences in the identification of non-contiguous
          // bins, or whether a grid is monotonic and non-overlapping [*]
          // (e.g. if a source expression contains both a FORTRAN
          // and C style model), so stick to this approach for now.
          //
          // [*] although these checks are currently commented out
          //
          if ( !PyArg_ParseTuple( args, (char*)"O&O&|O&",
                                  (converter)convert_to_contig_array< RealArray >,
                                  &pars,
                                  (converter)convert_to_contig_array< DoubleArray >,
                                  &xlo,
                                  (converter)convert_to_contig_array< DoubleArray >,
                                  &xhi ) )
            throw NoError("Error Parsing args");

          npy_intp npars = pars.get_size();
          if ( NumPars != npars ) {
            std::ostringstream err;
            err << "expected " << NumPars << " parameters, got " << npars;
            throw TypeError(err.str());
          }
          return;
        };
      };

      template <typename RealArray>
      class PyArgTuple4 : public PyArgTupleBase {
      public:
        PyArgTuple4( PyObject* args, npy_intp NumPars, RealArray& pars,
                     RealArray& fluxes, DoubleArray& xlo, DoubleArray& xhi ) :
          PyArgTupleBase( ) {
          //
          // The arguments are parsed as
          //   pars, fluxes, xlo
          //   pars, fluxes, xlo, xhi
          // where fluxes is the spectrum that is to be convolved
          // by the model.
          //
          if ( !PyArg_ParseTuple( args, (char*)"O&O&O&|O&",
                                  (converter)convert_to_contig_array< RealArray >,
                                  &pars,
                                  (converter)convert_to_contig_array< RealArray >,
                                  &fluxes,
                                  (converter)convert_to_contig_array< DoubleArray >,
                                  &xlo,
                                  (converter)convert_to_contig_array< DoubleArray >,
                                  &xhi ) )
            throw NoError("Error Parsing args");

          npy_intp npars = pars.get_size();
          if ( NumPars != npars ) {
            std::ostringstream err;
            err << "expected " << NumPars << " parameters, got " << npars;
            throw TypeError(err.str());
          }
          return;
        };
      };

      template<typename Real, typename RealArray>
      class xspecModelFctBase {
      public:

        virtual ~xspecModelFctBase( ) { }

        xspecModelFctBase( PyObject* arg, npy_intp numpars, bool hasnorm )
          : args(arg), ifl(1), NumPars(numpars), HasNorm(hasnorm) { }

        virtual void call_xspec( RealArray& result ) { }

        void eval( RealArray& result ) {

          PyArgTuple3<RealArray>( args, NumPars, pars, xlo, xhi );

          //
          // The grid to send to XSPEC (double precision).
          //
          create_grid(xlo, xhi, ear, gaps_index);

          nelem = int( xlo.get_size() );
          ngrid = ear.size();
          npts = ngrid - 1;

          // Number of bins to send to XSPEC
          int nout = ngrid;
          if (xhi) nout--;

          create_output(nout, result, error);

          call_xspec( result );

          int ngaps = gaps_index.size();
          if (ngaps > 0) {
            finalize_grid(nelem, result, gaps_index);
          }

          // Apply normalization if required
          if ( HasNorm )
            for (int i = 0; i < nelem; i++)
              result[i] *= pars[NumPars - 1];

        } // eval

        const char* get_err_msg( ) { return "XSPEC model evaluation failed"; }

      protected:

        PyObject* args;
        int ifl, ngrid, npts;
        std::vector<SherpaFloat> ear;
        RealArray pars, error;

      private:

        int NumPars, HasNorm, nelem;
        std::vector<int> gaps_index;
        DoubleArray xlo, xhi;

      }; // class xspecModelFctBase

      template<typename Real, typename RealArray, xsccCall XSpecFunc>
      class xspecModelFctC : public xspecModelFctBase<Real, RealArray>  {
      public:

        xspecModelFctC( PyObject* args, npy_intp NumPars, bool HasNorm )
          : xspecModelFctBase<Real, RealArray>( args, NumPars, HasNorm ) { }

        void call_xspec( RealArray& result ) {
          XSpecFunc( &this->ear[0], this->npts, &this->pars[0], this->ifl,
                     &result[0], &this->error[0], NULL );
          return;
        }

      }; // class xspecModelFctC

      template<typename Real, typename RealArray, xsf77Call XSpecFunc>
      class xspecModelFctF : public xspecModelFctBase<Real, RealArray>  {
      public:

        xspecModelFctF( PyObject* args, npy_intp NumPars, bool HasNorm )
          : xspecModelFctBase<Real, RealArray>( args, NumPars, HasNorm ) { }

        void call_xspec( RealArray& result ) {
          // convert to 32-byte float
          std::vector<float> fear(this->ngrid);
          CONVERTARRAY(this->ear, fear, this->ngrid);
          XSpecFunc( &fear[0], this->npts, &this->pars[0], this->ifl,
                     &result[0], &this->error[0] );
          return;
        }

      }; // class xspecModelFctF


      template<typename Real, typename RealArray, xsF77Call XSpecFunc>
      class xspecModelFctFD : public xspecModelFctBase<Real, RealArray>  {
      public:

        xspecModelFctFD( PyObject* args, npy_intp NumPars, bool HasNorm )
          : xspecModelFctBase<Real, RealArray>( args, NumPars, HasNorm ) { }

        void call_xspec( RealArray& result ) {
          XSpecFunc( &this->ear[0], this->npts, &this->pars[0], this->ifl,
                     &result[0], &this->error[0] );
          return;
        }

      }; // class xspecModelFctFD


      template<typename Real, typename RealArray>
      class xspecModelFctConBase {
      public:

        virtual ~xspecModelFctConBase( ) { }

        xspecModelFctConBase( PyObject* arg, npy_intp numpars )
          : args(arg), ifl(1), NumPars(numpars) { }

        virtual void call_xspec( RealArray& result ) { }

        void eval( RealArray& result ) {

          PyArgTuple4<RealArray>( args, NumPars, pars, fluxes, xlo, xhi );

          // XSpec traditionally refers to the input energy grid as ear.
          // std::vector<SherpaFloat> ear;
          create_contiguous_grid(xlo, xhi, ear);

          nelem = xlo.get_size();
          ngrid = ear.size();
          npts = ngrid - 1;

          // For now require the fluxes array to have the same
          // size as the input grid. If xhi is not given then
          // technically fluxes should be one less, but this is
          // likely to cause problems (as it doesn't match how
          // the rest of the interface works).
          if( nelem != fluxes.get_size() ) {
            std::ostringstream err;
            err << "flux array does not match the input grid: " << nelem
                << " and " << fluxes.get_size();
            throw ValueError(err.str());
          }

          // Number of bins to send to XSPEC
          int nout = ngrid;
          if (xhi) nout--;

          create_output(nout, result, error);

          // Copy over the flux array
          std::copy(&fluxes[0], &fluxes[nout], &result[0]);

          call_xspec( result );

        } // eval

        const char* get_err_msg( ) { return "XSPEC convolution model evaluation failed"; }

      protected:

        PyObject* args;
        int ifl, ngrid, npts;;
        std::vector<SherpaFloat> ear;
        RealArray pars, error;

      private:

        int NumPars, HasNorm, nelem;
        std::vector<int> gaps_index;
        DoubleArray xlo, xhi;
        RealArray fluxes;

      }; // class xspecModelFctConBase

      template<typename Real, typename RealArray, xsccCall XSpecFunc>
      class xspecModelFctConC : public xspecModelFctConBase<Real, RealArray>  {
      public:

        xspecModelFctConC( PyObject* arg, npy_intp numpars )
          : xspecModelFctConBase<Real, RealArray>( arg, numpars ) { }

        void call_xspec( RealArray& result ) {
          XSpecFunc( &this->ear[0], this->npts, &this->pars[0], this->ifl,
                     &result[0], &this->error[0], NULL );
          return;
        }

      }; // class xspecModelFctConC

      template<typename Real, typename RealArray, xsf77Call XSpecFunc>
      class xspecModelFctConF : public xspecModelFctConBase<Real, RealArray>  {
      public:

        xspecModelFctConF( PyObject* args, npy_intp NumPars )
          : xspecModelFctConBase<Real, RealArray>( args, NumPars ) { }

        void call_xspec( RealArray& result ) {
          // convert to 32-byte float
          std::vector<float> fear(this->ngrid);
          CONVERTARRAY(this->ear, fear, this->ngrid);
          XSpecFunc( &fear[0], this->npts, &this->pars[0], this->ifl,
                     &result[0], &this->error[0] );
          return;
        }

      }; // class xspecModelFctConF

      class xspecTableModelBase {
      public:

        virtual ~xspecTableModelBase( ) { }

        xspecTableModelBase( int arg ) : ifl(arg) { }

        virtual void call_xspec( std::vector<FloatArrayType>& fear,
                                 FloatArray& pars, npy_intp npars,
                                 char* filename, char*tabtype,
                                 FloatArray& result ) { }

        void eval( bool HasNorm, DoubleArray& xlo, DoubleArray& xhi,
                   npy_intp npars, FloatArray& pars, char* filename,
                   char* tabtype, FloatArray& result ) {

          //
          // The grid to send to XSPEC (double precision).
          //
          create_grid(xlo, xhi, ear, gaps_index);

          nelem = xlo.get_size();
          ngrid = ear.size();
          npts = ngrid - 1;
          ifl = 1;

          // convert to 32-byte float
          std::vector<FloatArrayType> fear(ngrid);
          CONVERTARRAY(ear, fear, ngrid);

          // Number of bins to send to XSPEC
          nout = ngrid;
          if (xhi) nout--;

          create_output(nout, result, error);

          // Swallow exceptions here

          call_xspec(fear, pars, npars, filename, tabtype, result);

          ngaps = gaps_index.size();
          if (ngaps > 0) {
            finalize_grid(nelem, result, gaps_index);
          }

        } // eval

        const char* get_err_msg( ) { return "XSPEC model evaluation failed"; }

      protected:
        int npts, ifl;
        FloatArray error;

      private:
        int nelem, nout, ngrid, ngaps;
        std::vector<SherpaFloat> ear;
        std::vector<int> gaps_index;

      }; // class xspecTableModelBase


      class xspecTableModelTabint : public xspecTableModelBase {
      public:
        xspecTableModelTabint( ) : xspecTableModelBase(1) { }
        void call_xspec( std::vector<FloatArrayType>& fear,
                         FloatArray& pars, npy_intp npars,
                         char* filename, char*tabtype,
                         FloatArray& result ) {
          tabint( &fear[0], this->npts, &pars[0], npars,
                  filename, this->ifl, tabtype,
                  &result[0], &this->error[0] );
        }
      }; // class xspecTableModelTabint

      //
      // In theory the following should work but have not tested, yet
      //
      // template<XSpecFuncVal XSpecFunc>
      // class xspecTableModel : public xspecTableModelBase {
      // public:
      //   xspecTableModel( ) : xspecTableModelBase(1) { }
      //   void call_xspec( std::vector<FloatArrayType>& fear,
      //                    FloatArray& pars, npy_intp npars,
      //                    int char* filename, char*tabtype,
      //                    FloatArray& result ) {
      //     XSpecFunc( &fear[0], this->npts, &pars[0], filename,
      //                this->ifl, &result[0], &this->error[0] );

      //   }

      // }; // class xspecTableModel


template <npy_intp NumPars, bool HasNorm, xsf77Call XSpecFunc>
PyObject* xspecmodelfct( PyObject* self, PyObject* args ) {

  xspecModelFctF<float, FloatArray, XSpecFunc> xspec_model =
    xspecModelFctF<float, FloatArray, XSpecFunc>( args, NumPars, HasNorm );
  try {
    FloatArray result;
    xspec_model.eval( result );
    return result.return_new_ref();
  } catch(const NoError& re) {
    return NULL;
  } catch(const ValueError& re) {
    PyErr_SetString( PyExc_ValueError, re.what() );
    return NULL;
  } catch(const TypeError& re) {
    PyErr_SetString( PyExc_TypeError, re.what() );
    return NULL;
  } catch(...) {
    // Even though the XSPEC model function is Fortran, it could call
    // C++ functions, so swallow exceptions here
    PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
    return NULL;
  }

}

template <npy_intp NumPars, bool HasNorm, xsF77Call XSpecFunc>
PyObject* xspecmodelfct_dbl( PyObject* self, PyObject* args ) {

  xspecModelFctFD<double, DoubleArray, XSpecFunc> xspec_model =
    xspecModelFctFD<double, DoubleArray, XSpecFunc>( args, NumPars, HasNorm );
  try {
    DoubleArray result;
    xspec_model.eval( result );
    return result.return_new_ref();
  } catch(const NoError& re) {
    return NULL;
  } catch(const ValueError& re) {
    PyErr_SetString( PyExc_ValueError, re.what() );
    return NULL;
  } catch(const TypeError& re) {
    PyErr_SetString( PyExc_TypeError, re.what() );
    return NULL;
  } catch(...) {
    // Even though the XSPEC model function is Fortran, it could call
    // C++ functions, so swallow exceptions here
    PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
    return NULL;
  }

}

template <npy_intp NumPars, bool HasNorm, xsccCall XSpecFunc>
PyObject* xspecmodelfct_C( PyObject* self, PyObject* args ) {

  xspecModelFctC<double, DoubleArray, XSpecFunc> xspec_model =
    xspecModelFctC<double, DoubleArray, XSpecFunc>( args, NumPars, HasNorm );
  try {
    DoubleArray result;
    xspec_model.eval( result );
    return result.return_new_ref();
  } catch(const NoError& re) {
    return NULL;
  } catch(const ValueError& re) {
    PyErr_SetString( PyExc_ValueError, re.what() );
    return NULL;
  } catch(const TypeError& re) {
    PyErr_SetString( PyExc_TypeError, re.what() );
    return NULL;
  } catch(...) {
    PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
    return NULL;
  }

}


// Handle convolution models, which are assumed to have C-style
// linkage.
//
// This template does not support non-contiguous grids.
template <npy_intp NumPars, xsccCall XSpecFunc>
PyObject* xspecmodelfct_con( PyObject* self, PyObject* args ) {

  xspecModelFctConC<double, DoubleArray, XSpecFunc> xspec_model =
    xspecModelFctConC<double, DoubleArray, XSpecFunc>( args, NumPars );
  try {
    DoubleArray result;
    xspec_model.eval(result);
    return result.return_new_ref();
  } catch(const NoError& re) {
    return NULL;
  } catch(const ValueError& re) {
    PyErr_SetString( PyExc_ValueError, re.what() );
    return NULL;
  } catch(const TypeError& re) {
    PyErr_SetString( PyExc_TypeError, re.what() );
    return NULL;
  } catch(...) {
    PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
    return NULL;
  }

}

// As there's only one FORTRAN convolution model, explicitly include
// F77 in the name (rather than have the FORTRAN interface be "default"
// version as it for the additive and multiplicative models).
//
template <npy_intp NumPars, xsf77Call XSpecFunc>
PyObject* xspecmodelfct_con_f77( PyObject* self, PyObject* args ) {

  xspecModelFctConF<float, FloatArray, XSpecFunc> xspec_model =
    xspecModelFctConF<float, FloatArray, XSpecFunc>( args, NumPars );
  try {
    FloatArray result;
    xspec_model.eval( result );
    return result.return_new_ref();
  } catch(const NoError& re) {
    return NULL;
  } catch(const ValueError& re) {
    PyErr_SetString( PyExc_ValueError, re.what() );
    return NULL;
  } catch(const TypeError& re) {
    PyErr_SetString( PyExc_TypeError, re.what() );
    return NULL;
  } catch(...) {
    // Even though the XSPEC model function is Fortran, it could call
    // C++ functions, so swallow exceptions here
    PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
    return NULL;
  }

}

// As of XSPEC 12.10.1, the table-model routines have been
// consolidated into one routine, so there is no need for
// a template. A templace could be used to allow compile-time
// specialization over additive versus multiplicative, but
// for now have a run-time check rather than multiple versions
// of this routine.
//
#ifdef XSPEC_12_10_1

PyObject* xspectablemodel( PyObject* self, PyObject* args, PyObject *kwds )
{

        PyArgTupleBase();

	FloatArray pars;
	DoubleArray xlo;
	DoubleArray xhi;
	char *filename, *tabtype;
	static char *kwlist[] = {(char*)"pars", (char*)"xlo", (char*)"xhi",
                                 (char*)"filename", (char*)"tabtype", NULL};

        // The grid arrays could be cast to FloatArray here, saving
        // conversion later on in this routine. However, that can then
        // lead to differences in the identification of non-contiguous
        // bins, or whether a grid is monotonic and non-overlapping [*]
        // (e.g. if a source expression contains both a FORTRAN
        // and C style model), so stick to this approach for now.
        //
        // [*] although these checks are currently commented out
        //
	if ( !PyArg_ParseTupleAndKeywords( args, kwds, (char*)"O&O&|O&ss",
                                           kwlist,
			(converter)convert_to_contig_array< FloatArray >,
			&pars,
			(converter)convert_to_contig_array< DoubleArray >,
			&xlo,
			(converter)convert_to_contig_array< DoubleArray >,
			&xhi,
                                           &filename,
                                           &tabtype) )
          return NULL;

        // This used to be specified at compile time, but with XSPEC 12.10.1
        // it has been changed to run time.
        bool HasNorm = strcmp(tabtype, "add") == 0;

        // Remove the final parameter if this is an additive model
        npy_intp npars = pars.get_size();
        if (HasNorm) { npars -= 1; }

        xspecTableModelTabint xspec_model = xspecTableModelTabint( );
        try {
          FloatArray result;
          xspec_model.eval(HasNorm, xlo, xhi, npars, pars, filename, tabtype, result);
          // Apply normalization if required (note that npars
          // has been reduced by 1 if HasNorm is true, so it is
          // correct to use npars and not 'npars - 1' here).
          //
          if ( HasNorm )
            for (int i = 0; i < xlo.get_size(); i++)
              result[i] *= pars[npars];
          return result.return_new_ref();
	} catch(const NoError& re) {
	  return NULL;
	} catch(const ValueError& re) {
	  PyErr_SetString( PyExc_ValueError, re.what() );
	  return NULL;
	} catch(const TypeError& re) {
	  PyErr_SetString( PyExc_TypeError, re.what() );
	  return NULL;
        } catch(...) {
          // Even though the XSPEC model function is Fortran, it could call
          // C++ functions, so swallow exceptions here
          PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
          return NULL;
        }

}

#else

typedef void (*XSpecFuncVal)( float* ear, int ne, float* param, const char* filenm, int ifl, float* photar, float* photer );

template <bool HasNorm, XSpecFuncVal XSpecFunc>
PyObject* xspectablemodel( PyObject* self, PyObject* args, PyObject *kwds ) {

        PyArgTupleBase();

	FloatArray pars;
	DoubleArray xlo;
	DoubleArray xhi;
	char *filename;
	static char *kwlist[] = {(char*)"pars", (char*)"xlo", (char*)"xhi",
			(char*)"filename", NULL};

        // The grid arrays could be cast to FloatArray here, saving
        // conversion later on in this routine. However, that can then
        // lead to differences in the identification of non-contiguous
        // bins, or whether a grid is monotonic and non-overlapping [*]
        // (e.g. if a source expression contains both a FORTRAN
        // and C style model), so stick to this approach for now.
        //
        // [*] although these checks are currently commented out
        //
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

        // xspecTableModel<XSpecFunc> xspec_model = xspecTableModel<XSpecFunc>( );
        // try {
        //   FloatArray result;
        //   xspec_model.eval(HasNorm, xlo, xhi, npars, pars, filename, NULL, result);
        //   // Apply normalization if required
        //   if ( HasNorm )
        //     for (int i = 0; i < xlo.get_size(); i++)
        //       result[i] *= pars[npars - 1];   // NOTE: NumPars not sent to template
        //   return result.return_new_ref();
        // } catch(std::runtime_error& re) {
        //   return NULL;
        // } catch(...) {
        //   // Even though the XSPEC model function is Fortran, it could call
        //   // C++ functions, so swallow exceptions here
        //   PyErr_SetString( PyExc_ValueError, xspec_model.get_err_msg() );
        //   return NULL;
        // }

        //
	// The grid to send to XSPEC (double precision).
	//
	std::vector<SherpaFloat> ear;
        std::vector<int> gaps_index;
	create_grid(xlo, xhi, ear, gaps_index);

	int nelem = int( xlo.get_size() );
	int ngrid = ear.size();
	int ifl = 1;

        // convert to 32-byte float
        std::vector<FloatArrayType> fear(ngrid);
	CONVERTARRAY(ear, fear, ngrid);

	// Number of bins to send to XSPEC
	int nout = ngrid;
	if (xhi) nout--;

	FloatArray result, error;
	create_output(nout, result, error);

	// Swallow exceptions here

        try {

          int npts = ngrid - 1;
          XSpecFunc( &fear[0], npts, &pars[0], filename, ifl,
                     &result[0], &error[0] );

	} catch(...) {

          PyErr_SetString( PyExc_ValueError,
                           (char*)"XSPEC model evaluation failed" );
          return NULL;

	}

	int ngaps = (int) gaps_index.size();
	if (ngaps > 0) {
	  finalize_grid(nelem, result, gaps_index);
	}

	// Apply normalization if required
	if ( HasNorm )
          for (int i = 0; i < nelem; i++)
            result[i] *= pars[npars - 1];  // NOTE: NumPars not sent to template

	return result.return_new_ref();

}

#endif

} } } /* namespace xspec, namespace astro, namespace sherpa */


#define _XSPECFCTSPEC(name, npars, has_norm) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct< npars, has_norm, \
				name##_ >))

#define XSPECMODELFCT(name, npars)  _XSPECFCTSPEC(name, npars, false)
#define XSPECMODELFCT_NORM(name, npars)  _XSPECFCTSPEC(name, npars, true)

// double precision
#define XSPECMODELFCT_DBL(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_dbl< npars, false, \
				name##_ >))

#define XSPECMODELFCT_C(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, false, name >))

#define XSPECMODELFCT_C_NORM(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_C< npars, true, name >))

#define XSPECMODELFCT_CON(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_con< npars, name >))

#define XSPECMODELFCT_CON_F77(name, npars) \
		FCTSPEC(name, (sherpa::astro::xspec::xspecmodelfct_con_f77< npars, name##_ >))


#ifdef XSPEC_12_10_1
#define XSPECTABLEMODEL        \
		{ (char*)"tabint", \
	(PyCFunction)((PyCFunctionWithKeywords)sherpa::astro::xspec::xspectablemodel), \
	METH_VARARGS|METH_KEYWORDS, \
	NULL }

#else

#define _XSPECTABLEMODELSPEC(name, has_norm) \
		{ (char*)#name, \
	(PyCFunction)((PyCFunctionWithKeywords)sherpa::astro::xspec::xspectablemodel< has_norm, name >), \
	METH_VARARGS|METH_KEYWORDS, \
	NULL }

#define XSPECTABLEMODEL(name) \
		_XSPECTABLEMODELSPEC(name, false)

#define XSPECTABLEMODEL_NORM(name) \
		_XSPECTABLEMODELSPEC(name, true)

#endif

#endif /* __sherpa_astro_xspec_extension_hh__ */
