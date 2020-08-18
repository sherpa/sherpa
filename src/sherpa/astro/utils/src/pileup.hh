// 
//  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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


#ifndef __sherpa_pileup_hh__
#define __sherpa_pileup_hh__

#include "PyWrapper.hh"

#define PILEUP_VERSION "0.1.1"

#define MAX_NUM_TERMS	30

int apply_pileup(unsigned int num_bins, const double *arf_source_ftime,
		 double *results, double *pileup_fractions, double *integral_ae,
		 double exposure_time,
		 unsigned int max_num_terms, unsigned int *num_terms,
		 const double *energ_lo, const double *energ_hi,
		 const double *specresp, double fracexpo, double frame_time,
		 double alpha, double g0, double num_regions, double psf_frac,
		 sherpa::usrfuncproto model_func, sherpa::PyWrapper* x);

#endif // __sherpa_pileup_hh__
