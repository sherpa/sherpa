//_C++_INSERT_SAO_COPYRIGHT_HERE_(2008)_
//_C++_INSERT_GPL_LICENSE_HERE_

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
