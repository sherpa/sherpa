//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#ifndef __sherpa_projection_hh__
#define __sherpa_projection_hh__
int projection(double* original_pars, const int op_size,
	       const double* pars_mins, const int mins_size,
	       const double* pars_maxs, const int maxs_size,
	       const double* pars_hardmins, const int hmins_size,
	       const double* pars_hardmaxs, const int hmaxs_size,
	       double* pars_elow, const int elow_size,
	       double* pars_ehi, const int ehi_size,
	       const double sigma, 
	       const double eps, 
	       const int maxiters,
	       const double remin,
	       double (*statfcn)(double*, int),
	       double (*fitfcn)(double (*statfcn)(double*, int),
				double*,double*,double*,int,int)) throw();
#endif
