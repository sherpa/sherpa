//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#ifndef __sherpa_info_matrix_hh__
#define __sherpa_info_matrix_hh__
int info_matrix(double* original_pars, const int op_size,
		const double* pars_mins, const int mins_size,
		const double* pars_maxs, const int maxs_size,
		const double* pars_hardmins, const int hmins_size,
		const double* pars_hardmaxs, const int hmaxs_size,
		double* info, const int info_rows, const int info_cols, 
		const double sigma, 
		const double eps, 
		const int maxiters,
		const double remin,
		double (*fcn)(double*, int)) throw();
#endif
