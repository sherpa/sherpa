// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_estutils_hh__
#define __sherpa_estutils_hh__
#include <math.h>
#ifdef __SUNPRO_CC
#include <sunmath.h>
#define NAN quiet_nan(0)
#endif

#define EST_SUCCESS    0
#define EST_FAILURE    1
#define EST_HARDMIN    2
#define EST_HARDMAX    3
#define EST_HARDMINMAX 4
#define EST_NEWMIN     5
#define EST_MAXITER    6
#define EST_HITNAN     7
//what about maxiters? nans?

// The intent of this structure is to return success, or
// an indicator of the kind of failure; and, if a failure,
// the parameter number for which the failure occurred.
struct est_return_code {
  int status;
  int par_number;
  int nfits;
};

int neville( int n, const double *x, const double *y, double xinterp,
	     double& answer ) throw();
int at_param_space_bound(const double par, 
			 const double hardmin, 
			 const double hardmax) throw();
void set_value(double* par, const double parmin, const double parmax, 
	       const double value) throw();
void set_value_from_step(double* par,const double parmin, 
			 const double parmax,const double initv,
			 const double f,const int upper) throw();
est_return_code get_onesided_interval(double* pars, 
				      const double* pars_mins,
				      const double* pars_maxs, 
				      const double* pars_hardmins,
				      const double* pars_hardmaxs, 
				      const int parnum,
				      const double min_stat, 
				      const double thr_stat,
				      const double sigma,
				      const double eps, 
				      const int maxiters,
				      const double remin,
				      const int upper, 
				      const int numpars,
				      double* par_bound,
				      double (*fcn)(double*, int)) throw();


est_return_code info_matrix(double* original_pars, const int op_size,
			    const double* pars_mins, const int mins_size,
			    const double* pars_maxs, const int maxs_size,
			    const double* pars_hardmins, const int hmins_size,
			    const double* pars_hardmaxs, const int hmaxs_size,
			    double* info, const int info_rows, 
			    const int info_cols, 
			    const double sigma, 
			    const double eps, 
			    const int maxiters,
			    const double remin,
			    double (*fcn)(double*, int)) throw();

est_return_code projection(double* original_pars, const int op_size,
			   const double* pars_mins, const int mins_size,
			   const double* pars_maxs, const int maxs_size,
			   const double* pars_hardmins, const int hmins_size,
			   const double* pars_hardmaxs, const int hmaxs_size,
			   double* pars_elow, const int elow_size,
			   double* pars_ehi, const int ehi_size,
			   int* pars_eflags, const int eflags_size,
			   const double sigma, 
			   const double eps, 
			   const double tol,
			   const int maxiters,
			   const double remin,
			   const int* parnums, const int parnumsize,
			   double (*statfcn)(double*, int),
			   double (*fitfcn)(double (*statfcn)(double*, int),
					    double*,double*,double*,
					    int,int)) throw();


#endif
