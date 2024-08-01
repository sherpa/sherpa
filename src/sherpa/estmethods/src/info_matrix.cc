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

#include <cstdlib>
#include <vector>
#include <float.h>

#include "estutils.hh"


// This function calculates the information matrix--*not* the 
// covariance matrix.  The calling function has to invert the
// information matrix to get the covariance matrix.

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
			    double (*fcn)(double*, int)) throw()
{
    int iter = 3;
    int i,j,k;
    int numpars = op_size;
    est_return_code status;
    status.status = EST_SUCCESS;
    status.par_number = -1;
    status.nfits = 0;  // Always zero for covariance, never refit

    // The confidence levels we get can't possibly correspond
    // to a sigma less than zero; i.e., if sigma == 0, then
    // the levels are equal to the parameter value itself.
    if (sigma < 0) {
      status.status = EST_FAILURE;
      return status;
    }
    
    double min_stat = fcn(original_pars, numpars);
    if (isnan(min_stat)) {
      status.status = EST_HITNAN;
      return status;
    }

    std::vector<double> h(iter);
    std::vector<double> d2f(iter);
    
    double delta_stat= pow(sigma,2.0);
    double thresh_stat = min_stat + delta_stat;

    std::vector<double> e(numpars);
    std::vector<double> pars(numpars);

    for (i = 0 ; i < numpars; i++)
      pars[i] = original_pars[i];
    
    double pb = 1.0;
    est_return_code s;
    for (i = 0; i < numpars; i++) {
      e[i] = 0.0;
      for ( j = 0 ; j < 2 ; j++ ) {
	s = get_onesided_interval(original_pars, pars_mins,
				  pars_maxs, pars_hardmins,
				  pars_hardmaxs, i,
				  min_stat, thresh_stat, sigma,
				  eps, maxiters, remin,
				  j, numpars,
				  &pb, fcn);
	if (s.status == EST_NEWMIN) {
	  return s;
	}
	if (s.status != EST_SUCCESS) {
	  // If failure for any reason, then just say that e[i]
	  // on this side is the magnitude of the parameter value,
	  // and try to keep going.
	  e[i] += fabs(0.0 - original_pars[i]) / 2.0;
	}
	else {
	  e[i] += fabs(pb-original_pars[i]) / 2.0;
	}
      }
    }
    
    double ratio = 0.707;
    double f1 = 0;
    double f2 = 0;
    
    for (i = 0; i < numpars; i++) {
      int found_nan = 0;
      for ( j = 0 ; j < iter ; j++ ) {
	h[j] = e[i] * pow(ratio,(double)(iter-(j+1)));
	pars[i] = original_pars[i] + h[j];
	// If we end up outside parameter limits, set to be
	// inside limits
	if (pars[i] < pars_hardmins[i]) {
	  pars[i] = pars_hardmins[i];
	}
	if (pars[i] > pars_hardmaxs[i]) {
	  pars[i] = pars_hardmaxs[i];
	}

	if (isnan(h[j]) || isnan(pars[i])) {
	  found_nan = 1;
	  break;
	  //status.status = EST_HITNAN;
	  //return status;
	}
	f1 = fcn(&pars[0], numpars);

	pars[i] = original_pars[i] - h[j];
	// If we end up outside parameter limits, set to be
	// inside limits
	if (pars[i] < pars_hardmins[i]) {
	  pars[i] = pars_hardmins[i];
	}
	if (pars[i] > pars_hardmaxs[i]) {
	  pars[i] = pars_hardmaxs[i];
	}

	if (isnan(h[j]) || isnan(pars[i])) {
	  found_nan = 1;
	  break;
	  //status.status = EST_HITNAN;
	  //return status;
	}
	f2 = fcn(&pars[0], numpars);

	d2f[j] = ( (2*min_stat) - (f1+f2) )/(h[j]*h[j]);
	if (isnan(d2f[j])) {
	  found_nan = 1;
	  break;
	  //status.status = EST_HITNAN;
	  //return status;
	}
      }

      if (found_nan == 1)
	info[i*numpars + i] = -DBL_MAX;
      else {
	if (neville(iter, &h[0], &d2f[0], 0.0, info[i*numpars + i]) !=
	    EXIT_SUCCESS)
	  info[i*numpars + i] = -DBL_MAX;
      }
      info[i*numpars + i] = -info[i*numpars + i];
      pars[i] = original_pars[i];
    }
    
    //
    // Now find the off-diagonal derivatives.
    // 
    for (i = 0; i < numpars; i++) {
      double p1 = original_pars[i];
      for ( j = i+1 ; j < numpars ; j++ ) {
	double p2 = original_pars[j];
	int found_nan = 0;
	for ( k = 0 ; k < iter ; k++ ) {
	  h[k] = pow(ratio,(double)(iter-(k+1)));
	  pars[i] = p1+h[k]/sqrt(info[i*numpars + i]);
	  pars[j] = p2+h[k]/sqrt(info[j*numpars + j]);

	  // If we end up outside parameter limits, set to be
	  // inside limits
	  if (pars[i] < pars_hardmins[i]) {
	    pars[i] = pars_hardmins[i];
	  }
	  if (pars[i] > pars_hardmaxs[i]) {
	    pars[i] = pars_hardmaxs[i];
	  }
	  if (pars[j] < pars_hardmins[j]) {
	    pars[j] = pars_hardmins[j];
	  }
	  if (pars[j] > pars_hardmaxs[j]) {
	    pars[j] = pars_hardmaxs[j];
	  }

	  if (isnan(h[k]) || isnan(pars[i]) || isnan(pars[j])) {
	    found_nan = 1;
	    break;
	    //status.status = EST_HITNAN;
	    //return status;
	  }
	  f1 = fcn(&pars[0], numpars);
	  
	  pars[i] = p1-h[k]/sqrt(info[i*numpars + i]);
	  pars[j] = p2-h[k]/sqrt(info[j*numpars + j]);
	  // If we end up outside parameter limits, set to be
	  // inside limits
	  if (pars[i] < pars_hardmins[i]) {
	    pars[i] = pars_hardmins[i];
	  }
	  if (pars[i] > pars_hardmaxs[i]) {
	    pars[i] = pars_hardmaxs[i];
	  }
	  if (pars[j] < pars_hardmins[j]) {
	    pars[j] = pars_hardmins[j];
	  }
	  if (pars[j] > pars_hardmaxs[j]) {
	    pars[j] = pars_hardmaxs[j];
	  }

	  if (isnan(h[k]) || isnan(pars[i]) || isnan(pars[j])) {
	    found_nan = 1;
	    break;
	    //status.status = EST_HITNAN;
	    //return status;
	  }
	  f2 = fcn(&pars[0], numpars);
	  d2f[k] = ( (2*min_stat) - (f1+f2) )/(h[k]*h[k]);
	  if (isnan(d2f[k])) {
	    found_nan = 1;
	    break;
	    //status.status = EST_HITNAN;
	    //return status;
	  }
	}

	if (found_nan == 1)
	  info[i*numpars + j] = -DBL_MAX;
	else {
	  if (neville(iter, &h[0], &d2f[0], 0.0, info[i*numpars + j]) !=
	      EXIT_SUCCESS)
	    info[i*numpars + j] = -DBL_MAX;
	}
	info[i*numpars + j] = -(info[i*numpars + j]+2)*
	  sqrt(info[i*numpars + i]*info[j*numpars + j])/2.;
	info[j*numpars + i] = info[i*numpars + j];
	pars[j] = p2;
      }
      pars[i] = p1;
    }
    
    //
    // If the statistic is chi**2, Cash, or Cstat, then the info matrix must be
    // divided by a factor of 2...
    // For now, the only statistics we have, so just leave this.  SMD 04/15/06
    //
    for ( i = 0 ; i < numpars ; i++ ) {
      for ( j = 0 ; j < numpars ; j++ ) {
	info[i*numpars + j] /= 2.;
      }
    }
    
    // Return the information matrix.  It needs to be inverted
    // by the user to get the covariance matrix.  SMD 04/26/06
    return status;
} 
