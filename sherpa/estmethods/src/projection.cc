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
#include <float.h>
#include <vector>
#include "estutils.hh"

static double minimize(double* pars, double* parmins, 
		double* parmaxs, int numpars, int parnum,
		double (*statfcn)(double*, int),
		double (*fitfcn)(double (*statfcn)(double*, int),
				 double*,double*,double*,int,int)) throw()
{
  double return_val = NAN;
  std::vector<double> orig_pars(numpars);

  for (int i = 0; i < numpars; i++)
    orig_pars[i] = pars[i];

  return_val = fitfcn(statfcn, pars, parmins, parmaxs, numpars, parnum);

  for (int i = 0; i < numpars; i++)
    pars[i] = orig_pars[i];

  return return_val;
}

static double make_projection(double* pars, const double* pars_hardmins,
		       const double* pars_hardmaxs, int numpars,
		       int parnum,double pstep,double chisq,double tol,
		       int* nfits,
		       double (*statfcn)(double*, int),
		       double (*fitfcn)(double (*statfcn)(double*, int),
					double*,double*,double*,int,int)) throw()
{
  // Problem:  In old Sherpa, this is where parameter parnum
  // was frozen, and the bounds were set to hard mins and
  // maxs.  How do we callback to fit method to accomplish
  // these things?
  std::vector<double> pars_new(numpars);
  std::vector<double> pars_newmins(numpars);
  std::vector<double> pars_newmaxs(numpars);

  for (int i = 0; i < numpars; i++) {
    pars_new[i] = pars[i];
    pars_newmins[i] = pars_hardmins[i];
    pars_newmaxs[i] = pars_hardmaxs[i];
  }
  // in minimize, set pars_newmins[parnum], pars_newmaxs[parnum]
  // to pnew, so that value is frozen there during new fit.
  
  double dpeps = 0.001;
  short itercount = 0;
  short ichange = 0;
  short istop = 0;
  double dp = 0;
  double fold = 0;
  double f = 0;
  double proj = 0.0;
  int hardbound = 0;

  if (pstep < 0.0)
    proj = NAN;
  else
    proj = NAN;

  do {
    if ( itercount ) fold = f;
    set_value(&pars_new[parnum], pars_newmins[parnum],
	      pars_newmaxs[parnum], pars[parnum]+exp(dp)*pstep);
    f = minimize(&pars_new[0], &pars_newmins[0], &pars_newmaxs[0],
		 numpars, parnum, statfcn, fitfcn);
    *nfits = *nfits + 1;
    if (isnan(f)) {
      return NAN;
    }
    if ( itercount ) {
      if ( fold < chisq && f >= chisq ) {
	istop = 1;
      } else if ( fold >= chisq && f < chisq ) {
	istop = 1;
      }
    }
    if ( pars_new[parnum] <= pars_hardmins[parnum] || 
	 pars_new[parnum] >= pars_hardmaxs[parnum] ) {
      if ( f < chisq ) {
	if ( pars_new[parnum] <= pars_hardmins[parnum] ) {
	  proj = NAN;
	  // Return code that lower projection bound not found?
	  // Or return -DBL_MAX to indicate a problem
	} else {
	  proj = NAN;
	  // Return code that upper projection bound not found?
	  // Or return DBL_MAX to indicate a problem
	}
	hardbound = 1;
      } else {
	ichange = 1;
      }
      break;
    } 
    if ( f == fold ) {
      // Return code that no change in stat value found?
      // Return hard limit to show that no limit was found
      // (i.e., if you are at the hard limit no way we found
      // reasonable bounds) -- sign of pstep indicates if
      // we are heading for lower bound (negative pstep) or
      // upper bound (positive pstep).
      if (pstep < 0.0)
	proj = NAN;
      else
	proj = NAN;
      hardbound = 1;
      break;
    }
    if ( !istop ) {
      if ( f < chisq ) {
	dp++;
	ichange = 1;
      } else {
	dp--;
	ichange = -1;
      }
    } 
    // If the exponential factor is getting smaller than 
    // DBL_EPSILON, then subtracting (factor * step size)
    // from parameter will bascially be zero -- no need then
    // to try exp(dp) even smaller than DBL_EPSILON
    if (exp(dp) < DBL_EPSILON)
      return NAN;
    itercount++;
  } while ( istop == 0 );
  istop = 0;
  double dplo = 0;
  double dphi = 0;
  double flo = 0;
  double fhi = 0;
  if ( ichange == 1 ) {
    dplo = dp - 1;
    dphi = dp;
    flo = fold;
    fhi = f;
  } else {
    dplo = dp;
    dphi = dp + 1;
    flo = f;
    fhi = fold;
  }
  double dpestold = dplo + (chisq-flo)*(dphi-dplo)/(fhi-flo);
  double dpest = 0;
  if ( !hardbound ) {
    do {
      dp = (dplo + dphi)/2.;
      set_value(&pars_new[parnum], pars_newmins[parnum],
		pars_newmaxs[parnum], pars[parnum]+exp(dp)*pstep);
      f = minimize(&pars_new[0], &pars_newmins[0], &pars_newmaxs[0],
		   numpars, parnum, statfcn, fitfcn);
      *nfits = *nfits + 1;
      if (isnan(f)) {
	return NAN;
      }
      if ( f < chisq && flo < chisq ) {
	flo = f;
	dplo = dp;
      } else if ( f < chisq && fhi < chisq ) {
	fhi = f;
	dphi = dp;
      } else if ( f >= chisq && flo < chisq ) {
	fhi = f;
	dphi = dp;
      } else if ( f >= chisq && fhi < chisq ) {
	flo = f;
	dplo = dp;
      }
      dpest = dplo + (chisq-flo)*(dphi-dplo)/(fhi-flo);
      if ( fabs(dpest-dpestold) < dpeps ) {
	if (tol > 0.0) {
	  set_value(&pars_new[parnum], pars_newmins[parnum],
		    pars_newmaxs[parnum], pars[parnum]+exp(dpest)*pstep);
	  f = minimize(&pars_new[0], &pars_newmins[0], &pars_newmaxs[0],
		       numpars, parnum, statfcn, fitfcn);
	  *nfits = *nfits + 1;
	  if (isnan(f)) {
	    return NAN;
	  }
	  if (fabs(f - chisq) < tol)
	    istop = 1;
	  else
	    return NAN;
	}
	else {
	  istop = 1;
	}
      } 
      else {
	dpestold = dpest;
      }
    } while ( istop == 0 );

    // Return size of error bar rather that upper or
    // lower bound
    // proj = pars[parnum]+exp(dpest)*pstep;
    proj = exp(dpest)*pstep;
  }
  
  return proj;
}

// This function use the projection method to estimate
// errors on best-fit parameter values.

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
			   const int maxiters, const double remin,
			   const int* parnums, const int parnumsize,
			   double (*statfcn)(double*, int),
			   double (*fitfcn)(double (*statfcn)(double*, int),
					    double*,double*,double*,int,
					    int)) throw()
{
  int i,j;
  int numpars = op_size;
  est_return_code status;
  status.status = EST_SUCCESS;
  status.par_number = -1;
  status.nfits = 0;  // For projection, return number of times 
                     // fit called as part of the status

  // The confidence levels we get can't possibly correspond
  // to a sigma less than zero; i.e., if sigma == 0, then
  // the levels are equal to the parameter value itself.
  if (!(sigma > 0)) {
    status.status = EST_FAILURE;
    return status;
  }

  double min_stat = statfcn(original_pars, numpars);
  if (isnan(min_stat)) {
    status.status = EST_HITNAN;
    return status;
  }

  double delta_stat= pow(sigma,2.0);
  double thresh_stat = min_stat + delta_stat;

  double pb = 1.0;
  est_return_code s;
  s.status = EST_SUCCESS;
  s.par_number = -1;

  for (i = 0; i < parnumsize; i++) {
    pars_eflags[i] = EST_SUCCESS;
    for ( j = 0 ; j < 2 ; j++ ) {
      s = get_onesided_interval(original_pars, pars_mins,
				pars_maxs, pars_hardmins,
				pars_hardmaxs, parnums[i],
				min_stat, thresh_stat, sigma,
				eps, maxiters, remin,
				j, numpars,
				&pb, statfcn);
      if (s.status == EST_NEWMIN) {
	status.status = s.status;
	status.par_number = s.par_number;
	return status;
      }
      //
      // If only one free parameter, then use bounds
      // determined by get_onesided_interval method.
      if (numpars == 1) {
	if (j == 0) {
	  if (s.status != EST_SUCCESS ||
	      !(pb > pars_hardmins[parnums[i]])) {
	    pars_elow[i] = NAN;
	    pars_eflags[i] = EST_HARDMIN;
	  }
	  else {
	    pars_elow[i] = pb - original_pars[parnums[i]];
	    pars_eflags[i] = EST_SUCCESS;
	  }
	}
	if (j == 1) {
	  if (s.status != EST_SUCCESS ||
	      !(pb < pars_hardmaxs[parnums[i]])) {
	    pars_ehi[i] = NAN;
	    if (pars_eflags[i] == EST_HARDMIN)
	      pars_eflags[i] = EST_HARDMINMAX;
	    else
	      pars_eflags[i] = EST_HARDMAX;
	  }
	  else {
	    pars_ehi[i] = pb - original_pars[parnums[i]];
	  }
	}
      }
      else if (numpars > 1) {
	double scale;
	if (j == 0) {
	  if (s.status != EST_SUCCESS) {
	    pars_elow[i] = NAN;
	    pars_eflags[i] = EST_HARDMIN;
	  }
	  else {
	    scale = pb - original_pars[parnums[i]];
	    
	    pars_elow[i] = make_projection(original_pars,
					   pars_hardmins,
					   pars_hardmaxs,
					   numpars, parnums[i],
					   scale,
					   thresh_stat,
					   tol,
					   &(status.nfits),
					   statfcn,
					   fitfcn);
	    if (isnan(pars_elow[i])) {
	      pars_elow[i] = NAN;
	      pars_eflags[i] = EST_HARDMIN;
	    }
	    else {
	      pars_eflags[i] = EST_SUCCESS;
	    }
	  }
	}
	if (j == 1) {
	  if (s.status != EST_SUCCESS) {
	    pars_ehi[i] = NAN;
	    if (pars_eflags[i] == EST_HARDMIN)
	      pars_eflags[i] = EST_HARDMINMAX;
	    else
	      pars_eflags[i] = EST_HARDMAX;
	  }
	  else {
	    scale = pb - original_pars[parnums[i]];

	    pars_ehi[i] = make_projection(original_pars,
					  pars_hardmins,
					  pars_hardmaxs,
					  numpars, parnums[i],
					  scale,
					  thresh_stat,
					  tol,
					  &(status.nfits),
					  statfcn,
					  fitfcn);

	    if (isnan(pars_ehi[i])) {
	      pars_ehi[i] = NAN;
	      if (pars_eflags[i] == EST_HARDMIN)
		pars_eflags[i] = EST_HARDMINMAX;
	      else
		pars_eflags[i] = EST_HARDMAX;
	    }
	  }
	}
      }
    }
  }
  
  return status;
} 
