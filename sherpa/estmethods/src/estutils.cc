//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#include <cstdlib>
#include <vector>

#include "estutils.hh"
#include "sherpa/utils.hh"
#include "fcmp.h"

// These are utility functions used by the Sherpa methods
// for estimating parameter uncertainties.  The most important
// function here is get_onesided_interval.  It is used by
// both the covariance and projection methods to set a scale
// for the region of parameter space explored by these methods.

// The neville interpolation code should live in utils dir but
// because the implementation plan the code cannot be modular.
int neville( int n, const double *x, const double *y, double xinterp,
	     double& answer ) throw() {

  std::vector<double> p( y, y+n );
  
  for ( int jj = 1; jj < n; jj++ )
    for ( int ii = n-1; ii >= jj; ii-- ) {
      double denom = x[ii] - x[ii-jj];
      if ( 0.0 == denom )
	return EXIT_FAILURE;
      p[ii] = ( (xinterp-x[ii-jj])*p[ii] - (xinterp-x[ii])*p[ii-1] ) / denom;
    }

  answer = p[n-1];
  return EXIT_SUCCESS;

}

// Often, these parameter estimiation routines need to check
// if they are butting up against the hard minimum or hard
// maximum.  This helper function checks for that.  Make sure
// hard minimum and hard maximum are passed to this function!

int at_param_space_bound(double* par, 
			 const double hardmin, 
			 const double hardmax) throw()
{
  if (!(*par > hardmin)) {
    *par = hardmin;
    return EST_HARDMIN;
  }
  else if (!(*par < hardmax)) {
    *par = hardmax;
    return EST_HARDMAX;
  }
  return EST_SUCCESS;
}

// Wrap up previous function, when we want to set value, and
// make sure it stays inside bounds, without returning error
// status that setting the value exceeded bounds.
void set_value(double* par, 
	       const double parmin, 
	       const double parmax,
	       const double value) throw()
{
  *par = value;
  at_param_space_bound(par, parmin, parmax);
}

// Whether we step above or below the current value depends on
// whether we are searching for the upper or lower estimated 
// error.  In this function, set the parameter value accordingly.

void set_value_from_step(double* par, const double parmin, 
			 const double parmax,const double initv,
			 const double f,const int upper) throw()
{
  if ( upper == 0 && f > 0. ){
    *par = initv - f;
  } else if ( upper == 0 && f < 0. ) {
    *par = initv + f;
  } else if ( upper ==  1 && f > 0. ) {
    *par = initv + f;
  } else if ( upper ==  1 && f < 0. ) {
    *par = initv - f;
  }
  set_value(par, parmin, parmax, *par);
}

double get_stat(double* new_min_stat, double* new_min_parval, 
		const int parnum, const double* pars, 
		const double* pars_mins, const double* pars_maxs,
		const int numpars, double (*fcn)(double*, int))
{
  double return_stat = fcn((double*)pars, numpars);
  if (return_stat < *new_min_stat) {
    if (!(pars[parnum] < pars_mins[parnum]) &&
	!(pars[parnum] > pars_maxs[parnum])) {
      *new_min_parval = pars[parnum];
      *new_min_stat = return_stat;
    }
  }
  return return_stat;
}

// This function sets a scale for parameter error estimate functions.
// If upper is true, set a scale for the upper error; if upper is
// false, set a scale for the lower error.
//
// Arrays pars, pars_mins, pars_maxs, pars_hardmins and pars_hardmaxs
// all need to be the same size--double arrays of length numpars.
//
// numpars is the number of parameters in these arrays; parnum is
// the *position* in these arrays.
//
// This function should be called twice for each parameter in pars.
// Thus if pars holds five parameters, this function will be called
// ten times--calculating an upper and lower scale for each parameter.
//
// par_bound is the address of a double.  par_bound holds either
// the upper (if upper is true) or lower (if upper is false) scale
// determined by this function.
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
				      double (*fcn)(double*, int)) throw()
{
  double frac;
  double deriv,cur_stat,new_stat;
  int at_boundary = EST_SUCCESS;
  double f = 1.;
  double diff = 0.;
  int iters = 0;
  est_return_code status;
  status.status = EST_SUCCESS;
  status.par_number = -1;
  status.nfits = 0;  // Should remain zero as fitting not called by
                     // this function

  if (sigma < 0 || eps < 0) {
    status.status = EST_FAILURE;
    return status;
  }

  double epshi = POW((sigma+eps)/2.,2.);
  double epslo = POW((sigma-eps)/2.,2.);
  
  double epsilon = FABS(epshi-epslo);
  
  double initv = pars[parnum];

  double new_min_stat = min_stat;
  double new_min_parval = initv;
  
  if ( initv != 0. ) {
    frac = 0.01 * initv;
  } else {
    if ( !upper ) {
      if ( pars_mins[parnum] < 0. ) {
	frac = (pars_mins[parnum])/100.;
      } else {
	at_boundary = at_param_space_bound(pars+parnum,
					   pars_hardmins[parnum],
					   pars_hardmaxs[parnum]);
	frac = 1.;
      }
    } else {
      if ( pars_maxs[parnum] > 0. ) {
	frac = (pars_maxs[parnum])/100.;
      } else {
	at_boundary = at_param_space_bound(pars+parnum,
					   pars_hardmins[parnum],
					   pars_hardmaxs[parnum]);
	frac = 1.;
      }
    }
  }
  
  if ( at_boundary == EST_SUCCESS ) {
    set_value_from_step(pars+parnum,pars_hardmins[parnum],
			pars_hardmaxs[parnum],initv,f*frac,upper);
    if ( initv == 0. ) {
      while (1) {
	new_stat = get_stat(&new_min_stat, &new_min_parval, 
			    parnum, pars, pars_mins, pars_maxs,
			    numpars, fcn);
	if ( new_stat > 1.2*thr_stat ) {
	  frac /= 10.;
	  set_value_from_step(pars+parnum,pars_hardmins[parnum],
			      pars_hardmaxs[parnum],initv,f*frac,upper);
	} else {
	  break;
	}
	if ( FABS(frac) < 1.e-30 ) break;
      }
    }
  }
  
  at_boundary = at_param_space_bound(pars+parnum, pars_hardmins[parnum],
				     pars_hardmaxs[parnum]);
  
  if ( at_boundary == EST_SUCCESS ) {
    new_stat = get_stat(&new_min_stat, &new_min_parval, 
			parnum, pars, pars_mins, pars_maxs,
			numpars, fcn);
    diff = new_stat - min_stat;
    while ( diff <= 0. ) {
      f = 1.3 * f;
      set_value_from_step(pars+parnum,pars_hardmins[parnum],
			  pars_hardmaxs[parnum],initv,f*frac,upper);
      at_boundary = at_param_space_bound(pars+parnum, pars_hardmins[parnum],
					 pars_hardmaxs[parnum]);
      if ( at_boundary != EST_SUCCESS ) break;
      new_stat = get_stat(&new_min_stat, &new_min_parval, 
			  parnum, pars, pars_mins, pars_maxs,
			  numpars, fcn);
      diff = new_stat - min_stat;
    }
  }
  
  if ( at_boundary == EST_SUCCESS ) {
    set_value_from_step(pars+parnum,pars_hardmins[parnum],
			pars_hardmaxs[parnum],initv,
			f*frac*sqrt(FABS(thr_stat-min_stat)/diff),
			upper);
    at_boundary = at_param_space_bound(pars+parnum, pars_hardmins[parnum],
				       pars_hardmaxs[parnum]);
  }
  
  //
  // Now iterate to the final solution.
  //
  // Note that in the following, (float)upper should not be
  // necessary because all that we are doing is estimating
  // the local slope (derivative) of the statistic, which
  // can be estimated either "to the left" or "to the right"
  // of the current point, regardless of its position relative
  // to the mode (i.e. we estimate "to the right" only).
  //
  if ( at_boundary == EST_SUCCESS ) {
    while ( iters < maxiters ) {
      cur_stat = get_stat(&new_min_stat, &new_min_parval, 
			  parnum, pars, pars_mins, pars_maxs,
			  numpars, fcn);
      
      if ( FABS(thr_stat - cur_stat) > epsilon ) {
	double step = FABS(pars[parnum]-initv)/10.;
	double test_bound_val = pars[parnum] + step;
	at_boundary = at_param_space_bound(&test_bound_val, 
					   pars_hardmins[parnum],
					   pars_hardmaxs[parnum]);
	if ( at_boundary != EST_SUCCESS ) break;
	pars[parnum] += step;
	new_stat = get_stat(&new_min_stat, &new_min_parval, 
			    parnum, pars, pars_mins, pars_maxs,
			    numpars, fcn);
	deriv = (new_stat - cur_stat)/step;
	test_bound_val = pars[parnum] +
	  (thr_stat-cur_stat)/deriv - step;
	
	if ( upper == 0 && test_bound_val >= initv ) {
	  test_bound_val = initv - step;
	} else if ( upper == 1 && test_bound_val <= initv ) {
	  test_bound_val = initv + step;
	}
	
	at_boundary = at_param_space_bound(&test_bound_val, 
					   pars_hardmins[parnum],
					   pars_hardmaxs[parnum]);
	if ( at_boundary != EST_SUCCESS ) break;
	pars[parnum] = test_bound_val;
      } else {
	break;
      }
      iters++;
    }
  }
  
  // Instead of optimizing model parameters here, in C scope, return
  // the EST_NEWMIN status code, fit in Python scope when the
  // exception is caught, and then call est_method again.  Have to do
  // that because the only function passed into this function is the
  // fit statistic function.  The difference between the old and new
  // minima has to be greater than the "remin" tolerance value.
  //
  // And it really only makes sense to proceed if the new "best value"
  // is within the soft limits; if you are at a soft limit, it's
  // pretty clear that parameter space is not going to constrain
  // your limits very well.
  //
  // If remin flag zero or less, don't indicate reminimization -- 
  // prevent triggering when new stat is *greater* than minimum (should
  // always be greater than minimum anyway).

  if (remin > 0.0 &&
      _sao_fcmp(min_stat, new_min_stat, remin) > 0 &&
      pars[parnum] < pars_maxs[parnum] &&
      pars[parnum] > pars_mins[parnum]) {
    pars[parnum] = new_min_parval;
    if (upper == 1)
      *par_bound = 1.0;
    else
      *par_bound = -1.0;
    new_stat = fcn(pars, numpars);
    status.status = EST_NEWMIN;
    return status;
  }

  if ( at_boundary == EST_SUCCESS && iters < maxiters ) {  // good
    *par_bound = pars[parnum];
    status.status = EST_SUCCESS;
  } else if ( at_boundary == EST_SUCCESS ) {                // not good
    if (upper == 1)
      *par_bound = pars_hardmaxs[parnum] - initv;
    else
      *par_bound = pars_hardmins[parnum] - initv;
    status.status = EST_MAXITER;
  } else {                                       // boundary hit
    if (upper == 1)
      *par_bound = pars_hardmaxs[parnum] - initv;
    else
      *par_bound = pars_hardmins[parnum] - initv;
    status.status = at_boundary;
    status.par_number = parnum;
  }

  pars[parnum] = initv;
  return status;
}
