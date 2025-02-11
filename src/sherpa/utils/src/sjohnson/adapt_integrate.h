#ifndef __adapt_integrate_h__
#define __adapt_integrate_h__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>

typedef double (*integrand) ( unsigned ndim, const double *x, void *);

int adapt_integrate(integrand f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, 
		    double reqAbsError, double reqRelError, 
		    double *val, double *estimated_error);

#ifdef __cplusplus
}
#endif

#endif /* __adapt_integrate_h__ */
