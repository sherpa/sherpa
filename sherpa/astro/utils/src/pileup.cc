//     Copyright (C) 2000 Massachusetts Institute of Technology 
 
//     Author:  John E. Davis <davis@space.mit.edu>

//     This program is free software; you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details. 

//     You should have received a copy of the GNU General Public License
//     along with this program; if not, write to the Free Software
//     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include "pileup.hh"

// Make this a power of 2 so that ffts are fast
#define NUM_POINTS      (1024*4)
#define XMALLOC(n,t) (t *)malloc((n)*sizeof(t))
#define ISIS_FREE(p) do {if (p) free ((char *)p); p = NULL;} while (0)

static double Max_Probability_Cutoff = 0.99999;

typedef struct {
  unsigned int num_terms;
  unsigned int max_num_terms;
  double num_frames;
  int has_psf_fraction;
  unsigned int num;
  double exposure_time;
  double frame_time;
  double *energies;
  double de;
  double *arf_s;		       // A(E)s(E)
  double *arf_time;
  double *results;
  double *arf_s_fft;
  double *arf_s_tmp;
  double *arf_s_fft_tmp;
  double *pileup_fractions;
  double integral_ae;
  double arf_frac_exposure;
} pileup_kernel_t;


// static void print_pileup_kernel( pileup_kernel_t* k ) {
//   printf( "num_terms=%d\n", k->num_terms );
//   printf( "max_num_terms=%d\n", k->max_num_terms );
//   printf( "num_frames=%g\n", k->num_frames );
//   printf( "has_psf_fraction=%d\n", k->has_psf_fraction );
//   printf( "num=%d\n", k->num );
//   printf( "exposure_time=%g\n", k->exposure_time );
//   printf( "de=%g\n", k->de );
//   printf( "frame_time=%g\n", k->frame_time );
//   printf( "integral_ae=%g\n", k->integral_ae );
//   printf( "arf_frac_exposure=%g\n", k->arf_frac_exposure );
//   printf( "\n" );
// }

#ifdef __cplusplus
extern "C" {
#endif

extern int JDMfftn (int ndim, int *dims, double * Re, double * Im,
                    int iSign, double scaling);
#ifdef __cplusplus
}
#endif

// These functions borrowed from jdmath
static unsigned int JDMbinary_search_d (double x, const double *xp,
					unsigned int n)
{
   unsigned int n0, n1, n2;
   
   n0 = 0;
   n1 = n;

   while (n1 > n0 + 1)
     {
	n2 = (n0 + n1) / 2;
	if (xp[n2] >= x) 
	  {
	     if (xp[n2] == x) return n2;
	     n1 = n2;
	  }
	else n0 = n2;
     }
   return n0;
}

// Here fft_s is a complex array with num*2*2 complex elements
static int setup_convolution_fft (double *s, unsigned int num, 
				  double *fft_s)
{
   double *right;
   unsigned int i;
   int num2;
   
   num2 = (int) num * 2;

   // Set the left part to 0 --- this serves as a pad
   memset ((char *) fft_s, 0, num2 * sizeof (double));
   right = fft_s + num2;

   for (i = 0; i < num; i++)
     {
       *right++ = s[i];	       // real
       *right++ = 0.0;	       // imag
     }

   return JDMfftn (1, &num2, fft_s, fft_s + 1, 2, -2);
}

static int do_convolution (double *fft_s, double *s, unsigned int num, 
			   double *fft_s_tmp)
{
   unsigned int num4;
   unsigned int i;
   int num2;

   (void) setup_convolution_fft (s, num, fft_s_tmp);

   // Multiply transforms
   num4 = 4 * num;
   i = 0;
   while (i < num4)
     {
	double re0, im0, re1, im1;
	unsigned int i1;
	
	i1 = i + 1;

	re0 = fft_s_tmp[i];
	re1 = fft_s[i];
	im0 = fft_s_tmp[i1];
	im1 = fft_s[i1];
	
	fft_s_tmp[i] = re0 * re1 - im0 * im1;
	fft_s_tmp[i1] = re0 * im1 + re1 * im0;

	i = i1 + 1;
     }

   num2 = (int) 2 * num;
   // Perform inverse and return left half of real part
   (void) JDMfftn (1, &num2, fft_s_tmp, fft_s_tmp + 1, -2, -2);
   for (i = 0; i < num; i++)
     s[i] = fft_s_tmp[2*i];
	
   return 0;
}

// integrate from: (flo, fhi, fy)
//               to: (tlo, thi) yielding => ty
static int rebin_histogram (double *fy, double *flo, double *fhi, int nf,
			    double *ty, double *tlo, double *thi, int nt)
{
   int f, t;

   f = 0;
   for (t = 0; t < nt; t++)
     {
        double t0, t1, s;

        t0 = tlo[t];
        t1 = thi[t];

        // Accumulate sum over 'from' bins
        // which overlap 'to' bin t

        s = 0.0;
        for ( ;f < nf; f++)
          {
             double f0, f1;
             double min_max, max_min;

             f0 = flo[f];
             f1 = fhi[f];

             if (t0 > f1)
               continue;
             if (f0 > t1)
               break;

             if (t0 > f0) max_min = t0;
             else max_min = f0;

             if (t1 < f1) min_max = t1;
             else min_max = f1;

             // prevent division by zero
             if (f0 == f1)
               return -1;

             s += fy[f] * (min_max - max_min) / (f1 - f0);

             // s += fy[f] * overlap_frac (f0, f1, t0, t1); 

             if (f1 > t1)
               break;
          }

        ty[t] = s;
     }

   return 0;
}

static int add_in_xspec (pileup_kernel_t *k, double *arf_s, double psf_frac,
			 double *results)
{
   unsigned int i;
   unsigned int num = k->num;

   psf_frac *= k->num_frames;
   for (i = 0; i < num; i++)
     results[i] += arf_s[i] * psf_frac;
   
   return 0;
}

static int perform_pileup (pileup_kernel_t *k,
			   double *arf_s,
			   double reg_size,
			   double g0,
			   double psf_frac,
			   double *gfactors,
			   double *pileup_dist,
			   double *results)
{
   unsigned int i;
   unsigned int num;
   unsigned int i_factorial;
   double fft_norm;
   double *fft_s;
   double *arf_s_fft_tmp;
   double *arf_s_tmp;
   double exp_factor;
   double integ_arf_s, integ_arf_s_n;
   double total_prob;

   num = k->num;
   fft_s = k->arf_s_fft;
   arf_s_tmp = k->arf_s_tmp;

   integ_arf_s = 0.0;
   psf_frac = psf_frac / reg_size;

   if (k->arf_frac_exposure > 0)
     psf_frac /= k->arf_frac_exposure;

   for (i = 0; i < num; i++)
     {
	arf_s_tmp[i] = results[i] = arf_s[i] * psf_frac;
	integ_arf_s += arf_s_tmp[i];
     }

   if (pileup_dist != NULL)
     {
	pileup_dist[0] = 0.0;
	pileup_dist[1] = integ_arf_s;
     }

   k->integral_ae = integ_arf_s/g0;
   k->num_terms = k->max_num_terms;

   if (integ_arf_s == 0.0)
     return 0;

   // Normalize by integ_arf_s to avoid possible floating point overflow. 
   // This will be corrected below.
   for (i = 0; i < num; i++)
     arf_s_tmp[i] /= integ_arf_s;

   exp_factor = exp (-k->integral_ae);

   (void) setup_convolution_fft (arf_s_tmp, num, fft_s);

   arf_s_fft_tmp = k->arf_s_fft_tmp;

   fft_norm = sqrt (2.0 * num);

   i_factorial = 1;
   integ_arf_s_n = integ_arf_s;
   total_prob = 1 + integ_arf_s;

   for (i = 2; i <= k->max_num_terms; i++)
     {
	unsigned int j;
	double norm_i;
	
	i_factorial *= i;
	integ_arf_s_n *= integ_arf_s;
	norm_i = integ_arf_s_n / i_factorial;
	total_prob += norm_i;

	(void) do_convolution (fft_s, arf_s_tmp, num, arf_s_fft_tmp);
	
	norm_i *= gfactors[i-2];

	for (j = 0; j < num; j++)
	  {
#ifndef HAVE_DJBFFT	     
	     arf_s_tmp[j] *= fft_norm;
#endif	     
	     results[j] += norm_i * arf_s_tmp[j];
	  }

	if (pileup_dist != NULL)
	  pileup_dist [i] = norm_i;

	if (total_prob * exp (-integ_arf_s) > Max_Probability_Cutoff)
	  {
	     k->num_terms = i;
	     break;
	  }
     }

   exp_factor *= k->num_frames * reg_size;

   // Apply correction to account for the number of effective frames
   exp_factor *= k->arf_frac_exposure;

   for (i = 0; i < num; i++)
     results [i] *= exp_factor;
   
   return 0;
}

static int convert_results (pileup_kernel_t *k, double* vals,
			    unsigned int num_bins,
			    double* energ_lo, double* energ_hi)
{
   unsigned int i;
   double *results;
   double *energies;
   unsigned int num;
   double *enlo, *enhi, *s_den;
   double de;


   num = k->num;
   results = k->results;
   energies = k->energies;

   if (NULL == (enlo = XMALLOC (3*num, double)))
     return -1;
   enhi = enlo + num;
   s_den = enhi + num;

   de = k->de;
   energies = k->energies;

   if (de >= energies[0])
     enlo[0] = (energies[0]/2.0);
   else
     enlo[0] = energies[0];
   
   for( i = 1; i < num; i++)
     enhi[num-1-i] = energies[num-1-i];
   
   for (i = 1; i < num; i++)
     enlo[i] = enhi[i-1];
   enhi[num-1] = (energies[num-1] + de);
   
   if (-1 == rebin_histogram (results, enlo, enhi, num,
			      vals, energ_lo, energ_hi, num_bins))
     {
	free ((char *)enlo);
	return -1;
     }
   free ((char *)enlo);

   return 0;
}

// By definition, S(y)|dy| = s(E)|dE|, where E = hc/y.
// So, s(E) = S(y) |dy/dE|
//          = S(y) hc/E^2
// We need the product, A(E)s(E)

static int convert_spectrum( pileup_kernel_t *k,
			     sherpa::usrfuncproto model_func,
			     sherpa::PyWrapper* x )
{
  
   double *arf_s, *energies, *arf, *enlo, *enhi, *s_den;
   unsigned int i, num;
   double de;

   if( NULL == (enlo = XMALLOC(k->num,double)))
     return -1;
   
   if( NULL == (enhi = XMALLOC(k->num,double)))
     return -1;

   if( NULL == (s_den = XMALLOC(k->num,double)))
     return -1;

   num = k->num;
   
   de = k->de;
   energies = k->energies;

   if (de >= energies[0])
    enlo[0] = (energies[0]/2.0);
   else
    enlo[0] = energies[0];
   
   for (i = 1; i < num; i++)       
     enhi[num-1-i] = energies[num-1-i];

   for (i = 1; i < num; i++)
     enlo[i] = enhi[i-1];   
   enhi[num-1] = (energies[num-1] + de);
   
   // evaluate finite grid using unconvolved source model so far
   if( EXIT_SUCCESS != model_func( &enlo[0], &enhi[0], &s_den[0], num, x ) )
     return -1;
   
   arf_s = k->arf_s;
   arf = k->arf_time;
   
   for( i = 0; i < num; i++ ) {
     // convolve source model with arf ( ARF specresp * frame_time ) 
     arf_s[i] = arf[i] * s_den[i];
     
     if (arf_s[i] < 0)
       arf_s[i] = 0;
   }

   ISIS_FREE(s_den);
   ISIS_FREE(enlo);
   ISIS_FREE(enhi);

   return 0;
}
   
static void init_coeffs (double alpha,
			 double *coeffs, unsigned int npiles)
{
   unsigned int i;

   alpha = fabs(alpha);
   for (i = 2; i <= npiles; i++)
     {
	coeffs[i - 2] = pow (alpha, (double)i-1);
     }
}

static int compute_kernel (pileup_kernel_t *k, double alpha, double g0,
			   double num_regions, double psf_frac,
			   double* results, unsigned int num_bins,
			   const double* energ_lo, const double* energ_hi,
			   sherpa::usrfuncproto model_func, sherpa::PyWrapper* x)
{
   double coeffs[MAX_NUM_TERMS+1];

   if (k == NULL)
     return -1;   

   // source model convolved with ARF in here
   if (-1 == convert_spectrum ( k, model_func, x ))
     return -1;
   
   init_coeffs (alpha, coeffs, k->max_num_terms);
      
   if (-1 == perform_pileup (k, k->arf_s, num_regions, g0, psf_frac,
			     coeffs, k->pileup_fractions,
			     k->results))
     return -1;

   // print_pileup_kernel( k );
   
   psf_frac = 1.0 - psf_frac;
   if (psf_frac > 0.0)
     {
	if (-1 == add_in_xspec (k, k->arf_s, psf_frac, k->results))
	  return -1;
     }
      
   if (-1 == convert_results (k, results, num_bins,
			      (double*) energ_lo, (double*) energ_hi ))
     return -1;
   
   return 0;
}

static void delete_kernel (pileup_kernel_t *k)
{
   if (NULL == k)
     return;

   if (k->arf_s != NULL) free (k->arf_s);
   if (k->arf_s_fft != NULL) free (k->arf_s_fft);
   if (k->arf_s_tmp != NULL) free (k->arf_s_tmp);
   if (k->arf_s_fft_tmp != NULL) free (k->arf_s_fft_tmp);
   if (k->arf_time != NULL) free(k->arf_time);
   if (k->energies != NULL) free(k->energies);
   if (k->results != NULL) free(k->results);
}


static pileup_kernel_t*
init_kernel(pileup_kernel_t *k, const double* arf_source,
	    double *pileup_fractions, double exposure_time,
	    unsigned int num_points, unsigned int max_num_terms,
	    const double *energ_lo, const double *energ_hi,
	    const double *specresp, double fracexpo, double frame_time)
{
  unsigned int i;
  double min_energy, max_energy;
  k->exposure_time = exposure_time;
  k->frame_time = frame_time;
  k->max_num_terms = max_num_terms;   
  k->num_frames = k->exposure_time / frame_time;
  k->num = NUM_POINTS;

//   k->num = num_points;
//   k->de = (max_energy - min_energy) / num_points;
  
//     It looks like John Davis hard-coded it like so:
    
//     max_energy = 15.0
//     min_energy = 0.0
//     NUM_POINTS = 1024*4
//     de = (max_energy - min_energy) / NUM_POINTS = 0.00366211 keV 
  
  k->de = (15.0 - 0.0) / NUM_POINTS;
  k->has_psf_fraction = 1;
  
  if (NULL == (k->energies = XMALLOC (NUM_POINTS, double)))
    return NULL;
  if ( NULL == (k->arf_time = XMALLOC (NUM_POINTS, double)))
    return NULL;
  if ( NULL == (k->arf_s = XMALLOC (NUM_POINTS, double)))
    return NULL;
  
  k->arf_frac_exposure = fracexpo;

//    arf_source_time = frame_time * arf_source
//    num_points is length of convolved model vals so far arf_source
//    NUM_POINTS is size of energy range 0.0 -> 15.0 kev: 4096

  k->energies[ 0 ] = 0.0001;  // MUST be non-zero
  k->arf_time[0] = 0.0;
  min_energy = energ_lo[0];
  max_energy = energ_hi[num_points-1];
  
  for(i = 1; i < k->num; i++) {
    double energ;
    
    k->energies[ i ] = energ = i * k->de;

    if( (energ >= max_energy) || ( energ < min_energy ) ) {
      k->arf_time[i] = 0.0;
      continue;
    }

    k->arf_time[i] = (specresp[JDMbinary_search_d(energ,energ_lo,num_points)]
		      * frame_time );
		      
    if (k->arf_time[i] < 0)
      k->arf_time[i] = 0;
  }

  if ( NULL == (k->results = XMALLOC (NUM_POINTS, double)))
    return NULL;
  k->pileup_fractions = pileup_fractions;
  
  if ((NULL == (k->arf_s_fft = XMALLOC (4*NUM_POINTS, double)))
      || (NULL == (k->arf_s_tmp = XMALLOC (NUM_POINTS, double)))
      || (NULL == (k->arf_s_fft_tmp = XMALLOC (4*NUM_POINTS, double))))
    {
      delete_kernel (k);
      return NULL;
    }
  
  return k;
}


int
apply_pileup(unsigned int num_bins, const double *arf_source,
	     double *results, double *pileup_fractions, double *integral_ae,
	     double exposure_time,
	     unsigned int max_num_terms, unsigned int *num_terms,
	     const double* energ_lo, const double* energ_hi,
	     const double* specresp, double fracexpo, double frame_time,
	     double alpha, double g0, double num_regions, double psf_frac,
	     sherpa::usrfuncproto model_func, sherpa::PyWrapper* x)
{
  pileup_kernel_t k;
  int rv = EXIT_FAILURE;

  if ((NULL == arf_source) || (NULL == results) ||
      (NULL == energ_lo) || (NULL == energ_hi) || (NULL == specresp) ||
      (NULL == pileup_fractions) ||
      (NULL == num_terms))
    return EXIT_FAILURE;
  
  if (NULL == init_kernel(&k, arf_source, pileup_fractions, exposure_time,
			  num_bins, max_num_terms, energ_lo, energ_hi,
			  specresp, fracexpo, frame_time))
    return EXIT_FAILURE;

  if (0 == compute_kernel(&k, alpha, g0, num_regions, psf_frac,
			  results, num_bins, energ_lo, energ_hi,
			  model_func, x)) {
    *integral_ae = k.integral_ae;
    *num_terms = k.num_terms;
    rv = EXIT_SUCCESS;
  }

  delete_kernel(&k);

  return rv;
}
