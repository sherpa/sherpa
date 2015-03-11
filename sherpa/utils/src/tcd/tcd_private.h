/*                                                                
**  Copyright (C) 1998-2007  Smithsonian Astrophysical Observatory 
*/                                                                

/*                                                                          */
/*  This program is free software; you can redistribute it and/or modify    */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 3 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful,         */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License along */
/*  with this program; if not, write to the Free Software Foundation, Inc., */
/*  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.             */
/*                                                                          */


/*H*****************************************************************
 * FILE NAME:  tcd_private.h
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * private header file for Transform, Convolution, Deconvolution (tcd)
 * library.  The function prototypes defined here should not be made
 * available to the general user.
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.        Date
 ----------       -----
 preClearCase     30March1998


 *
 H***************************************************************** */

#ifndef TCD_PRIVATE_H
#define TCD_PRIVATE_H

/* ------------------------------------------------------ */
/* Priviate functions */
/* ------------------------------------------------------ */



extern int tcdCheckAxes( long nAxes, long *lAxes );
extern int tcdCheckData( void *data, long nAxes, long *lAxes );



extern int tcdFFT1D( 
		    long        n, 
		    tcdComplex *data, 
		    long        step, 
		    short       dir 
		    );


typedef int (*Xform)(long, tcdComplex *, long, short );

extern int tcdGenXformDrive( 
			     Xform xform,
			     tcdComplex *data,
			     long nAxes,
			     long *lAxes,
			     short dir
			     );


extern int tcdKernelArray(
			  char *kernelArray,
			  char *kernelOrigin,
			  float **kernel,
			  long  *nAxes,
			  long  **lAxes,
			  long  **origin
			  );

extern int tcdDetLength(
			char *spec,   /* i: array string specification      */
			long  *lAxes, /* o: array of lengths                */
			long  level,  /* i: current dimension               */
			long  *off    /* o: offset in string thru recursion */
			);


extern int tcd_strNDarray(
			  tcdDATATYPE dtype,
			  char   *spec,  /* input array specification   */
			  void  *data,  /* ouptut data array           */
			  long  *nAxes,  /* output number of dimensions */
			  long **lAxes   /* length of data axes         */
			    );


extern int tcdKernelLib(
			tcdFUNTYPE *funType,
			float  *libParam,
			float **kernel,
			long   *nAxes,
			long  **lAxes,
			long  **origin
			);

extern int tcdGetKernelAxes(
                                tcdFUNTYPE libName,
                                float  *libParam,
                                float scale,
                                long  **lAxes,
                                long  **origin
                                );

extern int tcdGetLbgKernelAxes(
                                tcdFUNTYPE libName,
                                float  *libParam,
                                float scale,
                                long  **lAxes,
                                long  **origin
                                );

extern int tcdKernelLib_box(
			    float  *params,
			    float **kernel,
			    long   *nAxes,
			    long  **lAxes,
			    long  **origin
			    );


extern float tcdSignalEdge(
			   tcdDATATYPE dtype,
			   void *data,
			   long   nAxes,
			   long  *lAxes,
			   long  *atPixel,
			   float *params,
			   int *status
			   );

/* Valid ways to pad data */

enum tcdPadType
{
  tcdPAD2KERNEL,  /* Increase lengths by kernel size */
  tcdPAD2N,       /* Increase to power of 2^N        */
  tcdPAD2SPEC     /* Increase to specified size      */
};
typedef enum tcdPadType tcdPADTYPE;


/* >>>>>> This routine will be broken up into several other routines  <<<<<<<< */

extern int tcdPadData(
		      tcdPADTYPE   ptype,  /* i: how to pad data       */
		      long        *param,  /* i: parameters needed     */
 		      tcdDATATYPE  dtype,  /* i: input data type       */
		      void        *data,   /* i: pointer to data array */
		      long         nAxes,  /* i: number of data axes   */
		      long        *lAxes,  /* i: length of data axes   */
		      void       *output,  /* o: output data array     */
		      long       **nlAxes  /* o: length of output data array */
		      );



/* Valid kernel types functions */

enum tcdKernel
{
  tcdDATA,   /* input is a preallocated array (may need to pad) */
  tcdARRAY,  /* input is a text string defining array           */
  tcdLIB,    /* input is a call to built in library             */
  tcdFUN = tcdLIB,
  tcdAUXLIB  /* NOT IMPLEMENTED                                 */
};
typedef enum tcdKernel tcdKERNEL;


/* create or cast kernel data array (WILL BE WRAPPED BY SOME FUNCTIONS) */

extern int tcdBuildKernel(
			  tcdKERNEL  kType,   /* I: input kernel type            */
			  void      *kParam1, /* i: kernel specific params 1     */
			  void      *kParam2, /* i: kernel specific params 2     */
			  long       nAxes,   /* i: number dimension on input arr*/
			  float    **kernel,  /* o: output kernel                */
			  long     **lAxes,   /* o: output dim. of kernel        */
			  long     **kAxes    /* o: indices of origin in kernel  */
			  );



typedef float (*KernelData)(
			    float *p,
			    long   n, 
			    float *pos, 
			    int   *status
			    );

typedef int (*KernelSize)( 
			  float *p, 
			  long   n, 
			  long  *ax,  
			  long  *tor
			  );


extern int tcdKernelLibStep(
			    KernelData kernelData,
			    float *params,
			    float *kernel,
			    long   nAxes,
			    long  *lAxes,
			    long  *origin,
			    long   stepSize,
			    long  *atPixel,
			    float *atPos,
			    long   nAt,
			    int   *status
			    );


extern float tcdKernelLibEval(
			      KernelData kernelData,
			      float *params,
			      long   nAxes,
			      long   stepSize,
			      long  *atPixel,
			      float *atPos,
			      long   nAt,
			      int   *status
			      );


extern int tcdKernel_BOX_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_BOX_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_CONE_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_CONE_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_PYRAMID_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_PYRAMID_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_TOPHATo2D_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_TOPHATo2D_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_GAUS_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_GAUS_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_MEXHAT_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_MEXHAT_data( float *p, long n, float *pos, int *status);


extern int tcdKernel_BETA_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_BETA_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_POWER_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_POWER_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_EXP_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_EXP_data( float *p, long n, float *pos, int *status);

extern int tcdKernel_SINC_size( float *p, long n, long *ax, long *tor);

extern float tcdKernel_SINC_data( float *p, long n, float *pos, int *status);



extern float tcdKernelLoop_f(
                      tcdDATATYPE dtype,  /* i: input data type      */
                      void *data,         /* i: pointer to data      */
                      long   nAxes,       /* i: number of axes       */
                      long  *lAxes,       /* i: length of axes       */
                      long  *dOrigin,     /* i: origin of data       */
                      float *kernel,      /* i: pointer to kernel    */
                      long  *kAxes,       /* i: length of axes       */
                      long  *kOrigin,     /* i: origin of kernel     */
                      long  *atPixel,     /* i: current pixel location */
                      long  *atOffset,    /* i: current kernel offset  */
                      long  *signalPixel, /* i: current pixle location */
                      long   nAt,         /* i: current level of recursion */
                      float *params,      /* i: edge treatment parameters  */
                      float *kRenEdge,    /* o: integ of kernel over edge */
                      int   *status       /* o: error status */
                      );

extern int        tcdBuildAdaptiveKernel(long, tcdFUNTYPE, float*,
                                         float, float**,
                                         long**, long**, float*);

extern int        tcdBuildAdaptiveLbgKernel(long, tcdFUNTYPE, float*,
                                            float, float**,
                                            long**, long**, float*);

extern int        tcdAdaptiveSig(float*, long, long*, float*, float*, float*, 
                                 long, float*);

extern int        tcdNumPtsMatrix(long, long*, long*);

extern int        tcdCreateLongMatrix(long**, long, long*);

extern int        tcdSetLongMatrix(long*, long,long*, long);

extern int        tcdSetMatrix(float*, long,long*, float);

extern int        tcdSumMatrix(float*, long, float*);

extern int        tcdMultiplyMatrix(float*, float*, long, long*, float*);

extern int        tcdScaleMatrix(float*, float, long, long*, float*);

extern int        tcdSqrtMatrix(float*, long, float*);

extern int        tcdDivideMatrix(float*, float*, long, long*, float*);

extern int        tcdSubtractMatrix(float*, float*, long, long*, float*);

extern int        tcdSearchMask(float*, long, long*, float, float*,
                                long*, int*);

extern int        tcdSearchMatrix(float*, long, long*, float, long*,
                                  int*);

extern int        tcdAdaptivePrint(float*, float*, float*,
                                   float*, float*, long*,
                                   long, long*, long,
                                   float, float*);

extern int        tcdAdaptiveMaxMedMinMatrix(float*, long, float*, float*,
                                             float*);

extern int        tcd_sortLo2Hi_floats(const void*, const void *);

#endif
