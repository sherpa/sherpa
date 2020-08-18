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
 * FILE NAME:  tcd.h
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * public header file for Transform, Convolution, Deconvolution (tcd)
 * library.
 * 
 *
 * REVISION HISTORY:
 *
 * Ref. No.        Date
 ----------       -----
 preClearCase     30March1998


 *
 H***************************************************************** */

#ifdef __cplusplus
extern "C" {
#endif

#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifndef TCD_LIB
#define TCD_LIB

#define TCDPIXELTOOFFSETNULLORIG(nAxes, lAxes, pixel, offset) \
 { \
 int ii;       \
       offset = pixel[(nAxes)-1];      \
       \
       for (ii=(nAxes-2); ii>=0; ii-- )        \
       {       \
         offset = offset * lAxes[ii] + pixel[ii];      \
       }       \
 }

/* tcdlib error codes */


enum tcdErrorCode
{
  tcdSUCCESS,               /* success                                    */
  tcdERROR_ALLOC,           /* memory allocation error                    */
  tcdERROR_UNSUPORTTYPE,    /* unsupported data type                      */
  tcdERROR_NULLPTR,         /* user passed NULL pointer                   */
  tcdERROR_NAXES0,          /* user passed nAxes <= 0                     */
  tcdERROR_PADLTOLD,        /* req. padding < original image size         */
  tcdERROR_UNKWNPAD,        /* unknown padding specified                  */
  tcdERROR_LAXES0,          /* user specified lAxes[i] = 0                */
  tcdERROR_NULLSTRING,      /* user passed empty/NULL string              */
  tcdERROR_UNSUPORTNAXES,   /* kernel lib doesn't support specified nAxes */
  tcdERROR_NAXESMISMATCH,   /* kernel and image nAxes don't match         */
  tcdERROR_INCONSISTENT,    /* length of axes specified in string not same*/
  tcdERROR_NOTIMPLEMENTED,  /* functionality not implemented yet          */
  tcdERROR                  /* generic error                              */
};



/* tcdlib constants */

#define tcdPARAM_PAD 1

#define tcdFORWARD -1
#define tcdREVERSE  1

#define tcdFREE     1
#define tcdNOFREE   0


/* -----------------------------------------
   Define complex numbers 
   -----------------------------------------
*/

typedef struct { float r; float i; } tcdComplex;
typedef struct { double r; double i;} tcdDComplex;


/* --------------------------------------
   Define enumerated types
   --------------------------------------
*/


/* Valid data types in tcdlib */

enum tcdDataType
{
  tcdBYTE,     /* 1-byte unsigned integer      */
  tcdSHORT,    /* 2-byte signed integer        */
  tcdLONG,     /* 4-byte signed integer        */
  tcdFLOAT,    /* 4-byte floating point number */
  tcdDOUBLE,   /* 8-byte floating point number */
  tcdCOMPLEX,  /* 2 tcdFLOATs                  */
  tcdDCOMPLEX  /* 2 tcdDOUBLEs                 */
};
typedef enum tcdDataType tcdDATATYPE;



/* Valid transforms in tcdlib */

enum tcdTransform
{
  tcdFFT=0,  /* discrete foutier transform (FFT)             */
  tcdDFT=0,  /* discrete fourier transform                   */
  tcdDWT,    /* discrete wavelet transform (NOT IMPLEMENTED) */
  tcdDCT,    /* discrete cosine transform  (NOT IMPLEMENTED) */
  tcdLT      /* Laplace transform          (NOT IMPLEMENTED) */
};
typedef enum tcdTransform tcdTRANSFORM;




/* Valid treatment of data edges */

enum tcdSignalEdge
{
  tcdCONSTANT,  /* use a constant at edges              */
  tcdNEAREST,   /* extrapolate nearest data point       */
  tcdWRAP,      /* wrap around data array               */
  tcdMIRROR,    /* reflect data array                   */
  tcdRENORM,    /* renormalize kernel by area over edge */
  tcdNOOP       /* no-operation (NOT USED)              */
};
typedef enum tcdSignalEdge tcdSIGNALEDGE;


/* Valid predefined kernel functions */

enum tcdFunType
{
  tcdBOX,        /* N-D constant value                */
  tcdGAUS,       /* N-D Gaussian                      */
  tcdMEXHAT,     /* N-D Mexican Hat  */
  tcdTOPHATo2D,   /* 2D only Top Hat */
  tcdBETA,
  tcdEXP,
  tcdPOWER,
  tcdSINC,
  tcdCONE,
  tcdPYRAMID
};
typedef enum tcdFunType tcdFUNTYPE;


enum tcdConORCor
{
  tcdCONVOLVE = 1,
  tcdCORRELATE = -1
};
typedef enum tcdConORCor tcdConOrCor;




enum tcdScale
{
  tcdLINEAR,
  tcdLOG,
  tcdLN,
  tcdDB
};
typedef enum tcdScale tcdSCALE;


enum tcdConMethod
{
  tcdSLIDECONVOLVE,  
  tcdFFTCONVOLVE   
};
typedef enum tcdConMethod tcdCONMETHOD;

/* -----------------------------
   Function prototypes
   -----------------------------
*/ 



/* -----------     TRANSFORM  -------------  */


/* perform N-D transform of the initialized data array */

/* double precision */
extern int tcdTransformD(
			tcdTRANSFORM  tType,  /* i: transform type          */
			double        *params, /* i: transform parameters    */
			tcdDComplex   *data,  /* i/o: initialize data array */
			long          nAxes,  /* i: number of data axes     */
			long         *lAxes,  /* i: length of data axes     */
			long         *dOrigin /* i: origin of data axes     */
			);

/* initialize data array for transform.  The routine allocates the necessary 
   memory */

extern int tcdInitTransform(
			    tcdDATATYPE dtype, /* i: input data type          */
			    void *real,        /* i: input for real part      */
			    void *imag,        /* i: input for imaginary part */
			    long  nAxes,       /* i: number of axes           */
			    long *lAxes,       /* i: length of axes           */
			    tcdComplex **data  /* o: initialized data array   */
			    );


/* double precision */
extern int tcdInitTransformD(
			    tcdDATATYPE dtype, /* i: input data type          */
			    void *real,        /* i: input for real part      */
			    void *imag,        /* i: input for imaginary part */
			    long  nAxes,       /* i: number of axes           */
			    long *lAxes,       /* i: length of axes           */
			    tcdDComplex **data  /* o: initialized data array   */
			    );


/* deallocate memory allocated by tcdInitTransform */

extern int tcdFreeTransform(
			    tcdComplex **data  /* i/o: initialized data array */
			    );

/* double precision */
extern int tcdFreeTransformD(
			    tcdDComplex **data  /* i/o: initialized data array */
			    );

extern int tcdPowerSpectrum(
		     float *psparam,
		     tcdComplex *data,
		     long nAxes,
		     long *lAxes,
		     long *dOrigin,
		     float *output
		     );


/* ----------------  Kernel Routins  -------------------*/



/* parse user supplied string for kernel array specification and origin spec */

extern int tcdBuildKernelTxt(
			     char   *kernelArray,  /* i: string specifying kernel */
			     char   *kernelOrigin, /* i: string specifying origin */
			     float **kernel,       /* o: kernel data              */
			     long   *nAxes,        /* o: number of axes           */
			     long  **lAxes,        /* o: length of kernel axes    */
			     long  **origin        /* o: origin of kernel         */
			     );


/* build a kernel from one of the built in libraries using the libParams */

extern int tcdBuildKernelLib(
		       tcdFUNTYPE   libName,  /* i: which library to call         */
		       float       *libParam, /* i: parameters to pass to library */
		       float      **kernel,   /* o: kernel data                   */
		       long        *nAxes,    /* o: number of axes                */
		       long       **lAxes,    /* o: length of data axes           */
		       long       **origin    /* o: origin of kernel              */
		       );

/* deallocate memory allocated by tcdBuildKernel routines */

extern int tcdFreeKernel(
			 float **kernel,  /* i/o: kernel data    */ 
			 long  **laxes,   /* i/o: length of axes */
			 long  **kOrigin  /* i/o: origin of axes */
			 );





/* ----------------  Convolve Routins  -------------------*/



/* perform sliding cell convolution of two arrays (data and kernel). */

/* float */
extern int tcdSlideConvolve(
		    tcdDATATYPE dtype, /* i: input data data type     */
		    void  *data,       /* i: input data array         */
		    long   nAxes,      /* i: number of axes           */
		    long  *lAxes,      /* i: length of data axes      */
		    long  *dOrigin,    /* i: origin of data axes      */
		    float *kernel,     /* i: input kernel array       */
		    long  *kAxes,      /* i: length of kernel axes    */
		    long  *kOrigin,    /* i: origin of kernel axes    */
		    float *params,     /* i: edge treatment parameters*/
		    float *output      /* o: output data array        */
		    );

/* double */
extern int tcdSlideConvolveD(
		    tcdDATATYPE dtype, /* i: input data data type     */
		    void  *data,       /* i: input data array         */
		    long   nAxes,      /* i: number of axes           */
		    long  *lAxes,      /* i: length of data axes      */
		    long  *dOrigin,    /* i: origin of data axes      */
		    double *kernel,     /* i: input kernel array       */
		    long  *kAxes,      /* i: length of kernel axes    */
		    long  *kOrigin,    /* i: origin of kernel axes    */
		    double *params,    /* i: edge treatment parameters*/
		    double *output     /* o: output data array        */
		    );



/* perform convolution by using FFT's (currently only read data types) */

extern int tcdFFTConvolve(
			  tcdConOrCor nORr,   /* i: convolve or correlate*/
			  tcdDATATYPE dtype,  /* i: input data type      */
			  void   *data,       /* i: input data array     */
			  long    nAxes,      /* i: number of axes       */
			  long   *lAxes,      /* i: length of axes       */
			  long   *dOrigin,    /* i: origin of data array */
			  tcdDATATYPE ktype,  /* i: kernel data type     */
			  void   *kernel,     /* i: kernel data          */
			  long   *kAxes,      /* i: kernel axes          */
			  long   *kOrigin,    /* i: kernel origin        */
			  float **output,     /* o: output array         */
			  long  **newAxes,     /* o: new length array     */
			  tcdComplex **fftData, /* o: fft of data array    */
			  tcdComplex **fftKern  /* o: fft of kernal array  */
			  );

/* double precision */
extern int tcdFFTConvolveD(
			  tcdConOrCor nORr,   /* i: convolve or correlate*/
			  tcdDATATYPE dtype,  /* i: input data type      */
			  void   *data,       /* i: input data array     */
			  long    nAxes,      /* i: number of axes       */
			  long   *lAxes,      /* i: length of axes       */
			  long   *dOrigin,    /* i: origin of data array */
			  tcdDATATYPE ktype,  /* i: kernel data type     */
			  void   *kernel,     /* i: kernel data          */
			  long   *kAxes,      /* i: kernel axes          */
			  long   *kOrigin,    /* i: kernel origin        */
			  double **output,     /* o: output array         */
			  long  **newAxes,     /* o: new length array     */
			  tcdDComplex **fftData, /* o: fft of data array    */
			  tcdDComplex **fftKern  /* o: fft of kernal array  */
			  );


/* perform sliding cell convolution of two arrays ( data and kernel) where
   the kernel array varies as a function of pixel location.  The kernel is
   specified as function and is recomputed at every pixel location in the
   output data array.  Memory allocations are handled by the kernel function.
   The specified kernel function MUST accept the special signal params[2] = 1
   to mean to free the memory allocated to the kernel, axes, and origin arrays */


typedef int (*KernelFn)(                /* i: user defined kernel function*/
		       long    nAxes,   /* ii: = nAxes        */
		       long   *atPixel, /* ii: = pixel-origin */
		       float  *params,  /* ii: = params       */
		       float **data,    /* oo: kernel array   */
		       long  **lAxes,   /* oo: kernla length  */
		       long  **Origin   /* oo: kernel origin  */ 
		       );


extern int tcdSlideConvolveVarKernel(
				     tcdDATATYPE dtype,  /* i: input data type  */
				     void  *data,        /* i: input data array */
				     long   nAxes,       /* i: number of axes   */ 
				     long  *lAxes,       /* i: length of axes   */
				     long  *dOrigin,     /* i: origin of axes   */
				     KernelFn kernelfn,  /* i: kernel function  */
				     float *params,      /* i: parameters       */
				     float *output       /* o: output array */
				     );


/* Allocate memory for output data array */

extern int tcdInitConvolveOutput(
				 long    nAxes, /* i: number of axes        */
				 long   *lAxes, /* i: length of axes        */
				 float **output /* o: initialize data array */
				 );

extern int tcdInitConvolveOutputD(
				 long    nAxes, /* i: number of axes        */
				 long   *lAxes, /* i: length of axes        */
				 double **output /* o: initialize data array */
				 );

/* free memory allocated by tcdInitConvolveOutput */

extern int tcdFreeConvolveOutput(
				 float **output /* i/o: initialized data array */
				 );

extern int tcdFreeConvolveOutputD(
				 double **output /* i/o: initialized data array */
				 );




/*----------------MISC HELPER ROUTINES----------------------*/


/* convert a pixel location to an offset in the data array */

extern int tcdPixelToOffset(
			    long  nAxes,  /* i: number of data axes  */
			    long *lAxes,  /* i: length of data axes  */
			    long *origin, /* i: iorigin of data axes */
			    long *pixel,  /* i: pixel to convert     */
			    long *offset  /* o: offset into array    */
			    );

/* convert an offset in the data array into a pixel location */

extern int tcdOffsetToPixel(
			    long nAxes,   /* i: number of data axes */
			    long *lAxes,  /* i: length of data axes */
			    long *origin, /* i: origin of data */
			    long offset,  /* i: offset into data array */
			    long *pixel   /* o: returned pixel location */
			    );


/* Routines to pad data.  Padding is always done to "right" side of data */
/* This routine will allocate memory for the data array.  The returned
   data array has the same tcdDATATYPE as the input data array*/


/* Increase axes (lAxes) by addAxes amount, ie nlAxes[i]=lAxes[i]+addAxes[i]*/

extern int tcdPadDataWith(
			  tcdDATATYPE dtype, /* i: input data type        */
			  void  *data,       /* i: data array             */
			  long   nAxes,      /* i: num data axes          */
			  long  *lAxes,      /* i: length data axes       */
			  long  *addAxes,    /* i: increase data axes by  */
			  void  *output,     /* o: padded data array      */
			  long **nlAxes      /* o: length of padded array */
			  );

/* Pad data array to the amount specified in newAxes array.  newAxes must
   be >= lAxes */

extern int tcdPadDataSpec(
			  tcdDATATYPE dtype, /* i: input data type     */
			  void  *data,       /* i: data array          */
			  long   nAxes,      /* i: number data axes    */
			  long  *lAxes,      /* i: lengths of axes     */
			  long  *newAxes,    /* i: new lengths of axes */
			  void  *output      /* o: padded data array   */
			  );


/* Pad data array to next 2^N size.  The order parameter can be used to
   increase the size even more, ie 2^(N+order), ie 2^(N+1).  Typically order 
   will be == 0 */

extern int tcdPadData2N(
			tcdDATATYPE dtype, /* i: input data type        */
			void  *data,       /* i: data array             */
			long   nAxes,      /* i: number data axes       */
			long  *lAxes,      /* i: lengths of axes        */
			long   order,      /* i: which 2^(N+order)      */
			void  *output,     /* o: padded data array      */
			long **nlAxes      /* o: length of padded array */
			);



/* free memory allocated by tcdPadData routines */

extern int tcdFreePadData(
			  void *data,  /* i: padded data array            */
			  long **lAxes /* i: lengths of padded data array */
			  );




/* crop data array specifying corner to crop at*/

extern int tcdCropData(
                tcdDATATYPE  dType,   /* i: input data type           */
                void        *data,    /* i: input data array          */
                long         nAxes,   /* i: number of axes            */
                long        *lAxes,   /* i: length of axes            */
                long        *dOrigin, /* i: origin of axes            */
                long        *loLeft,  /* i: lower left corner r.t.o.  */
                long        *upRite,  /* i: upper right corner r.t.o. */
                void        *output,  /* o: output data array         */
                long       **oAxes,   /* o: output data array lengths */
                long       **oOrigin  /* o: output data array origin  */
                );


/* cast an array of 1 data type to an array of another data type. */
/* This routine allocates memory for the output data array        */
/* There is a special version to cast data to type tcdCOMPLEX     */

extern int tcdCastArray( 
			tcdDATATYPE inType,  /* i: input data type     */
			void *inData,        /* i: input data array    */
			long  nAxes,         /* i: number of data axes */
			long *lAxes,         /* i: length of data axes */
			tcdDATATYPE outType, /* i: output data type    */
			void *outData        /* o: output data array   */
			);

/* cast a 2 data arrays to the real & imaginary parts of a complex data array */
/* This routine allocates the memory for the output data array */

extern int tcdCastToComplex(
			    tcdDATATYPE dtype,/*i: input data types             */
			    void     *real,   /*i: pointer to real part of array */
			    void     *imag,   /*i: pointer to imag part of array */
			    long      nAxes,  /*i: number of axes                */
			    long     *lAxes,  /*i: length of axes                */
			    tcdComplex *data  /*o: complex data returned         */
			    );

extern int tcdCastToDComplex(
			    tcdDATATYPE dtype,/*i: input data types             */
			    void     *real,   /*i: pointer to real part of array */
			    void     *imag,   /*i: pointer to imag part of array */
			    long      nAxes,  /*i: number of axes                */
			    long     *lAxes,  /*i: length of axes                */
			    tcdDComplex *data  /*o: complex data returned         */
			    );


/* free memory allocated by tcdCast... routines */
/* NOT IMPLEMENTED */
extern int tcdFreeCastArray( 
			    void *data /* i/o: cast data array */
			    );



/* convert a complex data array, which can be interpreted as an array  of
   floats with real, imag, real, imag, real, imag, ..... into two arrays
   of floats with real, real, ... real, imag, imag, imag..... */
/* This is for easy output of the data array into two seperate peices */
/* This conversion is done IN PLACE which means that real and imag are
   merely pointers inside the data array.  The data array itself is mangled
   by this operation.*/
/* Do not free the real or imag arrays.  Only free the data array */

extern int tcdComplexReOrder(
			     tcdComplex *data,  /* i/o: input data array     */
			     long        nAxes, /* i: number of data axes    */
			     long       *lAxes, /* i: length of axes         */
			     float     **real,  /* o: real data array        */
			     float     **imag   /* o: imag data array        */
			     );





extern int tcdShiftArray(
			 tcdDATATYPE dType, /* i: input data type         */
			 void *data,        /* i/o: data array            */
			 long  nAxes,       /* i: number of axes          */
			 long *lAxes,       /* i: length of axes          */
			 long *dOrigin,     /* i: origin of axes          */
			 long *shiftBy      /* i: amount to shift axes by */
			 );

extern int tcdFlipArray(
			 tcdDATATYPE dType, /* i: input data type         */
			 void *data,        /* i/o: data array            */
			 long  nAxes,       /* i: number of axes          */
			 long *lAxes,       /* i: length of axes          */
			 long *dOrigin      /* i/0: origin of axes        */
			 );

/* -----------     CSMOOTH  -------------  */

/* Function prototypes */

extern int        tcdAdaptiveSmooth(tcdDATATYPE, tcdDATATYPE, tcdDATATYPE, 
                                    void*, float, float, float, 
                                    float, float, long, long*, long*,
                                    tcdFUNTYPE, int, int, int, int, float*,
                                    int, int, float, float,
                                    void*, void*, float*, float*, float*);

extern float *tcdLucy( tcdDATATYPE imtype,
             void *data,
             long nAxes,
             long *liAxes,
             long *oiOrigin,
             tcdDATATYPE ktype,
             void *kernel,
             long *lkAxes,
             long *okOrigin,
             long niter
             );

#ifdef __cplusplus
}
#endif


#endif
