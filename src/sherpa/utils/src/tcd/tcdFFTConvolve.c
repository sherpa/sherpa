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


#include "tcd.h"
#include "tcd_private.h"

int tcdPhaseShift( void *data, long nAxes, long *lAxes,
		   long *dOrigin, long *kAxes, long size );


/* double precision */
int tcdFFTConvolveD(
		   tcdConOrCor nORr,     /* i: convolve or correlate*/
		   tcdDATATYPE dtype,    /* i: input data type      */
		   void  *data,          /* i: input data array     */
		   long   nAxes,         /* i: number of axes       */
		   long  *lAxes,         /* i: length of axes       */
		   long  *dOrigin,       /* i: origin of data array */
		   tcdDATATYPE ktype,    /* i: kernel data type     */
		   void  *kernel,        /* i: kernel data          */
		   long  *kAxes,         /* i: kernel axes          */
		   long  *kOrigin,       /* i: kernel origin        */
		   double **output,       /* o: output array         */
		   long  **newAxes,      /* o: new length array     */
		   tcdDComplex **fftData, /* o: fft of data array    */
		   tcdDComplex **fftKern  /* o: fft of kernal array  */
		   )
{

  tcdDComplex *product  = NULL;

  void *data_p = NULL;
  void *kern_p = NULL;

  long ii;
  long nTotal = 1;

  int  padData = 0;
  int  padKern = 0;

  int status;
  double dxformParam[2] = { tcdFORWARD, 0 };



  /* check input data */
  if ( data != NULL )
    {
      status = tcdCheckAxes( nAxes, lAxes );
      if ( status != tcdSUCCESS ) return( status );
    }

  if ( kernel != NULL )
    {
      status = tcdCheckAxes( nAxes, kAxes );
      if ( status != tcdSUCCESS ) return( status );
    }

  /* check NULL pointers */

  if (( data == NULL ) && ( *fftData == NULL )) return(tcdERROR_NULLPTR);

  if (( kernel == NULL) && ( *fftKern == NULL )) return(tcdERROR_NULLPTR);

  if ( ( (data == NULL ) || (kernel == NULL) )  && 
       ( (*newAxes) == NULL ) ) return( tcdERROR_NULLPTR);

  /* determine if either data or kernel needs padded (pad to largest) */

  if ( ( data != NULL ) && ( kernel != NULL ) )
    {
      
      (*newAxes) = ( long *)calloc( nAxes, sizeof( long ));
      if ((*newAxes) == NULL ) return( tcdERROR_ALLOC );
      
      for ( ii=0;ii<nAxes;ii++) 
	{ 
	  (*newAxes)[ii] = ( lAxes[ii] > kAxes[ii] ) ? lAxes[ii] : kAxes[ii] ;
	}

    }


  nTotal = 1;
  for (ii=0;ii<nAxes;ii++)
    {
      if ( data != NULL ) 
	{
	  if ( (*newAxes)[ii] > lAxes[ii] ) padData = 1;
	}

      if ( kernel != NULL )
	{
	  if ( (*newAxes)[ii] > kAxes[ii] ) padKern = 1;
	}

      nTotal *= (*newAxes)[ii];

    }



  /* data routines */
  /* pad data if needed */

  if ( data != NULL )
    {

      if ( padData == 1)
	{
	  status = tcdPadDataSpec( dtype, data, nAxes, lAxes, (*newAxes), 
				   &data_p);
	  if ( status != tcdSUCCESS ) return( status );
	}
      else
	{
	  data_p = data;
	}

      /* copy to complex array */
      
      status = tcdInitTransformD( dtype, data_p, NULL, nAxes, (*newAxes),
				 fftData );
      if ( status != tcdSUCCESS ) return( status );
      if ( padData == 1 ) free( data_p );
      
      /* compute fft */
      
      status = tcdTransformD( tcdFFT, dxformParam, *fftData, nAxes, (*newAxes),
			     dOrigin );

    } /* end if data array */


  /* kernel routines */
  /* Pad kernel if needed */

  
  if ( kernel != NULL )
    {


      if ( padKern == 1)
	{
	  status = tcdPadDataSpec( ktype, kernel, nAxes, kAxes, (*newAxes),
				   &kern_p);
	  if ( status != tcdSUCCESS ) return( status );
	}
      else
	{
	  kern_p = kernel;
	}

      /* copy to complex array */
      
      status = tcdInitTransformD( ktype, kern_p, NULL, nAxes, (*newAxes),
				 fftKern );
      if ( status != tcdSUCCESS ) return( status );
      if ( padKern == 1 ) free( kern_p );
      
      /* compute fft */

      if (nAxes > 3)
	return tcdERROR_NOTIMPLEMENTED;

      if ((dOrigin[0] != kOrigin[0]) || 
	  ((nAxes == 2) && (dOrigin[1] != kOrigin[1])) ||
	  ((nAxes == 3) && (dOrigin[2] != kOrigin[2])))
	tcdPhaseShift(*fftKern, nAxes, *newAxes, kOrigin, kAxes,
		      sizeof(tcdDComplex));

      status = tcdTransformD( tcdFFT, dxformParam, *fftKern, nAxes, (*newAxes), 
			     kOrigin );


      if ( status != tcdSUCCESS ) return( status );

    } /* end if kernal specified */


  /* alloc memory for convol array */

  product = ( tcdDComplex *) calloc( nTotal, sizeof( tcdDComplex));
  if ( product == NULL ) return ( tcdERROR_ALLOC );

  /* multiply arrays */

  for (ii=0;ii< nTotal; ii++ )
    {
      product[ii].r = (*fftData)[ii].r * (*fftKern)[ii].r - nORr *
	              (*fftData)[ii].i * (*fftKern)[ii].i;
      product[ii].i = (*fftData)[ii].i * (*fftKern)[ii].r + nORr *
	              (*fftData)[ii].r * (*fftKern)[ii].i;
    }


  /* inverse fft */
  dxformParam[0] = tcdREVERSE ;
  
  if ( nORr == tcdCORRELATE )
    {
      status = tcdTransformD( tcdFFT, dxformParam, product, nAxes, (*newAxes),
			     NULL );
    }
  else
    {
      status = tcdTransformD( tcdFFT, dxformParam, product, nAxes, (*newAxes),
			     dOrigin );
    }

  if ( status != tcdSUCCESS) return ( status );

  /* cast array to real */
  status = tcdCastArray( tcdDCOMPLEX, product, nAxes, (*newAxes), tcdDOUBLE,
			 output );
  if ( status != tcdSUCCESS ) return( status );

  /* need to normalize */

  for (ii=0; ii< nTotal; ii++) (*output)[ii] *= nTotal;
  

  /* save data */
  
  free(product);

  return( tcdSUCCESS );

}

static int phase_shift_1d(void*, long, long*, long*, long*, long);
static int phase_shift_2d(void*, long, long*, long*, long*, long);
static int phase_shift_3d(void*, long, long*, long*, long*, long);

/*
  +----------------------------------------------
  +
  + Apply a phase shift to the data based on the direction of the
  + transform.
  +
  +----------------------------------------------
  */
int tcdPhaseShift( void *data, 
		   long nAxes, 
		   long *lAxes,
		   long *kOrigin, 
		   long *kAxes,
		   long size
		   )
{

  int  status, ii;

  long * kOriginSafe;

  kOriginSafe = (long *) malloc(sizeof(long)*nAxes);

  /* deal with origin values that may be out-of-range */
  for (ii=0; ii<nAxes; ii++)
    {
      kOriginSafe[ii] = kOrigin[ii] % kAxes[ii];
      if (kOriginSafe[ii] < 0)
	kOriginSafe[ii] += kAxes[ii];
    }

  switch (nAxes) {
  case 1: 
    status = phase_shift_1d(data, nAxes, lAxes, kOriginSafe, kAxes, size);
    break;
  case 2: 
    status = phase_shift_2d(data, nAxes, lAxes, kOriginSafe, kAxes, size);
    break;
  case 3: 
    status = phase_shift_3d(data, nAxes, lAxes, kOriginSafe, kAxes, size);
    break;
  default:
    status = tcdERROR;
  }

  free (kOriginSafe);

  return( status );

}

/* 2 dimensional phase shift - usual case */
static int phase_shift_2d ( void *data, 
			    long  nAxes, 
			    long *lAxes,
			    long *kOrigin, 
			    long *kAxes,
			    long size
		   )
{

  long jj;  /* loop varialbes */

  long nTotal;  /* number of data points */

  long	x_fore, x_aft;
  long	y_fore, y_aft;
  long	s_offset, d_offset;
  long	dummy[2]={0,0};
  long	pixel[2];
  char *shadow;

  nTotal = lAxes[0]*lAxes[1];

  shadow = malloc(nTotal*size);
  memset(shadow, 0, nTotal*size);

  x_fore = kOrigin[0];
  x_aft = kAxes[0] - x_fore;
  y_fore = kOrigin[1];
  y_aft = kAxes[1] - y_fore;

  /* 
     
  Padded Kernel
  
  . . . . . . . .
  . . . . . . . .
  . . . . . . . .
  . . . . . . . .
  . . . . . . . .
  . . . . . . . .
  K K . . . . . .
  K K . . . . . .
  
  
  
  Kernel Quadrants
  
  2 1
  3 4
  
  origin is at lower left of quandrant 1
  
  y_fore is height of quadrants 3, 4 (not including origin)
  y_aft is height of quadrants 1 and 2 (including origin)
  
  x_fore is width of quadrants 2, 3 (not including origin)
  x_aft is width of quadrants 1 and 4 (including origin)
  */
  
  /* start with 1 */
    for (jj=0; jj<y_aft; jj++)
      { 
	pixel[0] = kOrigin[0];
	pixel[1] = jj+kOrigin[1];
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
	pixel[0] = 0;
	pixel[1] = jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
	memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]), 
	       x_aft*size);
      }
  
  /* now 2 */
    for (jj=0; jj<y_aft; jj++)
      {
	pixel[0] = 0;
	pixel[1] = jj + kOrigin[1];
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
	pixel[0] = lAxes[0]-x_fore;
	pixel[1] = jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
	memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]), 
	       x_fore*size);
      }
  
  
  /* now 3 */
    for (jj=0; jj<y_fore; jj++)
      {
	pixel[0] = 0;
	pixel[1] = jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
	pixel[0] = lAxes[0]-x_fore;
	pixel[1] = lAxes[1]-y_fore+jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
	
	memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]), 
	       x_fore*size);
      }
  
  
  /* now 4 */
    for (jj=0; jj<y_fore; jj++)
      {
	pixel[0] = kOrigin[0];
	pixel[1] = jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
	pixel[0] = 0;
	pixel[1] = lAxes[1]-y_fore+jj;
	tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
	memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]), 
	       x_aft*size);
      }

  memcpy(data, shadow, nTotal*size);

  free(shadow);

  return tcdSUCCESS;
}
		  
/* 1-D phase shift */
static int phase_shift_1d ( void *data, 
			    long nAxes, 
			    long *lAxes,
			    long *kOrigin, 
			    long *kAxes,
			    long size
		   )
{

  long nTotal;  /* number of data points */

  long	x_fore, x_aft;
  long	s_offset, d_offset;
  long	dummy[1]={0};
  long	pixel[1];
  char *shadow;

  nTotal = lAxes[0];

  shadow = malloc(nTotal*size);
  memset(shadow, 0, nTotal*size);

  x_fore = kOrigin[0];
  x_aft = kAxes[0] - x_fore;

  /* 
     
  Padded Kernel
  
  K K . . . . . .
  
  
  
  Kernel Halves
  
  2 1

  
  origin is at left of half 1
  
  
  x_fore is width of half 2 (not including origin)
  x_aft is width of half 1 (including origin)
  */
  
  /* start with 1 */

    {
      pixel[0] = kOrigin[0];
      tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
      pixel[0] = 0;
      tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
      memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
	     x_aft*size);
    }
  
  /* now 2 */
    {
      pixel[0] = 0;
      tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &s_offset);
      pixel[0] = lAxes[0]-x_fore;
      tcdPixelToOffset(nAxes, lAxes, (long *) &dummy, (long *) &pixel, &d_offset);
      memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
	     x_fore*size);
    }

  memcpy(data, shadow, nTotal*size);

  free(shadow);

  return tcdSUCCESS;
}
		  
      

/* 3D pahse shift - the monster! 
   study 2D shift before takling this one
*/
static int phase_shift_3d (void *data, 
			   long  nAxes, 
			   long *lAxes,
			   long *kOrigin, 
			   long *kAxes,
			   long size
		   )
{

  long jj, kk;  /* loop varialbes */

  long nTotal;  /* number of data points */

  long	x_fore, x_aft;
  long	y_fore, y_aft;
  long	z_fore, z_aft;
  long	s_offset, d_offset;
  long	dummy[3]={0,0,0};
  long	pixel[3];
  char *shadow;

  nTotal = lAxes[0]*lAxes[1]*lAxes[2];

  shadow = malloc(nTotal*size);
  memset(shadow, 0, nTotal*size);

  x_fore = kOrigin[0];
  x_aft = kAxes[0] - x_fore;
  y_fore = kOrigin[1];
  y_aft = kAxes[1] - y_fore;
  z_fore = kOrigin[2];
  z_aft = kAxes[2] - z_fore;

  /* 
     
  Padded Kernel
  
      top                bottom
  . . . . . . . .    . . . . . . . .
  . . . . . . . .    . . . . . . . .
  . . . . . . . .    . . . . . . . .
  . . . . . . . .    . . . . . . . .
  . . . . . . . .    . . . . . . . .
  . . . . . . . .    . . . . . . . .
  K K . . . . . .    K K . . . . . .
  K K . . . . . .    K K . . . . . .
                     ^
                  (0,0,0)
  
  
  Kernel Octants
  
  2 1                6 5
  3 4                7 8
  
  origin is at lower left corner of octant 1
  
  x_fore is width of octants 2, 3, 6, 7 (not including origin)
  x_aft is width of octants 1, 4, 5, 8 (including origin)

  y_fore is height of octants 3, 4, 7, 8 (not including origin)
  y_aft is height of octants 1, 2, 5, 6 (including origin)
  
  z_fore is depth of octants 5, 6, 7, 8 (not including origin)
  z_aft is depth of octants 1, 2, 3, 4 (including origin)
  
  */
  
  /* Do top half of kernel to bottom half of cube */

  /* 3D start with octant 1 */
    for (jj=0; jj<y_aft; jj++)
      for (kk=0; kk<z_aft; kk++)
	{
	  pixel[0] = kOrigin[0];
	  pixel[1] = jj+kOrigin[1];
	  pixel[2] = kk+kOrigin[2];
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = 0;
	  pixel[1] = jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_aft*size);
	}
  
  /* 3D octant 2 */
    for (jj=0; jj<y_aft; jj++)
      for (kk=0; kk<z_aft; kk++)
	{
	  pixel[0] = 0;
	  pixel[1] = jj+kOrigin[1];
	  pixel[2] = kk+kOrigin[2];
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = lAxes[0]-x_fore;
	  pixel[1] = jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_fore*size);
      }
  
  
  /* 3D octant 3 */
    for (jj=0; jj<y_fore; jj++)
      for (kk=0; kk<z_aft; kk++)
	{
	  pixel[0] = 0;
	  pixel[1] = jj;
	  pixel[2] = kk+kOrigin[2];
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = lAxes[0]-x_fore;
	  pixel[1] = lAxes[1]-y_fore+jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_fore*size);
      }
  
  
  /* 3D octant 4 */
    for (jj=0; jj<y_fore; jj++)
      for (kk=0; kk<z_aft; kk++)
	{
	  pixel[0] = kOrigin[0];
	  pixel[1] = jj;
	  pixel[2] = kk+kOrigin[2];
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = 0;
	  pixel[1] = lAxes[1]-y_fore+jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_aft*size);
	}
  
  /* Do bottom half of kernel to top half of cube */

  /* 3D start with octant 5 */
    for (jj=0; jj<y_aft; jj++)
      for (kk=0; kk<z_fore; kk++)
	{
	  pixel[0] = kOrigin[0];
	  pixel[1] = jj+kOrigin[1];
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = 0;
	  pixel[1] = jj;
	  pixel[2] = lAxes[2]-z_fore+kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_aft*size);
	}
  
  /* 3D octant 6 */
    for (jj=0; jj<y_aft; jj++)
      for (kk=0; kk<z_fore; kk++)
	{
	  pixel[0] = 0;
	  pixel[1] = jj+kOrigin[1];
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = lAxes[0]-x_fore;
	  pixel[1] = jj;
	  pixel[2] = lAxes[2]-z_fore+kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_fore*size);
      }
  
  
  /* 3D octant 7 */
    for (jj=0; jj<y_fore; jj++)
      for (kk=0; kk<z_fore; kk++)
	{
	  pixel[0] = 0;
	  pixel[1] = jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = lAxes[0]-x_fore;
	  pixel[1] = lAxes[1]-y_fore+jj;
	  pixel[2] = lAxes[2]-z_fore+kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_fore*size);
      }
  
  
  /* 3D octant 8 */
    for (jj=0; jj<y_fore; jj++)
      for (kk=0; kk<z_fore; kk++)
	{
	  pixel[0] = kOrigin[0];
	  pixel[1] = jj;
	  pixel[2] = kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &s_offset);
	  pixel[0] = 0;
	  pixel[1] = lAxes[1]-y_fore+jj;
	  pixel[2] = lAxes[2]-z_fore+kk;
	  tcdPixelToOffset(nAxes, lAxes, (long*) &dummy, (long *) &pixel, &d_offset);
	  memcpy(&(shadow[d_offset*size]), &(((char *)data)[s_offset*size]),
		 x_aft*size);
	}

  memcpy(data, shadow, nTotal*size);

  free(shadow);

  return tcdSUCCESS;
}


