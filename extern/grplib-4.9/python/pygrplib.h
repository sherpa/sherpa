/*                                                                
**  Copyright (C) 2007  Smithsonian Astrophysical Observatory 
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


/* H*****************************************************************
 *
 * FILE NAME:  pygrplib.h
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:
 *
 * header info for pygrplib tool ( prototypes and macros )
 *
 *
 * REVISION HISTORY:
 *
 * Ref. No.         Date
   ----------       -----
   0.1              April2007 	File Created
H***************************************************************** */
/* Tool has been cleaned to Numpy API version 1.7 */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Python.h"
#include "numpy/arrayobject.h"  /* Used by NumPy */

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Macro definitions for pygrplib.c
 * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Macro to cast NumPy array data to a 1-d double array.
 */
#define IDATA(p) ((int *) ( PyArray_DATA(((PyArrayObject *)p)) ))
#define SDATA(p) ((short *) ( PyArray_DATA(((PyArrayObject *)p)) ))
#define DDATA(p) ((double *) ( PyArray_DATA(((PyArrayObject *)p)) ))


/* * * * * * * * * * * * * * * * * * * * * * * *
 * Function Declarations for pygrplib.c
 * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * This function returns the grouping and quality arrays that represent the
 * input data (countsArray) after it has been adaptively grouped so that
 * each group contains at least numCounts counts.
 */
static PyObject *grpAdaptive(PyObject *self,	/*i: Used by Python */
                             PyObject *args,	/*i: Python tuple of the arguments */
                             PyObject *kwds  /*i: Python tuple of keywords */
                             );

/*
 * This function returns the grouping and quality arrays that represent the
 * input data (countsArray) after it has been adaptively grouped so that the
 * signal to noise of each group is at least equal to the snr parameter. The
 * errorCol array gives the error for each element of the original array: if
 * it is not supplied then the error is taken to be the square root of the
 * element value.
 */
static PyObject *grpAdaptiveSnr(PyObject *self,	/*i: Used by Python */
                                PyObject *args,	/*i: Python tuple of the arguments */
                                PyObject *kwds  /*i: Python tuple of keywords */
                                );

/*
 * This function returns the grouping and quality arrays for a set of groups
 * defined by the low (binLowArray) and high (binHighArray) boundaries when
 * applied to the axis values of the data (axisArray).
 */
static PyObject *grpBin(PyObject *self,	/*i: Used by Python */
                        PyObject *args,	/*i: Python tuple of the arguments */
                        PyObject *kwds  /*i: Python tuple of keywords */
                        );

/*
 * This function allows you to calculate the grouping information needed to
 * group the input data (the axisArray array) to match the grouping of another
 * dataset
 */
static PyObject *grpBinFile(PyObject *self,	/*i: Used by Python */
                            PyObject *args,	/*i: Python tuple of the arguments */
                            PyObject *kwds  /*i: Python tuple of keywords */
                            );

/*
 * This function returns the grouping and quality arrays that represent an
 * array of numChans elements in which the groups are each grpWidth elements
 * wide.
 */
static PyObject *grpBinWidth(PyObject *self,	/*i: Used by Python */
                             PyObject *args,	/*i: Python tuple of the arguments */
                             PyObject *kwds  /*i: Python tuple of keywords */
                             );

/*
 * This function returnes the number of channels (i.e. elements) in each
 * group. The return value is an array whose length equals that of the input
 * data (the dataArray argument) and each element within a group contains the
 * same value.
 */
static PyObject *grpGetChansPerGroup(PyObject *self,	/*i: Used by Python */
                                     PyObject *args,	/*i: Python tuple of the arguments */
                                     PyObject *kwds  /*i: Python tuple of keywords */
                                     );

/*
 * This function applies the grouping information from the grouping parameter
 * to the dataArray parameter. The return value is an array whose length
 * equals that of the input data (the dataArray argument) and each element
 * within a group contains the same value.
 */
static PyObject *grpGetGroupSum(PyObject *self,	/*i: Used by Python */
                                PyObject *args,	/*i: Python tuple of the arguments */
                                PyObject *kwds  /*i: Python tuple of keywords */
                                );

/*
 * his function calculates which group each element in the input array
 * belongs to, where the groups are numbered from 1. The return value is
 * an array whose length equals that of the input data (the grouping argument)
 * and each element within a group contains the same value.
 */
static PyObject *grpGetGroupNum(PyObject *self,	/*i: Used by Python */
                                PyObject *args,	/*i: Python tuple of the arguments */
                                PyObject *kwds  /*i: Python tuple of keywords */
                                );

/*
 * In this routine, groups are created when the absolute value of the slope
 * of the input data (the axisArray and binArray arguments) is less than
 * the threshold value (the slope argument).
 */
static PyObject *grpMaxSlope(PyObject *self,	/*i: Used by Python */
                             PyObject *args,	/*i: Python tuple of the arguments */
                             PyObject *kwds  /*i: Python tuple of keywords */
                             );

/*
 * In this routine, groups are created when the absolute value of the slope of
 * the input data (the axisArray and binArray arguments) is more than the
 * threshold value (the slope argument).
 */
static PyObject *grpMinSlope(PyObject *self,	/*i: Used by Python */
                             PyObject *args,	/*i: Python tuple of the arguments */
                             PyObject *kwds  /*i: Python tuple of keywords */
                             );

/*
 * This function returns the grouping and quality arrays that represent an
 * array of numChans elements grouped into numGroups groups.
 */
static PyObject *grpNumBins(PyObject *self,	/*i: Used by Python */
                            PyObject *args,	/*i: Python tuple of the arguments */
                            PyObject *kwds  /*i: Python tuple of keywords */
                            );

/*
 * This function returns the grouping and quality arrays that represent
 * the input data (countsArray) after it has been grouped so that each
 * group contains at least numCounts counts. The optional parameters
 * maxLength and tabStops represent the maximum number of elements that
 * can be combined and an array representing those elements that should
 * be ignored respectively.
 */
static PyObject *grpNumCounts(PyObject *self,	/*i: Used by Python */
															PyObject *args,	/*i: Python tuple of the arguments */
															PyObject *kwds  /*i: Python tuple of keywords */
															);

/*
 * This function returns the grouping and quality arrays that represent
 * the input data (countsArray) after it has been grouped so that the signal
 * to noise of each group is at least equal to the snr parameter. The
 * optional parameters maxLength and tabStops represent the maximum number
 * of elements that can be combined into a group and an array representing
 * those elements that should be ignored respectively. The errorCol array
 * gives the error for each element of the original array: if it is not
 * supplied then the error is taken to be the square root of the element value.
 */
static PyObject *grpSnr(PyObject *self,	/*i: Used by Python */
												PyObject *args, /*i: Python tuple of the arguments */
												PyObject *kwds  /*i: Python tuple of keywords */
												);
