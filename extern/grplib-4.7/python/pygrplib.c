/*_C_INSERT_SAO_COPYRIGHT_HERE_(2007-20012)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

/* H*****************************************************************
 *
 * FILE NAME:  pygrplib.c
 *
 * DEVELOPMENT: tools
 *
 * DESCRIPTION:  Provides a c-extension for the Python language to
 * the grplib library of functions
 *
 *
 *
 * REVISION HISTORY:
 *
 * Ref. No.         Date
 ----------       -----
 0.1              April2007 	File Created
 H***************************************************************** */

#include "pygrplib.h"
#include "grplib.h"
#include "grp_priv.h"

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Utilities for Python Interaction
 * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Set up the methods table
 */
static PyMethodDef
    groupMethods[] =
        {
         /* Each entry in the groupMethods array is a PyMethodDef structure containing
          * 1) the Python name,
          * 2) the C-function that implements the function,
          * 3) flags indicating whether or not keywords are accepted for this function,
          * and 4) The docstring for the function.
          */
          {
           "grpAdaptive",
           (PyCFunction)grpAdaptive,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpAdaptive(countsArray, minCounts [,maxLength, tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpAdaptiveSnr",
           (PyCFunction)grpAdaptiveSnr,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpAdaptiveSnr(countsArray, snr [,maxLength, tabStops, errorCol]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpBin",
           (PyCFunction)grpBin,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpBin(dataArray, binLowArray, binHighArray [,tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpBinFile",
           (PyCFunction)grpBinFile,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpBinFile(dataArray, fDataArray, fGroupingCol, fQualCol [,tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpBinWidth",
           (PyCFunction)grpBinWidth,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpBinWidth(numChans, binWidth [, tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpGetChansPerGroup",
           (PyCFunction)grpGetChansPerGroup,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  chanspergrp = grpGetChansPerGroup(groupCol) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpGetGroupSum",
           (PyCFunction)grpGetGroupSum,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  grpdata = grpGetGroupSum(dataArray, groupCol) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpGetGroupNum",
           (PyCFunction)grpGetGroupNum,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  grpnum = grpGetGroupNum(groupCol) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpMaxSlope",
           (PyCFunction)grpMaxSlope,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpMaxSlope(dataArray, binArray, slope [,maxLength, tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpMinSlope",
           (PyCFunction)grpMinSlope,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpMinSlope(dataArray, binArray, slope [,maxLength, tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpNumBins",
           (PyCFunction)grpNumBins,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpNumBins(numChans, numBins [, tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpNumCounts",
           (PyCFunction)grpNumCounts,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpNumCounts(countsArray, numCounts [,maxLength,tabStops]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {
           "grpSnr",
           (PyCFunction)grpSnr,
           METH_VARARGS | METH_KEYWORDS,
           "Example:  (grouping, quality) = grpSnr(countsArray, snr [,maxLength, tabStops,errorCol]) \n\
     The keywords for the function are the parameters shown in 'Example.'"},
          {NULL, NULL, 0, NULL} /* Sentinel - marks the end of this structure */
        };/*end...groupMethods */

/*
 * Initialize the module
 * */
void initgroup(void);
void initgroup()
{
  (void)Py_InitModule("group", groupMethods);
  import_array(); /* Must be present for NumPy.  Called first after above line. */
}

/* Error Message Format Strings */
char groupmsg[1024];
#define GROUP_GENERIC_MSG        "%s() %s"
#define GROUP_TYPE_ERROR_MSG     "%s() Could not parse input arguments, please check input for correct type(s)"
#define GROUP_OPEN_FILE_MSG      "%s() Could not open parameter file"
#define GROUP_MEMORY_MSG         "%s() Could not allocate memory"
#define GROUP_CREATE_MSG         "%s() Could not create array object"
#define GROUP_DIFF_LENGTH_MSG    "%s() %s have differing length"
#define GROUP_INCOMP_TYPE_MSG    "%s() %s is an incompatible type"
#define GROUP_GT_ZERO_MSG        "%s() %s values must be > zero"

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Function Definitions
 * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * This function returns the grouping and quality arrays that represent the
 * input data (countsArray) after it has been adaptively grouped so that
 * each group contains at least numCounts counts.
 */
static PyObject *grpAdaptive(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpAdaptive";
  double maxLength = 0; /* number of elements that can be combined into a group */
  double minCounts = 0; /* how many counts to contain in each group */
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_countsArray = NULL; /* countsArray in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_countsArray = NULL;
  char *arr_tabStops = NULL;
  int stride_countsArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_countsArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"countsArray", "minCounts", "maxLength", "tabStops",
                           NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|dO!", kwlist,
                                   &PyArray_Type, &py_countsArray, &minCounts, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_countsArray->descr->type_num, NPY_DOUBLE))
    {
      py_countsArray
          = (PyArrayObject *)PyArray_Cast(py_countsArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_countsArray = PyArray_STRIDE(py_countsArray, 0);
      arr_countsArray = PyArray_BYTES(py_countsArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName,
              (char*)"The countsArray");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */
  if (minCounts <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"minCounts");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_countsArray->dimensions[0]; /* the number of channels is the size of the py_countsArray */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_countsArray = (double *)calloc(numChans, sizeof(double));
  if (c_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_countsArray[ii] = *((double*)((void*)(arr_countsArray
        + stride_countsArray * ii)));
  }

  /* Function called from grplib.c */
  isError = grp_do_adaptive(c_countsArray, numChans, minCounts, groupCol,
                            qualCol, c_tabStops, maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol);/* free the allocated memory */
  free(qualCol);
  free(c_countsArray);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpAdaptive*/

/*
 * This function returns the grouping and quality arrays that represent the
 * input data (countsArray) after it has been adaptively grouped so that the
 * signal to noise of each group is at least equal to the snr parameter. The
 * errorCol array gives the error for each element of the original array: if
 * it is not supplied then the error is taken to be the square root of the
 * element value.
 */
static PyObject *grpAdaptiveSnr(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpAdaptiveSnr";
  double snr = 0; /* signal to noise parameter */
  double maxLength = 0; /* number of elements that can be combined into a group */
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_countsArray = NULL; /* countsArray in c-style array */
  double *c_errorCol = NULL; /* errorCol in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  short useErrCols = 0; /* value indicating if a errorCol argument was passed to the function */
  int numErrs = 0; /* number of errors in errorCol */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_countsArray = NULL;
  char *arr_errCol = NULL;
  char *arr_tabStops = NULL;
  int stride_countsArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_errCol = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_countsArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_errorCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /* The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grp_do_snr in grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grp_do_snr in grplib.c */

  static char *kwlist[] = {"countsArray", "snr", "maxLength", "tabStops",
                           "errorCol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|dO!O!", kwlist,
                                   &PyArray_Type, &py_countsArray, &snr, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops,
                                   &PyArray_Type, &py_errorCol /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_countsArray->descr->type_num, NPY_DOUBLE))
    {
      py_countsArray
          = (PyArrayObject *)PyArray_Cast(py_countsArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_countsArray = PyArray_STRIDE(py_countsArray, 0);
      arr_countsArray = PyArray_BYTES(py_countsArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName,
              (char*)"The countsArray");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */
  if (snr <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  if (py_errorCol != NULL)
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_errorCol->descr->type_num, NPY_DOUBLE))
    {
      py_errorCol = (PyArrayObject *)PyArray_Cast(py_errorCol, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_errCol = PyArray_STRIDE(py_errorCol, 0);
      arr_errCol = PyArray_BYTES(py_errorCol);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The errorCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */

    useErrCols = 1; /* set value to true since we have a errorCol array */
    numErrs = py_errorCol->dimensions[0]; /* the number of tabs is the size of the py_errorCol */

    c_errorCol = (double *)calloc(numErrs, sizeof(double));
    if (c_errorCol == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numErrs; ii++)
    {
      c_errorCol[ii] = *((double*)((void*)(arr_errCol + stride_errCol * ii)));
    }
  }/*end if... */

  numChans = py_countsArray->dimensions[0]; /* the number of channels is the size of the py_countsArray */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */
  if (useErrCols && (numErrs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The errorCol and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!useErrCols)
  {
    c_errorCol = (double *)calloc(numChans, sizeof(double));
    if (c_errorCol == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  ii = 0;
  while (!useErrCols && (ii < numChans))
  {
    c_errorCol[ii] = 1.0; /*fill errorCol with 1's */
    ii++;
  }/*end while... */

  c_countsArray = (double *)calloc(numChans, sizeof(double));
  if (c_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_countsArray[ii] = *((double*)((void*)(arr_countsArray
        + stride_countsArray * ii)));
  }

  /* Function called from grplib.c */
  isError = grp_do_adaptive_snr(c_countsArray, numChans, snr, groupCol,
                                qualCol, c_tabStops, c_errorCol, useErrCols,
                                maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol);/* free the allocated memory */
  free(qualCol);
  free(c_countsArray);
  free(c_errorCol);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpAdaptiveSnr*/

/*
 * This function returns the grouping and quality arrays for a set of groups
 * defined by the low (binLowArray) and high (binHighArray) boundaries when
 * applied to the axis values of the data (axisArray).
 */
static PyObject *grpBin(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpBin";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numBins = -1; /* number of bins provided by binlow and binhigh */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  int colRealFlag = 0; /* value determining if dataArray values are ints or reals */
  double *c_dataArray = NULL; /* dataArray in c-style array */
  double *c_binLowArray = NULL; /* binLowArray in c-style array */
  double *c_binHighArray = NULL; /* binHighArray in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_dataArray = NULL;
  char *arr_binLowArray = NULL;
  char *arr_binHighArray = NULL;
  char *arr_tabStops = NULL;
  int stride_dataArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_binLowArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_binHighArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_dataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_binLowArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_binHighArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"dataArray", "binLowArray", "binHighArray",
                           "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!|O!", kwlist,
                                   &PyArray_Type, &py_dataArray, &PyArray_Type,
                                   &py_binLowArray, /* mandatory args */
                                   &PyArray_Type, &py_binHighArray, /* mandatory args */
                                   &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((py_dataArray == NULL) || (py_binLowArray == NULL) || (py_binHighArray
      == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_dataArray->descr->type_num, NPY_DOUBLE)
        && (PyArray_CanCastSafely(py_binLowArray->descr->type_num, NPY_DOUBLE))
        && (PyArray_CanCastSafely(py_binHighArray->descr->type_num, NPY_DOUBLE)))
    {
      /* determine if dataArray values are ints or reals before casting*/
      /* NOTE: int value for NPY_FLOAT  = 11
       * 			 int value for NPY_DOUBLE = 12
       * 			 int value for NPY_LONG   = 7 (which corresponds to c-type INT) */
      if (((py_dataArray->descr->type_num) == NPY_FLOAT)
          || ((py_dataArray->descr->type_num) == NPY_DOUBLE))
      {
        colRealFlag = 1;
      }/*end if... */

      py_dataArray = (PyArrayObject *)PyArray_Cast(py_dataArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_dataArray = PyArray_STRIDE(py_dataArray, 0);
      arr_dataArray = PyArray_BYTES(py_dataArray);

      py_binLowArray
          = (PyArrayObject *)PyArray_Cast(py_binLowArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_binLowArray = PyArray_STRIDE(py_binLowArray, 0);
      arr_binLowArray = PyArray_BYTES(py_binLowArray);

      py_binHighArray = (PyArrayObject *)PyArray_Cast(py_binHighArray,
                                                      NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_binHighArray = PyArray_STRIDE(py_binHighArray, 0);
      arr_binHighArray = PyArray_BYTES(py_binHighArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The Array");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_dataArray->dimensions[0];
  numBins = py_binLowArray->dimensions[0]; /* the number of channels is the size of the py_binLowArray */
  /* check to see if binlow has same size as binhigh */
  if (py_binLowArray->dimensions[0] != py_binHighArray->dimensions[0])
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"binLowArray and binHighArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_dataArray = (double *)calloc(numChans, sizeof(double));
  c_binLowArray = (double *)calloc(numChans, sizeof(double));
  c_binHighArray = (double *)calloc(numChans, sizeof(double));
  if (c_dataArray == NULL || c_binLowArray == NULL || c_binHighArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }

  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_binLowArray[ii] = *((double*)((void*)(arr_binLowArray
        + stride_binLowArray * ii)));
    c_dataArray[ii]
        = *((double*)((void*)(arr_dataArray + stride_dataArray * ii)));
    c_binHighArray[ii] = *((double*)((void*)(arr_binHighArray
        + stride_binHighArray * ii)));
  }

  /* Function called from grplib.c */
  isError = grp_do_bin(c_dataArray, numChans, c_binLowArray, c_binHighArray,
                       numBins, groupCol, qualCol, c_tabStops, NULL, 0,
                       colRealFlag);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /* free the allocated memory */
  free(qualCol);
  free(c_dataArray);
  free(c_binLowArray);
  free(c_binHighArray);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpBin*/

/*
 * This function allows you to calculate the grouping information needed to
 * group the input data (the axisArray array) to match the grouping of another
 * dataset
 */
static PyObject *grpBinFile(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpBinFile";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int fNumChans; /* number of channels in the file data */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  int colRealFlag = 0; /* value determining if dataArray values are ints or reals */
  double *c_dataArray = NULL; /* dataArray in c-style array */
  double *c_fDataArray = NULL; /* fDataArray in c-style array */
  short *c_fGroupCol = NULL; /* short int column of grouping data argument */
  short *c_fQualCol = NULL; /* short int column of quality data argument */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_dataArray = NULL;
  char *arr_fDataArray = NULL;
  char *arr_fGroupCol = NULL;
  char *arr_fQualCol = NULL;
  char *arr_tabStops = NULL;
  int stride_dataArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_fDataArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_fGroupCol = 0; /* used to find next value in non-contiguous arrays */
  int stride_fQualCol = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_dataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_fDataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_fGroupCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_fQualCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"dataArray", "fDataArray", "fGroupCol", "fQualCol",
                           "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!O!|O!", kwlist,
                                   &PyArray_Type, &py_dataArray, &PyArray_Type,
                                   &py_fDataArray, /* mandatory args */
                                   &PyArray_Type, &py_fGroupCol, &PyArray_Type,
                                   &py_fQualCol, /* mandatory args */
                                   &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((py_dataArray == NULL) || (py_fDataArray == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if ((PyArray_CanCastSafely(py_dataArray->descr->type_num, NPY_DOUBLE))
        && (PyArray_CanCastSafely(py_fDataArray->descr->type_num, NPY_DOUBLE)))
    {
      /* determine if dataArray values are ints or reals before casting */
      /* NOTE: int value for NPY_FLOAT  = 11
       * 			 int value for NPY_DOUBLE = 12
       * 			 int value for NPY_LONG   = 7 (which corresponds to c-type INT) */
      if (((py_dataArray->descr->type_num) == NPY_FLOAT)
          || ((py_dataArray->descr->type_num) == NPY_DOUBLE))
      {
        colRealFlag = 1;
      }/*end if... */

      py_dataArray = (PyArrayObject *)PyArray_Cast(py_dataArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_dataArray = PyArray_STRIDE(py_dataArray, 0);
      arr_dataArray = PyArray_BYTES(py_dataArray);

      py_fDataArray = (PyArrayObject *)PyArray_Cast(py_fDataArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_fDataArray = PyArray_STRIDE(py_fDataArray, 0);
      arr_fDataArray = PyArray_BYTES(py_fDataArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The Array");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */

  if ((py_fGroupCol == NULL) || (py_fQualCol == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    if ((py_fGroupCol->descr->type_num >= 17) || (py_fQualCol->descr->type_num
        >= 17))
    {/*types 17 and above include strings and other non-numerical values */
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName,
              (char*)"The groupCol or qualCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_fGroupCol = (PyArrayObject *)PyArray_Cast(py_fGroupCol, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_fGroupCol = PyArray_STRIDE(py_fGroupCol, 0);
    arr_fGroupCol = PyArray_BYTES(py_fGroupCol);

    py_fQualCol = (PyArrayObject *)PyArray_Cast(py_fQualCol, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_fQualCol = PyArray_STRIDE(py_fQualCol, 0);
    arr_fQualCol = PyArray_BYTES(py_fQualCol);
  }/*end else... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_dataArray->dimensions[0];
  fNumChans = py_fDataArray->dimensions[0]; /* the number of channels is the size of the py_fDataArray */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_dataArray = (double *)calloc(numChans, sizeof(double));
  c_fDataArray = (double *)calloc(numChans, sizeof(double));
  c_fGroupCol = (short *)calloc(numChans, sizeof(short));
  c_fQualCol = (short *)calloc(numChans, sizeof(short));
  if (c_dataArray == NULL || c_fDataArray == NULL || c_fGroupCol == NULL
      || c_fQualCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }

  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_dataArray[ii]
        = *((double*)((void*)(arr_dataArray + stride_dataArray * ii)));
    c_fDataArray[ii] = *((double*)((void*)(arr_fDataArray + stride_fDataArray
        * ii)));
    c_fGroupCol[ii]
        = *((short*)((void*)(arr_fGroupCol + stride_fGroupCol * ii)));
    c_fQualCol[ii] = *((short*)((void*)(arr_fQualCol + stride_fQualCol * ii)));
  }

  /* call the function from grplib.c */
  isError = grp_do_bin_file(c_dataArray, numChans, groupCol, qualCol,
                            c_tabStops, c_fDataArray, fNumChans, c_fGroupCol,
                            c_fQualCol, colRealFlag, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /* free the allocated memory */
  free(qualCol);
  free(c_dataArray);
  free(c_fDataArray);
  free(c_fGroupCol);
  free(c_fQualCol);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpBinFile*/

/*
 * This function returns the grouping and quality arrays that represent an
 * array of numChans elements in which the groups are each grpWidth elements
 * wide.
 */
static PyObject *grpBinWidth(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpBinWidth";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  long numChans = 0; /* number of channels in groupCol and qualCol */
  long binWidth = 0; /* number of bins */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_tabStops = NULL;
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"numChans", "binWidth", "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O!", kwlist, &numChans,
                                   &binWidth, /* mandatory args */
                                   &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);;
    return NULL;
  }/*end if... */
  if ((numChans <= 0) || (binWidth <= 0))
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_GENERIC_MSG, funcName,
            (char*)"The number of tab stops and number of channels "
              "specified in the argument list have different sizes");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* Function called from grplib.c */
  isError = grp_do_bin_width(numChans, binWidth, groupCol, qualCol, c_tabStops,
                             NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /*free allocated memory */
  free(qualCol);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpBinWidth*/

/*
 * This function returns the number of channels (i.e. elements) in each
 * group. The return value is an array whose length equals that of the input
 * data (the dataArray argument) and each element within a group contains the
 * same value.
 */
static PyObject *grpGetChansPerGroup(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpGetChansPerGroup";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  long numChans = 0; /* number of channels in groupCol and qualCol */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  long *chansPerGrpCol = NULL; /* used to store values when converting from c-style array to numpy array */
  double *groupData = NULL; /* value used to store values while creating output */
  short *c_groupCol = NULL; /* short int of c_groupCol */

  /* Handles case if array is not contiguous in memory */
  char *arr_groupCol = NULL;
  int stride_groupCol = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_groupCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *chansPerGrp = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"groupCol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, &PyArray_Type,
                                   &py_groupCol))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_groupCol == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    if (py_groupCol->descr->type_num >= 17)/*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The groupCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_groupCol = (PyArrayObject *)PyArray_Cast(py_groupCol, NPY_SHORT);

    /* Handles case if array is not contiguous in memory */
    stride_groupCol = PyArray_STRIDE(py_groupCol, 0);
    arr_groupCol = PyArray_BYTES(py_groupCol);
  }/*end else... */

  numChans = py_groupCol->dimensions[0];

  /* allocate memory for arrays */
  chansPerGrpCol = (long *)calloc(numChans, sizeof(long));
  if (chansPerGrpCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */

  c_groupCol = (short *)calloc(numChans, sizeof(short));
  if (c_groupCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */

  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_groupCol[ii] = *((short*)((void*)(arr_groupCol + stride_groupCol * ii)));
  }

  /* Function called from grplib.c */
  isError = set_chans_per_grp(c_groupCol, chansPerGrpCol, numChans);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  chansPerGrp = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if (NULL == chansPerGrp)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */

  groupData = DDATA(chansPerGrp);

  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = chansPerGrpCol[ii]; /*grab the data from groupCol and place in grouping data */
  }/*end for... */

  free(chansPerGrpCol); /*free allocated memory */
  free(c_groupCol);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("O", PyArray_Return(chansPerGrp));
}/*end...grpGetChansPerGroup*/

/*
 * This function applies the grouping information from the grouping parameter
 * to the dataArray parameter. The return value is an array whose length
 * equals that of the input data (the dataArray argument) and each element
 * within a group contains the same value.
 */
static PyObject *grpGetGroupSum(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpGetGroupSum";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  long numChans = 0; /* number of channels in groupCol and qualCol */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *grpDataCol = NULL; /* used to store values when converting from c-style array to numpy array */
  double *groupData = NULL; /* value used to store values while creating output */
  short *c_groupCol = NULL; /* short int of c_groupCol */
  double *c_dataArray = NULL; /*  The data array as a c-style array */

  /* Handles case if array is not contiguous in memory */
  char *arr_dataArray = NULL;
  char *arr_groupCol = NULL;
  int stride_groupCol = 0; /* used to find next value in non-contiguous arrays */
  int stride_dataArray = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_dataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_groupCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *grpData = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"dataArray", "groupCol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!", kwlist, &PyArray_Type,
                                   &py_dataArray, &PyArray_Type, &py_groupCol))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((py_dataArray == NULL) || (py_groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_dataArray->descr->type_num, NPY_DOUBLE))
    {
      py_dataArray = (PyArrayObject *)PyArray_Cast(py_dataArray, NPY_DOUBLE);
      stride_dataArray = PyArray_STRIDE(py_dataArray, 0);
      arr_dataArray = PyArray_BYTES(py_dataArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The dataArray");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */

    if (py_groupCol->descr->type_num >= 17)/*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The groupCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_groupCol = (PyArrayObject *)PyArray_Cast(py_groupCol, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_groupCol = PyArray_STRIDE(py_groupCol, 0);
    arr_groupCol = PyArray_BYTES(py_groupCol);
  }/*end else... */

  numChans = py_dataArray->dimensions[0];

  c_groupCol = (short *)calloc(numChans, sizeof(short));
  c_dataArray = (double *)calloc(numChans, sizeof(double));
  if (c_groupCol == NULL || c_dataArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_groupCol[ii] = *((short*)((void*)(arr_groupCol + stride_groupCol * ii)));
    c_dataArray[ii]
        = *((double*)((void*)(arr_dataArray + stride_dataArray * ii)));
  }

  /* allocate memory for arrays */
  grpDataCol = (double *)calloc(numChans, sizeof(double));
  if (grpDataCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */

  /* Function called from grplib.c */
  isError = set_grp_data(c_dataArray, c_groupCol, grpDataCol, numChans);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grpData = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if (NULL == grpData)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */

  groupData = DDATA(grpData);

  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = grpDataCol[ii]; /*grab the data from groupCol and place in grouping data */
  }/*end for... */

  free(grpDataCol);/*free allocated memory */
  free(c_dataArray);
  free(c_groupCol);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("O", PyArray_Return(grpData));
}/*end...grpGetGroupSum*/

/*
 * his function calculates which group each element in the input array
 * belongs to, where the groups are numbered from 1. The return value is
 * an array whose length equals that of the input data (the grouping argument)
 * and each element within a group contains the same value.
 */
static PyObject *grpGetGroupNum(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpGetGroupNum";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  long numChans = 0; /* number of channels in groupCol and qualCol */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  long *grpNumCol = NULL; /* used to store values when converting from c-style array to numpy array */
  double *groupData = NULL; /* value used to store values while creating output */
  short *c_groupCol = NULL; /* short int of c_groupCol */

  /* Handles case if array is not contiguous in memory */
  char *arr_groupCol = NULL;
  int stride_groupCol = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_groupCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *grpNum = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"groupCol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", kwlist, &PyArray_Type,
                                   &py_groupCol))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_groupCol == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    if ((py_groupCol->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The groupCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_groupCol = (PyArrayObject *)PyArray_Cast(py_groupCol, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_groupCol = PyArray_STRIDE(py_groupCol, 0);
    arr_groupCol = PyArray_BYTES(py_groupCol);
  }/*end else... */

  numChans = py_groupCol->dimensions[0];

  /* allocate memory for arrays */
  grpNumCol = (long *)calloc(numChans, sizeof(long));
  if (grpNumCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */

  c_groupCol = (short *)calloc(numChans, sizeof(short));
  if (c_groupCol == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_groupCol[ii] = *((short*)((void*)(arr_groupCol + stride_groupCol * ii)));
  }

  /* Function called from grplib.c */
  isError = set_grp_num(c_groupCol, grpNumCol, numChans);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grpNum = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if (NULL == grpNum)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */

  groupData = DDATA(grpNum);

  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = grpNumCol[ii]; /*grab the data from groupCol and place in grouping data */
  }/*end for... */

  free(grpNumCol); /*free allocated memory */
  free(c_groupCol);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("O", PyArray_Return(grpNum));
}/*end...grpGetGroupNum*/

/*
 * In this routine, groups are created when the absolute value of the slope
 * of the input data (the axisArray and binArray arguments) is less than
 * the threshold value (the slope argument).
 */
static PyObject *grpMaxSlope(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpMaxSlope";
  double maxLength = 0; /* number of elements that can be combined into a group */
  double slope = 0; /* the value of the slope to limit the grouping */
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numBins; /* number of bins in the binArray */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_dataArray = NULL; /* dataArray in c-style array */
  double *c_binArray = NULL; /* binArray in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_dataArray = NULL;
  char *arr_binArray = NULL;
  char *arr_tabStops = NULL;
  int stride_dataArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_binArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_dataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_binArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"dataArray", "binArray", "slope", "maxLength",
                           "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!d|dO!", kwlist,
                                   &PyArray_Type, &py_dataArray, &PyArray_Type,
                                   &py_binArray, &slope, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((py_dataArray == NULL) || (py_binArray == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_dataArray->descr->type_num, NPY_DOUBLE)
        && (PyArray_CanCastSafely(py_binArray->descr->type_num, NPY_DOUBLE)))
    {
      py_dataArray = (PyArrayObject *)PyArray_Cast(py_dataArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_dataArray = PyArray_STRIDE(py_dataArray, 0);
      arr_dataArray = PyArray_BYTES(py_dataArray);
      py_binArray = (PyArrayObject *)PyArray_Cast(py_binArray, NPY_DOUBLE);

      /* Handles case if array is not contiguous in memory */
      stride_binArray = PyArray_STRIDE(py_binArray, 0);
      arr_binArray = PyArray_BYTES(py_binArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The Array");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */
  if (slope <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    //c_tabStops  = SDATA(py_tabStops);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_dataArray->dimensions[0]; /* the number of channels is the size of the py_dataArray */
  numBins = py_binArray->dimensions[0]; /* the number of bins is the size of the py_binArray */

  if (numBins != numChans)
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The binArray and dataArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and dataArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_dataArray = (double *)calloc(numChans, sizeof(double));
  c_binArray = (double *)calloc(numChans, sizeof(double));
  if (c_dataArray == NULL || c_binArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_dataArray[ii]
        = *((double*)((void*)(arr_dataArray + stride_dataArray * ii)));
    c_binArray[ii] = *((double*)((void*)(arr_binArray + stride_binArray * ii)));
  }

  isError = grp_do_max_slope(c_dataArray, c_binArray, numChans, slope,
                             groupCol, qualCol, c_tabStops, maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /*free allocated memory */
  free(qualCol);
  free(c_tabStops);
  free(c_dataArray);
  free(c_binArray);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpMaxSlope*/

/*
 * In this routine, groups are created when the absolute value of the slope of
 * the input data (the axisArray and binArray arguments) is more than the
 * threshold value (the slope argument).
 */
static PyObject *grpMinSlope(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpMinSlope";
  double maxLength = 0; /* number of elements that can be combined into a group */
  double slope = 0; /* the value of the slope to limit the grouping */
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numBins; /* number of bins in the binArray */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_dataArray = NULL; /* dataArray in c-style array */
  double *c_binArray = NULL; /* binArray in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_dataArray = NULL;
  char *arr_binArray = NULL;
  char *arr_tabStops = NULL;
  int stride_dataArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_binArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_dataArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_binArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"dataArray", "binArray", "slope", "maxLength",
                           "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!d|dO!", kwlist,
                                   &PyArray_Type, &py_dataArray, &PyArray_Type,
                                   &py_binArray, &slope, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((py_dataArray == NULL) || (py_binArray == NULL))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_dataArray->descr->type_num, NPY_DOUBLE)
        && (PyArray_CanCastSafely(py_binArray->descr->type_num, NPY_DOUBLE)))
    {
      py_dataArray = (PyArrayObject *)PyArray_Cast(py_dataArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_dataArray = PyArray_STRIDE(py_dataArray, 0);
      arr_dataArray = PyArray_BYTES(py_dataArray);
      py_binArray = (PyArrayObject *)PyArray_Cast(py_binArray, NPY_DOUBLE);

      /* Handles case if array is not contiguous in memory */
      stride_binArray = PyArray_STRIDE(py_binArray, 0);
      arr_binArray = PyArray_BYTES(py_binArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The Array");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */
  if (slope <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    //c_tabStops  = SDATA(py_tabStops);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_dataArray->dimensions[0]; /* the number of channels is the size of the py_dataArray */
  numBins = py_binArray->dimensions[0]; /* the number of bins is the size of the py_binArray */

  if (numBins != numChans)
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The binArray and dataArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and dataArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_dataArray = (double *)calloc(numChans, sizeof(double));
  c_binArray = (double *)calloc(numChans, sizeof(double));
  if (c_dataArray == NULL || c_binArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_dataArray[ii]
        = *((double*)((void*)(arr_dataArray + stride_dataArray * ii)));
    c_binArray[ii] = *((double*)((void*)(arr_binArray + stride_binArray * ii)));
  }

  isError = grp_do_min_slope(c_dataArray, c_binArray, numChans, slope,
                             groupCol, qualCol, c_tabStops, maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /*free allocated memory */
  free(qualCol);
  free(c_tabStops);
  free(c_dataArray);
  free(c_binArray);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpMinSlope*/

/*
 * This function returns the grouping and quality arrays that represent an
 * array of numChans elements grouped into numGroups groups.
 */
static PyObject *grpNumBins(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpNumBins";
  int ii; /* loop variable */
  int isError; /* captures return value from grplib.c */
  long numChans = 0; /* number of channels in groupCol and qualCol */
  long numBins = 0; /* number of bins */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_tabStops = NULL;
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_tabStops = NULL; /*  The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grplib.c */

  static char *kwlist[] = {"numChans", "numBins", "tabStops", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O!", kwlist, &numChans,
                                   &numBins, /* mandatory args */
                                   &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if ((numChans <= 0) || (numBins <= 0))
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }
    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and numChans");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* Function called from grplib.c */
  isError = grp_do_num_bins(numChans, numBins, groupCol, qualCol, c_tabStops,
                            NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /*free allocated memory */
  free(qualCol);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpNumBins*/

/*
 * This function returns the grouping and quality arrays that represent
 * the input data (countsArray) after it has been grouped so that each
 * group contains at least numCounts counts. The optional parameters
 * maxLength and tabStops represent the maximum number of elements that
 * can be combined and an array representing those elements that should
 * be ignored respectively.
 */
static PyObject *grpNumCounts(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpNumCounts";
  double maxLength = 0; /* number of elements that can be combined into a group */
  double numCounts = 0; /* how many counts to contain in each group */
  int ii; /* loop variable */
  int isError; /* captures return value from grp_do_num_counts in grplib.c */
  int numChans; /* number of channels in groupCol and qualCol */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_countsArray = NULL; /* countsArray in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_countsArray = NULL;
  char *arr_tabStops = NULL;
  int stride_countsArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_countsArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *grouping = NULL; /* The result obtained from grp_do_num_counts in grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grp_do_num_counts in grplib.c */

  static char *kwlist[] = {"countsArray", "numCounts", "maxLength", "tabStops",
                           NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|dO!", kwlist,
                                   &PyArray_Type, &py_countsArray, &numCounts, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_countsArray->descr->type_num, NPY_DOUBLE))
    {
      py_countsArray
          = (PyArrayObject *)PyArray_Cast(py_countsArray, NPY_DOUBLE);

      /* Handles case if array is not contiguous in memory */
      stride_countsArray = PyArray_STRIDE(py_countsArray, 0);
      arr_countsArray = PyArray_BYTES(py_countsArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName,
              (char*)"The countsArray");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */

  if (numCounts <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }

    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);

    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  numChans = py_countsArray->dimensions[0]; /* the number of channels is the size of the py_countsArray */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  c_countsArray = (double *)calloc(numChans, sizeof(double));
  if (c_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }

  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_countsArray[ii] = *((double*)((void*)(arr_countsArray
        + stride_countsArray * ii)));
  }

  /* Function called from grplib.c */
  isError = grp_do_num_counts(c_countsArray, numChans, numCounts, groupCol,
                              qualCol, c_tabStops, maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /* free the allocated memory */
  free(qualCol);
  free(c_countsArray);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpNumCounts*/

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
static PyObject *grpSnr(PyObject *self, /*i: Used by Python */
PyObject *args, /*i: Python tuple of the arguments */
PyObject *kwds /*i: Python tuple of keywords */)
{
  char* funcName = "grpSnr";
  double snr = 0; /* signal to noise parameter */
  double maxLength = 0; /* number of elements that can be combined into a group */
  int ii = 0; /* loop variable */
  int isError = 0; /* captures return value from grplib.c */
  int numChans = 0; /* number of channels in groupCol and qualCol */
  int numTabs = 0; /* number of tabs in tabStop */
  npy_intp dims[1]; /* the dimensions of the arrays */
  npy_intp typenum; /* the typenum */
  double *c_countsArray = NULL; /* countsArray in c-style array */
  double *c_errorCol = NULL; /* errorCol in c-style array */
  double *groupData = NULL; /* used to store values when converting from c-style array to numpy array */
  double *qualData = NULL; /* used to store values when converting from c-style array to numpy array */
  short *groupCol = NULL; /* the GROUPING column */
  short *qualCol = NULL; /* the QUALITY column */
  short *c_tabStops = NULL; /* elements that should be ignored */
  short useErrCols = 0; /* value indicating if a errorCol argument was passed to the function */
  int numErrs = 0; /* number of errors in errorCol */
  int isTabStops = 0; /* a tabStop argument is passed in */

  /* Handles case if array is not contiguous in memory */
  char *arr_countsArray = NULL;
  char *arr_errCol = NULL;
  char *arr_tabStops = NULL;
  int stride_countsArray = 0; /* used to find next value in non-contiguous arrays */
  int stride_errCol = 0; /* used to find next value in non-contiguous arrays */
  int stride_tabStops = 0; /* used to find next value in non-contiguous arrays */

  PyArrayObject *py_countsArray = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_errorCol = NULL; /* Python Array Object that will be converted to a c-style array for processing */
  PyArrayObject *py_tabStops = NULL; /* The Python Object that will be turn into a numpy array Object */
  PyArrayObject *grouping = NULL; /* The result obtained from grp_do_snr in grplib.c */
  PyArrayObject *quality = NULL; /* The result obtained from grp_do_snr in grplib.c */

  static char *kwlist[] = {"countsArray", "snr", "maxLength", "tabStops",
                           "errorCol", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!d|dO!O!", kwlist,
                                   &PyArray_Type, &py_countsArray, &snr, /* mandatory args */
                                   &maxLength, &PyArray_Type, &py_tabStops,
                                   &PyArray_Type, &py_errorCol /* optional keyword args*/))
  {
    sprintf(groupmsg, GROUP_TYPE_ERROR_MSG, funcName);
    PyErr_SetString(PyExc_TypeError, groupmsg);
    return NULL;
  }/*end if... */
  if (py_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  else
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_countsArray->descr->type_num, NPY_DOUBLE))
    {
      py_countsArray
          = (PyArrayObject *)PyArray_Cast(py_countsArray, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_countsArray = PyArray_STRIDE(py_countsArray, 0);
      arr_countsArray = PyArray_BYTES(py_countsArray);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName,
              (char*)"The countsArray");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
  }/*end else... */
  if (snr <= 0)
  {
    sprintf(groupmsg, GROUP_GT_ZERO_MSG, funcName, (char*)"Scalar");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  if (py_tabStops != NULL)/* if a tabStop array is present */
  {
    if ((py_tabStops->descr->type_num) >= 17) /*types 17 and above include strings and other non-numerical values */
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The tabStops");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end if... */
    py_tabStops = (PyArrayObject *)PyArray_Cast(py_tabStops, NPY_SHORT);
    /* Handles case if array is not contiguous in memory */
    stride_tabStops = PyArray_STRIDE(py_tabStops, 0);
    arr_tabStops = PyArray_BYTES(py_tabStops);

    numTabs = py_tabStops->dimensions[0]; /* the number of tabs is the size of the py_tabStops */
    isTabStops = 1; /* set value to true since we have a tabStop array */

    c_tabStops = (short *)calloc(numTabs, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numTabs; ii++)
    {
      c_tabStops[ii]
          = *((short*)((void*)(arr_tabStops + stride_tabStops * ii)));
    }
  }/*end if... */

  if (py_errorCol != NULL)
  {
    /* Make sure the arrays are of correct type */
    if (PyArray_CanCastSafely(py_errorCol->descr->type_num, NPY_DOUBLE))
    {
      py_errorCol = (PyArrayObject *)PyArray_Cast(py_errorCol, NPY_DOUBLE);
      /* Handles case if array is not contiguous in memory */
      stride_errCol = PyArray_STRIDE(py_errorCol, 0);
      arr_errCol = PyArray_BYTES(py_errorCol);
    }/*end if... */
    else
    {
      sprintf(groupmsg, GROUP_INCOMP_TYPE_MSG, funcName, (char*)"The errorCol");
      PyErr_SetString(PyExc_TypeError, groupmsg);
      return NULL;
    }/*end else... */
    useErrCols = 1; /* set value to true since we have a errorCol array */
    numErrs = py_errorCol->dimensions[0]; /* the number of tabs is the size of the py_errorCol */

    c_errorCol = (double *)calloc(numErrs, sizeof(double));
    if (c_errorCol == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }

    /* loop through all of the coordinates and put the return value in c_value */
    for (ii = 0; ii < numErrs; ii++)
    {
      c_errorCol[ii] = *((double*)((void*)(arr_errCol + stride_errCol * ii)));
    }
  }/*end if... */

  numChans = py_countsArray->dimensions[0]; /* the number of channels is the size of the py_countsArray */
  if (isTabStops && (numTabs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The tabStops and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */
  if (useErrCols && (numErrs != numChans))
  {
    sprintf(groupmsg, GROUP_DIFF_LENGTH_MSG, funcName,
            (char*)"The errorCol and countsArray");
    PyErr_SetString(PyExc_ValueError, groupmsg);
    return NULL;
  }/*end if... */

  /* allocate memory for arrays */
  groupCol = (short *)calloc(numChans, sizeof(short));
  qualCol = (short *)calloc(numChans, sizeof(short));
  if ((qualCol == NULL) || (groupCol == NULL))
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }/*end if... */
  if (!useErrCols)
  {
    c_errorCol = (double *)calloc(numChans, sizeof(double));
    if (c_errorCol == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */
  if (!isTabStops)
  {
    c_tabStops = (short *)calloc(numChans, sizeof(short));
    if (c_tabStops == NULL)
    {
      sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
      PyErr_SetString(PyExc_MemoryError, groupmsg);
      return NULL;
    }/*end if... */
  }/*end if... */

  ii = 0;
  while (!useErrCols && (ii < numChans))
  {
    c_errorCol[ii] = 1.0; /*fill errorCol with 1's */
    ii++;
  }/*end while... */

  c_countsArray = (double *)calloc(numChans, sizeof(double));
  if (c_countsArray == NULL)
  {
    sprintf(groupmsg, GROUP_MEMORY_MSG, funcName);
    PyErr_SetString(PyExc_MemoryError, groupmsg);
    return NULL;
  }
  /* loop through all of the coordinates and put the return value in c_value */
  for (ii = 0; ii < numChans; ii++)
  {
    c_countsArray[ii] = *((double*)((void*)(arr_countsArray
        + stride_countsArray * ii)));
  }

  /* Function called from grplib.c */
  isError = grp_do_snr(c_countsArray, numChans, snr, groupCol, qualCol,
                       c_tabStops, c_errorCol, useErrCols, maxLength, NULL);

  dims[0] = numChans;
  typenum = NPY_DOUBLE;

  /* create the output arrays from the data returned in groupCol and qualCol */
  grouping = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  quality = (PyArrayObject *)PyArray_SimpleNew(1, dims, typenum);
  if ((NULL == grouping) || (NULL == quality))
  {
    sprintf(groupmsg, GROUP_CREATE_MSG, funcName);
    PyErr_SetString(PyExc_Exception, groupmsg);
    return NULL;
  }/*end if... */
  groupData = DDATA(grouping);
  qualData = DDATA(quality);
  for (ii = 0; ii < numChans; ii++)
  {
    groupData[ii] = groupCol[ii]; /*grab the data from groupCol and place in grouping data */
    qualData[ii] = qualCol[ii]; /*grab the data from qualCol and place in quality data */
  }/*end for... */

  free(groupCol); /* free the allocated memory */
  free(qualCol);
  free(c_countsArray);
  free(c_errorCol);
  free(c_tabStops);

  /* Return grouping and quality NumPy arrays */
  return Py_BuildValue("OO", PyArray_Return(grouping), PyArray_Return(quality));
}/*end...grpSnr*/
