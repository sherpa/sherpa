/*_C_INSERT_SAO_COPYRIGHT_HERE_(2002-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/


#include "grp_priv.h"
#include "grplib.h"


/*********************************************************************/
/************************  Utility functions  ************************/
/*********************************************************************/

/* Counts the number of tabstopped (non-zero) channels in array */
long count_bad_chans(short *tabStops, long numChans){
    
    long ii, counter = 0;

    for(ii = 0; ii < numChans; ii++){
        if(tabStops[ii]) counter++;
    }

    return(counter);
}

/* Rounds to the nearest integer */
long my_round(double num){
    return((long) (num + 0.5));
}

/* Counts the number of groups  - do bad groups "count"? */
long count_groups(short *groupCol, short *qualCol, long numChans){

    long ii, counter = 0;

    for(ii = 0; ii < numChans; ii++){
        if(groupCol[ii] == GRP_BEGIN) counter++;
    }

    return(counter);
}

/* Creates a "bad" group with quality = 2 */
void set_incomplete(short *groupCol, short *qualCol, long startChan,
                    long endChan){

    create_group(groupCol, startChan, endChan);
    set_quality(qualCol, (short) GRP_POOR, startChan, endChan);
}

/* Sets startChan's grouping to 1, and the rest to -1 */
void create_group(short *groupCol, long startChan, long endChan){

    long ii;

    for(ii = startChan; ii <= endChan; ii++){
        if(ii == startChan)
            groupCol[ii] = GRP_BEGIN;
        else
            groupCol[ii] = GRP_MIDDLE;
    }
}

/* Sets startChan-endChan inclusive with quality = qualVal */
void set_quality(short *qualCol, short qualVal, long startChan,
                 long endChan){

    long ii;

    for(ii = startChan; ii <= endChan; ii++){
        qualCol[ii] = qualVal;
    }
}

/* Marks channels in specified interval as "used" */
void mark_used(short *usedChans, long startChan, long endChan){

    long ii;

    for(ii = startChan; ii <= endChan; ii++){
        usedChans[ii] = (short) GRP_TRUE;
    }
}

/* Return the first channel with a data value greater than or equal
 * to the given value.
 */
long lower_bound(double value, double *dataCol, long numChans,
                 int isAscending, dsErrList *errList){

    long ii;

    if (isAscending) {
      for(ii = 0; ii < numChans; ii++){
        if(dataCol[ii] >= value) return(ii);
      }
    } else {
      for(ii = numChans-1; ii >= 0; ii--){
        if(dataCol[ii] >= value) 
	  return(ii);

      }
    }

    if(errList)
        dsErrAdd(errList, dsDMGROUPLOWERBOUNDERR, Individual,
                 Generic);
    else
        err_msg("ERROR: grp_priv.c:lower_bound(): No data greater "
                "than or equal to given value.\n");
    return(GRP_ERROR);
}

/* Return the last channel with a data value less than or equal
 * to the given value.
 */
long upper_bound(double value, double *dataCol, long numChans,
                 int isAscending, int isColReal, dsErrList *errList){

    long ii;

    /* if the data column is real (not int), the last bin contains
       data up to, but not including the upper limit.  */
    if (isColReal) {

      if (isAscending) {
	for(ii = (numChans - 1); ii >= 0; ii--){
	  if(dataCol[ii] < value) return(ii);
	}
      } else {
	for(ii = 0; ii < numChans; ii++){
	  if(dataCol[ii] < value) 
	    return(ii);
	}
      }

    /* if the data column is not real (i.e., int), the
       last bin contains data at the upper limit.  */
    } else {

      if (isAscending) {
	for(ii = (numChans - 1); ii >= 0; ii--){
	  if(dataCol[ii] <= value) return(ii);
	}
      } else {
	for(ii = 0; ii < numChans; ii++){
	  if(dataCol[ii] <= value) 
	    return(ii);
	}
      }

    }

    if(errList)
        dsErrAdd(errList, dsDMGROUPUPPERBOUNDERR, Individual,
                 Generic);
    else
      (isColReal) ?
	err_msg("ERROR: grp_priv.c:upper_bound(): No data less "
                "than given value.\n")
	  :
	err_msg("ERROR: grp_priv.c:upper_bound(): No data less "
		"than or equal to given value.\n");


    return(GRP_ERROR);
}

/* Check to see if any bins overlap */
int check_overlap(double *binLow, double *binHigh, long numBins){

    long ii, jj;
    
    for(ii = 0; ii < (numBins - 1); ii++){
        for(jj = (ii + 1); jj < numBins; jj++){

	  if((binLow[ii] < binHigh[jj]) && (binLow[ii] > binLow[jj]) &&
	     (fabs(binLow[ii]-binHigh[jj]) > binLow[ii]*FLT_EPSILON))
	    return(GRP_ERROR);	    
	  
	  if ((binHigh[ii] > binLow[jj]) && (binHigh[ii] < binHigh[jj]) &&
	     (fabs(binLow[jj]-binHigh[ii]) > binLow[jj]*FLT_EPSILON))
	    return(GRP_ERROR);
	    

        }
    }

    return(GRP_SUCCESS);
}

/* Check to make sure data is monotonically increasing */
int check_increasing(double *dataCol, long numChans){

    long ii;
    double tempNum;

    tempNum = dataCol[0];

    if(numChans == 1) return(GRP_SUCCESS);
    
    for(ii = 1; ii < numChans; ii++){
      if(dataCol[ii] <= tempNum)
	return(GRP_ERROR);
      else tempNum = dataCol[ii];
    }

    return(GRP_SUCCESS);
}

    
/* Check to make sure data is monotonically decreasing */
int check_decreasing(double *dataCol, long numChans){

    long ii;
    double tempNum;

    tempNum = dataCol[0];

    if(numChans == 1) return(GRP_SUCCESS);
    
    for(ii = 1; ii < numChans; ii++){
        if(dataCol[ii] > tempNum)
            return(GRP_ERROR);
        else tempNum = dataCol[ii];
    }

    return(GRP_SUCCESS);
}

