/*                                                                
**  Copyright (C) 2002-2007  Smithsonian Astrophysical Observatory 
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


 
#include "grplib.h"
#include "grp_priv.h"


/*********************************************************************/
/************************  Grouping Routines  ************************/
/*********************************************************************/

/* NUM_BINS
 * --------
 * Input parameters:
 *   numChans - number of channels in groupCol and qualCol
 *   numBins  - number of equal-sized bins desired
 *   groupCol - the GROUPING column
 *   qualCol  - the QUALITY column
 *   tabStops - array giving channels with tabs or stops
 */
int grp_do_num_bins(long numChans, long numBins, short *groupCol,
                    short *qualCol, short *tabStops,
                    dsErrList *errList){

    long numTabStops, binWidth, numGroups;
    short isError = GRP_FALSE;
    int returnVal = GRP_SUCCESS;

    /* Check for obviously bad inputs */
    if((numChans <= 0) || (numBins <= 0) || !groupCol || !qualCol
       || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    numTabStops = count_bad_chans(tabStops, numChans);
    binWidth = (double) ((double) numChans -
			 (double) numTabStops) / (double) numBins;

    if(!binWidth){
        if(errList)
            dsErrAdd(errList, dsDMGROUPZEROWIDTHERR, Accumulation,
                     Generic);
        else
            err_msg("WARNING: The calculated bin width rounds to "
                    "zero.\nIt will be reset to 1.");
        binWidth = 1;
        returnVal = GRP_WARNING;
    }
    
    isError = grp_do_bin_width(numChans, binWidth, groupCol,
                               qualCol, tabStops, errList);
    if(isError){
        if(isError == GRP_ERROR)
            return(GRP_ERROR);
        else returnVal = GRP_WARNING;
    }
    
    numGroups = count_groups(groupCol, qualCol, numChans);
    /* More groups produced than requested? */
    if(numGroups > numBins){
        if(errList)
            dsErrAdd(errList, dsDMGROUPEXTRAGROUPSERR, Accumulation,
                     Generic);
        else
            err_msg("WARNING: More groups produced than requested.\n");
        returnVal = GRP_WARNING;
    }
    /* Fewer groups produced than requested? */
    else if(numGroups < numBins){
        if(errList)
            dsErrAdd(errList, dsDMGROUPTOOFEWGROUPSERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: Fewer groups produced than requested.\n");
        return(GRP_ERROR);
    }
    
    return(returnVal);
}

/* NUM_CTS
 * -------
 * Input parameters:
 *   dataCol   - the column with the counts
 *   numChans  - number of channels in groupCol and qualCol
 *   numCounts - the minimum number of counts desired in each group
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   maxLength - maximum size of groups
 */
int grp_do_num_counts(double *dataCol, long numChans, double numCounts,
                      short *groupCol, short *qualCol, short *tabStops,
                      double maxLength, dsErrList *errList){

   long ii, counter = 0;
   double totCounts = 0.0;

   /* Check for obviously bad inputs */
   if(!dataCol || (numChans <= 0) || (numCounts <= 0) || !groupCol
      || !qualCol || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
   }
   
   if(maxLength <= 0.0)
       maxLength = MAX_BIN_LENGTH;

   for(ii = 0; ii < numChans; ii++){
       /* Are we in a tab or stop? */
       if(tabStops[ii]){
           if(counter != 0){
               set_incomplete(groupCol, qualCol, ii - counter, ii - 1);
               totCounts = 0.0;
               counter = 0;
           }
       }
       /* Are we at the end of the table? */
       else if(ii == (numChans - 1)){
           /* Does this complete a group? */
           if((totCounts + dataCol[ii] >= numCounts) ||
              ((counter + 1) >= maxLength)){
               groupCol[ii] = counter ? GRP_MIDDLE : GRP_BEGIN;
               qualCol[ii] = GRP_GOOD;
           }
           else{
               set_incomplete(groupCol, qualCol, ii - counter, ii);
           }
       }
       /* Are we at the end of a group or have reached maxLength? */
       else if((totCounts + dataCol[ii] >= numCounts) ||
               ((counter + 1) >= maxLength)){
           groupCol[ii] = counter ? GRP_MIDDLE : GRP_BEGIN;
           qualCol[ii] = GRP_GOOD;
           counter = 0;
           totCounts = 0.0;
       }
       /* Are we at the beginning of a group? */
       else if(counter == 0){
           groupCol[ii] = GRP_BEGIN;
           qualCol[ii] = GRP_GOOD;
           counter++;
           totCounts += dataCol[ii];
       }
       /* We must be in the middle of a group */
       else{
           groupCol[ii] = GRP_MIDDLE;
           qualCol[ii] = GRP_GOOD;
           counter++;
           totCounts += dataCol[ii];
       }
   } /* end for(ii) */
   
   return(GRP_SUCCESS);
}

/* BIN
 * ---
 * Input parameters:
 *   dataCol   - the column with the channel data
 *   numChans  - number of channels in groupCol and qualCol
 *   binLow    - array of lower group boundaries
 *   binHigh   - array of higher group boundaries
 *   numBins   - number of binLow-binHigh group pairs
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 */
int grp_do_bin(double *dataCol, long numChans, double *binLow,
               double *binHigh, long numBins, short *groupCol,
               short *qualCol, short *tabStops, dsErrList *errList,
               short partialBin, int isColReal){

    long ii, jj, binLength, counter, tempLow, tempHigh, chanHigh;
    short isComplete;
    double maxVal, dataHigh;
    int tmpVar = 0;
    int isAscending = 0;

    /* Check for obviously bad inputs */
    if(!dataCol || (numChans <= 0) || !binLow || !binHigh
       || (numBins < 0) || !groupCol || !qualCol || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    

    /* Check for monotonically increasing or decreasing  data */
    if (check_increasing(dataCol, numChans) == 0) {
      isAscending = 1;
    } else { 
      if ( check_decreasing(dataCol, numChans) != 0 ) {
	if(errList)
	  dsErrAdd(errList, dsDMGROUPBADDATAORDERERR, Accumulation,
		   Generic);
        else
	  err_msg("ERROR: Data column is not increasing/decreasing.\n");
        return(GRP_ERROR);
      }
    }

    if (isAscending) {
      dataHigh = dataCol[numChans - 1];
      maxVal = binHigh[numBins - 1];
    } else {
      dataHigh = dataCol[0];
      maxVal = binLow[0];
    }

    chanHigh = upper_bound(maxVal, dataCol, numChans, isAscending, isColReal, errList);

    /* Check if any bins overlap */
    if(check_overlap(binLow, binHigh, numBins)){
        if(errList)
            dsErrAdd(errList, dsDMGROUPOVERLAPBINSPECERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least two bins in binspec overlap.\n");
        return(GRP_ERROR);
    }
    
    /* Go through all binLow-binHigh pairs */
    for(ii = 0; ii < numBins; ii++){

        tempLow = lower_bound(binLow[ii], dataCol, numChans, isAscending, 
			      errList);
        tempHigh = upper_bound(binHigh[ii], dataCol, numChans, isAscending, 
			       isColReal, errList);

	if (!isAscending) {
	  tmpVar = tempLow;
	  tempLow = tempHigh;
	  tempHigh = tmpVar;
	}

        if((tempLow == GRP_ERROR) || (tempHigh == GRP_ERROR))
	  continue;
	/*            return(GRP_ERROR); */
        
        /* If there are no data within the pair, skip to next pair */
        if((tempHigh - tempLow) < 0)
            continue;

        binLength = tempHigh - tempLow;

        if(binLow[ii] > dataHigh){
            if(errList)
                dsErrAdd(errList, dsDMGROUPINVALIDBINERR, Accumulation,
                         Generic);
            else
                err_msg("ERROR: A bin boundary is invalid.\nMake "
                        "sure the binspec fits within the bounds "
                        "of the data.\n");
            return(GRP_ERROR);
        }

        isComplete = GRP_TRUE;
        
        /* Check for a complete group */
        for(jj = tempLow; jj <= tempHigh; jj++){
            if ( (jj > chanHigh) || tabStops[jj] ) {
                isComplete = GRP_FALSE;
                break;
            }
        }

        /* If it's the last group, and partialBin is TRUE, then it's
           not complete. */
        if( (partialBin && (ii == (numBins - 1))) ||
	    ((!isAscending) && (binLength < (numBins-1))) )
            isComplete = GRP_FALSE;
        
        counter = 0;
        
        /* Create group - good or bad */
        for(jj = tempLow; jj <= tempHigh; jj++){
            /* Are we in a tab or stop? */
            if(tabStops[jj]){
                counter = 0;
            }
            /* Are we at the end of the table? */
            else if(jj == (numChans - 1)){
                /* Is this a single-element group? */
                if(counter == 0){
                    if(isComplete){
                        groupCol[jj] = GRP_BEGIN;
                        qualCol[jj] = GRP_GOOD;
                        break;
                    }
                    else{
                        groupCol[jj] = GRP_BEGIN;
                        qualCol[jj] = GRP_POOR;
                        break;
                    }
                }
                /* Does this complete the group? */
                if(isComplete){
                    groupCol[jj] = GRP_MIDDLE;
                    qualCol[jj] = GRP_GOOD;
                    break;
                }
                else{
                    groupCol[jj] = GRP_MIDDLE;
                    qualCol[jj] = GRP_POOR;
                    break;
                }
            }
            /* Are we at the beginning of a group? */
            else if(counter == 0){
                if(isComplete){
                    groupCol[jj] = GRP_BEGIN;
                    qualCol[jj] = GRP_GOOD;
                    counter++;
                }
                else{
                    groupCol[jj] = GRP_BEGIN;
                    qualCol[jj] = GRP_POOR;
                    counter++;
                }
            }
            /* We must be in the middle of a group */
            else{
                if(isComplete){
                    groupCol[jj] = GRP_MIDDLE;
                    qualCol[jj] = GRP_GOOD;
                    counter++;
                }
                else{
                    groupCol[jj] = GRP_MIDDLE;
                    qualCol[jj] = GRP_POOR;
                    counter++;
                }
            }
        } /* end for(jj) */
    } /* end for(ii) */
    
    return(GRP_SUCCESS);
}

/* BIN_WIDTH
 * ---------
 * Input parameters:
 *   numChans - number of channels in groupCol and qualCol
 *   binWidth - size of bins desired
 *   groupCol - the GROUPING column
 *   qualCol  - the QUALITY column
 *   tabStops - array giving channels with tabs or stops
 */
int grp_do_bin_width(long numChans, long binWidth, short *groupCol,
                     short *qualCol, short *tabStops,
                     dsErrList *errList){

    long ii, counter = 0;

    /* Check for obviously bad inputs */
    if((numChans <= 0) || (binWidth <= 0) || !groupCol || !qualCol
       || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    for(ii = 0; ii < numChans; ii++){
        /* Are we in a tab or stop? */
        if(tabStops[ii]){
            if(counter != 0){
                set_incomplete(groupCol, qualCol, ii - counter,
                               ii - 1);
                counter = 0;
            }
        }
        /* Are we at the end of the table? */
        else if(ii == (numChans - 1)){
            /* Does this complete a group? */
            if((counter + 1) == binWidth){
                groupCol[ii] = counter ? GRP_MIDDLE : GRP_BEGIN;
                qualCol[ii] = GRP_GOOD;
            }
            else{
                set_incomplete(groupCol, qualCol, ii - counter, ii);
            }
        }
        /* Are we at the end of a group? */
        else if((counter + 1) == binWidth){
            groupCol[ii] = counter ? GRP_MIDDLE : GRP_BEGIN;
            qualCol[ii] = GRP_GOOD;
            counter = 0;
        }
        /* Are we at the beginning of a group? */
        else if(counter == 0){
            groupCol[ii] = 1;
            qualCol[ii] = GRP_GOOD;
            counter++;
        }
        /* We must be in the middle of a group */
        else{
            groupCol[ii] = GRP_MIDDLE;
            qualCol[ii] = GRP_GOOD;
            counter++;
        }
    } /* end for(ii) */
    
    return(GRP_SUCCESS);   
}

/* SNR
 * ---
 * Input parameters:
 *   dataCol   - the column with the data
 *   numChans  - number of channels in groupCol and qualCol
 *   snr       - the signal-to-noise ratio threshold
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   errorCol  - optional error column
 *   useErr    - if true, use errorCol data, else use counts
 *   maxLength - maximum size of groups
 */
int grp_do_snr(double *dataCol, long numChans, double snr,
               short *groupCol, short *qualCol, short *tabStops,
               double *errorCol, short useErr, double maxLength,
               dsErrList *errList){

    double runningSnr = 0.0;
    double runningSignal=0.0;
    double runningNoise=0.0;
    long ii, counter = 0;
    int returnVal = GRP_SUCCESS;

    /* Check for obviously bad inputs */
    if(!dataCol || (numChans <= 0) || (snr <= 0) || !groupCol
       || !qualCol || !tabStops || !errorCol){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    if(maxLength <= 0.0)
        maxLength = MAX_BIN_LENGTH;
   
    for(ii = 0; ii < numChans; ii++){

        if(tabStops[ii]){
            if(counter != 0){
                set_incomplete(groupCol, qualCol, ii - counter,
                               ii - 1);
                runningSnr = 0.0;
		runningSignal=0.0;
		runningNoise=0.0;
                counter = 0;
            }
        }
        /* Are we at the end of the table? */
        else if(ii == (numChans - 1)){
            /* Does this complete a group? */
            if(useErr){
                if(!errorCol[ii]){
                    if(errList)
                        dsErrAdd(errList, dsDMGROUPZEROERRORERR,
                                 Accumulation, Generic);
                    else
                        err_msg("WARNING: The supplied error column "
                                "contains zero-valued data.");
                    returnVal = GRP_WARNING;
                }
                else {
		  runningSignal += dataCol[ii];
		  runningNoise += errorCol[ii]*errorCol[ii];
		  runningSnr = (runningSignal / sqrt(runningNoise) );
		  runningSnr *= runningSnr;
		}
            }
            else if(dataCol[ii]) {
	      runningSignal += dataCol[ii];
	      runningSnr = runningSignal;
	    }

            if((sqrt(runningSnr) > snr) ||
               ((counter + 1) >= maxLength)){
                groupCol[ii] = counter ? GRP_MIDDLE : GRP_BEGIN;
                qualCol[ii] = GRP_GOOD;
		counter++;
            }
            else {
                set_incomplete(groupCol, qualCol, ii - counter, ii);
	    }
        }
        /* Have we reached maxLength? */
        else if((counter + 1) >= maxLength){
            groupCol[ii] = GRP_MIDDLE;
            qualCol[ii] = GRP_GOOD;
            counter = 0;
            runningSnr = 0.0;
	    runningSignal=0.0;
	    runningNoise=0.0;
        }
        /* Are we at the end of a group? */
        else if(sqrt(runningSnr) > snr){
            groupCol[ii] = GRP_MIDDLE;
            qualCol[ii] = GRP_GOOD;
            counter=0;
            runningSnr = 0.0;
	    runningSignal=0.0;
	    runningNoise=0.0;
        }

        /* Are we at the beginning of a group? */
	if(counter == 0){
            if(useErr){
                if(!errorCol[ii]){
                    if(errList)
                        dsErrAdd(errList, dsDMGROUPZEROERRORERR,
                                 Accumulation, Generic);
                    else
                        err_msg("WARNING: The supplied error column "
                                "contains zero-valued data.");
                    returnVal = GRP_WARNING;
                }
                else {
		  runningSignal += dataCol[ii];
		  runningNoise += errorCol[ii]*errorCol[ii];
		  runningSnr = (runningSignal/sqrt(runningNoise));
		  runningSnr *= runningSnr;
		}
            } /* no user error */
            else if(dataCol[ii]) {
	      runningSignal += dataCol[ii];
	      runningSnr = runningSignal;
	    }
            groupCol[ii] = GRP_BEGIN;
            qualCol[ii] = GRP_GOOD;
            counter++;
        } /* end else counter == 0 */
        /* We must be in the middle of a group */
        else if (ii != (numChans - 1)) {
            if(useErr){
                if(!errorCol[ii]){
                    if(errList)
                        dsErrAdd(errList, dsDMGROUPZEROERRORERR,
                                 Accumulation, Generic);
                    else
                        err_msg("WARNING: The supplied error column "
                                "contains zero-valued data.");
                    returnVal = GRP_WARNING;
                }
                else {
		  runningSignal += dataCol[ii];
		  runningNoise += errorCol[ii]*errorCol[ii];
		  runningSnr = (runningSignal/sqrt(runningNoise));
		  runningSnr *= runningSnr;
		}
	    }
            else if(dataCol[ii]) {
	      runningSignal += dataCol[ii];
	      runningSnr = runningSignal;
	    }
            groupCol[ii] = GRP_MIDDLE;
	    qualCol[ii] = GRP_GOOD;
            counter++;
        }
    } /* end for(ii) */
    
    return(returnVal);
}

/* ADAPTIVE
 * --------
 * Input parameters:
 *   dataCol   - the column with the data
 *   numChans  - number of channels in groupCol and qualCol
 *   minCounts - the minimum number of counts desired in each group
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   maxLength - maximum size of groups
 */
int grp_do_adaptive(double *dataCol, long numChans, double minCounts,
                    short *groupCol, short *qualCol, short *tabStops,
                    double maxLength, dsErrList *errList){
    
    short *usedChans;
    long ii, jj, tempLength, tempMax, curWidth = 0;
    long counter = 0;
    double groupCounts = 0.0;

    /* Check for obviously bad inputs */
    if(!dataCol || (numChans <= 0) || (minCounts <= 0) || !groupCol
       || !qualCol || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    if(maxLength <= 0.0)
        maxLength = MAX_BIN_LENGTH;
   
    /* Create and initialize used channel list */
    usedChans = (short *) calloc(numChans, sizeof(short));
    for(ii = 0; ii < numChans; ii++){
        if(tabStops[ii] || (qualCol[ii] != 0))
            usedChans[ii] = GRP_TRUE;
        else
            usedChans[ii] = GRP_FALSE;
    }
    
    /* Main loop through adaptive group sizes */
    while((curWidth + 1) <= maxLength){
        curWidth++;
        
        /* Determine maxLength each time as it might be limited */
        tempLength = 0;
        tempMax = 0;
        for(ii = 0; ii < numChans; ii++){
            if(!usedChans[ii]){
                tempLength++;
                if(tempLength > tempMax)
                    tempMax = tempLength;
            }
            else
                tempLength = 0;
        }
        if(tempMax < maxLength)
            maxLength = tempMax;
        
        /* Iterate over each row for each group size */
        for(ii = 0; ii < (numChans - curWidth); ii++){
            if(usedChans[ii]) continue;
            groupCounts = 0.0;
            /* Try to make groups of the current width */
            for(jj = 0; jj < curWidth; jj++){
                if(usedChans[ii + jj])
                    break;
                groupCounts += dataCol[ii + jj];
                if(jj == curWidth - 1){
                    if(groupCounts >= minCounts){
                        /* Enough counts - let's group it */
                        mark_used(usedChans, ii, ii + jj);
                        create_group(groupCol, ii, ii + jj);
                        set_quality(qualCol, GRP_GOOD, ii, ii + jj);
                    }
                }
            } /* end for(jj) */
        } /* end for(ii) */
    } /* end while() */
    
    /* Put unused channels into "bad" groups */
    for(ii = 0; ii < numChans; ii++){
        /* Are we in a used channel? */
        if(usedChans[ii]){
            if(counter != 0){
                set_incomplete(groupCol, qualCol, ii - counter,
                               ii - 1);
                counter = 0;
            }
        }
        /* Are we at the end of the table? */
        else if(ii == (numChans - 1)){
            /* Does this complete a group? */
            if(counter != 0)
                set_incomplete(groupCol, qualCol, ii - counter, ii);
            else
                set_incomplete(groupCol, qualCol, ii, ii);
        }
        /* Are we at the end of a group */
        else if(usedChans[ii + 1]){
            set_incomplete(groupCol, qualCol, ii - counter, ii);
            counter = 0;
        }
        /* Are we at the beginning of a group? */
        else{
            counter++;
        }
    } /* end for(ii) */
    
    free(usedChans);
    return(GRP_SUCCESS);
}

/* ADAPTIVE_SNR
 * ------------
 * Input parameters:
 *   dataCol   - the column with the data
 *   numChans  - number of channels in groupCol and qualCol
 *   snr       - the signal-to-noise ratio threshold
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   errorCol  - optional error column
 *   useErr    - if true, use errorCol data, else use counts
 *   maxLength - maximum size of groups
 */
int grp_do_adaptive_snr(double *dataCol, long numChans, double snr,
                        short *groupCol, short *qualCol,
                        short *tabStops, double *errorCol,
                        short useErr, double maxLength,
                        dsErrList *errList){

   long ii, jj, tempLength, tempMax, counter = 0;
   long curWidth = 0;
   double runningSnr    = 0.0;
   double runningSignal = 0.0;
   double runningNoise  = 0.0;
   short *usedChans;
   int returnVal = GRP_SUCCESS;

    /* Check for obviously bad inputs */
    if(!dataCol || (numChans <= 0) || (snr <= 0) || !groupCol
       || !qualCol || !tabStops || !errorCol){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    if(maxLength <= 0.0)
       maxLength = MAX_BIN_LENGTH;
   
   /* Create used channel list */
   usedChans = (short *) calloc(numChans, sizeof(short));
   for(ii = 0; ii < numChans; ii++){
      if(tabStops[ii] || (qualCol[ii] != 0))
         usedChans[ii] = GRP_TRUE;
      else
         usedChans[ii] = GRP_FALSE;
   }
   /* Main loop through adaptive group sizes */
   while((curWidth + 1) <= maxLength){
      curWidth++;

      /* Determine maxLength each time as it might be limited */
      tempLength = 0;
      tempMax = 0;
      for(ii = 0; ii < numChans; ii++){
         if(!usedChans[ii]){
            tempLength++;
            if(tempLength > tempMax)
               tempMax = tempLength;
         }
         else
            tempLength = 0;
      }
      if(tempMax < maxLength)
         maxLength = tempMax;

      /* Iterate over each row for each group size */
      for(ii = 0; ii < (numChans - curWidth); ii++){
         if(usedChans[ii]) continue;
         runningSnr = 0.0;
	 runningSignal = 0.0;
	 runningNoise = 0.0;
         /* Try to make groups of the current width */
         for(jj = 0; jj < curWidth; jj++){
            if(usedChans[ii + jj]){
               runningSnr = 0.0;
	       runningSignal = 0.0;
	       runningNoise = 0.0;
               break;
            }
            if(useErr){
                if(!errorCol[ii + jj]){
                    if(errList)
                        dsErrAdd(errList, dsDMGROUPZEROERRORERR,
                                 Accumulation, Generic);
                    else
                        err_msg("WARNING: The supplied error column "
                                "contains zero-valued data.");
                    returnVal = GRP_WARNING;
                }
                else {
                    runningSnr += pow((dataCol[ii + jj] /
                                       errorCol[ii + jj]), 2);
		    runningSignal += dataCol[ii+jj];
		    runningNoise += (errorCol[ii+jj]*errorCol[ii+jj]);
		    runningSnr = runningSignal/runningNoise;
		    runningSnr *= runningSnr;

		}
            }
            else if(dataCol[ii + jj]) {
	      runningSignal += dataCol[ii];
	      runningSnr = runningSignal;
	    }
            if(jj == (curWidth - 1)){
               if(sqrt(runningSnr) > snr){
                  /* Enough counts - let's group it */
                  mark_used(usedChans, ii, ii + jj);
                  create_group(groupCol, ii, ii + jj);
                  set_quality(qualCol, GRP_GOOD, ii, ii + jj);
               }
            }
         } /* end for(jj) */
      } /* end for(ii) */
   } /* end while() */

   /* Put unused channels into "bad" groups */
   for(ii = 0; ii < numChans; ii++){
      /* Are we in a used channel? */
      if(usedChans[ii]){
         if(counter != 0){
            set_incomplete(groupCol, qualCol, ii - counter, ii - 1);
            counter = 0;
         }
      }
      /* Are we at the end of the table? */
      else if(ii == (numChans - 1)){
         /* Does this complete a group? */
         if(counter != 0)
            set_incomplete(groupCol, qualCol, ii - counter, ii);
         else
            set_incomplete(groupCol, qualCol, ii, ii);
      }
      /* Are we at the end of a group */
      else if(usedChans[ii + 1]){
         set_incomplete(groupCol, qualCol, ii - counter, ii);
         counter = 0;
      }
      /* Are we at the beginning of a group? */
      else{
         counter++;
      }
   } /* end for(ii) */

   free(usedChans);
   return(GRP_SUCCESS);
}

/* BIN_FILE
 * --------
 * Input parameters:
 *   dataCol   - the column with the channel data
 *   numChans  - number of channels in groupCol and qualCol
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   fDataCol  - the data column from template file
 *   fNumChans - number of channels from template file
 *   fGroupCol - the grouping column from template file
 *   fQualCol  - the quality column from template file
 */
int grp_do_bin_file(double *dataCol, long numChans, short *groupCol,
                    short *qualCol, short *tabStops, double *fDataCol,
                    long fNumChans, short *fGroupCol, short *fQualCol, 
		    int isColReal, dsErrList *errList){

    short isError = GRP_FALSE;
    double *binLow, *binHigh, *tBinLow, *tBinHigh;
    long groupBegin = -1, tabBegin = -1;
    long ii, index = 0, tIndex = 0;

    /* Check for obviously bad inputs */
    if(!dataCol || (numChans <= 0) || !groupCol || !qualCol
       || !tabStops || !fDataCol || (fNumChans <= 0) || !fGroupCol
       || !fQualCol){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
      
    binLow = (double *) calloc(fNumChans, sizeof(double));
    binHigh = (double *) calloc(fNumChans, sizeof(double));
    tBinLow = (double *) calloc(fNumChans, sizeof(double));
    tBinHigh = (double *) calloc(fNumChans, sizeof(double));

    /* Fill in tabStops */
    for(ii = 0; ii < fNumChans; ii++){
        if((fQualCol[ii] == GRP_TABBED)
           || (fGroupCol[ii] == GRP_UNUSED))
            tabStops[ii] = GRP_TRUE;
        else
            tabStops[ii] = GRP_FALSE;
    }
    
    /* Need to determine binLow, binHigh, numBins, tabStops */
    for(ii = 0; ii < (fNumChans - 1); ii++){
        
        /* Are we starting a group? */
        if((groupBegin < 0) && (fGroupCol[ii] == GRP_BEGIN))
            groupBegin = ii;

        /* Are we starting a tab? */
        if((tabBegin < 0) && (fQualCol[ii] == GRP_TABBED))
            tabBegin = ii;
        
        /* Are we at the end of a group? */
        if((groupBegin >= 0) && ((fGroupCol[ii + 1] == GRP_BEGIN) ||
                          (fQualCol[ii + 1] == GRP_TABBED))){
            /* Make group */
            binLow[index] = fDataCol[groupBegin];
            binHigh[index] = fDataCol[ii + 1];
            groupBegin = -1;
            index++;
        }
        
        /* Are we at the end of a tab? */
        if((tabBegin >= 0) && (fQualCol[ii + 1] != GRP_TABBED)){
            /* Make tab */
            tBinLow[tIndex] = fDataCol[tabBegin];
            tBinHigh[tIndex] = fDataCol[ii + 1];
            tabBegin = -1;
            tIndex++;
        }
        
        /* Are we at the end of the table? Handle grouping. */
        if((ii == (fNumChans - 2)) && (fQualCol[ii] != GRP_TABBED)){
            if((groupBegin < 0) && (fGroupCol[ii + 1] == GRP_BEGIN)){
                binLow[index] = fDataCol[ii + 1];
                binHigh[index] = fDataCol[ii + 1];
                index++;
            }
            else if((groupBegin >= 0) &&
                    (fGroupCol[ii + 1] != GRP_BEGIN)){
                binLow[index] = fDataCol[groupBegin];
                binHigh[index] = fDataCol[ii + 1];
                index++;
            }
        }
        
        /* Are we at the end of the table? Handle tabs. */
        if(ii == (fNumChans - 2)){
            if((tabBegin < 0) && (fQualCol[ii + 1] == GRP_TABBED)){
                tBinLow[tIndex] = fDataCol[ii + 1];
                tBinHigh[tIndex] = fDataCol[ii + 1];
                tIndex++;
            }
            else if(tabBegin >= 0){
                tBinLow[tIndex] = fDataCol[tabBegin];
                tBinHigh[tIndex] = fDataCol[ii + 1];
                tIndex++;
            }
        }
    } /* end for(ii) */
   
    isError = grp_do_bin(dataCol, numChans, binLow, binHigh, 
			 index, groupCol, qualCol, tabStops, 
			 errList, 0, isColReal);
    if(isError)
        return(GRP_ERROR);

    isError = set_tabs(dataCol, groupCol, qualCol, numChans, tBinLow,
                       tBinHigh, tIndex, 1, isColReal, errList);
    if(isError){
        err_msg("Error setting tabs in method BIN_FILE\n");
        return(GRP_ERROR);
    }

    free(binLow);
    free(binHigh);
    free(tBinLow);
    free(tBinHigh);
    
    return(GRP_SUCCESS);
}

/* NONE
 * ----
 * Input parameters:
 *   numChans - number of channels in groupCol and qualCol
 *   groupCol - the GROUPING column
 *   qualCol  - the QUALITY column
 */
int grp_do_none(long numChans, short *groupCol, short *qualCol,
                dsErrList *errList){

    long ii;

    if((numChans <= 0) || !groupCol || !qualCol){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    /* Reset all columns */
    for(ii = 0; ii < numChans; ii++){
        groupCol[ii] = GRP_BEGIN;
        qualCol[ii] = GRP_GOOD;
    }
    
    return(GRP_SUCCESS);
}

/* MIN_SLOPE
 * ---------
 * Input parameters:
 *   dataCol   - the column representing the independent axis
 *   binCol    - the column with the data values
 *   numChans  - number of channels in groupCol and qualCol
 *   slope     - the minimum slope threshold
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   maxlength - maximum size of groups
 */
int grp_do_min_slope(double *dataCol, double *binCol, long numChans,
                     double slope, short *groupCol, short *qualCol,
                     short *tabStops, double maxlength,
                     dsErrList *errList){

    long ii, jj, counter = 0;
    double range = 0.0;
    double tempSlope = 0.0;
    short *usedChans;
    
    /* Check for obviously bad inputs */
    if(!dataCol || !binCol || (numChans < 2) || (slope <= 0)
       || !groupCol || !qualCol || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    if(maxlength <= 0.0)
        maxlength = MAX_BIN_LENGTH;
    
    /* Create and initialize used channel list */
    usedChans = (short *) calloc(numChans, sizeof(short));
    for(ii = 0; ii < numChans; ii++){
        if(tabStops[ii] || (qualCol[ii] != 0))
            usedChans[ii] = GRP_TRUE;
        else
            usedChans[ii] = GRP_FALSE;
    }
    
    ii = 0;
    jj = 1;
    while(ii < (numChans - 1)){
        
        /* Are we in a tab or stop? */
        if(tabStops[ii]){
            ii++;
            jj = (ii + 1);
        }
        else{
            while(jj < numChans){
                
                /* Calculate current slope */
                tempSlope = fabs(((binCol[jj] - binCol[ii])
                                  / (dataCol[jj] - dataCol[ii])));
                
                range = (dataCol[jj] - dataCol[ii]);
                
                /* Are we in a tab or stop? */
                if(tabStops[jj]){
                    ii++;
                    jj = (ii + 1);
                    break;
                }
                /* Are we at the end of the table? */
                else if(jj == (numChans - 1)){
                    /* Does this complete a group? */
                    if((tempSlope <= slope) ||
                       (range >= maxlength)){
                        mark_used(usedChans, ii, jj);
                        create_group(groupCol, ii, jj);
                        set_quality(qualCol, GRP_GOOD, ii, jj);
                        ii = jj;
                        break;
                    }
                    else{
                        ii++;
                        jj = (ii + 1);
                        break;
                    }
                }
                /* Are we at the end of a group or have reached *
                 * maxlength? */
                else if((tempSlope <= slope) ||
                        (range >= maxlength)){
                    mark_used(usedChans, ii, jj);
                    create_group(groupCol, ii, jj);
                    set_quality(qualCol, GRP_GOOD, ii, jj);
                    ii = (jj + 1);
                    jj = (ii + 1);
                    break;
                }
                /* Keep looking */
                else jj++;
                
            } /* end while(jj) */
        } /* end if */
    } /* end while(ii) */
    
    /* Put unused channels into "bad" groups */
    for(ii = 0; ii < numChans; ii++){
        /* Are we in a used channel? */
        if(usedChans[ii]){
            if(counter != 0){
                set_incomplete(groupCol, qualCol, ii - counter,
                               ii - 1);
                counter = 0;
            }
        }
        /* Are we at the end of the table? */
        else if(ii == (numChans - 1)){
            /* Does this complete a group? */
            if(counter != 0)
                set_incomplete(groupCol, qualCol, ii - counter, ii);
            else
                set_incomplete(groupCol, qualCol, ii, ii);
        }
        /* Are we at the end of a group */
        else if(usedChans[ii + 1]){
            set_incomplete(groupCol, qualCol, ii - counter, ii);
            counter = 0;
        }
        /* Are we at the beginning of a group? */
        else{
            counter++;
        }
    } /* end for(ii) */
    
    free(usedChans);
    return(GRP_SUCCESS);
}

/* MAX_SLOPE
 * ---------
 * Input parameters:
 *   dataCol   - the column representing the independent axis
 *   binCol    - the column with the data values
 *   numChans  - number of channels in groupCol and qualCol
 *   slope     - the maximum slope threshold
 *   groupCol  - the GROUPING column
 *   qualCol   - the QUALITY column
 *   tabStops  - array giving channels with tabs or stops
 *   maxlength - maximum size of groups
 */
int grp_do_max_slope(double *dataCol, double *binCol, long numChans,
                     double slope, short *groupCol, short *qualCol,
                     short *tabStops, double maxlength,
                     dsErrList *errList){
    
    long ii, jj, counter = 0;
    double range = 0.0;
    double tempSlope = 0.0;
    short *usedChans;
    
    /* Check for obviously bad inputs */
    if(!dataCol || !binCol || (numChans < 2) || (slope <= 0)
       || !groupCol || !qualCol || !tabStops){
        if(errList)
            dsErrAdd(errList, dsDMGROUPBADPARAMERR, Accumulation,
                     Generic);
        else
            err_msg("ERROR: At least one input parameter has an "
                    "invalid value.\n");
        return(GRP_ERROR);
    }
    
    if(maxlength <= 0.0)
        maxlength = MAX_BIN_LENGTH;
    
    /* Create and initialize used channel list */
    usedChans = (short *) calloc(numChans, sizeof(short));
    for(ii = 0; ii < numChans; ii++){
        if(tabStops[ii] || (qualCol[ii] != 0))
            usedChans[ii] = GRP_TRUE;
        else
            usedChans[ii] = GRP_FALSE;
    }
    
    ii = 0;
    jj = 1;
    while(ii < (numChans - 1)){
        
        /* Are we in a tab or stop? */
        if(tabStops[ii]){
            ii++;
            jj = (ii + 1);
        }
        else{
            while(jj < numChans){
                
                /* Calculate current slope */
                tempSlope = fabs(((binCol[jj] - binCol[ii])
                                  / (dataCol[jj] - dataCol[ii])));
                
                range = (dataCol[jj] - dataCol[ii]);
                
                /* Are we in a tab or stop? */
                if(tabStops[jj]){
                    ii++;
                    jj = (ii + 1);
                    break;
                }
                /* Are we at the end of the table? */
                else if(jj == (numChans - 1)){
                    /* Does this complete a group? */
                    if((tempSlope >= slope) ||
                       (range >= maxlength)){
                        mark_used(usedChans, ii, jj);
                        create_group(groupCol, ii, jj);
                        set_quality(qualCol, GRP_GOOD, ii, jj);
                        ii = jj;
                        break;
                    }
                    else{
                        ii++;
                        jj = (ii + 1);
                        break;
                    }
                }
                /* Are we at the end of a group or have reached
                 * maxlength? */
                else if((tempSlope >= slope) ||
                        (range >= maxlength)){
                    mark_used(usedChans, ii, jj);
                    create_group(groupCol, ii, jj);
                    set_quality(qualCol, GRP_GOOD, ii, jj);
                    ii = (jj + 1);
                    jj = (ii + 1);
                    break;
                }
                /* Keep looking */
                else jj++;
                
            } /* end while(jj) */
        } /* end if */
    } /* end while(ii) */
    
    /* Put unused channels into "bad" groups */
    for(ii = 0; ii < numChans; ii++){
        /* Are we in a used channel? */
        if(usedChans[ii]){
            if(counter != 0){
                set_incomplete(groupCol, qualCol, ii - counter,
                               ii - 1);
                counter = 0;
            }
        }
        /* Are we at the end of the table? */
        else if(ii == (numChans - 1)){
            /* Does this complete a group? */
            if(counter != 0)
                set_incomplete(groupCol, qualCol, ii - counter, ii);
            else
                set_incomplete(groupCol, qualCol, ii, ii);
        }
        /* Are we at the end of a group */
        else if(usedChans[ii + 1]){
            set_incomplete(groupCol, qualCol, ii - counter, ii);
            counter = 0;
        }
        /* Are we at the beginning of a group? */
        else{
            counter++;
        }
    } /* end for(ii) */
    
    free(usedChans);
    return(GRP_SUCCESS);
}

/*********************************************************************/
/************************  Utility Functions  ************************/
/*********************************************************************/

/* Create tabstops, given tab and stop boundaries */
int create_tabstops(long numChans, double *stopCol, double *tabCol,
		    int isStopColReal, int isTabColReal,
                    double *tBinLow, double *tBinHigh, long tNumBins,
                    double *sBinLow, double *sBinHigh, long sNumBins,
                    short *tabStops, int isAscending, 
		    dsErrList *errList){

    long ii, jj, tempLow, tempHigh, tmpVar;

    /* Initialize tabStops */
    for(ii = 0; ii < numChans; ii++)
        tabStops[ii] = GRP_FALSE;

    /* Go through stops */ 
    for(ii = 0; ii < sNumBins; ii++){
      tempLow = lower_bound(sBinLow[ii], stopCol, numChans, isAscending,
                              errList);

      tempHigh = upper_bound(sBinHigh[ii], stopCol, numChans, isAscending,
                             isStopColReal, errList);

	if (!isAscending) {
	  tmpVar = tempLow;
	  tempLow = tempHigh;
	  tempHigh = tmpVar;
	}

        if((tempLow == GRP_ERROR) || (tempHigh == GRP_ERROR))
            return(GRP_ERROR);

        for(jj = tempLow; jj <= tempHigh; jj++){
            if(jj >= numChans)
                continue;
            else
                tabStops[jj] = GRP_TRUE;
        }
    }
    
    /* Go through tabs */ 
    for(ii = 0; ii < tNumBins; ii++){
        tempLow = lower_bound(tBinLow[ii], tabCol, numChans, isAscending, 
			      errList);
        tempHigh = upper_bound(tBinHigh[ii], tabCol, numChans, isAscending,
                               isTabColReal, errList);

	if (!isAscending) {
	  tmpVar = tempLow;
	  tempLow = tempHigh;
	  tempHigh = tmpVar;
	}

        if((tempLow == GRP_ERROR) || (tempHigh == GRP_ERROR))
            return(GRP_ERROR);
        for(jj = tempLow; jj <= tempHigh; jj++){
            if(jj >= numChans)
                continue;
            else
                tabStops[jj] = GRP_TRUE;
        }
    }

    return(GRP_SUCCESS);
}

/* Set tabs in quality column given tabspec */
int set_tabs(double *dataCol, short *groupCol, short *qualCol,
             long numChans, double *binLow, double *binHigh,
             long numBins, int isAscending, 
	     int isColReal, dsErrList *errList){

    long ii, jj;
    double tempLow, tempHigh, tmpVar;

    /* Go through the binlow-binhigh pairs */
    for(ii = 0; ii < numBins; ii++){
        tempLow = lower_bound(binLow[ii], dataCol, numChans, 
			      isAscending, errList);
        tempHigh = upper_bound(binHigh[ii], dataCol, numChans, 
			       isAscending, isColReal, 
                               errList);

	if (!isAscending) {
	  tmpVar = tempLow;
	  tempLow = tempHigh;
	  tempHigh = tmpVar;
	}

        if((tempLow == GRP_ERROR) || (tempHigh == GRP_ERROR))
            return(GRP_ERROR);
        for(jj = tempLow; jj <= tempHigh; jj++){
            if(jj == tempLow)
                groupCol[jj] = GRP_BEGIN;
            else
                groupCol[jj] = GRP_MIDDLE;
            qualCol[jj] = GRP_TABBED;
        }
    }

    return(GRP_SUCCESS);
}

/* Set stops in grouping & quality columns given stopspec */
int set_stops(double *dataCol, short *groupCol, short *qualCol,
              long numChans, double *binLow, double *binHigh,
              long numBins, int isAscending, 
              int isColReal, dsErrList *errList){

    long ii, jj;
    double tempLow, tempHigh, tmpVar;

    /* Go through the binlow-binhigh pairs */
    for(ii = 0; ii < numBins; ii++){
        tempLow = lower_bound(binLow[ii], dataCol, numChans, 
			      isAscending, errList);
        tempHigh = upper_bound(binHigh[ii], dataCol, numChans, 
			       isAscending, isColReal, errList);
	if (!isAscending) {
	  tmpVar = tempLow;
	  tempLow = tempHigh;
	  tempHigh = tmpVar;
	}

        if((tempLow == GRP_ERROR) || (tempHigh == GRP_ERROR))
            return(GRP_ERROR);
        for(jj = tempLow; jj <= tempHigh; jj++){
            groupCol[jj] = GRP_BEGIN;
            qualCol[jj] = GRP_GOOD;
        }
    }

    return(GRP_SUCCESS);
}

/* Set GRP_DATA column given GROUPING and data column */
int set_grp_data(double *dataCol, short *groupCol, double *grpDataCol,
                 long numChans){

    long ii, jj, stopChan;
    double counter = 0.0;
    
    stopChan = (numChans - 1);
    
    for(ii = (numChans - 1); ii >= 0; ii--){
        if(groupCol[ii] == GRP_UNUSED){
	  grpDataCol[ii] = dataCol[ii];
            counter = 0.0;
            stopChan = (ii - 1);
        }
        else if((groupCol[ii] == GRP_BEGIN) || (ii == 0)){
            counter += dataCol[ii];
            for(jj = ii; jj <= stopChan; jj++){
                grpDataCol[jj] = counter;
            }
            counter = 0.0;
            stopChan = (ii - 1);
        }
        else counter += dataCol[ii];
    }
    
    return(GRP_SUCCESS);
}

/* Set CHANS_PER_GRP column given GROUPING column */
int set_chans_per_grp(short *groupCol, long *chansPerGrpCol,
                      long numChans){

    long ii, jj, stopChan;
    long counter = 1;

    stopChan = (numChans - 1);
    
    for(ii = (numChans - 1); ii >= 0; ii--){
      if (groupCol[ii] == GRP_UNUSED) {
	chansPerGrpCol[ii] = 1;
	stopChan = ii-1;
      } else {
        if((groupCol[ii] == GRP_BEGIN) || (ii == 0)){
	  for(jj = ii; jj <= stopChan; jj++)
	    chansPerGrpCol[jj] = counter;
	  counter = 1;
	  stopChan = (ii - 1);
        }
        else counter++;
      }
    }
    return(GRP_SUCCESS);
}

/* Set GRP_NUM column given GROUPING column */
int set_grp_num(short *groupCol, long *grpNumCol, long numChans){

    long ii;
    long counter = 0;

    for(ii = 0; ii < numChans; ii++){
      if (groupCol[ii] == GRP_UNUSED) {
	grpNumCol[ii] = 0;
      } else {
        if(groupCol[ii] == GRP_BEGIN){
            counter++;
            grpNumCol[ii] = counter;
        } else {
            grpNumCol[ii] = counter;
	}
      }
    }

    return(GRP_SUCCESS);
}

/* Set GRP_STAT_ERR column given GROUPING and STAT_ERR columns */
int set_grp_stat_err(double *grpStatErrCol, short *groupCol,
                     double *statErrCol, long numChans)
{
  
  
  long ii, jj, stopChan;
  double runningTotal = 0.0;
  
  stopChan = numChans -1;
  
  for (ii=0;ii<=stopChan;ii++) {
    
    switch ( groupCol[ii] ) {
      
    case GRP_UNUSED:
      grpStatErrCol[ii] = statErrCol[ii];
      break;
      
    case GRP_BEGIN:
      runningTotal = statErrCol[ii] * statErrCol[ii];
      jj = ii+1;
      while ( (jj<=stopChan) && (groupCol[jj] == GRP_MIDDLE )) {
	runningTotal += statErrCol[jj] * statErrCol[jj];
	jj++;
      }
      runningTotal = sqrt(runningTotal);
      grpStatErrCol[ii] = runningTotal;
      jj = ii+1;
      while ( (jj<=stopChan) && (groupCol[jj] == GRP_MIDDLE )) {
	grpStatErrCol[jj] = runningTotal;
	jj++;
      }
      ii=jj-1;
      
      break;
      
      
    case GRP_MIDDLE:
    default:
      /* Should never get here */
      return( GRP_ERROR );
    }
    
  }
  
  
  return(GRP_SUCCESS);
}

