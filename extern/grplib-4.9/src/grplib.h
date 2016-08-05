/*                                                                
**  Copyright (C) 1997-2007  Smithsonian Astrophysical Observatory 
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


#ifndef GRP_LIB_H
#define GRP_LIB_H


#ifndef STUB_ERROR_LIB
#include "liberr.h"
#else
#include "stub_dserr.h"
#endif


#include <math.h>
#include <stdlib.h>
#include <float.h>
 
/* Maximum bin length default */
#define MAX_BIN_LENGTH DBL_MAX

/* Return values */
enum return_val {GRP_ERROR = -1, GRP_SUCCESS, GRP_WARNING};

/* Quality column values */
enum quality_val {GRP_GOOD = 0, GRP_POOR = 2, GRP_TABBED = 5};

/* Grouping column values */
enum grouping_val {GRP_MIDDLE = -1, GRP_UNUSED, GRP_BEGIN};

/* Boolean values */
enum boolean_val {GRP_FALSE, GRP_TRUE};


/* Function declarations */

int grp_do_num_bins(long numChans, long numBins, short *groupCol,
                    short *qualCol, short *tabStops,
                    dsErrList *errList);

int grp_do_num_counts(double *dataCol, long numChans, double numCounts,
                      short *groupCol, short *qualCol, short *tabStops,
                      double maxLength, dsErrList *errList);

int grp_do_bin(double *dataCol, long numChans, double *binLow,
               double *binHigh, long numBins, short *groupCol,
               short *qualCol, short *tabStops, dsErrList *errList,
               short partialBin, int isColReal);

int grp_do_bin_width(long numChans, long binWidth, short *groupCol,
                     short *qualCol, short *tabStops,
                     dsErrList *errList);

int grp_do_snr(double *dataCol, long numChans, double snr,
               short *groupCol, short *qualCol, short *tabStops,
               double *errorCol, short useErr, double maxLength,
               dsErrList *errList);

int grp_do_adaptive(double *dataCol, long numChans, double minCounts,
                    short *groupCol, short *qualCol, short *tabStops,
                    double maxLength, dsErrList *errList);

int grp_do_adaptive_snr(double *dataCol, long numChans, double snr,
                        short *groupCol, short *qualCol,
                        short *tabStops, double *errorCol,
                        short useErr, double maxLength,
                        dsErrList *errList);

int grp_do_min_slope(double *dataCol, double *binCol, long numChans,
                     double slope, short *groupCol, short *qualCol,
                     short *tabStops, double maxlength,
                     dsErrList *errList);

int grp_do_max_slope(double *dataCol, double *binCol, long numChans,
                     double slope, short *groupCol, short *qualCol,
                     short *tabStops, double maxlength,
                     dsErrList *errList);

int grp_do_bin_file(double *dataCol, long numChans, short *groupCol,
                     short *qualCol, short *tabStops,
                     double *fDataCol, long fNumChans,
                     short *fGroupCol, short *fQualCol,
                     int isColReal, dsErrList *errList);

int grp_do_none(long numChans, short *groupCol, short *qualCol,
                dsErrList *errList);

int create_tabstops(long numChans, double *stopCol, double *tabCol,
		    int isStopColReal, int isTabColReal,
                    double *tBinLow, double *tBinHigh, long tNumBins,
                    double *sBinLow, double *sBinHigh, long sNumBins,
                    short *tabStops, int isAscending,
		    dsErrList *errList);

int set_tabs(double *dataCol, short *groupCol, short *qualCol,
             long numChans, double *binLow, double *binHigh,
             long numBins, int isAscending, int isColReal,
	     dsErrList *errList);

int set_stops(double *dataCol, short *groupCol, short *qualCol,
              long numChans, double *binLow, double *binHigh,
              long numBins, int isAscending, int isColReal,
	      dsErrList *errList);

int set_grp_data(double *dataCol, short *groupCol, double *grpDataCol,
                 long numChans);

int set_chans_per_grp(short *groupCol, long *chansPerGrpCol,
                      long numChans);

int set_grp_num(short *groupCol, long *grpNumCol, long numChans);

int set_grp_stat_err(double *grpStatErrCol, short *groupCol,
                     double *statErrCol, long numChans);
int check_increasing(double *dataCol, long numChans);

int check_decreasing(double *dataCol, long numChans);

#endif
