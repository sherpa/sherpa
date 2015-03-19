/*_C_INSERT_SAO_COPYRIGHT_HERE_(1997-2007)_*/
/*_C_INSERT_GPL_LICENSE_HERE_*/

#ifndef GRP_PRIV_H
#define GRP_PRIV_H

#ifndef STUB_ERROR_LIB
#include "liberr.h"
#else
#include "stub_dserr.h"
#endif



/* Function declarations */

long count_bad_chans(short *tabStops, long numChans);

long my_round(double num);

long count_groups(short *groupCol, short *qualCol, long numChans);

void set_incomplete(short *groupCol, short *qualCol, long startChan,
                    long endChan);

void create_group(short *groupCol, long startChan, long endChan);

void set_quality(short *qualCol, short qualVal, long startChan,
                 long endChan);

void mark_used(short *usedChans, long startChan, long endChan);

long lower_bound(double value, double *dataCol, long numChans,
                 int isAscending, dsErrList *errList);

long upper_bound(double value, double *dataCol, long numChans,
                 int isAscending, int isColReal, dsErrList *errList);

int check_overlap(double *binLow, double *binHigh, long numBins);


#endif
