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
