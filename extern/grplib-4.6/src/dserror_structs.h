/*                                                                
**  Copyright (C) 1998,2007  Smithsonian Astrophysical Observatory 
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

#ifndef _DSERROR_STRUCTS_H
#define _DSERROR_STRUCTS_H

/*****************************************************************************/
/*  			DATA STRUCTURES REQUIRED			     */
/*****************************************************************************/

typedef enum dserrbool {dsErrFalse, dsErrTrue} dsErrBool;

typedef enum dserrtype {Accumulation, Individual, Both} dsErrType;

typedef enum dserrseverity {None, Warning, Fatal} dsErrSeverity;

typedef enum dserrmsgtype {Custom, Generic} dsErrMsgType;

typedef enum dserrdatatypes{dsErrShort, dsErrInteger, dsErrLong, dsErrDouble, dsErrLongDouble, dsErrChar, dsErrString, dsErrBad} dsErrDataTypes;

#define dsERRLENGTH  1024
#define dsERRSEVNONE 0
#define dsERRSEVWARNING 1
#define dsERRSEVFATAL 2

typedef long dsErrCode;

/*typedef char dsErrMsg[dsERRLENGTH];*/
typedef char* dsErrMsg;

/* this is a bitmap that specifies which group's errors to include in the
   hash map when the library is initialized.  It must be setup and passed in
   to the initialization function. */
typedef unsigned short dsErrGroup;

/* Error structure, this contains the basic elements that define a given
   error */
typedef struct dserr{
  dsErrCode error_code_t;
  dsErrSeverity error_sev_e;
  dsErrMsg error_msg_t;
} dsErr;

/* Whereas dsErr defines an error in the hypothetical sense, this structure 
   defines an error in a concrete sense - it stores information pertaining
   to a specific instance of an error */
typedef struct dserrinstance{
  dsErr error_t;
  long count;
  dsErrType error_type_e;
} dsErrInstance;

/* error node in the error list */
typedef struct dserrnode{
  dsErrInstance error_instance_t;
  struct dserrnode *next_p;
  struct dserrnode *prev_p;
} dsErrNode;

/* error list */
typedef struct dserrlist{
  dsErrNode *head_p;
  dsErrNode *tail_p;
  long size;
  long contains_fatal;
  long contains_warning;
} dsErrList;

/* memory node in memory management list */
typedef struct dserrmemnode{
  void *data_p;
  struct dserrmemnode *next_p;
} dsErrMemNode;

/* memory management list */
typedef struct dserrmemlist{
  dsErrMemNode *head_p;
} dsErrMemList;

#endif



