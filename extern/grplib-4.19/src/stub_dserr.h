

#ifdef STUB_ERROR_LIB

#ifndef INIT_ERR_LIB_H

typedef struct {
  short zoo;
} dsErrList;

typedef long dsErrCode;

typedef enum dserrbool {dsErrFalse, dsErrTrue} dsErrBool;


extern void      err_msg(const char *, ...);

extern dsErrCode dsErrAdd(dsErrList *error_list_p, ... );


enum { Individual, Accumulation, Generic, 
       dsDMGROUPBADDATAORDERERR,
       dsDMGROUPBADPARAMERR,
       dsDMGROUPEXTRAGROUPSERR,
       dsDMGROUPINVALIDBINERR,
       dsDMGROUPLOWERBOUNDERR,
       dsDMGROUPOVERLAPBINSPECERR,
       dsDMGROUPTOOFEWGROUPSERR,
       dsDMGROUPUPPERBOUNDERR,
       dsDMGROUPZEROERRORERR,
       dsDMGROUPZEROWIDTHERR
};


extern long dsErrGetErrorCt(dsErrList *error_list_p);

extern long dsErrGetNumOccur(dsErrList *error_list_p, ... );

extern dsErrBool dsErrRemoveAllCode(dsErrList *error_list_p, ... ); 


#define INIT_ERR_LIB_H

#endif

#endif

