

#ifdef STUB_ERROR_LIB

#ifndef INIT_ERR_LIB_H

typedef struct {
  short zoo;
} dsErrList;

#define dsErrCode int

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


#define INIT_ERR_LIB_H

#endif

#endif

