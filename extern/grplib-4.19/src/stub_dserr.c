
#ifdef STUB_ERROR_LIB

#include "stub_dserr.h"

void err_msg( const char *blah, ... )
{
  return;
}


dsErrCode dsErrAdd( dsErrList *blah, ... )
{
  return(0);
}


long dsErrGetErrorCt( dsErrList *blah )
{
  return(0);
}

long dsErrGetNumOccur( dsErrList *blah, ... )
{
  return(0);
}

dsErrBool dsErrRemoveN( dsErrList *blah, ... )
{
  return(0);
}

dsErrBool dsErrRemoveAllCode( dsErrList *error_list_p, ... )
{
  return(0);
} 


#endif

