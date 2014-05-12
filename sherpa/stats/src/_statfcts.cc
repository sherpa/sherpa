//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
#include "sherpa/stat_extension.hh"
#include "sherpa/stats.hh"


static PyMethodDef StatFcts[] = {

  STATERRFCT( calc_chi2gehrels_errors ),
  STATERRFCT( calc_chi2constvar_errors ),
  STATERRFCT( calc_chi2datavar_errors ),
  STATERRFCT( calc_chi2xspecvar_errors ),

  STATFCT( calc_cash_stat ),
  STATFCT( calc_cstat_stat ),
  STATFCT( calc_chi2_stat ),
  STATFCT( calc_chi2modvar_stat ),
  STATFCT( calc_lsq_stat ),

  { NULL, NULL, 0, NULL }

};


SHERPAMOD(_statfcts, StatFcts)
