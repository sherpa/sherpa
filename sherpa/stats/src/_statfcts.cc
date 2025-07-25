//
//  Copyright (C) 2007, 2015, 2016, 2025
//  Smithsonian Astrophysical Observatory
//
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//


#include "sherpa/stat_extension.hh"
#include "sherpa/stats.hh"

static PyMethodDef StatFcts[] = {

  STATERRFCT( calc_chi2gehrels_errors ),
  STATERRFCT( calc_chi2constvar_errors ),
  STATERRFCT( calc_chi2datavar_errors ),
  STATERRFCT( calc_chi2xspecvar_errors ),

  STATFCT( calc_chi2_stat ),
  STATFCT( calc_chi2modvar_stat ),
  STATFCT( calc_lsq_stat ),

  LKLHD_STATFCT( calc_cash_stat ),
  LKLHD_STATFCT( calc_cstat_stat ),
  WSTATFCT( calc_wstat_stat ),

  { NULL, NULL, 0, NULL }

};


SHERPAMOD(_statfcts, StatFcts)
