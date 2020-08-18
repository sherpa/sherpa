// 
//  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

#ifndef __sherpa_projection_hh__
#define __sherpa_projection_hh__
int projection(double* original_pars, const int op_size,
	       const double* pars_mins, const int mins_size,
	       const double* pars_maxs, const int maxs_size,
	       const double* pars_hardmins, const int hmins_size,
	       const double* pars_hardmaxs, const int hmaxs_size,
	       double* pars_elow, const int elow_size,
	       double* pars_ehi, const int ehi_size,
	       const double sigma, 
	       const double eps, 
	       const int maxiters,
	       const double remin,
	       double (*statfcn)(double*, int),
	       double (*fitfcn)(double (*statfcn)(double*, int),
				double*,double*,double*,int,int)) throw();
#endif
