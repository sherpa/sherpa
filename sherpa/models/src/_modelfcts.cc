// 
//  Copyright (C) 2010  Smithsonian Astrophysical Observatory
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

#include "sherpa/model_extension.hh"
#include "sherpa/models.hh"

extern "C" {
  void init_modelfcts();
}

static PyMethodDef ModelFcts[] = {

  MODELFCT1D( box1d, 3  ),
  MODELFCT1D( const1d, 1 ),
  MODELFCT1D( cos, 3 ),
  MODELFCT1D( delta1d, 2 ),
  MODELFCT1D( erf, 3 ),
  MODELFCT1D( erfc, 3 ),
  MODELFCT1D( exp, 3 ),
  MODELFCT1D( exp10, 3 ),
  MODELFCT1D( gauss1d, 3 ),
  MODELFCT1D( log, 3 ),
  MODELFCT1D( log10, 3 ),
  MODELFCT1D( ngauss1d, 3 ),
  MODELFCT1D_NOINT( poisson, 2 ),
  MODELFCT1D( poly1d, 10 ),
  MODELFCT1D_NOINT( logparabola, 4 ),
  MODELFCT1D( powlaw, 3 ),
  MODELFCT1D( sin, 3 ),
  MODELFCT1D( sqrt, 2 ),
  MODELFCT1D( stephi1d, 2 ),
  MODELFCT1D( steplo1d, 2 ),
  MODELFCT1D( tan, 3 ),

  MODELFCT2D( box2d, 5 ),
  MODELFCT2D( const2d, 1 ),
  MODELFCT2D( delta2d, 3 ),
  MODELFCT2D_NOINT( gauss2d, 6 ),
  MODELFCT2D_NOINT( sigmagauss2d, 6 ),
  MODELFCT2D_NOINT( ngauss2d, 6 ),
  MODELFCT2D( poly2d, 9 ),
  
  PY_MODELFCT1D_INT((char*)"integrate1d",
		 (char*)"integrate user functions\n\nExample:\n int_array = integrate1d(func, param_array, xlo_array, xhi_array)" ),

  { NULL, NULL, 0, NULL }

};


SHERPAMODELMOD(_modelfcts, ModelFcts)
