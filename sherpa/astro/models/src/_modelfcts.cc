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

#define _MODELFCTPTR(name) \
  sherpa::astro::models::name< SherpaFloat, SherpaFloatArray >

#include "sherpa/model_extension.hh"
#include "sherpa/astro/models.hh"

extern "C" {
  void init_modelfcts();
}


static PyMethodDef ModelFcts[] = {

  MODELFCT1D_NOINT( atten, 3 ),
  MODELFCT1D_NOINT( bbody, 3 ),
  MODELFCT1D_NOINT( bbodyfreq, 2 ),
  MODELFCT1D_NOINT( beta1d, 4 ),
  MODELFCT1D( bpl1d, 5 ),
  MODELFCT1D_NOINT( dered, 2 ),
  MODELFCT1D_NOINT( edge, 3 ),
  MODELFCT1D( linebroad, 3 ),
  MODELFCT1D( lorentz1d, 3 ),
  MODELFCT1D_NOINT( nbeta1d, 4 ),
  MODELFCT1D( schechter, 3 ),

  MODELFCT2D_NOINT( beta2d, 7 ),
  MODELFCT2D_NOINT( devau, 6 ),
  MODELFCT2D_NOINT( sersic, 7 ),
  MODELFCT2D_NOINT( hr, 6 ),
  MODELFCT2D_NOINT( lorentz2d, 6 ),

  { NULL, NULL, 0, NULL }

};


SHERPAMODELMOD(_modelfcts, ModelFcts)
