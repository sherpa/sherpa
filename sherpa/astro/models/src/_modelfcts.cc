//_C++_INSERT_SAO_COPYRIGHT_HERE_(2007)_
//_C++_INSERT_GPL_LICENSE_HERE_
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
