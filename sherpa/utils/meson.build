# sherpa/utils

sources = [
  '__init__.py',
  'akima.py',
  'err.py',
  'formatting.py',
  'logging.py',
  'numeric_types.py',
  'parallel.py',
  'random.py',
  'testing.py',
]

py.install_sources(
  sources,
  subdir: 'sherpa/utils'
)

# extension: _utils
#
utils_sources = [
  'src/cephes/const.c',
  'src/cephes/fabs.c',
  'src/cephes/isnan.c',
  'src/cephes/mtherr.c',
  'src/cephes/polevl.c',
  'src/cephes/ndtri.c',
  'src/cephes/gamma.c',
  'src/cephes/igam.c',
  'src/cephes/igami.c',
  'src/cephes/incbet.c',
  'src/cephes/incbi.c',
  'src/sjohnson/Faddeeva.cc',
  'src/_utils.cc'
]

utils_include = [
  incdir_numpy,
  '../include',
  'src/cephes',
  'src/gsl',
  'src/sjohnson'
]

py.extension_module(
  '_utils',
  utils_sources,
  subdir: 'sherpa/utils',
  include_directories: utils_include,
  install: true
)

# extension: integration
#
integration_sources = [
  'src/gsl/err.c',
  'src/gsl/error.c',
  'src/gsl/stream.c',
  'src/gsl/strerror.c',
  'src/gsl/message.c',
  'src/gsl/qng.c',
  'src/sjohnson/adapt_integrate.c',
  'src/integration.cc'
]

integration_include = [
  '../include',
  'src/',
  'src/gsl',
  'src/sjohnson'
]

py.extension_module(
  'integration',
  integration_sources,
  subdir: 'sherpa/utils',
  include_directories: integration_include,
  install: true
)

# extension: _psf
#
psf_sources = [
  'src/tcd/tcdCastArray.c',
  'src/tcd/tcdError.c',
  'src/tcd/tcdFFTConvolve.c',
  'src/tcd/tcdInitConvolveOut.c',
  'src/tcd/tcdInitTransform.c',
  'src/tcd/tcdPadData.c',
  'src/tcd/tcdPixelArith.c',
  'src/tcd/tcdTransform.c',
  'src/_psf.cc'
]

psf_include = [
  incdir_numpy,
  '/home/dburke/micromamba/envs/temp/include',  # TODO
  '../include',
  'src/tcd'
]

psf_link = [
  '-L/home/dburke/micromamba/envs/temp/lib',  # TODO
  '-lfftw3'
]

# TODO: how to pass in the FFTW3 location?
#
py.extension_module(
  '_psf',
  psf_sources,
  subdir: 'sherpa/utils',
  include_directories: psf_include,
  link_args: psf_link,
  install: true
)

subdir('tests')
