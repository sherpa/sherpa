
sources = [
  'tcd.h',
  'tcd_private.h',
  'tcdCastArray.c',
  'tcdError.c',
  'tcdFFTConvolve.c',
  'tcdInitConvolveOut.c',
  'tcdInitTransform.c',
  'tcdPadData.c',
  'tcdPixelArith.c',
  'tcdTransform.c',
  ]

tcd = static_library(
  'tcd',
  sources,
  dependencies: [fftw3_dep],
  install: false
)

tcd_dep = declare_dependency(
  link_with: tcd,
  include_directories: include_directories('.'),
)
