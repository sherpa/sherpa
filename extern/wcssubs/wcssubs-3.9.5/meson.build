
sources = [
  'cel.c',
  'distort.c',
  'fitsfile.c',
  'hget.c',
  'hput.c',
  'iget.c',
  'imhfile.c',
  'imio.c',
  'lin.c',
  'platepos.c',
  'poly.c',
  'proj.c',
  'sph.c',
  'tnxpos.c',
  'zpxpos.c',
  'wcs.c',
  'wcsinit.c',
  'wcscon.c',
  'wcslib.c',
  'wcstrig.c',
  'worldpos.c',
  'dateutil.c',
  'fileutil.c',
  'dsspos.c',
  'wcslib.h',
  'wcs.h',
  'fitshead.h',
  'fitsfile.h'
]

wcssubs = static_library(
  'wcssubs',
  sources,
  install: false
)

wcssubs_dep = declare_dependency(
  link_with: wcssubs,
  include_directories: include_directories('.')
)
