
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

# Hide compiler warnings
extra_args = [
  '-Wno-deprecated-non-prototype',
  '-Wno-format',
  '-Wno-misleading-indentation',
]

if cc.get_id() == 'gcc'
  extra_args += [
    '-Wno-format-overflow',
    '-Wno-stringop-truncation',
    '-Wno-stringop-overflow',
    '-Wno-unused-result',
    '-Wno-maybe-uninitialized'
  ]
elif cc.get_id() == 'clang'
  extra_args += [
    '-Wno-array-parameter',
    '-Wno-absolute-value',
    '-Wno-sometimes-uninitialized',
    '-Wno-constant-logical-operand'
  ]
endif

wcssubs = static_library(
  'wcssubs',
  sources,
  c_args: extra_args,
  install: false
)

wcssubs_dep = declare_dependency(
  link_with: wcssubs,
  include_directories: include_directories('.')
)
