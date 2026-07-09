
lgen = generator(
  flex,
  output: '@BASENAME@.c',
  arguments: ['-o', '@OUTPUT@', '@INPUT@']
)

bison_args = [
  '-d',
  '--name-prefix=regYY',
  '@INPUT@',
  '--defines=@OUTPUT1@',
  '--output=@OUTPUT0@'
]

pgen = generator(
  bison,
  output: ['@BASENAME@.c', '@BASENAME@.h'],
  arguments: bison_args
)

reglexer = lgen.process('reglexer.l')
regparser = pgen.process('regparser.y')

sources = [
  'region_area.c',
  'region_convert.c',
  'region_create.c',
  'region_extent.c',
  'region_eval.c',
  'region_logic.c',
  'region_repr.c',
  'region_shape.c',
  'regShape_utils.c',
  'regAnnulus.c',
  'regField.c',
  'regPolygon.c',
  'regSector.c',
  'regMask.c',
  'regCircle.c',
  'regEllipse.c',
  'regBox.c',
  'regPoint.c',
  'regRectangle.c',
  'regPie.c',
  'cxcregion.h',
  'region_priv.h',
  'regpar.h',
  reglexer,
  regparser
  ]

# How many of these are always needed, versus only for linux or macOS?
#
extra_args = [
  '-D_POSIX_C_SOURCE',  # fileno
  '-D_DEFAULT_SOURCE',  # str[n]casecmp, glibc
  '-D_BSD_SOURCE',      # str[n]casecmp, glibc <= 2.19
]

if host_machine.system() == 'darwin'
  extra_args += [
    '-D_DARWIN_C_SOURCE',   # snprintf
  ]
endif

region = static_library(
  'region',
  sources,
  c_args: extra_args,
  install: false
)

region_dep = declare_dependency(
  link_with: region,
  include_directories: include_directories('.')
)
