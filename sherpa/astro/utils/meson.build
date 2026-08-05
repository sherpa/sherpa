# sherpa/astro/utils

sources = [
  '__init__.py',
  'smoke.py',
  'xspec.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/astro/utils'
)

# extension: _wcs (optional)
#
if wcssubs_dep.found()
  wcs_sources = [
    'src/_wcs.cc'
  ]

  wcs_include = [
    incdir_numpy,
    '../../include'
  ]

  py.extension_module(
    '_wcs',
    wcs_sources,
    subdir: 'sherpa/astro/utils',
    include_directories: wcs_include,
    dependencies: wcssubs_dep,
    install: true
  )
endif


# extension: _region (optional)
#
if region_dep.found()
  region_sources = [
    'src/_region.cc'
  ]

  region_include = [
    incdir_numpy,
    '../../include'
  ]

  py.extension_module(
    '_region',
    region_sources,
    subdir: 'sherpa/astro/utils',
    include_directories: region_include,
    install: true,
    dependencies: region_dep,
    cpp_args: reg_macros
  )
endif

# extension: _pileup
#
pileup_sources = [
  'src/fftn.c',
  'src/pileup.cc',
  'src/_pileup.cc'
]

pileup_include = [
  incdir_numpy,
  '../../include'
]

py.extension_module(
  '_pileup',
  pileup_sources,
  subdir: 'sherpa/astro/utils',
  include_directories: pileup_include,
  install: true
)

# extension: _utils
#
utils_sources = [
  'src/_utils.cc'
]

utils_include = [
  incdir_numpy,
  '../../include'
]

py.extension_module(
  '_utils',
  utils_sources,
  subdir: 'sherpa/astro/utils',
  include_directories: utils_include,
  install: true
)

subdir('tests')
