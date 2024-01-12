# sherpa/optmethods

sources = [
  '__init__.py',
  'ncoresde.py',
  'ncoresnm.py',
  'optfcts.py',
  'opt.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/optmethods'
)

# extension: _saoopt
#
saoopt_sources = [
  'src/Simplex.cc',
  'src/_saoopt.cc'
]

saoopt_include = [
  incdir_numpy,
  '../include',
  '../utils/src/gsl'
]

py.extension_module(
  '_saoopt',
  saoopt_sources,
  subdir: 'sherpa/optmethods',
  include_directories: saoopt_include,
  install: true
)

# extension: _tstoptfct
#
tstoptfct_sources = [
  'tests/_tstoptfct.cc'
]

tstoptfct_include = [
  incdir_numpy,
  '../include'
]

py.extension_module(
  '_tstoptfct',
  tstoptfct_sources,
  subdir: 'sherpa/optmethods',
  include_directories: tstoptfct_include,
  install: true
)

subdir('tests')
