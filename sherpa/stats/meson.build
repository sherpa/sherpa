# sherpa/stats

sources = [
  '__init__.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/stats'
)

# extension: _statfcts
#
statfcts_sources = [
  'src/_statfcts.cc'
]

statfcts_include = [
  incdir_numpy,
  '../include',
]

py.extension_module(
  '_statfcts',
  statfcts_sources,
  subdir: 'sherpa/stats',
  include_directories: statfcts_include,
  install: true
)

subdir('tests')
