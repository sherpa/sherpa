# sherpa/estmethods

sources = [
  '__init__.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/estmethods'
)

# extension: _est_funcs
#
est_funcs_sources = [
  'src/estutils.cc',
  'src/info_matrix.cc',
  'src/projection.cc',
  'src/estwrappers.cc'
]

est_funcs_include = [
  incdir_numpy,
  '../include'
]

py.extension_module(
  '_est_funcs',
  est_funcs_sources,
  subdir: 'sherpa/estmethods',
  include_directories: est_funcs_include,
  install: true
)

subdir('tests')
