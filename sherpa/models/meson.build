# sherpa/models

sources = [
  '__init__.py',
  'basic.py',
  'model.py',
  'op.py',
  'parameter.py',
  'regrid.py',
  'template.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/models'
)

# extension: _modelfcts
#
modelfcts_sources = [
  'src/_modelfcts.cc'
]

modelfcts_include = [
  incdir_numpy,
  '../include',
]

py.extension_module(
  '_modelfcts',
  modelfcts_sources,
  subdir: 'sherpa/models',
  include_directories: modelfcts_include,
  install: true
)

subdir('tests')
