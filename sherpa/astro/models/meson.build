# sherpa/astro/models

sources = [
  '__init__.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/astro/models'
)

# extension: _modelfcts
#
modelfcts_sources = [
  '../../utils/src/sjohnson/Faddeeva.cc',
  'src/_modelfcts.cc'
]

modelfcts_include = [
  incdir_numpy,
  '../../include',
  '../../utils/src/sjohnson'
]

py.extension_module(
  '_modelfcts',
  modelfcts_sources,
  subdir: 'sherpa/astro/models',
  include_directories: modelfcts_include,
  install: true
)

subdir('tests')
