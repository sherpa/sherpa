# sherpa/astro/xspec

sources = [
  '__init__.py',
  'utils.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/astro/xspec'
)

# extension: _xspec (optional; currently unbuilt)
#
xspec_sources = [
  'src/_xspec.cc'
]

xspec_include = [
  incdir_numpy,
  '../../include',
]

#py.extension_module(
#  '_xspec',
#  xspec_sources,
#  subdir: 'sherpa/astro/xspec',
#  include_directories: xspec_include,
#  install: true
#)

subdir('tests')
