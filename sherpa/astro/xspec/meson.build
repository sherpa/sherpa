# sherpa/astro/xspec

sources = [
  '__init__.py',
  'utils.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa/astro/xspec'
)

# extension: _xspec (optional)
#
if xspec_dep.found()
  xspec_sources = [
    'src/_xspec.cc'
  ]

  xspec_include = [
    incdir_numpy,
    '../../include',
  ]

  # On macOS catching the XSPEC exceptions can fail with the default
  # -fvisibility=hidden option, so switch the setting just for this
  # module (do for all platforms just in case).
  #
  py.extension_module(
    '_xspec',
    xspec_sources,
    cpp_args: xspec_macros + ['-fvisibility=default'],
    subdir: 'sherpa/astro/xspec',
    include_directories: xspec_include,
    dependencies: xspec_dep,
    install: true
  )
endif

subdir('tests')
