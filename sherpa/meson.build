# sherpa

# Generate version.py for sdist
# TODO: add support for --meson-dist
meson.add_dist_script(
   ['../scripts/get_sherpa_version', '--meson-dist',
    '--module', 'sherpa/_version.py']
)
if not fs.exists('_version.py')
  generate_version = custom_target(
    'generate-version',
    install: true,
    build_always_stale: true,
    build_by_default: true,
    output: '_version.py',
    input: '../scripts/get_sherpa_version',
    command: [py, '@INPUT@', '--module', '@OUTPUT@'],
    install_dir: py.get_install_dir() / 'sherpa',
    install_tag: 'python-runtime'
  )
else
  # When building from sdist, _version.py exists and should be included
  py.install_sources(
    ['_version.py'],
    subdir : 'sherpa'
  )
endif

sources = [
  '__init__.py',
  'all.py',
  'conftest.py',
  'data.py',
  'fit.py',
  'instrument.py',
  'io.py',
  'logposterior.py'
]

py.install_sources(
  sources,
  subdir: 'sherpa',
)

data = [
  'sherpa.rc',
  'sherpa-standalone.rc'
]

install_data(
  data,
  install_dir: install_dir / 'sherpa'
)

subdir('astro')
subdir('estmethods')
subdir('image')
subdir('models')
subdir('optmethods')
subdir('plot')
subdir('sim')
subdir('static')
subdir('stats')
subdir('ui')
subdir('utils')

subdir('include')

subdir('tests')
