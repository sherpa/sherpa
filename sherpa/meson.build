# sherpa

sources = [
  '__init__.py',
  '_version.py',        # this should be auto-generated
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
