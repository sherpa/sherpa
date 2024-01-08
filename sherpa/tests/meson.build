# sherpa/tests/

sources = [
  'test_data_notebook.py',
  'test_data.py',
  'test_fit_notebook.py',
  'test_fit_unit.py',
  'test_instrument.py',
  'test_sherpa.py',
  'test_stack.py',
  'test_tests.py'
]

data = [
  'a.lis',
  'b.lis'
]

py.install_sources(
  sources,
  subdir: 'sherpa/tests'
)

install_data(
  data,
  install_dir: install_dir / 'sherpa/tests'
)
