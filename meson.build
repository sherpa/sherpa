# Top-level build for Sherpa
#
project(
  'sherpa',
  ['cpp', 'c'],
  # scipy claims we can not include the git commit hash here
  version: '4.16.1.dev0',
  # scipy (January 2024) uses this lower limit, so let's use it too
  meson_version: '>= 1.1.0',
  default_options: [
    # what are good minimums?
    'c_std=c99',
    'cpp_std=c++11'
  ]
)

py = import('python').find_installation(pure: false)

# For when we want to install non-Python files
#
install_dir = py.get_install_dir()

incdir_numpy = run_command(py,
    [
      '-c',
      'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
    ],
    check: true
  ).stdout().strip()

subdir('sherpa')
