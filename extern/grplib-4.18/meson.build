# grplib
# - this is always going to build without liberr support

project('grplib',
        ['c'],
	version: '4.8',
	default_options: [
	  'c_std=c11',
	  'python.allow_limited_api=false',
	]
	)

# If a limited API is used, what is the version to use (this needs to
# be updated to match the minimum-supported Python version in pyproject.toml).
#
limited_api_version = '3.11'

py = import('python').find_installation(pure: false)

incdir_numpy = run_command(py,
    [
      '-c',
      'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
    ],
    check: true
  ).stdout().strip()

subdir('src')
subdir('python')
