# grplib
# - this is always going to build without liberr support

project('grplib',
        ['c'],
	version: '4.8',
	default_options: [ 'c_std=c11' ]
	)

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
