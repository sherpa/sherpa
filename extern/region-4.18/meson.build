# region

project('region',
        ['c'],
	version: '4.14',
	default_options: [ 'c_std=c11' ]
	)

py = import('python').find_installation(pure: false)

bison = find_program('bison', required: false)
if not bison.found()
  # This path has not been tested
  bison = find_program('byacc', required: true)
endif

flex = find_program('flex')

subdir('src')
