# region

project('region',
        ['c'],
	version: '4.14',
	default_options: [ 'c_std=c11' ]
	)

py = import('python').find_installation(pure: false)

bison = find_program('bison', required: false)
bison_no_deprecated = false
if bison.found()
  if find_program('bison', required: false, version: '> 2.3').found()
    bison_no_deprecated = true
  endif
else
  bison = find_program('byacc', required: true)
endif

flex = find_program('flex')

subdir('src')
