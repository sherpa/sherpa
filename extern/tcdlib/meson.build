# tcd

project('tcd',
        ['c'],
	version: '4',
	default_options: [ 'c_std=c11' ]
	)

py = import('python').find_installation(pure: false)

fftw3_dep = dependency('fftw3', required : true)

subdir('src')
