# fftw
# - I do not want to convert this to meson ...
#

project('fftw',
        ['c'],
	version: '3.3.10',
	default_options: [
	  'c_std=c11',
	  'default_library=shared'
	],
	# only need 0.56.0 but follow the sherpa package
        meson_version: '>= 1.6.0',
	)

mod = import('unstable-external_project')

copts = [
]
p = mod.add_project('configure',
  configure_options: copts
)

fftw3_dep = p.dependency('fftw3')
meson.override_dependency('fftw3', fftw3_dep)
