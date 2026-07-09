# Top-level build for Sherpa
#
project(
  'sherpa',
  ['cpp', 'c'],

  # From https://github.com/sherpa/sherpa/pull/1948#issuecomment-2692558307
  subproject_dir: 'extern',

  # Follow NumPy and SciPy and modify the version from pyproject.toml
  version: run_command(
    ['scripts/get_sherpa_version'],
    check: true
  ).stdout().strip(),

  meson_version: '>= 1.6.0',
  default_options: [
    'c_std=c11',
    'cpp_std=c++11'
  ]
)

fs = import('fs')

cc = meson.get_compiler('c')

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

# Rely on the subproject installing into the "site-packages"
# directory, which it would be nice not to (see issue gh-50).
# For now the grplib_proj and stklib_proj variables are unused
# (outside of just forcing the code to be built).
#
if get_option('build-group')
  grplib_proj = subproject('grplib-4.18')
endif

if get_option('build-stk')
  stklib_proj = subproject('stklib-4.18')
endif

# Unlike group and stk, the _region module requires a library which
# can either be built as a subproject or taken from an external
# location (e.g. CIAO). So we can build the dependency, take if from
# CIAO, or not build _region. When integrating this dependency we have
# two ways to build the code (taking advantage of other CIAO code or
# not); rather than trying to support all possible options this extra
# question is only checked when we do not build the dependency (since
# the assumption is that if the dependency is being built then the
# CIAO code is not available to link to).
#
if get_option('build-region')
  region_dep = subproject('region-4.18').get_variable('region_dep')
  reg_macros = []

else
  region_prefix = get_option('region-prefix')
  if region_prefix == ''
    region_dep = disabler()
  else
    # Is is likely that region-use-cxc-parser is always going to be
    # true at this point, but allow it to be optional.
    #
    reg_libdir = region_prefix / 'lib'
    reg_incdir = include_directories(region_prefix / 'include')
    reg_dep = cc.find_library(
      'region',
      dirs: reg_libdir,
      has_headers: ['cxcregion.h'],
      header_include_directories: reg_incdir,
      required: false,
      disabler: true
    )

    if get_option('region-use-cxc-parser')
      ascdm_dep = cc.find_library(
        'ascdm',
        dirs: reg_libdir,
        has_headers: ['ascdm.h'],
        header_include_directories: reg_incdir,
        required: true,
      )
      reg_links = [reg_dep, ascdm_dep]
      reg_macros = ['-DUSE_CXCDM_PARSER']
    else
      reg_links = [reg_dep]
      reg_macros = []
    endif

    region_dep = declare_dependency(
      dependencies: reg_links,
      include_directories: reg_incdir
    )
  endif
endif

# Similar to region, the wcs module requires a library
# which can either be built as a subproject or taken from an
# external location (e.g. CIAO).
#
if get_option('build-wcssubs')
  wcssubs_dep = subproject('wcssubs').get_variable('wcssubs_dep')

else
  wcssubs_prefix = get_option('wcssubs-prefix')
  if wcssubs_prefix == ''
    wcssubs_dep = disabler()
  else
    wcs_libdir = wcssubs_prefix / 'lib'
    wcs_incdir = include_directories(wcssubs_prefix / 'include')
    wcs_dep = cc.find_library(
      'wcs',
      dirs: wcs_libdir,
      has_headers: ['wcs.h'],
      header_include_directories: wcs_incdir,
      required: false,
      disabler: true
    )

    wcssubs_dep = declare_dependency(
      dependencies: [wcs_dep],
      include_directories: wcs_incdir
    )
  endif
endif

# The idea is to make FFTW3 an optional dependency, but for now
# require it is available.
#
# Sherpa 4.18.0 and earlier used a built-in copy of FFTW3 if one was
# not present on the system but this has been dropped.
#
fftw3_dep = dependency('fftw3', required : true)

# Configure the XSPEC module. Rather than have a separate
# build-xspec option (to indicate the module should be built)
# and xspec-prefix (which is needed when building the module),
# just use xspec-prefix (leave it as empty not to build).
#
# We could just use $HEADAS rather than require a location, but the
# CXC xspec-modelsonly conda package does interesting things with
# HEADAS which means it can't be easily used.
#
xspec_prefix = get_option('xspec-prefix')
if xspec_prefix != ''
  xspec_incdir = include_directories(xspec_prefix / 'include')
  xspec_libdir = xspec_prefix / 'lib'

  # Find the XSPEC version. There is the possibility that a user
  # wants to specify the version, rather than pick it up directly,
  # but for now do not support this usage.
  #
  xspec_version = run_command(py,
    ['scripts/get_xspec_version'],
    check: true
  ).stdout().strip()

  message('Found XSPEC version:', xspec_version)

  xspec_macros = run_command(py,
    ['scripts/make_xspec_macros', xspec_version],
    check: true
  ).stdout().strip().split()

  message('XSPEC macros:', xspec_macros)

  # Find the libraries
  xspec_libs = get_option('xspec-libraries')
  xspec_lib_deps = []
  foreach xlib : xspec_libs
    xlib_dep = cc.find_library(
      xlib,
      dirs: xspec_libdir,
      header_include_directories: xspec_incdir,
      required: true,
    )
    xspec_lib_deps += xlib_dep
  endforeach

  xspec_dep = declare_dependency(
    dependencies: xspec_lib_deps,
    include_directories: [xspec_incdir],
  )

else
  xspec_dep = disabler()
endif

# Summarize the options.
#
message('-----------------------')
message('Build summary:')

deps = {
  'region': region_dep,
  'wcssubs': wcssubs_dep,
  'fftw3': fftw3_dep,      # this is currently a requirement
  'xspec': xspec_dep
}

foreach label, dep : deps
  if dep.found()
    message('  ', label, ':   true')
  else
    message('  ', label, ':   false')
  endif
endforeach
message('-----------------------')

subdir('sherpa')
