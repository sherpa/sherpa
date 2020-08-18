#
#  Copyright (C) 2014, 2016, 2017, 2018, 2020
#     Smithsonian Astrophysical Observatory
#
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#


from numpy.distutils.core import Extension

# Include directory for Sherpa headers
sherpa_inc = ['src/sherpa/include', 'src/sherpa/utils/src']

header_deps = {
    'myArray': (),
    'array': (),
    'constants': (),
    'extension': ('array',),
    'integration': (),
    'model_extension': ('extension', 'integration'),
    'models': ('constants', 'utils'),
    'stat_extension': ('extension',),
    'stats': ('utils',),
    'utils': ('constants','extension'),
    'astro/models': ('constants', 'utils'),
    'astro/utils': (),
    'astro/xspec_extension': ('extension',),
    }

def get_deps(deps):
    deps = set(deps)
    alldeps = set()

    while deps:
        dep = deps.pop()
        if dep not in alldeps:
            alldeps.add(dep)
            deps.update(header_deps[dep])

    return [sherpa_inc[0] + '/sherpa/' + d + '.hh' for d in alldeps]

####
# EXTENSIONS WITH EXTERNAL DEPENDENCIES
####

def build_psf_ext(library_dirs, include_dirs, libraries):
    return Extension('sherpa.utils._psf',
             ['src/sherpa/utils/src/tcd/tcdCastArray.c',
              'src/sherpa/utils/src/tcd/tcdError.c',
              'src/sherpa/utils/src/tcd/tcdFFTConvolve.c',
              'src/sherpa/utils/src/tcd/tcdInitConvolveOut.c',
              'src/sherpa/utils/src/tcd/tcdInitTransform.c',
              'src/sherpa/utils/src/tcd/tcdPadData.c',
              'src/sherpa/utils/src/tcd/tcdPixelArith.c',
              'src/sherpa/utils/src/tcd/tcdTransform.c',
              'src/sherpa/utils/src/_psf.cc'],
             sherpa_inc + ['src/sherpa/utils/src/tcd'] + include_dirs,
             library_dirs=library_dirs,
             libraries=libraries,
             depends=(get_deps(['extension', 'utils'])+
                      ['src/sherpa/utils/src/tcd/tcd.h',]))

def build_wcs_ext(library_dirs, include_dirs, libraries):
    return Extension('sherpa.astro.utils._wcs',
                 ['src/sherpa/astro/utils/src/_wcs.cc'],
                 sherpa_inc + include_dirs,
                 library_dirs=library_dirs,
                 libraries=libraries,
                 depends=get_deps(['extension']))

def build_region_ext(library_dirs, include_dirs, libraries):
    return Extension('sherpa.astro.utils._region',
                 ['src/sherpa/astro/utils/src/_region.cc'],
                 sherpa_inc + include_dirs,
                 library_dirs=library_dirs,
                 libraries=(libraries),
                 depends=get_deps(['extension']))

def build_xspec_ext(library_dirs, include_dirs, libraries, define_macros=None):
    return Extension('sherpa.astro.xspec._xspec',
                  ['src/sherpa/astro/xspec/src/_xspec.cc'],
                  sherpa_inc + include_dirs,
                  library_dirs=library_dirs,
                  runtime_library_dirs=library_dirs,
                  libraries=libraries,
                  define_macros=define_macros,
                  depends=(get_deps(['astro/xspec_extension'])))

def build_ext(name, library_dirs, include_dirs, libraries, **kwargs):
    func = globals().get('build_'+name+'_ext')
    return func(library_dirs, include_dirs, libraries, **kwargs)


def build_lib_arrays(command, libname):
    library_dirs = getattr(command, libname+'_lib_dirs').split(' ')
    include_dirs = getattr(command, libname+'_include_dirs').split(' ')
    libraries = getattr(command, libname+'_libraries').split(' ')
    return [library_dirs, include_dirs, libraries]

###
# Static Extensions
###

estmethods = Extension('sherpa.estmethods._est_funcs',
              ['src/sherpa/estmethods/src/estutils.cc',
               'src/sherpa/estmethods/src/info_matrix.cc',
               'src/sherpa/estmethods/src/projection.cc',
               'src/sherpa/estmethods/src/estwrappers.cc'],
              (sherpa_inc + ['src/sherpa/utils/src/gsl']),
              depends=(get_deps(['extension', 'utils']) +
                       ['src/sherpa/estmethods/src/estutils.hh',
                        'src/sherpa/estmethods/src/info_matrix.hh',
                        'src/sherpa/estmethods/src/projection.hh',
                        'src/sherpa/utils/src/gsl/fcmp.h']))


utils = Extension('sherpa.utils._utils',
              ['src/sherpa/utils/src/cephes/const.c',
               'src/sherpa/utils/src/cephes/fabs.c',
               'src/sherpa/utils/src/cephes/isnan.c',
               'src/sherpa/utils/src/cephes/mtherr.c',
               'src/sherpa/utils/src/cephes/polevl.c',
               'src/sherpa/utils/src/cephes/ndtri.c',
               'src/sherpa/utils/src/cephes/gamma.c',
               'src/sherpa/utils/src/cephes/igam.c',
               'src/sherpa/utils/src/cephes/igami.c',
               'src/sherpa/utils/src/cephes/incbet.c',
               'src/sherpa/utils/src/cephes/incbi.c',
               'src/sherpa/utils/src/_utils.cc'],
              sherpa_inc + ['src/sherpa/utils/src/cephes',
                            'src/sherpa/utils/src/gsl'],
              depends=(get_deps(['extension', 'utils'])+
                       ['src/sherpa/utils/src/gsl/fcmp.h',
                        'src/sherpa/utils/src/cephes/cephes.h']))

modelfcts = Extension('sherpa.models._modelfcts',
              ['src/sherpa/models/src/_modelfcts.cc'],
              sherpa_inc,
              depends=get_deps(['model_extension', 'models']))

saoopt = Extension('sherpa.optmethods._saoopt',
              ['src/sherpa/optmethods/src/_saoopt.cc',
               'src/sherpa/optmethods/src/Simplex.cc'],
              sherpa_inc + ['src/sherpa/utils/src/gsl'],
              depends=(get_deps(['myArray', 'extension']) +
                       ['src/sherpa/include/sherpa/fcmp.hh',
                        'src/sherpa/include/sherpa/MersenneTwister.h',
                        'src/sherpa/include/sherpa/functor.hh',
                        'src/sherpa/optmethods/src/DifEvo.hh',
                        'src/sherpa/optmethods/src/DifEvo.cc',
                        'src/sherpa/optmethods/src/NelderMead.hh',
                        'src/sherpa/optmethods/src/NelderMead.cc',
                        'src/sherpa/optmethods/src/Opt.hh',
                        'src/sherpa/optmethods/src/PyWrapper.hh',
                        'src/sherpa/optmethods/src/RanOpt.hh',
                        'src/sherpa/optmethods/src/Simplex.hh',
                        'src/sherpa/optmethods/src/Simplex.cc',
                        'src/sherpa/optmethods/src/minpack/LevMar.hh',
                        'src/sherpa/optmethods/src/minpack/LevMar.cc',
                        'src/sherpa/optmethods/src/minim.hh']))

tstoptfct = Extension('sherpa.optmethods._tstoptfct',
              ['src/sherpa/optmethods/tests/_tstoptfct.cc'],
              sherpa_inc,
              depends=(get_deps(['extension']) +
                       ['src/sherpa/include/sherpa/fcmp.hh',
                        'src/sherpa/include/sherpa/MersenneTwister.h',
                        'src/sherpa/include/sherpa/functor.hh',
                        'src/sherpa/optmethods/tests/tstopt.hh',
                        'src/sherpa/optmethods/tests/tstoptfct.hh',
                        'src/sherpa/optmethods/src/DifEvo.hh',
                        'src/sherpa/optmethods/src/DifEvo.cc',
                        'src/sherpa/optmethods/src/NelderMead.hh',
                        'src/sherpa/optmethods/src/NelderMead.cc',
                        'src/sherpa/optmethods/src/Opt.hh',
                        'src/sherpa/optmethods/src/PyWrapper.hh',
                        'src/sherpa/optmethods/src/RanOpt.hh',
                        'src/sherpa/optmethods/src/Simplex.hh',
                        'src/sherpa/optmethods/src/Simplex.cc',
                        'src/sherpa/optmethods/src/minpack/LevMar.hh',
                        'src/sherpa/optmethods/src/minpack/LevMar.cc']))

statfcts = Extension('sherpa.stats._statfcts',
              ['src/sherpa/stats/src/_statfcts.cc'],
              sherpa_inc,
              depends=get_deps(['stat_extension', 'stats']))

integration = Extension('sherpa.utils.integration',
              ['src/sherpa/utils/src/gsl/err.c',
               'src/sherpa/utils/src/gsl/error.c',
               'src/sherpa/utils/src/gsl/stream.c',
               'src/sherpa/utils/src/gsl/strerror.c',
               'src/sherpa/utils/src/gsl/message.c',
               'src/sherpa/utils/src/gsl/qng.c',
               'src/sherpa/utils/src/adapt_integrate.c',
               'src/sherpa/utils/src/integration.cc'],
              sherpa_inc + ['src/sherpa/utils/src',
                            'src/sherpa/utils/src/gsl'],
              depends=(get_deps(['integration'])+
                       ['src/sherpa/utils/src/adapt_integrate.h',
                        'src/sherpa/utils/src/gsl/gsl_integration.h']))

astro_modelfcts = Extension('sherpa.astro.models._modelfcts',
              ['src/sherpa/astro/models/src/_modelfcts.cc'],
              sherpa_inc,
              depends=get_deps(['model_extension', 'astro/models']))

pileup = Extension('sherpa.astro.utils._pileup',
              ['src/sherpa/astro/utils/src/fftn.c',
               'src/sherpa/astro/utils/src/_pileup.cc',
               'src/sherpa/astro/utils/src/pileup.cc'],
              sherpa_inc + ['src/sherpa/astro/utils/src'],
              depends=(get_deps(['extension']) +
                       ['src/sherpa/astro/utils/src/pileup.hh',
                        'src/sherpa/astro/utils/src/PyWrapper.hh',
                        'src/sherpa/astro/utils/src/fftn.inc']))

astro_utils = Extension('sherpa.astro.utils._utils',
              ['src/sherpa/astro/utils/src/_utils.cc'],
              (sherpa_inc + ['src/sherpa/utils/src/gsl']),
              depends=(get_deps(['extension', 'utils', 'astro/utils'])+
                       ['src/sherpa/utils/src/gsl/fcmp.h']))


static_ext_modules = [
                   estmethods,
                   utils,
                   modelfcts,
                   saoopt,
                   tstoptfct,
                   statfcts,
                   integration,
                   astro_modelfcts,
                   pileup,
                   astro_utils,
                ]
