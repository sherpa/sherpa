#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2014)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

#import setuptools
from numpy.distutils.core import setup, Extension

from setup_helpers import get_deps, sherpa_inc, build, clean, install, sdist


# FIXME THIS SECTION IS TEMPORARY AND PART OF IT WILL PROBABLY BE MOVED TO SETUP.CFG
extern_inc='build/include'
extern_lib='build/lib'

# COMMON METADATA
meta = {
        'name' : 'sherpa',
        'version' : '4.6.2',
        'author' : 'Smithsonian Astrophysical Observatory / Chandra X-Ray Center',
        'author_email' : 'cxchelp@head.cfa.harvard.edu',
        'url' : 'http://cxc.harvard.edu/sherpa/',
        'description' : 'Modeling and fitting package for scientific data analysis',
        'license' : 'GNU GPL v3',
        'long_description' : 'Modeling and fitting package for scientific data analysis',
        'platforms' : 'Linux, Mac OS X',
        'install_requires' : ['numpy', 'pyfits', 'matplotlib'],

        }


build_xspec = False
# END TEMPORARY CONFIGURATION SECTION

packages=['sherpa',
                'sherpa.estmethods',
                'sherpa.image',
                'sherpa.models',
                'sherpa.optmethods',
                'sherpa.plot',
                'sherpa.sim',
                'sherpa.stats',
                'sherpa.ui',
                'sherpa.utils',
                'sherpa.astro',
                'sherpa.astro.io',
                'sherpa.astro.models',
                'sherpa.astro.optical',
                'sherpa.astro.sim',
                'sherpa.astro.ui',
                'sherpa.astro.utils',
                ]

meta['packages'] = packages

package_data={'sherpa': ['include/sherpa/*.hh',
                               'include/sherpa/astro/*.hh',
                               'tests/test_*.py'],
                    'sherpa.estmethods': ['tests/test_*.py'],
                    'sherpa.image': ['tests/test_*.py'],
                    'sherpa.models': ['tests/test_*.py'],
                    'sherpa.optmethods': ['tests/test_*.py'],
                    'sherpa.plot': ['tests/test_*.py'],
                    'sherpa.sim': ['tests/test_*.py'],
                    'sherpa.stats': ['tests/test_*.py'],
                    'sherpa.ui': ['tests/test_*.py'],
                    'sherpa.utils': ['tests/test_*.py'],
                    'sherpa.astro': ['tests/test_*.py'],
                    'sherpa.astro.io': ['tests/test_*.py'],
                    'sherpa.astro.models': ['tests/test_*.py'],
                    'sherpa.astro.optical': ['tests/test_*.py'],
                    'sherpa.astro.sim': ['tests/test_*.py'],
                    'sherpa.astro.ui': ['tests/test_*.py'],
                    'sherpa.astro.utils': ['tests/test_*.py'],
                    }

meta['package_data'] = package_data

data_files=[('sherpa', ['sherpa/sherpa.rc']),
            ('', ['build/lib/group.so']),
            ]

meta['data_files'] = data_files

if build_xspec:
    packages.append('sherpa.astro.xspec')
    package_data['sherpa.astro.xspec': ['tests/test_*.py']]

#####
#Extensions
#####

clibs = [

    ('sherpa',
     {'sources': ['sherpa/utils/src/gsl/fcmp.c'],
      'sourceDir' : 'sherpa/utils',
      'libs' : [],
      'libdirs' : [],
      'include_dirs': ['sherpa/utils/src'],
      'headerExportDir' : [],
      })
    ]

meta['libraries'] = clibs

estmethods = Extension('sherpa.estmethods._est_funcs',
              ['sherpa/estmethods/src/estutils.cc',
               'sherpa/estmethods/src/info_matrix.cc',
               'sherpa/estmethods/src/projection.cc',
               'sherpa/estmethods/src/estwrappers.cc'],
              (sherpa_inc + ['sherpa/utils/src/gsl']),
              libraries=(['sherpa']),
              depends=(get_deps(['extension', 'utils']) +
                       ['sherpa/estmethods/src/estutils.hh',
                        'sherpa/estmethods/src/info_matrix.hh',
                        'sherpa/estmethods/src/projection.hh',
                        'sherpa/utils/src/gsl/fcmp.h']))


utils = Extension('sherpa.utils._utils',
              ['sherpa/utils/src/cephes/const.c',
               'sherpa/utils/src/cephes/fabs.c',
               'sherpa/utils/src/cephes/isnan.c',
               'sherpa/utils/src/cephes/mtherr.c',
               'sherpa/utils/src/cephes/polevl.c',
               'sherpa/utils/src/cephes/ndtri.c',
               'sherpa/utils/src/cephes/gamma.c',
               'sherpa/utils/src/cephes/igam.c',
               'sherpa/utils/src/cephes/igami.c',
               'sherpa/utils/src/cephes/incbet.c',
               'sherpa/utils/src/cephes/incbi.c',
               'sherpa/utils/src/_utils.cc'],
              sherpa_inc + ['sherpa/utils/src/cephes',
                            'sherpa/utils/src/gsl'],
              libraries=(['sherpa']),
              depends=(get_deps(['extension', 'utils'])+
                       ['sherpa/utils/src/gsl/fcmp.h',
                        'sherpa/utils/src/cephes/cephes.h']))

modelfcts = Extension('sherpa.models._modelfcts',
              ['sherpa/models/src/_modelfcts.cc'],
              sherpa_inc,
              depends=get_deps(['model_extension', 'models']))

saoopt = Extension('sherpa.optmethods._saoopt',
              ['sherpa/optmethods/src/_saoopt.cc',
               'sherpa/optmethods/src/Simplex.cc'],
              sherpa_inc + ['sherpa/utils/src/gsl'],
              libraries=(['sherpa']),
              depends=(get_deps(['myArray', 'extension']) +
                       ['sherpa/include/sherpa/fcmp.hh',
                        'sherpa/include/sherpa/MersenneTwister.h',
                        'sherpa/include/sherpa/functor.hh',
                        'sherpa/optmethods/src/DifEvo.hh',
                        'sherpa/optmethods/src/DifEvo.cc',
                        'sherpa/optmethods/src/NelderMead.hh',
                        'sherpa/optmethods/src/NelderMead.cc',
                        'sherpa/optmethods/src/Opt.hh',
                        'sherpa/optmethods/src/PyWrapper.hh',
                        'sherpa/optmethods/src/RanOpt.hh',
                        'sherpa/optmethods/src/Simplex.hh',
                        'sherpa/optmethods/src/Simplex.cc',
                        'sherpa/optmethods/src/minpack/LevMar.hh',
                        'sherpa/optmethods/src/minpack/LevMar.cc']))

tstoptfct = Extension('sherpa.optmethods._tstoptfct',
              ['sherpa/optmethods/tests/_tstoptfct.cc'],
              sherpa_inc,
              depends=(get_deps(['extension']) +
                       ['sherpa/include/sherpa/fcmp.hh',
                        'sherpa/include/sherpa/MersenneTwister.h',
                        'sherpa/include/sherpa/functor.hh',
                        'sherpa/optmethods/src/DifEvo.hh',
                        'sherpa/optmethods/src/DifEvo.cc',
                        'sherpa/optmethods/src/NelderMead.hh',
                        'sherpa/optmethods/src/NelderMead.cc',
                        'sherpa/optmethods/src/Opt.hh',
                        'sherpa/optmethods/src/PyWrapper.hh',
                        'sherpa/optmethods/src/RanOpt.hh',
                        'sherpa/optmethods/src/Simplex.hh',
                        'sherpa/optmethods/src/Simplex.cc',
                        'sherpa/optmethods/src/minpack/LevMar.hh',
                        'sherpa/optmethods/src/minpack/LevMar.cc']))

statfcts = Extension('sherpa.stats._statfcts',
              ['sherpa/stats/src/_statfcts.cc'],
              sherpa_inc,
              depends=get_deps(['stat_extension', 'stats']))

pykdtree = Extension('sherpa.utils._pykdtree',
              ['sherpa/utils/src/_pykdtree.cc'],
              sherpa_inc + ['sherpa/utils/src'],
              depends=(get_deps([]) +
                       ['sherpa/utils/src/kdtree++/allocator.hpp',
                        'sherpa/utils/src/kdtree++/function.hpp',
                        'sherpa/utils/src/kdtree++/iterator.hpp',
                        'sherpa/utils/src/kdtree++/kdtree.hpp',
                        'sherpa/utils/src/kdtree++/region.hpp',
                        'sherpa/utils/src/kdtree++/node.hpp']))

integration = Extension('sherpa.utils.integration',
              ['sherpa/utils/src/gsl/err.c',
               'sherpa/utils/src/gsl/error.c',
               'sherpa/utils/src/gsl/stream.c',
               'sherpa/utils/src/gsl/strerror.c',
               'sherpa/utils/src/gsl/message.c',
               'sherpa/utils/src/gsl/qng.c',
               'sherpa/utils/src/adapt_integrate.c',
               'sherpa/utils/src/integration.cc'],
              sherpa_inc + ['sherpa/utils/src',
                            'sherpa/utils/src/gsl'],
              depends=(get_deps(['integration'])+
                       ['sherpa/utils/src/adapt_integrate.h',
                        'sherpa/utils/src/gsl/gsl_integration.h']))

astro_modelfcts = Extension('sherpa.astro.models._modelfcts',
              ['sherpa/astro/models/src/_modelfcts.cc'],
              sherpa_inc,
              depends=get_deps(['model_extension', 'astro/models']))

pileup = Extension('sherpa.astro.utils._pileup',
              ['sherpa/astro/utils/src/fftn.c',
               'sherpa/astro/utils/src/_pileup.cc',
               'sherpa/astro/utils/src/pileup.cc'],
              sherpa_inc + ['sherpa/astro/utils/src'],
              depends=(get_deps(['extension']) +
                       ['sherpa/astro/utils/src/pileup.hh',
                        'sherpa/astro/utils/src/PyWrapper.hh',
                        'sherpa/astro/utils/src/fftn.inc']))

astro_utils = Extension('sherpa.astro.utils._utils',
              ['sherpa/astro/utils/src/_utils.cc'],
              (sherpa_inc + ['sherpa/utils/src/gsl']),
              libraries=(['sherpa']),
              depends=(get_deps(['extension', 'utils', 'astro/utils'])+
                       ['sherpa/utils/src/gsl/fcmp.h']))

####
# EXTENSIONS WITH EXTERNAL DEPENDENCIES
####

psf = Extension('sherpa.utils._psf',
              ['sherpa/utils/src/tcd/tcdCastArray.c',
               'sherpa/utils/src/tcd/tcdError.c',
               'sherpa/utils/src/tcd/tcdFFTConvolve.c',
               'sherpa/utils/src/tcd/tcdInitConvolveOut.c',
               'sherpa/utils/src/tcd/tcdInitTransform.c',
               'sherpa/utils/src/tcd/tcdPadData.c',
               'sherpa/utils/src/tcd/tcdPixelArith.c',
               'sherpa/utils/src/tcd/tcdTransform.c',
               'sherpa/utils/src/_psf.cc'],
              sherpa_inc + ['sherpa/utils/src/tcd', extern_inc],
              library_dirs=[extern_lib],
              libraries=['fftw3'],
              depends=(get_deps(['extension', 'utils'])+
                       ['sherpa/utils/src/tcd/tcd.h',]))

wcs = Extension('sherpa.astro.utils._wcs',
                  ['sherpa/astro/utils/src/_wcs.cc'],
                  sherpa_inc + [extern_inc],
                  library_dirs=[extern_lib],
                  libraries=['wcs'],
                  depends=get_deps(['extension']))

region = Extension('sherpa.astro.utils._region',
                  ['sherpa/astro/utils/src/_region.cc'],
                  sherpa_inc + [extern_inc],
                  library_dirs=[extern_lib],
                  libraries=(['region', 'wcs']),
                  depends=get_deps(['extension']))


####
# FORTRAN EXTENSIONS
####
minpack = Extension('sherpa.optmethods._minpack',
              ['sherpa/optmethods/src/minpack/_minpack.pyf',
               'sherpa/optmethods/src/minpack/covar.f',
               'sherpa/optmethods/src/minpack/lmdif.f',
               'sherpa/optmethods/src/minpack/mylmdif.f'])

minim =  Extension('sherpa.optmethods._minim',
              ['sherpa/optmethods/src/_minim.pyf',
               'sherpa/optmethods/src/minim.f',
               'sherpa/optmethods/src/syminv.f'])

meta['ext_modules'] = [estmethods,
               utils,
               modelfcts,
               saoopt,
               tstoptfct,
               statfcts,
               pykdtree,
               integration,
               astro_modelfcts,
               pileup,
               astro_utils,
               psf,
               wcs,
               region,
               minpack,
               minim,

               ]

meta['cmdclass'] = {
                    'build': build,
                    'clean' : clean,
                    'install' : install,
                    'sdist' : sdist,
                    }

setup(**meta)