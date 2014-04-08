#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

###
# Check that numpy is installed and with a good version
###
try:
    import imp
    imp.find_module('numpy')
except ImportError:
    import sys
    print >>sys.stderr, (
            "You need to install NUMPY in order to build Sherpa\n"
            "Other dependencies will be automatically installed\n"
            "Please install NUMPY (e.g. pip install numpy) and try again."
            )
    sys.exit(2)

import setuptools
from numpy.distutils.command.build import build as _build
from numpy.distutils.command.install import install as _install
from sdist import sdist as _sdist
from numpy.distutils.core import Command, Extension, setup as _setup
from distutils.command.clean import clean as _clean
from subprocess import call
from multiprocessing import cpu_count
import os


# Include directory for Sherpa headers
sherpa_inc = ['sherpa/include', 'sherpa/utils/src']

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
        next = deps.pop()
        if next not in alldeps:
            alldeps.add(next)
            deps.update(header_deps[next])

    return [sherpa_inc[0] + '/sherpa/' + d + '.hh' for d in alldeps]

####
# EXTENSIONS WITH EXTERNAL DEPENDENCIES
####

def build_psf_ext(library_dirs, include_dirs, libraries):
     return Extension('sherpa.utils._psf',
              ['sherpa/utils/src/tcd/tcdCastArray.c',
               'sherpa/utils/src/tcd/tcdError.c',
               'sherpa/utils/src/tcd/tcdFFTConvolve.c',
               'sherpa/utils/src/tcd/tcdInitConvolveOut.c',
               'sherpa/utils/src/tcd/tcdInitTransform.c',
               'sherpa/utils/src/tcd/tcdPadData.c',
               'sherpa/utils/src/tcd/tcdPixelArith.c',
               'sherpa/utils/src/tcd/tcdTransform.c',
               'sherpa/utils/src/_psf.cc'],
              sherpa_inc + ['sherpa/utils/src/tcd'] + include_dirs,
              library_dirs=library_dirs,
              libraries=libraries,
              depends=(get_deps(['extension', 'utils'])+
                       ['sherpa/utils/src/tcd/tcd.h',]))

def build_wcs_ext(library_dirs, include_dirs, libraries):
     return Extension('sherpa.astro.utils._wcs',
                  ['sherpa/astro/utils/src/_wcs.cc'],
                  sherpa_inc + include_dirs,
                  library_dirs=library_dirs,
                  libraries=libraries,
                  depends=get_deps(['extension']))

def build_region_ext(library_dirs, include_dirs, libraries):
     return Extension('sherpa.astro.utils._region',
                  ['sherpa/astro/utils/src/_region.cc'],
                  sherpa_inc + include_dirs,
                  library_dirs=library_dirs,
                  libraries=(libraries),
                  depends=get_deps(['extension']))

###
# Static Extensions
###

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
# FORTRAN EXTENSIONS
####
minpack = Extension('sherpa.optmethods._minpack',
              ['sherpa/optmethods/src/minpack/_minpack.pyf',
               'sherpa/optmethods/src/minpack/covar.f',
               'sherpa/optmethods/src/minpack/lmdif.f',
               'sherpa/optmethods/src/minpack/mylmdif.f',
               ])

minim =  Extension('sherpa.optmethods._minim',
              ['sherpa/optmethods/src/_minim.pyf',
               'sherpa/optmethods/src/minim.f',
               'sherpa/optmethods/src/syminv.f'])

####
### GROUP
####

group = Extension('group',
              ['extern/grplib-4.6/python/pygrplib.c'],
              ['extern/grplib-4.6/src'],
              library_dirs=['extern/grplib-4.6/src/.libs/'],
              libraries=['grp'],
              depends=['extern/grplib-4.6/python/pygrplib.h']
             )

def clean_sherpa():
    prefix = os.getcwd()
    os.chdir('extern')
    call(['make', 'uninstall'])
    call(['make', 'distclean'])
    try:
        os.remove('built')
    except:
        pass
    os.chdir(prefix)

def build_sherpa(configure):
    prefix=os.getcwd()
    os.chdir('extern')
    out = call(configure)
    if out != 0: exit(out)
    cflags = '-fPIC'
    out = call(['make', 'CFLAGS='+cflags,'-j'+str(cpu_count()+1), 'install'])
    if out != 0: exit(out)
    open('built', 'w').close()
    os.chdir(prefix)

class build(_build):
    def run(self):
        if not os.path.exists('extern/built'):
            configure = self.get_finalized_command('sherpa_config', True).build_configure()
            self.warn('built configure string' + str(configure))
            build_sherpa(configure)
        _build.run(self)

class install(_install):
    def run(self):
        if not os.path.exists('extern/built'):
            configure = self.get_finalized_command('sherpa_config', True).build_configure()
            self.warn('built configure string' + str(configure))
            build_sherpa(configure)
        _install.run(self)

class clean(_clean):
    def run(self):
        _clean.run(self)
        clean_sherpa()

class sdist(_sdist):
    def run(self):
        clean_sherpa()
        self.get_finalized_command('sherpa_config', True).build_configure()
        _sdist.run(self)



def setup(*args, **kwargs):

    kwargs['libraries'] = [

        ('sherpa',
         {'sources': ['sherpa/utils/src/gsl/fcmp.c'],
          'sourceDir' : 'sherpa/utils',
          'libs' : [],
          'libdirs' : [],
          'include_dirs': ['sherpa/utils/src'],
          'headerExportDir' : [],
          })
        ]



    ext_modules = [group,
                   estmethods,
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
                   minpack,
                   minim,
                   ]

    packages = kwargs['packages']
    package_data = kwargs['package_data']

    class sherpa_config(Command):
        description = "Configure Sherpa build options. If in doubt, ignore this command and stick to defaults. See setup.cfg for more information."
        user_options = [
                        ('fftw', None, "Whether Sherpa should build the embedded fftw3 library ('shipped', default)"),
                        ('fftw-include-dirs', None, "Where the fftw3 headers are located, if fftw is 'shipped'"),
                        ('fftw-lib-dirs', None, "Where the fftw3 libraries are located, if fftw is 'shipped'"),
                        ('fftw-libraries', None, "Name of the libraries that should be linked as fftw3"),
                        ('region', None, "Whether Sherpa should build the embedded region library ('shipped', default)"),
                        ('region-include-dirs', None, "Where the region headers are located, if region is 'shipped'"),
                        ('region-lib-dirs', None, "Where the region libraries are located, if region is 'shipped'"),
                        ('region-libraries', None, "Name of the libraries that should be linked as region"),
                        ('wcs-include-dirs', None, "Where the wcs subroutines headers are located"),
                        ('wcs-lib-dirs', None, "Where the wcs subroutines libraries are located"),
                        ('wcs-libraries', None, "Name of the libraries that should be linked as wcs"),
                        ('with-xspec', None, "Whether sherpa must build the XSPEC module (default False)")
                        ]

        def initialize_options(self):
            self.fftw=None
            self.fftw_include_dirs='build/include'
            self.fftw_lib_dirs='build/lib'
            self.fftw_libraries='fftw3'
            self.region=None
            self.region_include_dirs='build/include'
            self.region_lib_dirs='build/lib'
            self.region_libraries='region'
            self.wcs_include_dirs='build/include'
            self.wcs_lib_dirs='build/lib'
            self.wcs_libraries='wcs'
            self.with_xspec=False

        def finalize_options(self):
            pass

        def build_configure(self):
            configure = ['./configure','--disable-shared','--prefix='+os.getcwd()+'/build']
            if self.fftw != 'local':
                configure.append('--enable-fftw')
            ext_modules.append(build_psf_ext(*self._build_lib_arrays('fftw')))
            ext_modules.append(build_wcs_ext(*self._build_lib_arrays('wcs')))
            ld1, inc1, l1 = self._build_lib_arrays('wcs')
            if self.region != 'local':
                configure.append('--enable-region')
            ld2, inc2, l2 = self._build_lib_arrays('region')
            ld, inc, l = (ld1+ld2, inc1+inc2, l1+l2)
            ext_modules.append(build_region_ext(ld, inc, l))

            if self.with_xspec:
                packages.append('sherpa.astro.xspec')
                package_data['sherpa.astro.xspec': ['tests/test_*.py']]

            return configure

        def _build_lib_arrays(self, libname):
            library_dirs = getattr(self, libname+'_lib_dirs').split(' ')
            include_dirs = getattr(self, libname+'_include_dirs').split(' ')
            libraries = getattr(self, libname+'_libraries').split(' ')
            return [library_dirs, include_dirs, libraries]

    kwargs['ext_modules'] = ext_modules

    kwargs['cmdclass'] = {
                    'build': build,
                    'clean' : clean,
                    'install' : install,
                    'sdist' : sdist,
                    'sherpa_config' : sherpa_config,
                    }

#     kwargs['entry_points'] = {
#                              "distutils.setup.keywords": [
#                                                           "install_requires       = setuptools.dist:check_requirements",
#                              ]
#                              }

    return _setup(*args, **kwargs)
