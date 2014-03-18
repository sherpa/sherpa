#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

import sys
import os
import platform
from numpy.distutils.core import setup, Extension

###############################################################################
#
# Local configuration
#
###############################################################################


#
# Set defaults
#

conf = {
    'cfitsio_lib': 'cfitsio',
    'cfitsio_library_dir': None,
    'xspec_library_dir': None,
    'reg_library_dir': None,
    'reg_include_dir': None,
    'fftw_library_dir' : '/usr/lib',
    'fftw_include_dir' : '/usr/include',
    'wcs_library_dir' : None,
    'wcs_include_dir' : None,
    'fortran_lib' : None,
    'fortran_library_dir' : None
    }

#
# To store f2py temp files somewhere they can be easily removed,
# create the "build" directory if it doesn't exist, and set the
# tempfile module to use "./build" as the directory to store
# temp files.  (This covers request not to leave temp files in
# /tmp, where they would otherwise just accumulate over the
# course of many builds.)
#

import tempfile
tempfile.tempdir='./build'
try:
    os.mkdir("./build",0775)
except OSError, e:
    # This error means the directory already exists,
    # so pass
    if (e.errno == 17):
        pass
    else:
        raise e

#
# Process 'var=value' command-line arguments, replacing the defaults with
# the specified values
#

while ((len(sys.argv) > 1) and (sys.argv[1].count('=') == 1) and
       (not sys.argv[1].startswith('--'))):
    var, val = sys.argv.pop(1).split('=')
    
    if var not in conf:
        sys.stderr.write("error: '%s' is not a valid configuration variable\n"
                         % var)
        sys.exit(1)

    conf[var] = val


# The Fortran compiler class respects environment variables
# FDEBUG and FOPT, to set the debug and optimization flags respectively
# So no need to examine those env. variabales here, the relevant
# FCompiler class will pick them up.
#
# For the C and C++ compilers, we DO have to a bit of work here.
# distutils.sysconfig._config_vars will stick in spurious -g and
# -O3 flags, even if you want to set these flags for yourself!
#
# So, here I will examine if environment variables DFLAGS or OFLAGS
# have been set.  These are values set in CIAO Makefiles; but even if
# you call this setup file directly, you can set these env. variables
# by hand, in your shell, and that will govern whether you use -g or
# -O*here (or -xO* on Solaris).
#
# If env. variables DFLAGS or OCFLAGS doesn't exist, then defaults from
# the C and C++ compiler classes will be used.  SMD 10/30/2008
if (os.environ.has_key('DFLAGS') == True):
    if (os.environ['DFLAGS'] == ''):
        import distutils.sysconfig
        old_str = distutils.sysconfig._config_vars['CFLAGS']
        distutils.sysconfig._config_vars['CFLAGS'] = old_str.replace('-g','')
        old_str = distutils.sysconfig._config_vars['OPT']
        distutils.sysconfig._config_vars['OPT'] = old_str.replace('-g','')

if (os.environ.has_key('OCFLAGS') == True):
    oflag = os.environ['OCFLAGS']
    import distutils.sysconfig
    old_str = distutils.sysconfig._config_vars['CFLAGS']
    distutils.sysconfig._config_vars['CFLAGS'] = old_str.replace('-O ','').replace('-O0','').replace('-O1','').replace('-O2','').replace('-O3','').replace('-O4','').replace('-O5','') + ' ' + oflag
    old_str = distutils.sysconfig._config_vars['OPT']
    distutils.sysconfig._config_vars['OPT'] = old_str.replace('-O ','').replace('-O0','').replace('-O1','').replace('-O2','').replace('-O3','').replace('-O4','').replace('-O5','') + ' ' + oflag


###############################################################################
#
# Platform-specific build options
#
###############################################################################
needf90 = False

#
# Link against libCstd when building C++ extensions on Solaris
#

cpp_libs = []

if platform.system() == 'SunOS':
    cpp_libs.extend(['Cstd', 'sunmath'])


#
# Link against libfui and libmtsk when building with Sun's Fortran compiler
#

from numpy.distutils.fcompiler.sun import SunFCompiler
SunFCompiler.old_get_libraries = SunFCompiler.get_libraries

def get_libraries(self):
    libs = self.old_get_libraries()
    for l in ['fui', 'mtsk']:
        if l not in libs:
            libs.append(l)
    return libs

SunFCompiler.get_libraries = get_libraries

#
# On Darwin, g77 is not in /usr/bin; the user needs to have
# installed it in /sw/bin with fink.  But then the GnuFCompiler
# class needs to have /sw/lib added as a library path, so we
# can get the g2c library.
#
if platform.system() == 'Darwin':
    needf90 = True 
    from numpy.distutils.fcompiler.gnu import GnuFCompiler
    GnuFCompiler.old_get_libraries = GnuFCompiler.get_libraries
    GnuFCompiler.old_get_library_dirs = GnuFCompiler.get_library_dirs

    # Have also found that on OS X, keeping cc_dynamic on the
    # link line prevents linkage--fortunately it appears that
    # none of our code requires cc_dynamic anyway.
    def get_libraries_gnuf(self):
        from numpy.distutils.fcompiler.gnu import Gnu95FCompiler
        libs = self.old_get_libraries()
        if ('cc_dynamic' in libs and
            isinstance(self, Gnu95FCompiler) == False):
            libs.remove('cc_dynamic')
        return libs

    def get_libraries_g95(self):
        libs = self.old_get_libraries()
        new_libs = ['SystemStubs']
        if (conf['fortran_lib'] != None):
            new_libs = [conf['fortran_lib'], 'SystemStubs']
        for l in new_libs:
            if l not in libs:
                libs.append(l)
        return libs

    # Need to add path to library for g2c, even though
    # g2c itself *is* already in list of libraries
    def get_library_dirs(self):
        dirs = self.old_get_library_dirs()
        if (conf['fortran_library_dir'] != None and
            conf['fortran_library_dir'] != './'):
            if conf['fortran_library_dir'] not in dirs:
                dirs.append(conf['fortran_library_dir'])
        else:
            # Try likely paths if fortran_library_dir points nowhere
            if '/sw/lib' not in dirs:
                dirs.append('/sw/lib')
            if '/opt/local/lib' not in dirs:
                dirs.append('/opt/local/lib') 
        return dirs

    GnuFCompiler.get_libraries = get_libraries_gnuf
    GnuFCompiler.get_library_dirs = get_library_dirs

    # Block gfortran from adding superfluous -arch flags
    def _universal_flags(self, cmd):
        # These need to come from CCEXTRA_ARGS, without setting
        # LDFLAGS
        return ['-m64']

    from numpy.distutils.fcompiler.gnu import Gnu95FCompiler
    Gnu95FCompiler._universal_flags = _universal_flags


    # If on the Intel Mac, must use g95 to compile Fortran
    # Use gcc as linker, because g95 just isn't helping us make
    # bundles on OS X
    from numpy.distutils.fcompiler.g95 import G95FCompiler
    G95FCompiler.old_get_library_dirs = G95FCompiler.get_library_dirs
    G95FCompiler.old_get_libraries = G95FCompiler.get_libraries
    G95FCompiler.get_library_dirs = get_library_dirs
    G95FCompiler.get_libraries = get_libraries_g95
    G95FCompiler.executables['linker_so'] = ["gcc","-undefined dynamic_lookup -bundle"]
else:
    # On any other platform, still need to modify G95Compiler to
    # add f95 to library path, for XSPEC module build
    def get_libraries_g95(self):
        libs = self.old_get_libraries()
        if (conf['fortran_lib'] != None):
            if conf['fortran_lib'] not in libs:
                libs.append(conf['fortran_lib'])
        return libs

    def get_library_dirs(self):
        dirs = self.old_get_library_dirs()
        if (conf['fortran_library_dir'] != None and
            conf['fortran_library_dir'] != './'):
            if conf['fortran_library_dir'] not in dirs:
                dirs.append(conf['fortran_library_dir'])
        return dirs

    from numpy.distutils.fcompiler.g95 import G95FCompiler
    G95FCompiler.old_get_libraries = G95FCompiler.get_libraries
    G95FCompiler.get_libraries = get_libraries_g95
    G95FCompiler.old_get_library_dirs = G95FCompiler.get_library_dirs
    G95FCompiler.get_library_dirs = get_library_dirs

###############################################################################
#
# Handling of header interdependencies
#
###############################################################################

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



###############################################################################
#
# Extension modules
#
###############################################################################

    
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

#
# Standard modules
#

extension_modules = [
    # sherpa.estmethods._est_funcs
    Extension('sherpa.estmethods._est_funcs',
              ['sherpa/estmethods/src/estutils.cc',
               'sherpa/estmethods/src/info_matrix.cc',
               'sherpa/estmethods/src/projection.cc',
               'sherpa/estmethods/src/estwrappers.cc'],
              (sherpa_inc + ['sherpa/utils/src/gsl']),
              libraries=(cpp_libs + ['sherpa']),
              depends=(get_deps(['extension', 'utils']) +
                       ['sherpa/estmethods/src/estutils.hh',
                        'sherpa/estmethods/src/info_matrix.hh',
                        'sherpa/estmethods/src/projection.hh',
                        'sherpa/utils/src/gsl/fcmp.h'])),

    # sherpa.models._modelfcts
    Extension('sherpa.models._modelfcts',
              ['sherpa/models/src/_modelfcts.cc'],
              sherpa_inc,
              libraries=cpp_libs,
              depends=get_deps(['model_extension', 'models'])),

##################################optmethods###################################
    
    # sherpa.optmethods._minpack
    Extension('sherpa.optmethods._minpack',
              ['sherpa/optmethods/src/minpack/_minpack.pyf',
               'sherpa/optmethods/src/minpack/covar.f',               
               'sherpa/optmethods/src/minpack/lmdif.f',
               'sherpa/optmethods/src/minpack/mylmdif.f']),

    # simplex <==> neldermead
    # leave minim here for debug purpose only
    # sherpa.optmethods._minim
    Extension('sherpa.optmethods._minim',
              ['sherpa/optmethods/src/_minim.pyf',
               'sherpa/optmethods/src/minim.f',
               'sherpa/optmethods/src/syminv.f']),

###############################################################################

##    # sherpa.optmethods._chokkan
##     Extension('sherpa.optmethods._chokkan',
##               ['sherpa/optmethods/src/chokkan/_chokkan.cc',
##                'sherpa/optmethods/src/chokkan/lbfgs.c'],
##               sherpa_inc + ['sherpa/optmethods/src/chokkan'],
##               depends=(get_deps(['myArray', 'extension']) +
##                        ['sherpa/optmethods/src/chokkan/lbfgs.h'])),

##     # powell: bobyqa & newuoa 
##     Extension('sherpa.optmethods._powell',
##               ['sherpa/optmethods/src/powell/_powell.pyf',
##                'sherpa/optmethods/src/powell/altmov.f',
##                'sherpa/optmethods/src/powell/bigden.f',
##                'sherpa/optmethods/src/powell/biglag.f',
##                'sherpa/optmethods/src/powell/bobyqa.f',
##                'sherpa/optmethods/src/powell/bobyqb.f',
##                'sherpa/optmethods/src/powell/bupdate.f',
##                'sherpa/optmethods/src/powell/newuoa.f',
##                'sherpa/optmethods/src/powell/newuob.f',
##                'sherpa/optmethods/src/powell/prelim.f',
##                'sherpa/optmethods/src/powell/rescue.f',
##                'sherpa/optmethods/src/powell/trsapp.f',
##                'sherpa/optmethods/src/powell/trsbox.f',
##                'sherpa/optmethods/src/powell/update.f']),

##     # sherpa.optmethods._odrpack
##     Extension('sherpa.optmethods._odrpack',
##               ['sherpa/optmethods/src/odrpack/_odrpack.pyf',
##                'sherpa/optmethods/src/odrpack/real_precision.f95',
##                'sherpa/optmethods/src/odrpack/myodr.f95',
##                'sherpa/optmethods/src/odrpack/lpkbls.f95',
##                'sherpa/optmethods/src/odrpack/odr.f95']),

##     # sherpa.optmethods._port
##     Extension('sherpa.optmethods._port',
##               ['sherpa/optmethods/src/port/_port.pyf',
##                'sherpa/optmethods/src/port/port.f',
##                'sherpa/optmethods/src/port/myport.f',
##                'sherpa/optmethods/src/port/dpptri.f']),

##    # sherpa.optmethods._stogo
##     Extension('sherpa.optmethods._stogo',
##              ['sherpa/optmethods/src/StoGo/_stogo.cc',
##               'sherpa/optmethods/src/StoGo/linalg.cc',
##               'sherpa/optmethods/src/StoGo/tools.cc'],
##               depends=(get_deps(['myArray', 'extension']) +
##                        ['sherpa/include/sherpa/functor.hh',
##                        'sherpa/optmethods/src/StoGo/linag.h',
##                        'sherpa/optmethods/src/StoGo/tools.h',
##                        'sherpa/optmethods/src/StoGo/Local.hh',
##                        'sherpa/optmethods/src/StoGo/Global.hh'])),

 
###############################################################################

    # sherpa.optmethods._saoopt
    Extension('sherpa.optmethods._saoopt',
              ['sherpa/optmethods/src/_saoopt.cc',
               'sherpa/optmethods/src/Simplex.cc'],
              sherpa_inc + ['sherpa/utils/src/gsl'],
              libraries=(cpp_libs + ['sherpa']),              
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
                        'sherpa/optmethods/src/minpack/LevMar.cc'])),

    # sherpa.optmethods._tstoptfct
    Extension('sherpa.optmethods._tstoptfct',
              ['sherpa/optmethods/tests/_tstoptfct.cc'],
              sherpa_inc,
              libraries=cpp_libs,
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
                        'sherpa/optmethods/src/minpack/LevMar.cc'])),

##################################optmethods###################################

    # sherpa.stats._statfcts
    Extension('sherpa.stats._statfcts',
              ['sherpa/stats/src/_statfcts.cc'],
              sherpa_inc,
              libraries=cpp_libs,
              depends=get_deps(['stat_extension', 'stats'])),

    # sherpa.utils._utils
    Extension('sherpa.utils._utils',
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
              libraries=(cpp_libs + ['sherpa']),
              depends=(get_deps(['extension', 'utils'])+
                       ['sherpa/utils/src/gsl/fcmp.h',
                        'sherpa/utils/src/cephes/cephes.h'])),


    # sherpa.utils._pykdtree
    Extension('sherpa.utils._pykdtree',
              ['sherpa/utils/src/_pykdtree.cc'],
              sherpa_inc + ['sherpa/utils/src'],
              libraries=cpp_libs,
              depends=(get_deps([]) +
                       ['sherpa/utils/src/kdtree++/allocator.hpp',
                        'sherpa/utils/src/kdtree++/function.hpp',
                        'sherpa/utils/src/kdtree++/iterator.hpp',
                        'sherpa/utils/src/kdtree++/kdtree.hpp',
                        'sherpa/utils/src/kdtree++/region.hpp',
                        'sherpa/utils/src/kdtree++/node.hpp'])),

    # sherpa.utils._psf
    Extension('sherpa.utils._psf',
              ['sherpa/utils/src/tcd/tcdCastArray.c',
               'sherpa/utils/src/tcd/tcdError.c',
               'sherpa/utils/src/tcd/tcdFFTConvolve.c',
               'sherpa/utils/src/tcd/tcdInitConvolveOut.c',
               'sherpa/utils/src/tcd/tcdInitTransform.c',
               'sherpa/utils/src/tcd/tcdPadData.c',
               'sherpa/utils/src/tcd/tcdPixelArith.c',
               'sherpa/utils/src/tcd/tcdTransform.c',
               'sherpa/utils/src/_psf.cc'],
              sherpa_inc + ['sherpa/utils/src/tcd', conf['fftw_include_dir']],
              library_dirs=[conf['fftw_library_dir']],
              libraries=(cpp_libs + ['fftw3']),
              depends=(get_deps(['extension', 'utils'])+
                       ['sherpa/utils/src/tcd/tcd.h'])),
    
    # sherpa.utils.integration
    Extension('sherpa.utils.integration',
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
              libraries=cpp_libs,
              depends=(get_deps(['integration'])+
                       ['sherpa/utils/src/adapt_integrate.h',
                        'sherpa/utils/src/gsl/gsl_integration.h'])),

    # sherpa.astro.models._modelfcts
    Extension('sherpa.astro.models._modelfcts',
              ['sherpa/astro/models/src/_modelfcts.cc'],
              sherpa_inc,
              libraries=cpp_libs,
              depends=get_deps(['model_extension', 'astro/models'])),

    # sherpa.astro.utils._pileup
    Extension('sherpa.astro.utils._pileup',
              ['sherpa/astro/utils/src/fftn.c',
               'sherpa/astro/utils/src/_pileup.cc',
               'sherpa/astro/utils/src/pileup.cc'],
              sherpa_inc + ['sherpa/astro/utils/src'],
              libraries=cpp_libs,
              depends=(get_deps(['extension']) +
                       ['sherpa/astro/utils/src/pileup.hh',
                        'sherpa/astro/utils/src/PyWrapper.hh',
                        'sherpa/astro/utils/src/fftn.inc'])),

    # sherpa.astro.utils._utils
    Extension('sherpa.astro.utils._utils',
              ['sherpa/astro/utils/src/_utils.cc'],
              (sherpa_inc + ['sherpa/utils/src/gsl']),
              libraries=(cpp_libs + ['sherpa']),
              depends=(get_deps(['extension', 'utils', 'astro/utils'])+
                       ['sherpa/utils/src/gsl/fcmp.h'])),
    
    ]

#
# WCS module (optional)
#

if (conf['wcs_library_dir'] is not None and
    conf['wcs_include_dir'] is not None):
    extension_modules.append(
        # sherpa.astro.utils._wcs
        Extension('sherpa.astro.utils._wcs',
                  ['sherpa/astro/utils/src/_wcs.cc'],
                  sherpa_inc + [conf['wcs_include_dir']],
                  library_dirs=[conf['wcs_library_dir']],
                  libraries=(cpp_libs + ['wcs']),
                  depends=get_deps(['extension'])),
        
        )

#
# Region module (optional)
#

if (conf['reg_library_dir'] is not None and
    conf['reg_include_dir'] is not None and 
    conf['wcs_library_dir'] is not None and
    conf['wcs_include_dir'] is not None):

    reg_lib_dirs = [conf['reg_library_dir'], conf['wcs_library_dir']]
    if conf['cfitsio_library_dir'] is not None:
        reg_lib_dirs.append(conf['cfitsio_library_dir'])

    extension_modules.append(
        # sherpa.astro.utils._region
        Extension('sherpa.astro.utils._region',
                  ['sherpa/astro/utils/src/_region.cc'],
                  sherpa_inc + [conf['reg_include_dir'], conf['wcs_include_dir']],
                  library_dirs=reg_lib_dirs,
                  libraries=(cpp_libs +
                             ['region', 'ascdm', conf['cfitsio_lib'], 'wcs']),
                  depends=get_deps(['extension'])),
        )

#
# XSPEC module (optional)
#

if conf['xspec_library_dir'] is not None:
    import numpy.distutils.fcompiler
    fc = numpy.distutils.fcompiler.new_fcompiler(requiref90=needf90)
    fc.customize()
    
    xspec_libs = (cpp_libs +
                  ['XSFunctions', 'XSModel', 'XSUtil', 'XS',
                   'CCfits', conf['cfitsio_lib']] +
                  fc.get_libraries())

    xspec_library_dirs = (fc.get_library_dirs() +
                          [conf['xspec_library_dir']])
        
    if conf['cfitsio_library_dir'] is not None:
        xspec_library_dirs.append(conf['cfitsio_library_dir'])

    extension_modules.append(
        # sherpa.astro.xspec._xspec
        Extension('sherpa.astro.xspec._xspec',
                  ['sherpa/astro/xspec/src/_xspec.cc'],
                  sherpa_inc,
                  library_dirs=xspec_library_dirs,
                  runtime_library_dirs=xspec_library_dirs,
                  libraries=xspec_libs,
                  depends=(get_deps(['astro/xspec_extension'])))
        )



###############################################################################
#
# Run setup
#
###############################################################################

# CIAO 4.6 release, Sherpa package 1
setup(name='sherpa',
      version='4.6.1',
      author='Smithsonian Astrophysical Observatory / Chandra X-Ray Center',
      author_email='cxchelp@head.cfa.harvard.edu',
      url='http://cxc.harvard.edu/sherpa/',
      description='Modeling and fitting package for scientific data analysis',
      license='GNU GPL v3',
      long_description='Modeling and fitting package for scientific data analysis',
      platforms='Linux, Mac OS X, Solaris',
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
                'sherpa.astro.xspec'],
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
                    'sherpa.astro.xspec': ['tests/test_*.py']},
      libraries=clibs,
      ext_modules=extension_modules,
      data_files=[('sherpa', ['sherpa/sherpa.rc']),
                  ('sherpa', ['sherpa/00-sherpa_startup.py']),
                  ('sherpa', ['sherpa/ipython_config.py'])],
      )
