#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

###############################################################################
#
# Handling of header interdependencies
#
###############################################################################

from numpy.distutils.command.build import build as _build
from numpy.distutils.command.install import install as _install
from numpy.distutils.command.sdist import sdist as _sdist
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

def build_sherpa():
    prefix=os.getcwd()
    os.chdir('extern')
    call(['./configure','--disable-shared','--prefix='+prefix+'/build'])
    out = call(['make', 'CFLAGS=-fPIC','-j'+str(cpu_count()+1), 'install'])
    if out != 0: exit(out)
    open('built', 'w').close()
    os.chdir(prefix)

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

class build(_build):
    def run(self):
        if not os.path.exists('extern/built'):
            build_sherpa()
        _build.run(self)

class clean(_clean):
    def run(self):
        _clean.run(self)
        clean_sherpa()

class install(_install):
    def run(self):
        if not os.path.exists('extern/built'):
            build_sherpa()
        _install.run(self)

class sdist(_sdist):
    def run(self):
        clean_sherpa()
        _sdist.run(self)
