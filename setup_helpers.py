#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

###############################################################################
#
# Handling of header interdependencies
#
###############################################################################

from numpy.distutils.command.build import build
from numpy.distutils.command.install import install
from distutils.command.clean import clean
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
    call(['make', 'CFLAGS=-fPIC','-j'+str(cpu_count()+1), 'install'])
    os.chdir(prefix)

def clean_sherpa():
    prefix = os.getcwd()
    os.chdir('extern')
    call(['make', 'distclean'])
    os.chdir(prefix)

class sherpa_build(build):
    def run(self):
        build_sherpa()
        build.run(self)

class sherpa_clean(clean):
    def run(self):
        clean.run(self)
        clean_sherpa()

class sherpa_install(install):
    def run(self):
        build_sherpa()
        install.run(self)
