#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

###############################################################################
#
# Handling of header interdependencies
#
###############################################################################

from numpy.distutils.command.build import build
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

class sherpa_build(build):
    def run(self):
        prefix=os.getcwd()
        os.chdir('extern')
        call(['./configure','--disable-shared','--prefix='+prefix+'/build'])
        call(['make', '-j'+str(cpu_count()+1), 'install'])
        os.chdir(prefix)
        build.run(self)

class sherpa_clean(clean):
    def run(self):
        clean.run(self)
        prefix = os.getcwd()
        os.chdir('extern')
        call(['make', 'distclean'])
        os.chdir(prefix)
