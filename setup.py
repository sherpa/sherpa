# 
#  Copyright (C) 2014  Smithsonian Astrophysical Observatory
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


# ##
# Check that numpy is installed and with a good version
###
try:
    import imp

    imp.find_module('numpy')
except ImportError:
    import sys

    print >> sys.stderr, (
        "You need to install NUMPY in order to build Sherpa\n"
        "Other dependencies will be automatically installed\n"
        "Please install NUMPY (e.g. pip install numpy) and try again."
    )
    sys.exit(2)

import os

try:
    import setuptools
except:
    import sys

    print >> sys.stderr, (
        "WARNING\n"
        "Could not import setuptools.\n"
        "This might lead to an incomplete installation\n"
    )
from numpy.distutils.core import setup

from helpers.extensions import static_ext_modules
#from helpers import commands

import versioneer

versioneer.VCS='git'
versioneer.versionfile_source = 'sherpa/_version.py'
versioneer.versionfile_build = 'sherpa/_version.py'
versioneer.tag_prefix = ''
versioneer.parentdir_prefix = 'sherpa-'

meta = dict(name='sherpa',
            version=versioneer.get_version(),
            author='Smithsonian Astrophysical Observatory / Chandra X-Ray Center',
            author_email='cxchelp@head.cfa.harvard.edu',
            url='http://cxc.harvard.edu/sherpa/',
            description='Modeling and fitting package for scientific data analysis',
            license='GNU GPL v3',
            long_description=open('README.md', 'rt').read(),
            platforms='Linux, Mac OS X',
            install_requires=['numpy',],
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
                      'sherpa.astro.datastack',
                      'sherpa.astro.io',
                      'sherpa.astro.models',
                      'sherpa.astro.optical',
                      'sherpa.astro.sim',
                      'sherpa.astro.ui',
                      'sherpa.astro.utils',
                      'sherpa.astro.xspec',
            ],
            package_data={'sherpa': ['include/sherpa/*.hh',
                                     'include/sherpa/astro/*.hh',
                                     'tests/*'],
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
                          'sherpa.astro.datastack': ['tests/data/*', 'tests/*.py'],
                          'sherpa.astro.io': ['tests/test_*.py'],
                          'sherpa.astro.models': ['tests/test_*.py'],
                          'sherpa.astro.optical': ['tests/test_*.py'],
                          'sherpa.astro.sim': ['tests/test_*.py'],
                          'sherpa.astro.ui': ['tests/test_*.py'],
                          'sherpa.astro.utils': ['tests/test_*.py'],
                          },
            data_files=[('sherpa',
		    ['sherpa/sherpa.rc', 'sherpa/sherpa-standalone.rc']), ],
            ext_modules=static_ext_modules, cmdclass=versioneer.get_cmdclass(),
            entry_points={
                'console_scripts': [
                    'sherpa_test = sherpa:clitest',
                ]
            },
            classifiers=[
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                'Programming Language :: C',
                'Programming Language :: Fortran',
                'Programming Language :: Python :: 2.7',
                'Programming Language :: Python :: Implementation :: CPython',
                'Topic :: Scientific/Engineering :: Astronomy',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            )

setup(**meta)

