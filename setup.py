#
#  Copyright (C) 2014, 2017, 2018, 2019, 2020, 2021, 2022
# Smithsonian Astrophysical Observatory
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

import sys

# python_requires will stop pip, but also let users who are using
# 'python setup.py develop'.
#
if sys.version_info < (3, 7):
    sys.stderr.write("Sherpa 4.14 (and later) requires Python 3.7 or later.\n\n")
    sys.stderr.write("Please use Sherpa 4.13.1 if you need to use Python 3.6\n")
    sys.exit(1)

try:
    import setuptools
except:
    print((
        "WARNING\n"
        "Could not import setuptools.\n"
        "This might lead to an incomplete installation\n"
    ), file=sys.stderr)

# The module used to try to find the numpy module and error out if
# that could not be found, but it seems simpler to just error out
# here.
#
try:
    from numpy.distutils import core

except ImportError:
    print((
        "You need to install NUMPY in order to build Sherpa\n"
        "Other dependencies will be automatically installed\n"
        "Please install NUMPY (e.g. pip install numpy) and try again."
    ), file=sys.stderr)
    sys.exit(2)

from helpers import commands as sherpa_commands
from helpers.extensions import static_ext_modules

import versioneer

commands = versioneer.get_cmdclass(sherpa_commands)

meta = dict(name='sherpa',
            version=versioneer.get_version(),
            author='Smithsonian Astrophysical Observatory / Chandra X-Ray Center',
            author_email='cxchelp@head.cfa.harvard.edu',
            url='http://cxc.harvard.edu/sherpa/',
            description='Modeling and fitting package for scientific data analysis',
            license='GNU GPL v3',
            long_description=open('README.md', 'rt').read(),
            long_description_content_type='text/markdown',
            platforms='Linux, Mac OS X',
            python_requires='~=3.7',
            install_requires=['numpy'],
            tests_require=['pytest-xvfb', 'pytest>=5.0,!=5.2.3'],
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
            package_data={# Ordered to match 'find sherpa -name tests | sort'
                          'sherpa.astro.datastack': ['tests/data/*', 'tests/*.py'],
                          'sherpa.astro.io': ['tests/test_*.py'],
                          'sherpa.astro.models': ['tests/test_*.py'],
                          'sherpa.astro.optical': ['tests/test_*.py'],
                          'sherpa.astro.sim': ['tests/test_*.py'],
                          'sherpa.astro': ['tests/test_*.py'],
                          'sherpa.astro.ui': ['tests/data/*', 'tests/test_*.py'],
                          'sherpa.astro.utils': ['tests/test_*.py'],
                          'sherpa.astro.xspec': ['tests/test_*.py'],
                          'sherpa.estmethods': ['tests/test_*.py'],
                          'sherpa.image': ['tests/test_*.py'],
                          'sherpa.models': ['tests/test_*.py'],
                          'sherpa.optmethods': ['tests/*.hh', 'tests/_tstoptfct.cc', 'tests/test_*.py'],
                          'sherpa.plot': ['tests/test_*.py'],
                          'sherpa.sim': ['tests/test_*.py'],
                          'sherpa.stats': ['tests/test_*.py'],
                          'sherpa': ['include/sherpa/*.hh',
                                     'include/sherpa/astro/*.hh',
                                     'tests/*',
                                     'static/css/*'],
                          'sherpa.ui': ['tests/test_*.py'],
                          'sherpa.utils': ['tests/test_*.py'],
                          },
            data_files=[('sherpa',
                         ['sherpa/sherpa.rc', 'sherpa/sherpa-standalone.rc']), ],
            ext_modules=static_ext_modules,
            cmdclass=commands,
            entry_points={
                'console_scripts': [
                    'sherpa_test = sherpa:clitest',
                    'sherpa_smoke = sherpa:_smoke_cli',
                ]
            },
            classifiers=[
                'Intended Audience :: Science/Research',
                'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                'Programming Language :: C',
                'Programming Language :: Python :: 3.7',
                'Programming Language :: Python :: 3.8',
                'Programming Language :: Python :: 3.9',
                'Programming Language :: Python :: 3.10',
                'Programming Language :: Python :: Implementation :: CPython',
                'Topic :: Scientific/Engineering :: Astronomy',
                'Topic :: Scientific/Engineering :: Physics'
                ],
            )

core.setup(**meta)
