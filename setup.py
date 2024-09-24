#
# Copyright (C) 2014, 2017 - 2024
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

import os
import sys

# This is done before we load in any non-core modules to avoid people
# installing software that they can not use.
#
if sys.version_info < (3, 10):
    sys.stderr.write("Sherpa 4.17.0 (and later) requires Python 3.10 or later.\n\n")
    sys.stderr.write("Please use Sherpa 4.16.1 if you need to use Python 3.9\n")
    sys.exit(1)

# We need to import setuptools so that 'python setup.py develop'
# works, but it isn't needed for 'pip install .'. Is this still true?
# I am leaving in as I can imagine it might require a newer setuptools
# than we currently support.
#
from setuptools import setup

# How do we use local modules like helpers? This is a hack based on
# discussions around
# https://github.com/python-versioneer/python-versioneer/issues/193
#
sys.path.append(os.path.dirname(__file__))

from helpers import commands as sherpa_commands
from helpers.extensions import static_ext_modules

import versioneer

# First provide helpful messages if contributors try and run legacy commands.

HELP = {
    "test": """
Note: tests are no-longer run using 'python setup.py test'. Instead
you will need to use pytest. For example:

    pip install -e .[test]
    pytest

will ensure that pytest and pytest-xvfb are installed before running
the tests.
""",

    "develop": """
Note: 'python setup.py develop' is no-longer supported. Please use
either of

    pip install -e .
    pip install -e .[test]
""",

    "install": """
Note: 'python setup.py install' is no-longer supported. Please use

    pip install .
""",

    "build_sphinx": """
Note: the documentation is now built by saying

    cd docs
    make html
"""
}

for opt, help in HELP.items():
    if opt in sys.argv:
        print(help)
        sys.exit(1)


setup(version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(sherpa_commands),
      ext_modules=static_ext_modules)
