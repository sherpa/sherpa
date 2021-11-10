#
# Copyright (C) 2014, 2017, 2018, 2019, 2020, 2021, 2022, 2023
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

# Ideally we would rely on the configuration in setup.cfg to avoid
# this situation, but as we may have users who still use old systems
# add an error message. As Sherpa releases occur the minimum-supported
# Python version is going to increase (beyond 3.7) but it is not really
# worth spending a lot of time on the error message here.
#
# This is done before we load in any non-core modules to avoid people
# installing software that they can not use.
#
if sys.version_info < (3, 9):
    sys.stderr.write("Sherpa 4.15.1 (and later) requires Python 3.9 or later.\n\n")
    sys.stderr.write("Please use Sherpa 4.15.0 if you need to use Python 3.8\n")
    sys.exit(1)

# We need to import setuptools so that 'python setup.py develop'
# works, but it isn't needed for 'pip install .'. Is this still true?
# I am leaving in as I can imagine it might require a newer setuptools
# than we currently support.
#
import setuptools
from numpy.distutils.core import setup

# How do we use local modules like helpers? This is a hack based on
# discussions around
# https://github.com/python-versioneer/python-versioneer/issues/193
#
sys.path.append(os.path.dirname(__file__))

from helpers import commands as sherpa_commands
from helpers.extensions import static_ext_modules

import versioneer

# First provide helpful messages if contributors try and run legacy commands.

TEST_HELP = """
Note: tests are no-longer run using 'python setup.py test'. Instead
you will need to use pytest. For example:

    pip install -e .[test]
    pytest

will ensure that pytest and pytest-xvfb are installed before running
the tests.
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)


setup(version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(sherpa_commands),
      ext_modules=static_ext_modules)
