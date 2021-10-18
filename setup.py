#
# Copyright (C) 2014, 2017, 2018, 2019, 2020, 2021, 2022
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

# How do we use local modules like helpers? This is a hack based on
# discussions around
# https://github.com/python-versioneer/python-versioneer/issues/193
#
sys.path.append(os.path.dirname(__file__))

from helpers import commands as sherpa_commands
from helpers.extensions import static_ext_modules

import versioneer

commands = versioneer.get_cmdclass(sherpa_commands)

meta = dict(version=versioneer.get_version(),
            ext_modules=static_ext_modules,
            cmdclass=commands
            )

core.setup(**meta)
