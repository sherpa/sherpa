#
#  Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

"""Access to the Sherpa test data directory.

This is used when running Sherpa tests.
"""

import logging
import os

import pkg_resources


warning = logging.getLogger(__name__).warning


# This code has been moved out of sherpa.utils.testing to make
# it easy to import (sherpa.utils requires sherpa.utils._utils
# which seems to cause problems in certain situations).
#


# I am not convinced the conversion from __file__ access to
# pkg_resources is correct here. Using __file__ is frowned upon,
# hence the change, and as we can not assume Python 3.7 or later
# at the moment we can not use importlib_resources.
#
def _get_datadir():
    """Setup the location of the Sherpa test data, if installed."""

    try:
        import sherpatest
        datadir = pkg_resources.resource_filename('sherpatest', '')

    except ImportError:
        try:
            import sherpa
            datadir = pkg_resources.resource_filename('sherpa', '')
            datadir = os.path.join(datadir, os.pardir,
                                   'sherpa-test-data', 'sherpatest')
        except ImportError:
            # neither sherpatest nor sherpa can be found, falling back to None
            datadir = None

    # Check the directory exists
    if datadir is None or not os.path.exists(datadir) or not os.listdir(datadir):
        return None

    return datadir


DATADIR = _get_datadir()


def get_datadir():
    """Return the location of the Sherpa test data.

    The data directory is determined either by the existance of the
    sherpatest module, sherpa-test-directory submodule, or
    by being set manually with set_datadir.
    """

    return DATADIR


def set_datadir(dname):
    """Set the location of the Sherpa test data.

    """

    from sherpa.utils.testing import SherpaTestCase

    global DATADIR
    if not os.path.exists(dname) or not os.listdir(dname):
        warning("Unable to set the datadir to {}".format(dname))
        dname = None

    DATADIR = dname
    SherpaTestCase.datadir = dname
