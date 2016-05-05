#
# Copyright (C) 2015  Smithsonian Astrophysical Observatory
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
from sherpa.utils import _get_datadir
import pytest
import os

try:
    import mock as mock_module
except ImportError:
    from unittest import mock as mock_module


class _TestData():
    datadir = _get_datadir()

    skip_message = "required test data missing"

    def make_path(self, *segments):
        return os.path.join(self.datadir, *segments)


_testdata = _TestData()


@pytest.fixture(scope="function")
def testdata(request):
    if _testdata.datadir is None:
        pytest.skip(_testdata.skip_message)
    return _testdata


@pytest.fixture(scope="session")
def mock(request):
    return mock_module
