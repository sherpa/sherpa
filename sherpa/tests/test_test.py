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
from sherpa import conftest
import pytest

try:
    import mock
except ImportError:
    from unittest import mock


@pytest.fixture
def datadir_test(request):
    """
    inject test path
    """
    conftest._testdata.datadir = "test_path"


@pytest.fixture
def datadir_none(request):
    """
    simulate scenario where datadir is None
    """
    datadir = conftest._testdata.datadir
    msg = conftest._testdata.skip_message

    def fin():
        conftest._testdata.datadir = datadir
        conftest._testdata.skip_message = msg
    request.addfinalizer(fin)

    conftest._testdata.datadir = None
    conftest._testdata.skip_message = "testing tests get skipped when no data dir is available"


def test_testdata(datadir_test, testdata):
    assert testdata.datadir == "test_path"


def test_make_path(datadir_test, testdata):
    assert testdata.make_path("one", "two") == "test_path/one/two"


def test_make_path_none(datadir_none, testdata):
    # This test should not even be run
    assert testdata.make_path("one", "two") == "foo"


def test_skip(datadir_none, testdata):
    assert False  # This test should not even be run
