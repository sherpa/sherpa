#
#  Copyright (C) 2022
#  MIT
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
'''Pytest fixtures for running doctests

Pytest looks for conftest.py files that define fixtures only in the
directory tree that is being tested. When we run doctests in the
documentation only, it will therefore not find that fixtures that are
defined in sherpa/conftest.py, thus this file copies over the fixtures we
need.
'''
import pytest
from sherpa.utils.testing import get_datadir

@pytest.fixture(autouse=True)
def add_sherpa_test_data_dir(doctest_namespace):
    '''Define `data_dir` for doctests

    We make this an autouse=True fixture, because that means that
    the variable `data_dir` is available in every file that we may
    want to test without any special markup to that file. There is a
    risk that this is too magical and could confuse developers
    in the future, but it seems more important that this keeps any
    extra markup out of those files and keep them clean.
    '''
    # We could error out here, but we only want the tests that
    # use it to fail.
    #
    path = get_datadir()
    if path is None:
        pytest.skip("sherpa-test-data not found")

    if not path.endswith('/'):
        path += '/'

    doctest_namespace["data_dir"] = path
