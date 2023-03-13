#
#  Copyright (C) 2017, 2018, 2019, 2023
#  Smithsonian Astrophysical Observatory
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

import importlib
import sherpa
import pytest

SKIP_MESSAGE = "This test should have been skipped"


def test_decorator_exception_when_pytest_is_absent(monkeypatch):
    """
    Test that when pytest is not installed (simulated with monkeypatch) the decorators raise an exception when
    a test is run.
    """
    import sys
    monkeypatch.setitem(sys.modules, 'pytest', None)
    importlib.reload(sherpa.utils.testing)

    from sherpa.utils import testing

    @testing.requires_data
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_pylab
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_fits
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_group
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_stk
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_ds9
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_xspec
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)

    @testing.requires_package("non.existent.package")
    def foo():
        raise RuntimeError(SKIP_MESSAGE)
    assert_decorated_function(foo)


def assert_decorated_function(function):
    """
    Pytest cannot be assumed to be installed by the regular user, unlike unittest, which is part of Python's
    standard library. The tests ensures that the decorators can still be used but if the tests are run they throw
    and exception prompting users to install pytest, in those cases where pytest is not installed automatically.
    """

    with pytest.raises(ImportError) as excinfo:
        function()

    assert sherpa.utils.testing.PYTEST_MISSING_MESSAGE == str(excinfo.value)
