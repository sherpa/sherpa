#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

import pytest

# Whilelist of known warnings. One can associate different warning messages to the same warning class
known_warnings = {
    DeprecationWarning: ['unorderable dtypes; returning scalar but in the future this will be an error',
                         "Non-string object detected for the array ordering. Please pass in 'C', 'F', 'A', or 'K' instead"]
}


@pytest.fixture(scope="function", autouse=True)
def capture_all_warnings(request, recwarn):
    """
    This fixture will run automatically before and after every test function is executed.
    It uses pytest's infrastructure to get all recorded warnings and match them against the while list. If an
    unknown warning is found, then the fixture finalizer will fail the specific test function.

    In the verbose pytest report the test function will show twice if an unknown warning is captured: one with the
    actual result of the test and one with an ERROR. The warning will be shown as part of the stderr stream.

    Parameters
    ----------
    request standard injected service for pytest fixtures
    recwarn injected pytest service for accessing recorded warnings

    """
    def fin():
        warnings = [w for w in recwarn.list
                    if type(w.message) not in known_warnings
                    or str(w.message) not in known_warnings[type(w.message)]]

        assert 0 == len(warnings)

    request.addfinalizer(fin)
