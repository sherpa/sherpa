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

from sherpa import smoke
import pytest
import mock


def test_success():
    try:
        smoke()
    except SystemExit:
        pytest.fail("smoke test should have passed")


def test_failure():
    with pytest.raises(SystemExit) as cm:
        smoke(require_failure=True)

    assert "Test failures were detected" == str(cm.value)


@mock.patch.dict('sys.modules', astropy=None)
def test_fits_failure():
    with pytest.raises(SystemExit) as cm:
        smoke(fits="astropy")

    assert "Requested fits as astropy but module not found" == str(cm.value)


@mock.patch.dict('sys.modules', values={"sherpa.astro.xspec": None})
def test_xspec_failure():
    with pytest.raises(SystemExit) as cm:
        smoke(xspec=True)

    assert "Requested xspec as xspec but module not found" == str(cm.value)