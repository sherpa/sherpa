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

# Basic tests of sherpa.astro.ui.utils routines.
#
# This supplements test_astro_ui_ui.py (in that tests are being moved
# out of a single file) but uses pytest instead of unittest.
#

import os
import tempfile

import pytest

# from numpy.testing import assert_almost_equal

from sherpa.astro import ui
from sherpa.utils import requires_data, requires_fits

tmpdir = tempfile.gettempdir()


# Note that the logic in load_table_model is quite convoluted,
# since it has a lot of fall-through cases, so these checks are
# used in case of any refactoring.
#

@requires_fits
@pytest.mark.skipif(not(os.path.isdir(tmpdir)),
                    reason='temp directory does not exist')
def test_load_table_model_fails_with_dir():
    """Check that load_table_model fails with invalid input: directory

    The temporary directory is used for this (the test is skipped if
    it does not exist).
    """

    ui.clean()
    assert ui.list_model_components() == []
    with pytest.raises(IOError):
        ui.load_table_model('tmpdir', tmpdir)

    assert ui.list_model_components() == []


@requires_fits
@pytest.mark.skipif(not(os.path.exists('/dev/null')),
                    reason='/dev/null does not exist')
def test_load_table_model_fails_with_dev_null():
    """Check that load_table_model fails with invalid input: /dev/null

    This simulates an empty file (and relies on the system
    containing a /dev/null file that reads in 0 bytes).
    """

    ui.clean()
    assert ui.list_model_components() == []
    # There is currently no explicit test for no data, so an
    # internal error is raised.
    with pytest.raises(ValueError):
        ui.load_table_model('devnull', '/dev/null')

    assert ui.list_model_components() == []


@requires_data
def test_load_table_model_fails_with_text_column(make_data_path):
    """Check that load_table_model fails with invalid input: text column

    The first column is text (and an ASCII file) so it is
    expected to fail.
    """

    # Check that this file hasn't been changed (as I am re-using it for
    # this test)
    infile = make_data_path('table.txt')
    assert os.path.isfile(infile)

    ui.clean()
    assert ui.list_model_components() == []
    with pytest.raises(ValueError) as exc:
        ui.load_table_model('stringcol', infile)

    emsg = 'The file {} could not be loaded, probably '.format(infile) + \
           'because it contained spurious data and/or strings'
    assert str(exc.value) == emsg

    assert ui.list_model_components() == []
