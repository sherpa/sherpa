#
#  Copyright (C) 2016, 2018, 2021, 2022
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

# Basic tests of sherpa.astro.ui.utils routines.
#
# This supplements test_astro_ui_ui.py (in that tests are being moved
# out of a single file)
#

import logging
import os

import pytest

from sherpa.astro import ui
from sherpa.utils.err import ArgumentTypeErr, DataErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec


# Note that the logic in load_table_model is quite convoluted,
# since it has a lot of fall-through cases, so these checks are
# used in case of any refactoring. Since load_table_model and
# load_xstable_model are similar in these regards the tests are
# used for both.
#
# The tests are manually repeated, rather than sending in the
# load function as a parameter, as I was unable to get the latter
# to work well.
#

@requires_fits
def test_load_table_model_fails_with_dir(tmp_path):
    """Check that the function fails with invalid input: directory

    The temporary directory is used for this.
    """

    tmpdir = tmp_path / 'load_table_model'
    tmpdir.mkdir()

    ui.clean()
    assert ui.list_model_components() == []
    with pytest.raises(IOError):
        ui.load_table_model('tmpdir', str(tmpdir))

    assert ui.list_model_components() == []


@requires_fits
@requires_xspec
def test_load_xstable_model_fails_with_dir(tmp_path):
    """Check that the function fails with invalid input: directory

    The temporary directory is used for this.
    """

    tmpdir = tmp_path / 'load_xstable_model'
    tmpdir.mkdir()

    ui.clean()
    assert ui.list_model_components() == []
    with pytest.raises(IOError):
        ui.load_xstable_model('tmpdir', str(tmpdir))

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

    # The error depends on the load function
    with pytest.raises(ValueError):
        ui.load_table_model('devnull', '/dev/null')

    assert ui.list_model_components() == []


@requires_fits
@requires_xspec
@pytest.mark.skipif(not(os.path.exists('/dev/null')),
                    reason='/dev/null does not exist')
def test_load_xstable_model_fails_with_dev_null():
    """Check that load_table_model fails with invalid input: /dev/null

    This simulates an empty file (and relies on the system
    containing a /dev/null file that reads in 0 bytes).
    """

    ui.clean()
    assert ui.list_model_components() == []

    # The error depends on the load function
    with pytest.raises(IOError):
        ui.load_xstable_model('devnull', '/dev/null')

    assert ui.list_model_components() == []


@requires_fits
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

    # The error depends on the load function.
    with pytest.raises(Exception):
        ui.load_table_model('stringcol', infile)

    assert ui.list_model_components() == []


@requires_fits
@requires_xspec
@requires_data
def test_load_xstable_model_fails_with_text_column(make_data_path):
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

    # The error depends on the load function.
    with pytest.raises(Exception):
        ui.load_xstable_model('stringcol', infile)

    assert ui.list_model_components() == []


class NonIterableObject:
    """Something that tuple(..) of will error out on"""

    pass


@pytest.mark.parametrize("func",
                         [ui.notice2d_id, ui.ignore2d_id,
                          ui.notice2d_image, ui.ignore2d_image])
def test_spatial_filter_errors_out_invalid_id(func):
    """Just check we create the expected error message.

    Somewhat contrived.
    """

    ids = NonIterableObject()
    with pytest.raises(ArgumentTypeErr,
                       match="'ids' must be an identifier or list of identifiers"):
        func(ids)


@requires_data
@requires_fits
def test_1160_original(make_data_path, clean_astro_ui, caplog):
    """Close to the originally-reported bug"""

    with SherpaVerbosity("ERROR"):
        ui.load_pha(make_data_path("3c273.pi"))

    assert len(caplog.records) == 0
    ui.notice(0.5, 6)
    assert len(caplog.records) == 1
    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 0.00146:14.9504 -> 0.4672:6.57 Energy (keV)"

    # Use data.get_filter rather than ui.get_filter as can
    # control the accuracy.
    #
    data = ui.get_data()

    fstr = "0.4672:6.5700"
    assert data.get_filter(format="%.4f") == fstr

    ui.set_grouping([1] * 1024)
    assert len(caplog.records) == 2

    assert data.get_filter(format="%.4f") == fstr

    r = caplog.record_tuples[1]
    assert r[0] == "sherpa.ui.utils"
    assert r[1] == logging.INFO
    assert r[2] == "dataset 1: 0.4672:6.57 Energy (keV) (unchanged)"


def test_get_dep_1160(clean_astro_ui):
    """Bizarre behavior was observed trying to fix #1160 so test it.

    Hold on: this is because the ui get_dep normalizes by the
    bin width (and it's likely that we have other tests of
    this behavior), but I'm going to leave this here.
    """

    grping = [1, -1, 1, -1, 1]
    vals = [1, 2, 3, 4, 5]
    gvals = [3, 7, 5]
    dvals = [1.5, 3.5, 5]  # gvals divided by bin width

    ui.load_arrays(1, [1, 2, 3, 4, 5], vals, ui.DataPHA)
    ui.set_grouping(grping)
    ui.group()

    # These really should be a given, but test just in case
    pha = ui.get_data()
    assert pha.grouped
    assert pha.grouping == pytest.approx(grping)

    assert pha.get_dep(filter=False) == pytest.approx(vals)
    assert pha.get_dep(filter=True) == pytest.approx(gvals)

    # Now the actual checks. For PHA data ui.get_dep actually
    # calls pha.get_y with the pha._rate field set to False,
    # and not pha.get_dep.
    #
    assert ui.get_dep(filter=False) == pytest.approx(dvals)
    assert ui.get_dep(filter=True) == pytest.approx(dvals)
