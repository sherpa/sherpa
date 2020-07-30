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

"""Test handling of fitting a background to PHA data.

There are a few tests scattered around that do this, but this
attempts to provide a solid test base. There is likely to be
overlap in the WSTAT tests, for one.
"""

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.astro.data import DataPHA
from sherpa.utils.testing import requires_data, requires_fits, requires_group


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_filter(id, make_data_path, clean_astro_ui, hide_logging):
    """Can I change filtering of the background?"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    sdata = ui.get_data(id)
    bdata = ui.get_bkg(id)

    # Need an argument for notice/ignore_id
    iid = 1 if id is None else id

    # The idea here is to check that we can have
    # different filters for the source and background
    # regions.
    #
    assert sdata.get_dep(filter=True).size == 46
    assert bdata.get_dep(filter=True).size == 46

    ui.ignore_id(iid, None, 0.5)
    ui.ignore_id(iid, 7, None)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 40

    # need to reset the background filter
    ui.notice_id(iid, None, None, bkg_id=1)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 46

    ui.ignore_id(iid, None, 0.3, bkg_id=1)
    ui.ignore_id(iid, 8, None, bkg_id=1)

    assert sdata.get_dep(filter=True).size == 40
    assert bdata.get_dep(filter=True).size == 42


@requires_data
@requires_fits
@requires_group
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_group(id, make_data_path, clean_astro_ui, hide_logging):
    """Can I change grouping of the background?"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    sdata = ui.get_data(id)
    bdata = ui.get_bkg(id)

    dtype = np.dtype('>i2')  # why not np.dtype(np.int16)?
    assert sdata.grouping.dtype == dtype
    assert bdata.grouping.dtype == dtype

    assert sdata.quality.dtype == dtype
    assert bdata.quality.dtype == dtype

    assert sdata.grouped
    assert bdata.grouped

    # Check default grouping (taken from source)
    sgrouping = sdata.grouping.copy()
    assert sgrouping.sum() == -932
    assert (bdata.grouping == sgrouping).all()

    # Change the grouping of the background daga
    #
    ui.group_counts(id, 5, bkg_id=1)
    assert (sdata.grouping == sgrouping).all()
    assert bdata.grouping.sum() == -952

    # check we can ungroup just the background
    ui.ungroup(id, bkg_id=1)
    assert sdata.grouped
    assert not bdata.grouped

    assert bdata.grouping.sum() == -952


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_models(id, make_data_path, clean_astro_ui, hide_logging):
    """Is the model set correctly for the background"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    # Check the source models
    #
    ui.set_source(id, ui.powlaw1d.pl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl)

    smdl = ui.get_source(id)
    assert smdl.name == 'powlaw1d.pl'

    bmdl = ui.get_bkg_source(id)
    assert bmdl.name == 'powlaw1d.bpl'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((38564.608926889 * (powlaw1d.pl + 0.134921 * (powlaw1d.bpl)))))'

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf((38564.608926889 * powlaw1d.bpl)))'


@requires_data
@requires_fits
@pytest.mark.parametrize("id", [None, 1, "bgnd"])
def test_setup_pha1_file_models_two(id, make_data_path, clean_astro_ui, hide_logging):
    """Is the model set correctly for two background components"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(id, infile)

    # Tweak scaling values to make output easier to check.
    # The ARF is changed as well, to a different value, to
    # check what value is used in the model expression.
    #
    sdata = ui.get_data(id)
    sdata.exposure = 100
    sdata.get_arf().exposure = 150
    sdata.backscal = 0.1
    sdata.areascal = 0.8

    # scale factor to correct to source is
    #   = (100 / 1000) * (0.1 / 0.4) * (0.8 / 0.4)
    #     exp = 0.1   back = 0.25  area = 2
    #   = 0.025 if exclude the area scaling
    #   = 0.05  if include the area scaling
    #
    bdata1 = ui.get_bkg(id, bkg_id=1)
    bdata1.exposure = 1000
    # bdata1.get_arf().exposure = 1500  no ARF
    bdata1.backscal = 0.4
    bdata1.areascal = 0.4

    # add a second dataset (copy of the background)
    #
    # scale factor
    #   = (100 / 2000) * (0.1 / 0.8) * (0.8 / 0.5)
    #     exp = 0.05   back = 0.125  area = 1.6
    #   = 0.00625 if exclude the area scaling
    #   = 0.01    if include the area scaling
    #
    bdata2 = ui.unpack_pha(make_data_path('3c273_bg.pi'))
    bdata2.exposure = 2000
    # bdata2.get_arf().exposure = 2500  no ARF
    bdata2.backscal = 0.8
    bdata2.areascal = 0.5
    ui.set_bkg(id, bdata2, bkg_id=2)

    # Check the source models:
    #  - first with one
    #  - then with two
    #
    ui.set_source(id, ui.powlaw1d.pl)
    ui.set_bkg_source(id, ui.powlaw1d.bpl)

    smdl = ui.get_source(id)
    assert smdl.name == 'powlaw1d.pl'

    bmdl = ui.get_bkg_source(id)
    assert bmdl.name == 'powlaw1d.bpl'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((100 * (powlaw1d.pl + 0.03 * (powlaw1d.bpl)))))'

    bmdl = ui.get_bkg_model(id)
    assert bmdl.name == 'apply_rmf(apply_arf((1000 * powlaw1d.bpl)))'

    ui.set_bkg_source(id, ui.polynom1d.bpl2, bkg_id=2)

    bmdl = ui.get_bkg_model(id, bkg_id=1)
    assert bmdl.name == 'apply_rmf(apply_arf((1000 * powlaw1d.bpl)))'

    bmdl = ui.get_bkg_model(id, bkg_id=2)
    assert bmdl.name == 'apply_rmf(apply_arf((2000 * polynom1d.bpl2)))'

    smdl = ui.get_model(id)
    assert smdl.name == 'apply_rmf(apply_arf((100 * (powlaw1d.pl + 0.03 * (powlaw1d.bpl + polynom1d.bpl2)))))'
