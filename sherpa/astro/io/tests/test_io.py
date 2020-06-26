#
#  Copyright (C) 2016, 2017, 2018, 2020  Smithsonian Astrophysical Observatory
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

import numpy as np
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.astro import ui

from tempfile import NamedTemporaryFile

@pytest.fixture
def setup(make_data_path):
    ui.clean()


@requires_data
@requires_fits
@requires_xspec
def test_mod_fits(make_data_path):
    tablemodelfile = make_data_path("xspec-tablemodel-RCS.mod")
    ui.load_table_model("tmod", tablemodelfile)
    tmod = ui.get_model_component("tmod")
    assert tmod.name == "xstablemodel.tmod"

@requires_fits
def test_warnings_are_gone_arrays():
    ui.load_arrays(1, [1, 2, 3], [4, 5, 6])
    #  We now have logic in conftest.py to catch white-listed warnings and fail on unexpected ones.
    #  We just need to make any warnings bubble up, here and in the following test.
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=True, clobber=True)
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=False, clobber=True)

@requires_fits
@requires_data
def test_warnings_are_gone_pha(make_data_path):
    pha = make_data_path("3c273.pi")
    ui.load_pha(pha)
    with NamedTemporaryFile() as f:
        ui.save_data(1, f.name, ascii=False, clobber=True)

