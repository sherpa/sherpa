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

        
@requires_fits
@requires_data
def test_scaling_staterr(make_data_path):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files wehre made with
    dmtcalc sherpa-test-data/sherpatest/source.pi sherpa-test-data/sherpatest/source_edit.pi expression="stat_err=stat_err/exposure"
    dmcopy "sherpa-test-data/sherpatest/source_edit.pi[cols channel,pi,RATE=count_rate,stat_err]" sherpa-test-data/sherpatest/source_rate.pi clobber=yes

    '''
    ui.load_pha("phacounts", make_data_path("source.pi"), use_errors=True)
    ui.load_pha("pharate", make_data_path("source_rate.pi"), use_errors=True)
    assert np.allclose(ui.get_data("phacounts").counts,
                        ui.get_data("pharate").counts)
    assert np.allclose(ui.get_data("phacounts").staterror,
                        ui.get_data("pharate").staterror)
