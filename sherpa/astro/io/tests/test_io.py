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


def assert_staterr(use_errors):
    assert np.all(ui.get_data("phacounts").counts == 
                  pytest.approx(ui.get_data("pharate").counts))
    if use_errors is True: 
        assert np.all(ui.get_data("phacounts").staterror == 
                      pytest.approx(ui.get_data("pharate").staterror))
    else:
        assert ui.get_data("phacounts").staterror is None
        assert ui.get_data("pharate").staterror is None
    
    for n in ['phacounts', 'pharate']:
        ui.group_bins(n, 16)
    ui.set_analysis('energy')
    ui.ignore(None, 3.)
    if use_errors is True:
        assert np.all(ui.get_data("phacounts").get_error() == 
                  pytest.approx(ui.get_data("pharate").get_error()))
    else:
        assert ui.get_data("phacounts").get_error() is None
        assert ui.get_data("pharate").get_error() is None
    
        
@requires_fits
@requires_data
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc source.pi source_edit.pi expression="stat_err=stat_err/exposure"
    dmcopy "source_edit.pi[cols channel,pi,RATE=count_rate,stat_err]" source_rate.pi clobber=yes
    # Then fix units to counts/s manually in the header
    '''
    ui.load_pha("phacounts", make_data_path("source.pi"),
                use_errors=use_errors)
    ui.load_pha("pharate", make_data_path("source_rate.pi"),
                use_errors=use_errors)
    for n in ['phacounts', 'pharate']:
        ui.load_arf(n, make_data_path("source.arf"))
        ui.load_rmf(n, make_data_path("source.rmf"))

    assert_staterr(use_errors)


@requires_fits
@requires_data
@pytest.mark.parametrize("use_errors", [True, False])
def test_scaling_staterr_pha2(make_data_path, use_errors):
    '''Regression test for https://github.com/sherpa/sherpa/issues/800

    Notes
    -----
    Test files are made with (in sherpa-test-data/sherpatest)
    dmtcalc 3c120_pha2.gz 3c120_pha2_edit.gz expression="RATE=(float)counts/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit.gz 3c120_pha2_edit2.gz expression="stat_err=stat_err/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit2.gz 3c120_pha2_edit3.gz expression="backgroun_up=(float)background_up/exposure" clobber=yes
    dmtcalc 3c120_pha2_edit3.gz 3c120_pha2_edit4.gz expression="backgroun_down=(float)background_down/exposure" clobber=yes

    dmcopy "3c120_pha2_edit4.gz[cols -COUNTS,-BACKGROUND_UP,-BACKGROUND_DOWN][cols *,BACKGROUND_UP=backgroun_up,BACKGROUND_DOWN=backgroun_down][col -COUNTS]" 3c120_pha2_rate.gz clobber=yes
    Then fix units manually.
    '''
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "phacounts")
    ui.load_pha(make_data_path("3c120_pha2_rate.gz"),
                use_errors=use_errors)
    ui.copy_data(9, "pharate")
    for n in ['phacounts', 'pharate']:
            ui.load_arf(n, make_data_path("3c120_meg_-1.arf.gz"))
            ui.load_rmf(n, make_data_path("3c120_meg_-1.rmf.gz"))
    assert_staterr(use_errors)
