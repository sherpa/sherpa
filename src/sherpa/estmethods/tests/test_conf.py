#
#  Copyright (C) 2019, 2020, 2021  Smithsonian Astrophysical Observatory
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

import numpy as np

import pytest

from sherpa.fit import Fit
from sherpa.data import Data1D
from sherpa.models.basic import Polynom1D
from sherpa.estmethods import Confidence
from sherpa import ui


@pytest.fixture
def setUp(hide_logging):

    x = [-13, -5, -3, 2, 7, 12]
    y = np.asarray([102.3, 16.7, -0.6, -6.7, -9.9, 33.2])
    err = np.ones(6) * 5
    data = Data1D('tst', x, y, err)
    mdl = Polynom1D('mdl')
    return data, mdl


def cmp_results(result):

    parnames = ('mdl.c0', 'mdl.c1', 'mdl.c2')
    parvals = np.array([-9.384889507344322, -2.4169154937347925,
                        0.47827334261018023])
    parmins = np.array([-2.917507940156449, -0.250889317129555,
                        -0.03126766429871808])
    parmaxes = np.array([2.91750794015645, 0.250889317129555,
                         0.03126766429871808])

    for got, expected in zip(result.parnames, parnames):
        assert got == expected

    assert len(result.parnames) == len(parnames)

    assert result.parvals == pytest.approx(parvals)
    assert result.parmins == pytest.approx(parmins)
    assert result.parmaxes == pytest.approx(parmaxes)


@pytest.mark.parametrize('thaw_c1', [True, False])
def test_low_level(thaw_c1, setUp):
    data, mdl = setUp
    if thaw_c1:
        mdl.c1.thaw()

    mdl.c2.thaw()
    f = Fit(data, mdl, estmethod=Confidence())
    mdl.c2 = 1
    f.fit()

    if not thaw_c1:
        mdl.c1.thaw()
        f.fit()

    result = f.est_errors()
    cmp_results(result)


@pytest.mark.parametrize('thaw_c1', [True, False])
def tst_ui(thaw_c1, setUp, clean_ui):
    data, mdl = setUp

    ui.load_arrays(1, data.x, data.y, data.staterror)
    ui.set_source(1, ui.polynom1d.mdl)
    if thaw_c1:
        ui.thaw(mdl.c1)

    ui.thaw(mdl.c2)
    mdl.c2 = 1
    ui.fit()

    if not thaw_c1:
        ui.thaw(mdl.c1)
        ui.fit()

    ui.conf()
    result = ui.get_conf_results()
    cmp_results(result)
