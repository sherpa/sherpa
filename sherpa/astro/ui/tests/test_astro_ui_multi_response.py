#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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

"""
Test the handling of "multiple" responses - that is, those handled by
load_multi_arfs and load_multi_rmfs in the UI layer.

Note that the load_multi_xxx commands are just simple wrappers
around load_arf/rmf, so this is aiming at checking we are passing
the right data down to the lower-level elements, and that they
are working correctly once done so.

"""

import numpy as np

import pytest

from sherpa.astro import ui
from sherpa.astro.instrument import create_arf


def get_egrid():
    """The standard energy grid / dataset for the tests."""

    egrid = np.arange(0.1, 10, 0.01)

    nchans = egrid.size - 1
    assert nchans == 989

    chans = np.arange(1, nchans + 1)
    y = np.zeros_like(chans)
    dset = ui.DataPHA('faked', chans, y)

    return egrid[:-1].copy(), egrid[1:].copy(), dset


def expected_arf2():
    """What is the composite ARF expected to be"""

    nchans = 989
    y1 = np.ones(nchans) * 1.1
    y2 = np.ones(nchans) * 0.7

    chans = np.arange(1, nchans + 1)

    idx, = np.where(chans < 490)
    y1[~idx] = 0.1
    y2[idx] = 0

    idx, = np.where((chans > 591) & (chans < 611))
    y1[idx] = 0

    return y1, y2


# The ARF is split into two components that overlap.
# The exposure times are assumed to be the same (not sure
# if we support different exposure times).
#
# ARF1 overlaps all of ARF2, except for a small
# energy range
#
def make_arf2():
    """Create two ARFs that overlap"""

    elo, ehi, dset = get_egrid()
    y1, y2 = expected_arf2()

    arf1 = create_arf(elo, ehi, specresp=y1)
    arf2 = create_arf(elo, ehi, specresp=y2)
    return arf1, arf2, dset


@pytest.mark.parametrize("id", [None, 1, "three"])
def test_set_multi_arfs(id):

    arf1, arf2, dset = make_arf2()

    ui.clean()

    ui.set_data(id, dset)
    d = ui.get_data(id=id)
    assert d.response_ids == []

    ui.set_arf(id=id, arf=arf1, resp_id=1)
    assert d.response_ids == [1]

    ui.set_arf(id=id, arf=arf2, resp_id=2)
    assert d.response_ids == [1, 2]


def test_eval_multi_arf():

    arf1, arf2, dset = make_arf2()

    ui.clean()

    ui.set_data(1, dset)
    ui.set_arf(id=1, arf=arf1, resp_id=1)
    ui.set_arf(id=1, arf=arf2, resp_id=2)

    mdl = ui.create_model_component('const1d', 'mdl')
    mdl.c0 = 4
    ui.set_source(mdl)

    # The analysis setting appears to depend on how the
    # data is set up. This is just to check it is energy.
    #
    d = ui.get_data()
    assert d.units == 'energy'

    # Easiest way to evaluate the model is to grab the
    # data from plot_source / plot_model
    #
    # The source doesn't care about how the instrument is set
    # up.
    #
    splot = ui.get_source_plot()
    assert (splot.y == 4).all()

    # The model plot does care about the instrument
    #
    yarf1, yarf2 = expected_arf2()
    expected = (yarf1 + yarf2) * 4

    mplot = ui.get_model_plot()
    assert mplot.y == pytest.approx(expected)

    # spot checks (just to validate the expected arf code)
    #
    # ymid represents the region where ARF1 is 0
    #
    yfirst = 4 * 1.1
    ylast = 4 * (0.1 + 0.7)
    ymid = 4 * 0.7

    assert mplot.y[0] == pytest.approx(yfirst)
    assert mplot.y[-1] == pytest.approx(ylast)
    assert mplot.y[600] == pytest.approx(ymid)
