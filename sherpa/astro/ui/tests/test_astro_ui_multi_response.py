#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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
from sherpa.astro.data import DataRMF


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


def make_rmf(elo, ehi, offset):
    """offset is the number of channels anove/below"""

    assert elo.size == ehi.size
    nchans = elo.size

    if offset > 0:
        nokay = nchans - offset

        n_grp = np.ones(nchans, dtype=np.int16)
        n_chan = np.ones(nchans, dtype=np.int16)
        n_chan[offset] = 2
        n_chan[offset + 1:] = 3

        f_chan = np.ones(nchans, dtype=np.int16)
        f_chan[offset] = 1
        f_chan[offset + 1:] = np.arange(1, nokay, dtype=np.int16)

        matrix = np.asarray([1] * offset + [0.85, 0.15] + [0.3, 0.6, 0.1] * (nokay - 1))

    elif offset < 0:
        offset = -offset
        nokay = nchans - offset

        n_grp = np.ones(nchans, dtype=np.int16)
        n_chan = np.ones(nchans, dtype=np.int16)
        n_chan[:nokay - 1] = 3
        n_chan[nokay - 1] = 2

        f_chan = nchans * np.ones(nchans, dtype=np.int16)
        f_chan[:nokay - 1] = np.arange(offset, nchans - 1, dtype=np.int16)
        f_chan[nokay - 1] = nchans - 1

        matrix = np.asarray([0.3, 0.6, 0.1] * (nokay - 1) +
                            [0.35, 0.65] +
                            [1] * offset)

    else:
        n_grp = np.ones(nchans, dtype=np.int16)
        n_chan = 3 * np.ones(nchans, dtype=np.int16)
        n_chan[0] = 2
        n_chan[-1] = 2
        f_chan = np.arange(0, nchans, dtype=np.int16)
        f_chan[0] = 1
        matrix = np.asarray([0.85, 0.15] * 1 +
                            [0.3, 0.6, 0.1] * (nchans - 2) +
                            [0.35, 0.65] * 1)

    return DataRMF("dummy", detchans=nchans,
                   energ_lo=elo, energ_hi=ehi,
                   n_grp=n_grp, n_chan=n_chan,
                   f_chan=f_chan, matrix=matrix,
                   e_min=elo, e_max=ehi)


# The RMF has two components that are offset from each
# other.
#
def make_rmf2():
    """Create two RMFs that overlap"""

    elo, ehi, dset = get_egrid()

    rmf1 = make_rmf(elo, ehi, 20)
    rmf2 = make_rmf(elo, ehi, -10)
    return rmf1, rmf2, dset


@pytest.mark.parametrize("id", [None, 1, "three"])
def test_set_multi_arfs(id, clean_astro_ui):

    arf1, arf2, dset = make_arf2()

    ui.set_data(id, dset)
    d = ui.get_data(id=id)
    assert d.response_ids == []

    ui.set_arf(id=id, arf=arf1, resp_id=1)
    assert d.response_ids == [1]

    ui.set_arf(id=id, arf=arf2, resp_id=2)
    assert d.response_ids == [1, 2]


@pytest.mark.parametrize("id", [None, 1, "three"])
def test_set_multi_arfs_reorder(id, clean_astro_ui):

    arf1, arf2, dset = make_arf2()

    ui.set_data(id, dset)
    d = ui.get_data(id=id)
    assert d.response_ids == []

    ui.set_arf(id=id, arf=arf2, resp_id=2)
    assert d.response_ids == [2]

    ui.set_arf(id=id, arf=arf1, resp_id=1)
    assert d.response_ids == [2, 1]


@pytest.mark.parametrize("id", [None, 1, "three"])
def test_set_multi_rmfs(id, clean_astro_ui):

    rmf1, rmf2, dset = make_rmf2()

    ui.set_data(id, dset)
    d = ui.get_data(id=id)
    assert d.response_ids == []

    ui.set_rmf(id=id, rmf=rmf1, resp_id=1)
    assert d.response_ids == [1]

    ui.set_rmf(id=id, rmf=rmf2, resp_id=2)
    assert d.response_ids == [1, 2]


@pytest.mark.parametrize("id", [None, 1, "three"])
def test_set_multi_rmfs_reorder(id, clean_astro_ui):

    rmf1, rmf2, dset = make_rmf2()

    ui.set_data(id, dset)
    d = ui.get_data(id=id)
    assert d.response_ids == []

    ui.set_rmf(id=id, rmf=rmf2, resp_id=2)
    assert d.response_ids == [2]

    ui.set_rmf(id=id, rmf=rmf1, resp_id=1)
    assert d.response_ids == [2, 1]


def check_eval_multi_arf():
    """Test that the data is handled correctly

    For use by test_eval_multi_arf.
    """

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


def test_eval_multi_arf(clean_astro_ui):
    """See also test_eval_multi_arf_reorder"""

    arf1, arf2, dset = make_arf2()
    ui.set_data(1, dset)
    ui.set_arf(id=1, arf=arf1, resp_id=1)
    ui.set_arf(id=1, arf=arf2, resp_id=2)

    check_eval_multi_arf()


def test_eval_multi_arf_reorder(clean_astro_ui):
    """Change the order of setting the ARFs

    Should be the same as test_eval_multi_arf
    """

    arf1, arf2, dset = make_arf2()
    ui.set_data(1, dset)
    ui.set_arf(id=1, arf=arf2, resp_id=2)
    ui.set_arf(id=1, arf=arf1, resp_id=1)

    check_eval_multi_arf()


def check_eval_multi_rmf():
    """Test that the data is handled correctly

    For use by test_eval_multi_rmf.
    """

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

    # The model plot does care about the instrument. There should
    # be two equal responses, offset from each other. As this is
    # done for each energy they should "cancel out" (apart from
    # doubling the signal)  apart from the start/end bins
    #
    # I haven't been able to convince myself I understand the handling
    # at the start/end of the RMF, so I am using this just as a
    # regression test.
    #
    mplot = ui.get_model_plot()

    range_is_8 = np.arange(11, 968)
    assert mplot.y[range_is_8] == pytest.approx(8)

    # handle "overflow" bins
    assert mplot.y[0] == pytest.approx(84.6)
    assert mplot.y[-1] == pytest.approx(43.0)

    # handle bins below and above the offsets
    range_is_low = np.arange(1, 11)
    is_low = [4.2, 4., 4., 4., 4., 4., 4., 4., 5.2, 7.6]
    assert mplot.y[range_is_low] == pytest.approx(is_low)

    range_is_high = np.arange(968, 988)
    is_high = [6.8, 4.4, 4., 4., 4., 4., 4., 4., 4., 4.,
               4., 4., 4., 4., 4., 4., 4., 4., 4., 4.2]
    assert mplot.y[range_is_high] == pytest.approx(is_high)


def test_eval_multi_rmf(clean_astro_ui):
    """See also test_eval_multi_rmf_reorder"""

    rmf1, rmf2, dset = make_rmf2()
    ui.set_data(1, dset)
    ui.set_rmf(id=1, rmf=rmf1, resp_id=1)
    ui.set_rmf(id=1, rmf=rmf2, resp_id=2)

    check_eval_multi_rmf()


def test_eval_multi_rmf_reorder(clean_astro_ui):
    """Change the order of setting the RMFs

    Should be the same as test_eval_multi_rmf
    """

    rmf1, rmf2, dset = make_rmf2()
    ui.set_data(1, dset)
    ui.set_rmf(id=1, rmf=rmf2, resp_id=2)
    ui.set_rmf(id=1, rmf=rmf1, resp_id=1)

    check_eval_multi_rmf()
