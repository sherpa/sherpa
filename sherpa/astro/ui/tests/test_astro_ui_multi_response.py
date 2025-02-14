#
#  Copyright (C) 2019, 2020, 2022 - 2024
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

"""
Test the handling of "multiple" responses - that is, those handled by
load_multi_arfs and load_multi_rmfs in the UI layer. There is also
now a check of handling multi-matrix RMF files.

Note that the load_multi_xxx commands are just simple wrappers
around load_arf/rmf, so this is aiming at checking we are passing
the right data down to the lower-level elements, and that they
are working correctly once done so.

"""

from typing import Sequence

import numpy as np

import pytest

from sherpa.astro import io, ui
from sherpa.astro.data import DataRMF
from sherpa.astro.instrument import RMF1D, create_arf, create_delta_rmf
from sherpa.astro.io.types import HeaderItem, Header, BlockList
from sherpa.utils.err import PlotErr
from sherpa.utils.testing import requires_data, requires_fits


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

    assert mplot.y[11:968] == pytest.approx(8)

    # handle "overflow" bins
    assert mplot.y[0] == pytest.approx(84.6)
    assert mplot.y[-1] == pytest.approx(43.0)

    # handle bins below and above the offsets
    range_is_low = np.arange(1, 11)
    is_low = [4.2, 4., 4., 4., 4., 4., 4., 4., 5.2, 7.6]
    assert mplot.y[range_is_low] == pytest.approx(is_low)

    is_high = [6.8, 4.4, 4., 4., 4., 4., 4., 4., 4., 4.,
               4., 4., 4., 4., 4., 4., 4., 4., 4., 4.2]
    assert mplot.y[968:988] == pytest.approx(is_high)


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


def check_eval_multi_arfrmf():
    """Test that the data is handled correctly

    For use by test_eval_multi_arfrmf.
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

    # Comparison to the "truth" is harder than the previous checks
    # so just hard-code it.
    #
    y = ui.get_model_plot().y

    assert y[0] == pytest.approx(93.06)
    assert y[1] == pytest.approx(4.62)
    assert y[2:479] == pytest.approx(4.4)
    assert y[479] == pytest.approx(3.2)
    assert y[480] == pytest.approx(0.8)
    assert y[481:498] == pytest.approx(0.4)
    assert y[498] == pytest.approx(1.24)
    assert y[499] == pytest.approx(2.92)
    assert y[500:570] == pytest.approx(3.2)
    assert y[570] == pytest.approx(3.08)
    assert y[571] == pytest.approx(2.84)
    assert y[572:589] == pytest.approx(2.8)
    assert y[589] == pytest.approx(2.92)
    assert y[590] == pytest.approx(3.16)
    assert y[591:968] == pytest.approx(3.2)
    assert y[968] == pytest.approx(3.08)
    assert y[969] == pytest.approx(2.84)
    assert y[970:987] == pytest.approx(2.8)
    assert y[987] == pytest.approx(2.94)
    assert y[988] == pytest.approx(30.1)


def test_eval_multi_arfrmf(clean_astro_ui):
    """See also test_eval_multi_arfrmf_reorder"""

    arf1, arf2, dset = make_arf2()
    rmf1, rmf2, dset = make_rmf2()

    ui.set_data(1, dset)
    ui.set_arf(id=1, arf=arf1, resp_id=1)
    ui.set_arf(id=1, arf=arf2, resp_id=2)

    ui.set_rmf(id=1, rmf=rmf1, resp_id=1)
    ui.set_rmf(id=1, rmf=rmf2, resp_id=2)

    check_eval_multi_arfrmf()


def test_eval_multi_arfrmf_reorder(clean_astro_ui):
    """See also test_eval_multi_arfrmf"""

    arf1, arf2, dset = make_arf2()
    rmf1, rmf2, dset = make_rmf2()

    ui.set_data(1, dset)

    ui.set_rmf(id=1, rmf=rmf2, resp_id=2)
    ui.set_rmf(id=1, rmf=rmf1, resp_id=1)

    ui.set_arf(id=1, arf=arf2, resp_id=2)
    ui.set_arf(id=1, arf=arf1, resp_id=1)

    check_eval_multi_arfrmf()


def setup_order_plot(make_data_path):
    """Set up a faked dataset with multiple orders."""

    pha = make_data_path('3c273.pi')
    ui.load_pha(pha)

    # It has already loaded in one response
    arf = ui.get_arf(resp_id=1)
    arf.specresp *= 0.5

    for order, scale in enumerate([0.4, 0.25], 2):
        ui.load_arf(make_data_path('3c273.arf'), resp_id=order)
        ui.load_rmf(make_data_path('3c273.rmf'), resp_id=order)

        arf = ui.get_arf(resp_id=order)
        arf.specresp *= scale

    ui.set_source(ui.powlaw1d.pl)

    ui.notice(0.5, 7)
    ui.ignore(3, 4)


@requires_data
@requires_fits
def test_get_order_plot_multi(make_data_path, clean_astro_ui):
    """Rather than fake data, use a known dataset.

    Here we pretend we have three orders but with the same
    response (except that the ARF is 0.5, 0.4, 0.25 of the
    normal ARF).
    """

    setup_order_plot(make_data_path)

    fplot = ui.get_fit_plot()
    oplot = ui.get_order_plot()

    # The idea is to compare the X range of plot_fit to plot_order
    # (but accessed just using the plot objects rather than creating
    # an actual plot).
    #
    # First some safety checks
    assert fplot.dataplot.xlo == pytest.approx(fplot.modelplot.xlo)
    assert fplot.dataplot.xhi == pytest.approx(fplot.modelplot.xhi)

    assert len(oplot.xlo) == 3
    assert len(oplot.xhi) == 3
    assert len(oplot.y) == 3

    assert oplot.xlo[1] == pytest.approx(oplot.xlo[0])
    assert oplot.xlo[2] == pytest.approx(oplot.xlo[0])

    assert oplot.xhi[1] == pytest.approx(oplot.xhi[0])
    assert oplot.xhi[2] == pytest.approx(oplot.xhi[0])

    # We know the y values are 0.5, 0.4, 0.25 times the original arf
    # so we can compare them.
    #
    assert oplot.y[1] == pytest.approx(oplot.y[0] * 0.4 / 0.5)
    assert oplot.y[2] == pytest.approx(oplot.y[0] * 0.25 / 0.5)

    xlo = oplot.xlo[0]
    xhi = oplot.xhi[0]
    assert len(xlo) == 564
    assert xlo[0] == pytest.approx(0.46720001101493835)
    assert xhi[-1] == pytest.approx(9.869600296020508)

    # The model plot is technically drawn the same way as the order plot
    # (ungrouped) but it uses different code (sherpa.astro.plot.ModelHistogram)
    # so let's compare.
    #
    mplot = ui.get_model_plot()
    assert mplot.xlo[0] == pytest.approx(0.46720001101493835)
    assert mplot.xhi[-1] == pytest.approx(9.869600296020508)

    # Also compare to the fit plot (which is grouped)
    #
    assert fplot.modelplot.xlo[0] == pytest.approx(0.46720001101493835)
    assert fplot.modelplot.xhi[-1] == pytest.approx(9.869600296020508)

    # How does the overall model plot y values compare to the three
    # orders?
    #
    # Unfortunately the selected groups are different so can not
    # directly compare.
    #
    # y = oplot.y[0] + oplot.y[1] + oplot.y[2]
    # assert y == pytest.approx(mplot.y)


@requires_data
@requires_fits
def test_plot_order_multi(make_data_path, clean_astro_ui):
    """Rather than fake data, use a known dataset.

    Most of the testing is done by test_get_order_plot_multi;
    this just checks we can make an actual plot.
    """

    setup_order_plot(make_data_path)

    # All we are doing is checking we can call plot_order.
    #
    ui.plot_order()


@requires_data
@requires_fits
@pytest.mark.parametrize("orders", [9, [9], [1, 2, 9]])
def test_get_plot_order_knows_bad_response_multiple(orders, make_data_path, clean_astro_ui):
    """Check we error out when we have multiple orders"""

    setup_order_plot(make_data_path)
    with pytest.raises(PlotErr, match="'9' is not a valid order"):
        ui.get_order_plot(orders=orders)


@requires_data
@requires_fits
@pytest.mark.parametrize("orders", [2, [2], [1, 2]])
def test_get_plot_order_knows_bad_response_single(orders, make_data_path, clean_astro_ui):
    """Historically the order-plot code would not error out when there was
    only one response. So check on this behavior as a regression test.

    There is an argument to say that maybe we should just ignore the
    order argument when we know there's only one response, but this
    does not seem particularly sensible.

    """

    pha = make_data_path('3c273.pi')
    ui.load_pha(pha)
    ui.set_source(ui.powlaw1d.pl)

    # This currently works, no matter the order value.
    #
    ui.get_order_plot(orders=orders)


def create_multi_matrix_rmf(outfile: str
                            ) -> tuple[np.ndarray, np.ndarray,
                                       np.ndarray, np.ndarray,
                                       Sequence[float]]:
    """Use set_hdus to write out a multi-matrix RMF."""

    # First block is "perfect" and the second adds a "blob"
    # component.
    #
    egrid = np.linspace(0.1, 2.1, 21)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf1 = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi,
                            name="perfect")
    rmf1.header["TESTBLCK"] = 1

    # reduce the "energy" resolution of the second block.
    #
    e4grid = np.linspace(0.1, 2.1, 6)
    e4lo = e4grid[:-1]
    e4hi = e4grid[1:]
    matrix = np.zeros((5, 4), dtype=np.float32)
    blur = [0.1, 0.3, 0.4, 0.2]
    for i in range(5):
        matrix[i] = blur

    rmf2 = DataRMF("blob", detchans=rmf1.detchans,
                   energ_lo=e4lo, energ_hi=e4hi,
                   n_grp=np.ones(5, dtype=np.int32),
                   f_chan=np.asarray([2, 4, 6, 8, 10], dtype=np.int32),
                   n_chan=4 * np.ones(5, dtype=np.int32),
                   matrix=matrix.flatten(),
                   e_min=elo, e_max=ehi,
                   header={"TESTBLCK": 2})

    # Extract the data and then reconstruct a "multi" RMF.
    #
    blocks1 = io._pack_rmf(rmf1)
    blocks2 = io._pack_rmf(rmf2)

    header = Header([HeaderItem(name="TESTKEY", value=12)])
    blocks = [blocks1.blocks[0], blocks2.blocks[0], blocks1.blocks[1]]
    blist = BlockList(blocks=blocks, header=header)

    io.backend.set_hdus(outfile, blist)
    return (elo, ehi, e4lo, e4hi, blur)


@requires_fits
def test_load_multi_matrix_rmf(tmp_path, clean_astro_ui, caplog):
    """Can we load in a multi-matrix RMF?

    We have to create the file on disk since we do not have one to
    add to sherpa-test-data at this time. Based on
    sherpa/astri/io/tests/test_io_response.py::test_read_multi_matrix_rmf

    """

    outpath = tmp_path / "multi.rmf"
    outfile = str(outpath)
    elo, ehi, e4lo, e4hi, blur = create_multi_matrix_rmf(outfile)

    ui.dataspace1d(1, 20, 1, dstype=ui.DataPHA)

    # What happens if we try to read this in? We could use unpack_rmf
    # but use load_rmf to check how that works
    #
    assert len(caplog.record_tuples) == 0
    ui.load_rmf(outfile)
    assert len(caplog.record_tuples) == 0

    rmf = ui.get_rmf()

    # Note that the RMF is wrapped within the RMF1D class here,
    # unlike sherpa/astri/io/tests/test_io_response.py.
    #
    assert isinstance(rmf, RMF1D)

    # This appears like the "high-res" version. I do not want to
    # test access to the internals of the DataMultiRMF class here.
    #
    assert rmf.detchans == 20
    assert rmf.offset == 1

    assert rmf.energ_lo == pytest.approx(elo)
    assert rmf.energ_hi == pytest.approx(ehi)
    assert rmf.e_min == pytest.approx(elo)
    assert rmf.e_max == pytest.approx(ehi)
    assert rmf.n_grp == pytest.approx([1] * 20)

    # Set up the ARF.
    #
    yarf = np.full(20, 0.8)
    arf = create_arf(elo, ehi, specresp=yarf)
    ui.set_arf(arf)

    # When you apply it you get the combined result.
    #
    mdl = ui.create_model_component("gauss1d", "mdl")
    mdl.pos = 1.1
    mdl.fwhm = 0.8
    mdl.ampl = 1e4

    ui.set_source(mdl)

    # The "perfect" response.
    #
    expected_perfect = mdl(elo, ehi)

    # Create a 2D array to represent the blurry matrix
    #
    blurry_matrix = np.zeros((5, 20))
    for idx, fchan in enumerate([2, 4, 6, 8, 10]):
        blurry_matrix[idx, fchan - 1:fchan + 3] = blur

    expected_blurry = mdl(e4lo, e4hi) @ blurry_matrix

    # We include both matrices, and then include the "ARF".
    #
    expected = 0.8 * (expected_perfect + expected_blurry)

    # Correct for the bin widths.
    #
    yplot = expected / (ehi - elo)

    mplot = ui.get_model_plot()
    assert mplot.xlo == pytest.approx(elo)
    assert mplot.xhi == pytest.approx(ehi)
    assert mplot.y == pytest.approx(yplot)


@requires_fits
def test_save_multi_matrix_rmf(tmp_path, clean_astro_ui, caplog):
    """Can we write out a multi-matrix RMF?

    As we don't currently have one "in the wild" we manually create
    one. This replicates the logic from test_load_multi_matrix_rmf.

    This is a regression test since we currently only write out one of
    the matrices.

    """

    outfile1 = str(tmp_path / "multi1.rmf")
    elo, ehi, e4lo, e4hi, blur = create_multi_matrix_rmf(outfile1)

    yarf = np.full(20, 0.8)
    arf = create_arf(elo, ehi, specresp=yarf)

    mrmf = ui.unpack_rmf(outfile1)

    ui.dataspace1d(1, 20, 1, dstype=ui.DataPHA)
    ui.set_rmf(mrmf)
    ui.set_arf(arf)

    outfile2 = str(tmp_path / "multi2.rmf")
    assert len(caplog.record_tuples) == 0
    ui.save_rmf(outfile2)
    assert len(caplog.record_tuples) == 0

    # Try to validate the on-disk format using UI routines (so this is a
    # check of how much we can round-trip data).
    #
    ormf = ui.unpack_rmf(outfile2)
    assert isinstance(ormf, RMF1D)

    ui.dataspace1d(1, 20, 1, dstype=ui.DataPHA, id=2)
    ui.set_rmf(2, ormf)
    ui.set_arf(2, arf)

    mdl = ui.create_model_component("gauss1d", "g1")
    mdl.pos = 1.1
    mdl.fwhm = 0.8
    mdl.ampl = 1e4

    ui.set_source(1, mdl)
    ui.set_source(2, mdl)

    # The two plots should be different, as m1plot has the combination
    # response
    #
    # - high-res is "perfect"
    # - low-res is blurry
    #
    # and m2plot should just be the ideal plot. This is not needed here,
    # but once we can write out multi-matrix RMF files we can change to
    # testing equality of m1 and m2.
    #
    m1 = ui.get_model_plot(id=1, recalc=True).y.copy()
    m2 = ui.get_model_plot(id=2, recalc=True).y.copy()

    assert m2 != pytest.approx(m1)

    # Check that this is the "perfect" RMF being used. We need to account
    # for bin widths in the plot, and the ARF.
    #
    expected = 0.8 * mdl(elo, ehi) / (ehi - elo)
    assert m2 == pytest.approx(expected)
