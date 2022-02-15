#
#  Copyright (C) 2017, 2018, 2020, 2021, 2022
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

import numpy as np
import pytest

from sherpa.astro.data import DataPHA, DataIMG, DataIMGInt
from sherpa.astro import ui
from sherpa.astro import utils
from sherpa.astro.utils import do_group, filter_resp, range_overlap_1dint
from sherpa.data import Data1D, Data1DInt, Data2D, Data2DInt
from sherpa.utils.err import IOErr
from sherpa.utils.testing import requires_data, requires_fits


# See https://github.com/sherpa/sherpa/issues/405
def test_filter_resp_nochans():
    """What happens if no channels?

    It could be an error, or empty arrays could be returned.
    """

    # fake up an RMF
    ngrp = [1, 1]
    f_chan = [1, 2]
    n_chan = [1, 1]
    matrix = [1.0, 1.0]
    offset = 1

    with pytest.raises(ValueError) as excinfo:
        filter_resp([], ngrp, f_chan, n_chan, matrix, offset)

    emsg = "There are no noticed channels"
    assert str(excinfo.value) == emsg


# See https://github.com/sherpa/sherpa/issues/405
@requires_data
@requires_fits
def test_rmf_filter_no_chans(make_data_path):
    """This is not really a unit test but placed here as it
    is related to test_filter_resp_nochans.
    """

    rmffile = make_data_path('3c273.rmf')
    rmf = ui.unpack_rmf(rmffile)
    with pytest.raises(ValueError) as excinfo:
        rmf.notice([])

    emsg = "There are no noticed channels"
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("lo, hi, expected",
                         [(None, None, [1] * 5),
                          (0, 100, [1] * 5),
                          (20, 90, [1] * 5),
                          (1, 10, None),
                          (1, 20, None),
                          (90, 100, None),
                          (100, 110, None),
                          (10, 19, None),
                          (10, 10, None),
                          (91, 91, None),
                          # edge handling for densities is different at
                          # low and high ends
                          (20, 20, [1, 0, 0, 0, 0]),
                          (30, 30, [0, 1, 0, 0, 0]),
                          (44, 44, [0, 0, 1, 0, 0]),
                          (60, 60, [0, 0, 0, 1, 0]),
                          (80, 80, [0, 0, 0, 0, 1]),
                          (90, 0, None),
                          # limits fall within each bin
                          (21, 21, [1, 0, 0, 0, 0]),
                          (32, 32, [0, 1, 0, 0, 0]),
                          (54, 54, [0, 0, 1, 0, 0]),
                          (79, 79, [0, 0, 0, 1, 0]),
                          (85, 85, [0, 0, 0, 0, 1]),
                          # ranges - single bin (exact)
                          (20, 30, [1, 0, 0, 0, 0]),
                          (30, 44, [0, 1, 0, 0, 0]),
                          (44, 60, [0, 0, 1, 0, 0]),
                          (60, 80, [0, 0, 0, 1, 0]),
                          (80, 90, [0, 0, 0, 0, 1]),
                          # ranges - single bin (subset)
                          (21, 29, [0.8, 0, 0, 0, 0]),
                          (30, 37, [0, 0.5, 0, 0, 0]),
                          (51, 60, [0, 0, 0.5625, 0, 0]),
                          (68, 72, [0, 0, 0, 0.2, 0]),
                          (81, 82, [0, 0, 0, 0, 0.1]),
                          # ranges - partial overlap, all within grid
                          (27, 84, [0.3, 1, 1, 1, 0.4]),
                          (32, 78, [0, 12 / 14, 1, 0.9, 0]),
                          (20, 24, [0.4, 0, 0, 0, 0]),
                          (83, 90, [0, 0, 0, 0, 0.7]),
                          # partial overlap, but lo or hi past grid
                          (0, 26, [0.6, 0, 0, 0, 0]),
                          (0, 30, [1.0, 0, 0, 0, 0]),
                          (0, 37, [1.0, 0.5, 0, 0, 0]),
                          (0, 44, [1.0, 1.0, 0, 0, 0]),
                          (26, 100, [0.4, 1, 1, 1, 1]),
                          (30, 100, [0, 1, 1, 1, 1]),
                          (44, 100, [0, 0, 1, 1, 1]),
                          # partial overlap, starting/ending on grid
                          (20, 61, [1.0, 1.0, 1.0, 0.05, 0]),
                          (20, 80, [1.0, 1.0, 1.0, 1.0, 0]),
                          (20, 81, [1.0, 1.0, 1.0, 1.0, 0.1]),
                          (20, 50, [1, 1, 0.375, 0, 0]),
                          (75, 90, [0, 0, 0, 0.25, 1]),
                          (46, 90, [0, 0, 0.875, 1, 1])
                          ])
@pytest.mark.parametrize("reverse", [False, True])
def test_range_overlap_1dint_ascending(lo, hi, expected, reverse):

    grid = np.asarray([20, 30, 44, 60, 80, 90])

    # wavelength grids are in descending order, but first/second
    # elements of axes are still in low/high order.
    #
    if reverse:
        grid = grid[::-1]
        axes = (grid[1:], grid[:-1])
        if expected is not None:
            expected = expected[::-1]
    else:
        axes = (grid[:-1], grid[1:])

    got = range_overlap_1dint(axes, lo, hi)
    if expected is None:
        assert got is None
        return

    assert got == pytest.approx(expected)


def test_do_group_invalid_scheme():
    """Check we error out and not segfault if the name is unknown."""

    with pytest.raises(ValueError) as ve:
        do_group([1, 2, 3], [1, 1, 1], 'foo')

    assert str(ve.value) == 'unsupported group function: foo'


def test_do_group_check_lengths():
    """Check we error out if lengths do not match."""

    with pytest.raises(TypeError) as te:
        do_group([1, 2, 3], [1, 1], 'sum')

    assert str(te.value) == 'input array sizes do not match, data: 3 vs group: 2'


@pytest.mark.parametrize("func,expected",
                         [("sum", [3, 12, 6]),
                          ("_sum_sq", np.sqrt([5, 50, 36])),
                          ("_min", [1, 3, 6]),
                          ("_max", [2, 5, 6]),
                          ("_middle", [1.5, 4, 6]),
                          ("_make_groups", [1, 2, 3])
                         ])
def test_do_group_check_func(func, expected):
    """Test grouping functions on simple input."""

    ans = do_group([1, 2, 3, 4, 5, 6], [1, -1, 1, -1, -1, 1], func)
    assert ans == pytest.approx(expected)


# Note that the grouping code expects the data to be >= 0 so
# support for -1 in the following is basically "an extra" (i.e. we
# could decide to error out if values < 0 are included, which
# would mean the test would fail, but for now it is supported).
#
@pytest.mark.parametrize("func", ["sum", "_min", "_max", "_middle"])
def test_do_group_nothing_func(func):
    """Special case: no grouping"""

    ans = do_group([6, 4, 2, 1, -1, 3], [1, 1, 1, 1, 1, 1], func)
    assert ans == pytest.approx([6, 4, 2, 1, -1, 3])


def test_do_group_nothing_sumsq():
    """Special case: no grouping"""

    ans = do_group([6, 4, 2, 1, -1, 3], [1, 1, 1, 1, 1, 1], "_sum_sq")
    assert ans == pytest.approx([6, 4, 2, 1, 1, 3])


def test_do_group_nothing_group():
    """Special case: no grouping

    The _make_groups uses the first element in the input as the starting
    term, adding 1 to each element.
    """

    ans = do_group([6, 4, 2, 1, -1, 3], [1, 1, 1, 1, 1, 1], "_make_groups")
    assert ans == pytest.approx([6, 7, 8, 9, 10, 11])


@pytest.mark.parametrize("func,expected",
                         [('sum', 15),
                          ('_sum_sq', np.sqrt(67)),
                          ('_min', -1), ('_max', 6),
                          ('_middle', 2.5), ('_make_groups', 6)])
def test_do_group_single_group(func, expected):
    """There's only one group."""

    ans = do_group([6, 4, 2, 1, -1, 3], [1, -1, -1, -1, -1, -1], func)
    assert ans == pytest.approx([expected])


def make_data(data_class):
    """Create a test data object of the given class.

    Using a string means it is easier to support the various PHA
    "types" - eg basic, grouping, grouping+quality.

    """

    x0 = np.asarray([1, 3, 7, 12])
    y = np.asarray([2, 3, 4, 5])
    if data_class == "1d":
        return Data1D('x1', x0, y)

    if data_class == "1dint":
        return Data1DInt('xint1', x0, np.asarray([3, 5, 8, 15]), y)

    chans = np.arange(1, 5)
    if data_class == "pha":
        return DataPHA('pha', chans, y)

    # We want to provide PHA tests that check out the grouping and
    # quality handling (but it is not worth trying all different
    # variants), so we have "grp" for grouping and no quality [*], and
    # "qual" for grouping and quality.
    #
    # [*] by which I mean we have not called ignore_bad, not that
    # there is no quality array.
    #
    grp = np.asarray([1, -1, 1, 1])
    qual = np.asarray([0, 0, 2, 0])
    pha = DataPHA('pha', chans, y, grouping=grp, quality=qual)
    if data_class == "grp":
        return pha

    if data_class == "qual":
        pha.ignore_bad()
        return pha

    x0 = np.asarray([1, 2, 3] * 2)
    x1 = np.asarray([1, 1, 1, 2, 2, 2])
    y = np.asarray([2, 3, 4, 5, 6, 7])
    if data_class == "2d":
        return Data2D('x2', x0, x1, y, shape=(2, 3))

    if data_class ==  "2dint":
        return Data2DInt('xint2', x0, x1, x0 + 1, x1 + 1, y, shape=(2, 3))

    if data_class == "img":
        return DataIMG('img', x0, x1, y, shape=(2, 3))

    if data_class == "imgint":
        return DataIMGInt('imgi', x0, x1, x0 + 1, x1 + 1, y, shape=(2, 3))

    assert False


@pytest.mark.parametrize("data_class", ["1d", "1dint", "pha", "grp", "qual"])
def test_calc_data_sum_invalid_range(data_class):
    """lo > hi"""

    data = make_data(data_class)
    with pytest.raises(IOErr) as err:
        utils.calc_data_sum(data, 10, 2)

    # TODO: This error message is not great.
    #
    assert str(err.value) == "the energy range is not consistent, 10 !< 2"


@pytest.mark.parametrize("data_class", ["1d", "1dint", "pha", "grp", "qual"])
def test_calc_data_sum_no_range(data_class):
    """Call calc_data_sum(data)

    The following comment holds for the remaining calc_data_sum tests:
    we are primarily interested in using these as regression tests. In
    most places we can explain what is going on based on the behavior
    of the notice methods, but there may be places where it is less
    clear.

    One such place is the handling of channels that have been removed
    by a quality filter (i.e. a DataPHA object on which ignore_bad has
    been called that has quality != 0 channels). The code currently
    doesn't care about the quality filter, and this is what we test
    against, but we may decide to change this behavior at some point.

    As well as the test of calc_data_sum the tests often call a method
    on the data object before calling the method and then compare it
    to the same call after calling the calc routine. This is to check
    that we have restored the filter expression correctly.

    """

    data = make_data(data_class)
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data) == 14
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("name", ["lo", "hi"])
@pytest.mark.parametrize("limit,expected", [(0, 0),
                                            (1, 2),
                                            (2, 0),
                                            (3, 3),
                                            (4, 0),
                                            (6, 0),
                                            (7, 4),
                                            (11, 0),
                                            (12, 5),
                                            (13, 0),
                                            (15, 0)])
def test_calc_data_sum_only_one_limit_1d(name, limit, expected):
    """Call calc_data_sum(data, lo or hi): Data1D

    The code is written so that if either lo or hi is given and the
    other is None then we behave the same, hence the use of the name
    parameter for the test. The code is written as a regression test,
    rather than from first principles. The 1D/1DInt/PHA cases are
    subtly different, hence the different tests.

    """

    data = make_data("1d")
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data, **{name: limit}) == expected
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("name", ["lo", "hi"])
@pytest.mark.parametrize("limit,expected", [(0, 0),
                                            (1, 0),
                                            (2, 2),
                                            (3, 0),
                                            (4, 3),
                                            (6, 0),
                                            (7, 0),
                                            (11, 0),
                                            (12, 0),
                                            (13, 5),
                                            (15, 0)])
def test_calc_data_sum_only_one_limit_1dint(name, limit, expected):
    """Call calc_data_sum(data, lo or hi): Data1DInt"""

    data = make_data("1dint")
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data, **{name: limit}) == expected
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("name", ["lo", "hi"])
@pytest.mark.parametrize("limit,expected", [(0, 0),
                                            (1, 2),
                                            (2, 3),
                                            (3, 4),
                                            (4, 5),
                                            (5, 0),
                                            (6, 0),
                                            (15, 0)])
def test_calc_data_sum_only_one_limit_pha(name, limit, expected):
    """Call calc_data_sum(data, lo or hi): DataPHA"""

    data = make_data("pha")
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data, **{name: limit}) == expected
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("data_class", ["grp", "qual"])
@pytest.mark.parametrize("name", ["lo", "hi"])
@pytest.mark.parametrize("limit,expected", [(0, 0),
                                            (1, 5),
                                            (2, 5),
                                            (3, 4),
                                            (4, 5),
                                            (5, 0),
                                            (6, 0),
                                            (15, 0)])
def test_calc_data_sum_only_one_limit_pha_grouped(name, limit, expected, data_class):
    """Call calc_data_sum(data, lo or hi): DataPHA (group/quality)"""

    data = make_data(data_class)
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data, **{name: limit}) == expected
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("frange,expected", [((0, 20), 14),
                                             ((0, 10), 9),
                                             ((2, 3), 3),
                                             ((2, 10), 7),
                                             ((3, 7), 7),
                                             ((4, 6), 0),
                                             ((3, 12), 12),
                                             ((1, 13), 14),
                                             ((20, 22), 0)])
def test_calc_data_sum_filtered_1d(frange, expected):
    """Filter out everything and then check with both lo/hi set: Data1D

    For the 1D case the filter ranges should be inclusive. We do not
    test too many cases, as the assumption is that the existing
    notice/ignore tests cover these more thoroughly.

    """

    data = make_data("1d")
    data.ignore()
    assert not data.mask
    assert utils.calc_data_sum(data, *frange) == expected
    assert not data.mask


@pytest.mark.parametrize("frange,expected", [((0, 20), 14),
                                             ((0, 10), 9),
                                             ((2, 3), 2),
                                             ((2, 10), 9),
                                             ((3, 7), 3),
                                             ((4, 6), 3),
                                             ((3, 12), 7),
                                             ((1, 13), 14),
                                             ((20, 22), 0)])
def test_calc_data_sum_filtered_1dint(frange, expected):
    """Filter out everything and then check with both lo/hi set: Data1DInt

    For the 1DInt case the filter ranges should be inclusive for the
    lower limit and exclusive for the upper limit. We do not test too
    many cases, as the assumption is that the existing notice/ignore
    tests cover these more thoroughly.

    """

    data = make_data("1dint")
    data.ignore()
    assert not data.mask
    assert utils.calc_data_sum(data, *frange) == expected
    assert not data.mask


@pytest.mark.parametrize("frange,expected", [((0, 10), 14),
                                             ((0, 4), 14),
                                             ((2, 3), 7),
                                             ((2, 10), 12),
                                             ((1, 1), 2),
                                             ((2, 2), 3),
                                             ((3, 3), 4),
                                             ((4, 4), 5),
                                             ((3, 4), 9),
                                             ((3, 5), 9),
                                             ((3, 12), 9),
                                             ((1, 13), 14),
                                             ((20, 22), 0)])
def test_calc_data_sum_filtered_pha(frange, expected):
    """Filter out everything and then check with both lo/hi set: DataPHA

    We test a few more cases than the 1D/1DInt cases because we want
    to include tests of the bad-quality channel case (not relevant
    here with the "pha" test but it is relevant for the "grp" and
    "qual" data obects we test below in
    test_calc_data_sum_filtered_pha_grouped.

    """

    data = make_data("pha")
    data.ignore()
    assert not data.mask
    assert utils.calc_data_sum(data, *frange) == expected
    assert not data.mask


@pytest.mark.parametrize("data_class", ["grp", "qual"])
@pytest.mark.parametrize("frange,expected", [((0, 10), 14),
                                             ((0, 4), 14),
                                             ((2, 3), 9),
                                             ((2, 10), 14),
                                             ((1, 1), 5),
                                             ((2, 2), 5),
                                             ((3, 3), 4),  # this is technically filtered-out by bad quality
                                             ((4, 4), 5),
                                             ((3, 4), 9),
                                             ((3, 5), 9),
                                             ((3, 12), 9),
                                             ((1, 13), 14),
                                             ((20, 22), 0)])
def test_calc_data_sum_filtered_pha_grouped(frange, expected, data_class):
    """Filter out everything and then check with both lo/hi set: DataPHA (grouped/quality)"""

    data = make_data(data_class)
    data.ignore()
    assert not data.mask
    assert utils.calc_data_sum(data, *frange) == expected
    assert not data.mask


@pytest.mark.parametrize("data_class", ["2d", "img", "2dint", "imgint"])
def test_calc_data_sum_no_range_2d(data_class):
    """What happens when data is not 1D?

    It looks like we still sum up the data.
    """

    data = make_data(data_class)
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum(data) == 27
    assert data.get_dep(filter=True) == pytest.approx(orig)


@pytest.mark.parametrize("data_class", ["1d", "1dint",
                                        "pha", "grp", "qual",
                                        "2d", "2dint"])
def test_calc_data_sum2d_no_range_1d(data_class):
    """What happens when data is not an IMG class

    Note that this includes 2d/2dint classes as they all fail.

    """

    data = make_data(data_class)
    with pytest.raises(AttributeError) as err:
        utils.calc_data_sum2d(data)

    assert str(err.value).endswith(" object has no attribute 'notice2d'")


@pytest.mark.parametrize("data_class", ["img", "imgint"])
def test_calc_data_sum2d_no_range_2d(data_class):
    """Call calc_data_sum2d(data)

    The following comment holds for the remaining calc_data_sum2d tests:
    we are primarily interested in using these as regression tests. In
    most places we can explain what is going on based on the behavior
    of the notice2d methods, but there may be places where
    it is less clear.

    As well as the test of calc_data_sum2d the tests often call a method
    on the data object before calling the method and then compare it
    to the same call after calling the calc routine. This is to check
    that we have restored the filter expression correctly.

    """

    data = make_data(data_class)
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum2d(data) == 27
    assert data.get_dep(filter=True) == pytest.approx(orig)


# XFAIL: We can not use DataIMGInt here because of issue #1379
@pytest.mark.parametrize("data_class", ["img", pytest.param("imgnit", marks=pytest.mark.xfail)])
def test_calc_data_sum2d_filtered_2d(data_class):
    """Call calc_data_sum2d(data, region)"""

    data = make_data(data_class)
    data.notice2d("rect(0, 0, 2, 2)", ignore=True)
    assert np.iterable(data.mask)
    omask = data.mask.copy()
    orig = data.get_dep(filter=True).copy()
    assert utils.calc_data_sum2d(data, "rect(0, 0, 2, 3)") == 16
    assert data.mask == pytest.approx(omask)
    assert data.get_dep(filter=True) == pytest.approx(orig)
