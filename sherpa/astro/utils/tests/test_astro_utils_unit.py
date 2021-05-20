#
#  Copyright (C) 2017, 2018, 2020, 2021
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

from sherpa.astro import ui
from sherpa.astro.utils import do_group, filter_resp, range_overlap_1dint
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
    """Check we error out if lengths do not match."""

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
    term, addnig 1 to each element.
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
