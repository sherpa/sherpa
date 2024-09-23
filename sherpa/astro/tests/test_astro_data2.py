#
#  Copyright (C) 2020 - 2024
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

"""Continued testing of sherpa.astro.data."""

import logging
import pickle
import re
import warnings

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataIMG, DataIMGInt, DataPHA, DataRMF
from sherpa.astro.instrument import create_delta_rmf, matrix_to_rmf
from sherpa.astro import io
from sherpa.astro.io.wcs import WCS
from sherpa.data import Data2D, Data2DInt
from sherpa.models import Const1D, Delta2D, Polynom1D, Polynom2D
from sherpa.stats._statfcts import calc_chi2datavar_errors
from sherpa.utils import dataspace2d
from sherpa.utils.err import DataErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_group, requires_region, requires_wcs


def test_can_not_group_ungrouped():
    """Does setting the grouping setting fail with no data?"""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert not pha.grouped
    with pytest.raises(DataErr,
                       match="data set 'name' does not specify grouping flags"):
        pha.grouped = True


def test_pha_get_indep_when_all_filtered():
    """Regression test."""

    pha = DataPHA("pha", [1, 2, 3], [9, 7, 8])
    pha.ignore()

    # This is with 'pha.mask = False'
    #
    indep = pha.get_indep()
    assert len(indep) == 1
    assert indep[0] == pytest.approx([])


def test_pha_get_indep_when_all_filtered():
    """Regression test."""

    pha = DataPHA("pha", [1, 2, 3], [9, 7, 8])
    pha.mask = [False] * 3

    indep = pha.get_indep()
    assert len(indep) == 1
    assert indep[0] == pytest.approx([])


def test_pha_get_indep_when_partial_filtered():
    """Regression test."""

    pha = DataPHA("pha", [1, 2, 3], [9, 7, 8])
    pha.mask = [True, False, True]

    indep = pha.get_indep()
    assert len(indep) == 1
    assert indep[0] == pytest.approx([1, 3])


def test_get_mask_is_none():
    """This is a regression test.

    The test name is no-longer valid since the return value is now, as
    of 4.17.0, [True,True,True] and not None.

    """
    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 3)


def test_get_mask_is_none_when_all_filtered():
    """This is a regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.ignore()
    assert pha.mask is False
    assert pha.get_mask() is None


def test_get_noticed_channels_no_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.get_filter() == "1:3"
    assert pha.get_noticed_channels() == pytest.approx([1, 2, 3])
    assert pha.mask is True


def test_get_noticed_channels_no_filter_manual():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.mask = [True, True, True]
    assert pha.mask == pytest.approx([1, 1, 1])
    assert pha.get_filter() == "1:3"
    assert pha.get_noticed_channels() == pytest.approx([1, 2, 3])


def test_get_noticed_channels_partial_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.ignore(2, 2)
    assert pha.mask == pytest.approx([1, 0, 1])
    assert pha.get_filter() == "1,3"
    assert pha.get_noticed_channels() == pytest.approx([1, 3])


def test_get_noticed_channels_removed_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.ignore()
    assert pha.mask is False
    assert pha.get_filter() == ""
    assert pha.get_noticed_channels() == pytest.approx([])


def test_get_noticed_channels_removed_filter2():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.mask = [False, False, False]
    assert pha.mask == pytest.approx([0, 0, 0])
    assert pha.get_filter() == ""
    assert pha.get_noticed_channels() == pytest.approx([])


def test_get_noticed_expr_no_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.get_noticed_expr() == "1-3"


def test_get_noticed_expr_no_filter_manual():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.mask = [True, True, True]
    assert pha.get_noticed_expr() == "1-3"


def test_get_noticed_expr_partial_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.ignore(2, 2)
    assert pha.get_noticed_expr() == "1,3"


def test_get_noticed_expr_removed_filter():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.ignore()
    assert pha.get_noticed_expr() == "No noticed channels"


def test_get_noticed_expr_removed_filter2():
    """Regression test."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha.mask = [False, False, False]
    assert pha.get_noticed_expr() == "No noticed channels"


# Historically we have needed to use ndarray and not lists, so check.
@pytest.mark.parametrize("chans",
                         [np.asarray([1, 2, 3]),
                          [1, 2, 3]
                          ])
def test_get_filter_expr_channel(chans):
    """Check get_filter_expr is called"""

    pha = DataPHA('name', chans, [1, 1, 1])
    assert pha.get_filter_expr() == '1-3 Channel'

    pha.ignore(None, 1)
    assert pha.get_filter_expr() == '2-3 Channel'


# Historically we have needed to use ndarray and not lists, so check.
@pytest.mark.parametrize("chans",
                         [np.asarray([1, 2, 3]),
                          [1, 2, 3]
                          ])
def test_get_filter_is_empty(chans):

    pha = DataPHA('name', chans, [1, 1, 1])
    assert pha.get_filter() == '1:3'
    assert pha.mask is True

    pha.ignore()
    assert pha.get_filter() == ''
    assert pha.mask is False

    pha.mask = [False, False, False]
    assert pha.get_filter() == ''

def test_pha_get_noticed_channels_when_empty():
    """A regression test."""

    empty = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="^The size of 'empty' has not been set$"):
        empty.get_noticed_channels()


def test_pha_filter_when_empty(caplog):
    """A regression test."""

    empty = DataPHA("empty", None, None)
    assert len(caplog.record_tuples) == 0

    with SherpaVerbosity("INFO"):
        empty.ignore(hi=3)

    assert len(caplog.records) == 1

    emsg = "Skipping dataset empty: The size of 'empty' has not been set"
    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.INFO
    assert r[2] == emsg


def test_pha_get_mask_when_empty(caplog):
    """A regression test."""

    empty = DataPHA("empty", None, None)
    with pytest.raises(DataErr,
                       match="^The size of 'empty' has not been set$"):
        empty.get_mask()


def test_pha_channel_empty_remains_empty():
    """A regression test."""

    empty = DataPHA("empty", None, None)
    assert empty.channel is None
    empty.channel = None
    assert empty.channel is None


@pytest.mark.parametrize("chans", [None, [1, 2, 3]])
def test_pha_group_when_empty(chans):
    """A regression test."""

    empty = DataPHA("empty", chans, None)
    with pytest.raises(DataErr,
                       match="^data set 'empty' can not be grouped as channel or counts is not set"):
        empty.group_counts(5)


@pytest.mark.parametrize("chtype,expected,args",
                         [("channel", '1:10', []),
                          ("channel", '', [(False, 1, 10)]),
                          ("channel", '2:9', [(True, 2, 9)]),
                          ("channel", '2:3,7:9', [(True, 2, 9), (False, 4, 6)]),
                          ("channel", '1:4,7:9', [(True, 2, 9), (False, 4, 6), (True, 0, 4)]),
                          ("channel", '2:3,5:10', [(True, 2, 9), (False, 4, 6), (True, 5, 13)]),
                          ("channel", '', [(True, 2, 9), (False, 4, 6), (True, 5, 13), (False, 0, 13)]),
                          ("channel", '1:10', [(True, 2, 9), (False, 4, 6), (True, 0, 13)]),
                          # None checks
                          ("channel", '1:3', [(True, None, 3)]),
                          ("channel", '4:10', [(False, None, 3)]),
                          ("channel", '5:10', [(True, 5, None)]),
                          ("channel", '1:4', [(False, 5, None)]),
                          ("channel", '1:3,5:10', [(True, 5, None), (True, None, 3)]),
                          ("channel", '4', [(False, 5, None), (False, None, 3)]),
                          # a few checks of non-integer channel limits (we don't explicitly
                          # say what this means so just check we know what it does)
                          # These are no-longer valid
                          # ("channel", '3:7', [(True, 2.8, 7.9)]),
                          # ("channel", '3:7', [(True, 2.1, 7.2)]),
                          # ("channel", '1:2,8:10', [(False, 2.8, 7.9)]),
                          # energy
                          ("energy", '0.2:2.2', []),
                          ("energy", '', [(False, 0.3, 2.1)]),
                          ("energy", '', [(False, 0, 3)]),
                          ("energy", '0.4:2.0', [(True, 0.51, 1.98)]),
                          ("energy", '0.4:1.2,1.6:2.0', [(True, 0.51, 1.98), (False, 1.24, 1.51)]),
                          ("energy", '0.2:1.4,1.6:2.0', [(True, 0.51, 1.98), (False, 1.24, 1.51), (True, 0.001, 1.32)]),
                          ("energy", '0.4:1.2,1.4:2.2', [(True, 0.51, 1.98), (False, 1.24, 1.51), (True, 1.46, 12.2)]),
                          ("energy", '', [(True, 0.51, 1.98), (False, 1.24, 1.51), (True, 1.46, 12.2), (False, 0.01, 13)]),
                          ("energy", '0.2:2.2', [(True, 0.51, 1.98), (False, 1.24, 1.51), (True, 0.01, 13)]),
                          # None checks
                          ("energy", '0.2:0.8', [(True, None, 0.65)]),
                          ("energy", '0.8:2.2', [(False, None, 0.65)]),
                          ("energy", '0.8:2.2', [(True, 0.95, None)]),
                          ("energy", '0.2:0.8', [(False, 0.95, None)]),
                          ("energy", '0.2:0.8,1.0:2.2', [(True, 1.05, None), (True, None, 0.65)]),
                          ("energy", '0.2:2.2', [(True, 0.95, None), (True, None, 0.65)]),
                          ("energy", '0.8:1.0', [(False, 1.05, None), (False, None, 0.65)]),
                          ("energy", '', [(False, 0.95, None), (False, None, 0.65)]),
                          # wavelength
                          ("wave", '5.6:62.0', []),
                          ("wave", '', [(False, 1, 70)]),
                          ("wave", '6.2:31.0', [(True, 6.5, 25)]),
                          ("wave", '6.2:8.9,12.4:31.0', [(True, 6.5, 25), (False, 9.1, 12)]),
                          ("wave", '5.6:10.3,12.4:31.0', [(True, 6.5, 25), (False, 9.1, 12), (True, 1, 10)]),
                          ("wave", '6.2:8.9,10.3:62.0', [(True, 6.5, 25), (False, 9.1, 12), (True, 12, 70)]),
                          ("wave", '5.6:62.0', [(True, 6.5, 25), (False, 9.1, 12), (True, 1, 70)]),
                          # None checks
                          ("wave", '5.6:10.3', [(True, None, 9.1)]),
                          ("wave", '10.3:62.0', [(False, None, 9.1)]),
                          ("wave", '10.3:62.0', [(True, 12.0, None)]),
                          ("wave", '5.6:10.3', [(False, 12.0, None)]),
                          ("wave", '5.6:10.3,12.4:62.0', [(True, 12.5, None), (True, None, 9.1)]),
                          ("wave", '5.6:62.0', [(True, 12.0, None), (True, None, 9.1)]),
                          ("wave", '10.3:12.4', [(False, 12.5, None), (False, None, 9.1)]),
                          ("wave", '', [(False, 12.0, None), (False, None, 9.1)]),
                         ])
def test_pha_get_filter_checks_ungrouped(chtype, expected, args):
    """Check we get the filter we expect

    chtype is channel, energy, or wavelength
    expected is the expected response
    args is a list of 3-tuples of (flag, loval, hival) where
    flag is True for notice and False for ignore; they define
    the filter to apply
    """

    chans = np.arange(1, 11, dtype=int)
    counts = np.ones(10, dtype=int)
    pha = DataPHA('data', chans, counts)

    # Use an ARF to create a channel to energy mapping
    # The 0.2-2.2 keV range maps to 5.636-61.992 Angstrom
    #
    egrid = 0.2 * np.arange(1, 12)
    arf = DataARF('arf', egrid[:-1], egrid[1:], np.ones(10))
    pha.set_arf(arf)

    pha.units = chtype
    for (flag, lo, hi) in args:
        if flag:
            pha.notice(lo, hi)
        else:
            pha.ignore(lo, hi)

    assert pha.get_filter(format='%.1f') == expected


def test_pha_get_xerr_all_bad_channel_no_group():
    """get_xerr handles all bad values [channel]

    It's not obvious what it is meant to be doing here.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  quality=[2, 2, 2])

    assert pha.get_xerr() == pytest.approx([0.5, 0.5, 0.5])

    pha.ignore_bad()
    assert pha.get_filter() == ''
    assert pha.get_xerr() == pytest.approx([0.5, 0.5, 0.5])


def test_pha_get_xerr_all_bad_channel_group():
    """get_xerr handles all bad values [channel]

    The behavior with grouping is different, presumably because
    we assume we have grouping when we have a quality array.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, 1, 1],
                  quality=[2, 2, 2])

    assert pha.get_xerr() == pytest.approx([0.5, 0.5, 0.5])

    assert pha.grouped
    pha.ignore_bad()
    assert pha.get_filter() == ''
    assert pha.get_xerr() == pytest.approx([])


def test_pha_get_xerr_all_bad_energy_no_group():
    """get_xerr handles all bad values [energy]

    It's not obvious what it is meant to be doing here.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  quality=[2, 2, 2])

    ebins = np.asarray([3.0, 5., 8.0, 12.0])
    rlo = ebins[:-1]
    rhi = ebins[1:]
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)
    pha.set_rmf(rmf)
    pha.units = 'energy'

    assert pha.get_xerr() == pytest.approx([1.0, 1.5, 2.0])

    pha.ignore_bad()
    assert pha.get_filter() == ''
    assert pha.get_xerr() == pytest.approx([1.0, 1.5, 2.0])


def test_pha_get_xerr_all_bad_energy_group():
    """get_xerr handles all bad values [energy]

    The behavior with grouping is different, presumably because
    we assume we have grouping when we have a quality array.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, 1, 1],
                  quality=[2, 2, 2])

    ebins = np.asarray([3.0, 5., 8.0, 12.0])
    rlo = ebins[:-1]
    rhi = ebins[1:]
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)
    pha.set_rmf(rmf)
    pha.units = 'energy'

    assert pha.get_xerr() == pytest.approx([1.0, 1.5, 2.0])

    assert pha.grouped
    pha.ignore_bad()

    # Should this error out or not?
    assert pha.get_filter() == ''
    # with pytest.raises(DataErr,
    #                    match="mask excludes all data"):
    #     pha.get_filter()

    assert pha.get_xerr() == pytest.approx([])


@pytest.mark.parametrize("ignore", [False, True])
@pytest.mark.parametrize("lbl,lo,hi", [('lo', 1.5, 2.5),
                                       ('lo', 1.5, 2),
                                       ('hi', 1, 2.5)])
def test_pha_channel_limits_are_integers(ignore, lbl, lo, hi):
    """Ensure channels are integers."""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, -1, 1])

    func = pha.ignore if ignore else pha.notice
    with pytest.raises(DataErr,
                       match=f"unknown {lbl} argument: 'must be an integer channel value'"):
        func(lo, hi)


def test_288_a():
    """The issue from #288 which was working"""

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    assert pha.mask is True
    pha.ignore(3, 4)

    # I use approx because it gives a nice answer, even though
    # I want eqiuality not approximation in this test. Fortunately
    # with bools the use of approx is okay (it can tell the
    # difference between 0 and 1, aka False and True).
    #
    assert pha.mask == pytest.approx([True, False, True])


def test_288_a_energy():
    """The issue from #288 which was working

    test_288_a but with a response so we test energy filters
    """

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    rlo = channels
    rhi = channels + 1
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    assert pha.mask is True
    pha.ignore(3, 4)

    # I use approx because it gives a nice answer, even though
    # I want equality not approximation in this test. Fortunately
    # with bools the use of approx is okay (it can tell the
    # difference between 0 and 1, aka False and True).
    #
    assert pha.mask == pytest.approx([True, False, True])


def test_288_b():
    """The issue from #288 which was failing

    We now error out with a non-integer channel
    """

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    assert pha.mask is True
    with pytest.raises(DataErr,
                       match="unknown lo argument: 'must be an integer channel value'"):
        pha.ignore(3.1, 4)


def test_288_b_energy():
    """The issue from #288 which was failing

    test_288_b but with a response so we test energy filters
    """

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    rlo = channels
    rhi = channels + 1
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    assert pha.mask is True
    pha.ignore(3.1, 4)

    assert pha.mask == pytest.approx([True, False, True])


@requires_group
def test_grouping_non_numpy():
    """Historically the group* calls would fail oddly if y is not numpy

    TypeError: grpNumCounts() Could not parse input arguments, please check input for correct type(s)

    This has now been addressed but the test has been left in.
    """

    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [0, 0, 0, 2, 1, 1, 0, 0, 0, 0]

    pha = DataPHA('416', x, y)
    pha.group_counts(3)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
    assert pha.quality == pytest.approx(quality)


@requires_group
def test_416_a():
    """The first test case from issue #416

    This used to use channels but it has been changed to add an RMF so
    we can filter in energy space, as it is not clear what non-integer
    channels should mean.

    """

    # if y is not a numpy array then group_counts errors out
    # with a strange error. Another reason why DataPHA needs
    # to validate input
    #
    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 0, 0])

    pha = DataPHA('416', x, y)

    rmf = create_delta_rmf(x, x + 1, e_min=x, e_max=x + 1,
                           name='416')
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    pha.notice(4.5, 6.5)

    # There are two ways to get the mask:
    #   - pha.mask returns the grouped mask (if the data is
    #     grouped)
    #   - pha.get_mask() always returns the ungrouped mask
    #
    mask_ungrouped = [False] * 3 + [True] * 3 + [False] * 4
    mask_grouped = [False] * 3 + [True] * 2 + [False] * 4

    assert pha.mask == pytest.approx(mask_ungrouped)
    assert pha.get_mask() == pytest.approx(mask_ungrouped)
    assert pha.grouping is None
    assert pha.quality is None

    # The grouping is done only for the noticed data range.
    pha.group_counts(3)

    assert pha.mask == pytest.approx(mask_grouped)
    assert pha.get_mask() == pytest.approx(mask_ungrouped)

    # Check we get the expected grouping: the first 3 and last 4
    # channels are excluded by the notice call above, so they are 0,
    # which means we only have 3 channels with data.
    #
    grouping = [0] * 3 + [1, -1, 1] + [0] * 4
    assert pha.grouping == pytest.approx(grouping)

    # As with grouoping, the first 3 and last 4 channels are not
    # changed, so have a quality of 0. For the remaining three
    # channels we know the first group is okay (this corresponds to
    # [1, -1] from the grouping array above, correspondinf to counts
    # [2, 1]) but the second group (the second [1] in grouping,
    # corresponding to counts of [1]) does not meet the grouping
    # criteria (it sums to 1!) and so has a quality of 2.
    #
    quality = [0] * 3 + [0, 0, 2] + [0] * 4
    assert pha.quality == pytest.approx(quality)

    dep = pha.get_dep(filter=True)
    assert dep == pytest.approx([3, 1])


@requires_group
def test_416_c():
    """The third test case from issue #416

    This used to use channels but it has been changed to add an RMF so
    we can filter in energy space, as it is not clear what non-integer
    channels should mean.

    """

    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 0, 0])

    pha = DataPHA('416', x, y)

    rmf = create_delta_rmf(x, x + 1, e_min=x, e_max=x + 1,
                           name='416')
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    # When using channels this used notice(3.5, 6.5)
    # but using energy space we need to use a different
    # range to match the ones the original channel filter
    # used.
    #
    pha.notice(4.5, 6.5)

    # this should be ~pha.mask
    tabstops = [True] * 3 + [False] * 3 + [True] * 4
    assert ~pha.mask == pytest.approx(tabstops)

    assert pha.grouping is None
    assert pha.quality is None
    pha.group_counts(3, tabStops=~pha.mask)

    grouping = [0] * 3 + [1, -1, 1] + [0] * 4
    assert pha.grouping == pytest.approx(grouping)

    # the second grouped bin has a quality of 2 as
    # it only contains 1 count
    quality = np.zeros(10, dtype=int)
    quality[5] = 2
    assert pha.quality == pytest.approx(quality)

    pha.ignore_bad()

    assert pha.grouping == pytest.approx(grouping)
    assert pha.quality == pytest.approx(quality)

    dep = pha.get_dep(filter=False)
    assert dep == pytest.approx(y)

    # It is not at all obvious why we get 8 bins returned
    # here. The ignore_bad has removed any existing
    # filters, but why do we get 8, not 10, values?
    # Well, one bin has been removed (quality=2)
    # and two bins have merged into 1. Hence the 8.
    #
    dep = pha.get_dep(filter=True)
    exp = np.zeros(8)
    exp[3] = 3
    assert dep == pytest.approx(exp)


@pytest.fixture
def make_test_image():
    """A simple image

    Note that normally you'd have logical axes of 1:31,
    1:21 here and then a WCS, but I've decided to number
    the axes differently (in physical units) as there is
    no requirement that the logical units are 1:nx/ny.
    """

    x1, x0 = np.mgrid[3830:3850, 4245:4275]

    # What is the ordering of shape? At the moment going for
    # NumPy style (ny, nx), but there is no credible documentation
    # (any documentation was added to describe the behavior we
    # have now).
    #
    shape = x0.shape

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.ones(x0.size)
    return DataIMG('d', x0, x1, y, shape=shape)


@pytest.fixture
def make_test_image_sky():
    """A simple image with just sky WCS

    """

    crval = [2000.5, -5000.5]
    cdelt = [2.0, 4.0]
    crpix = [-2.0, 3.0]
    sky = WCS("physical", "LINEAR", crval=crval, crpix=crpix, cdelt=cdelt)

    # logical: x=1, 2, 3
    #          y=1, 2
    #
    x1, x0 = np.mgrid[1:3, 1:4]
    shape = x0.shape

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.ones(x0.size)
    return DataIMG('sky-ey', x0, x1, y, shape=shape, sky=sky)


# This is a regression test - the values were calculated by
# Sherpa/WCS and not from first principles.
#
WORLD_X0 = np.asarray([30.10151131, 30., 29.89848869, 30.10154255, 30., 29.89845745])
WORLD_X1 = np.asarray([ 9.89998487, 9.9000001, 9.89998487, 9.99998461, 10., 9.99998461])


@pytest.fixture
def make_test_image_world():
    """A simple image with just world WCS

    """

    crval = [30, 10]
    cdelt = [-0.1, 0.1]
    crpix = [2.0, 2.0]
    eqpos = WCS("world", "WCS", crval=crval, crpix=crpix, cdelt=cdelt)

    # logical: x=1, 2, 3
    #          y=1, 2
    #
    x1, x0 = np.mgrid[1:3, 1:4]
    shape = x0.shape

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.ones(x0.size)
    return DataIMG('world-ey', x0, x1, y, shape=shape, eqpos=eqpos)


@pytest.fixture
def make_test_pha():
    """A simple PHA"""

    chans = np.asarray([1, 2, 3, 4], dtype=np.int16)
    counts = np.asarray([1, 2, 0, 3], dtype=np.int16)
    return DataPHA('p', chans, counts)


@pytest.fixture
def make_grouped_pha():
    """A simple PHA with grouping

    Note that this does have a quality=2 bin that is
    ignored.
    """

    chans = np.asarray([1, 2, 3, 4, 5], dtype=np.int16)
    counts = np.asarray([1, 2, 0, 3, 12], dtype=np.int16)
    grp = np.asarray([1, -1, -1, 1, 1], dtype=np.int16)
    qual = np.asarray([0, 0, 0, 0, 2], dtype=np.int16)
    pha = DataPHA('grp', chans, counts,
                  grouping=grp, quality=qual)
    pha.ignore_bad()
    return pha


@pytest.fixture
def make_quality_pha():
    """A simple PHA with grouping/quality data

    This is different to make_grouped_data as

    - it is more focussed on quality issues so has more
      channels, and more "bad" ones
    - does not call ignore_bad.

    """

    chans = np.arange(1, 10, dtype=np.int16)
    counts = np.asarray([1, 2, 0, 3, 12, 2, 9, 8, 7], dtype=np.int16)

    # The first group extends across the first set of "5" bad
    # channels, as does the last group (but that also includes a "5"
    # channel).
    grp = np.asarray([1, -1, -1, -1, 1, 1, 1, -1, -1], dtype=np.int16)
    qual = np.asarray([0, 5, 5, 0, 0, 0, 2, 2, 5], dtype=np.int16)
    return DataPHA('grp', chans, counts, grouping=grp, quality=qual)


def test_img_get_img(make_test_image):
    img = make_test_image
    ival = img.get_img()
    assert ival.shape == (20, 30)
    assert ival == pytest.approx(np.ones(20 * 30).reshape((20, 30)))


def test_img_get_img_filter_none1(make_test_image):
    """get_img when all the data has been filtered: mask is False

    It is not obvious what is meant to happen here (the docs suggest
    the filter is ignored but there is some handling of filters), so
    this should be treated as a regression test. See issue #1447

    """
    img = make_test_image
    img.notice2d(ignore=True)

    # safety check to ensure all the data has been ignored
    assert img.mask is False

    shape = (20, 30)
    expected = np.ones(shape)

    ival = img.get_img()
    assert ival.shape == shape
    assert ival == pytest.approx(expected)


@requires_region
def test_img_get_img_filter_none2(make_test_image):
    """get_img when all the data has been filtered: mask is array of False"""

    img = make_test_image
    img.notice2d("rect(0,0,10000,10000)", ignore=True)

    # safety check to ensure all the data has been ignored
    assert np.iterable(img.mask)
    assert not np.any(img.mask)

    shape = (20, 30)
    ival = img.get_img()
    assert ival.shape == shape
    assert not np.any(np.isfinite(ival))


@requires_region
def test_img_get_img_filter_some(make_test_image):
    """get_img when some of the data has been filtered.

    Unlike filtering out all the data, this does filter the response.

    """
    img = make_test_image
    # use a shape that's easy to filter
    img.notice2d("rect(4250, 3840,4256,3842)")

    # safety check to ensure that a subset of the data has been masked out
    assert np.iterable(img.mask)
    assert np.any(img.mask)
    assert not np.all(img.mask)

    # It looks like RECT is inclusive for low and high edges.
    shape = (20, 30)
    expected = np.zeros(20 * 30) * np.nan
    idx = np.hstack((np.arange(305, 312), np.arange(335, 342), np.arange(365, 372)))
    expected[idx] = 1
    expected.resize(shape)

    ival = img.get_img()
    assert ival.shape == shape

    # pytest.approx follows IEEE so nan != nan, hence we
    # have to filter out the values we expect.
    #
    good = np.isfinite(expected)
    assert np.isfinite(ival) == pytest.approx(good)
    assert ival[good] == pytest.approx(expected[good])


def image_callable(x0, x1):
    """Check that we call the routine correctly (DataIMG/get_img)"""

    assert len(x0) == 20 * 30
    assert len(x1) == 20 * 30
    assert x0[0] == pytest.approx(4245)
    assert x1[0] == pytest.approx(3830)
    assert x0[-1] == pytest.approx(4274)
    assert x1[-1] == pytest.approx(3849)
    return np.ones(x0.size) + 2


def image_callable_filtered(x0, x1):
    """Check that we call the routine correctly (DataIMG/get_img)"""

    assert len(x0) == 21
    assert len(x1) == 21
    assert x0[0] == pytest.approx(4250)
    assert x1[0] == pytest.approx(3840)
    assert x0[-1] == pytest.approx(4256)
    assert x1[-1] == pytest.approx(3842)
    return np.ones(x0.size) + 2


def image_callable_filtered2(x0, x1):
    """Check that we call the routine correctly (DataIMG/get_img)"""

    assert len(x0) == 11
    assert len(x1) == 11
    assert x0[0] == pytest.approx(4247)
    assert x1[0] == pytest.approx(3831)
    assert x0[-1] == pytest.approx(4248)
    assert x1[-1] == pytest.approx(3834)
    return np.ones(x0.size) + 2


def image_callable_none(x0, x1):
    """Check that we call the routine correctly (DataIMG/get_img)"""

    assert len(x0) == 0
    assert len(x1) == 0
    return np.asarray([])


def test_img_get_img_model(make_test_image):
    """What happens when we give a callable function to get_img?

    The idea is that it will be a model, but all we need is
    a callable.

    """
    img = make_test_image
    ival, mval = img.get_img(image_callable)

    shape = (20, 30)
    expected1 = np.ones(shape)
    expected2 = np.ones(shape) * 3

    # The data
    assert ival.shape == shape
    assert ival == pytest.approx(expected1)

    # The callable
    assert mval.shape == shape
    assert mval == pytest.approx(expected2)


def test_img_get_img_model_filter_none1(make_test_image):
    """See test_img_get_img_filter_none1. Issue #1447"""

    img = make_test_image
    img.notice2d(ignore=True)
    with pytest.raises(DataErr,
                       match="mask excludes all data"):
        img.get_img(image_callable)


@requires_region
def test_img_get_img_model_filter_none2(make_test_image):
    """See test_img_get_img_filter_none2. Issue #1447"""

    img = make_test_image
    img.notice2d("rect(2000,3000,7000,5000)", ignore=True)
    ival, mval = img.get_img(image_callable_none)

    shape = (20, 30)
    assert ival.shape == shape
    assert mval.shape == shape

    assert not np.any(np.isfinite(ival))
    assert not np.any(np.isfinite(mval))


@requires_region
def test_img_get_img_model_filter_some(make_test_image):
    """get_img with a callable and having a filter"""

    img = make_test_image
    # use a shape that's easy to filter
    img.notice2d("rect(4250, 3840,4256,3842)")

    ival, mval = img.get_img(image_callable_filtered)

    shape = (20, 30)
    idx = np.hstack((np.arange(305, 312), np.arange(335, 342), np.arange(365, 372)))

    expected1 = np.zeros(20 * 30) * np.nan
    expected1[idx] = 1
    expected1.resize(shape)

    expected2 = np.zeros(20 * 30) * np.nan
    expected2[idx] = 3
    expected2.resize(shape)

    assert ival.shape == shape
    assert mval.shape == shape

    # pytest.approx follows IEEE so nan != nan, hence we
    # have to filter out the values we expect.
    #
    good = np.isfinite(expected1)
    assert np.isfinite(ival) == pytest.approx(good)
    assert np.isfinite(mval) == pytest.approx(good)
    assert ival[good] == pytest.approx(expected1[good])
    assert mval[good] == pytest.approx(expected2[good])


@requires_region
def test_img_get_img_model_filter_some2(make_test_image):
    """test_img_get_img_model_filter_some but with a non-rectangular filter

    We have been using a filter that is rectangular, so matches the
    grid. Let's see what happens if the filter like a circle so that
    the bounding box does not match the filter.

    """

    img = make_test_image
    img.notice2d("circle(4247.8, 3832.1, 2)")

    # check
    assert img.mask.sum() == 11

    ival, mval = img.get_img(image_callable_filtered2)

    shape = (20, 30)
    idx = np.asarray([32, 33, 34, 61, 62, 63, 64, 92, 93, 94, 123])

    expected1 = np.zeros(20 * 30) * np.nan
    expected1[idx] = 1
    expected1.resize(shape)

    expected2 = np.zeros(20 * 30) * np.nan
    expected2[idx] = 3
    expected2.resize(shape)

    assert ival.shape == shape
    assert mval.shape == shape

    # pytest.approx follows IEEE so nan != nan, hence we
    # have to filter out the values we expect.
    #
    good = np.isfinite(expected1)
    assert np.isfinite(ival) == pytest.approx(good)
    assert np.isfinite(mval) == pytest.approx(good)
    assert ival[good] == pytest.approx(expected1[good])
    assert mval[good] == pytest.approx(expected2[good])


def test_img_can_not_set_coord(make_test_image):
    """The coord attribute is not writeable.

    It used to be, but now we require the user to change
    it with the set_coord method.
    """
    d = make_test_image

    # This dataset does not have a physical system, but we do not get
    # a DataErr but an AttributeError. The text message depends on
    # Python version (e.g. 3.9, 3.10, and 3.11 all have different
    # messages) so we do not check the message, just the class.
    #
    with pytest.raises(AttributeError):
        d.coord = "physical"


def test_img_set_coord_invalid(make_test_image):
    """An invalid coord setting"""
    d = make_test_image
    assert d.coord == 'logical'

    emsg = "unknown coordinates: 'bob'\n"
    emsg += "Valid options: logical, image, physical, world, wcs"
    with pytest.raises(DataErr,
                       match=emsg):
        d.set_coord('bob')

    assert d.coord == 'logical'


@pytest.mark.parametrize('coord,expected',
                         [('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_set_coord_notset(coord, expected, make_test_image):
    """A valid coord setting but we don't have the data"""

    d = make_test_image
    with pytest.raises(DataErr,
                       match=f"data set 'd' does not contain a {expected} coordinate system"):
        d.set_coord(coord)

    assert d.coord == 'logical'


@pytest.mark.parametrize('coord,expected',
                         [('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_get_coord_notset(coord, expected, make_test_image):
    """Check get_physical/world fail when there's no WCS"""
    d = make_test_image

    meth = getattr(d, f"get_{coord}")
    with pytest.raises(DataErr,
                       match=f"data set 'd' does not contain a {expected} coordinate system"):
        meth()


def test_img_set_coord_image(make_test_image):
    """Can set to image though"""
    d = make_test_image
    assert d.coord == 'logical'

    d.set_coord('image')
    assert d.coord == 'logical'


def test_img_get_coord_image(make_test_image):
    """Can call get_image though"""
    d = make_test_image

    cs = d.get_image()

    x1, x0 = np.mgrid[3830:3850, 4245:4275]
    x0 = x0.flatten()
    x1 = x1.flatten()

    assert cs[0] == pytest.approx(x0)
    assert cs[1] == pytest.approx(x1)
    assert len(cs) == 2


@pytest.fixture
def read_test_image(make_data_path):
    from sherpa.astro.io import read_image
    filename = 'acisf07999_000N001_r0035_regevt3_srcimg.fits'
    d = read_image(make_data_path(filename))
    d.name = 'test.img'
    return d


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord,expected',
                         [('logical', 'logical'),
                          ('image', 'logical'),
                          ('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_file_set_coord(coord, expected, read_test_image):
    """call set_coord with an image with WCS"""
    d = read_test_image
    assert d.coord == 'logical'
    d.set_coord(coord)
    assert d.coord == expected


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'image', 'physical', 'world', 'wcs'])
def test_img_file_get_logical(coord, read_test_image):
    """get_logical when coord is set"""
    d = read_test_image
    d.set_coord(coord)

    yexp, xexp = np.mgrid[1:378, 1:170]
    xexp = xexp.flatten()
    yexp = yexp.flatten()

    x, y = d.get_logical()
    assert x == pytest.approx(xexp)
    assert y == pytest.approx(yexp)


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'image', 'physical', 'world', 'wcs'])
def test_img_file_get_physical(coord, read_test_image):
    """get_physical when coord is set"""
    d = read_test_image
    d.set_coord(coord)

    yexp, xexp = np.mgrid[4333.1298828125:4710:1, 3062.3100585938:3231:1]
    xexp = xexp.flatten()
    yexp = yexp.flatten()

    x, y = d.get_physical()
    assert x == pytest.approx(xexp)
    assert y == pytest.approx(yexp)


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'image', 'physical', 'world', 'wcs'])
def test_img_file_get_world(coord, read_test_image):
    """get_world when coord is set"""
    d = read_test_image
    d.set_coord(coord)

    # Since the pixel size isn't guaranteed to be constant
    # just check the corners. Note that this is not a
    # check from first principles: it just checks that we
    # get the same answer as previous calls to this routine.
    #
    x, y = d.get_world()

    # BL
    assert x[0] == pytest.approx(150.02683651815326)
    assert y[0] == pytest.approx(2.6402818651328728)

    # TR
    assert x[-1] == pytest.approx(150.00385708212673)
    assert y[-1] == pytest.approx(2.6916707654223244)

    # BR
    assert x[168] == pytest.approx(150.00385224075313)
    assert y[168] == pytest.approx(2.640284264823834)

    # TL
    assert x[169 * 377 - 169] == pytest.approx(150.0268422985145)
    assert y[169 * 377 - 169] == pytest.approx(2.691668318963721)


def test_img_get_axes_logical(make_test_image):
    """Does get_axes work?"""
    d = make_test_image
    x, y = d.get_axes()

    assert x == pytest.approx(np.arange(1, 31, 1))
    assert y == pytest.approx(np.arange(1, 21, 1))


@requires_fits
@requires_data
@pytest.mark.parametrize('coord', ['logical', 'image'])
def test_img_file_get_axes_logical(coord, read_test_image):
    """get_axes when coord is set: logical"""
    d = read_test_image
    d.set_coord(coord)
    x, y = d.get_axes()

    assert x == pytest.approx(np.arange(1, 170, 1))
    assert y == pytest.approx(np.arange(1, 378, 1))


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['physical'])
def test_img_file_get_axes_physical(coord, read_test_image):
    """get_axes when coord is set: physical"""
    d = read_test_image
    d.set_coord(coord)
    x, y = d.get_axes()

    assert x == pytest.approx(np.arange(3062.3100585938, 3231, 1))
    assert y == pytest.approx(np.arange(4333.1298828125, 4710, 1))


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['world', 'wcs'])
def test_img_file_get_axes_world(coord, read_test_image):
    """get_axes when coord is set: world"""
    d = read_test_image
    d.set_coord(coord)
    x, y = d.get_axes()

    assert x.size == 169
    assert y.size == 377

    # This is an interesting combination of the corners from
    # test_img_file_get_world
    assert x[0] == pytest.approx(150.02683651815326)
    assert y[0] == pytest.approx(2.6402818651328728)
    assert x[-1] == pytest.approx(150.00385224075313)
    assert y[-1] == pytest.approx(2.691668318963721)


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord,expected',
                         [('logical', 'x0'),
                          ('image', 'x0'),
                          ('physical', 'x0 (pixels)'),
                          ('world', 'RA (deg)'),
                          ('wcs', 'RA (deg)')])
def test_img_file_get_x0label(coord, expected, read_test_image):
    """get_x0label"""
    d = read_test_image
    d.set_coord(coord)
    assert d.get_x0label() == expected


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'physical', 'world'])
@pytest.mark.parametrize('label', ['', 'not a label', 'x0'])
def test_img_file_set_x0label(coord, label, read_test_image):
    """set_x0label"""
    d = read_test_image
    d.set_coord(coord)
    d.set_x0label(label)
    assert d.get_x0label() == label


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord,expected',
                         [('logical', 'x1'),
                          ('image', 'x1'),
                          ('physical', 'x1 (pixels)'),
                          ('world', 'DEC (deg)'),
                          ('wcs', 'DEC (deg)')])
def test_img_file_get_x1label(coord, expected, read_test_image):
    """get_x1label"""
    d = read_test_image
    d.set_coord(coord)
    assert d.get_x1label() == expected


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'physical', 'world'])
@pytest.mark.parametrize('label', ['', 'not a label', 'x0'])
def test_img_file_set_x1label(coord, label, read_test_image):
    """set_x1label"""
    d = read_test_image
    d.set_coord(coord)
    d.set_x1label(label)
    assert d.get_x1label() == label


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'physical', 'world'])
def test_img_file_get_ylabel(coord, read_test_image):
    """get_ylabel"""
    d = read_test_image
    d.set_coord(coord)
    assert d.get_ylabel() == 'y'


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize('coord', ['logical', 'physical', 'world'])
@pytest.mark.parametrize('label', ['', 'not a label', 'y'])
def test_img_file_set_ylabel(coord, label, read_test_image):
    """set_ylabel"""
    d = read_test_image
    d.set_coord(coord)
    d.set_ylabel(label)
    assert d.get_ylabel() == label


def test_img_get_bounding_mask_nofilter(make_test_image):
    """get_bounding_mask with no filter"""
    d = make_test_image
    ans = d.get_bounding_mask()
    assert len(ans) == 2
    assert ans[0]
    assert ans[1] is None


def test_img_get_bounding_mask_nodata(make_test_image):
    """get_bounding_mask with all data masked"""
    d = make_test_image
    d.notice2d(ignore=True)
    ans = d.get_bounding_mask()
    assert len(ans) == 2
    assert not ans[0]
    assert ans[1] is None


@requires_region
def test_img_get_bounding_mask_filtered(make_test_image):
    """get_bounding_mask with data partially filtered"""
    d = make_test_image
    d.notice2d('ellipse(4260,3840,3,2,0)')
    ans = d.get_bounding_mask()

    mask = np.zeros(5 * 7, dtype=bool)
    for i in [3,  8,  9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24,
              25, 26, 31]:
        mask[i] = True

    assert len(ans) == 2
    assert ans[0] == pytest.approx(mask)
    assert ans[1] == (5, 7)


@requires_region
def test_img_get_filter(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape)
    assert d.get_filter() == shape.capitalize()


@requires_region
def test_img_get_filter_exclude(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape, ignore=True)

    expected = 'Field()&!' + shape.capitalize()
    assert d.get_filter() == expected


@requires_region
def test_img_get_filter_none(make_test_image):
    """Simple get_filter check on an image: no data"""
    d = make_test_image

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape)
    d.notice2d(ignore=True)

    # It's not clear what the filter should be here
    assert d.get_filter() == ''


@requires_region
def test_img_get_filter_combined(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape1 = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape1)

    shape2 = 'rect(4258,3830,4264,3841)'
    d.notice2d(shape2)

    shape2 = shape2.replace('rect', 'rectangle')
    shape = shape1.capitalize() + '|' + shape2.capitalize()
    assert d.get_filter() == shape


@requires_region
def test_img_get_filter_excluded(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape1 = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape1)

    shape2 = 'rect(4258,3830,4264,3841)'
    d.notice2d(shape2, ignore=True)

    shape2 = shape2.replace('rect', 'rectangle')
    shape = shape1.capitalize() + '&!' + shape2.capitalize()
    assert d.get_filter() == shape


def check_ignore_ignore(d):
    """Check removing the shapes works as expected.

    Tests must use requires_region
    """

    from sherpa.astro.utils._region import Region

    shape1 = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape1, ignore=True)

    mask1 = ~Region(shape1).mask(d.x0, d.x1).astype(bool)
    assert d.mask == pytest.approx(mask1)

    expected = 'Field()&!' + shape1.capitalize()
    assert d.get_filter() == expected

    shape2 = 'rect(4258,3830,4264,3841)'
    d.notice2d(shape2, ignore=True)

    mask2 = ~Region(shape2).mask(d.x0, d.x1).astype(bool)
    assert d.mask == pytest.approx(mask1 & mask2)

    shape2 = shape2.replace('rect', 'rectangle')
    expected = 'Field()&!' + shape1.capitalize() + '&!' + shape2.capitalize()
    assert d.get_filter() == expected


@requires_region
def test_img_get_filter_included_excluded(make_test_image):
    """Simple get_filter check on an image.

    Just to match test_img_get_filter_excluded_excluded.
    """
    d = make_test_image
    check_ignore_ignore(d)


@requires_region
def test_img_get_filter_excluded_excluded(make_test_image):
    """Simple get_filter check on an image.

    Here we want to check the behavior when d.mask is False.  I am not
    sure this makes sense, but this is done to show the current
    behavior. Note that d.notice2d(ignore=True) is meant to ignore all
    points but it (currently) doesn't add a region since we know from
    the mask all the points are ignored and so there's no need to add
    a "no region" filter: if you have !field() and then union s1 we
    would get

        !field()|s

    but this is the same as

        s

    Instead if we do !field().subtract(s) then it's the same as
    !field(). There probably is something we could improve here.

    """
    d = make_test_image

    assert d.mask
    d.notice2d(ignore=True)
    assert not d.mask

    # It is not at all obvious to me that we should get the
    # same results as test_img_get_filter_included_excluded,
    # as we start with ignoring all points.
    #
    # However, this is just to check the existing behavior,
    # which was not changed in #968.
    #
    check_ignore_ignore(d)


@requires_region
def test_img_get_filter_compare_filtering(make_test_image):
    """Check calling notice2d(ignore=True) with 2 shapes is same as once.

    """
    d = make_test_image

    shape1 = 'ellipse(4260,3840,3,2,0)'
    shape2 = 'rect(4258,3830,4264,3841)'
    d.notice2d(shape1, ignore=True)
    d.notice2d(shape2, ignore=True)
    assert d._region is not None

    maska = d.mask.copy()

    d.notice2d()
    assert d._region is None
    assert d.mask is True

    exc = f"field()-{shape1}-{shape2}"
    d.notice2d(exc)

    maskb = d.mask.copy()
    assert maskb == pytest.approx(maska)

    # just check we have some True and False values
    assert maska.min() == 0
    assert maska.max() == 1


def test_pha_size_grouped(make_grouped_pha):
    """Regression test: what is the size of this?

    Is it always DETCHANS or does it change? Test what we say.
    """

    pha = make_grouped_pha
    assert pha.quality_filter is not None
    assert pha.size == 5
    assert pha.quality_filter.sum() == 4


def test_pha_size_quality(make_quality_pha):
    """Regression test: what is the size of this?

    Is it always DETCHANS or does it change? Test what we say.
    """

    pha = make_quality_pha
    pha.ignore_bad()
    assert pha.quality_filter is not None
    assert pha.size == 9
    assert pha.quality_filter.sum() == 4


def test_grouped_quality_filter_expr(make_grouped_pha):
    """What is the filter expression?

    This is a regression test.
    """

    pha = make_grouped_pha
    pha.get_filter() == "1:9"  # does not exclude bad channel


def test_quality_quality_filter_expr(make_quality_pha):
    """What is the filter expression?

    This is a regression test.
    """

    pha = make_quality_pha
    pha.ignore_bad()

    # Note that the first group covers channels 1-4, but channels 2 and 3
    # are excluded, so this could be written as "1,4-..." but then
    # this loses the fact that the first group is 1 and 4, (so it
    # can be thought of as being correct).
    #
    pha.get_filter() == "1:9"  # does not exclude bad channels at end


def test_pha_quality_all_bad_basic_checks():
    """Regression test

    Note this only sets the quality field, not the grouping field.

    """

    pha = DataPHA("q", [1, 2, 3, 4], [9, 0, 1, 64])
    fvals = [12, 2, 7, 8]
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 4)
    assert pha.get_filter() == "1:4"
    assert pha.get_x() == pytest.approx([1, 2, 3, 4])
    assert pha.apply_filter(fvals) == pytest.approx(fvals)
    assert pha.apply_grouping(fvals) == pytest.approx(fvals)

    pha.quality = [2, 2, 2, 5]
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 4)
    assert pha.get_filter() == "1:4"
    assert pha.get_x() == pytest.approx([1, 2, 3, 4])
    assert pha.apply_filter(fvals) == pytest.approx(fvals)
    assert pha.apply_grouping(fvals) == pytest.approx(fvals)

    pha.ignore_bad()
    assert pha.mask == pytest.approx([False] * 4)
    assert pha.get_mask() == pytest.approx([False] * 4)
    assert pha.get_filter() == ""
    assert pha.get_x() == pytest.approx([1, 2, 3, 4])
    assert pha.apply_filter(fvals) == pytest.approx([])
    assert pha.apply_grouping(fvals) == pytest.approx(fvals)


@pytest.mark.parametrize("qual,fexpr,mask,counts",
                         [([2, 0, 0, 0], "1:4", [0, 1, 1, 1], [1, 64]),
                          ([0, 0, 0, 2], "1:4", [1, 1, 1, 0], [9, 1]),
                          ([0, 2, 2, 0], "1:4", [1, 0, 0, 1], [9, 64])
                         ])
def test_pha_quality_bad_range_checks(qual, fexpr, mask, counts):
    """Regression test when a group is all bad quality.

    We want to test start, end, and middle of the channels.
    """

    pha = DataPHA("q", [1, 2, 3, 4], [9, 0, 1, 64], quality=qual,
                  grouping=[1, 1, -1, 1])

    pha.ignore_bad()
    assert pha.get_filter() == fexpr
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(mask)
    assert pha.get_dep(filter=True) == pytest.approx(counts)


def test_pha_quality_change_mask(make_quality_pha):
    """A regression test."""

    pha = make_quality_pha
    pha.ignore_bad()
    assert pha.mask is True
    pha.mask = [1, 1, 0]
    assert pha.mask == pytest.approx([True, True, False])


def test_pha_quality_change_mask_ungrouped(make_quality_pha):
    """A regression test."""

    pha = make_quality_pha
    pha.ignore_bad()
    pha.ungroup()
    assert pha.mask is True
    pha.mask = [1, 1, 0, 1, 1, 0, 0, 1, 1]
    assert pha.mask == pytest.approx([True, True, False, True, True, False, False, True, True])


def test_pha_quality_change_mask_fullsize(make_quality_pha):
    """A regression test.

    Can we give the mask the "full" size (i.e. detchans)?
    No.

    """

    pha = make_quality_pha
    pha.ignore_bad()
    with pytest.raises(DataErr,
                       match="^size mismatch between grouped data and mask: 3 vs 9$"):
        pha.mask = [1, 1, 0, 1, 1, 0, 0, 1, 1]


def test_pha_quality_change_mask_wrong_size(make_quality_pha):
    """A regression test."""

    pha = make_quality_pha
    pha.ignore_bad()
    with pytest.raises(DataErr,
                       match="^size mismatch between grouped data and mask: 3 vs 5$"):
        pha.mask = [1, 1, 0, 1, 1]


def test_pha_quality_change_mask_ungrouped_wrong_size(make_quality_pha):
    """A regression test."""

    pha = make_quality_pha
    pha.ignore_bad()
    pha.ungroup()
    with pytest.raises(DataErr,
                       match="^size mismatch between independent axis and mask: 9 vs 5$"):
        pha.mask = [1, 1, 0, 1, 1]


def test_pha_change_channels(make_test_pha):
    """What happens if we change the channel/count values?

    We have several ways of getting the independent and dependent
    axes.
    """
    pha = make_test_pha

    channels = [1, 2, 3, 4]
    counts = [1, 2, 0, 3]
    assert np.all(pha.channel == channels)
    assert np.all(pha.counts == counts)

    assert len(pha.get_indep()) == 1
    assert np.all(pha.get_indep()[0] == channels)
    assert np.all(pha.get_dep() == counts)

    assert len(pha.indep) == 1
    assert np.all(pha.indep[0] == channels)
    assert np.all(pha.dep == counts)

    channels2 = [2, 3, 4, 5]
    counts2 = [20, 30, 20, 10]
    pha.channel = channels2
    pha.counts = counts2
    assert np.all(pha.channel == channels2)
    assert np.all(pha.counts == counts2)

    assert len(pha.get_indep()) == 1
    assert np.all(pha.get_indep()[0] == channels2)
    assert np.all(pha.get_dep() == counts2)

    assert len(pha.indep) == 1
    assert np.all(pha.indep[0] == channels2)
    assert np.all(pha.dep == counts2)


def test_pha_add_channels(make_test_pha):
    """What happens if we increase the number of channels/counts?

    Extends test_pha_change_channels
    """
    pha = make_test_pha

    channels2 = np.arange(1, 6, dtype=int)
    with pytest.raises(DataErr,
                       match="independent axis can not change size: 4 to 5"):
        pha.channel = channels2


def test_pha_remove_channels(make_test_pha):
    """What happens if we decrease the number of channels/counts?

    Extends test_pha_change_channels
    """
    pha = make_test_pha

    channels2 = np.arange(1, 4, dtype=int)
    with pytest.raises(DataErr,
                       match="independent axis can not change size: 4 to 3"):
        pha.channel = channels2


@pytest.mark.parametrize("requested,expected",
                         [("bin", "channel"), ("Bin", "channel"),
                          ("channel", "channel"), ("ChannelS", "channel"),
                          ("chan", "channel"),
                          ("energy", "energy"), ("ENERGY", "energy"),
                          ("Energies", "energy"),
                          ("WAVE", "wavelength"), ("wavelength", "wavelength"),
                          ("Wavelengths", "wavelength"),
                          ("chan This Is Wrong", "channel"),  # should this be an error?
                          ("WAVEY GRAVY", "wavelength")  # should this be an error?
                          ])
def test_pha_valid_units(requested, expected, make_test_pha):
    """Check we can set the units field of a PHA object"""
    pha = make_test_pha
    pha.units = requested
    assert pha.units == expected


@pytest.mark.parametrize("invalid", ["Bins", "BINNING", "wavy", "kev", "angstrom"])
def test_pha_invalid_units(invalid, make_test_pha):
    """Check we can not set units to an invalid value"""
    pha = make_test_pha
    with pytest.raises(DataErr,
                       match=f"unknown quantity: '{invalid}'"):
        pha.units = invalid


@pytest.mark.parametrize("invalid", ["RATE", "COUNTS", "rates", "count", "count-rate"])
def test_pha_analysis_type_invalid(invalid, make_test_pha):
    pha = make_test_pha
    with pytest.raises(DataErr,
                       match=f"unknown plot type '{invalid}', choose 'rate' or 'counts'"):
        pha.set_analysis("channel", type=invalid)


def test_pha_analysis_plot_fac_valid(make_test_pha):
    """Historically we've allowed 2.0 as an argument, so check it still works"""
    pha = make_test_pha
    assert pha.plot_fac == 0
    pha.plot_fac = 2.0
    assert pha.plot_fac == 2


@pytest.mark.parametrize("invalid", ["1", 2.01, 0.5, complex(1)])
def test_pha_analysis_plot_fac_invalid(invalid, make_test_pha):
    pha = make_test_pha
    # Need to protect the '(1+0j)' brackets.
    #
    emsg = re.escape(f"unknown plot_fac setting: '{invalid}'")
    with pytest.raises(DataErr,
                       match=emsg):
        pha.plot_fac = invalid


@pytest.mark.parametrize("invalid", ["1", 2.01, 0.5, complex(1)])
def test_pha_analysis_factor_invalid(invalid, make_test_pha):
    pha = make_test_pha
    # Need to protect the '(1+0j)' brackets.
    #
    emsg = re.escape(f"unknown factor setting: '{invalid}'")
    with pytest.raises(DataErr,
                       match=emsg):
        pha.set_analysis("channel", factor=invalid)


def test_pha_get_specresp_no_response(make_test_pha):
    pha = make_test_pha
    assert pha.get_specresp() is None


def test_pha_ignore_bad_no_quality(make_test_pha):
    pha = make_test_pha
    assert pha.quality is None
    with pytest.raises(DataErr,
                       match="data set 'p' does not specify quality flags"):
        pha.ignore_bad()


def test_pha_quality_noticed_channels_no_filter(make_quality_pha):
    """Regression test."""

    pha = make_quality_pha
    chans0 = pha.get_noticed_channels()
    assert chans0 == pytest.approx(np.arange(1, 10))

    pha.ignore_bad()

    chans1 = pha.get_noticed_channels()
    assert chans1 == pytest.approx([1, 4, 5, 6])


def test_pha_quality_noticed_channels_with_filter(make_quality_pha):
    """Regression test."""

    pha = make_quality_pha
    assert pha.grouped
    pha.ignore_bad()

    # The bins are 1-4 (although 2 and 3 are excluded), 5, 6, 7-9 (all
    # excluded). Ignoring at hi=2 makes this interesting, as the first
    # group should then be removed.
    #
    pha.ignore(hi=2)

    chans1 = pha.get_noticed_channels()
    assert chans1 == pytest.approx([5, 6])

    pha.notice(lo=2)
    chans2 = pha.get_noticed_channels()
    assert chans2 == pytest.approx([1, 4, 5, 6])


def test_pha_quality_ignore_bad_clear_filter(make_quality_pha):
    """Regression test."""

    pha = make_quality_pha

    assert pha.get_filter() == "1:9"
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 9)
    assert pha.quality_filter is None

    # channels 2,3 and 7-9 are "bad"
    pha.ignore(hi=3)

    assert pha.get_filter() == "5:9"
    assert pha.mask == pytest.approx([False] + [True] * 3)
    assert pha.get_mask() == pytest.approx([False] * 4 + [True] * 5)
    assert pha.quality_filter is None

    # This resets the previous filters
    pha.ignore_bad()

    qflags = [True] * 1 + [False] * 2 + [True] * 3 + [False] * 3
    assert pha.get_filter() == "1:9"
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(qflags)
    assert pha.quality_filter == pytest.approx(qflags)

    pha.ignore(hi=3)

    assert pha.get_filter() == "5:6"
    assert pha.mask == pytest.approx([False] + [True] * 2)
    assert pha.get_mask() == pytest.approx([False] * 2 + [True] * 2)
    assert pha.quality_filter == pytest.approx(qflags)

    pha.ignore(lo=2, hi=4)

    assert pha.get_filter() == "5:6"
    assert pha.mask == pytest.approx([False] + [True] * 2)
    assert pha.get_mask() == pytest.approx([False] * 2 + [True] * 2)
    assert pha.quality_filter == pytest.approx(qflags)

    # This removes the quality filter!
    pha.notice()

    assert pha.get_filter() == "1:9"
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 9)
    assert pha.quality_filter is None


def test_pha_grouping_changed_no_filter_1160(make_test_pha):
    """What happens when the grouping is changed?

    See also test_pha_grouping_changed_filter_1160
    """

    pha = make_test_pha
    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([1, 2, 0, 3])

    # grouping set but not grouped
    pha.grouping = [1, 1, 1, 1]
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx([1, 2, 0, 3])

    # now grouped
    pha.grouped = True
    d3 = pha.get_dep(filter=True)
    assert d3 == pytest.approx([1, 2, 0, 3])

    pha.grouping = [1, 1, -1, 1]
    d4 = pha.get_dep(filter=True)
    assert d4 == pytest.approx([1, 2, 3])


def test_pha_grouping_changed_filter_1160(make_test_pha):
    """What happens when the grouping is changed?

    See also test_pha_grouping_changed_no_filter_1160
    """

    pha = make_test_pha
    pha.notice(2, 5)

    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([2, 0, 3])

    # grouping set but not grouped
    pha.grouping = [1, 1, 1, 1]
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx([2, 0, 3])

    # now grouped
    pha.grouped = True
    d3 = pha.get_dep(filter=True)
    assert d3 == pytest.approx([2, 0, 3])

    pha.grouping = [1, 1, -1, 1]
    d4 = pha.get_dep(filter=True)
    assert d4 == pytest.approx([2, 3])


def test_pha_grouping_changed_1160_grped_no_filter(make_grouped_pha):
    """Test based on work on #1160

    This is probably no different to
    test_pha_grouping_changed_no_filter_1160 and
    test_pha_grouping_changed_filter_1160 but separated out
    as more tests here are probably useful.
    """

    # Do we care about adding a response?
    pha = make_grouped_pha

    # why does this not understand the "bad quality" filter?
    ofilter = "1:5"
    assert pha.get_filter() == ofilter

    # Change the grouping
    pha.grouping = [1] * 5

    # Although no grouping, we still have the bad filter in place
    assert pha.get_dep(filter=False) == pytest.approx([1, 2, 0, 3, 12])
    assert pha.get_dep(filter=True) == pytest.approx([1, 2, 0, 3])

    assert pha.get_filter() == ofilter


def test_pha_grouping_changed_1160_grped_with_filter(make_grouped_pha):
    """Test based on work on #1160

    See test_pha_grouping_changed_1160_grped_no_filter
    """

    pha = make_grouped_pha

    # Can not say
    #   pha.notice(2, 4)
    # as we have already done a filter, so this would not
    # act to ignore the first channel.
    #
    # The ignore(lo=5) line is not needed as it is already excluded by
    # a bad-quality channel, but users can say this, so check the the
    # response.
    #
    pha.ignore(hi=1)
    pha.ignore(lo=5)

    # Dropping channel 1 means the first group gets dropped, so we
    # only have channel 4 left.
    #
    assert pha.get_filter() == "4"

    assert pha.get_dep(filter=False) == pytest.approx([1, 2, 0, 3, 12])
    assert pha.get_dep(filter=True) == pytest.approx([3])

    # Change the grouping; it would be nice if it could have
    # recognized the requested range was > 1 and <= 5 but the current
    # code does not support this.
    #
    pha.grouping = [1] * 5

    assert pha.get_dep(filter=False) == pytest.approx([1, 2, 0, 3, 12])
    assert pha.get_dep(filter=True) == pytest.approx([3])

    assert pha.get_filter() == "4"


def test_pha_grouping_changed_1160_ungrped_with_filter(make_grouped_pha):
    """Test based on work on #1160

    A version of
    test_pha_grouping_changed_1160_grped_with_filter
    but the data is not grouped, even though the grouping
    is set/changed.

    """

    pha = make_grouped_pha

    # Apply the filter whilst still grouped
    pha.ignore(hi=1)
    pha.ignore(lo=5)

    pha.ungroup()

    # The filtering does not change because of the ungroup call,
    # although we might like it too.
    assert pha.get_filter() == "4"

    assert pha.get_dep(filter=False) == pytest.approx([1, 2, 0, 3, 12])
    assert pha.get_dep(filter=True) == pytest.approx([3])

    # Change the grouping
    pha.grouping = [1] * 5

    assert pha.get_dep(filter=False) == pytest.approx([1, 2, 0, 3, 12])
    assert pha.get_dep(filter=True) == pytest.approx([3])

    assert pha.get_filter() == "4"


@requires_fits
@requires_data
def test_1160(make_data_path):
    """Use the dataset we reported this with just as an extra check

    It is slightly different to the other #1160 tests above because

    a) we do not have a non-zero quality bin
    b) we have an instrument response (this should not really
       matter here)

    """

    from sherpa.astro.io import read_pha
    pha = read_pha(make_data_path("3c273.pi"))

    fexpr = "0.47:6.57"

    pha.notice(0.5, 6)
    assert pha.get_dep(filter=False).shape == (1024, )
    assert pha.get_dep(filter=True).shape == (41, )
    assert pha.mask.shape == (46, )
    assert pha.mask.sum() == 41
    assert pha.get_filter(format="%.2f") == fexpr

    pha.grouping = [1] * 1024
    assert pha.get_dep(filter=False).shape == (1024, )
    assert pha.get_dep(filter=True).shape == (418, )
    assert pha.mask.shape == (1024, )
    assert pha.mask.sum() == 418
    assert pha.get_filter(format="%.2f") == fexpr


def test_pha_remove_grouping(make_test_pha):
    """Check we can remove the grouping array.

    See issue #1659
    """

    pha = make_test_pha
    assert pha.grouping is None
    assert not pha.grouped

    no_data = [1, 2, 0, 3]
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx(no_data)

    pha.grouping = [1, -1, 1, -1]
    assert not pha.grouped
    pha.grouped = True
    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([3, 3])

    # Can we remove the grouping column?
    pha.grouping = None

    # Check the get_dep behavior before grouped as this is currently
    # causes a TypeError rather than not changing a boolean flag
    # variable.
    #
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx(no_data)

    # This thinks that pha.grouped is still set
    assert not pha.grouped


@pytest.mark.xfail
@pytest.mark.parametrize("grouping", [True, [1, 1], np.ones(10)])
def test_pha_grouping_size(grouping, make_test_pha):
    """Check we error out if grouping has the wrong size"""

    pha = make_test_pha
    with pytest.raises(DataErr) as de:
        pha.grouping = grouping

    assert str(de.value) == 'size mismatch between channel and grouping'


def test_pha_remove_quality(make_test_pha):
    """Check we can remove the quality array."""

    pha = make_test_pha
    assert pha.quality is None

    no_data = [1, 2, 0, 3]
    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx(no_data)

    pha.quality = [0, 0, 0, 2]
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx(no_data)

    pha.quality = None
    d3 = pha.get_dep(filter=True)
    assert d3 == pytest.approx(no_data)


@pytest.mark.xfail
def test_pha_remove_quality_bad(make_test_pha):
    """Check we can remove the quality array after calling ignore_bad

    Here we ensure we have a "bad" value that will be
    marked bad by ignore_bad.

    What is the expected behavior after removing the
    quality array? See #1427
    """

    pha = make_test_pha
    assert pha.quality is None

    no_data = [1, 2, 0, 3]

    pha.quality = [0, 0, 0, 2]
    pha.ignore_bad()
    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([1, 2, 0])

    # At the moment d2 == [1, 2, 0] so the quality filter remains
    pha.quality = None
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx(no_data)


def test_pha_quality_bad_filter(make_test_pha, caplog):
    """What is the filter expression when ignore bad + filter

    Also check the screen output (there is none for the PHA case,
    unlike the UI version).

    """

    pha = make_test_pha
    assert pha.get_filter() == "1:4"

    assert len(caplog.record_tuples) == 0
    with SherpaVerbosity("INFO"):
        pha.ignore(hi=1)

    assert pha.get_filter() == "2:4"
    assert len(caplog.record_tuples) == 0

    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([2, 0, 3])

    pha.quality = [0, 0, 0, 2]
    with SherpaVerbosity("INFO"):
        pha.ignore_bad()

    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx([2, 0])
    assert pha.get_filter() == "2:3"
    assert len(caplog.record_tuples) == 0


def test_pha_quality_bad_filter2(make_quality_pha, caplog):
    """A different set of bad channels and groups to test_pha_quality_bad_filter

    This also does not include a filter, since this messes everything up at
    this time.
    """

    pha = make_quality_pha
    assert pha.get_filter() == "1:9"
    assert pha.grouped

    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([6, 12, 2, 24])

    with SherpaVerbosity("INFO"):
        pha.ignore_bad()

    assert len(caplog.record_tuples) == 0

    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx([4, 12, 2])
    assert pha.get_filter() == "1:9"

    assert len(caplog.record_tuples) == 0


@pytest.mark.xfail
def test_pha_quality_bad_filter_remove(make_test_pha):
    """test_pha_quality_bad_filter then remove the quality array

    What is the expected behavior after removing the
    quality array? See #1427
    """

    pha = make_test_pha
    pha.ignore(hi=1)
    pha.quality = [0, 0, 0, 2]
    pha.ignore_bad()

    # At the moment the filter still includes the quality filter
    pha.quality = None
    assert pha.get_filter() == "2:4"


@pytest.mark.parametrize("field,expected",
                         [("channel", [1, 2, 3, 4, 5, 6, 7, 8, 9]),
                          ("counts", [1, 2, 0, 3, 12, 2, 9, 8, 7]),
                          ("grouping", [1, -1, -1, -1, 1, 1, 1, -1, -1]),
                          ("quality", [0, 5, 5, 0, 0, 0, 2, 2, 5])
                          ])
def test_pha_quality_bad_field(field, expected, make_quality_pha):
    """After ignore_bad what does the field return?"""

    pha = make_quality_pha
    assert getattr(pha, field) == pytest.approx(expected)

    pha.ignore_bad()
    assert getattr(pha, field) == pytest.approx(expected)


def test_pha_quality_bad_mask(make_quality_pha):
    """What does the mask look like?"""

    pha = make_quality_pha
    assert pha.mask is True

    pha.ignore_bad()
    assert pha.mask is True


def test_pha_quality_bad_mask_grouped(make_quality_pha):
    """What does the mask look like?"""

    pha = make_quality_pha
    pha.ignore_bad()
    pha.group()
    assert pha.mask is True


def test_pha_quality_bad_get_mask(make_quality_pha):
    """What does the get_mask() look like?"""

    pha = make_quality_pha
    assert pha.get_mask() == pytest.approx([1] * 9)

    pha.ignore_bad()
    assert pha.get_mask() == pytest.approx([1, 0, 0, 1, 1, 1, 0, 0, 0])


def test_pha_no_quality_ignore_bad(make_test_pha):
    """What happens if call ignore_bad and no quality data"""

    pha = make_test_pha
    with pytest.raises(DataErr,
                       match="^data set 'p' does not specify quality flags$"):
        pha.ignore_bad()


@requires_group
def test_pha_change_quality_values(caplog):
    """What happens if we change the quality column?

    This is a regression test as it is likely we should change the filter,
    but we have not thought through the consequences. See also #1427
    """

    pha = DataPHA('ex', [1, 2, 3, 4, 5, 6, 7], [1, 2, 1, 0, 2, 2, 1])
    pha.group_counts(5)
    assert pha.quality == pytest.approx([0, 0, 0, 0, 0, 2, 2])
    assert pha.get_dep(filter=True) == pytest.approx([6, 3])
    assert pha.get_filter() == '1:7'

    assert pha.quality_filter is None
    assert len(caplog.records) == 0

    pha.ignore_bad()
    assert len(caplog.records) == 0

    assert pha.quality_filter == pytest.approx([True] * 5 + [False] * 2)
    assert pha.get_dep(filter=True) == pytest.approx([6])
    assert pha.get_filter() == '1:7'

    # With no tabStops set it uses ~pha.get_mask() which in this case
    # is [False] * 5 + [True] * 2,
    #
    pha.group_counts(4)
    assert len(caplog.records) == 0

    assert pha.quality == pytest.approx([0, 0, 0, 2, 2, 0, 0])

    # Should quality filter be reset?
    assert pha.quality_filter == pytest.approx([True] * 5 + [False] * 2)
    assert pha.get_dep(filter=True) == pytest.approx([4, 2])
    assert pha.get_filter() == '1:7'


def test_pha_group_adapt_check():
    """Regression test.

    This was found when investigating ignore_bad issues and felt to
    be worth a test to check how this code behaves.
    """

    counts = [4, 2, 3, 1, 5, 6, 7]
    pha = DataPHA("ex", [1, 2, 3, 4, 5, 6, 7], counts)
    pha.group_adapt(6)

    # The grouping may change if the adaptive scheme changes.
    assert pha.grouping == pytest.approx([1, -1, 1, 1, -1, 1, 1])

    # The group library behaves oddly (the last element being 2).
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 2])

    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 3, 6, 6, 7])


def test_pha_group_ignore_bad_then_filter(caplog):
    """Regression test."""

    counts = [4, 2, 3, 1, 5, 6, 7]
    pha = DataPHA("ex", [1, 2, 3, 4, 5, 6, 7], counts)

    # The equivalent of pha.group_adapt(6) with CIAO 4.17 but this
    # may change, so set the data manually. Compare to
    # test_pha_group_adapt_check
    #
    # pha.group_adapt(6)
    pha.grouping = [1, -1, 1, 1, -1, 1, 1]
    pha.quality = [0, 0, 2, 0, 0, 0, 2]
    pha.group()
    assert len(caplog.records) == 0

    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 7)
    assert pha.get_filter() == '1:7'
    assert pha.quality_filter is None

    pha.ignore_bad()
    assert len(caplog.records) == 0

    qual_mask = [True] * 2 + [False] + [True] * 3 + [False]
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(qual_mask)
    assert pha.get_filter() == '1:7'
    assert pha.quality_filter == pytest.approx(qual_mask)
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 2])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 6, 6])

    pha.ignore(4, 5)
    assert len(caplog.records) == 0

    assert pha.mask == pytest.approx([True, False, True])
    assert pha.get_mask() == pytest.approx([True] * 2 + [False] * 2 + [True])
    assert pha.get_filter() == '1:2,6'
    assert pha.quality_filter == pytest.approx(qual_mask)
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 2])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 6])


def test_pha_group_ignore_bad_then_group(caplog):
    """Regression test."""

    counts = [4, 2, 3, 1, 5, 6, 7]
    pha = DataPHA("ex", [1, 2, 3, 4, 5, 6, 7], counts)
    pha.group_adapt(6)
    pha.ignore_bad()
    assert len(caplog.records) == 0

    qual_mask = [True] * 2 + [False] + [True] * 3 + [False]
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(qual_mask)
    assert pha.quality_filter == pytest.approx(qual_mask)

    # Change the grouping. What happens with the existing "bad
    # quality" data?
    #
    pha.group_counts(4)
    assert len(caplog.records) == 0

    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(qual_mask)
    assert pha.get_filter() == '1:7'
    assert pha.quality_filter == pytest.approx(qual_mask)
    assert pha.quality == pytest.approx([0, 2, 0, 0, 0, 0, 0])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([4, 2, 6, 6])

    # Shouldn't this be a no-op. It isn't because the group call
    # didn't change the quality_filter array, so it now changes what
    # are the good/bad channels.
    #
    pha.ignore_bad()
    assert len(caplog.records) == 0

    qual_mask = [True] + [False] + [True] * 5
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(qual_mask)
    assert pha.get_filter() == '1:7'
    assert pha.quality_filter == pytest.approx(qual_mask)
    assert pha.quality == pytest.approx([0, 2, 0, 0, 0, 0, 0])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([4, 3, 6, 6, 7])


def test_pha_filter_ignore_bad_filter(caplog):
    """A regression test.

    Mix filtering, ignore_bad, and more filtering.
    """

    counts = np.asarray([4, 2, 3, 1, 5, 6, 7])
    pha = DataPHA("ex", [1, 2, 3, 4, 5, 6, 7], counts)

    pha.ignore(lo=4, hi=4)
    assert len(caplog.records) == 0

    data_mask = [True] * 3 + [False] + [True] * 3
    assert pha.mask == pytest.approx(data_mask)
    assert pha.get_mask() == pytest.approx(data_mask)
    assert pha.get_filter() == '1:3,5:7'
    assert pha.quality_filter is None
    assert pha.quality is None
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx(counts[[0, 1, 2, 4, 5, 6]])

    pha.group_counts(5)
    assert len(caplog.records) == 0

    assert pha.mask == pytest.approx([True] * 2 + [False] + [True] * 3)
    assert pha.get_mask() == pytest.approx(data_mask)
    assert pha.get_filter() == '1:3,5:7'
    assert pha.quality_filter is None
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 0])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 3, 5, 6, 7])

    pha.ignore_bad()
    assert len(caplog.records) == 1

    r = caplog.records[-1]
    assert r.name == "sherpa.astro.data"
    assert r.levelname == "WARNING"
    assert r.getMessage() == "filtering grouped data with quality flags, previous filters deleted"

    new_mask = [True] * 2 + [False] + [True] * 4
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx(new_mask)
    assert pha.get_filter() == '1:7'
    assert pha.quality_filter == pytest.approx(new_mask)
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 0])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 1, 5, 6, 7])

    pha.ignore(lo=2, hi=2)
    assert len(caplog.records) == 1

    assert pha.mask == pytest.approx([False] + [True] * 4)
    assert pha.get_mask() == pytest.approx([False] * 2 + [True] * 4)
    assert pha.get_filter() == '4:7'
    assert pha.quality_filter == pytest.approx(new_mask)
    assert pha.quality == pytest.approx([0, 0, 2, 0, 0, 0, 0])
    assert pha.get_dep(filter=False) == pytest.approx(counts)
    assert pha.get_dep(filter=True) == pytest.approx([1, 5, 6, 7])


@pytest.mark.parametrize("field", ["grouping", "quality"])
def test_pha_change_xxx_non_integer_value(field, make_test_pha):
    """What happens if send grouping/quality values that can not be converted to an array?"""

    pha = make_test_pha
    invalid = [None, "x", {}, set()]
    with pytest.raises(DataErr,
                       match="Array must be a sequence of integers or None"):
        setattr(pha, field, invalid)


def test_pha_change_grouping_type(make_test_pha):
    """Check the grouping column is converted to int"""
    pha = make_test_pha
    grp = np.asarray([1.0, -1.0, -1.0, 1.0])
    pha.grouping = grp

    # Since integer values can do an exact check
    assert (pha.grouping == np.asarray([1, -1, -1, 1])).all()
    assert pha.grouping.dtype == np.int16


def test_pha_change_quality_type(make_test_pha):
    """Check the quality column is converted to int"""
    pha = make_test_pha
    # technically negative numbers are allowed
    qual = np.asarray([0.0, 2.0, 5.0, -1.0])
    pha.quality = qual

    # Since integer values can do an exact check
    assert (pha.quality == np.asarray([0, 2, 5, -1])).all()
    assert pha.quality.dtype == np.int16


@pytest.mark.parametrize("label", ["grouping", "quality"])
def test_pha_change_grouping_rounding(label, make_test_pha):
    """What happens with non-integer values?

    Unlike test_pha_change_grouping/quality_type we can more-easily
    use the same input array, which makes it easier to test both
    columns with the same routine. It is actually unclear what
    we should do with input like this - should we error out,
    silently truncate, or perhaps warn the user. For the moment
    test assuming silent truncation.

    """

    pha = make_test_pha
    vals = np.asarray([0.5, 1.5, -0.5, 0.9])
    setattr(pha, label, vals)

    got = getattr(pha, label)
    assert (got == np.asarray([0, 1, 0, 0])).all()


@requires_group
def test_pha_ignore_bad_group_quality(caplog):
    """Check handling of ignore_bad when quality and grouping set.

    This used to be called test_416_b but has been expanded to
    check a few more things. See also
    test_pha_ignore_bad_quality which is meant to
    be the same but with an ungrouped dataset (so the
    results won't quite match).

    """

    # The energy range matches the channel values to make
    # things easier.
    #
    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 1, 0])

    pha = DataPHA('416', x, y)

    rmf = create_delta_rmf(x, x + 1, e_min=x, e_max=x + 1,
                           name='416')
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    assert pha.get_filter(format="%.1f") == "1.0:11.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(1, 11))

    # No grouping or filtering yet
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y)

    assert not pha.grouped

    # After this we have
    # - two groups, channels 1-5 and 6-10
    # - the first group has quality=0, the second quality=2
    # - the noticed range is channels 3-7 before grouping
    #   which becomes 1-11 after grouping (i.e. all points)
    #
    pha.notice(3.5, 6.5)
    assert pha.get_filter(format="%.1f") == "3.0:7.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(3, 7))

    omask = [False] * 2 + [True] * 4 + [False] * 4
    assert pha.mask == pytest.approx(omask)
    assert pha.get_mask() == pytest.approx(omask)

    # Only filtering
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y[2:6])

    # tabStops is not set, so uses the current mask
    pha.group_counts(3)
    assert pha.get_filter(format="%.1f") == "3.0:7.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(3, 7))

    # Grouped and filtered
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([3, 1])

    assert pha.mask == pytest.approx([False] * 2 + [True] * 2 + [False] * 4)
    assert pha.get_mask() == pytest.approx(omask)

    grouping = [0, 0, 1, -1, -1,  1, 0, 0, 0, 0]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 0, 0, 0, 0]
    assert pha.quality == pytest.approx(quality)
    assert pha.quality_filter is None

    assert pha.grouped

    # By calling ignore_bad we have
    # - removed the channels with quality=2, which is
    #   channels 6-10
    # - removed the noticed range
    #
    assert len(caplog.record_tuples) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        pha.ignore_bad()

    # check captured log
    #
    emsg = 'filtering grouped data with quality flags, previous filters deleted'
    assert caplog.record_tuples == [
        ('sherpa.astro.data', logging.WARNING, emsg)
        ]

    assert pha.grouped

    # We have reverted the energy filter, so the mask attribute
    # is back to a boolean.
    #
    assert type(pha.mask) is bool
    assert pha.mask is True

    # However, get_mask reflects the quality filter, so is all True
    # except for the 6th element.
    #
    single_bad = [True] * 5 + [False] + [True] * 4
    assert pha.get_mask() == pytest.approx(single_bad)

    # What about the quality fields?
    #
    assert pha.quality == pytest.approx(quality)
    assert pha.quality_filter == pytest.approx(single_bad)

    # Saying all that though, the filter expression does not
    # know we are ignoring channels 6-10.
    #
    # TODO: This is likely a bug.
    #
    assert pha.get_filter(format="%.1f") == "1.0:11.0"
    assert pha.get_noticed_channels() == pytest.approx([1, 2, 3, 4, 5, 7, 8, 9, 10])

    # Grouped and quality-filtered (even though get_filter
    # returns 1:11 here).
    #
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([0, 0, 3, 0, 0, 1, 0])

    # check there have been no more messages.
    #
    assert len(caplog.record_tuples) == 1


@pytest.mark.parametrize("groupit", [False, True])
def test_pha_ignore_bad_quality(groupit, caplog):
    """Check handling of ignore_bad when quality set but no grouping.

    See test_pha_ignore_bad_group_quality. The case when
    the quality array is not set is handled earlier by
    test_pha_ignore_bad_no_quality

    The groupit flag is used to ensure the results are
    the same if the data has no grouping data at all
    (False) or has grouping but is not used (True).

    """

    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 1, 0])

    pha = DataPHA('416', x, y)

    rmf = create_delta_rmf(x, x + 1, e_min=x, e_max=x + 1,
                           name='416')
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    grps = np.asarray([1, -1, -1, -1, -1] * 2)
    if groupit:
        pha.grouping = grps

    assert not pha.grouped

    assert pha.get_filter(format="%.1f") == "1.0:11.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(1, 11))

    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y)

    # After this we have
    # - the noticed range is channels 3-7
    #
    pha.notice(3.5, 6.5)
    assert pha.get_filter(format="%.1f") == "3.0:7.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(3, 7))

    # Only filtering
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y[2:6])

    mask = [False] * 2 + [True] * 4 + [False] * 4
    assert pha.mask == pytest.approx(mask)
    assert pha.get_mask() == pytest.approx(mask)

    if groupit:
        assert pha.grouping == pytest.approx(grps)
    else:
        assert pha.grouping is None

    assert not pha.grouped
    assert pha.quality is None
    assert pha.quality_filter is None

    # Now apply quality filtering without grouping. We choose
    # the same quality range as test_pha_grouped_filtered_quality_warns
    #
    quality = [0] * 5 + [2] * 5
    pha.quality = quality
    assert pha.quality == pytest.approx(quality)
    assert pha.quality_filter is None

    # By calling ignore_bad we have
    # - removed the channels with quality=2, which is
    #   channels 6-10
    #
    assert len(caplog.record_tuples) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        pha.ignore_bad()

    assert not pha.grouped

    # check captured log; at the moment this DOES NOT warn the
    # user about the filter being removed.
    #
    assert len(caplog.record_tuples) == 0

    # The mask changed (the channel=6 value is now filtered out).
    #
    mask2 = [False] * 2 + [True] * 3 + [False] * 5
    assert pha.mask == pytest.approx(mask2)
    assert pha.get_mask() == pytest.approx(mask2)

    # What about the quality fields?
    #
    assert pha.quality == pytest.approx(quality)
    assert pha.quality_filter is None

    # The filter expression has changed to reflect the quality filter;
    # this is unlike the grouped version above.
    #
    assert pha.get_filter(format="%.1f") == "3.0:6.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(3, 6))

    # noticed and quality-filtered.
    #
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y[2:5])

    # check there have been no more messages.
    #
    assert len(caplog.record_tuples) == 0


@requires_group
def test_361():
    """Check issue #361

    This is also tested in test_filter_bad_notice_361 in
    sherpa/astro/ui/tests/test_filtering.py using the UI
    interface.
    """

    # energy ranges are
    #   0.1-0.2, 0.2-0.3, ..., 1.0-1.1
    # and when grouped we get
    #   0.1-0.3, 0.3-0.5, 0.5-0.7, 0.7-0.9, 0.9-1.1
    # with counts
    #   12, 6, 11, 8, 3
    # and then the quality array knocks out the
    #   0.5-0.7 group (11 counts).
    #
    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([5, 7, 2, 4, 6, 5, 8, 0, 1, 2])
    grp = np.asarray([1, -1] * 5)
    qual = np.zeros(10)
    qual[4:6] = 2

    pha = DataPHA('361', x, y,
                  grouping=grp, quality=qual)

    elo = x * 0.1
    ehi = (x + 1) * 0.1
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi,
                           name='4361')
    pha.set_arf(rmf)
    pha.set_analysis('energy')

    assert pha.grouped
    assert pha.get_dep() == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([12, 6, 11, 8, 3])
    assert pha.get_noticed_channels() == pytest.approx(np.arange(1, 11))

    pha.ignore_bad()
    assert pha.get_dep() == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([12, 6, 8, 3])
    assert pha.get_noticed_channels() == pytest.approx([1, 2, 3, 4, 7, 8, 9, 10])

    pha.notice(0.35, 0.8)
    assert pha.get_dep() == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([6, 8])

    # The issue in #361 seems to come from evaluating an array
    # of the expected length as created by the model. We can
    # be more-direct here and check the problematic call.
    #
    assert pha.get_noticed_channels() == pytest.approx([3, 4, 7, 8])


def test_grouped_pha_get_dep(make_grouped_pha):
    """Quality filtering and grouping is applied: get_dep

    As noted in issue #1438 it's not obvious what get_y is meant to
    return. It is not the same as get_dep as there's post-processing.
    So just test the current behavior.

    """
    pha = make_grouped_pha

    # grouped counts are [3, 3, 12]
    # channel widths are [3, 1, 1]
    # which gives [1, 3, 12]
    # but the last group is marked bad by quality,
    # so we expect [1, 3]
    #
    assert pha.get_dep() == pytest.approx([1, 3])


def test_grouped_pha_get_y(make_grouped_pha):
    """Quality filtering and grouping is applied: get_y

    As noted in issue #1438 it's not obvious what get_y is meant to
    return. It is not the same as get_dep as there's post-processing.
    So just test the current behavior.

    """
    pha = make_grouped_pha

    # grouped counts are [3, 3, 12]
    # channel widths are [3, 1, 1]
    # which gives [1, 3, 12]
    # but the last group is marked bad by quality,
    # so we expect [1, 3]
    #
    assert pha.get_y() == pytest.approx([1, 3])


def test_quality_pha_get_dep(make_quality_pha):
    """Regression test for quality + filtering."""

    pha = make_quality_pha

    # ungrouped counts no quality are [1, 2, 0, 3, 12, 2, 9, 8, 7]
    # ungrouped counts with quality are [1, 3, 12, 2]
    #
    # grouped counts no quality are [6, 12, 2, 24]
    # grouped counts with quality are [4, 12, 2]
    #
    all_counts = [1, 2, 0, 3, 12, 2, 9, 8, 7]
    assert pha.get_dep(filter=False) == pytest.approx(all_counts)
    assert pha.get_dep(filter=True) == pytest.approx([6, 12, 2, 24])

    pha.ignore_bad()
    assert pha.get_dep(filter=False) == pytest.approx(all_counts)
    assert pha.get_dep(filter=True) == pytest.approx([4, 12, 2])


def test_quality_pha_get_y(make_quality_pha):
    """Regression test for quality + filtering."""

    pha = make_quality_pha

    # bin widths are
    #    no quality  [4, 1, 1, 3]
    #       quality  [2, 1, 1]
    #
    # grouped counts no quality are [6, 12, 2, 24]
    # grouped counts with quality are [4, 12, 2]
    #
    expected1 = np.asarray([1.5, 12, 2, 8])
    assert pha.get_y(filter=False) == pytest.approx(expected1)
    assert pha.get_y(filter=True) == pytest.approx(expected1)

    pha.ignore_bad()
    # expected2 = np.asarray([2, 12, 2])
    expected2 = np.asarray([1, 12, 2])  # TODO: why do we get this?
    assert pha.get_y(filter=False) == pytest.approx(expected2)
    assert pha.get_y(filter=True) == pytest.approx(expected2)


def test_quality_pha_get_indep(make_quality_pha):
    """Regression test for quality + filtering.

    This ignores the channel settings.
    """

    pha = make_quality_pha
    assert pha.units == "channel"

    chans = np.arange(1, 10)

    for filt in [False, True]:
        indep = pha.get_indep(filter=filt)
        assert len(indep) == 1
        assert indep[0] == pytest.approx(chans)

    pha.ignore_bad()

    indep = pha.get_indep(filter=False)
    assert len(indep) == 1
    assert indep[0] == pytest.approx(chans)

    indep = pha.get_indep(filter=True)
    assert len(indep) == 1
    assert indep[0] == pytest.approx([1, 4, 5, 6])


def test_grouped_pha_mask(make_grouped_pha):
    """What is the default mask setting?"""
    pha = make_grouped_pha
    assert np.isscalar(pha.mask)
    assert pha.mask is True


def test_grouped_pha_get_mask(make_grouped_pha):
    """What is the default get_mask value?"""
    pha = make_grouped_pha
    assert pha.get_mask() == pytest.approx([True] * 4 + [False])


def test_quality_pha_mask(make_quality_pha):
    """What is the default mask setting?"""
    pha = make_quality_pha
    assert np.isscalar(pha.mask)
    assert pha.mask is True

    pha.ignore_bad()
    assert np.isscalar(pha.mask)
    assert pha.mask is True


def test_quality_pha_get_mask(make_quality_pha):
    """What is the default get_mask value?"""
    pha = make_quality_pha
    assert pha.get_mask() == pytest.approx([True] * 9)

    pha.ignore_bad()
    # This is a regression test
    assert pha.get_mask() == pytest.approx([True] + [False] * 2 + [True] * 3 + [False] * 3)


@pytest.mark.parametrize("field,expected",
                         [("channel", np.arange(1, 10)),
                          ("counts", [1, 2, 0, 3, 12, 2, 9, 8, 7]),
                          ("grouping", [1, -1, -1, -1, 1, 1, 1, -1, -1]),
                          ("quality", [0, 5, 5, 0, 0, 0, 2, 2, 5]),
                          ("backscal", [0.2, 99, 98, 0.4, 0.5, 0.6, 2, 3, 4]),
                          ("areascal", 0.9)
                          ])
def test_quality_pha_fields(field, expected, make_quality_pha):
    """Regression test: do we get the expected fields?"""

    pha = make_quality_pha

    # fake in a backscal array but scalar areascal
    pha.areascal = 0.9
    pha.backscal = [0.2, 99, 98, 0.4, 0.5, 0.6, 2, 3, 4]

    pha.ignore_bad()

    assert getattr(pha, field) == pytest.approx(expected)


# Should this really return "1:4" as the fifth channel has been
# excluded? At the moment check the current behavior.
#
def test_grouped_pha_get_filter(make_grouped_pha):
    """What is the default get_filter value?"""
    pha = make_grouped_pha
    assert pha.get_filter() == "1:5"


def test_grouped_pha_set_filter(make_grouped_pha):
    """What happens with a simple filter?"""
    pha = make_grouped_pha
    pha.ignore(hi=2)
    assert pha.get_filter() == "4"


# What should get_dep(filter=False) return here? Should it
# include the quality=2 filtered bin (12) or not? At the
# moment it does, so we test against this behavior, but it
# might be something we want to change.
#
@pytest.mark.parametrize("filter,expected",
                         [(False, [1, 2, 0, 3, 12]),
                          (True, [3, 3])])
def test_grouped_pha_get_dep(filter, expected, make_grouped_pha):
    """Quality filtering and grouping is applied: get_dep"""
    pha = make_grouped_pha
    assert pha.get_dep(filter=filter) == pytest.approx(expected)


@pytest.mark.parametrize("filter,expected",
                         [(False, [1, 2, 0, 3, 12]),
                          (True, [3])])
def test_grouped_pha_filter_get_dep(filter, expected, make_grouped_pha):
    """What happens after a simple filter?

    We do this because the behavior of test_grouped_get_filter
    and test_grouped_pha_get_dep has been unclear.
    """
    pha = make_grouped_pha
    pha.ignore(hi=2)
    assert pha.get_dep(filter=filter) == pytest.approx(expected)


def setup_pha_quality_example():
    """A PHA dataset for some tests."""

    chans = np.arange(1, 11, dtype=np.int16)
    pha = DataPHA("x", chans, chans)
    pha.grouping = [1, -1, -1, 1, 1, 1, -1, 1, 1, 1]
    pha.quality = [0, 5, 0, 0, 0, 2, 2, 0, 0, 5]
    pha.exposure = 5.0

    elo = chans * 0.1
    ehi = elo + 0.1
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)
    pha.set_rmf(rmf)

    # The groups are irrelevant to the model evaluation. We do
    # care about the energy bins
    #
    #   channel   1   2   3   4   5   6   7   8   9   10
    #     elo    0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    #     ehi    0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1
    #    bad?         *               ?   ?           *
    #
    # where
    #    * indicates a "very bad" bin (qual=1 or 5)
    #    ? indicates a "slightly bad" bin (qual=2)
    #
    # The groups are complicated for the first group, since it
    # contanis a bad channel.
    #
    #  group        1     2   3     4     5   6    7
    #  channels  1,2*,3   4   5   6?-7?   8   9   10
    #
    # As of 4.17.0 there's no difference in handing the
    # "severity" of the quality setting.
    #
    pha.set_analysis("energy")
    return pha, elo, ehi


def test_pha_get_indep_grouped_quality():
    """A regression test.

    This is to check the behavior with ignore_bad.

    See also test_pha_eval_model_to_fit_grouped_quality

    """

    pha, _, _ = setup_pha_quality_example()

    expected = np.arange(1, 11, dtype=np.int16)

    # Checks: no filtering
    #
    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected)

    # Checks: what does the bad filter do?
    #
    pha.ignore_bad()
    expected1 = expected[[0, 2, 3, 4, 7, 8]]

    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected1)

    # This range is in the "bad quality" region so it should not
    # change the result.
    #
    pha.ignore(0.6, 0.8)

    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected1)

    # Subset the data
    #
    pha.ignore(hi=0.25)  # this is a bad-quality bin
    pha.ignore(lo=0.95)
    expected2 = expected[[2, 3, 4, 7]]

    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected2)

    # Re-allow the low-energy range, but it should
    # still drop the 0.2-0.3 keV bin. It does not.
    #
    pha.notice(hi=0.4)
    expected3 = expected[[0, 1, 2, 3, 4, 7]]

    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected3)

    # Drop all filters. It should not drop the quality filter, but it
    # currently does.
    #
    pha.notice()

    i1, = pha.get_indep(filter=False)
    i2, = pha.get_indep(filter=True)

    assert i1 == pytest.approx(expected)
    assert i2 == pytest.approx(expected)


def test_pha_eval_model_to_fit_grouped_quality():
    """A regression test.

    This is to check the behavior with ignore_bad.

    See also test_pha_get_indep_grouped_quality
    """

    pha, elo, ehi = setup_pha_quality_example()

    mdl = Polynom1D()
    mdl.c0 = 0.1
    mdl.c1 = 1.2

    # Need to scale by the exposure time for the expected values.
    #
    expected = 5 * mdl(elo, ehi)
    assert (expected > 0).all()  # just check

    resp = pha.get_full_response()
    full_mdl = resp(mdl)

    # Checks: no filtering
    #
    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected)

    # Checks: what does the bad filter do?
    #
    pha.ignore_bad()
    expected1 = expected[[0, 2, 3, 4, 7, 8]]

    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected1)

    # This range is in the "bad quality" region so it should not
    # change the result.
    #
    pha.ignore(0.6, 0.8)

    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected1)

    # Subset the data
    #
    pha.ignore(hi=0.25)  # this is a bad-quality bin
    pha.ignore(lo=0.95)
    expected2 = expected[[2, 3, 4, 7]]

    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected2)

    # Re-allow the low-energy range, but it should
    # still drop the 0.2-0.3 keV bin. It does not.
    #
    pha.notice(hi=0.4)
    expected3 = expected[[0, 1, 2, 3, 4, 7]]

    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected3)

    # Drop all filters. It should not drop the quality filter, but it
    # currently does.
    #
    pha.notice()

    m1 = pha.eval_model(full_mdl)
    m2 = pha.eval_model_to_fit(full_mdl)

    assert m1 == pytest.approx(expected)
    assert m2 == pytest.approx(expected)


def test_grouped_pha_set_y_invalid_size(make_grouped_pha):
    """What happens if change a grouped PHA counts/y setting?

    See also test_grouped_pha_set_related_invalid_size which
    is essentially the same but for other fields.
    """
    pha = make_grouped_pha

    # Pick an array that matches the grouped/filtered data size.
    # This is "not actionable" as we can't work out how to change
    # the counts channels, so it should error.
    #
    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and y: 5 vs 2"):
        pha.set_dep([2, 3])


@pytest.mark.parametrize("related", ["staterror", "syserror",
                                     "y", "counts",
                                     "backscal", "areascal",
                                     "grouping", "quality"])
def test_grouped_pha_set_related_invalid_size(related, make_grouped_pha):
    """Can we set the value to a 2-element array?"""
    pha = make_grouped_pha

    # Pick an array that matches the grouped/filtered data size.
    # This is "not actionable" as we can't work out how to change
    # the counts channels, so it should error.
    #
    # Handle y/counts alias here
    related = "y" if related == "counts" else related
    emsg = f"size mismatch between independent axis and {related}: 5 vs 2"

    with pytest.raises(DataErr,
                       match=emsg):
        setattr(pha, related, [2, 3])


@pytest.mark.parametrize("column", ["staterror", "syserror",
                                    "y", "counts",
                                    "backscal", "areascal",
                                    "grouping", "quality"])
def test_pha_check_related_fields_correct_size(column, make_grouped_pha):
    """Can we set the value to a 2-element array?"""

    d = DataPHA('example', None, None)
    setattr(d, column, np.asarray([2, 10, 3]))

    with pytest.raises(DataErr,
                       match="independent axis can not change size: 3 to 4"):
        d.indep = (np.asarray([2, 3, 4, 5]), )


@pytest.mark.parametrize("label", ["filter", "grouping"])
def test_pha_no_group_apply_xxx_invalid_size(label, make_test_pha):
    """Check apply_filter/grouping tests the data length: no quality/group

    Issue #1439 points out that quality handling creates different results.

    """
    pha = make_test_pha

    func = getattr(pha, f"apply_{label}")
    elabel = "filtered data" if label == "filter" else "data"
    msg = f"size mismatch between {elabel} and array: 4 vs 2"
    with pytest.raises(DataErr,
                       match=f"^{msg}$"):
        func([1, 2])


@pytest.mark.parametrize("vals", [[1], [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_no_group_filtered_apply_filter_invalid_size(vals, make_test_pha):
    """Check apply_filter tests the data length: no quality/group, filtered

    This behaves differently to the apply_grouping case
    """

    pha = make_test_pha
    pha.ignore(hi=2)

    # safety check to make sure we've excluded points
    assert not np.all(pha.mask)
    assert np.any(pha.mask)

    with pytest.raises(DataErr,
                       match="^size mismatch between filtered data and array: 2 vs [18]$"):
        pha.apply_filter(vals)


@pytest.mark.parametrize("vals", [[1], [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_no_group_filtered_apply_grouping_invalid_size(vals, make_test_pha):
    """Check apply_grouping tests the data length: no quality/group, filtered

    This behaves differently to the apply_filter case
    """

    pha = make_test_pha
    pha.ignore(hi=2)

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 4 vs [18]$"):
        pha.apply_grouping(vals)


@pytest.mark.parametrize("label", ["filter", "grouping"])
@pytest.mark.parametrize("vals", [[1], (2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_zero_quality_apply_xxx_invalid_size(label, vals, make_test_pha):
    """Check apply_filter/grouping tests the data length: quality set to 0

    We can not use make_grouped_pha and then set the quality array to
    0's as that does not remove the quality setting (see issue #1427)
    so we replicate most of make_grouped_pha with make_test_pha.

    """
    pha = make_test_pha
    pha.grouping = [1, -1, -1, 1]
    pha.quality = [0] * 4
    pha.group()

    func = getattr(pha, f"apply_{label}")
    elabel = "filtered data" if label == "filter" else "data"
    msg = f"size mismatch between {elabel} and array: 4 vs [128]"
    with pytest.raises(DataErr,
                       match=f"^{msg}$"):
        func(vals)


@pytest.mark.parametrize("vals", [(2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_zero_quality_filtered_apply_filter_invalid_size(vals, make_test_pha):
    """Check apply_filter tests the data length: quality set to 0, filtered"""

    pha = make_test_pha
    pha.grouping = [1, -1, -1, 1]
    pha.quality = [0] * 4
    pha.group()

    pha.ignore(hi=2)

    # safety check to make sure we've excluded points
    assert not np.all(pha.mask)
    assert np.any(pha.mask)

    with pytest.raises(DataErr,
                       match="^size mismatch between filtered data and array: 1 vs [28]$"):
        pha.apply_filter(vals)


@pytest.mark.parametrize("vals", [[1], (2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_zero_quality_filtered_apply_grouping_invalid_size(vals, make_test_pha):
    """Check apply_grouping tests the data length: quality set to 0, filtered"""

    pha = make_test_pha
    pha.grouping = [1, -1, -1, 1]
    pha.quality = [0] * 4
    pha.group()

    pha.ignore(hi=2)

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 4 vs [128]$"):
        pha.apply_grouping(vals)


@pytest.mark.parametrize("vals", [(2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_quality_apply_filter_invalid_size(vals, make_grouped_pha):
    """Check apply_filter tests the data length: with quality set"""

    pha = make_grouped_pha

    with pytest.raises(DataErr,
                       match="^size mismatch between filtered data and array: 4 vs [28]$"):
        pha.apply_filter(vals)


@pytest.mark.parametrize("vals", [(2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_quality_filtered_apply_filter_invalid_size(vals, make_grouped_pha):
    """Check apply_filter tests the data length: with quality set, filtered"""

    pha = make_grouped_pha
    pha.ignore(hi=1)

    # safety check to make sure we've excluded points
    assert pha.mask == pytest.approx([False, True])
    assert pha.get_mask() == pytest.approx([False, False, False, True])

    with pytest.raises(DataErr,
                       match="^size mismatch between filtered data and array: 1 vs [28]$"):
        pha.apply_filter(vals)


@pytest.mark.parametrize("vals", [pytest.param([42], marks=pytest.mark.xfail), [10, 20, 35, 42, 55]])
def test_pha_quality_filtered_apply_filter_match_filter(vals, make_grouped_pha):
    """What happens if the array has the correct size?"""

    pha = make_grouped_pha
    pha.ignore(hi=1)

    # XFAIL: when the array matches the filtered data there's a problem
    # matching to the ungrouped data.
    got = pha.apply_filter(vals)
    assert got == pytest.approx([42])


@pytest.mark.parametrize("vals", [[1], (2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_quality_apply_grouping_invalid_size(vals, make_grouped_pha):
    """Check apply_grouping tests the data length: with quality set"""

    pha = make_grouped_pha

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 5 vs [128]$"):
        pha.apply_grouping(vals)


@pytest.mark.parametrize("vals", [[1], (2, 3), [1, 2, 3, 4, 5, 6, 7, 8]])
def test_pha_quality_filtered_apply_grouping_invalid_size(vals, make_grouped_pha):
    """Check apply_grouping tests the data length: with quality set, filtered"""

    pha = make_grouped_pha
    pha.ignore(hi=1)

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 5 vs [128]+$"):
        pha.apply_grouping(vals)


def test_pha_quality_apply_grouping_size_matches_detchans(make_grouped_pha):
    """Regression test: send in DETCHANS values to group.

    Is this an error or not? There's arguments for both, so test what we do.
    """

    pha = make_grouped_pha
    pha.ignore(hi=1)
    got = pha.apply_grouping([10, 12, 2, 4, 84])
    assert got == pytest.approx([24, 4])


def test_pha_quality_apply_grouping_size_matches_quality(make_grouped_pha):
    """Regression test: send in"good" values to group.

    Is this an error or not? There's arguments for both, so test what we do.
    """

    pha = make_grouped_pha
    pha.ignore(hi=1)
    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 5 vs 4$"):
        pha.apply_grouping([10, 12, 2, 4])


def test_pha_apply_filter_check():
    """Check that apply_filter works as expected.

    We go through a number of stages - e.g.

      - no filter or group
      - only group
      - group and filter

    """

    chans = np.arange(1, 21)
    counts = np.ones(20)
    data = DataPHA("ex", chans, counts)

    all_vals = np.arange(1, 21)
    filt_vals = np.arange(5, 17)

    expected = np.arange(1, 21)
    got = data.apply_filter(all_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)

    grouping = np.asarray([1, -1] * 10)
    data.grouping = grouping
    data.quality = [0] * 20

    assert not data.grouped
    data.group()
    assert data.grouped

    expected = np.asarray([3, 7, 11, 15, 19, 23, 27, 31, 35, 39])
    got = data.apply_filter(all_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)

    # This ignores the first two groups, channels 1-2 and 3-4,
    # and the last two groups, channels 17-18 and 19-20.
    # Note that channel 17 is ignored even though not explicitly
    # asked because of the use of ignore.
    #
    data.ignore(hi=4)
    data.ignore(lo=18)

    expected = np.asarray([11, 15, 19, 23, 27, 31])
    got = data.apply_filter(all_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)

    # Now the data has been filtered we can check what happens when
    # the input argument has less channels in.
    #
    got = data.apply_filter(filt_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)

    # Remove the grouping
    #
    data.ungroup()
    assert not data.grouped

    # Note that we still send in vals=arange(5, 17)
    #
    expected = filt_vals.copy()
    got = data.apply_filter(filt_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)

    # We can send in the full array too
    got = data.apply_filter(all_vals, groupfunc=np.sum)
    assert got == pytest.approx(expected)



@pytest.mark.parametrize("backscal", [0, -2.1e-10, -12])
def test_datapha_negative_scale_check(backscal):
    """Check of < 0 condition in DataPHA._check_scale.

    The current test relies on us not restricting the
    backscal value to > 0. Perhaps we should do that and
    then we may be able to remove the code being tested
    in _check_scale.
    """

    # Create a PHA with no ARF or RMF
    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  exposure=10.0, backscal=0.2)

    assert pha.get_backscal() == pytest.approx(0.2)

    pha.backscal = backscal
    assert pha.get_backscal() == pytest.approx(1.0)


def test_datapha_apply_grouping_quality_filter_length_check():
    """Check we get an error for this case."""

    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
    grouping = np.asarray([1, -1, 1, 1])
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  grouping=grouping)

    assert pha.grouped

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and quality_filter: 4 vs 5"):
        pha.quality_filter = np.asarray([1, 1, 1, 1, 1], dtype=bool)


def test_datapha_apply_grouping_quality_filter_scalar():
    """A regression test.

    Check we get an error for this case.
    """

    channels = np.arange(1, 5, dtype=np.int16)
    counts = np.asarray([10, 5, 12, 7], dtype=np.int16)
    grouping = np.asarray([1, -1, 1, 1])
    pha = DataPHA('test-pha', channel=channels, counts=counts,
                  grouping=grouping)

    assert pha.grouped

    # TODO: The error message here is not ideal.
    #
    with pytest.raises(DataErr,
                       match="^Array must be a sequence or None$"):
        pha.quality_filter = True


@requires_fits
@requires_data
def test_xmmrgs_notice(make_data_path):
    """Test that notice and ignore works on XMMRGS dataset, which is
    ordered in increasing wavelength, not energy"""
    from sherpa.astro.io import read_pha, read_rmf
    dat = read_pha(make_data_path('xmmrgs/P0112880201R1S004SRSPEC1003.FTZ'))
    rmf = read_rmf(make_data_path('xmmrgs/P0112880201R1S004RSPMAT1003.FTZ'))
    dat.set_rmf(rmf)
    dat.units = 'wave'
    dat.notice(18.8, 19.2)
    assert len(dat.get_dep(filter=True)) == 41
    assert dat.get_filter(format='%.2f') == '18.80:19.21'

    dat.ignore(10, 19.)
    assert len(dat.get_dep(filter=True)) == 20
    assert dat.get_filter(format='%.2f') == '19.01:19.21'


def test_pickle_image_filter_none(make_test_image):
    """Check we can pickle/unpickle without a region filter.

    This test assumes we have region support, but we do not
    currently have any test builds without it so do not
    bother skipping.

    """

    d = make_test_image
    assert d._region is None

    d2 = pickle.loads(pickle.dumps(d))
    assert d2._region is None


@requires_region
@pytest.mark.parametrize("ignore,region,expected",
                         [(False, 'circle(4255, 3840, 20)', 'Circle(4255,3840,20)'),
                          (True, 'circle(4255, 3840, 20)', 'Field()&!Circle(4255,3840,20)'),
                          (False, 'circle(4255, 3840, 20) - field()', 'Circle(4255,3840,20)&!Field()'),
                          (True, 'circle(4255, 3840, 20) - field()', 'Field()&!Circle(4255,3840,20)|Field()'),
                          ])
def test_pickle_image_filter(ignore, region, expected, make_test_image):
    """Check we can pickle/unpickle with a region filter."""

    from sherpa.astro.utils._region import Region

    d = make_test_image
    d.notice2d(region, ignore=ignore)
    assert isinstance(d._region, Region)
    assert str(d._region) == expected

    d2 = pickle.loads(pickle.dumps(d))
    assert isinstance(d2._region, Region)
    assert str(d2._region) == expected


def test_img_sky_create(make_test_image_sky):
    d = make_test_image_sky
    assert d.sky is not None
    assert d.eqpos is None


def test_img_world_create(make_test_image_world):
    d = make_test_image_world
    assert d.sky is None
    assert d.eqpos is not None


def test_img_sky_show(make_test_image_sky):
    d = make_test_image_sky
    out = str(d).split("\n")
    assert out[0] == "name      = sky-ey"
    assert out[1] == "x0        = Int64[6]"
    assert out[2] == "x1        = Int64[6]"
    assert out[3] == "y         = Float64[6]"
    assert out[4] == "shape     = (2, 3)"
    assert out[5] == "staterror = None"
    assert out[6] == "syserror  = None"
    assert out[7] == "sky       = physical"
    assert out[8] == " crval    = [ 2000.5,-5000.5]"
    assert out[9] == " crpix    = [-2., 3.]"
    assert out[10] == " cdelt    = [2.,4.]"
    assert out[11] == "eqpos     = None"
    assert out[12] == "coord     = logical"
    assert len(out) == 13


def test_img_world_show(make_test_image_world):
    d = make_test_image_world
    out = str(d).split("\n")
    assert out[0] == "name      = world-ey"
    assert out[1] == "x0        = Int64[6]"
    assert out[2] == "x1        = Int64[6]"
    assert out[3] == "y         = Float64[6]"
    assert out[4] == "shape     = (2, 3)"
    assert out[5] == "staterror = None"
    assert out[6] == "syserror  = None"
    assert out[7] == "sky       = None"
    assert out[8] == "eqpos     = world"
    assert out[9] == " crval    = [30.,10.]"
    assert out[10] == " crpix    = [2.,2.]"
    assert out[11] == " cdelt    = [-0.1, 0.1]"
    assert out[12] == " crota    = 0"
    assert out[13] == " epoch    = 2000"
    assert out[14] == " equinox  = 2000"
    assert out[15] == "coord     = logical"
    assert len(out) == 16


@requires_wcs
def test_img_sky_pickle(make_test_image_sky):
    """Very basic test of pickling"""
    d = make_test_image_sky
    d.set_coord("physical")

    d2 = pickle.loads(pickle.dumps(d))
    assert d2.coord == "physical"
    assert d2.eqpos is None

    # We don't have an easy way to check for WCS equivalence
    # so just rely on string representation.
    #
    assert str(d2.sky) == str(d.sky)

    # check the independent axes are converted
    assert (d2.x0 == d.x0).all()
    assert (d2.x1 == d.x1).all()


@requires_wcs
def test_img_world_pickle(make_test_image_world):
    """Very basic test of pickling"""
    d = make_test_image_world
    d.set_coord("wcs")

    d2 = pickle.loads(pickle.dumps(d))
    assert d2.coord == "world"
    assert d2.sky is None

    # We don't have an easy way to check for WCS equivalence
    # so just rely on string representation.
    #
    assert str(d2.sky) == str(d.sky)

    # check the independent axes are converted
    assert (d2.x0 == d.x0).all()
    assert (d2.x1 == d.x1).all()


@requires_wcs  # not for all, but easiest this way
@pytest.mark.parametrize("path", [[],
                                  ["logical"],
                                  ["physical", "logical", "physical", "logical", "physical", "logical"]])
def test_img_sky_logical(path, make_test_image_sky):
    """The logical axes are as expected. Inspired by issue 1380."""
    d = make_test_image_sky
    for coord in path:
        d.set_coord(coord)

    x1, x0 = np.mgrid[1:3, 1:4]
    assert (d.x0 == x0.flatten()).all()
    assert (d.x1 == x1.flatten()).all()


@requires_wcs  # not for all, but easiest this way
@pytest.mark.parametrize("path", [[],
                                  ["logical"],
                                  ["world", "logical", "world", "logical", "world", "logical"]])
def test_img_world_logical(path, make_test_image_world):
    """The logical axes are as expected. Inspired by issue 1380."""
    d = make_test_image_world
    for coord in path:
        d.set_coord(coord)

    x1, x0 = np.mgrid[1:3, 1:4]
    assert (d.x0 == x0.flatten()).all()
    assert (d.x1 == x1.flatten()).all()


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["physical"],
                                  ["logical", "physical", "logical", "physical", "logical"]])
def test_img_sky_physical(path, make_test_image_sky):
    """The physical axes are as expected. Inspired by issue 1380."""
    d = make_test_image_sky
    for coord in path:
        d.set_coord(coord)

    d.set_coord("physical")
    x1, x0 = np.mgrid[1:3, 1:4]
    x0 = (x0 + 2.0) * 2.0 + 2000.5
    x1 = (x1 - 3.0) * 4.0 - 5000.5

    assert (d.x0 == x0.flatten()).all()
    assert (d.x1 == x1.flatten()).all()


def test_img_world_physical(make_test_image_world):
    """The physical axes are not defined."""
    d = make_test_image_world
    with pytest.raises(DataErr,
                       match="data set 'world-ey' does not contain a physical coordinate system"):
        d.set_coord("physical")


def test_img_sky_world(make_test_image_sky):
    """The world axes are not defined."""
    d = make_test_image_sky
    with pytest.raises(DataErr,
                       match="data set 'sky-ey' does not contain a world coordinate system"):
        d.set_coord("world")


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["logical"],
                                  ["world", "logical", "world", "logical", "world", "logical"]])
def test_img_world_world(path, make_test_image_world):
    """The world axes are as expected. Inspired by issue 1380."""
    d = make_test_image_world
    for coord in path:
        d.set_coord(coord)

    d.set_coord("world")
    assert d.x0 == pytest.approx(WORLD_X0)
    assert d.x1 == pytest.approx(WORLD_X1)


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["physical"],
                                  ["logical", "physical", "logical", "physical"]])
def test_img_sky_get_logical(path, make_test_image_sky):
    """Check get_logical works"""
    d = make_test_image_sky
    for coord in path:
        d.set_coord(coord)

    x1, x0 = np.mgrid[1:3, 1:4]
    a, b = d.get_logical()
    assert (a == x0.flatten()).all()
    assert (b == x1.flatten()).all()


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["world"],
                                  ["logical", "world", "logical", "world"]])
def test_img_world_get_logical(path, make_test_image_world):
    """Check get_logical works"""
    d = make_test_image_world
    for coord in path:
        d.set_coord(coord)

    x1, x0 = np.mgrid[1:3, 1:4]
    a, b = d.get_logical()
    assert (a == x0.flatten()).all()
    assert (b == x1.flatten()).all()


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["logical"],
                                  ["physical", "logical", "physical", "logical", "physical", "logical"]])
def test_img_sky_get_physical(path, make_test_image_sky):
    """Check get_physical works"""
    d = make_test_image_sky
    for coord in path:
        d.set_coord(coord)

    x1, x0 = np.mgrid[1:3, 1:4]
    x0 = (x0 + 2.0) * 2.0 + 2000.5
    x1 = (x1 - 3.0) * 4.0 - 5000.5

    a, b = d.get_physical()
    assert (a == x0.flatten()).all()
    assert (b == x1.flatten()).all()


def test_img_world_get_physical(make_test_image_world):
    """Check get_physical errors out"""
    d = make_test_image_world
    with pytest.raises(DataErr,
                       match="data set 'world-ey' does not contain a physical coordinate system"):
        d.get_physical()


def test_img_sky_get_world(make_test_image_sky):
    """Check get_world errors out"""
    d = make_test_image_sky
    with pytest.raises(DataErr,
                       match="data set 'sky-ey' does not contain a world coordinate system"):
        d.get_world()


@requires_wcs
@pytest.mark.parametrize("path", [[],
                                  ["logical"],
                                  ["world", "logical", "world", "logical"]])
def test_img_world_get_world(path, make_test_image_world):
    """Check get_world works"""
    d = make_test_image_world
    for coord in path:
        d.set_coord(coord)

    a, b = d.get_world()
    assert a == pytest.approx(WORLD_X0)
    assert b == pytest.approx(WORLD_X1)


@requires_wcs
@requires_region
def test_img_sky_can_filter(make_test_image_sky):
    """Check we can filter the image using physical coordinates"""
    data = make_test_image_sky
    assert data.coord == "logical"
    assert data.mask
    assert data.get_filter() == ""

    data.set_coord("physical")
    assert data.coord == "physical"
    assert data.mask
    assert data.get_filter() == ""

    data.notice2d("rect(2009,-5006,2011,-5000)", ignore=True)
    assert data.mask == pytest.approx([1, 1, 1, 1, 1, 0])
    assert data.get_filter() == "Field()&!Rectangle(2009,-5006,2011,-5000)"


@requires_wcs
@requires_region
def test_img_sky_can_filter_change_coords(make_test_image_sky, caplog):
    """What happens to a filter after changing coordinates?

    This is a regression test.
    """
    data = make_test_image_sky

    data.set_coord("physical")
    data.notice2d("rect(2009,-5006,2011,-5000)", ignore=True)

    assert len(caplog.records) == 0
    with caplog.at_level(logging.INFO, logger='sherpa'):
        data.set_coord("image")

    assert len(caplog.records) == 1

    assert data.coord == "logical"
    assert data.mask
    assert data.get_filter() == ""

    r = caplog.record_tuples[0]
    assert r[0] == "sherpa.astro.data"
    assert r[1] == logging.WARN
    assert r[2] == "Region filter has been removed from 'sky-ey'"


def test_arf_checks_energy_length():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = np.arange(2, 9)
    dummy = []

    with pytest.raises(ValueError,
                       match="The energy arrays must have the same size, not 4 and 7"):
        DataARF("dummy", elo, ehi, dummy)


def test_rmf_checks_energy_length():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = np.arange(2, 9)
    dummy = []

    with pytest.raises(ValueError,
                       match="The energy arrays must have the same size, not 4 and 7"):
        DataRMF("dummy", 1024, elo, ehi, dummy, dummy, dummy, dummy)


@pytest.mark.parametrize("startchan,na,nb,nc",
                         [(0, 0, 2, 17),  # channel 0 is not defined, what should it do?
                          (1, 0, 3, 16),
                          (2, 0, 4, 15),
                          (3, 1, 4, 14),
                          (4, 2, 4, 13),
                          (9, 7, 4, 8),
                          (12, 10, 4, 5),
                          (16, 14, 4, 1),
                          (17, 15, 4, 0),
                          (18, 16, 3, 0), # channel 20 is not defined, what should it do?
                          (19, 17, 2, 0), # channel 20,21 is not defined, what should it do?
                          (20, 18, 1, 0), # channel 20,21,22 is not defined, what should it do?
                          (21, 19, 0, 0), # channel 21,22,23 is not defined, what should it do?
                          ])
def test_rmf_simple_filter_check(startchan, na, nb, nc):
    """Check we behave as expected for a very simple filter.

    nc is not needed since it's 19 - na - nb, but specify it
    This assumes offset=1.

    It is not at all clear why the filter selects 4 channels once away
    from the edges. Is it a small bug in the code?

    """

    egrid = np.arange(0.1, 2.1, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    # This is a "perfect" response.
    #
    mvals = np.linspace(1, 19, 19)
    assert rmf.apply_rmf(mvals) == pytest.approx(mvals)

    selected = rmf.notice([startchan, startchan + 1, startchan + 2])
    expected = [False] * na + [True] * nb + [False] * nc
    assert selected == pytest.approx(expected)

    # Drop everything but the selected values.
    #
    expected = mvals.copy()
    expected[~selected] = 0
    assert rmf.apply_rmf(mvals[selected]) == pytest.approx(expected)


def test_rmf_complains_about_filter_mismatch():
    """Check we error out if, after filtering, we are given the wrong size."""

    egrid = np.arange(0.1, 2.1, 0.1)
    elo = egrid[:-1]
    ehi = egrid[1:]
    rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

    rmf.notice([9, 10, 11])
    msg = "Mismatched filter between ARF and RMF or PHA and RMF"
    with pytest.raises(TypeError, match=f"^{msg}$"):
        rmf.apply_rmf([1] * 19)


def test_rmf_invalid_offset():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = elo + 1
    dummy = []

    with pytest.raises(ValueError,
                       match="offset must be >=0, not -1"):
        DataRMF("dummy", 1024, elo, ehi, dummy, dummy, dummy, dummy, offset=-1)


def test_rmf_invalid_offset_not_integer():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = elo + 1
    dummy = []

    with pytest.raises(ValueError,
                       match="^offset must be an integer, not 0.5$"):
        DataRMF("dummy", 1024, elo, ehi, dummy, dummy, dummy, dummy, offset=0.5)


@pytest.mark.parametrize("offset", [0, 1, 2, 5])
def test_rmf_offset_check_square(offset, caplog):
    """What happens if offset is set?

    This uses a blurry / made-up RMF but with a 1 to 1 mapping of
    channel and energy grids.

    """

    # The channel range is offset, offset + 1, offset + 2, offset + 3.
    # The response is not "perfect".
    #
    elo = np.asarray([1.1, 1.2, 1.4, 1.5])
    ehi = np.asarray([1.2, 1.4, 1.5, 2.0])
    rmf = DataRMF("dummy", 4, elo, ehi,
                  n_grp=np.asarray([1, 1, 1, 1]),
                  f_chan=np.asarray([offset, offset, offset + 1, offset + 2]),
                  n_chan=np.asarray([1, 2, 2, 2]),
                  matrix=np.asarray([1.0, 0.6, 0.4, 0.5, 0.5, 0.2, 0.8]),
                  e_min=elo, e_max=ehi, offset=offset)

    assert len(caplog.record_tuples) == 0

    # The model, evaluated on the bins, returns
    #    mdl.c0 * bin-width
    # So, for this situation
    #    [0.2, 0.4, 0.2, 1.0]
    #
    mdl = Const1D()
    mdl.c0 = 2
    mvals = mdl(elo, ehi)
    assert mvals == pytest.approx([0.2, 0.4, 0.2, 1.0])

    # The RMF is slightly blury:
    # - [1.0, 0.0, 0.0, 0.0]
    # - [0.6, 0.4, 0.0, 0.0]
    # - [0.0, 0.5, 0.5, 0.0]
    # - [0.0, 0.0, 0.2, 0.8]
    #
    expected = [1.0 * 0.2 + 0.6 * 0.4,
                0.4 * 0.4 + 0.5 * 0.2,
                0.5 * 0.2 + 0.5 * 0.4,
                0.8 * 1.0]

    got = rmf.apply_rmf(mvals)
    assert got == pytest.approx(expected)

    # Now filter the response and try again. Select just the
    # third bin, which should select the last three "energy"
    # bins.
    #
    nchans = [offset + 2]

    selected = rmf.notice(nchans)
    assert selected == pytest.approx([False, True, True, True])

    expected2 = [0.0 * 0.2 + 0.6 * 0.4,
                 0.4 * 0.4 + 0.5 * 0.2,
                 0.5 * 0.2 + 0.5 * 0.4,
                 0.8 * 1.0]

    got2 = rmf.apply_rmf(mvals[selected])
    assert got2 == pytest.approx(expected2)


@pytest.mark.parametrize("offset", [0, 1, 2, 5])
def test_rmf_offset_check_rectangular(offset):
    """What happens if offset is set?

    Here the RMF grid has more bins than channels, and it's more
    constrained than the one_to_one case which might cover different
    code paths in the notice call.

    """

    # 10 channels and 20 energy bins.
    #
    ematrix = np.linspace(0.1, 21 * 0.1, 21)
    ebounds = np.linspace(0.05, 2.2, 11)
    elo = ematrix[:-1]
    ehi = ematrix[1:]
    e_min = ebounds[:-1]
    e_max = ebounds[1:]

    # Create the full matrix.
    #
    full_matrix = np.zeros((20, 10))
    for idx in range(2, 18):
        n = idx // 2
        full_matrix[idx][n:n + 2] = [0.6, 0.4]

    full_matrix[1][1:3] = [0.8, 0.2]
    full_matrix[18][8:10] = [0.2, 0.8]
    full_matrix[19][-1] = 1.0

    # Special case the multi-group rows.
    #
    full_matrix[0] = [0.0, 0.1, 0.1, 0.0, 0.4, 0.4, 0.0, 0.0, 0.0, 0.0]
    full_matrix[2] = [0.0, 0.1, 0.1, 0.0, 0.0, 0.0, 0.4, 0.4, 0.0, 0.0]
    full_matrix[12] = [0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0, 0.4, 0.4]
    full_matrix[18] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.4, 0.4]

    # Cmopress the matrix into the form that DataRMF needs. This step
    # assumes that offset=1.
    #
    n_grp, f_chan, n_chan, matrix = matrix_to_rmf(full_matrix)

    # Correct for the offset.
    #
    f_chan += (offset - 1)

    rmf = DataRMF("dummy", 10, elo, ehi, n_grp=n_grp, f_chan=f_chan,
                  n_chan=n_chan, matrix=matrix, e_min=e_min,
                  e_max=e_max, offset=offset)

    # Have different values for each bin: 0.5, 0.6, ..., 2.4
    #
    mvals = np.linspace(0.5, 2.4, 20)

    # Since we have the full matrix we can calculate it rather than
    # work it out manually.
    #
    expected = mvals @ full_matrix

    got = rmf.apply_rmf(mvals)
    assert got == pytest.approx(expected)

    # Now filter the response and try again:
    #
    # We check the RMF convolution still works after each notice call,
    # even if nothing has been ignored, just to try and check all the
    # corner cases.
    #
    # - drop the first channel
    nchans1 = offset + np.arange(1, 10)

    selected1 = rmf.notice(nchans1)
    assert selected1.all()
    expected1 = expected

    got1 = rmf.apply_rmf(mvals[selected1])
    assert got1 == pytest.approx(expected1)

    # - drop the last channel
    nchans2 = offset + np.arange(0, 9)

    selected2 = rmf.notice(nchans2)
    assert selected2 == pytest.approx([True] * 19 + [False])

    selected2 = rmf.notice(offset + np.arange(0, 9))
    assert selected2 == pytest.approx([True] * 19 + [False])

    expected2 = mvals[selected2] @ full_matrix[selected2, :]
    got2 = rmf.apply_rmf(mvals[selected2])
    assert got2 == pytest.approx(expected2)

    # - drop a range of channels (start and end)
    #
    # This is a regression test as there's no documentation on
    # filter_resp to indicate what it should return in this case.
    # DJB thinks it should have returned
    #
    #  T F T F F F T T T T T T T T F F F F T F
    #
    # rather than
    #
    #  T F T F T T T T T T T T T T F F F F T F
    #
    # [or it could not fik=kter out those ranges within the
    # start/end points].
    #
    nchans3 = offset + np.asarray([4, 5, 6])

    selected3 = rmf.notice(nchans3)
    assert selected3 == pytest.approx([True, False] * 2 + [True] * 10 + [False] * 4 + [True, False])

    # It is not clear what the RMF application does here.
    #
    # What did I naively expect?
    #
    # expected3 = mvals[selected3] @ full_matrix[selected3, :]
    # assert expected3 == pytest.approx([0, 0.12, 1.26, 2.14, 2.91, 3.54, 2.83, 1, 1.6, 1.6])
    #
    # How is this calculated?
    expected3 = [0, 0, 1.14, 2.14, 2.91, 3.54, 2.83, 1, 0, 0]

    got3 = rmf.apply_rmf(mvals[selected3])
    assert got3 == pytest.approx(expected3)

    # Check we get back the original values.
    #
    assert rmf.notice(None) is None

    got_last = rmf.apply_rmf(mvals)
    assert got_last == pytest.approx(expected)


@pytest.mark.parametrize("offset", [0, 1, 2, 5])
def test_rmf_offset_check_basics(offset):
    """Check out one of the tests from
    sherpa/astro/utils/tests/test_astro_utils.py::test_filter_resp_basics

    """

    fm = np.asarray([[0.0, 0.0, 0.0, 0.0, 0.0],
                     [0.0, 1.0, 0.0, 0.0, 0.0],
                     [0.0, 1.2, 1.8, 0.0, 0.0],
                     [0.0, 0.0, 2.0, 0.0, 0.0],
                     [0.0, 0.0, 3.3, 3.7, 0.0],
                     [0.0, 0.0, 0.0, 4.0, 0.0]])

    n_grp, f_chan, n_chan, matrix = matrix_to_rmf(fm, startchan=offset)

    # Correct for the offset.
    #
    delta = offset - 1

    # Make up energy ranges
    #
    ematrix = np.asarray([0.1, 0.2, 0.25, 0.35, 0.4, 0.5, 0.6])
    elo = ematrix[:-1]
    ehi = ematrix[1:]

    ebounds = np.asarray([0.05, 0.15, 0.25, 0.35, 0.45, 0.55])
    emin = ebounds[:-1]
    emax = ebounds[1:]

    rmf = DataRMF("x", 5, n_grp=n_grp, f_chan=f_chan, n_chan=n_chan,
                  matrix=matrix, offset=offset, energ_lo=elo,
                  energ_hi=ehi, e_min=emin, e_max=emax)

    mvals = np.asarray([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])

    # No filter:
    expected = mvals @ fm
    got = rmf.apply_rmf(mvals)
    assert got == pytest.approx(expected)

    # Now apply channel selections; need to correct for the offset.
    #
    # First channel:
    nchans = np.asarray([1]) + delta

    selected = rmf.notice(nchans)
    expected = mvals[selected] @ fm[selected, :]
    # Check that we get "no model" out
    assert expected == pytest.approx(np.zeros(5))
    got = rmf.apply_rmf(mvals[selected])
    assert got == pytest.approx(expected)

    # Second channel:
    nchans = np.asarray([2]) + delta

    selected = rmf.notice(nchans)
    expected = mvals[selected] @ fm[selected, :]
    got = rmf.apply_rmf(mvals[selected])
    assert got == pytest.approx(expected)

    # Add an explicit check here, rather than one that just checks
    # it is internally consistent.
    #
    assert got == pytest.approx([0, 0.56, 0.54, 0, 0])

    # Third channel:
    nchans = np.asarray([3]) + delta

    selected = rmf.notice(nchans)
    expected = mvals[selected] @ fm[selected, :]
    got = rmf.apply_rmf(mvals[selected])
    assert got == pytest.approx(expected)

    # Fourth channel:
    nchans = np.asarray([4]) + delta

    selected = rmf.notice(nchans)
    expected = mvals[selected] @ fm[selected, :]
    got = rmf.apply_rmf(mvals[selected])
    assert got == pytest.approx(expected)

    # Fifth channel:
    nchans = np.asarray([5]) + delta

    selected = rmf.notice(nchans)
    expected = mvals[selected] @ fm[selected, :]
    got = rmf.apply_rmf(mvals[selected])
    assert got == pytest.approx(expected)

    # [2,3] and [1,2,3] should give the same results
    # as channel 1 doesn't match anything
    #
    nchans1 = np.asarray([2,3]) + delta
    nchans2 = np.asarray([1,2,3]) + delta

    selected1 = rmf.notice(nchans1)
    got1 = rmf.apply_rmf(mvals[selected1])
    selected2 = rmf.notice(nchans2)
    got2 = rmf.apply_rmf(mvals[selected2])

    assert selected2 == pytest.approx(selected1)
    assert got2 == pytest.approx(got1)

    # Add an explicit check here, rather than one that just checks
    # it is internally consistent.
    #
    assert got1 == pytest.approx([0, 0.56, 2.99, 1.85, 0])


@pytest.mark.parametrize("subtract", [True, False])
def test_pha_no_bkg(subtract):
    """Just check we error out

    Given the way the code works, it errors out both ways.
    """

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr,
                       match="data set 'dummy' does not have any associated backgrounds"):
        pha.subtracted = subtract


@pytest.mark.parametrize("attr", ["response", "background"])
def test_pha_xxx_ids_invalid_not_an_iterable(attr):
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr,
                       match=f"{attr} ids 'None' does not appear to be an array"):
        setattr(pha, f"{attr}_ids", None)


@pytest.mark.parametrize("attr", ["response", "background"])
def test_pha_xxx_ids_invalid_not_known(attr):
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr,
                       match=re.escape(f"3 is not a valid {attr} id in []")):
        setattr(pha, f"{attr}_ids", [3])


def test_pha_set_analysis_rate_invalid():
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr,
                       match="unknown plot type 'None', choose 'rate' or 'counts'"):
        pha.set_analysis("channel", type=None)


def test_pha_get_ylabel_yfac0():
    """This does not depend on the backend"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    assert pha.plot_fac == 0
    assert pha.get_ylabel() == 'Counts/channel'



def test_pha_get_ylabel_yfac1(all_plot_backends):
    """Basic check

    The label depends on the backend, so we just want the dummy
    backend used here.
    """

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    pha.plot_fac = 1

    ylabel = pha.get_ylabel()
    assert ylabel.startswith('Counts/channel X Channel')
    assert "1" in ylabel


@requires_fits
@requires_data
def test_1209_rsp(make_data_path):
    """Do we pick up the header keywords from a RSP matrix.

    This is related to issue #1209
    """

    # We could set up channels and counts, but let's not.
    #
    d = DataPHA("dummy", None, None)
    assert d.header["TELESCOP"] == "none"
    assert d.header["INSTRUME"] == "none"
    assert d.header["FILTER"] == "none"

    infile = make_data_path("xmmrgs/P0112880201R1S004RSPMAT1003.FTZ")
    rsp = io.read_rmf(infile)
    d.set_rmf(rsp)

    assert d.header["TELESCOP"] == "XMM"
    assert d.header["INSTRUME"] == "RGS1"
    assert d.header["FILTER"] == "NONE"


@requires_fits
@requires_data
@pytest.mark.parametrize("mode,fexpr",
                         [(["arf", "rmf"], ""),
                          (["arf"], ""),
                          (["rmf"], "Medium")])
def test_1209_response(mode, fexpr, make_data_path):
    """Do we pick up the header keywords from ARF and/or RMF

    This is related to issue #1209

    We use a non-Chandra dataset for the responses just to
    ensure we understand other missions. Note that SWIFT and ROSAT
    are tested in test_astro_data_xxx_unit.py so we want something
    other than those two.
    """

    # We could set up channels and counts, but let's not.
    #
    d = DataPHA("dummy", None, None)
    assert d.header["TELESCOP"] == "none"
    assert d.header["INSTRUME"] == "none"
    assert d.header["FILTER"] == "none"

    # We hide the warnings about ENERG_LO being 0 in the input files
    # as we are not testing this here.
    #
    with warnings.catch_warnings(record=True):

        if "arf" in mode:
            infile = make_data_path("MNLup_2138_0670580101_EMOS1_S001_spec.arf")
            arf = io.read_arf(infile)
            d.set_arf(arf)

        if "rmf" in mode:
            infile = make_data_path("MNLup_2138_0670580101_EMOS1_S001_spec.rmf")
            rmf = io.read_rmf(infile)
            d.set_rmf(rmf)

    assert d.header["TELESCOP"] == "XMM"
    assert d.header["INSTRUME"] == "EMOS1"

    # The FILTER setting:
    #   is ""       in the ARF
    #      "Medium" in the RMF
    # so the output depends on the selected response.
    #
    # We could work it out, but specify it as an input to the test.
    # It turns out to be a good test that we see different behavior
    # depending on the loaded data!
    #
    assert d.header["FILTER"] == fexpr


@requires_fits
@requires_data
def test_1209_background(make_data_path):
    """Do we pick up the header keywords from the background?

    This is related to issue #1209

    We use a non-Chandra dataset.
    """

    # We need to set up the channels array to match the background.
    #
    d = DataPHA("dummy", np.arange(0, 800, dtype=np.int16), None)
    assert d.header["TELESCOP"] == "none"
    assert d.header["INSTRUME"] == "none"
    assert d.header["FILTER"] == "none"

    infile = make_data_path("MNLup_2138_0670580101_EMOS1_S001_specbg.fits")
    bkg = io.read_pha(infile)
    d.set_background(bkg)

    assert d.header["TELESCOP"] == "XMM"
    assert d.header["INSTRUME"] == "EMOS1"
    assert d.header["FILTER"] == "Medium"


@pytest.fixture
def make_dataimgint():
    """Create a simple IMG Int data set."""

    # a 1 by 2 grid.
    #
    x1, x0 = np.mgrid[10:12, -5:-4]
    shape = x0.shape
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.asarray([10, 5])

    return DataIMGInt("ival", x0, x1, x0 + 1, x1 + 1,
                      y, shape=shape)


def test_dataimgint_create(make_dataimgint):
    """Check we can create a basic integrated image data set.

    See issue #1379
    """

    x0 = np.asarray([-5, -5])
    x1 = np.asarray([10, 11])

    img = make_dataimgint

    assert (img.dep == [10, 5]).all()

    assert len(img.indep) == 4
    assert (img.indep[0] == x0).all()
    assert (img.indep[1] == x1).all()
    assert (img.indep[2] == (x0 + 1)).all()
    assert (img.indep[3] == (x1 + 1)).all()

    assert img.header == {}


def test_dataimgint_show(make_dataimgint):
    """Check we can show a basic integrated image data set.

    See issue #1379
    """

    img = make_dataimgint

    # This fails because there's problems getting x0 and x0lo
    # attributes.
    #
    out = str(img).split("\n")

    # Do we expect the x0/x1 output or x0lo/../x1hi
    # output? For the moment just test what we do return.
    #
    assert out[0] == "name      = ival"
    assert out[1] == "x0        = Float64[2]"
    assert out[2] == "x1        = Float64[2]"
    assert out[3] == "y         = Int64[2]"
    assert out[4] == "shape     = (2, 1)"
    assert out[5] == "staterror = None"
    assert out[6] == "syserror  = None"
    assert out[7] == "sky       = None"
    assert out[8] == "eqpos     = None"
    assert out[9] == "coord     = logical"
    assert len(out) == 10


def test_dataimgint_x0lo(make_dataimgint):
    assert make_dataimgint.x0lo == pytest.approx([-5, -5])


def test_dataimgint_x1lo(make_dataimgint):
    assert make_dataimgint.x1lo == pytest.approx([10, 11])


def test_dataimgint_x0hi(make_dataimgint):
    assert make_dataimgint.x0hi == pytest.approx([-4, -4])


def test_dataimgint_x1hi(make_dataimgint):
    assert make_dataimgint.x1hi == pytest.approx([11, 12])


def test_dataimgint_get_x0(make_dataimgint):
    x0 = np.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_dataimgint.get_x0() == x).all()


def test_dataimgint_x0(make_dataimgint):
    x0 = np.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_dataimgint.x0 == x).all()


def test_dataimgint_get_x1(make_dataimgint):
    x1 = np.asarray([10, 11])
    x = (x1 + x1 + 1) / 2

    assert (make_dataimgint.get_x1() == x).all()


def test_dataimgint_x1(make_dataimgint):
    x1 = np.asarray([10, 11])
    x = (x1 + x1 + 1) / 2

    assert (make_dataimgint.x1 == x).all()


def test_dataimgint_get_y(make_dataimgint):
    assert (make_dataimgint.get_y() == [10, 5]).all()


def test_dataimgint_y(make_dataimgint):
    assert (make_dataimgint.y == [10, 5]).all()


def test_dataimgint_get_dep(make_dataimgint):
    assert (make_dataimgint.get_dep() == [10, 5]).all()


def test_dataimgint_get_x0label(make_dataimgint):
    assert make_dataimgint.get_x0label() == "x0"


def test_dataimgint_get_x1label(make_dataimgint):
    assert make_dataimgint.get_x1label() == "x1"


def test_dataimgint_get_ylabel(make_dataimgint):
    assert make_dataimgint.get_ylabel() == "y"


def test_dataimgint_get_axes(make_dataimgint):
    """This copies the Data2DInt case but is different"""
    axes = make_dataimgint.get_axes()
    assert len(axes) == 4

    # What are these values? They are not the input values
    # to DataIMGInt.
    #
    assert (axes[0] == [-0.5]).all()
    assert (axes[1] == [-0.5, 0.5]).all()
    assert (axes[2] == [0.5]).all()
    assert (axes[3] == [0.5, 1.5]).all()


@pytest.mark.xfail
def test_dataimgint_notice(make_dataimgint):
    """basic notice call

    It is not entirely clear whether we expect the
    notice call to work here when notice2d is present.
    """
    img = make_dataimgint

    # The mask attribute can be True, False, or a ndarray. Fortunately
    # using an ndarray as a truthy value throws a ValueError.
    #
    assert img.mask

    # Data is defined on x0=-5, x1=10,11
    # so this excludes the second point.
    #
    img.notice(x1lo=10, x1hi=11)
    assert (img.mask == np.asarray([True, False])).all()


@pytest.mark.xfail
def test_dataimgint_ignore(make_dataimgint):
    """basic ignore call"""
    img = make_dataimgint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert (img.mask == np.asarray([False, True])).all()


def test_dataimgint_ignore_get_filter(make_dataimgint):
    """What exactly does get_filter return here?

    The current behavior does not look sensible.
    """
    img = make_dataimgint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert img.get_filter() == ''


def test_dataimgint_ignore_get_filter_expr(make_dataimgint):
    """What exactly does get_filter_expr return here?

    The current behavior does not look sensible.
    """
    img = make_dataimgint

    assert img.mask
    img.notice(x1lo=10, x1hi=11, ignore=True)
    assert img.get_filter_expr() == ''


# given how notice test above fails, how is this working?
def test_dataimgint_notice_get_x0(make_dataimgint):
    """basic notice call + get_x0"""
    img = make_dataimgint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_x0() == np.asarray([-4.5, -4.5])).all()
    assert (img.get_x0(True) == np.asarray([-4.5])).all()


@pytest.mark.xfail
def test_dataimgint_notice_get_x1(make_dataimgint):
    """basic notice call + get_x1"""
    img = make_dataimgint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_x1() == np.asarray([10.5, 11.5])).all()
    assert (img.get_x1(True) == np.asarray([10.5])).all()


@pytest.mark.xfail
def test_dataimgint_notice_get_y(make_dataimgint):
    """basic notice call + get_y"""
    img = make_dataimgint
    img.notice(x1lo=10, x1hi=11)
    assert (img.get_y() == np.asarray([10, 5])).all()
    assert (img.get_y(True) == np.asarray([10])).all()


@requires_region
def test_dataimgint_notice2d(make_dataimgint):
    """basic notice2d call.

    Given that we only have two items the testing is not
    going to be extensive.
    """
    img = make_dataimgint

    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.mask == np.asarray([True, False])).all()


@requires_region
def test_dataimgint_ignore2d(make_dataimgint):
    """basic ignore2d call.

    Given that we only have two items the testing is not
    going to be extensive.
    """
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)", ignore=True)
    assert (img.mask == np.asarray([False, True])).all()


@requires_region
def test_dataimgint_notice2d_get_filter(make_dataimgint):
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert img.get_filter() == 'Rectangle(-100,10,100,11)'


@requires_region
def test_dataimgint_notice2d_get_filter_expr(make_dataimgint):
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert img.get_filter_expr() == 'Rectangle(-100,10,100,11)'


@requires_region
def test_dataimgint_notice2d_get_x0(make_dataimgint):
    """basic notice2d call + get_x0"""
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.get_x0() == np.asarray([-4.5, -4.5])).all()
    assert (img.get_x0(True) == np.asarray([-4.5])).all()


@requires_region
def test_dataimgint_notice2d_get_x1(make_dataimgint):
    """basic notice2d call + get_x1"""
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.get_x1() == np.asarray([10.5, 11.5])).all()
    assert (img.get_x1(True) == np.asarray([10.5])).all()


@requires_region
def test_dataimgint_notice2d_get_y(make_dataimgint):
    """basic notice2d call + get_y"""
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.get_y() == np.asarray([10, 5])).all()
    assert (img.get_y(True) == np.asarray([10])).all()


def test_dataimgint_get_dims(make_dataimgint):
    assert make_dataimgint.get_dims() == (1, 2)


def test_dataimgint_get_img(make_dataimgint):
    img = make_dataimgint
    ival = img.get_img()
    assert ival.shape == (2, 1)
    assert (ival == np.asarray([[10], [5]])).all()


def test_dataimgint_get_img_model_no_filter(make_dataimgint):
    """Check we can evaluate a model

    The Data2DInt case also adds a filter to check that the routine
    ignores this filter, but as we currently don't understand the
    filtering we skip this step.

    """
    img = make_dataimgint

    # This model evaluates
    #   mdl.c + mdl.cx1 * x0 + mdl.cy1 * x1
    #
    # which becomes, because we use the middle of the bin
    #
    #   10 + 1 * (-4.5) + 10 * (10.5, 11.5)
    #   = (110.5, 120.5)
    #
    mdl = Polynom2D()
    mdl.c = 10
    mdl.cy1 = 10
    mdl.cx1 = 1

    ivals = img.get_img(mdl)
    assert len(ivals) == 2
    assert ivals[0].shape == (2, 1)
    assert ivals[1].shape == (2, 1)
    assert (ivals[0] == np.asarray([[10], [5]])).all()
    assert (ivals[1] == np.asarray([[110.5], [120.5]])).all()


def test_dataimgint_get_max_pos(make_dataimgint):
    assert make_dataimgint.get_max_pos() == (-4.5, 10.5)


def test_dataimgint_get_bounding_mask(make_dataimgint):
    assert make_dataimgint.get_bounding_mask() == (True, None)


@pytest.mark.parametrize("method",
                         ["get_error",
                          "get_imgerr",
                          "get_staterror",
                          "get_syserror",
                          "get_yerr"
                         ])
def test_dataimgint_method_is_none(method, make_dataimgint):
    """Check those methods that return None"""
    func = getattr(make_dataimgint, method)
    assert func() is None


@pytest.mark.parametrize("attribute",
                         ["eqpos",
                          "sky",
                          "staterror",
                          "syserror"
                         ])
def test_dataimgint_attribute_is_none(attribute, make_dataimgint):
    """Check those attributes that return None"""
    attr = getattr(make_dataimgint, attribute)
    assert attr is None


def test_dataimgint_no_sky(make_dataimgint):
    """Basic check (rely on base class to check all the combinations)."""

    with pytest.raises(DataErr,
                       match="data set 'ival' does not contain a physical coordinate system"):
        make_dataimgint.get_physical()


@requires_wcs
def test_dataimgint_sky(make_dataimgint):
    """We can convert coordinates.

    We assume the base class tests are good here, so this is a
    minimal check.
    """

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])

    # The "logical" coordinates are
    #  lo = [-5, 10], [-5, 11]
    #  hi = [-4, 11], [-4, 12]
    #
    # so these get converted to
    #
    #   new = (orig - crpix) * cdelt + crval
    #
    # which is
    #  lo = [87.5, 125.5], [87.5, 127.5]
    #  hi = [89.5, 127.5], [89.5, 129.5]
    #
    x0 = np.asarray([87.5, 87.5])
    x1 = np.asarray([125.5, 127.5])

    sky = img.get_physical()

    assert len(sky) == 4
    assert sky[0] == pytest.approx(x0)
    assert sky[1] == pytest.approx(x1)
    assert sky[2] == pytest.approx(x0 + 2)
    assert sky[3] == pytest.approx(x1 + 2)


def test_dataimgint_sky_coords_unchanged(make_dataimgint):
    """Just because sky is set we don't change axis data."""

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])

    x1 = np.asarray([10, 11])
    x = (x1 + x1 + 1) / 2
    assert img.get_x1() == pytest.approx(x)


@requires_wcs
def test_dataimgint_set_sky(make_dataimgint):
    """We can change to the SKY coordinate system"""

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])

    assert img.coord == "logical"
    img.set_coord("physical")
    assert img.coord == "physical"


@requires_wcs
def test_dataimgint_set_sky_x0hi(make_dataimgint):
    """x0hi is changed

    We don't check all attributes.
    """

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])

    img.set_coord("physical")
    x0 = np.asarray([87.5, 87.5])
    x = x0 + 2
    assert img.x0hi == pytest.approx(x)


@requires_wcs
def test_dataimgint_set_sky_get_x1(make_dataimgint):
    """get_x1 is changed

    We don't check all accessors.
    """

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])

    img.set_coord("physical")
    x1 = np.asarray([125.5, 127.5])
    x = (x1 + x1 + 2) / 2
    assert img.get_x1() == pytest.approx(x)


@requires_wcs
def test_dataimgint_sky_coords_reset(make_dataimgint):
    """We can get back to the logical units

    We only check one of the values.
    """

    img = make_dataimgint
    img.sky= WCS("sky", "LINEAR",
                 crval=[100.5, 110.5],
                 crpix=[1.5, 2.5],
                 cdelt=[2, 2])
    img.set_coord("physical")
    img.set_coord("logical")

    x1 = np.asarray([10, 11])
    x = (x1 + x1 + 1) / 2
    assert img.get_x1() == pytest.approx(x)


@pytest.mark.parametrize("dclass", [Data2D, DataIMG])
def test_1379_evaluation_unintegrated(dclass):
    """Check that delta2d does not evaluate (i.e. only 0's).

    This is based on the code that lead to showing #1379.
    """

    x0, x1, y, shape = dataspace2d([10,15])

    data = dclass("temp", x0, x1, y, shape=shape)

    # It is important that xpos/ypos is not set to either an integer
    # value or to half-pixel (as this is used for the bin edges in the
    # integrated case).
    #
    mdl = Delta2D("mdl")
    mdl.xpos = 4.3
    mdl.ypos = 8.9
    mdl.ampl = 100
    assert mdl.integrate  # ensure it's integrates

    out = data.eval_model(mdl)
    assert len(out) == len(y)
    assert set(out) == {0.0}


@pytest.mark.parametrize("dclass", [Data2DInt, DataIMGInt])
def test_1379_evaluation_integrated(dclass):
    """Check that delta2d does get evaluate at some point.

    This is based on the code that lead to showing #1379.
    """

    x0, x1, y, shape = dataspace2d([10,15])

    data = dclass("temp", x0 - 0.5, x1 - 0.5,
                  x0 + 0.5, x1 + 0.5, y,
                  shape=shape)

    mdl = Delta2D("mdl")
    mdl.xpos = 4.3
    mdl.ypos = 8.9
    mdl.ampl = 100
    assert mdl.integrate  # ensure it's integrates

    out = data.eval_model(mdl)
    assert len(out) == len(y)
    assert set(out) == {0.0, 100.0}
    assert out.sum() == 100.0
    assert out.argmax() == 83

    # An internal check that we are actually selecting the correct
    # pixel.
    #
    assert x0[83] == 4.0
    assert x1[83] == 9.0


@pytest.mark.parametrize("dclass", [Data2DInt, DataIMGInt])
def test_1379_evaluation_model_not_integrated(dclass):
    """If the delta2D model is not integrated all bets are off.

    This is based on the code that lead to showing #1379.
    """

    x0, x1, y, shape = dataspace2d([10,15])

    data = dclass("temp", x0 - 0.5, x1 - 0.5,
                  x0 + 0.5, x1 + 0.5, y,
                  shape=shape)

    mdl = Delta2D("mdl")
    mdl.xpos = 4.3
    mdl.ypos = 8.9
    mdl.ampl = 100
    mdl.integrate = False

    # As the integrate flag is False it behaves like the Data2D/DataIMG
    # case (since the xpos/ypos value is chosen not to fall on a
    # bin edge).
    #
    out = data.eval_model(mdl)
    assert len(out) == len(y)
    assert set(out) == {0.0}


@requires_fits
@requires_data
@requires_wcs
@pytest.mark.parametrize("coord", ["logical", "image", "physical", "world", "wcs"])
def test_1380_data(coord, make_data_path):
    """The contour data should ideally remain the same.

    See also sherpa/astro/ui/tests/test_astro_ui_plot.py::test_1380_plot

    This is the origin of the problem.
    """

    infile = make_data_path("image2.fits")
    img = io.read_image(infile)

    assert isinstance(img, DataIMG)
    assert img.coord == "logical"

    (x0_1, x1_1, y_1, xl_1, yl_1) = img.to_contour()

    # We do not check the output of this call. It is important
    # to call set_coord rather than change the coord attribute.
    #
    img.set_coord(coord)
    img.to_contour()

    # We do check that we get back the same data as we
    # originally did.
    #
    img.set_coord("logical")
    (x0_3, x1_3, y_3, xl_3, yl_3) = img.to_contour()

    assert xl_3 == xl_1
    assert yl_3 == yl_1
    assert (y_3 == y_1).all()
    assert (x0_3 == x0_1).all()
    assert (x1_3 == x1_1).all()


@requires_fits
@requires_data
@requires_wcs
def test_1380_pickle(make_data_path):
    """Can we pickle and restore an image?

    The fix for 1380 added new data that is pickled, so just
    check it works. Technically this should work but the
    state handling has had to be tweaked to allow old state
    files to be read in, so this just checks that new data
    is not affected by this. We don't have any "old" state
    files lying around that we can use here.

    There are a number of existing image pickle tests but
    they don't check the coordinate settings used here.
    """

    infile = make_data_path("image2.fits")
    img = io.read_image(infile)
    img.set_coord("physical")

    x0_1, x1_1 = img.indep

    img2 = pickle.loads(pickle.dumps(img))

    assert img2.coord == "physical"

    x0_2, x1_2 = img2.indep

    # this test should not need pytest.approx
    assert (x0_2 == x0_1).all()
    assert (x1_2 == x1_1).all()

    img2.set_coord("logical")
    assert img.coord == "physical"
    assert img2.coord == "logical"

    img.set_coord("logical")

    x0_3, x1_3 = img.indep
    x0_4, x1_4 = img2.indep

    assert (x0_4 == x0_3).all()
    assert (x1_4 == x1_3).all()

    assert (x0_3 != x0_1).all()

    # This is an internal check and may get changed if the
    # implementation changes.
    #
    assert img._orig_indep_axis[0] == "logical"
    assert img2._orig_indep_axis[0] == "logical"


def test_image_apply_filter_invalid_size(make_test_image):
    """Does an image error out if the filter is sent an invalid size?

    Test related to issue #1439 which is an issue with the DataPHA class.
    """

    data = make_test_image

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 600 vs 2$"):
        data.apply_filter([1, 2])


def test_image_filtered_apply_filter_invalid_size(make_test_image):
    """Does an image error out if the filter is sent an invalid size after a filter?"""

    data = make_test_image

    # "Fake" a filter (this is a perfectly-valid way to set up a filter,
    # at least at the time this code was written).
    #
    data.mask = np.ones(data.y.size, dtype=bool)
    data.mask[0] = False

    with pytest.raises(DataErr,
                       match="^size mismatch between data and array: 600 vs 2$"):
        data.apply_filter([1, 2])


def test_pha_subtract_bkg_no_staterror():
    """Check what happens with no staterror function for the background.

    The idea is that the data has a statistical error column but the
    background does not. In this case we can't calculate an error. The
    code currently returns None rather than raising an error.

    """

    chans = np.arange(1, 5)
    counts = np.asarray([10, 9, 3, 7])
    errs = np.asarray([3, 2, 1, 2])
    data = DataPHA("ex", chans, counts, errs)

    bcounts = np.asarray([2, 1, 2, 4])
    bkg = DataPHA("bkg", chans, bcounts)

    data.set_background(bkg)
    data.subtract()

    assert bkg.staterror is None
    assert data.staterror is not None

    # The data.get_staterror call is the code being tested here, but
    # the other asserts are being made to check that nothing has
    # changed.
    #
    assert bkg.get_staterror() is None
    assert data.get_staterror() is None


def test_pha_subtract_bkg_filter_false():
    """Check what happens with background and filter=False

    Looks like this has not been tested, so add an explicit check.

    """

    # Note that the background has a different set of groups to the
    # data, but it is over-ridden when computing the
    # background-subtracted values.
    #
    chans = np.arange(1, 5)
    counts = np.asarray([10, 9, 3, 7])
    grps = np.asarray([1, 1, -1, 1])
    data = DataPHA("ex", chans, counts, grouping=grps)

    bcounts = np.asarray([2, 0, 2, 4])
    bgrps = np.asarray([1, -1, -1, 1])
    bkg = DataPHA("bkg", chans, bcounts, grouping=bgrps)

    data.set_background(bkg)
    data.subtract()

    bgot = bkg.get_staterror(filter=False, staterrfunc=np.sqrt)
    assert bgot == pytest.approx([2, 2])

    expected = np.sqrt(np.asarray([10, 12, 7]) + np.asarray([2, 2, 4]))
    got = data.get_staterror(filter=False, staterrfunc=np.sqrt)
    assert got == pytest.approx(expected)


def test_pha_subtract_bkg_filter_cih2datavar():
    """Check what happens with background and chi2datavar

    Looks like this has not been tested, so add an explicit check.

    Follows test_pha_subtract_bkg_filter_false but uses a different
    error function.

    """

    chans = np.arange(1, 5)
    counts = np.asarray([10, 9, 3, 7])
    data = DataPHA("ex", chans, counts)

    bcounts = np.asarray([2, 0, 2, 4])
    bkg = DataPHA("bkg", chans, bcounts)

    data.set_background(bkg)
    data.subtract()

    bgot = bkg.get_staterror(staterrfunc=calc_chi2datavar_errors)
    assert bgot == pytest.approx([np.sqrt(2), 0, np.sqrt(2), np.sqrt(4)])

    expected = np.sqrt(np.asarray([10, 9, 3, 7]) + np.asarray([2, 0, 2, 4]))
    got = data.get_staterror(staterrfunc=calc_chi2datavar_errors)
    assert got == pytest.approx(expected)


@requires_group
@pytest.mark.parametrize("opt,arg",
                         [("bins", 2), ("width", 2),
                          ("counts", 3), ("adapt", 3),
                          pytest.param("snr", 3, marks=pytest.mark.xfail),
                          ("adapt_snr", 3)])
def test_1643_group_options(opt, arg):
    """Check the behavior noted in #1646 and if it's been fixed
    yet or not (it comes from the group module which is not
    part of Sherpa).
    """

    pha = DataPHA("grp", np.arange(1, 12),
                  [1, 1, 1, 1, 8, 2, 6, 4, 1, 1, 1])

    pha.ignore(hi=4)
    pha.ignore(lo=9)

    meth = getattr(pha, f"group_{opt}")

    assert pha.get_filter("5:8")
    meth(arg, tabStops=~pha.mask)
    assert pha.get_filter("5:8")

    # excluded channels are 0
    assert pha.grouping[0:4] == pytest.approx([0] * 4)
    assert pha.quality[0:4] == pytest.approx([0] * 4)

    assert pha.grouping[8:] == pytest.approx([0] * 3)
    assert pha.quality[8:] == pytest.approx([0] * 3)


def test_pha_group_background(caplog):
    """Do we group the background?"""

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = [1, 1, -1]

    src.set_background(bkg)

    assert not src.grouped
    assert not bkg.grouped

    src.group()

    assert src.grouped
    assert bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_group_background_not_set(caplog):
    """Do we group the background, but grouping not set?"""

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = None

    src.set_background(bkg)

    assert not src.grouped
    assert bkg.grouped  # TODO: why is this set?

    src.group()

    assert src.grouped
    assert bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_ungroup_background(caplog):
    """Do we ungroup the background?"""

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = [1, 1, -1]
    src.grouped = True
    bkg.grouped = True

    src.set_background(bkg)

    assert src.grouped
    assert bkg.grouped

    src.ungroup()

    assert not src.grouped
    assert not bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_ungroup_background_not_set(caplog):
    """Do we ungroup the background, but grouping not set?"""

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = None
    src.grouped = True

    with pytest.raises(DataErr,
                       match="data set 'bkg' does not specify grouping flags"):
        bkg.grouped = True

    src.set_background(bkg)

    assert src.grouped
    assert bkg.grouped

    src.ungroup()

    assert not src.grouped
    assert not bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_ungroup_background_after(caplog):
    """Set the background before setting the grouping

    The set_background method does some extra work, so let's
    see what happens if we don't do that.
    """

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.set_background(bkg)

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = [1, 1, -1]
    src.grouped = True
    bkg.grouped = True

    assert src.grouped
    assert bkg.grouped

    src.ungroup()

    assert not src.grouped
    assert not bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_ungroup_background_after_not_set(caplog):
    """Set the background before setting the grouping

    The set_background method does some extra work, so let's
    see what happens if we don't do that.
    """

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.set_background(bkg)

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = None
    src.grouped = True

    with pytest.raises(DataErr,
                       match="data set 'bkg' does not specify grouping flags"):
        bkg.grouped = True

    assert src.grouped
    assert not bkg.grouped

    src.ungroup()

    assert not src.grouped
    assert not bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_ungroup_background_remove(caplog):
    """Can we remove the grouping after grouping?

    This is to try and check a corner case.
    """

    src = DataPHA("src", [1, 2, 3], [1, 2, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [0, 1, 1])

    src.set_background(bkg)

    src.grouping = [1, 1, 1]
    src.quality = [0, 0, 0]
    bkg.grouping = [1, 1, -1]
    bkg.quality = [0, 0, 0]
    src.grouped = True
    bkg.grouped = True

    # TODO: should this remove the grouping flag?
    bkg.grouping = None

    src.ungroup()

    assert not src.grouped
    assert not bkg.grouped

    assert len(caplog.record_tuples) == 0


def test_pha_check_background_ids_basic():
    """Check background_ids can be used.

    This is a basic run through to check the behavior.
    """

    pha = DataPHA("src", [1, 2, 3], [1, 1, 1])
    b1 = DataPHA("b1", [1, 2, 3], [1, 1, 1])
    b2 = DataPHA("b2", [1, 2, 3], [1, 1, 1])

    assert len(pha.background_ids) == 0

    pha.set_background(b1)
    assert pha.background_ids == pytest.approx([1])

    pha.set_background(b2, id="up")
    assert pha.background_ids == pytest.approx([1, "up"])

    pha.delete_background(1)
    assert pha.background_ids == pytest.approx(["up"])

    # Check we can delete a background by setting background_ids
    #
    pha.background_ids = []
    assert pha.background_ids == []

    # Does the order matter?
    #
    pha.set_background(b2, id="up")
    pha.set_background(b1)
    assert pha.background_ids == pytest.approx(["up", 1])

    # Remove one element.
    #
    pha.background_ids = [1]
    assert pha.background_ids == pytest.approx([1])

    # We can do the following, which is technically a no-op but may
    # change some internal state. This also tests using an
    # iterable-that-is-not-a-list.
    #
    pha.background_ids = {1}
    assert pha.background_ids == pytest.approx([1])

    # We can change to a currently-unused background as long as we've
    # used it before.
    #
    pha.background_ids = ["up"]
    assert pha.background_ids == pytest.approx(["up"])

    pha.background_ids = [1, "up"]
    assert pha.background_ids == pytest.approx([1, "up"])

    # We can not set an unknown background.
    #
    with pytest.raises(DataErr,
                       match=r"^foo is not a valid background id in \['up', 1\]$"):
        pha.background_ids = ["foo"]

    # And it hasn't changed.
    #
    assert pha.background_ids == pytest.approx([1, "up"])


@pytest.mark.parametrize("bkg_id", [None, 1, "foo"])
def test_pha_delete_unknown_background(bkg_id):
    """What happens if we delete an unknown background?"""

    pha = DataPHA("src", [1, 2, 3], [1, 1, 1])
    b1 = DataPHA("b1", [1, 2, 3], [1, 1, 1])
    b2 = DataPHA("b2", [1, 2, 3], [1, 1, 1])

    pha.set_background(b1, id="up")
    pha.set_background(b1, id="down")

    # This is treated as a no-op
    pha.delete_background(bkg_id)

    assert pha.background_ids == pytest.approx(["up", "down"])


def test_pha_check_response_ids_basic():
    """Check response_ids can be used.

    This is a basic run through to check the behavior.
    """

    pha = DataPHA("src", [1, 2, 3], [1, 1, 1])

    elo = np.arange(2, 5)
    ehi = elo + 1
    rmf1 = create_delta_rmf(elo, ehi, name="rmf1")
    rmf2 = create_delta_rmf(elo, ehi, name="rmf2")

    pha.set_response(rmf=rmf1)
    assert pha.response_ids == pytest.approx([1])

    pha.set_response(rmf=rmf2, id="up")
    assert pha.response_ids == pytest.approx([1, "up"])

    pha.delete_response(1)
    assert pha.response_ids == pytest.approx(["up"])

    # Check we can delete a response by setting response_ids
    #
    pha.response_ids = []
    assert pha.response_ids == []

    # Does the order matter?
    #
    pha.set_response(rmf=rmf2, id="up")
    pha.set_response(rmf=rmf1)
    assert pha.response_ids == pytest.approx(["up", 1])

    # Remove one element.
    #
    pha.response_ids = [1]
    assert pha.response_ids == pytest.approx([1])

    # We can do the following, which is technically a no-op but may
    # change some internal state. This also tests using an
    # iterable-that-is-not-a-list.
    #
    pha.response_ids = {1}
    assert pha.response_ids == pytest.approx([1])

    # We can change to a currently-unused response as long as we've
    # used it before.
    #
    pha.response_ids = ["up"]
    assert pha.response_ids == pytest.approx(["up"])

    pha.response_ids = [1, "up"]
    assert pha.response_ids == pytest.approx([1, "up"])

    # We can not set an unknown response.
    #
    with pytest.raises(DataErr,
                       match=r"^foo is not a valid response id in \['up', 1\]$"):
        pha.response_ids = ["foo"]

    # And it hasn't changed.
    #
    assert pha.response_ids == pytest.approx([1, "up"])


def test_pha_delete_missing_background_is_a_noop():
    """Check this call returns.

    """

    pha = DataPHA("src", [1, 2, 3], [1, 1, 1])
    bkg = DataPHA("bkg", [1, 2, 3], [1, 1, 1])
    pha.set_background(bkg)
    pha.subtract()

    assert pha.subtracted
    assert pha.background_ids == pytest.approx([1])
    assert not bkg.subtracted
    assert bkg.background_ids == []

    # deleting a non-existent background is a no-op:
    #  - dataset with no backgrounds
    #  - dataset with a different-background to the requested
    #
    bkg.delete_background()
    pha.delete_background(2)

    # Minimal check that "nothing has happened".
    #
    assert pha.subtracted
    assert pha.background_ids == pytest.approx([1])
    assert not bkg.subtracted
    assert bkg.background_ids == []


def test_img_checks_coord_nonsense():
    """What happens when coord is set to a nonsense value?"""

    with pytest.raises(DataErr,
                       match="^unknown coordinates: 'nonsense'"):
        DataIMG("ex", [1, 2, 1, 2], [1, 1, 2, 2], [1, 2, 3, 4],
                coord="nonsense")


@pytest.mark.parametrize("coord", ["physical", "world", "wcs"])
def test_img_checks_coord_no_transform(coord):
    """What happens when coord is set to xxx but no transform?"""

    with pytest.raises(DataErr,
                       match="^data set 'ex' does not contain a .* coordinate system$"):
        DataIMG("ex", [1, 2, 1, 2], [1, 1, 2, 2], [1, 2, 3, 4],
                coord=coord)


@requires_group
@pytest.mark.parametrize("asarray", [True, False])
def test_group_xxx_tabtops_not_ndarray(asarray):
    """What happens if tabStops is not a ndarray?"""

    pha = DataPHA("test", [1, 2, 3, 4, 5], [2, 3, 4, 5, 6])
    tabstops = [1, 1, 0, 0, 1]
    if asarray:
        tabstops = np.asarray(tabstops)

    # This should only group channels 3 and 4.
    pha.group_width(2, tabStops=tabstops)

    assert pha.get_y() == pytest.approx([2, 3, 4.5, 6])
    assert pha.mask is True
    assert pha.get_mask() == pytest.approx([True] * 5)


@requires_group
@pytest.mark.parametrize("asarray", [True, False])
@pytest.mark.parametrize("nelem", [4, 6])
def test_group_xxx_tabtops_wrong_size(asarray, nelem):
    """What happens if tabStops is not a ndarray?"""

    pha = DataPHA("test", [1, 2, 3, 4, 5], [2, 3, 4, 5, 6])
    tabstops = [0] * nelem
    if asarray:
        tabstops = np.asarray(tabstops)

    emsg = r"^grpBinWidth\(\) The number of tab stops and number of channels specified in the argument list have different sizes$"
    with pytest.raises(ValueError, match=emsg):
        pha.group_width(2, tabStops=tabstops)


@requires_group
def test_group_xxx_tabstops_already_grouped():
    """Check what happens if tabStops is sent ~pha.mask when already grouped."""

    pha = DataPHA("grp", [1, 2, 3, 4, 5, 6], [12, 2, 9, 2, 4, 5])
    pha.mask = [1, 0, 1, 1, 1, 0]
    assert pha.get_y(filter=True) == pytest.approx([12, 9, 2, 4])

    pha.grouping = [1, 1, 1, -1, 1, 1]
    pha.grouped = True
    assert pha.get_y(filter=True) == pytest.approx([12, 5.5, 4])

    tstops = ~pha.mask
    assert tstops == pytest.approx([False, True, False, False, True])

    # Apply the mask as the tabStops (after inversion) where
    # len(tstops) < nchannel but does match the number of groups.
    #
    pha.group_counts(8, tabStops=tstops)

    assert pha.grouping == pytest.approx([1, 0, 1, 1, -1, 0])
    assert pha.quality == pytest.approx([0, 0, 0, 2, 2, 0])
    assert pha.get_y(filter=True) == pytest.approx([12, 9, 3])


@requires_wcs
@requires_region
def test_dataimg_axis_ordering():
    """What are the x0/x1 axes meant to be?

    See also test_dataimg_axis_ordering_1880
    """

    # nx=3, ny=2
    #
    x1, x0 = np.mgrid[1:3, 1:4]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(6) * 10 + 10

    sky = WCS("physical", "LINEAR", [100, 200], [1, 1], [10, 10])
    eqpos = WCS("world", "WCS", [30, 50], [100, 200], [-0.1, 0.1])
    orig = DataIMG("faked", x0, x1, y, shape=(2, 3), sky=sky,
                   eqpos=eqpos)

    orig.set_coord("physical")

    orig.notice2d("circle(110, 210, 6)", True)

    assert orig.get_filter() == "Field()&!Circle(110,210,6)"
    assert orig.get_dep(filter=True) == pytest.approx([10, 20, 30, 40, 60])
    a0, a1 = orig.get_indep()
    assert a0 == pytest.approx([100, 110, 120] * 2)
    assert a1 == pytest.approx([200, 200, 200, 210, 210, 210])


@requires_wcs
@requires_region
def test_dataimg_axis_ordering_1880():
    """What are the x0/x1 axes meant to be? See issues #1789 #1880

    See also test_dataimg_axis_ordering. This is a regression test to
    catch if we ever decide to update the DataIMG code.

    """

    # nx=2, ny=3
    #
    x1, x0 = np.mgrid[1:4, 1:3]
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = np.arange(6) * 10 + 10

    sky = WCS("physical", "LINEAR", [100, 200], [1, 1], [10, 10])
    eqpos = WCS("world", "WCS", [30, 50], [100, 200], [-0.1, 0.1])
    orig = DataIMG("faked", x0, x1, y, shape=(2, 3), sky=sky,
                   eqpos=eqpos)

    orig.set_coord("physical")

    # This should remove the pixel with value 50 but it actually
    # removes 40. This is an issue with how DataIMG requires the x0/x1
    # arrays (I think) but for this test we just test the existing
    # behavior. See issues #1880 and #1789.
    #
    orig.notice2d("circle(110, 210, 6)", True)

    assert orig.get_filter() == "Field()&!Circle(110,210,6)"
    assert orig.get_dep(filter=True) == pytest.approx([10, 20, 30, 50, 60])
    a0, a1 = orig.get_indep()
    assert a0 == pytest.approx([100, 110] * 3)
    assert a1 == pytest.approx([200, 200, 210, 210, 220, 220])


def test_rmf_get_indep():
    """Check this routine."""

    ebins = np.asarray([3.0, 5., 8.0, 12.0])
    rlo = ebins[:-1]
    rhi = ebins[1:]
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)

    xlo, xhi = rmf.get_indep()
    assert xlo == pytest.approx(rlo)
    assert xhi == pytest.approx(rhi)


def test_rmf_get_dep_simple():
    """Check this routine."""

    ebins = np.asarray([3.0, 5., 8.0, 12.0])
    rlo = ebins[:-1]
    rhi = ebins[1:]
    rmf = create_delta_rmf(rlo, rhi, e_min=rlo, e_max=rhi)

    # This is an "ideal" response.
    #
    y = rmf.get_dep()
    assert y == pytest.approx([1, 1, 1])


def test_rmf_get_dep_complex():
    """Check this routine."""

    ebins = np.asarray([1.1, 1.2, 1.5, 2.0, 3.0, 6.0])
    rlo = ebins[:-1]
    rhi = ebins[1:]
    rmf = DataRMF("x", 5, rlo, rhi,
                  [1, 1, 1, 1, 1],  # n_grp
                  [2, 3, 4, 4, 3],  # f_chan
                  [2, 1, 1, 1, 2],  # n_chan
                  [0.4, 0.6, 1.0, 1.0, 1.0, 0.2, 0.8],  # matrix
                  e_min=rlo, e_max=rhi)

    # This RMF only fills in channels 2 to 4.
    #
    y = rmf.get_dep()
    assert y == pytest.approx([0, 0.4, 1.8, 2.8, 0])
