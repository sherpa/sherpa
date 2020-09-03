#
#  Copyright (C) 2020
#        Smithsonian Astrophysical Observatory
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

import numpy as np

import pytest

from sherpa.astro.data import DataPHA
from sherpa.utils.err import DataErr


def test_can_not_group_ungrouped():
    """Does setting the grouping setting fail with no data?"""

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert not pha.grouped
    with pytest.raises(DataErr) as exc:
        pha.grouped = True

    assert str(exc.value) == "data set 'name' does not specify grouping flags"


def test_get_mask_is_none():

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.mask is True
    assert pha.get_mask() is None


def test_get_filter_expr_channel():
    """Check get_filter_expr is called"""

    pha = DataPHA('name', np.asarray([1, 2, 3]), [1, 1, 1])
    assert pha.get_filter_expr() == '1-3 Channel'

    pha.ignore(None, 1)
    assert pha.get_filter_expr() == '2-3 Channel'


def test_get_filter_is_empty():

    # Need to send in numpy arrays otherwise the code fails as it
    # assumes a numpy array. This should be addressed upstream.
    #
    # pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    pha = DataPHA('name', np.asarray([1, 2, 3]), [1, 1, 1])
    assert pha.get_filter() == '1:3'
    pha.ignore()
    assert pha.get_filter() == 'No noticed bins'


@pytest.mark.xfail
def test_need_numpy_channels():
    """We need to convert the  input to NumPy arrrays.

    This should be done in the constructor to DataPHA, but
    for now error out if not set.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.get_filter() == '1:3'
    # This errors out with a TypeError on 'elo = self.channel - 0.5'
    pha.ignore()
    assert pha.get_filter() == 'No noticed bins'


@pytest.mark.parametrize("chan", [0, -1, 4])
def test_error_on_invalid_channel_ungrouped(chan):
    """Does channel access fail when outside the bounds?

    For ungrouped data it currently does not, but just
    acts as an identity function.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha._to_channel(chan) == chan
    assert pha._from_channel(chan) == chan


@pytest.mark.parametrize("chan,exp1,exp2",
                         [(0, 0, 0),
                          (-1, -1, -1)])
def test_error_on_invalid_channel_grouped(chan, exp1, exp2):
    """Does channel access fail when outside the bounds?

    It is not clear what _from_channel is doing here, so
    just check the responses.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, -1, 1])
    assert pha.grouped
    assert pha._to_channel(chan) == exp1
    assert pha._from_channel(chan) == exp2


@pytest.mark.parametrize("chan", [-2, 4])
def test_error_on_invalid_channel_grouped2(chan):
    """Does channel access fail when outside the bounds?
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, -1, 1])
    assert pha.grouped

    assert pha._from_channel(chan) == chan


def test_288_a():
    """The issue from #288 which was working"""

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    assert pha.mask
    pha.ignore(3, 4)

    # I use approx because it gives a nice answer, even though
    # I want eqiuality not approximation in this test. Fortunately
    # with bools the use of approx is okay (it can tell the
    # difference between 0 and 1, aka False and True).
    #
    assert pha.mask == pytest.approx([True, False, True])


def test_288_b():
    """The issue from #288 which fails"""

    channels = np.arange(1, 6)
    counts = np.asarray([5, 5, 10, 10, 2])
    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    pha = DataPHA('x', channels, counts, grouping=grouping)

    assert pha.mask
    pha.ignore(3.1, 4)

    assert pha.mask == pytest.approx([True, True, True])


@pytest.mark.xfail
def test_grouping_non_numpy():
    """Historically the group* calls will fail oddly if y is not numpy

    TypeError: grpNumCounts() Could not parse input arguments, please check input for correct type(s)
    """

    x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    y = [0, 0, 0, 2, 1, 1, 0, 0, 0, 0]

    pha = DataPHA('416', x, y)
    pha.group_counts(3)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
    assert pha.quality == pytest.approx(quality)


def test_416_a():
    """The first test case from issue #416"""

    # if y is not a numpy array then group_counts errors out
    # with a strange error. Another reason why DataPHA needs
    # to validate input
    #
    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 0, 0])

    pha = DataPHA('416', x, y)
    pha.notice(3.5, 6.5)

    mask = [False, False, False, True, True, True, False, False, False, False]
    assert pha.mask == pytest.approx(mask)

    pha.group_counts(3)

    # We have a simplified mask
    mask = [False, False]
    assert pha.mask == pytest.approx(mask)

    # the "full" mask can be retrieved with get_mask
    mask = [False] * 10
    assert pha.get_mask() == pytest.approx(mask)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
    assert pha.quality == pytest.approx(quality)

    dep = pha.get_dep(filter=True)
    assert dep == pytest.approx([])


def test_416_b(caplog):
    """The second test case from issue #416

    This is to make sure this hasn't changed.
    """

    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 0, 0])

    pha = DataPHA('416', x, y)
    pha.notice(3.5, 6.5)
    pha.group_counts(3)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        pha.ignore_bad()

    # It's not obvious why this has switched to a boolean
    assert pha.mask

    # Mask is also interesting (currently just reporting
    # this behavior)
    mask = [True] * 5 + [False] * 5
    assert pha.get_mask() == pytest.approx(mask)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
    assert pha.quality == pytest.approx(quality)

    dep = pha.get_dep(filter=True)
    assert dep == pytest.approx([3])

    # check captured log
    #
    emsg = 'filtering grouped data with quality flags, previous filters deleted'
    assert caplog.record_tuples == [
        ('sherpa.astro.data', logging.WARNING, emsg)
        ]


def test_416_c():
    """The third test case from issue #416
    """

    x = np.asarray([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
    y = np.asarray([0, 0, 0, 2, 1, 1, 0, 0, 0, 0])

    pha = DataPHA('416', x, y)
    pha.notice(3.5, 6.5)

    # this should be ~pha.mask
    tabstops = [True] * 3 + [False] * 3 + [True] * 4
    assert ~pha.mask == pytest.approx(tabstops)

    pha.group_counts(3, tabStops=~pha.mask)
    pha.ignore_bad()

    grouping = [0] * 3 + [1, -1, 1] + [0] * 4
    assert pha.grouping == pytest.approx(grouping)

    # the second grouped bin has a quality of 2 as
    # it only contains 1 count
    quality = np.zeros(10, dtype=np.int)
    quality[5] = 2
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
