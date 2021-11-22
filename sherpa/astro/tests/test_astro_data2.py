#
#  Copyright (C) 2020, 2021
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

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataIMG, DataPHA, DataRMF
from sherpa.astro.instrument import create_delta_rmf
from sherpa.astro.utils._region import Region
from sherpa.utils.err import DataErr
from sherpa.utils.testing import requires_data, requires_fits, requires_group


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


def test_need_numpy_channels():
    """We didn't used to convert channels to a NumPy array which broke
    this logic - the ignore line would error out due to an operation on
    self.channel
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha.get_filter() == '1:3'

    pha.ignore()
    assert pha.get_filter() == 'No noticed bins'


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


@pytest.mark.parametrize("chan", [0, -1, 4])
def test_error_on_invalid_channel_ungrouped(chan):
    """Does channel access fail when outside the bounds?

    For ungrouped data it currently does not, but just
    acts as an identity function.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1])
    assert pha._from_channel(chan) == chan


@pytest.mark.parametrize("chan,exp1,exp2",
                         [(0, 1, 3),
                          (-1, 1, 1)])
def test_error_on_invalid_channel_grouped(chan, exp1, exp2):
    """Does channel access fail when outside the bounds?

    It is not clear what _from_channel is doing here, so
    just check the responses.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, -1, 1])
    assert pha.grouped
    assert pha._from_channel(chan) == exp2


@pytest.mark.parametrize("chan", [-2, 4])
def test_error_on_invalid_channel_grouped2(chan):
    """Does channel access fail when outside the bounds?

    This one does error out in _from_channel.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, -1, 1])
    assert pha.grouped
    with pytest.raises(DataErr) as exc:
        pha._from_channel(chan)

    # The error message is wrong
    # assert str(exc.value) == 'invalid group number: {}'.format(chan)
    assert str(exc.value) == 'invalid group number: {}'.format(chan - 1)


def test_pha_get_xerr_all_bad_channel_no_group():
    """get_xerr handles all bad values [channel]

    It's not obvious what it is meant to be doing here.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  quality=[2, 2, 2])

    assert pha.get_xerr() == pytest.approx([1, 1, 1])

    pha.ignore_bad()
    assert pha.get_filter() == ''
    assert pha.get_xerr() == pytest.approx([1, 1, 1])


def test_pha_get_xerr_all_bad_channel_group():
    """get_xerr handles all bad values [channel]

    The behavior with grouping is different, presumably because
    we assume we have grouping when we have a quality array.
    """

    pha = DataPHA('name', [1, 2, 3], [1, 1, 1],
                  grouping=[1, 1, 1],
                  quality=[2, 2, 2])

    assert pha.get_xerr() == pytest.approx([1, 1, 1])

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

    assert pha.get_xerr() == pytest.approx([2.0, 3.0, 4.0])

    pha.ignore_bad()
    assert pha.get_filter() == ''
    assert pha.get_xerr() == pytest.approx([2.0, 3.0, 4.0])


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

    assert pha.get_xerr() == pytest.approx([2.0, 3.0, 4.0])

    assert pha.grouped
    pha.ignore_bad()

    # Should this error out or not?
    assert pha.get_filter() == ''
    # with pytest.raises(DataErr) as de:
    #     pha.get_filter()

    # assert str(de.value) == 'mask excludes all data'

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
    with pytest.raises(DataErr) as exc:
        func(lo, hi)

    assert str(exc.value) == f"unknown {lbl} argument: 'must be an integer channel value'"


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

    assert pha.mask
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

    assert pha.mask
    with pytest.raises(DataErr) as de:
        pha.ignore(3.1, 4)

    assert str(de.value) == "unknown lo argument: 'must be an integer channel value'"


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

    assert pha.mask
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

    mask = [False, False, False, True, True, True, False, False, False, False]
    assert pha.mask == pytest.approx(mask)

    pha.group_counts(3)

    # We have a simplified mask
    mask = [True, True]
    assert pha.mask == pytest.approx(mask)

    # the "full" mask can be retrieved with get_mask
    mask = [True] * 10
    assert pha.get_mask() == pytest.approx(mask)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
    assert pha.quality == pytest.approx(quality)

    dep = pha.get_dep(filter=True)
    assert dep == pytest.approx([3, 1])


@requires_group
def test_416_b(caplog):
    """The second test case from issue #416

    This is to make sure this hasn't changed.

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

    pha.group_counts(3, tabStops=~pha.mask)
    pha.ignore_bad()

    grouping = [0] * 3 + [1, -1, 1] + [0] * 4
    assert pha.grouping == pytest.approx(grouping)

    # the second grouped bin has a quality of 2 as
    # it only contains 1 count
    quality = np.zeros(10, dtype=int)
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
def make_test_pha():
    """A simple PHA"""

    chans = np.asarray([1, 2, 3, 4], dtype=np.int16)
    counts = np.asarray([1, 2, 0, 3], dtype=np.int16)
    return DataPHA('p', chans, counts)


def test_img_set_coord_invalid(make_test_image):
    """An invalid coord setting"""
    d = make_test_image
    assert d.coord == 'logical'

    with pytest.raises(DataErr) as exc:
        d.set_coord('bob')

    emsg = "unknown coordinates: 'bob'\n"
    emsg += "Valid options: logical, image, physical, world, wcs"
    assert str(exc.value) == emsg

    assert d.coord == 'logical'


@pytest.mark.parametrize('coord,expected',
                         [('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_set_coord_notset(coord, expected, make_test_image):
    """A valid coord setting but we don't have the data"""

    d = make_test_image
    with pytest.raises(DataErr) as exc:
        d.set_coord(coord)

    emsg = "data set 'd' does not contain a {} coordinate system".format(expected)
    assert str(exc.value) == emsg

    assert d.coord == 'logical'


@pytest.mark.parametrize('coord,expected',
                         [('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_get_coord_notset(coord, expected, make_test_image):
    """Check get_physical/world fail when there's no WCS"""
    d = make_test_image

    meth = getattr(d, 'get_{}'.format(coord))
    with pytest.raises(DataErr) as exc:
        meth()

    emsg = "data set 'd' does not contain a {} coordinate system".format(expected)
    assert str(exc.value) == emsg


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
@pytest.mark.parametrize('coord,expected',
                         [('logical', 'x0'),
                          ('image', 'x0'),
                          ('physical', 'x0 (pixels)'),
                          ('world', 'RA (deg)'),
                          ('wcs', 'RA (deg)')])
def test_img_file_get_xlabel(coord, expected, read_test_image):
    """get_x0label"""
    d = read_test_image
    d.set_coord(coord)
    assert d.get_x0label() == expected


@requires_fits
@requires_data
@pytest.mark.parametrize('coord,expected',
                         [('logical', 'x1'),
                          ('image', 'x1'),
                          ('physical', 'x1 (pixels)'),
                          ('world', 'DEC (deg)'),
                          ('wcs', 'DEC (deg)')])
def test_img_file_get_ylabel(coord, expected, read_test_image):
    """get_x1label"""
    d = read_test_image
    d.set_coord(coord)
    assert d.get_x1label() == expected


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


def test_img_get_bounding_mask_filtered(make_test_image):
    """get_bounding_mask with data partially filtered"""
    d = make_test_image
    d.notice2d('ellipse(4260,3840,3,2,0)')
    ans = d.get_bounding_mask()
    print(np.where(ans[0]))

    mask = np.zeros(5 * 7, dtype=bool)
    for i in [3,  8,  9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24,
              25, 26, 31]:
        mask[i] = True

    assert len(ans) == 2
    assert ans[0] == pytest.approx(mask)
    assert ans[1] == (5, 7)


def test_img_get_filter(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape)
    assert d.get_filter() == shape.capitalize()


def test_img_get_filter_exclude(make_test_image):
    """Simple get_filter check on an image."""
    d = make_test_image
    assert d.get_filter() == ''

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape, ignore=True)

    expected = 'Field()&!' + shape.capitalize()
    assert d.get_filter() == expected


def test_img_get_filter_none(make_test_image):
    """Simple get_filter check on an image: no data"""
    d = make_test_image

    shape = 'ellipse(4260,3840,3,2,0)'
    d.notice2d(shape)
    d.notice2d(ignore=True)

    # It's not clear what the filter should be here
    assert d.get_filter() == ''


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
    """Check removing the shapes works as expected."""

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


def test_img_get_filter_included_excluded(make_test_image):
    """Simple get_filter check on an image.

    Just to match test_img_get_filter_excluded_excluded.
    """
    d = make_test_image
    check_ignore_ignore(d)


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

    exc = "field()-{}-{}".format(shape1, shape2)
    d.notice2d(exc)

    maskb = d.mask.copy()
    assert maskb == pytest.approx(maska)

    # just check we have some True and False values
    assert maska.min() == 0
    assert maska.max() == 1


@pytest.mark.parametrize("requested,expected",
                         [("bin", "channel"), ("Bin", "channel"),
                          ("channel", "channel"), ("ChannelS", "channel"),
                          ("chan", "channel"),
                          ("energy", "energy"), ("ENERGY", "energy"),
                          ("Energies", "energy"),
                          ("WAVE", "wavelength"), ("wavelength", "wavelength"),
                          ("Wavelengths", "wavelength"),
                          ("chan This Is Wrong", "channel"),  # should this be an error?
                          ("WAVEY GRAVY", "wavelength")  # shouls this be an error?
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
    with pytest.raises(DataErr) as de:
        pha.units = invalid

    assert str(de.value) == f"unknown quantity: '{invalid}'"


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


@pytest.mark.xfail
def test_pha_grouping_changed_filter_1160(make_test_pha):
    """What happens when the grouping is changed?

    See also test_pha_grouping_changed_filter_1160
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


@requires_fits
@requires_data
def test_xmmrgs_notice(make_data_path):
    '''Test that notice and ignore works on XMMRGS dataset, which is
    ordered in increasing wavelength, not energy'''
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


@pytest.mark.parametrize("ignore,region,expected",
                         [(False, 'circle(4255, 3840, 20)', 'Circle(4255,3840,20)'),
                          (True, 'circle(4255, 3840, 20)', 'Field()&!Circle(4255,3840,20)'),
                          (False, 'circle(4255, 3840, 20) - field()', 'Circle(4255,3840,20)&!Field()'),
                          (True, 'circle(4255, 3840, 20) - field()', 'Field()&!Circle(4255,3840,20)|Field()'),
                          ])
def test_pickle_image_filter(ignore, region, expected, make_test_image):
    """Check we can pickle/unpickle with a region filter.

    This test assumes we have region support, but we do not
    currently have any test builds without it so do not
    bother skipping.

    """

    d = make_test_image
    d.notice2d(region, ignore=ignore)
    assert isinstance(d._region, Region)
    assert str(d._region) == expected

    d2 = pickle.loads(pickle.dumps(d))
    assert isinstance(d2._region, Region)
    assert str(d2._region) == expected


def test_arf_checks_energy_length():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = np.arange(2, 9)
    dummy = []

    with pytest.raises(ValueError) as ve:
        DataARF("dummy", elo, ehi, dummy)

    assert str(ve.value) == "The energy arrays must have the same size, not 4 and 7"


def test_rmf_checks_energy_length():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = np.arange(2, 9)
    dummy = []

    with pytest.raises(ValueError) as ve:
        DataRMF("dummy", 1024, elo, ehi, dummy, dummy, dummy, dummy)

    assert str(ve.value) == "The energy arrays must have the same size, not 4 and 7"


def test_rmf_invalid_offset():
    """Just check we error out"""

    elo = np.arange(1, 5)
    ehi = elo + 1
    dummy = []

    with pytest.raises(ValueError) as ve:
        DataRMF("dummy", 1024, elo, ehi, dummy, dummy, dummy, dummy, offset=-1)

    assert str(ve.value) == "offset must be >=0, not -1"


@pytest.mark.parametrize("subtract", [True, False])
def test_pha_no_bkg(subtract):
    """Just check we error out

    Given the way the code works, it errors out both ways.
    """

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        pha.subtracted = subtract

    assert str(de.value) == "data set 'dummy' does not have any associated backgrounds"


@pytest.mark.parametrize("attr", ["response", "background"])
def test_pha_xxx_ids_invalid_not_an_iterable(attr):
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        setattr(pha, f"{attr}_ids", None)

    assert str(de.value) == f"{attr} ids 'None' does not appear to be an array"


@pytest.mark.parametrize("attr", ["response", "background"])
def test_pha_xxx_ids_invalid_not_known(attr):
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        setattr(pha, f"{attr}_ids", [3])

    # The error message could be better (use list to remove the dict_keys)
    # but it is not a high priority.
    #
    assert str(de.value) == f"3 is not a valid {attr} id in dict_keys([])"


def test_pha_set_analysis_rate_invalid():
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        pha.set_analysis("channel", type=None)

    assert str(de.value) == "unknown plot type 'none', choose 'rate' or 'counts'"


def test_pha_ignore_bad_no_quality():
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        pha.ignore_bad()

    assert str(de.value) == "data set 'dummy' does not specify quality flags"
