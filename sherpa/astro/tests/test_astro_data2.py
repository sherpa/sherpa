#
#  Copyright (C) 2020, 2021, 2022
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
import warnings

import numpy as np

import pytest

from sherpa.astro.data import DataARF, DataIMG, DataIMGInt, DataPHA, DataRMF
from sherpa.astro.instrument import create_delta_rmf
from sherpa.astro import io
from sherpa.astro.io.wcs import WCS
from sherpa.astro.utils._region import Region
from sherpa.data import Data2D, Data2DInt
from sherpa.models import Delta2D, Polynom2D
from sherpa.plot import backend, dummy_backend
from sherpa.utils import dataspace2d
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
    # assert str(exc.value) == f"invalid group number: {chan}"
    assert str(exc.value) == f"invalid group number: {chan - 1}"


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

    emsg = f"data set 'd' does not contain a {expected} coordinate system"
    assert str(exc.value) == emsg

    assert d.coord == 'logical'


@pytest.mark.parametrize('coord,expected',
                         [('physical', 'physical'),
                          ('world', 'world'),
                          ('wcs', 'world')])
def test_img_get_coord_notset(coord, expected, make_test_image):
    """Check get_physical/world fail when there's no WCS"""
    d = make_test_image

    meth = getattr(d, f"get_{coord}")
    with pytest.raises(DataErr) as exc:
        meth()

    emsg = f"data set 'd' does not contain a {expected} coordinate system"
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

    exc = f"field()-{shape1}-{shape2}"
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
    with pytest.raises(DataErr) as de:
        pha.units = invalid

    assert str(de.value) == f"unknown quantity: '{invalid}'"


@pytest.mark.parametrize("invalid", ["RATE", "COUNTS", "rates", "count", "count-rate"])
def test_pha_analysis_type_invalid(invalid, make_test_pha):
    pha = make_test_pha
    with pytest.raises(DataErr) as err:
        pha.set_analysis("channel", type=invalid)

    assert str(err.value) == f"unknown plot type '{invalid}', choose 'rate' or 'counts'"


def test_pha_analysis_plot_fac_valid(make_test_pha):
    """Historically we've allowed 2.0 as an argument, so check it still works"""
    pha = make_test_pha
    assert pha.plot_fac == 0
    pha.plot_fac = 2.0
    assert pha.plot_fac == 2


@pytest.mark.parametrize("invalid", ["1", 2.01, 0.5, complex(1)])
def test_pha_analysis_plot_fac_invalid(invalid, make_test_pha):
    pha = make_test_pha
    with pytest.raises(DataErr) as err:
        pha.plot_fac = invalid

    assert str(err.value) == f"unknown plot_fac setting: '{invalid}'"


@pytest.mark.parametrize("invalid", ["1", 2.01, 0.5, complex(1)])
def test_pha_analysis_factor_invalid(invalid, make_test_pha):
    pha = make_test_pha
    with pytest.raises(DataErr) as err:
        pha.set_analysis("channel", factor=invalid)

    assert str(err.value) == f"unknown factor setting: '{invalid}'"


def test_pha_get_spectresp_no_response(make_test_pha):
    pha = make_test_pha
    assert pha.get_specresp() is None


def test_pha_ignore_bad_no_quality(make_test_pha):
    pha = make_test_pha
    assert pha.quality is None
    with pytest.raises(DataErr) as err:
        pha.ignore_bad()

    assert str(err.value) == "data set 'p' does not specify quality flags"


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

    # This currently raises a DataErr due to
    # 'size mismatch between mask and data array'
    pha.grouping = [1, 1, -1, 1]
    d4 = pha.get_dep(filter=True)
    assert d4 == pytest.approx([2, 3])


@pytest.mark.xfail
def test_pha_remove_grouping(make_test_pha):
    """Check we can remove the grouping array."""

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
    # This thinks that pha.grouped is still set
    assert not pha.grouped
    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx(no_data)


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


def test_pha_quality_bad_filter(make_test_pha):
    """What is the filter expression when ignore bad + filter"""

    pha = make_test_pha
    assert pha.get_filter() == "1:4"

    pha.ignore(hi=1)
    assert pha.get_filter() == "2:4"

    d1 = pha.get_dep(filter=True)
    assert d1 == pytest.approx([2, 0, 3])

    pha.quality = [0, 0, 0, 2]
    pha.ignore_bad()

    d2 = pha.get_dep(filter=True)
    assert d2 == pytest.approx([2, 0])
    assert pha.get_filter() == "2:3"


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


@pytest.mark.xfail
@pytest.mark.parametrize("label", ["grouping", "quality"])
@pytest.mark.parametrize("vals", [True, [1, 1], np.ones(10)])
def test_pha_column_size(label, vals, make_test_pha):
    """Check we error out if column=label has the wrong size"""

    pha = make_test_pha
    with pytest.raises(DataErr) as de:
        # This does not throw an error
        setattr(pha, label, vals)

    assert str(de.value) == f"size mismatch between channel and {label}"


@pytest.mark.xfail
def test_pha_change_grouping_type(make_test_pha):
    """Check the grouping column is converted to int"""
    pha = make_test_pha
    grp = np.asarray([1.0, -1.0, -1.0, 1.0])
    pha.grouping = grp

    # Since integer values can do an exact check
    assert (pha.grouping == np.asarray([1, -1, -1, 1])).all()
    # At the moment the array is not converted to an int type
    assert pha.grouping.dtype == np.int16


@pytest.mark.xfail
def test_pha_change_quality_type(make_test_pha):
    """Check the quality column is converted to int"""
    pha = make_test_pha
    # technically negative numbers are allowed
    qual = np.asarray([0.0, 2.0, 5.0, -1.0])
    pha.quality = qual

    # Since integer values can do an exact check
    assert (pha.quality == np.asarray([0, 2, 5, -1])).all()
    # At the moment the array is not converted to an int type
    assert pha.quality.dtype == np.int16


@pytest.mark.xfail
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

    # At the moment there is no conversion so the return
    # value matches the input (i.e. vals).
    #
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

    # Only filtering
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx(y[2:6])

    pha.group_counts(3)
    assert pha.get_filter(format="%.1f") == "1.0:11.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(1, 11))

    # Grouped and filtered
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([3, 2])

    assert pha.mask == pytest.approx([True] * 2)
    assert pha.get_mask() == pytest.approx([True] * 10)

    grouping = [1, -1, -1, -1, -1,  1, -1, -1, -1, -1.]
    assert pha.grouping == pytest.approx(grouping)

    quality = [0, 0, 0, 0, 0, 2, 2, 2, 2, 2]
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
    assert pha.mask

    # However, get_mask reflects the quality filter, so is 5 True
    # followed by 5 False.
    #
    mask = [True] * 5 + [False] * 5
    assert pha.get_mask() == pytest.approx(mask)

    # What about the quality fields?
    #
    assert pha.quality == pytest.approx(quality)
    assert pha.quality_filter == pytest.approx([True] * 5 + [False] * 5)

    # Saying all that though, the filter expression does not
    # know we are ignoring channels 6-10.
    #
    # TODO: This is likely a bug.
    #
    assert pha.get_filter(format="%.1f") == "1.0:11.0"
    assert pha.get_noticed_channels() == pytest.approx(np.arange(1, 6))

    # Grouped and quality-filtered (even though get_filter
    # returns 1:11 here).
    #
    assert pha.get_dep(filter=False) == pytest.approx(y)
    assert pha.get_dep(filter=True) == pytest.approx([3])

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

    assert str(de.value) == "unknown plot type 'None', choose 'rate' or 'counts'"


def test_pha_ignore_bad_no_quality():
    """Just check we error out"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    with pytest.raises(DataErr) as de:
        pha.ignore_bad()

    assert str(de.value) == "data set 'dummy' does not specify quality flags"


def test_pha_get_ylabel_yfac0():
    """This does not depend on the backend"""

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    assert pha.plot_fac == 0
    assert pha.get_ylabel() == 'Counts/channel'


@pytest.mark.parametrize("override_plot_backend", [dummy_backend])
def test_pha_get_ylabel_yfac1(override_plot_backend):
    """Basic check

    The label depends on the backend, so we just want the dummy
    backend used here. **UNFORTUNATELY** - either because the
    override_plot_backend fixture is not well written, the current
    approach to setting up the plot backend does not handle it being
    swapped out (e.g. see #1191), or a combination of the two - the
    test doesn't work well if there is a non-dummy backend loaded.

    """

    chans = np.arange(1, 4)
    counts = np.ones_like(chans)
    pha = DataPHA("dummy", chans, counts)

    pha.plot_fac = 1

    # This is ugly - hopefully #1191 will fix this
    #
    ylabel = pha.get_ylabel()
    if backend.__name__.endswith('.dummy_backend'):
        assert ylabel == 'Counts/channel X Channel^1'
    else:
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

    # We could set up channels and counts, but let's not.
    #
    d = DataPHA("dummy", None, None)
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


@pytest.mark.xfail
def test_dataimgint_show(make_dataimgint):
    """Check we can show a basic integrated image data set.

    See issue #1379
    """

    img = make_dataimgint

    # This fails because there's problems getting x0 and x0lo
    # attributes.
    #
    out = str(img).split("\n")

    # I don't know what the correct output should be, so the
    # x0/x0lo order could be wrong.
    #
    assert out[0] == "name      = ival"
    assert out[1] == "x0lo      = Int64[2]"
    assert out[2] == "x1lo      = Int64[2]"
    assert out[3] == "x0hi      = Int64[2]"
    assert out[4] == "x1hi      = Int64[2]"
    assert out[5] == "y         = Int64[2]"
    assert out[6] == "shape     = (2, 1)"
    assert out[7] == "staterror = None"
    assert out[8] == "syserror  = None"
    assert out[9] == "sky       = None"
    assert out[10] == "eqpos     = None"
    assert out[11] == "coord     = logical"
    assert len(out) == 12


@pytest.mark.xfail
def test_dataimgint_x0lo(make_dataimgint):
    assert make_dataimgint.x0lo == pytest.approx([-5, -5])


@pytest.mark.xfail
def test_dataimgint_x1lo(make_dataimgint):
    assert make_dataimgint.x1lo == pytest.approx([10, 11])


@pytest.mark.xfail
def test_dataimgint_x0hi(make_dataimgint):
    assert make_dataimgint.x0hi == pytest.approx([-4, -4])


@pytest.mark.xfail
def test_dataimgint_x1hi(make_dataimgint):
    assert make_dataimgint.x1hi == pytest.approx([11, 12])


@pytest.mark.xfail
def test_dataimgint_get_x0(make_dataimgint):
    """This copies the Data2DInt case but fails"""
    x0 = np.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_dataimgint.get_x0() == x).all()


@pytest.mark.xfail
def test_dataimgint_x0(make_dataimgint):
    """This copies the Data2DInt case but fails"""
    x0 = np.asarray([-5, -5])
    x = (x0 + x0 + 1) / 2

    assert (make_dataimgint.x0 == x).all()


@pytest.mark.xfail
def test_dataimgint_get_x1(make_dataimgint):
    """This copies the Data2DInt case but fails"""
    x1 = np.asarray([10, 11])
    x = (x1 + x1 + 1) / 2

    assert (make_dataimgint.get_x1() == x).all()


@pytest.mark.xfail
def test_dataimgint_x1(make_dataimgint):
    """This copies the Data2DInt case but fails"""
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


@pytest.mark.xfail
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


@pytest.mark.xfail
def test_dataimgint_notice2d(make_dataimgint):
    """basic notice2d call.

    Given that we only have two items the testing is not
    going to be extensive.
    """
    img = make_dataimgint

    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.mask == np.asarray([True, False])).all()


@pytest.mark.xfail
def test_dataimgint_ignore2d(make_dataimgint):
    """basic ignore2d call.

    Given that we only have two items the testing is not
    going to be extensive.
    """
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)", ignore=True)
    assert (img.mask == np.asarray([False, True])).all()


@pytest.mark.xfail
def test_dataimgint_notice2d_get_filter(make_dataimgint):
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert img.get_filter() == 'Rectangle(-100,10,100,11)'


@pytest.mark.xfail
def test_dataimgint_notice2d_get_filter_expr(make_dataimgint):
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert img.get_filter_expr() == 'Rectangle(-100,10,100,11)'


@pytest.mark.xfail
def test_dataimgint_notice2d_get_x0(make_dataimgint):
    """basic notice2d call + get_x0"""
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.get_x0() == np.asarray([-4.5, -4.5])).all()
    assert (img.get_x0(True) == np.asarray([-4.5])).all()


@pytest.mark.xfail
def test_dataimgint_notice2d_get_x1(make_dataimgint):
    """basic notice2d call + get_x1"""
    img = make_dataimgint
    img.notice2d("rect(-100, 10, 100, 11)")
    assert (img.get_x1() == np.asarray([10.5, 11.5])).all()
    assert (img.get_x1(True) == np.asarray([10.5])).all()


@pytest.mark.xfail
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


@pytest.mark.xfail
def test_dataimgint_get_max_pos(make_dataimgint):
    """No idea wy this fails"""
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

    with pytest.raises(DataErr) as de:
        make_dataimgint.get_physical()

    assert str(de.value) == "data set 'ival' does not contain a physical coordinate system"


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


@pytest.mark.xfail
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


@pytest.mark.xfail
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


@pytest.mark.xfail
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
