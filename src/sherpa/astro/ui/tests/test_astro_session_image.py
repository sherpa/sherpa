#
#  Copyright (C) 2022
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

"""Astro-specific tests of the image functionality for Session

Notes
-----
A number of tests send data to an external DS9 process (which will be
started if needed) and then check the values displayed by the process.
These tests are highly-likely to fail if the tests are run in
parallel, such as with

    pytest -n auto

as there is no support in our DS9 interface to either create a new DS9
instance per worker, or to "lock" a DS9 so that the tests can be run
in serial.

"""

import numpy as np

import pytest

from sherpa.astro.ui.utils import Session as AstroSession

from sherpa.astro.data import DataIMG

from sherpa.utils.err import ArgumentTypeErr
from sherpa.utils.testing import requires_ds9


def example_data():
    """Create an example data set.

    It's not entirely clear whether this is a valid object, as we do
    not clearly describe what the units of the intependent axis are
    meant, or allowed, to be: e.g.  must they by logical units with
    any coord system added via WCS, or can they be anything, as used
    here? This means that some of the results (e.g. the filters used)
    may need to be updated if we ever resolve this.

    """

    x1, x0 = np.mgrid[-4:5, 6:15]

    shp = x1.shape

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / (1 + np.sqrt((x0 - 11)**2 + (x1 + 2)**2))

    return DataIMG('example', x0, x1, y, shape=shp)


@requires_ds9
@pytest.mark.parametrize("idval", [None, "foo"])
@pytest.mark.parametrize("funcname", ["notice", "ignore"])
def test_image_xxx2d_not_sent_an_image(funcname, idval):
    """Check we error out

    This test takes a surprising amount of time as we don't
    error out before setting up DS9, so we try to reduce the
    number of parameters.

    If we address #1583 the funcname parameter can be removed, as we
    would expect the image_data call to fail.

    """

    idarg = 1 if idval is None else idval
    s = AstroSession()
    s.load_arrays(idarg, [1, 2, 3], [1, 2, 3])

    # Should this fail? see #1583
    if idval is None:
        s.image_data()
    else:
        s.image_data(idval)

    func = getattr(s, f"{funcname}2d_image")
    with pytest.raises(ArgumentTypeErr,
                       match="^'img' must be a image data set$"):
        if idval is None:
            func()
        else:
            func(idval)


@requires_ds9
@pytest.mark.parametrize("funcname", ["notice", "ignore"])
def test_xxx2d_image_no_region(funcname):
    """What happens when no region has been supplied?"""

    s = AstroSession()

    # This is a bit awkward as we have to send a region to
    # DS9 with XPA so we can check we can read it back
    # with notice2d_image.
    #
    d = example_data()
    s.set_data(d)
    s.image_data()

    func = getattr(s, f"{funcname}2d_image")
    with pytest.raises(ValueError,
                       match="^unable to parse region string: ''$"):
        func()


def send_filters_to_ds9():
    """Filters used by test_xxx2d_image"""

    from sherpa.image import backend
    backend.xpaset("regions", data="circle(10,0,3)")
    backend.xpaset("regions", data="box(8,4,8,6,30)")


FILTER_NTOT = 81
FILTER_NSET = 47
FILTER_A = "Circle(10,0,3)"
FILTER_B = "RotBox(8,4,8,6,30)"


@requires_ds9
def test_notice2d_image():
    """A valid region"""

    s = AstroSession()

    # This is a bit awkward as we have to send a region to
    # DS9 with XPA so we can check we can read it back
    # with notice2d_image.
    #
    s.set_data(example_data())
    s.image_data()

    d = s.get_data()
    assert s.get_filter() == ""
    assert d.mask is True

    send_filters_to_ds9()

    s.notice2d_image()
    assert s.get_filter() == f"{FILTER_A}|{FILTER_B}"

    assert len(d.mask) == FILTER_NTOT
    assert d.mask.sum() == FILTER_NSET


@requires_ds9
def test_ignore2d_image():
    """A valid region

    Unlike test_notice2d_image this uses an explicit
    dataset identifier.
    """

    s = AstroSession()

    # This is a bit awkward as we have to send a region to
    # DS9 with XPA so we can check we can read it back
    # with notice2d_image.
    #
    s.set_data("ex", example_data())
    s.image_data("ex")

    d = s.get_data("ex")
    assert s.get_filter("ex") == ""
    assert d.mask is True

    send_filters_to_ds9()

    s.ignore2d_image("ex")
    assert s.get_filter("ex") == f"Field()&!{FILTER_A}&!{FILTER_B}"

    # 47 is from test_notice2d_image
    assert len(d.mask) == FILTER_NTOT
    assert d.mask.sum() == (FILTER_NTOT - FILTER_NSET)
