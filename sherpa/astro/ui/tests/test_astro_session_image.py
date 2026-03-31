#
#  Copyright (C) 2022, 2025
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
from sherpa.utils.testing import requires_ds9, requires_region, requires_wcs


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


@requires_region
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


@requires_region
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


@requires_region
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


@requires_ds9
def test_get_data_image_dataimg():
    """Check that the OBJECT keyword is used"""

    from sherpa.image import DataImage

    s = AstroSession()
    d = example_data()
    d.header = {"OBJECT": "A big blob", "NOT_OBJ": "foo"}
    s.set_data(d)

    obj = s.get_data_image()
    assert isinstance(obj, DataImage)
    assert obj.name == 'A_big_blob'
    assert obj.eqpos is None
    assert obj.sky is None

    assert d.y.ndim == 1
    y = d.y.reshape((9, 9))
    assert obj.y == pytest.approx(y)

    # sanity check
    assert y[2, 5] == pytest.approx(100.0)


@requires_ds9
def test_get_data_image_dataimg_repeated():
    """If OBJECT is set and then not set, what happens?"""

    s = AstroSession()
    d = example_data()
    d.header = {"OBJECT": "A big blob", "NOT_OBJ": "foo"}
    s.set_data(d)

    # This changes the name field (see test_get_data_image_dataimg)
    _ = s.get_data_image()

    d = example_data()
    d.header = {"NOT_OBJ": "foo"}
    s.set_data(d)

    # What is the name field now?
    obj = s.get_data_image()
    assert obj.name == 'Data'


@requires_ds9
@requires_wcs
def test_send_wcs(tmp_path, check_str):
    """Check we can send WCS information to DS9.

    Validating this is awkward and it has revealed some subtle issues
    with how DS9 interacts with its preference settings (e.g. ICRS vs
    FK5) that can result in subtle changes. Hence the need to be
    explicit.

    """

    from sherpa.astro.io.wcs import WCS

    s = AstroSession()
    s.dataspace2d([4, 3])

    # Set the WCS to something that is easy to check when returned
    # from DS9.
    #
    d = s.get_data()
    d.sky = WCS(type="LINEAR", name="physical", crval=[2000.0, 3000.0],
                crpix=[0.5, 0.5], cdelt=[2.0, 2.0])

    d.eqpos = WCS(type="WCS", name="world", crval=[50.5, -23.4],
                  crpix=[2000.0, 3000.0], cdelt=[-0.036, 0.036],
                  crota=0.0, epoch=2000.0, equinox=2000.0)

    # This needs to call the image methd, not just prepare_image
    s.image_data()

    resp = s.image_xpaget("wcs")
    assert resp == "wcs\n"

    # Is this the best way to check the WCS info that DS9 has?
    # Ensure we use ICRS, since the results are slightly different
    # if the FK5 system is in use.
    #
    s.image_xpaset("wcs sky icrs")
    wcsfile = tmp_path / "test.wcs"
    s.image_xpaset(f"wcs save {wcsfile}")

    # How much does the output here depend on system info (e.g
    # numerical precision) and is the order fixed across versions of
    # DS9?
    #
    cts = wcsfile.read_text()

    lines = ["WCSAXES =                    2 / Number of WCS axes",
             "CRPIX1  =                  0.5 / Reference pixel on axis 1",
             "CRPIX2  =                  0.5 / Reference pixel on axis 2",
             "CRVAL1  =                 50.5 / Value at ref. pixel on axis 1",
             "CRVAL2  =                -23.4 / Value at ref. pixel on axis 2",
             "CTYPE1  = 'RA---TAN'           / Type of co-ordinate on axis 1",
             "CTYPE2  = 'DEC--TAN'           / Type of co-ordinate on axis 2",
             "CDELT1  =               -0.072 / Pixel size on axis 1",
             "CDELT2  =                0.072 / Pixel size on axis 2",
             "RADESYS = 'ICRS    '           / Reference frame for RA/DEC values",
             # "EQUINOX =               2000.0 / [yr] Epoch of reference equinox",
             "WCSAXESP=                    2 / Number of WCS axes",
             "WCSNAMEP= 'PHYSICAL'           / Reference name for the coord. frame",
             "CRPIX1P =                  1.0 / Reference pixel on axis 1",
             "CRPIX2P =                  1.0 / Reference pixel on axis 2",
             "CRVAL1P =               2001.0 / Value at ref. pixel on axis 1",
             "CRVAL2P =               3001.0 / Value at ref. pixel on axis 2",
             "CTYPE1P = 'x       '           / Type of co-ordinate on axis 1",
             "CTYPE2P = 'y       '           / Type of co-ordinate on axis 2",
             "CDELT1P =                  2.0 / Pixel size on axis 1",
             "CDELT2P =                  2.0 / Pixel size on axis 2"
             ]

    # No idea why the two blank/empty lines
    expected = [f"{l:80s}" for l in lines] + ["", ""]
    check_str(cts, expected)
