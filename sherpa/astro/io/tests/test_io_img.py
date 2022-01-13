#
#  Copyright (C) 2021
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

from sherpa.astro.data import DataIMG
from sherpa.astro import io
from sherpa.utils.testing import requires_data, requires_fits


def backend_is(name):
    """Are we using the given backend?"""
    return io.backend.__name__ == f"sherpa.astro.io.{name}_backend"


@requires_fits
@requires_data
def test_image_write_basic(make_data_path, tmp_path):
    """Check we can write out an image file as a FITS file
    and keep useful data.
    """

    def check_header(obj):
        # check header for a selected set of keywords
        hdr = obj.header

        # selected OGIP keywords
        if not backend_is("crates"):
            assert hdr["HDUNAME"] == "EVENTS_IMAGE"

        assert hdr["HDUCLASS"] == "OGIP"
        assert hdr["HDUCLAS1"] == "EVENTS"
        assert hdr["HDUCLAS2"] == "ACCEPTED"
        assert "HDUCLAS3" not in hdr
        assert hdr["HDUVERS"] == "1.0.0"

        # a few header keywords to check they are handled,
        # including data types (string, float, integer,
        # logical).
        #
        assert hdr["OBJECT"] == "CSC"
        assert hdr["INSTRUME"] == "ACIS"
        assert hdr["GRATING"] == "NONE"
        assert hdr["DETNAM"] == "ACIS-0"
        assert hdr["RAND_SKY"] == pytest.approx(0.5)
        assert hdr["RAND_PI"] == pytest.approx(1.0)
        assert hdr["DATE-OBS"] == "2006-11-25T09:25:17"
        assert isinstance(hdr["CLOCKAPP"], (bool, np.bool_))
        assert hdr["CLOCKAPP"]
        assert hdr["TIMEZERO"] == 0

        assert "EQUINOX" not in hdr

    def check_data(obj):
        """Basic checks of the data"""

        assert obj.shape == (377, 169)
        assert obj.staterror is None
        assert obj.syserror is None

        assert isinstance(obj.sky, io.wcs.WCS)
        assert isinstance(obj.eqpos, io.wcs.WCS)

        assert obj.sky.name == "physical"
        assert obj.sky.type == "LINEAR"
        assert obj.sky.crval == pytest.approx([3061.8101, 4332.6299])
        assert obj.sky.crpix == pytest.approx([0.5, 0.5])
        assert obj.sky.cdelt == pytest.approx([1.0, 1.0])

        assert obj.eqpos.name == "world"
        assert obj.eqpos.type == "WCS"
        assert obj.eqpos.crval == pytest.approx([149.8853, 2.60795])
        assert obj.eqpos.crpix == pytest.approx([4096.5, 4096.5])
        cd = 0.492 / 3600
        assert obj.eqpos.cdelt == pytest.approx([-cd, cd])
        assert obj.eqpos.crota == pytest.approx(0)
        assert obj.eqpos.epoch == pytest.approx(2000.0)
        assert obj.eqpos.equinox == pytest.approx(2000.0)

        assert obj.coord == "logical"

        assert obj.x0.shape == (63713, )
        assert obj.x1.shape == (63713, )
        assert obj.y.shape == (63713, )

        yexp, xexp = np.mgrid[1:378, 1:170]
        xexp = xexp.flatten()
        yexp = yexp.flatten()

        assert obj.x0 == pytest.approx(xexp)
        assert obj.x1 == pytest.approx(yexp)

        vals = obj.y.reshape((377, 169))
        expected = np.asarray([0, 1, 0, 0, 0, 0, 1, 0, 0, 2])
        assert vals[184:194, 83] == pytest.approx(expected)

        # What do we expect for the data types? It is backend-specific
        #
        assert obj.x0.dtype == np.dtype("float64")
        assert obj.x1.dtype == np.dtype("float64")
        if backend_is("crates"):
            assert obj.y.dtype == np.dtype("float64")

        elif backend_is("pyfits"):
            assert obj.y.dtype == np.dtype(">i2")

        else:
            pytest.fail("Unrecognized IO backend")

    infile = make_data_path("acisf07999_000N001_r0035_regevt3_srcimg.fits")
    indata = io.read_image(infile)
    assert isinstance(indata, DataIMG)
    assert indata.name.endswith("/acisf07999_000N001_r0035_regevt3_srcimg.fits")
    check_header(indata)
    check_data(indata)

    # check we can write it out - be explicit with all options
    #
    outfile = tmp_path / "test.img"
    outfile = str(outfile)  # our IO routines don't recognize paths
    io.write_image(outfile, indata, ascii=False, clobber=False)

    outdata = io.read_image(outfile)
    assert isinstance(outdata, DataIMG)
    assert outdata.name.endswith("/test.img")
    check_header(outdata)
    check_data(outdata)
