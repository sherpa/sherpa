#
#  Copyright (C) 2007, 2015, 2016, 2017, 2018, 2019, 2023
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

import tempfile
import os
import sherpa
from sherpa.image import Image, DataImage, ModelImage, RatioImage, \
    ResidImage

from sherpa.utils.testing import requires_ds9

# Create a rectangular array for the tests just to ensure that
# there are no issues with Fortran/C order.
#
_ny = 10
_nx = 13


class Data():
    def __init__(self):
        self.name = None
        self.y = np.arange(0, _ny * _nx).reshape(_ny, _nx) / 2.0
        self.eqpos = None
        self.sky = None

    def get_img(self, model=None):
        if model is not None:
            return (self.y, self.y)
        else:
            return self.y


data = Data()


# The order returned by the xpaget call below is neither C or Fortran
# style. For instance,
#
# >>> y = np.asarray([1,2,3,10,20,30]).reshape(2,3)
# >>> im = Image()
# >>> im.image(y)
# >>> got = im.xpaget("data image 1 1 3 2 yes")
# >>> got
# '3\n20\n30\n1\n2\n10\n'
# >>> im.xpaget("data image 1 1 3 2 no")
# '3,1 = 3\n2,2 = 20\n3,2 = 30\n1,1 = 1\n2,1 = 2\n1,2 = 10\n'
#
def get_arr_from_imager(im, yexp):
    """Return data from image object, using the yexp value for shape/type"""

    assert yexp.ndim == 2

    (ny, nx) = yexp.shape
    dtype = yexp.dtype.type

    def proc(s):
        """Convert 'i,j = z' to a tuple (i, j, z)"""
        # no error checking
        toks = s.split('=')
        z = dtype(toks[1])
        toks = toks[0].split(',')
        i = int(toks[0])
        j = int(toks[1])
        return (i, j, z)

    # There is almost-certainly a better way to do this, but for
    # now be explicit (the array sizes are not expected to be large,
    # so it doesn't need to be efficient).
    out = np.zeros((ny, nx), dtype=dtype)
    d = im.xpaget("data image 1 1 {} {} no".format(nx, ny))
    for l in d.split('\n'):
        if l.strip() == '':
            continue
        i, j, z = proc(l)
        out[j - 1, i - 1] = z

    return out


@requires_ds9
def test_ds9():
    ctor = sherpa.image.ds9_backend.DS9.DS9Win
    im = ctor(sherpa.image.ds9_backend.DS9._DefTemplate, False)
    im.doOpen()
    im.showArray(data.y)
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    assert data_out == pytest.approx(data.y)


@requires_ds9
def test_image():
    im = Image()
    im.image(data.y)
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    assert data_out == pytest.approx(data.y)


@requires_ds9
def test_data_image():
    im = DataImage()
    im.prepare_image(data)
    im.image()
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    assert data_out == pytest.approx(data.y)


@requires_ds9
def test_model_image():
    im = ModelImage()
    im.prepare_image(data, 1)
    im.image()
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    assert data_out == pytest.approx(data.y)


@requires_ds9
def test_ratio_image():
    im = RatioImage()
    im.prepare_image(data, 1)
    im.image()
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    # All values but the first will be 1, because the first
    # model pixel will be zero, and therefore the ratio function
    # reassigns the ratio there to be one.
    expval = np.ones(data.y.shape)
    expval[0, 0] = 0
    assert data_out == pytest.approx(expval)


@requires_ds9
def test_resid_image():
    im = ResidImage()
    im.prepare_image(data, 1)
    im.image()
    data_out = get_arr_from_imager(im, data.y)
    im.xpaset("quit")
    # Return value is all zeros
    assert data_out == pytest.approx(np.zeros_like(data.y))


@requires_ds9
def test_connection_with_x_file():
    """Check that the connection works even if there is a
    file called x (this checks that the xpaset call is properly
    escaped when 'xpaset sherpa [BITPIX=..,x=..,y=..,]' is
    called.
    """

    origdir = os.getcwd()
    dname = tempfile.mkdtemp()
    try:
        os.chdir(dname)

        ofile = 'x'
        with open(ofile, 'w') as fh:
            fh.write('')

        im = Image()
        im.image(data.y)
        data_out = get_arr_from_imager(im, data.y)

    finally:
        os.unlink(ofile)
        os.chdir(origdir)
        os.rmdir(dname)

    im.xpaset("quit")
    assert data_out == pytest.approx(data.y)


@requires_ds9
@pytest.mark.parametrize('coordsys', ['', 'image', 'physical'])
def test_image_getregion(coordsys):
    """Test issue #319.

    Check that image_getregion returns a string, not a byte string.
    Since the data set has no coordinate systems, the image and
    physical systems return the same region. Actually, at present
    it does *not* directly check the image_getregion call, since
    that is in the sherpa.ui module. This is a test of the lower-level
    code, and it's not clear if this is the best idea.
    """

    # This is not ideal.
    from sherpa.image import ds9_backend

    im = ds9_backend.imager
    im.doOpen()
    im.showArray(data.y)

    # Use XPA to set a region in the imager
    imshape = 'image; circle 8.5 7.0 0.8'
    im.xpaset('regions format ds9')
    im.xpaset('regions', data=imshape)

    rval = ds9_backend.get_region(coordsys)

    im.xpaset("quit")

    # extract the coordinates to allow for numeric testing.
    #
    assert rval.startswith('circle(')
    assert rval.endswith(');')
    toks = rval[7:-2].split(',')
    assert len(toks) == 3

    assert float(toks[0]) == pytest.approx(8.5)
    assert float(toks[1]) == pytest.approx(7.0)
    assert float(toks[2]) == pytest.approx(0.8)
