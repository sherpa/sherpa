#
#  Copyright (C) 2007, 2015, 2016  Smithsonian Astrophysical Observatory
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
from numpy.testing import assert_allclose

import tempfile

import os
import sherpa
from sherpa.image import Image, DataImage, ModelImage, RatioImage, \
    ResidImage

from sherpa.utils import SherpaTestCase, requires_ds9


# Create a rectangular array for the tests just to ensure that
# there are no issues with Fortran/C order.
#
_ny = 10
_nx = 13


class Data(object):
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
        toks = s.split(b'=')
        z = dtype(toks[1])
        toks = toks[0].split(b',')
        i = int(toks[0])
        j = int(toks[1])
        return (i, j, z)

    # There is almost-certainly a better way to do this, but for
    # now be explicit (the array sizes are not expected to be large,
    # so it doesn't need to be efficient).
    out = np.zeros((ny, nx), dtype=dtype)
    d = im.xpaget("data image 1 1 {} {} no".format(nx, ny))
    for l in d.split(b'\n'):
        if l.strip() == b'':
            continue
        i, j, z = proc(l)
        out[j - 1, i - 1] = z

    return out


_atol = 0.0
_rtol = 1.0e-6


@requires_ds9
class test_image(SherpaTestCase):
    def test_ds9(self):
        ctor = sherpa.image.ds9_backend.DS9.DS9Win
        im = ctor(sherpa.image.ds9_backend.DS9._DefTemplate, False)
        im.doOpen()
        im.showArray(data.y)
        data_out = get_arr_from_imager(im, data.y)
        im.xpaset("quit")
        assert_allclose(data.y, data_out, atol=_atol, rtol=_rtol)

    def test_image(self):
        im = Image()
        im.image(data.y)
        data_out = get_arr_from_imager(im, data.y)
        im.xpaset("quit")
        assert_allclose(data.y, data_out, atol=_atol, rtol=_rtol)

    def test_data_image(self):
        im = DataImage()
        im.prepare_image(data)
        im.image()
        data_out = get_arr_from_imager(im, data.y)
        im.xpaset("quit")
        assert_allclose(data.y, data_out, atol=_atol, rtol=_rtol)

    def test_model_image(self):
        im = ModelImage()
        im.prepare_image(data, 1)
        im.image()
        data_out = get_arr_from_imager(im, data.y)
        im.xpaset("quit")
        assert_allclose(data.y, data_out, atol=_atol, rtol=_rtol)

    def test_ratio_image(self):
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
        assert_allclose(expval, data_out, atol=_atol, rtol=_rtol)

    def test_resid_image(self):
        im = ResidImage()
        im.prepare_image(data, 1)
        im.image()
        data_out = get_arr_from_imager(im, data.y)
        im.xpaset("quit")
        # Return value is all zeros
        assert_allclose(data.y * 0, data_out, atol=_atol, rtol=_rtol)

    def test_connection_with_x_file(self):
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

        assert_allclose(data.y, data_out, atol=_atol, rtol=_rtol)
