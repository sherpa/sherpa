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

import unittest
import numpy
# from numpy.testing import assert_allclose

import tempfile

import os
import sherpa
from sherpa.image import *
from sherpa.utils import SherpaTestCase, requires_ds9

# Create a 10x10 array for the tests.
class Data(object):
    def __init__(self):
        self.name = None
        self.y = numpy.arange(0,(10*10)/2,0.5)
        self.y = self.y.reshape(10,10)
        self.eqpos = None
        self.sky = None

    def get_img(self,model=None):
        if model is not None:
            return (self.y, self.y)
        else:
            return self.y

data = Data()

def get_arr_from_imager(im):
    # DS9 returns data as unordered string
    # Turn it into a 10x10 array
    data_out = im.xpaget("data image 1 1 10 10 yes")
    data_out = data_out.split()
    data_out = numpy.array(data_out).reshape(10,10)
    data_out = numpy.float_(data_out)
    return data_out

# TODO: use numpy.testing.assert_allclose (or related) so that a
#       per-pixel check can be done, rather than an aggregated one.
#       Hmmm, apparently the test is doing something different than
#       I think it is, as the values aren't an exact match. So what
#       exactly is get_arr_from_imager returning?

@requires_ds9
class test_image(SherpaTestCase):
    def test_ds9(self):
        im = sherpa.image.ds9_backend.DS9.DS9Win(sherpa.image.ds9_backend.DS9._DefTemplate, False)
        im.doOpen()
        im.showArray(data.y)
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        self.assertEqualWithinTol((data.y - data_out).sum(), 0.0, 1e-4)

    def test_image(self):
        im = Image()
        im.image(data.y)
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        self.assertEqualWithinTol((data.y - data_out).sum(), 0.0, 1e-4)

    def test_data_image(self):
        im = DataImage()
        im.prepare_image(data)
        im.image()
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        self.assertEqualWithinTol((data.y - data_out).sum(), 0.0, 1e-4)

    def test_model_image(self):
        im = ModelImage()
        im.prepare_image(data, 1)
        im.image()
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        self.assertEqualWithinTol((data.y - data_out).sum(), 0.0, 1e-4)

    def test_ratio_image(self):
        im = RatioImage()
        im.prepare_image(data, 1)
        im.image()
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        # The sum is 99, because the first model pixel
        # will be zero, and therefore the ratio function
        # reassigns the ratio there to be one.
        self.assertEqualWithinTol(data_out.sum(), 99.0, 1e-4)

    def test_resid_image(self):
        im = ResidImage()
        im.prepare_image(data, 1)
        im.image()
        data_out = get_arr_from_imager(im)
        im.xpaset("quit")
        self.assertEqualWithinTol(data_out.sum(), 0.0, 1e-4)

    def test_connection_with_x_file(self):
        """Check that the connection works even if there is a
        file called x (this checks that the xpaset call is properly
        escaped).

        """

        origdir = os.getcwd()
        dname = tempfile.mkdtemp()
        try:
            os.chdir(dname)

            ofile = 'x'
            open(ofile, 'w').write('')

            im = Image()
            im.image(data.y)
            data_out = get_arr_from_imager(im)

        finally:
            os.unlink(ofile)
            os.chdir(origdir)
            os.rmdir(dname)

        im.xpaset("quit")

        # assert_allclose(data.y, data_out, atol=0.0, rtol=1e-6)
        self.assertEqualWithinTol((data.y - data_out).sum(), 0.0, 1e-4)
