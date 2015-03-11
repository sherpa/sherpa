# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

import numpy
import os
from sherpa.image import *
import sherpa.image.DS9
from sherpa.utils import SherpaTestCase

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

class test_image(SherpaTestCase):
    if (os.environ.has_key("DISPLAY") == True):
        def test_ds9(self):
            im = sherpa.image.DS9.DS9Win(sherpa.image.DS9._DefTemplate, False)
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
    else:
        pass
