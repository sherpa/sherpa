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

import os
import numpy
from sherpa.utils import NoNewAttributesAfterInit, bool_cast

import logging
warning = logging.getLogger(__name__).warning
backend = None

try:
    import ds9_backend as backend

except Exception, e:
    # if DS9 is not found for some reason, like inside gdb
    # give a useful warning and fall back on dummy_backend of noops
    warning("imaging routines will not be available, \n" +
            "failed to import sherpa.image.ds9_backend due to \n'%s: %s'" %
            (type(e).__name__, str(e)))
    import dummy_backend as backend


__all__ = ('Image', 'DataImage', 'ModelImage', 'RatioImage',
           'ResidImage', 'PSFImage', 'PSFKernelImage', 'SourceImage',
           'ComponentModelImage', 'ComponentSourceImage')

__metaclass__ = type


class Image(NoNewAttributesAfterInit):

    def __init__(self):
        NoNewAttributesAfterInit.__init__(self)
        
    def close():
        backend.close()
    close = staticmethod(close)
    
    def delete_frames():
        backend.delete_frames()
    delete_frames = staticmethod(delete_frames)
    
    def get_region(coord):
        return backend.get_region(coord)
    get_region = staticmethod(get_region)
    
    def image(self, array, shape=None, newframe=False, tile=False):
        newframe = bool_cast(newframe)
        tile = bool_cast(tile)
        if shape is None:
            backend.image(array, newframe, tile)
        else:
            backend.image(array.reshape(shape), newframe, tile)
    
    def open():
        backend.open()
    open = staticmethod(open)

    def set_wcs(self, keys):
        backend.wcs( keys )
    
    def set_region(reg, coord):
        backend.set_region(reg, coord)
    set_region = staticmethod(set_region)
    
    def xpaget(arg):
        return backend.xpaget(arg)
    xpaget = staticmethod(xpaget)
    
    def xpaset(arg, data=None):
        return backend.xpaset(arg, data=None)
    xpaset = staticmethod(xpaset)

class DataImage(Image):
    """Image data.

    Attributes
    ----------
    name : str
    y : array_like
       The image data (pixel valuess) as a 2D array.
    eqpos :
       Coordinate transform to the "world" system.
    sky :
       Coordinate transform to the "physical" system.

    """

    def __init__(self):
        self.y = None
        self.eqpos = None
        self.sky = None
        self.name = 'Data'
        Image.__init__(self)

    def __str__(self):
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)
        return (('name   = %s\n' % self.name)+
                ('y      = %s\n' % y)+
                ('eqpos  = %s\n' % self.eqpos)+
                ('sky    = %s\n' % self.sky))
    
    def prepare_image(self, data):
        self.y = data.get_img()
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)
        header = getattr(data, 'header', None)
        if header is not None:
            obj = header.get('OBJECT')
            if obj is not None:
                self.name = str(obj).replace(" ", "_")


    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))


class ModelImage(Image):

    def __init__(self):
        self.name = 'Model'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

    def __str__(self):
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)
        return (('name   = %s\n' % self.name)+
                ('y      = %s\n' % y)+
                ('eqpos  = %s\n' % self.eqpos)+
                ('sky    = %s\n' % self.sky))
    
    def prepare_image(self, data, model):
        self.y = data.get_img(model)
        self.y = self.y[1]
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)

    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))


class SourceImage(ModelImage):
    def __init__(self):
        ModelImage.__init__(self)
        self.name = 'Source'

    def prepare_image(self, data, model):
        #self.y = data.get_img(model)
        #self.y = self.y[1]

        self.y = data.eval_model(model)
        data._check_shape()
        self.y = self.y.reshape(*data.shape)

        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class RatioImage(Image):

    def __init__(self):
        self.name = 'Ratio'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

    def __str__(self):
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)
        return (('name   = %s\n' % self.name)+
                ('y      = %s\n' % y)+
                ('eqpos  = %s\n' % self.eqpos)+
                ('sky    = %s\n' % self.sky))

    def _calc_ratio(self, ylist):
        data = numpy.array(ylist[0])
        model = numpy.asarray(ylist[1])
        bad = numpy.where(model == 0.0)
        data[bad] = 0.0
        model[bad] = 1.0
        return (data / model)

    def prepare_image(self, data, model):
        self.y = data.get_img(model)
        self.y = self._calc_ratio(self.y)
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)

    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))

class ResidImage(Image):

    def __init__(self):
        self.name = 'Residual'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

    def __str__(self):
        y = self.y
        if self.y is not None:
            y = numpy.array2string(self.y, separator=',', precision=4, suppress_small=False)
        return (('name   = %s\n' % self.name)+
                ('y      = %s\n' % y)+
                ('eqpos  = %s\n' % self.eqpos)+
                ('sky    = %s\n' % self.sky))
    
    def _calc_resid(self, ylist):
        return ylist[0] - ylist[1]
            
    def prepare_image(self, data, model):
        self.y = data.get_img(model)
        self.y = self._calc_resid(self.y)
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)

    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))


class PSFImage(DataImage):

    def prepare_image(self, psf, data=None):
        psfdata = psf.get_kernel(data, False)
        DataImage.prepare_image(self, psfdata)
        self.name = psf.kernel.name


class PSFKernelImage(DataImage):

    def prepare_image(self, psf, data=None):
        psfdata = psf.get_kernel(data)
        DataImage.prepare_image(self, psfdata)
        self.name = 'PSF_Kernel'


class ComponentSourceImage(ModelImage):

    def prepare_image(self, data, model):
        ModelImage.prepare_image(self, data, model)
        #self.name = "Source component '%s'" % model.name
        self.name = "Source_component"

class ComponentModelImage(ModelImage):

    def prepare_image(self, data, model):
        ModelImage.prepare_image(self, data, model)
        #self.name = "Model component '%s'" % model.name
        self.name = "Model_component"
