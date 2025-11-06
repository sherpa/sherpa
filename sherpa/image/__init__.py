#
#  Copyright (C) 2007, 2016, 2021, 2024, 2026
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

"""Support image display with an external display tool.

At present the only supported application is DS9 [DS9]_, which is
connected to via XPA [XPA]_.

References
----------

..  [DS9] SAOImageDS9, "An image display and visualization tool for astronomical data", https://ds9.si.edu/

..  [XPA] "The XPA Messaging System", https://github.com/ericmandel/xpa

"""

import logging

import numpy as np

from sherpa.utils import NoNewAttributesAfterInit, bool_cast, display_fields

warning = logging.getLogger(__name__).warning

try:
    from . import ds9_backend as backend

except Exception as e:
    # if DS9 is not found for some reason, like inside gdb
    # give a useful warning and fall back on dummy_backend of noops
    warning("imaging routines will not be available, \n"
            "failed to import sherpa.image.ds9_backend due to \n"
            "'%s: %s'", type(e).__name__, str(e))
    from . import dummy_backend as backend


__all__ = ('Image', 'DataImage', 'ModelImage', 'RatioImage',
           'ResidImage', 'PSFImage', 'PSFKernelImage', 'SourceImage',
           'ComponentModelImage', 'ComponentSourceImage')


# As with the Plot and Contour classes, the base Image class works
# with explicit arrays but the derived classes work with Sherpa
# objects (Data and Model), and extract the pixel values from
# them. This means that the image method ends up causing issues for
# type checkers, as the derived classes have a different
# signature. There are also issues with the prepare_image call,
# although this is not defined for the base Image class.
#
class Image(NoNewAttributesAfterInit):
    """Base class for sending image data to an external viewer."""

    _fields: list[str] = []
    """The fields to include in the string output.

    Names that end in ! are treated as scalars, otherwise they are
    passed through NumPy's array2string.
    """

    def __str__(self) -> str:
        return display_fields(self, self._fields)

    @staticmethod
    def close():
        """Stop the image viewer."""
        backend.close()

    @staticmethod
    def delete_frames():
        """Delete all the frames open in the image viewer."""
        backend.delete_frames()

    @staticmethod
    def get_region(coord):
        """Return the region defined in the image viewer.

        Parameters
        ----------
        coord : str
           The name of the coordinate system (the empty string means
           to use the current system).

        Returns
        -------
        region : str
           The region, or regions, or the empty string.

        """
        return backend.get_region(coord)

    def image(self,
              array,
              shape=None,
              newframe=False,
              tile=False
              ):
        """Send the data to the image viewer to display.

        Parameters
        ----------
        array
           The pixel values
        shape
           The shape of the data (optional).
        newframe
           Should the pixels be displayed in a new frame?
        tile
           Should the display be tiled?

        """
        newframe = bool_cast(newframe)
        tile = bool_cast(tile)
        if shape is None:
            vals = array
        else:
            vals = array.reshape(shape)

        backend.image(vals, newframe, tile)

    @staticmethod
    def open():
        """Start the image viewer."""
        backend.open()

    def set_wcs(self, keys):
        """Send the WCS informatiom to the image viewer.

        Parameters
        ----------
        keys
           The eqpos and sky transforms, and the name of the display.

        """
        backend.wcs(keys)

    @staticmethod
    def set_region(reg, coord):
        """Set the region to display in the image viewer.

        Parameters
        ----------
        reg : str
           The region to display.
        coord : str
           The name of the coordinate system (the empty string means
           to use the current system).

        """
        backend.set_region(reg, coord)

    @staticmethod
    def xpaget(arg):
        """Return the result of an XPA call to the image viewer.

        Send a query to the image viewer.

        Parameters
        ----------
        arg : str
           A command to send to the image viewer via XPA.

        Returns
        -------
        returnval : str

        """
        return backend.xpaget(arg)

    @staticmethod
    def xpaset(arg, data=None):
        """Return the result of an XPA call to the image viewer.

        Send a command to the image viewer.

        Parameters
        ----------
        arg : str
           A command to send to the image viewer via XPA.
        data : optional
           The data for the command.

        """
        backend.xpaset(arg, data=None)


class DataImage(Image):
    """Image data.

    Attributes
    ----------
    name : str
    y : array_like
       The image data (pixel values) as a 2D array.
    eqpos :
       Coordinate transform to the "world" system.
    sky :
       Coordinate transform to the "physical" system.

    """

    _fields: list[str] = ["name!", "y", "eqpos!", "sky!"]
    """The fields to include in the string output.

    Names that end in ! are displayed directly, otherwise they are
    passed through NumPy's array2string.

    """

    def __init__(self):
        self.y = None
        self.eqpos = None
        self.sky = None
        self.name = 'Data'
        Image.__init__(self)

    def prepare_image(self, data):
        self.y = data.get_img()
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)
        header = getattr(data, 'header', None)

        # Clear out any previous version.
        self.name = "Data"

        if header is None:
            return

        obj = header.get('OBJECT')
        if obj is not None:
            self.name = str(obj).replace(" ", "_")

    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))


class ModelImage(Image):
    """Model data."""

    _fields: list[str] = ["name!", "y", "eqpos!", "sky!"]
    """The fields to include in the string output.

    Names that end in ! are displayed directly, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.name = 'Model'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

    def prepare_image(self, data, model):
        self.y = data.get_img(model)
        self.y = self.y[1]
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)

    def image(self, shape=None, newframe=False, tile=False):
        Image.image(self, self.y, shape, newframe, tile)
        Image.set_wcs(self, (self.eqpos, self.sky, self.name))


class SourceImage(ModelImage):
    """The source model (before convolution) data."""

    def __init__(self):
        ModelImage.__init__(self)
        self.name = 'Source'

    def prepare_image(self, data, model):
        # self.y = data.get_img(model)
        # self.y = self.y[1]

        self.y = data.eval_model(model)
        data._check_shape()
        self.y = self.y.reshape(*data.shape)

        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class RatioImage(Image):
    """The data divide by the model."""

    _fields: list[str] = ["name!", "y", "eqpos!", "sky!"]
    """The fields to include in the string output.

    Names that end in ! are displayed directly, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.name = 'Ratio'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

    def _calc_ratio(self, ylist):
        data = np.array(ylist[0])
        model = np.asarray(ylist[1])
        bad = np.where(model == 0.0)
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
    """The data - model image."""

    _fields: list[str] = ["name!", "y", "eqpos!", "sky!"]
    """The fields to include in the string output.

    Names that end in ! are displayed directly, otherwise they are
    passed through NumPy's array2string.
    """

    def __init__(self):
        self.name = 'Residual'
        self.y = None
        self.eqpos = None
        self.sky = None
        Image.__init__(self)

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
    """The PSF image."""

    def prepare_image(self, psf, data=None):
        psfdata = psf.get_kernel(data, False)
        DataImage.prepare_image(self, psfdata)
        self.name = psf.kernel.name


class PSFKernelImage(DataImage):
    """The PSF kernel image."""

    def prepare_image(self, psf, data=None):
        psfdata = psf.get_kernel(data)
        DataImage.prepare_image(self, psfdata)
        self.name = 'PSF_Kernel'


class ComponentSourceImage(ModelImage):
    """The unconvolved source component."""

    def prepare_image(self, data, model):
        ModelImage.prepare_image(self, data, model)
        # self.name = "Source component '%s'" % model.name
        self.name = "Source_component"


class ComponentModelImage(ModelImage):
    """The model component."""

    def prepare_image(self, data, model):
        ModelImage.prepare_image(self, data, model)
        # self.name = "Model component '%s'" % model.name
        self.name = "Model_component"
