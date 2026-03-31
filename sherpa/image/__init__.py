#
#  Copyright (C) 2007, 2016, 2021, 2024-2025
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
from typing import cast

import numpy as np

from sherpa.astro.io.wcs import WCS
from sherpa.data import Data2D
from sherpa.instrument import PSFModel
from sherpa.models.model import Model
from sherpa.utils import NoNewAttributesAfterInit, bool_cast


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


__all__ = ('Image', 'BaseImage', 'DataImage', 'ModelImage', 'RatioImage',
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

    @staticmethod
    def close() -> None:
        """Stop the image viewer."""
        backend.close()

    @staticmethod
    def delete_frames() -> None:
        """Delete all the frames open in the image viewer."""
        backend.delete_frames()

    @staticmethod
    def get_region(coord: str) -> str:
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

    # This version could be a staticmethod but derived classes can not be.
    def image(self,
              array: np.ndarray,
              shape: tuple[int, ...] | None = None,
              newframe: bool = False,
              tile: bool = False
              ) -> None:
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
    def open() -> None:
        """Start the image viewer."""
        backend.open()

    @staticmethod
    def set_wcs(keys: tuple[WCS | None, WCS | None, str]) -> None:
        """Send the WCS informatiom to the image viewer.

        Parameters
        ----------
        keys
           The eqpos and sky transforms, and the name of the display.

        """
        backend.wcs(keys)

    @staticmethod
    def set_region(reg: str, coord: str) -> None:
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
    def xpaget(arg: str) -> str:
        """Query the image viewer via XPA.

        Retrieve the results of a query to the image viewer.

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
    def xpaset(arg: str, data: str | bytes | None = None) -> None:
        """Send the image viewer a command via XPA.

        Send a command to the image viewer.

        Parameters
        ----------
        arg : str
           A command to send to the image viewer via XPA.
        data : optional
           The data for the command.

        """
        backend.xpaset(arg, data=None)


# This is intended as an internal class. However it is exposed to
# users in an attempt to let them understand the class structure when
# viewing the documentation.
#
class BaseImage(Image):
    """Store the image data to display.

    This is used to separate the base `Image` support from the derived
    user classes like `DataImage` and `ModelImage`.

    """

    name: str = "undefined"
    """The name of the image"""

    def __init__(self) -> None:
        self.y: np.ndarray | None = None
        """The pixel values to display (as a 2D array) or None."""

        self.eqpos: WCS | None = None
        """Optional coordinate transform to the "world" system."""

        self.sky: WCS | None = None
        """Optional coordinate transform to the "physical" system."""

        super().__init__()

    def __str__(self) -> str:
        y = self.y
        if self.y is not None:
            y = np.array2string(self.y, separator=',', precision=4,
                                suppress_small=False)
        return (f'name   = {self.name}\n'
                f'y      = {y}\n'
                f'eqpos  = {self.eqpos}\n'
                f'sky    = {self.sky}\n')

    # As this class is not derived from ABCMeta we can not mark this
    # with @abstractmethod. The reason for defining it here is that it
    # makes it easier to use BaseImage as a type to indicate an image
    # object that provides prepare_image.
    #
    def prepare_image(self, data: Data2D, *args, **kwargs) -> None:
        """Extract and store the pixel values to display."""
        raise NotImplementedError()

    def image(self,
              shape: tuple[int, ...] | None = None,
              newframe: bool = False,
              tile: bool = False
              ) -> None:
        """Send the data to the image viewer to display.

        Parameters
        ----------
        shape
           The shape of the data (optional).
        newframe
           Should the pixels be displayed in a new frame?
        tile
           Should the display be tiled?

        """

        if self.y is None:
            raise DS9Err("prepare_image has not been called")

        super().image(self.y, shape, newframe, tile)
        self.set_wcs((self.eqpos, self.sky, self.name))


class DataImage(BaseImage):
    """Image data."""

    name: str = "Data"

    def prepare_image(self, data: Data2D) -> None:
        """Extract and store the pixel values to display."""

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


class ModelImage(BaseImage):
    """Model data."""

    name: str = "Model"

    def prepare_image(self, data: Data2D, model: Model) -> None:
        """Extract and store the pixel values to display."""

        y = data.get_img(model)
        self.y = y[1]
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class SourceImage(ModelImage):
    """The source model (before convolution) data."""

    name: str = "Source"

    def prepare_image(self, data: Data2D, model: Model) -> None:
        """Extract and store the pixel values to display."""

        # _check_shape ensures that data.shape is not None,
        # which implies that data.eval_model(model) will not
        # be None.
        #
        data._check_shape()
        y = cast(np.ndarray, data.eval_model(model))
        self.y = y.reshape(*cast(tuple[int, ...], data.shape))

        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class RatioImage(BaseImage):
    """The data divide by the model."""

    name: str = "Ratio"

    def _calc_ratio(self,
                    ylist: tuple[np.ndarray, np.ndarray]
                    ) -> np.ndarray:
        data = np.array(ylist[0])
        model = np.asarray(ylist[1])
        bad = np.where(model == 0.0)
        data[bad] = 0.0
        model[bad] = 1.0
        return (data / model)

    def prepare_image(self, data: Data2D, model: Model) -> None:
        """Extract and store the pixel values to display."""

        y = data.get_img(model)
        self.y = self._calc_ratio(y)
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class ResidImage(BaseImage):
    """The data - model image."""

    name: str = "Residual"

    def _calc_resid(self,
                    ylist: tuple[np.ndarray, np.ndarray]
                    ) -> np.ndarray:
        return ylist[0] - ylist[1]

    def prepare_image(self, data: Data2D, model: Model) -> None:
        """Extract and store the pixel values to display."""

        y = data.get_img(model)
        self.y = self._calc_resid(y)
        self.eqpos = getattr(data, 'eqpos', None)
        self.sky = getattr(data, 'sky', None)


class PSFImage(DataImage):
    """The PSF image."""

    def prepare_image(self, psf: PSFModel, data=None) -> None:
        """Extract and store the pixel values to display."""

        psfdata = psf.get_kernel(data, False)
        super().prepare_image(psfdata)
        self.name = psf.kernel.name


class PSFKernelImage(DataImage):
    """The PSF kernel image."""

    name: str = "PSF_Kernel"

    def prepare_image(self, psf: PSFModel, data=None) -> None:
        """Extract and store the pixel values to display."""

        psfdata = psf.get_kernel(data)
        super().prepare_image(psfdata)
        # What is the best name here?
        self.name = "PSF_Kernel"


class ComponentSourceImage(SourceImage):
    """The unconvolved source component."""

    name: str = "Source_component"


class ComponentModelImage(ModelImage):
    """The model component."""

    name: str = "Model_component"
