#
#  Copyright (C) 2008, 2020-2021, 2024, 2026
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

"""Handle basic WCS support.

This only supports a limited number of 2D transformations.

"""

import logging
from typing import overload

import numpy as np

from sherpa.utils import NoNewAttributesAfterInit, arr2str
from sherpa.utils.types import ArrayType


warning = logging.getLogger(__name__).warning

try:
    from sherpa.astro.utils._wcs import pix2world, world2pix  # type: ignore
except ImportError:
    warning('WCS support is not available')

    def pix2world(*args, **kwargs):
        "Pixel to World Coordinates"
        raise RuntimeError("No WCS support")

    def world2pix(*args, **kwargs):
        "World to Pixel Coordinates"
        raise RuntimeError("No WCS support")


class WCS(NoNewAttributesAfterInit):
    """Represent a World Coordinate System transformation.

    This does not support all possible astronomical transforms.

    Parameters
    ----------
    name
       The name of the transformation (expected to be "world" or
       "physical").
    type
       The transform type (expected to be one of "LINEAR", "WCS",
       or "TAN-P").
    crval
       The reference point in the output coordinate system (a
       two-element sequence).
    crpix
       The reference point in the input coordinate system (a
       two-element sequence).
    cdelt
       The width and height of a pixel in the output coordinate
       system (a two-element sequence).
    crota
       The rotation of the transformation.
    epoch
       The epoch of the transformation (not always used).
    equinox
       The equinox of the transformation (not always used).

    """

    def __init__(self,
                 name: str,
                 type: str,
                 crval: ArrayType,
                 crpix: ArrayType,
                 cdelt: ArrayType,
                 crota: float = 0.0,
                 epoch: float = 2000.0,
                 equinox: float = 2000.0
                 ) -> None:
        self.name = name
        self.type = type
        self.crval = np.asarray(crval, dtype=float)
        self.crpix = np.asarray(crpix, dtype=float)
        self.cdelt = np.asarray(cdelt, dtype=float)
        self.crota = crota
        self.epoch = epoch
        self.equinox = equinox
        super().__init__()

    def __repr__(self):
        return f"<{type(self).__name__} Coordinate instance '{self.name}'>"

    def __str__(self):
        val = [self.name,
               f' crval    = {arr2str(self.crval)}',
               f' crpix    = {arr2str(self.crpix)}',
               f' cdelt    = {arr2str(self.cdelt)}']

        if self.type == 'WCS':
            val.append(f' crota    = {self.crota:g}')
            val.append(f' epoch    = {self.epoch:g}')
            val.append(f' equinox  = {self.equinox:g}')

        return '\n'.join(val)

    @overload
    def apply(self,
              x0: float,
              x1: float
              ) -> tuple[np.float64, np.float64]:
        ...

    @overload
    def apply(self,
              x0: ArrayType,
              x1: ArrayType
              ) -> tuple[np.ndarray, np.ndarray]:
        ...

    def apply(self,
              x0: float | ArrayType,
              x1: float | ArrayType
              ) -> tuple[np.float64 | np.ndarray, np.float64 | np.ndarray]:
        """Convert the input coordinates to the output system.

        Parameters
        ----------
        x0
          Coordinate or coordinates of the first axis.
        x1
          Coordinate or coordinates of the second axis (must match the
          size of x0).

        Returns
        -------
        retval
           The coordinate (or coordinates) of the points in the output
           system as the two-element tuple (c0, c1). If x0 and x1 are
           scalars then c0 and c1 will be scalars otherwise they will
           be arrays.

        Examples
        --------

        >>> crval = [200, -100]
        >>> crpix = [5, 10]
        >>> cdelt = [2, 4]
        >>> c1 = WCS("physical", "LINEAR", crval, crpix, cdelt)
        >>> x, y = c1.apply(5, 10)
        >>> float(x), float(y)
        (200.0, -100.0)
        >>> c1.apply([5, 10], [10, 2])
        (array([200., 210.]), array([-100., -132.]))

        """

        return pix2world(self.type, x0, x1,
                         self.crpix, self.crval, self.cdelt,
                         self.crota, self.equinox, self.epoch)

    @overload
    def invert(self,
               x0: float,
               x1: float
               ) -> tuple[np.float64, np.float64]:
        ...

    @overload
    def invert(self,
               x0: ArrayType,
               x1: ArrayType
               ) -> tuple[np.ndarray, np.ndarray]:
        ...

    def invert(self,
               x0: float | ArrayType,
               x1: float | ArrayType
               ) -> tuple[np.float64 | np.ndarray, np.float64 | np.ndarray]:
        """Convert the output coordinates to the inpyt system.

        Parameters
        ----------
        x0
          Coordinate or coordinates of the first axis.
        x1
          Coordinate or coordinates of the second axis (must match the
          size of x0).

        Returns
        -------
        retval
           The coordinate (or coordinates) of the points in the input
           system as the two-element tuple (c0, c1). If x0 and x1 are
           scalars then c0 and c1 will be scalars otherwise they will
           be arrays.

        Examples
        --------

        >>> crval = [200, -100]
        >>> crpix = [5, 10]
        >>> cdelt = [2, 4]
        >>> c1 = WCS("physical", "LINEAR", crval, crpix, cdelt)
        >>> x, y = c1.invert(200, -100)
        >>> float(x), float(y)
        (5.0, 10.0)
        >>> c1.invert([200, 210], [-100, -132])
        (array([ 5., 10.]), array([10.,  2.]))

        """

        return world2pix(self.type, x0, x1,
                         self.crpix, self.crval, self.cdelt,
                         self.crota, self.equinox, self.epoch)
