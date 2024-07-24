#
#  Copyright (C) 2008, 2020, 2021, 2024
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

import numpy as np

from sherpa.utils import NoNewAttributesAfterInit

warning = logging.getLogger(__name__).warning

try:
    from sherpa.astro.utils._wcs import pix2world, world2pix  # type: ignore
except ImportError:
    warning('WCS support is not available')

    def pix2world(*args, **kwargs):
        raise RuntimeError("No WCS support")

    def world2pix(*args, **kwargs):
        raise RuntimeError("No WCS support")


class WCS(NoNewAttributesAfterInit):

    def __init__(self, name, type, crval, crpix, cdelt,
                 crota=0.0, epoch=2000.0, equinox=2000.0):
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
        def tostr(vals):
            return np.array2string(vals, separator=',', precision=4,
                                   suppress_small=False)

        val = [self.name,
               f' crval    = {tostr(self.crval)}',
               f' crpix    = {tostr(self.crpix)}',
               f' cdelt    = {tostr(self.cdelt)}']

        if self.type == 'WCS':
            val.append(f' crota    = {self.crota:g}')
            val.append(f' epoch    = {self.epoch:g}')
            val.append(f' equinox  = {self.equinox:g}')

        return '\n'.join(val)

    def apply(self, x0, x1):
        x0, x1 = pix2world(self.type, x0, x1,
                           self.crpix, self.crval, self.cdelt,
                           self.crota, self.equinox, self.epoch)
        return (x0, x1)

    def invert(self, x0, x1):
        x0, x1 = world2pix(self.type, x0, x1,
                           self.crpix, self.crval, self.cdelt,
                           self.crota, self.equinox, self.epoch)
        return (x0, x1)
