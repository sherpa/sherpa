# 
#  Copyright (C) 2008  Smithsonian Astrophysical Observatory
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
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.astro.utils._wcs import pix2world, world2pix

class WCS(NoNewAttributesAfterInit):

    def __init__(self, name, type, crval, crpix, cdelt,
                 crota=0.0, epoch=2000.0, equinox=2000.0):
        self.name = name
        self.type = type
        self.crval = numpy.asarray(crval, dtype=float)
        self.crpix = numpy.asarray(crpix, dtype=float)
        self.cdelt = numpy.asarray(cdelt, dtype=float)
        self.crota = crota
        self.epoch = epoch
        self.equinox = equinox
        NoNewAttributesAfterInit.__init__(self)

    def __repr__(self):
        return ("<%s Coordinate instance '%s'>" %
                (type(self).__name__, self.name))

    def __str__(self):
        val = [self.name,
               ' crval    = %s' % numpy.array2string(self.crval, separator=',', precision=4, suppress_small=False),
               ' crpix    = %s' % numpy.array2string(self.crpix, separator=',', precision=4, suppress_small=False),
               ' cdelt    = %s' % numpy.array2string(self.cdelt, separator=',', precision=4, suppress_small=False)]

        if self.type == 'WCS':
            val.append(' crota    = %g' % self.crota)
            val.append(' epoch    = %g' % self.epoch)
            val.append(' equinox  = %g' % self.equinox)

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
