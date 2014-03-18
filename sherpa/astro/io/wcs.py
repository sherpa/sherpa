#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2008)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
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
