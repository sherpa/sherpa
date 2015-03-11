# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
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
from sherpa.models.parameter import Parameter, tinyval
from sherpa.models.model import ArithmeticModel, modelCacher1d
from sherpa.utils.err import ModelErr
from sherpa.utils import SherpaFloat, sao_fcmp

_tol = numpy.finfo(numpy.float32).eps

# Optical Models for SED Analysis
#
# This Sherpa Python module contains optical models for fitting to SEDs.
# These models are Python versions of models found in the Specview
# application for analyzing spectra and SEDs.  These models are meant
# to be used in conjunction with Specview to serve the VAO SED project.
#
# These models work in wavelength space (Angstroms).
#

__all__ = ('AbsorptionEdge', 'AccretionDisk', 'AbsorptionGaussian', 'AbsorptionLorentz', 'EmissionLorentz', 'OpticalGaussian', 'EmissionGaussian', 'AbsorptionVoigt', 'BlackBody', 'Bremsstrahlung', 'BrokenPowerlaw', 'CCM', 'LogAbsorption', 'LogEmission', 'Polynomial', 'Powerlaw', 'Recombination', 'EmissionVoigt', 'XGal', 'FM', 'LMC', 'SM', 'SMC', 'Seaton')

# The speed of light in km/s
c_km = 2.99792458e+5

# Helper function, to do a particular sort of crude interpolation
# from tables supplied for certain extinction curves.  No general
# applicability outside this function, so we do not make it available
# to code outside this module.  Ported from the Specview file
# AbstractExtinction.java for Iris use.  5/26/11 SMD
def _extinct_interp(xtable, etable, x):
    out = numpy.zeros_like(x)
    last = len(xtable) - 1

    for i in xrange(len(x)):
        xval = x[i]
        if (xval <= xtable[0]):
            out[i] = etable[0]
        elif (xval >= xtable[last]):
            out[i] = etable[last]
        else:
            index = 0
            for j in xrange(last+1):
                if (xval < xtable[j]):
                    index = j
                    break
            x1 = xtable[index-1]
            x2 = xtable[index]
            e1 = etable[index-1]
            e2 = etable[index]

            out[i] = ((e2 - e1) / (x2 - x1)) * (xval - x1) + e1

    return out
    
# This model sets in edge (in Angstroms) beyond which absorption
# is a significant feature to the spectrum or SED.
class AbsorptionEdge(ArithmeticModel):

    def __init__(self, name='absorptionedge'):
        self.edgew = Parameter(name, 'edgew', 5000., tinyval, frozen=True, units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)
        self.index = Parameter(name, 'index', 3.0, alwaysfrozen=True,
                               hidden=True)

        ArithmeticModel.__init__(self, name,
                                 (self.edgew, self.tau, self.index))

    # We can turn on model caching with this commented-out feature,
    # if we find we need it.
    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.ones_like(x)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s edgew cannot be zero' % self.name)

        idx = (x <= p[0])
        y[idx] = numpy.exp(-(p[1]*numpy.power(x[idx]/p[0], p[2])))
        return y

# This model is an accretion disk continuum function.
class AccretionDisk(ArithmeticModel):

    def __init__(self, name='accretiondisk'):

        self.ref = Parameter(name, 'ref', 5000., frozen=True, units='angstroms')
        self.beta = Parameter(name, 'beta', 0.5, -10, 10)
        self.ampl = Parameter(name, 'ampl', 1.)
        self.norm = Parameter(name, 'norm', 20000.0, tinyval, alwaysfrozen=True,
                              hidden=True)

        ArithmeticModel.__init__(self, name,
                                 (self.ref, self.beta, self.ampl, self.norm))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 in x:
            raise ValueError('model evaluation failed, ' +
                             'x cannot be zero')

        if p[3] == 0.0:
            raise ValueError('model evaluation failed, ' +
                             'norm cannot be zero')

        return p[2]*numpy.power(x/p[3], -p[1])*numpy.exp(-p[0]/x)


# This model calculates a Gaussian function expressed in
# equivalent width, and models absorption due to this Gaussian.
class AbsorptionGaussian(ArithmeticModel):
    """Absorption Gaussian function expressed in equivalent width."""

    def __init__(self, name='absorptiongaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.ewidth = Parameter(name, 'ewidth', 1.)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.ewidth, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044 ?
        delta = numpy.abs((x - p[1]) / sigma)
        ampl  = p[2] / sigma / 2.50662828  # document this constant

        idx = (delta < self.limit)
        y[idx] = 1.0 - ampl * numpy.exp(- delta[idx] * delta[idx] / 2.0)

        return y

# This model calculates a Lorentzian function expressed in
# equivalent width, and models absorption due to this Lorentzian.
class AbsorptionLorentz(ArithmeticModel):
    """Absorption Lorentz function expressed in equivalent width."""

    def __init__(self, name='absorptionlorentz'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.ewidth = Parameter(name, 'ewidth', 1.)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ewidth))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = (1.0 / x - 1.0 / p[1]) * p[1] * c_km / p[0]
        y = 1.0 + 4.0 * y * y
        y *= 1.571 * p[0] * p[1] / c_km
        y = 1.0 - p[2] / y
        return y

# This model computes a Lorentzian profile for emission features.
class EmissionLorentz(ArithmeticModel):
    """Emission Lorentz function expressed in equivalent width."""

    def __init__(self, name='emissionlorentz'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.kurt = Parameter(name, 'kurt', 2., frozen=True)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.flux, self.kurt))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        sigma = p[0] * p[1] / c_km
        arg = numpy.power(numpy.abs(x - p[1]), p[3]) + sigma/2.0 * sigma/2.0

        arg[arg<1.0e-15] = 1.0e-15
        
        return p[2] * sigma / arg / (numpy.pi*2)

# This model computes an absorption Gaussian feature expressed in
# optical depth.
class OpticalGaussian(ArithmeticModel):
    """Absorption Gaussian function expressed in optical depth."""

    def __init__(self, name='opticalgaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.tau, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044 ?
        delta = numpy.abs((x - p[1]) / sigma)

        idx = (delta < self.limit)
        y[idx] = numpy.exp(-p[2] * numpy.exp(- delta[idx] * delta[idx] / 2.0))

        return y

# This model computes a Gaussian profile for emission features.
class EmissionGaussian(ArithmeticModel):
    """Emission Gaussian function."""

    def __init__(self, name='emissiongaussian'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.skew = Parameter(name, 'skew', 1., tinyval, frozen=True)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.flux,
                                  self.skew, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        if 0.0 == p[3]:
            raise ValueError('model evaluation failed, ' +
                             '%s skew cannot be zero' % self.name)

        y = numpy.zeros_like(x)
        sigma = p[1] * p[0] / 705951.5     # = 2.9979e5 / 2.354820044
        delta = numpy.abs((x - p[1]) / sigma)
        idx = (delta < self.limit)

        arg = - delta * delta / 2.0
        if sao_fcmp(p[3], 1.0, _tol) == 0:
            y[idx] = p[2] * numpy.exp(arg[idx]) / sigma / 2.50662828

        else:
            left = (arg <= p[1])
            arg[left] = numpy.exp(arg[left])
            right = ~left
            arg[right] = numpy.exp(arg[right] / p[3] / p[3])
            y[idx] = 2.0 * p[2] * arg[idx] / sigma / 2.50662828 / (1.0 + p[3])

        return y

# This model computes absorption as a Voigt function -- i.e., with
# a Gaussian core and Lorentzian wings.
class AbsorptionVoigt(ArithmeticModel):
    """Absorption Voigt function expressed in equivalent width."""

    def __init__(self, name='absorptionvoigt'):
        self.center = Parameter(name, 'center', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ew = Parameter(name, 'ew', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")
        self.lg = Parameter(name, 'lg', 1., tinyval, hard_min=tinyval)

        # Create core and wings from Gaussian and Lorentz
        self._core = AbsorptionGaussian()
        self._wings = AbsorptionLorentz()
        
        ArithmeticModel.__init__(self, name, (self.center, self.ew, self.fwhm, self.lg))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        # Combining two absorption components means
        # multiplying them, it appears (at least according
        # to addAbsorption in NarrowBandFunction.java)

        core_pars = numpy.array([ p[2], p[0], p[1] / 2.0, 4.0 ])
        wing_pars = numpy.array([ p[2] * p[3], p[0], p[1] / 2.0 ])
        return self._core.calc(core_pars, x) * self._wings.calc(wing_pars, x)

# This model computes continuum emission as a blackbody function.
class BlackBody(ArithmeticModel):
    """Blackbody model."""

    def __init__(self, name='blackbody'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")

        self._argmin = 1.0e-3
        self._argmax = 1000.0

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        c1 = 1.438786e8
        efactor = c1 / p[2]
        if ((efactor / p[0]) > self._argmax):
            # raise error exp too big
            raise ValueError('model evaluation failed, either temperature or reference wavelength too small')

        numer = p[1] * numpy.power(p[0], 5.0) * (numpy.exp(efactor / p[0]) - 1.0)

        y = numpy.zeros_like(x)
        x0 = numpy.where(x > 0.0)[0]
        if (len(x0) > 0):
            arg = numpy.zeros_like(x)
            arg[x0] = efactor / x[x0]
            denon = numpy.zeros_like(x)
            denon[x0] = numpy.power(x[x0], 5)
            argmin_slice = numpy.where(arg < self._argmin)[0]
            if (len(argmin_slice) > 0): 
                denon[argmin_slice] *= arg[argmin_slice] * (1.0 + 0.5 * arg[argmin_slice])
            arg = numpy.where(arg > self._argmax, self._argmax, arg)

            arg_slice = numpy.where(arg >= self._argmin)[0]
            if (len(arg_slice) > 0):
                denon[arg_slice] *= numpy.exp(arg[arg_slice]) - 1.0

            y[x0] = numer / denon[x0]

        return y

# This model computes continuum emission with the bremsstrahlung function.
class Bremsstrahlung(ArithmeticModel):
    """Bremsstrahlung model."""

    def __init__(self, name='bremsstrahlung'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        return p[1] * numpy.power((p[0] / x), 2) * numpy.exp(-1.438779e8 / x / p[2])

# This model computes continuum emission with a broken power-law;
# that is, the power-law index changes after a break at a particular
# wavelength.
class BrokenPowerlaw(ArithmeticModel):
    """Broken power-law model."""
    
    def __init__(self, name='brokenpowerlaw'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.index1 = Parameter(name, 'index1', 0.1, -10.0, 10.0)
        self.index2 = Parameter(name, 'index2', -0.1, -10.0, 10.0)

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.index1, self.index2))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        arg = x / p[0]
        arg = p[1] * (numpy.where(arg > 1.0, numpy.power(arg, p[3]), numpy.power(arg, p[2])))
        return arg

# This model computes extinction using the function published by
# Cardelli, Clayton, and Mathis
# (ApJ, 1989, vol 345, pp 245)
class CCM(ArithmeticModel):
    """Cardelli, Clayton, and Mathis extinction curve."""
    
    def __init__(self, name='ccm'):
        self.ebv = Parameter(name, 'ebv', 0.5)
        self.r = Parameter(name, 'r', 3.2)

        ArithmeticModel.__init__(self, name, (self.ebv, self.r))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.zeros_like(x)
        y2 = numpy.zeros_like(x)
        y3 = numpy.zeros_like(x)

        a = numpy.zeros_like(x)
        b = numpy.zeros_like(x)

        x = 1000.0 / (x / 10.0)

        # Infrared wavelengths
        xp = numpy.zeros_like(x)
        ir_slice = numpy.where((x >= 0.3) & (x <= 1.1))[0]
        if (len(ir_slice) > 0):
            xp[ir_slice] = numpy.power(x[ir_slice], 1.61)
            a[ir_slice] = 0.574 * xp[ir_slice]
            b[ir_slice] = -0.527 * xp[ir_slice]

        # Optical
        opt_slice = numpy.where((x > 1.1) & (x <= 3.3))[0]
        if (len(opt_slice) > 0):
            y[opt_slice] = x[opt_slice] - 1.82

            a[opt_slice] = 1.0 + 0.17699 * y[opt_slice] - 0.50477 * y[opt_slice] * y[opt_slice] - 0.02427 * numpy.power(y[opt_slice],3) + 0.72085 * numpy.power(y[opt_slice],4) + 0.01979 * numpy.power(y[opt_slice],5) - 0.77530 * numpy.power(y[opt_slice],6) + 0.32999 * numpy.power(y[opt_slice],7)

            b[opt_slice] = 0.0 + 1.41338 * y[opt_slice] + 2.28305 * y[opt_slice] * y[opt_slice] + 1.07233 * numpy.power(y[opt_slice],3) - 5.38434 * numpy.power(y[opt_slice],4) - 0.62551 * numpy.power(y[opt_slice],5) + 5.30260 * numpy.power(y[opt_slice],6) - 2.09002 * numpy.power(y[opt_slice],7)

        # Near-UV
        nuv_slice = numpy.where((x > 3.3) & (x <= 8.0))[0]
        if (len(nuv_slice) > 0):
            a[nuv_slice] = 0.0
            b[nuv_slice] = 0.0

            nuv_slice2 = numpy.where((x >= 5.9) & (x <= 8.0))[0]
            if (len(nuv_slice2) > 0):
                y[nuv_slice2] = x[nuv_slice2] - 5.9
                y2[nuv_slice2] = y[nuv_slice2] * y[nuv_slice2]
                y3[nuv_slice2] = y2[nuv_slice2] * y[nuv_slice2]

                a[nuv_slice2] = -0.04473 * y2[nuv_slice2] - 0.009779 * y3[nuv_slice2]
                b[nuv_slice2] = 0.21300 * y2[nuv_slice2] + .120700 * y3[nuv_slice2]

            a[nuv_slice] = a[nuv_slice] + 1.752 - 0.316 * x[nuv_slice] - 0.104 / (0.341 + numpy.power((x[nuv_slice]-4.67),2))

            b[nuv_slice] = b[nuv_slice] - 3.090 + 1.825 * x[nuv_slice] + 1.206 / (0.263 + numpy.power((x[nuv_slice]-4.62),2))

        # Far-UV
        fuv_slice = numpy.where((x > 8.0) & (x <= 20.0))[0]
        if (len(fuv_slice) > 0):
            y[fuv_slice] = x[fuv_slice] - 8.0
            y2[fuv_slice] = y[fuv_slice] * y[fuv_slice]
            y3[fuv_slice] = y2[fuv_slice] * y[fuv_slice]

            a[fuv_slice] = -1.073 - 0.628 * y[fuv_slice] + 0.137 * y2[fuv_slice] - 0.070 * y3[fuv_slice]
            b[fuv_slice] = 13.670 + 4.257 * y[fuv_slice] - 0.420 * y2[fuv_slice] + 0.374 * y3[fuv_slice]

        # Final extinction curve
        aext = p[1] * a + b
        return numpy.power(10.0,(-0.4 * p[0] * aext))

# This model computes absorption using a Gaussian function expressed
# in optical depth, and using the log of the FWHM.
class LogAbsorption(ArithmeticModel):
    """Log of absorption Gaussian function expressed in optical depth."""

    def __init__(self, name='logabsorption'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, 
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.tau = Parameter(name, 'tau', 0.5)

        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos,
                                              self.tau))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        y = numpy.ones_like(x)

        alpha = 0.69314718 / numpy.log(1.0 + p[0] / 2.9979e5 / 2.0)
        if (alpha <= 1.0):
            alpha = 1.0001

        y = numpy.where(x >= p[1], p[2] * numpy.power((x / p[1]), -alpha), p[2] * numpy.power((x / p[1]), alpha))
        
        return numpy.exp(-y)

# This model computes emission using a Gaussian function expressed
# in optical depth, and using the log of the FWHM.
class LogEmission(ArithmeticModel):
    """Log of the emission Gaussian function."""

    def __init__(self, name='logemission'):

        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval,
                              units="km/s")
        self.pos = Parameter(name, 'pos', 5000., tinyval, frozen=True, units='angstroms')
        self.flux = Parameter(name, 'flux', 1.)
        self.skew = Parameter(name, 'skew', 1., tinyval, frozen=True)
        self.limit = Parameter(name, 'limit', 4., alwaysfrozen=True,
                               hidden=True )

        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.flux,
                                  self.skew, self.limit))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s fwhm cannot be zero' % self.name)

        if 0.0 == p[1]:
            raise ValueError('model evaluation failed, ' +
                             '%s pos cannot be zero' % self.name)

        if 0.0 == p[3]:
            raise ValueError('model evaluation failed, ' +
                             '%s skew cannot be zero' % self.name)

        arg = 0.69314718 / numpy.log (1.0 + p[0] / 2.9979e5 / 2.0);
        if (arg <= 1.0):
            arg = 1.0001

        fmax = (arg - 1.0) * p[2] / p[1] / 2.0

        if (p[3] == 1.0):
            return numpy.where(x >= p[1], fmax * numpy.power((x / p[1]), -arg), fmax * numpy.power((x / p[1]), arg)) 

        arg1 = 0.69314718 / numpy.log (1.0 + p[3] * p[0] / 2.9979e5 / 2.0)
        fmax = (arg - 1.0) * p[2] / p[1] / (1.0 + (arg - 1.0) / (arg1 - 1.0))

        return numpy.where(x <= p[1], fmax * numpy.power((x / p[1]), arg), fmax * numpy.power((x / p[1]), -arg1))

# This model computes continuum emission as a polynomial,
# y = c0 + c1 * (x - offset) + c2 * (x - offset)^2 + c3 * (x - offset)^3 + c4 * (x - offset)^4 + c5 * (x - offset)^5
class Polynomial(ArithmeticModel):
    """Polynomial model."""
    
    def __init__(self, name='polynomial'):
        pars = []
        
        for i in xrange(6):
            pars.append(Parameter(name, 'c%d' % i, 0, frozen=True))
        pars[0].val = 1
        pars[0].frozen = False
        for p in pars:
            setattr(self, p.name, p)

        self.offset = Parameter(name, 'offset', 0, frozen=True)
        pars.append(self.offset)
        ArithmeticModel.__init__(self, name, pars)

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        y = numpy.zeros_like(x)
        xtemp = x - p[6]
        y += p[5]
        for i in [4,3,2,1,0]:
            y = y*xtemp + p[i]

        return y

# This model computes continuum emission using a power-law.
class Powerlaw(ArithmeticModel):
    """Power-law model."""
    
    def __init__(self, name='powerlaw'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.index = Parameter(name, 'index', -0.5, -10.0, 10.0)

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.index))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        arg = x / p[0]
        arg = p[1] * numpy.power(arg, p[2])

        return arg

# This model computes the continuum with an optically thin
# recombination function.
class Recombination(ArithmeticModel):
    """Optically thin recombination model."""
    
    def __init__(self, name='recombination'):
        self.refer = Parameter(name, 'refer', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.ampl = Parameter(name, 'ampl', 1., tinyval, hard_min=tinyval, units="angstroms")
        self.temperature = Parameter(name, 'temperature', 3000., tinyval, hard_min=tinyval, units="Kelvin")
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")

        ArithmeticModel.__init__(self, name, (self.refer, self.ampl, self.temperature, self.fwhm))

    #@modelCacher1d 
    def calc(self, p, x, xhi=None, **kwargs):
        if 0.0 == p[0]:
            raise ValueError('model evaluation failed, ' +
                             '%s refer cannot be zero' % self.name)
        
        x = numpy.asarray(x, dtype=SherpaFloat)
        sigma = p[0] * p[3] / 705951.5; # = 2.9979e5 / 2.354820044
        delta = 1.440e8 * (1.0 / x - 1.0 / p[0]) / p[2]

        return numpy.where(delta < 0.0, p[1] * numpy.exp(-numpy.power((x - p[0]), 2.0) / numpy.power(sigma, 2.0) / 2.0), p[1] * numpy.power((p[0] / x), 2.0) * numpy.exp(-delta))

# This model computes emission as a Voigt function -- i.e., with
# a Gaussian core and Lorentzian wings.
class EmissionVoigt(ArithmeticModel):
    """Emission Voigt function."""

    def __init__(self, name='emissionvoigt'):
        self.center = Parameter(name, 'center', 5000., tinyval, hard_min=tinyval, frozen=True, units="angstroms")
        self.flux = Parameter(name, 'flux', 1.)
        self.fwhm = Parameter(name, 'fwhm', 100., tinyval, hard_min=tinyval, units="km/s")
        self.lg = Parameter(name, 'lg', 1., tinyval, hard_min=tinyval, frozen=True)

        # Create core and wings from Gaussian and Lorentz
        self._core = EmissionGaussian()
        self._wings = EmissionLorentz()
        
        ArithmeticModel.__init__(self, name, (self.center, self.flux, self.fwhm, self.lg))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        # Combining two emission components means
        # adding them, it appears (at least according
        # to addEmission in NarrowBandFunction.java)

        core_pars = numpy.array([ p[2], p[0], p[1], 1.0])
        wing_pars = numpy.array([ p[3]*p[2], p[0], p[1], 2.0])
        return self._core.calc(core_pars, x) + self._wings.calc(wing_pars, x)

# This model computes the extragalactic extinction function of
# Calzetti, Kinney and Storchi-Bergmann, 1994, ApJ, 429, 582
class XGal(ArithmeticModel):
    """Extragalactic extinction function of Calzetti, Kinney and Storchi-Bergmann"""
    
    def __init__(self, name='xgal'):
        self.ebv = Parameter(name, 'ebv', 0.5)

        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)

        if 0.0 in x:
            raise ValueError('model evaluation failed, ' +
                             'x cannot be zero')

        x = 1000.0 / x

        # Formula from paper with zero point moved to (x = 0)
        ext = ((0.011 * x - 0.198) * x + 1.509) * x

        # Normalize the result according to Kailash Sahu's calculations
        ext *= 2.43

        return numpy.power(10.0,(-0.4 * p[0] * ext))

# This model computes the extinction curve for wavelengths for UV
# spectra below 3200 A.  Fitzpatrick and Massa 1988 extinction curve
# with Drude UV bump.
# See Fitzpatrick and Massa (ApJ, 1988, vol. 328, p. 734)

class FM(ArithmeticModel):
    """Extinction curve for UV spectra below 3200 A (Fitzpatrick and Massa 1988)"""
    
    def __init__(self, name='fm'):
        self.ebv = Parameter(name, 'ebv', 0.5)  # E(B-V)
        self.x0 = Parameter(name, 'x0', 4.6)    # Position of Drude bump
        self.width = Parameter(name, 'width', 0.06) # Width of Drude bump 
        self.c1 = Parameter(name, 'c1', 0.2)
        self.c2 = Parameter(name, 'c2', 0.1)
        self.c3 = Parameter(name, 'c3', 0.02)
        self.c4 = Parameter(name, 'c4', 0.1)
        ArithmeticModel.__init__(self, name, (self.ebv, self.x0,
                                              self.width, self.c1, self.c2,
                                              self.c3, self.c4))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        x = 10000. / x
        av = 3.14

        dru = x * x / ( p[2] * p[2] * x * x + numpy.power(( x - p[1] ),2.0))
        dx = x - 5.9
        fuv = 0.5392 * dx * dx + 0.0564 * numpy.power(dx,3.0)
        ext = av + p[3] + p[4] * x + p[5] * dru
        ext = numpy.where(x > 5.9, p[6] * fuv + ext, ext)

        return numpy.power(10.0,( -0.4 * ext * p[0] ))


# This model computes the extinction curve using the
#  LMC extinction curve from Howarth 1983 MNRAS, 203, 301
class LMC(ArithmeticModel):
    """Howarth's LMC extinction curve (Howarth 1983 MNRAS, 203, 301)"""
    
    def __init__(self, name='lmc'):
        self.ebv = Parameter(name, 'ebv', 0.5)
        self._R = 3.1
        
        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        # convert from wavelength in Angstroms to 1/microns
        x = 10000.0 / x
        
        extmag = numpy.zeros_like(x)

        # Infrared - extend optical results linearly to 0 at 1/lambda = 0
        slice1 = numpy.where(x <= 1.83)[0]
        slice2 = numpy.where((x > 1.83) & (x <= 2.75))[0]
        slice3 = numpy.where(x > 2.75)[0]
        x = numpy.where(x > 10.96, 10.96, x)

        if (len(slice1) > 0):
            extmag[slice1] = ((1.86 - 0.48 * x[slice1]) * x[slice1] - 0.1) * x[slice1]

        if (len(slice2) > 0):
            extmag[slice2] = self._R + 2.04 * (x[slice2] - 1.83) + 0.094 * (x[slice2] - 1.83) * (x[slice2] - 1.83)

        # continue out to lambda = 912 A
        if (len(slice3) > 0):
            extmag[slice3] = self._R - 0.236 + 0.462 * x[slice3] + 0.105 * x[slice3] * x[slice3] + 0.454 / ((x [slice3]- 4.557) * (x[slice3] - 4.557) + 0.293 )

        return numpy.power(10.0,( -0.4 * extmag * p[0] ))

# This model computes the galactic interstellar extinction function
# from Savage & Mathis (1979, ARA&A, 17, 73-111)
class SM(ArithmeticModel):
    """Galactic interstellar extinction function from Savage & Mathis (1979, ARA&A, 17, 73-111)"""
    def __init__(self, name='sm'):
        self.ebv = Parameter(name, 'ebv', 0.5)

        self._xtab = numpy.array([0.00,  0.29,  0.45,  0.80,  1.11,  1.43,
                                  1.82,  2.27,  2.50,  2.91,  3.65,  4.00,
                                  4.17,  4.35,  4.57,  4.76,  5.00,  5.26,
                                  5.56,  5.88,  6.25,  6.71,  7.18,  8.00,
                                  8.50,  9.00,  9.50,  10.00])
        self._extab = numpy.array([0.00,  0.16,  0.38,  0.87,  1.50,  2.32,
                                   3.10,  4.10,  4.40,  4.90,  6.20,  7.29,
                                   8.00,  8.87,  9.67,  9.33,  8.62,  8.00,
                                   7.75,  7.87,  8.12,  8.15,  8.49,  9.65,
                                   10.55, 11.55, 12.90, 14.40])
                                  
        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        # convert from wavelength in Angstroms to 1/microns
        x = 10000.0 / x

        extmag = numpy.zeros_like(x)
        extmag = _extinct_interp(self._xtab, self._extab, x)
        return numpy.power(10.0,( -0.4 * extmag * p[0]))

# This model computes the SMC extinction function of
# Prevot et al. 1984, A&A, 132, 389-392
class SMC(ArithmeticModel):
    """Prevot et.al. 1984 extinction curve for the SMC"""
    def __init__(self, name='smc'):
        self.ebv = Parameter(name, 'ebv', 0.5)

        self._xtab = numpy.array([0.00,  0.29,  0.45,  0.80,  1.11,  1.43,
                                  1.82,  2.35,  2.70,  3.22,  3.34,  3.46,
                                  3.60,  3.75,  3.92,  4.09,  4.28,  4.50,
                                  4.73,  5.00,  5.24,  5.38,  5.52,  5.70,
                                  5.88,  6.07,  6.27,  6.48,  6.72,  6.98,
                                  7.23,  7.52,  7.84])
        self._extab = numpy.array([-3.10, -2.94, -2.72, -2.23, -1.60, -0.78,
                                   0.00,  1.00,  1.67,  2.29,  2.65,  3.00,
                                   3.15,  3.49,  3.91,  4.24,  4.53,  5.30,
                                   5.85,  6.38,  6.76,  6.90,  7.17,  7.71,
                                   8.01,  8.49,  9.06,  9.28,  9.84, 10.80,
                                   11.51, 12.52, 13.54])
                                  
        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        # convert from wavelength in Angstroms to 1/microns
        x = 10000.0 / x

        extmag = numpy.zeros_like(x)
        extmag = _extinct_interp(self._xtab, self._extab, x)
        return numpy.power(10.0,( -0.4 * extmag * p[0]))

# This model computes Seaton's interstellar extinction function.
# The formulae are based on an adopted value of R = 3.20.
#
# This function implements Seaton's function as originally
# implemented in STScI's Synphot program.
#
# For wavelengths > 3704 Angstrom, the function interpolates
# linearly in 1/lambda in Seaton's table 3. For wavelengths
# < 3704 Angstrom, the class uses the formulae from Seaton's
# table 2. The formulae match at the endpoints of their respective
# intervals. There is a mismatch of 0.009 mag/ebmv at nu=2.7
# (lambda=3704 Angstrom). Seaton's tabulated value of 1.44 mag at
# 1/lambda = 1.1 may be in error; 1.64 seems more consistent with
# his other values.
#
# Wavelength range allowed is 0.1 to 1.0 microns; outside this
# range, the class extrapolates the function.
#
# References:
#
# lambda < 1000		same as lambda = 1000.
# 1000 < lambda < 3704	Seaton (1979) MNRAS 187,73p.
# 3704 < lambda < 10,000	Nandy (1975) A+A 44, 195. (corrected to R=3.2)
# 10000 < lambda		    extrapolate linearly in 1/lambda
    
class Seaton(ArithmeticModel):
    """Seaton's interstellar extinction function (1979 MNRAS)"""
    def __init__(self, name='seaton'):
        self.ebv = Parameter(name, 'ebv', 0.5)

        self._xtab = numpy.array([0.0,  1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
                                  1.6, 1.7, 1.8, 1.9, 2.0, 2.1,
                                  2.2, 2.3, 2.4, 2.5, 2.6, 2.7])
        self._extab = numpy.array([0.0,   1.36, 1.64, 1.84, 2.04, 2.24, 2.44,
                                   2.66, 2.88, 3.14, 3.36, 3.56, 3.77,
                                   3.96, 4.15, 4.26, 4.40, 4.52, 4.64])
                                  
        ArithmeticModel.__init__(self, name, (self.ebv,))

    #@modelCacher1d
    def calc(self, p, x, xhi=None, **kwargs):
        x = numpy.asarray(x, dtype=SherpaFloat)
        # convert from wavelength in Angstroms to 1/microns
        x = 10000.0 / x

        extmag = numpy.zeros_like(x)

        ir_slice = numpy.where(x <= 1.0)[0]
        opt_slice = numpy.where((x > 1.0) & (x <2.7))[0]
        uv1_slice = numpy.where((x >= 2.7) &  (x < 3.65))[0]
        uv2_slice = numpy.where((x >= 3.65) & (x < 7.14))[0]
        uv3_slice = numpy.where((x >= 7.14) &  (x <= 10.0))[0]
        uv_extra_slice = numpy.where(x > 10.0)[0]

        # Infrared - extend optical results linearly to 0 at 1/lambda = 0
        if (len(ir_slice) > 0):
            extmag[ir_slice] = self._extab[1] * x[ir_slice] * x[ir_slice]

        # Optical - interpolate in Seaton's table 3
        if (len(opt_slice) > 0):
            extmag[opt_slice] = _extinct_interp(self._xtab, self._extab, x[opt_slice])

        # UV - use analytic formulae from Seaton's table 2
        if (len(uv1_slice) > 0):
            extmag[uv1_slice] = 1.56 + 1.048 * x[uv1_slice] + 1.01 / (( x[uv1_slice] - 4.6 ) * ( x[uv1_slice] - 4.6 ) + 0.280)

        # UV again
        if (len(uv2_slice) > 0):
            extmag[uv2_slice] = 2.29 + 0.848 * x[uv2_slice] + 1.01 / (( x[uv2_slice] - 4.6 ) * ( x[uv2_slice] - 4.6 ) + 0.280)

        # and more UV still
        if (len(uv3_slice) > 0):
            extmag[uv3_slice] = 16.17 + x[uv3_slice] * ( -3.20 + 0.2975 * x[uv3_slice] )

        # Extrapolate beyond 1/lambda = 10.0
        if (len(uv_extra_slice) > 0):
            x[uv_extra_slice] = numpy.where(x[uv_extra_slice] < 50.0, x[uv_extra_slice], 50.0)
            extmag[uv_extra_slice] = 16.17 + x[uv_extra_slice] * ( -3.20 + 0.2975 * x[uv_extra_slice] )

        return numpy.power(10.0,(-0.4 * extmag * p[0]))
