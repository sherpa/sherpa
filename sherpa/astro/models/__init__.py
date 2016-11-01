# coding: utf-8
#
#  Copyright (C) 2007, 2016  Smithsonian Astrophysical Observatory
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
from sherpa.astro.utils import apply_pileup
from sherpa.utils.err import ModelErr
from sherpa.utils import *
import sherpa.astro.models._modelfcts

import six

__all__ = ('Atten', 'BBody', 'BBodyFreq', 'Beta1D', 'BPL1D', 'Dered', 'Edge',
           'LineBroad', 'Lorentz1D', 'NormBeta1D', 'Schechter',
           'Beta2D', 'DeVaucouleurs2D', 'HubbleReynolds', 'Lorentz2D',
           'JDPileup', 'MultiResponseSumModel', 'Sersic2D', 'Disk2D',
           'Shell2D')


class Atten(ArithmeticModel):
    """Model the attenuation by the Inter-Stellar Medium (ISM).

    This model calculates the transmission of the interstellar medium
    using the description of the ISM absorption of [1]_. It includes
    neutral He autoionization features. Between 1.2398 and 43.655
    Angstroms (i.e. in the 0.28-10 keV range) the model also accounts
    for metals as described in [2]_. It should only be used when the
    independent axis has units of Angstroms.

    Attributes
    ----------
    hcol
        The column density of HI in atoms cm^-2.
    heiRatio
        The ratio of the HeI to HI column densities.
    heiiRatio
        The ratio of the HeII to HI column densities.

    Notes
    -----
    The code uses the best available photoionization cross-sections to
    date from the atomic data literature and combines them in an
    arbitrary mixture of the three ionic species: HI, HeI, and HeII.

    This model provided courtesy of Pat Jelinsky.

    The grid version is evaluated by numerically intgerating the
    function over each bin using a non-adaptive Gauss-Kronrod scheme
    suited for smooth functions [3]_, falling over to a simple
    trapezoid scheme if this fails.

    References
    ----------

    .. [1] Rumph, Bowyer, & Vennes 1994, AJ 107, 2108
           http://adsabs.harvard.edu/abs/1994AJ....107.2108R

    .. [2] Morrison & McCammon 1983, ApJ 270, 119
           http://adsabs.harvard.edu/abs/1983ApJ...270..119M

    .. [3] https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html

    """

    def __init__(self, name='atten'):
        self.hcol = Parameter(name, 'hcol', 1e+20, 1e+17, 1e+24, tinyval)
        self.heiRatio = Parameter(name, 'heiRatio', 0.1, 0, 1)
        self.heiiRatio = Parameter(name, 'heiiRatio', 0.01, 0, 1)
        ArithmeticModel.__init__(self, name,
                                 (self.hcol, self.heiRatio, self.heiiRatio))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        # atten should act like xsphabs, never integrate.
        kwargs['integrate']=False
        return _modelfcts.atten(*args, **kwargs)


class BBody(ArithmeticModel):
    """A one-dimensional Blackbody model.

    This model can be used when the independent axis is in energy
    or wavelength space.

    Attributes
    ----------
    space
        This parameter is not fit (``alwaysfrozen`` is set), and
        should be set to either 0, when the independent axis is
        energy with units of keV, or 1 when the axis is wavelength
        with units of Angstrom.
    kT
        The temperature if the blackbody, in keV.
    ampl
        The amplitude of the blackbody component.

    See Also
    --------
    BBodyFreq

    Notes
    -----
    The blackbody emission is calculated as a function of energy using
    the expression::

        f(E) = ampl * E^2 / (exp(E / kT) - 1)

    where E is the photon energy, and kT is the blackbody temperature
    (both in keV). The amplitude, ampl, is related to the ratio of
    source radius to distance by::

        ampl = (2 * pi / (c^2 * h^3)) (R / d)^2
             = 9.884 x 10^31 (R / d)^2

    with Planck's constant (h) specified in keV-s and the speed of
    light (c) specified in cm/s, and with R and d representing the
    radius of, and distance to, the source respectively.

    There are two conditions when the above equation is not used:

    - if E/kt < 10^-4 then f(E) = ampl * E * kT

    - if E/kT > 60, f(E) = 0.

    """

    def __init__(self, name='bbody'):
        self.space = Parameter(name, 'space', 0, 0, 1, 0, 1,
                               '0 - energy | 1 - wave', alwaysfrozen=True)
        self.kT = Parameter(name, 'kT', 1, 0.1, 1000, 0, 1e+10, 'keV')
        self.ampl = Parameter(name, 'ampl', 1, 1e-20, 1e+20, 1e-20, 1e+20)
        ArithmeticModel.__init__(self, name, (self.space, self.kT, self.ampl))


    def guess(self, dep, *args, **kwargs):
        if args[0][0] > args[0][-1]:
            self.space.val = 1
        Emax = get_peak(dep, *args)
        tMax = Emax/1.594
        kt = {'val':tMax,
              'min':tMax/_guess_ampl_scale,
              'max':tMax*_guess_ampl_scale}
        param_apply_limits(kt, self.kt, **kwargs)

        norm = guess_amplitude(dep, *args)
        modampl = norm['val']*(numpy.exp(Emax/tMax)-1)/numpy.square(Emax)
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bbody(*args, **kwargs)


class BBodyFreq(ArithmeticModel):
    """A one-dimensional Blackbody model (frequency).

    This model can be used when the independent axis is in frequency
    space.

    Attributes
    ----------
    T
        The temperature if the blackbody, in Kelvin.
    ampl
        The amplitude of the blackbody component.

    See Also
    --------
    BBody

    Notes
    -----
    The blackbody emission is calculated as a function of frequency
    (v) using Wien's law (hv >> kT):


        f(v) = ampl * 2 * h * v^3 * exp(-h * v / (k * T)) / c^2

    where T is the blackbody temperature in Kelvin, h is Planck's
    constant, k is Boltzmann's constant, and c the speed of light.
    """

    def __init__(self, name='bbodyfreq'):
        self.T = Parameter(name, 'T', 1e+06, 1000, 1e+10, 1000, 1e+10)
        self.ampl = Parameter(name, 'ampl', 1, 0, hard_min=0)
        ArithmeticModel.__init__(self, name, (self.T, self.ampl))


    def guess(self, dep, *args, **kwargs):
        vmax = get_peak(dep, *args)
        tMax = vmax/5.88e+10
        t = {'val':tMax,
             'min':tMax/_guess_ampl_scale,
             'max':tMax*_guess_ampl_scale}
        norm = guess_amplitude(dep, *args)
        c_cm = 2.99792458e+10
        h_erg = 6.6260693e-27
        factor = numpy.exp(2.82)*numpy.square(c_cm)/h_erg/2.
        modampl = norm['val']*factor/numpy.power(vmax,3.)
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)
        param_apply_limits(t, self.t, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bbodyfreq(*args, **kwargs)


class Beta1D(ArithmeticModel):
    """One-dimensional beta model function.

    The beta model is a Lorentz model with a varying power law.

    Attributes
    ----------
    r0
        The core radius.
    beta
        This parameter controls the slope of the profile at large
        radii.
    xpos
        The reference point of the profile. This is frozen by default.
    ampl
        The amplitude refers to the maximum value of the model, at
        x = xpos.

    See Also
    --------
    Beta2D, Lorentz1D, NormBeta1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * (1 + ((x - xpos) / r0)^2)^(0.5 - 3 * beta)

    The grid version is evaluated by numerically intgerating the
    function over each bin using a non-adaptive Gauss-Kronrod scheme
    suited for smooth functions [1]_, falling over to a simple
    trapezoid scheme if this fails.

    References
    ----------

    .. [1] https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html

    """

    def __init__(self, name='beta1d'):
        self.r0 = Parameter(name, 'r0', 1, tinyval, hard_min=tinyval)
        self.beta = Parameter(name, 'beta', 1, 1e-05, 10, 1e-05, 10)
        self.xpos = Parameter(name, 'xpos', 0, 0, frozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.beta, self.xpos, self.ampl))

    def get_center(self):
        return (self.xpos.val,)

    def set_center(self, xpos, *args, **kwargs):
        self.xpos.set(xpos)

    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        param_apply_limits(pos, self.xpos, **kwargs)

        ref = guess_reference(self.r0.min, self.r0.max, *args)
        param_apply_limits(ref, self.r0, **kwargs)

        norm = guess_amplitude_at_ref(self.r0.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)


    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.beta1d(*args, **kwargs)


class BPL1D(ArithmeticModel):
    """One-dimensional broken power-law function.

    Attributes
    ----------
    gamma1
       The slope of the power law below the break point.
    gamma2
       The slope of the power law above the break point.
    eb
       The position of the break point.
    ref
       The reference position for the amplitude.
    ampl
       The amplitude, defined with respect to the reference position.

    See Also
    --------
    Beta1D, Lorentz1D, NormBeta1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * (x / ref)^(-gamma1)   x <= eb

             = ampl' * (x / ref)^(-gamma2)  otherwise

       ampl' = ampl (eb / ref)^(gamma2 - gamma1)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='bpl1d'):
        self.gamma1 = Parameter(name, 'gamma1', 0, -10, 10)
        self.gamma2 = Parameter(name, 'gamma2', 0, -10, 10)
        self.eb = Parameter(name, 'eb', 100, 0, 1000, 0)
        self.ref = Parameter(name, 'ref', 1, frozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 1e-20)
        ArithmeticModel.__init__(self, name,
                                 (self.gamma1, self.gamma2,
                                  self.eb, self.ref, self.ampl))


    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        param_apply_limits(pos, self.eb, **kwargs)

        ref = guess_reference(self.ref.min, self.ref.max, *args)
        param_apply_limits(ref, self.ref, **kwargs)

        norm = guess_amplitude_at_ref(self.ref.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.bpl1d(*args, **kwargs)


# TODO: what are the units of the independent axis: Angstrom?

class Dered(ArithmeticModel):
    """A de-reddening model.

    The integrate flag of this model should be set to False when
    used with an integrated grid.

    Attributes
    ----------
    rv
        The ratio of total to selective extinction.
    nhgal
        The absorbing column density of H_gal in units of 10^20 cm^-2.

    Notes
    -----
    This dereddening model uses the analytic formula for the mean
    extension law described in [1]_::

        A(lambda) = E(B-V) * (a * rv + b)
                  = 1.086 tau(lambda)

    where tau(lambda) is the wavelength-dependent optical depth::

        I(lambda) = I(0) * exp(-tau(lambda))

    and a and b are computed using wavelength-dependent formulae,
    which are not reproduced here, for the wavelength range 1000
    Angstroms to 3.3 microns. The relationship between the color
    excess and the column density (nhgal) is [2]_::

        E(B-V) = nhgal / 58.0

    where the units of nhgal is 10^20 cm^-2. The value of the ratio
    of total to selective extinction, rv, is initially set to 3.1, the
    standard value for the diffuse ISM. The final model form is::

        I(lambda) = I(0) exp(-nhgal * (a * ev + b) / (58.0 * 1.086)

    This model provided courtesy of Karl Forster.

    References
    ----------

    .. [1] Cardelli, Clayton, & Mathis 1989, ApJ 345, 245
           http://adsabs.harvard.edu/abs/1989ApJ...345..245C

    .. [2] Bohlin, Savage, & Drake 1978, ApJ 224, 132
           http://adsabs.harvard.edu/abs/1978ApJ...224..132B

    """

    def __init__(self, name='dered'):
        self.rv = Parameter(name, 'rv', 10, 1e-10, 1000, tinyval)
        self.nhgal = Parameter(name, 'nhgal', 1e-07, 1e-07, 100000)
        ArithmeticModel.__init__(self, name, (self.rv, self.nhgal))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.dered(*args, **kwargs)


class Edge(ArithmeticModel):
    """Photoabsorption edge model.

    This model can be used when the independent axis is in energy
    or wavelength space.

    Attributes
    ----------
    space
        This parameter is not fit (``alwaysfrozen`` is set), and
        should be set to either 0, when the independent axis is
        energy with units of keV, or 1 when the axis is wavelength
        with units of Angstrom.
    thresh
        The edge position (in energy or wavelength units matching
        the data grid).
    abs
        The absorption coefficient.

    Notes
    -----
    A phenomenological photoabsorption edge model as a function of
    energy::

        f(x) = exp(-abs * (x / thresh)^-3)   if x >= thresh

             = 1.0                           otherwise

    or, as a function of wavelength::

        f(x) = exp(-abs * (x / thresh)^3)    if x <= thresh

             = 1.0                           otherwise

    """

    def __init__(self, name='edge'):
        self.space = Parameter(name, 'space', 0, 0, 1, 0, 1,
                               '0 - energy | 1 - wave')
        self.thresh = Parameter(name, 'thresh', 1, 0, hard_min=0)
        self.abs = Parameter(name, 'abs', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.space, self.thresh, self.abs))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.edge(*args, **kwargs)


# DOC-NOTE:
#    The equation is different to the ahelp file, but matches the code
#
# DOC-TODO:
#    add note about integrated version
#    convince myself that it is doing something useful.
#

class LineBroad(ArithmeticModel):
    """A one-dimensional line-broadening profile.

    Attributes
    ----------
    ampl
        The amplitude of the line.
    rest
        The rest wavelength.
    vsini
        The rotation velocity (v sin(i)), in km/s.

    Notes
    -----
    The model is::

        f(lambda) = 2 * ampl * c * sqrt(x) / (pi * rest * vsini)

                x = 1 - ((lambda - rest) * c / (rest * vsini))^2

    """

    def __init__(self, name='linebroad'):
        self.ampl = Parameter(name, 'ampl', 1, 0, hard_min=0)
        self.rest = Parameter(name, 'rest', 1000, tinyval, hard_min=tinyval)
        self.vsini = Parameter(name, 'vsini', tinyval, tinyval,
                               hard_min=tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.rest, self.vsini))


    def guess(self, dep, *args, **kwargs):
        ref = guess_reference(self.rest.min, self.rest.max, *args)
        param_apply_limits(ref, self.rest, **kwargs)

        norm = guess_amplitude_at_ref(self.rest.val, dep, *args)
        fwhm = get_fwhm(dep, *args)
        c_km = 2.99792458e+5
        if self.rest.val != 0:
            vsini = 2.*fwhm*c_km/numpy.sqrt(3.)/self.rest.val
            vs = {'val': vsini,
                  'min': vsini/_guess_ampl_scale,
                  'max': vsini*_guess_ampl_scale }
            param_apply_limits(vs, self.vsini, **kwargs)

        modampl = norm['val']*numpy.pi*self.vsini.val*self.rest.val/2./c_km
        mod = {'val': modampl,
               'min': modampl/_guess_ampl_scale,
               'max': modampl*_guess_ampl_scale}
        param_apply_limits(mod, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.linebroad(*args, **kwargs)


# DOC-NOTE: for some reason the division in the equation in the notes
#           section confuses sphinx (it thinks it is a section title).
#
class Lorentz1D(ArithmeticModel):
    """One-dimensional normalized Lorentz model function.

    Attributes
    ----------
    fwhm
        The full-width half maximum of the line.
    pos
        The center of the line.
    ampl
        The amplitude refers to the integral of the model.

    See Also
    --------
    Beta1D, NormBeta1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) =                A * fwhm
               --------------------------------------
               2 * pi * (0.25 * fwhm^2 + (x - pos)^2)

           A = ampl / integral f(x) dx

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='lorentz1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, 0, hard_min=0)
        self.pos = Parameter(name, 'pos', 1)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

        norm = guess_amplitude(dep, *args)
        if fwhm != 10:
            aprime = norm['val']*self.fwhm.val*numpy.pi/2.
            ampl = {'val': aprime,
                    'min': aprime/_guess_ampl_scale,
                    'max': aprime*_guess_ampl_scale}
            param_apply_limits(ampl, self.ampl, **kwargs)
        else:
            param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.lorentz1d(*args, **kwargs)


class NormBeta1D(ArithmeticModel):
    """One-dimensional normalized beta model function.

    This is the same model as the ``Beta1D`` model but with a
    different slope parameter and normalisation.

    Attributes
    ----------
    pos
        The center of the line.
    w
        The line width.
    alpha
        The slope of the profile at large radii.
    ampl
        The amplitude refers to the integral of the model.

    See Also
    --------
    Beta1D, Lorentz1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = A * (1 + ((x - pos) / w)^2)^(-alpha)

           A = ampl / integral f(x) dx

    The grid version is evaluated by numerically intgerating the
    function over each bin using a non-adaptive Gauss-Kronrod scheme
    suited for smooth functions [1]_, falling over to a simple
    trapezoid scheme if this fails.

    References
    ----------

    .. [1] https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html

    """

    def __init__(self, name='normbeta1d'):
        self.pos = Parameter(name, 'pos', 0)
        self.width = Parameter(name, 'width', 1, tinyval, hard_min=tinyval)
        self.index = Parameter(name, 'index', 2.5, 0.5, 1000, 0.5)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.pos, self.width, self.index, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        ampl = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        norm = (fwhm['val']*numpy.sqrt(numpy.pi)*
                numpy.exp(lgam(self.index.val-0.5)-lgam(self.index.val)))
        for key in ampl.keys():
            ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.nbeta1d(*args, **kwargs)


class Schechter(ArithmeticModel):
    """One-dimensional Schecter model function.

    This model is for integrated data grids only.

    Attributes
    ----------
    alpha
        The slope of the power-law component.
    ref
        The reference position, which controls the switch from the
        power-law to the exponential.
    norm
        The normalisation of the model.

    Notes
    -----
    The functional form of the model for grids is::

        f(xlo,xhi) = norm * (xlo / ref)^alpha
                          * exp(-xlo / ref)
                          * (xhi - xlo) / ref

    and for points the model is::

        f(x) = 0

    """

    def __init__(self, name='schechter'):
        self.alpha = Parameter(name, 'alpha', 10)
        self.ref = Parameter(name, 'ref', 1)
        self.norm = Parameter(name, 'norm', 1)
        ArithmeticModel.__init__(self, name, (self.alpha, self.ref, self.norm))


    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.norm, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        if not self.integrate:
            raise ModelErr('alwaysint', self.name)
        return _modelfcts.schechter(*args, **kwargs)


class Beta2D(ArithmeticModel):
    """Two-dimensional beta model function.

    The beta model is a Lorentz model with a varying power law.

    Attributes
    ----------
    r0
        The core radius.
    xpos
        The center of the model on the x0 axis.
    ypos
        The center of the model on the x1 axis.
    ellip
        The ellipticity of the model.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.
    alpha
        This parameter controls the slope of the profile at large
        radii.

    See Also
    --------
    Beta1D, DeVaucouleurs2D, HubbleReynolds, Lorentz2D, Sersic2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl * (1 + r(x0,x1)^2)^(-alpha)

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                  r0^2 * (1-ellip)^2

        xoff(x0,x1) = (x0 - xpos) * cos(theta) + (x1 - ypos) * sin(theta)

        yoff(x0,x1) = (x1 - ypos) * cos(theta) - (x0 - xpos) * sin(theta)

    The grid version is evaluated by adaptive multidimensional
    integration scheme on hypercubes using cubature rules, based
    on code from HIntLib ([1]_) and GSL ([2]_).

    References
    ----------

    .. [1] HIntLib - High-dimensional Integration Library
           http://mint.sbg.ac.at/HIntLib/

    .. [2] GSL - GNU Scientific Library
           http://www.gnu.org/software/gsl/

    """

    def __init__(self, name='beta2d'):
        self.r0 = Parameter(name, 'r0', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians', True)
        self.ampl = Parameter(name, 'ampl', 1)
        self.alpha = Parameter(name, 'alpha', 1, -10, 10)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl, self.alpha))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.beta2d(*args, **kwargs)


class DeVaucouleurs2D(ArithmeticModel):
    """Two-dimensional de Vaucouleurs model.

    This is a formulation of the R^(1/4) law introduced by [1]_. It
    is a special case of the ``Sersic2D`` model with ``n=4``,
    as described in [2]_, [3]_, and [4]_.

    Attributes
    ----------
    r0
        The core radius.
    xpos
        The center of the model on the x0 axis.
    ypos
        The center of the model on the x1 axis.
    ellip
        The ellipticity of the model.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Beta2D, HubbleReynolds, Lorentz2D, Sersic2D

    Notes
    -----
    The model used is the same as the ``Sersic2D`` model with ``n=4``.

    References
    ----------

    .. [1] de Vaucouleurs G., 1948, Ann. d’Astroph. 11, 247
           http://adsabs.harvard.edu/abs/1948AnAp...11..247D

    .. [2] http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html

    .. [3] Graham, A. & Driver, S., 2005, PASA, 22, 118
           http://adsabs.harvard.edu/abs/2005PASA...22..118G

    .. [4] Ciotti, L. & Bertin, G., A&A, 1999, 352, 447-451
           http://adsabs.harvard.edu/abs/1999A%26A...352..447C

    """

    def __init__(self, name='devaucouleurs2d'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.devau(*args, **kwargs)


class HubbleReynolds(ArithmeticModel):
    """Two-dimensional Hubble-Reynolds model.

    Attributes
    ----------
    r0
        The core radius.
    xpos
        The center of the model on the x0 axis.
    ypos
        The center of the model on the x1 axis.
    ellip
        The ellipticity of the model.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Beta2D, DeVaucouleurs2D, Lorentz2D, Sersic2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl / (1 + r(x0,x1))^2

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                  r0^2 * (1-ellip)^2

        xoff(x0,x1) = (x0 - xpos) * cos(theta) + (x1 - ypos) * sin(theta)

        yoff(x0,x1) = (x1 - ypos) * cos(theta) - (x0 - xpos) * sin(theta)

    The grid version is evaluated by adaptive multidimensional
    integration scheme on hypercubes using cubature rules, based
    on code from HIntLib ([1]_) and GSL ([2]_).

    References
    ----------

    .. [1] HIntLib - High-dimensional Integration Library
           http://mint.sbg.ac.at/HIntLib/

    .. [2] GSL - GNU Scientific Library
           http://www.gnu.org/software/gsl/

    """


    def __init__(self, name='hubblereynolds'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.hr(*args, **kwargs)


class Lorentz2D(ArithmeticModel):
    """Two-dimensional un-normalised Lorentz function.

    Attributes
    ----------
    fwhm
        The full-width half maximum.
    xpos
        The center of the model on the x0 axis.
    ypos
        The center of the model on the x1 axis.
    ellip
        The ellipticity of the model.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Beta1D, DeVaucouleurs2D, HubbleReynolds, Lorentz1D, Sersic2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl / (1 + 4 * r(x0,x1)^2)

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                 fwhm^2 * (1-ellip)^2

        xoff(x0,x1) = (x0 - xpos) * cos(theta) + (x1 - ypos) * sin(theta)

        yoff(x0,x1) = (x1 - ypos) * cos(theta) - (x0 - xpos) * sin(theta)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='lorentz2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians',frozen=True)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.fwhm, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.lorentz2d(*args, **kwargs)


class JDPileup(ArithmeticModel):
    """A CCD pileup model for the ACIS detectors on Chandra.

    This model is based on the work by John Davis ([1]_). It is
    intended only for modeling the pileup in one-dimensional spectra
    obtained in imaging mode (i.e. no grating data), but can be used
    with the zeroth-order spectrum of a grating data set.

    Attributes
    ----------
    alpha
        The alpha parameter parameterizes "grade migration" in the
        detector and represents the fraction of piled-up events that
        result in a good grade.
    g0
        The probabilty of assigning a grade of zero. This should
        remain frozen.
    f
        The fraction of flux falling into the pileup region. This
        should remain frozen.
    n
        The number of detection cells. This parameter can not be fit.
    ftime
        The frame time in seconds (as given by the ``EXPTIME``
        keyword of the event file). This parameter can not be fit.
    fracexp
        The fractional exposure that the source experienced while
        dithering on the chip (as given by the ``FRACEXPO`` keyword
        in the ARF file). This parameter can not be fit.
    nterms
        The maximum number of photons considered for pileup in a
        single readout frame. This should not be changed from its
        default value of 30. This parameter can not be fit.

    Notes
    -----
    As the pileup model is inherently non-linear, it is strongly
    advised that multiple optimization methods are used to thoroughly
    investigate the search space for the model.

    The alpha parameter should vary with photon energy and detector
    position, but for simplicity it is treated as independent of
    energy and location.

    An example of using this model to fit a Chandra spectrum is
    provided in [2]_.

    References
    ----------

    .. [1] Davis, J, 2001, ApJ, 562, 575-582.
           http://adsabs.harvard.edu/abs/2001ApJ...562..575D

    .. [2] "Using A Pileup Model",
           http://cxc.harvard.edu/sherpa/threads/pileup/

    """

    def __init__(self, name='jdpileup'):
        self.alpha = Parameter(name, 'alpha', 0.5, 0, 1, 0, 1)
        self.g0 = Parameter(name, 'g0', 1, tinyval, 1, tinyval, 1, frozen=True)
        self.f = Parameter(name, 'f', 0.95, 0.9, 1, 0, 1)
        self.n = Parameter(name, 'n', 1, tinyval, 100, tinyval, 2048,
                           alwaysfrozen=True)
        self.ftime = Parameter(name, 'ftime', 3.241, tinyval, 5, tinyval, 100,
                               'sec', alwaysfrozen=True)
        self.fracexp = Parameter(name, 'fracexp', 0.987, 0, 1, 0, 1,
                                 alwaysfrozen=True)
        self.nterms = Parameter(name, 'nterms', 30, 1, 100, 1, 100,
                                alwaysfrozen=True)
        self._results = None
        ArithmeticModel.__init__(self, name,
                                 (self.alpha, self.g0, self.f, self.n,
                                  self.ftime, self.fracexp, self.nterms))

    def __str__(self):
        s = ArithmeticModel.__str__(self)

        if self._results is not None:
            pileup_fractions, integral_ae, num_terms = self._results

            sum = 0.0
            if num_terms > 0:
                sum = pileup_fractions[1]

            sum_piled = pileup_fractions[2:num_terms + 1].sum()
            sum += sum_piled

            pn = numpy.exp(-integral_ae)

            s += '\n\n'

            for i in six.moves.xrange(1, num_terms + 1):
                pn *= integral_ae / float(i)
                s += '   %d: %g  %g\n' % (i, pn,
                                          pileup_fractions[i] / sum)

            s += '   *** pileup fraction: %g' % (sum_piled / sum)

        return s

    def calc(self, p, arf_source, exposure_time, min_energy, max_energy,
             specresp, model):
        (alpha, g0, psf_frac, num_regions, frame_time, fracexpo,
         max_num_terms) = p

        out = apply_pileup(arf_source, exposure_time, int(max_num_terms),
                           min_energy, max_energy, specresp, fracexpo,
                           frame_time, alpha, g0, num_regions, psf_frac, model)
        self._results = out[1:]
        return out[0]


class MultiResponseSumModel(ArithmeticModel):
    pass

class Sersic2D(ArithmeticModel):
    """Two-dimensional Sersic model.

    This is a generalization of the ``DeVaucouleurs2D`` model,
    in which the exponent ``n`` can vary ([1]_, [2]_, and [3]_).

    Attributes
    ----------
    r0
        The core radius.
    xpos
        The center of the model on the x0 axis.
    ypos
        The center of the model on the x1 axis.
    ellip
        The ellipticity of the model.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.
    n
        The Sersic index (n=4 replicates the ``DeVaucouleurs2D``
        model).

    See Also
    --------
    Beta2D, DeVaucouleurs2D, HubbleReynolds, Lorentz2D

    Notes
    -----
    The functional form of the model for points is can be
    expressed as the following::

        f(x0,x1) = ampl * exp(-b(n) * (r(x0,x1)^(1/n) - 1))

            b(n) = 2 * n - 1 / 3 + 4 / (405 * n) + 46 / (25515 * n^2)

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                  r0^2 * (1-ellip)^2

        xoff(x0,x1) = (x0 - xpos) * cos(theta) + (x1 - ypos) * sin(theta)

        yoff(x0,x1) = (x1 - ypos) * cos(theta) - (x0 - xpos) * sin(theta)

    The grid version is evaluated by adaptive multidimensional
    integration scheme on hypercubes using cubature rules, based
    on code from HIntLib ([4]_) and GSL ([5]_).

    References
    ----------

    .. [1] http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html

    .. [2] Graham, A. & Driver, S., 2005, PASA, 22, 118
           http://adsabs.harvard.edu/abs/2005PASA...22..118G

    .. [3] Ciotti, L. & Bertin, G., A&A, 1999, 352, 447-451
           http://adsabs.harvard.edu/abs/1999A%26A...352..447C

    .. [4] HIntLib - High-dimensional Integration Library
           http://mint.sbg.ac.at/HIntLib/

    .. [5] GSL - GNU Scientific Library
           http://www.gnu.org/software/gsl/

    """

    def __init__(self, name='sersic2d'):
        self.r0 = Parameter(name, 'r0', 10, 0, hard_min=0)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999)
        self.theta = Parameter(name, 'theta', 0, -2*numpy.pi, 2*numpy.pi, -2*numpy.pi,
                               4*numpy.pi, 'radians')
        self.ampl = Parameter(name, 'ampl', 1)
        self.n = Parameter(name,'n', 1, .1, 10, 0.01, 100, frozen=True )
        ArithmeticModel.__init__(self, name,
                                 (self.r0, self.xpos, self.ypos, self.ellip,
                                  self.theta, self.ampl, self.n))
        self.cache = 0

    def get_center(self):
        return (self.xpos.val, self.ypos.val)

    def set_center(self, xpos, ypos, *args, **kwargs):
        self.xpos.set(xpos)
        self.ypos.set(ypos)

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        rad = guess_radius(*args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(rad, self.r0, **kwargs)


    def calc(self, *args, **kwargs):
        kwargs['integrate']=bool_cast(self.integrate)
        return _modelfcts.sersic(*args, **kwargs)

### disk2d and shell2d models
### Contributed by Christoph Deil of the HESS project
### Added to CIAO 4.6 for Dec. 2013 release SMD

class Disk2D(ArithmeticModel):
    """Two-dimensional uniform disk model.

    Attributes
    ----------
    xpos
        The center of the disk on the x0 axis.
    ypos
        The center of the disk on the x1 axis.
    ampl
        The amplitude of the signal within the disk.
    r0
        The radius of the disk.

    See Also
    --------
    Shell2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl  if (x0 - xpos)^2 + (x1 - ypos)^2 <= r0^2

                 = 0     otherwise

    For grids, the x0lo,x1lo values are used in the above equation.

    This model was provided by Christoph Deil.
    """

    def __init__(self, name='disk2d'):
        self.xpos = Parameter(name, 'xpos', 0) # p[0]
        self.ypos = Parameter(name, 'ypos', 0) # p[1]
        self.ampl = Parameter(name, 'ampl', 1) # p[2]
        self.r0 = Parameter(name, 'r0', 1, 0) # p[3]
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl, self.r0))

    def calc(self, p, x, y, *args, **kwargs):
        # Compute radii
        r2 = (x - p[0]) ** 2 + (y - p[1]) ** 2

        # Return ampl when r2 <= r0 else return 0
        return numpy.select([r2 <= p[3] ** 2], [p[2]])

# DOC-NOTE: TODO finish the functional form description

class Shell2D(ArithmeticModel):
    """A homogenous spherical 3D shell projected onto 2D.

    Attributes
    ----------
    xpos
        The center of the shell on the x0 axis.
    ypos
        The center of the shell on the x1 axis.
    ampl
        The amplitude.
    r0
        The line-of-sight distance.
    width
        The width of the shell.

    See Also
    --------
    Disk2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl * r(x0,x1)

                 = 0     otherwise

        r(x0,x1) = ?

    For grids, the x0lo,x1lo values are used in the above equation.

    It is not correct for very large shells on the sky.

    This model was provided by Christoph Deil.
    """

    def __init__(self, name='shell2d'):
        self.xpos = Parameter(name, 'xpos', 0) # p[0]
        self.ypos = Parameter(name, 'ypos', 0) # p[1]
        self.ampl = Parameter(name, 'ampl', 1) # p[2]
        self.r0 = Parameter(name, 'r0', 1, 0) # p[3]
        self.width = Parameter(name, 'width', 0.1, 0)
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl, self.r0, self.width))

    def calc(self, p, x, y, *args, **kwargs):
        """Homogeneously emitting spherical shell,
        projected along the z-direction
        (this is not 100% correct for very large shells on the sky)."""
        (xpos, ypos, ampl, r_0, width) = p

        r2 = (x - xpos) * (x - xpos) + (y - ypos) * (y - ypos)
        r_out = r_0 + width
        r_in2, r_out2 = r_0 * r_0, r_out * r_out
        # r_in3, r_out3 = r_in * r_in2, r_out * r_out2
        # We only call abs() in sqrt() to avoid warning messages.
        sphere_out = numpy.sqrt(numpy.abs(r_out2 - r2))
        sphere_in = numpy.sqrt(numpy.abs(r_in2 - r2))
        # Note: for r > r_out 'numpy.select' fills automatically zeros!
        non_normalized = numpy.select([r2 <= r_in2, r2 <= r_out2],
                                      [sphere_out - sphere_in, sphere_out])
        # Computed with Mathematica:
        integral = 2 * numpy.pi / 3 * (r_out ** 3 - r_0 ** 3)
        # integral = 1
        return ampl * integral * non_normalized
