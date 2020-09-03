#
#  Copyright (C) 2010, 2016, 2018, 2019, 2020
#       Smithsonian Astrophysical Observatory
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

from sherpa.models import Parameter, ArithmeticModel
from .parameter import Parameter, tinyval
from .model import ArithmeticModel, modelCacher1d, CompositeModel, \
    ArithmeticFunctionModel, RegriddableModel2D, RegriddableModel1D
from sherpa.utils.err import ModelErr
from sherpa.utils import bool_cast, get_position, guess_amplitude, \
    guess_amplitude_at_ref, \
    guess_amplitude2d, guess_bounds, guess_fwhm, guess_position, \
    guess_reference, interpolate, linear_interp, param_apply_limits, \
    sao_fcmp

from . import _modelfcts

import logging
warning = logging.getLogger(__name__).warning


__all__ = ('Box1D', 'Const1D', 'Cos', 'Delta1D', 'Erf', 'Erfc', 'Exp', 'Exp10',
           'Gauss1D', 'Log', 'Log10', 'LogParabola', 'NormGauss1D', 'Poisson',
           'Polynom1D', 'PowLaw1D', 'Scale1D', 'Sin', 'Sqrt', 'StepHi1D',
           'StepLo1D', 'Tan', 'Voigt',
           'Box2D', 'Const2D', 'Delta2D', 'Gauss2D', 'SigmaGauss2D',
           'NormGauss2D', 'Polynom2D', 'Scale2D', 'UserModel', 'TableModel',
           'Integrate1D')

DBL_EPSILON = numpy.finfo(numpy.float).eps


class Voigt(RegriddableModel1D):
    """Voigt Profile
    http://publikationen.badw.de/de/003395768
    """

    def __init__(self, name='voigt'):
        # fwhm_g = 2 * alpha
        # fwhm_l - 2 * gamma        
        self.fwhm_g = Parameter(name, 'fwhm_g', 0.05, min=0.0)
        self.fwhm_l = Parameter(name, 'fwhm_l', 0.05, min=0.0)
        self.pos = Parameter(name, 'pos', 0.0)
        self.ampl  = Parameter(name, 'ampl', 1.0)
        ArithmeticModel.__init__(self, name,
                                 (self.fwhm_g, self.fwhm_l, self.pos, self.ampl))
        return

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.wofz(*args, **kwargs)

    
class Box1D(RegriddableModel1D):
    """One-dimensional box function.

    The model is flat between ``xlow`` and ``xhi`` (both limits are
    inclusive), where it is set to the ``ampl`` parameter. Outside
    this range the model is zero.

    Attributes
    ----------
    xlow
        Coordinate of the lower cut-off.
    xhi
        Coordinate of the upper cut-off.
    ampl
        The amplitude of the box.

    See Also
    --------
    Box2D, Const1D, Delta1D, StepLo1D, StepHi1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl if xlow <= x <= xhi
             = 0       otherwise

    and for an integrated grid it is::

        f(lo,hi) = ampl         if lo >= xlow and hi <= xhi
                 = 0            if hi <= xlow or lo >= xhi
                 = ampl * g     where g is the fraction of lo,hi
                                that falls within xlo,xhi

    This behavior is different to how the amplitude is handled in
    other models, such as ``Const1D``.
    """

    def __init__(self, name='box1d'):
        self.xlow = Parameter(name, 'xlow', 0)
        self.xhi = Parameter(name, 'xhi', 0)
        self.ampl = Parameter(name, 'ampl', 1, -1, 1)
        ArithmeticModel.__init__(self, name, (self.xlow, self.xhi, self.ampl))

    def guess(self, dep, *args, **kwargs):
        lo, hi = guess_bounds(args[0])
        norm = guess_amplitude(dep, *args)
        param_apply_limits(lo, self.xlow, **kwargs)
        param_apply_limits(hi, self.xhi, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.box1d(*args, **kwargs)


class Const(ArithmeticModel):
    def __init__(self, name='const'):
        self.c0 = Parameter(name, 'c0', 1)
        ArithmeticModel.__init__(self, name, (self.c0,))

    def guess(self, dep, *args, **kwargs):
        min = dep.min()
        max = dep.max()
        ylo = 0
        if numpy.abs(min - 0) > DBL_EPSILON:
            ylo = min / 100.
            if min < 0: ylo = 0
        yhi = 0
        if numpy.abs(max - 0) > DBL_EPSILON:
            yhi = -10 * max
            if max > 0: yhi = 100 * max
        param_apply_limits({'val': (max + min) / 2.,
                            'min': ylo, 'max': yhi},
                           self.c0, **kwargs)


class Const1D(RegriddableModel1D, Const):
    """A constant model for one-dimensional data.

    Attributes
    ----------
    c0
        The amplitude of the model.

    See Also
    --------
    Box1D, Const2D, Delta1D, Scale1D, StepLo1D, StepHi1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl

    and for an integrated grid it is::

        f(xlo,xhi) = ampl * (xhi - xlo)

    """
    def __init__(self, name='const1d'):
        Const.__init__(self, name)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.const1d(*args, **kwargs)


class Cos(RegriddableModel1D):
    """One-dimensional cosine function.

    Attributes
    ----------
    period
        The period of the cosine, in units of the independent axis.
    offset
        The offset (related to the phase) of the cosine.
    ampl
        The amplitude of the cosine.

    See Also
    --------
    Sin, Tan

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * cos (2 * pi * (x - offset) / period)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='cos'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.cos(*args, **kwargs)


class Delta1D(RegriddableModel1D):
    """One-dimensional delta function.

    The delta function model is only non-zero at a single point
    (or bin for integrated grids).

    Attributes
    ----------
    pos
        The position of the signal.
    ampl
        The amplitude.

    See Also
    --------
    Box1D, Const1D, Delta2D, StepLo1D, StepHi1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl if x == pos
             = 0       otherwise

    and for an integrated grid it is::

        f(lo,hi) = ampl         if lo <= pos <= hi
                 = 0            otherwise

    This behavior is different to how the amplitude is handled in
    other models, such as ``Const1D``.
    """

    def __init__(self, name='delta1d'):
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.pos, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.delta1d(*args, **kwargs)


class Erf(RegriddableModel1D):
    """One-dimensional error function.

    The function is described at [1]_.

    Attributes
    ----------
    ampl
        The amplitude of the model.
    offset
        The offset of the model.
    sigma
        The scaling factor.

    See Also
    --------
    Erfc

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * erf((x - offset) / sigma)

        erf(y) = (2 / sqrt(pi)) Int_0^y exp(-t^2) dt

    and for an integrated grid it is the integral of this over
    the bin.

    References
    ----------

    .. [1] http://mathworld.wolfram.com/Erf.html

    """

    def __init__(self, name='erf'):
        self.ampl = Parameter(name, 'ampl', 1, 0)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.sigma = Parameter(name, 'sigma', 1, 1e-10, 10, tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.offset, self.sigma))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.erf(*args, **kwargs)


class Erfc(RegriddableModel1D):
    """One-dimensional complementary error function.

    The function is described at [1]_.

    Attributes
    ----------
    ampl
        The amplitude of the model.
    offset
        The offset of the model.
    sigma
        The scaling factor.

    See Also
    --------
    Erf

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * erfc((x - offset) / sigma)

        erfc(y) = (2 / sqrt(pi)) Int_y^infinity exp(-t^2) dt

    and for an integrated grid it is the integral of this over
    the bin.

    References
    ----------

    .. [1] http://mathworld.wolfram.com/Erfc.html

    """

    def __init__(self, name='erfc'):
        self.ampl = Parameter(name, 'ampl', 1, 0)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.sigma = Parameter(name, 'sigma', 1, 1e-10, 10, tinyval)
        ArithmeticModel.__init__(self, name,
                                 (self.ampl, self.offset, self.sigma))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.erfc(*args, **kwargs)


class Exp(RegriddableModel1D):
    """One-dimensional exponential function.

    Attributes
    ----------
    offset
        The offset of the model.
    coeff
        The scaling factor.
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp10, Log, Log10, LogParabola, Sqrt

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * exp(coeff * (x - offset))

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='exp'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.exp(*args, **kwargs)


class Exp10(RegriddableModel1D):
    """One-dimensional exponential function, base 10.

    Attributes
    ----------
    offset
        The offset of the model.
    coeff
        The scaling factor.
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp, Log, Log10, LogParabola, Sqrt

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * 10^(coeff * (x - offset))

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='exp10'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.exp10(*args, **kwargs)


class Gauss1D(RegriddableModel1D):
    """One-dimensional gaussian function.

    Attributes
    ----------
    fwhm
        The Full-Width Half Maximum of the gaussian. It is related to
        the sigma value by: FWHM = sqrt(8 * log(2)) * sigma.
    pos
        The center of the gaussian.
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Gauss2D, NormGauss1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * exp(-4 * log(2) * (x - pos)^2 / fwhm^2)

    and for an integrated grid it is the integral of this over
    the bin.

    Examples
    --------

    Compare the gaussian and normalized gaussian models:

    >>> m1 = sherpa.models.basic.Gauss1D()
    >>> m2 = sherpa.models.basic.NormGauss1D()
    >>> m1.pos, m2.pos = 10, 10
    >>> m1.ampl, m2.ampl = 10, 10
    >>> m1.fwhm, m2.fwhm = 5, 5
    >>> m1(10)
    10.0
    >>> m2(10)
    1.8788745573993026
    >>> m1.fwhm, m2.fwhm = 1, 1
    >>> m1(10)
    10.0
    >>> m2(10)
    9.394372786996513

    The normalised version will sum to the amplitude when given
    an integrated grid - i.e. both low and high edges rather than
    points - that covers all the signal (and with a bin size a lot
    smaller than the FWHM):

    >>> m1.fwhm, m2.fwhm = 12.2, 12.2
    >>> grid = np.arange(-90, 110, 0.01)
    >>> glo, ghi = grid[:-1], grid[1:]
    >>> m1(glo, ghi).sum()
    129.86497637060958
    >>> m2(glo, ghi).sum()
    10.000000000000002

    """

    def __init__(self, name='gauss1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.gauss1d(*args, **kwargs)


class Log(RegriddableModel1D):
    """One-dimensional natural logarithm function.

    Attributes
    ----------
    offset
        The offset of the model.
    coeff
        The scaling factor.
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp, Exp10, Log10, LogParabola, Sqrt

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * log (coeff * (x - offset))

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='log'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.log(*args, **kwargs)


class Log10(RegriddableModel1D):
    """One-dimensional logarithm function, base 10.

    Attributes
    ----------
    offset
        The offset of the model.
    coeff
        The scaling factor.
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp, Exp10, Log, LogParabola, Sqrt

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * log_10 (coeff * (x - offset))

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='log10'):
        self.offset = Parameter(name, 'offset', 0)
        self.coeff = Parameter(name, 'coeff', -1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name,
                                 (self.offset, self.coeff, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.log10(*args, **kwargs)


class LogParabola(RegriddableModel1D):
    """One-dimensional log-parabolic function.

    Attributes
    ----------
    ref
        The reference point for the normalization.
    c1
        The power-law index (gamma).
    c2
        The curvature of the parabola (beta).
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp, Exp10, Log, Log10, Sqrt

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * (x / ref) ^ (-c1 - c2 * log_10 (x / ref))

    The grid version is evaluated by numerically intgerating the
    function over each bin using a non-adaptive Gauss-Kronrod scheme
    suited for smooth functions [1]_, falling over to a simple
    trapezoid scheme if this fails.

    References
    ----------

    .. [1] https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html

    """

    def __init__(self, name='logparabola'):
        self.ref = Parameter(name, 'ref', 1, alwaysfrozen=True)
        self.c1 = Parameter(name, 'c1', 1)
        self.c2 = Parameter(name, 'c2', 1)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.ref, self.c1,
                                              self.c2, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.logparabola(*args, **kwargs)


_gfactor = numpy.sqrt(numpy.pi / (4 * numpy.log(2)))


class NormGauss1D(RegriddableModel1D):
    """One-dimensional normalised gaussian function.

    Attributes
    ----------
    fwhm
        The Full-Width Half Maximum of the gaussian. It is related to
        the sigma value by: FWHM = sqrt(8 * log(2)) * sigma.
    pos
        The center of the gaussian.
    ampl
        The amplitude refers to the integral of the model over the
        range -infinity to infinity.


    See Also
    --------
    Gauss1D, NormGauss2D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * exp(-4 * log(2) * (x - pos)^2 / fwhm^2)
               ----------------------------------------------
                       sqrt(pi / (4 * log(2))) * fwhm

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='normgauss1d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.pos = Parameter(name, 'pos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.fwhm, self.pos, self.ampl))

    def get_center(self):
        return (self.pos.val,)

    def set_center(self, pos, *args, **kwargs):
        self.pos.set(pos)

    def guess(self, dep, *args, **kwargs):
        ampl = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(pos, self.pos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

        # Apply normalization factor to guessed amplitude
        norm = numpy.sqrt(numpy.pi / _gfactor) * self.fwhm.val
        for key in ampl.keys():
            if ampl[key] is not None:
                ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.ngauss1d(*args, **kwargs)


class Poisson(RegriddableModel1D):
    """One-dimensional Poisson function.

    A model expressing the ratio of two Poisson distributions of mean
    mu, one for which the random variable is x, and the other for
    which the random variable is equal to mu itself.

    Attributes
    ----------
    mean
        The mean of the first distribution.
    ampl
        The amplitude of the model.

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * mean! exp((x - mean) * log(mean)) / x!

    The grid version is evaluated by numerically intgerating the
    function over each bin using a non-adaptive Gauss-Kronrod scheme
    suited for smooth functions [1]_, falling over to a simple
    trapezoid scheme if this fails.

    References
    ----------

    .. [1] https://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html

    """

    def __init__(self, name='poisson'):
        self.mean = Parameter(name, 'mean', 1, 1e-05, hard_min=tinyval)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.mean, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        pos = get_position(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(pos, self.mean, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.poisson(*args, **kwargs)


class Polynom1D(RegriddableModel1D):
    """One-dimensional polynomial function of order 8.

    The maximum order of the polynomial is 8. The default setting has
    all parameters frozen except for ``c0``, which means that the
    model acts as a constant.

    Attributes
    ----------
    c0
        The constant term.
    c1
        The amplitude of the (x-offset) term.
    c2
        The amplitude of the (x-offset)^2 term.
    c3
        The amplitude of the (x-offset)^3 term.
    c4
        The amplitude of the (x-offset)^4 term.
    c5
        The amplitude of the (x-offset)^5 term.
    c6
        The amplitude of the (x-offset)^6 term.
    c7
        The amplitude of the (x-offset)^7 term.
    c8
        The amplitude of the (x-offset)^8 term.
    offset
        There is a degeneracy between ``c0`` and ``offset``, so it
        is recommended that at least one of these remains frozen.


    See Also
    --------
    PowLaw1D, Polynom2D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = sum_(i=0)^(i=8) c_i * (x - offset)^i

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='polynom1d'):
        pars = []

        for i in range(9):
            pars.append(Parameter(name, 'c%d' % i, 0, frozen=True))
        pars[0].val = 1
        pars[0].frozen = False
        for p in pars:
            setattr(self, p.name, p)

        self.offset = Parameter(name, 'offset', 0, frozen=True)
        pars.append(self.offset)

        ArithmeticModel.__init__(self, name, pars)

    def guess(self, dep, *args, **kwargs):
        xmin = args[0].min()
        xmax = args[0].max()
        ymin = dep.min()
        ymax = dep.max()
        dydx = (ymax - ymin) / (xmax - xmin)
        dydx2 = (ymax - ymin) / ((xmax - xmin) * (xmax - xmin))

        xlo = 0
        if numpy.abs(xmin - 0) >= DBL_EPSILON:
            xlo = -xmin
            if xmin < 0:
                xlo = xmin

        xhi = 0
        if numpy.abs(xmax - 0) >= DBL_EPSILON:
            xhi = -xmax
            if xmax > 0:
                xhi = xmax

        ylo = 0
        if numpy.abs(ymin - 0) >= DBL_EPSILON:
            ylo = -ymin
            if ymin < 0:
                ylo = ymin

        yhi = 0
        if numpy.abs(ymax - 0) >= DBL_EPSILON:
            yhi = -ymax
            if ymax > 0:
                yhi = ymax

        c0 = {'val': (ymax + ymin) / 2.0, 'min': ylo, 'max': yhi}
        c1 = {'val': 0.0, 'min': -100 * dydx, 'max': 100 * dydx}
        c2 = {'val': 0.0, 'min': -100 * dydx2, 'max': 100 * dydx2}
        c3 = {'val': 0.0, 'min': ylo, 'max': yhi}
        off = {'val': 0.0, 'min': xlo, 'max': xhi}

        param_apply_limits(c0, self.c0, **kwargs)
        param_apply_limits(c1, self.c1, **kwargs)
        param_apply_limits(c2, self.c2, **kwargs)
        param_apply_limits(c3, self.c3, **kwargs)
        param_apply_limits(c3, self.c4, **kwargs)
        param_apply_limits(c3, self.c5, **kwargs)
        param_apply_limits(c3, self.c6, **kwargs)
        param_apply_limits(c3, self.c7, **kwargs)
        param_apply_limits(c3, self.c8, **kwargs)
        param_apply_limits(off, self.offset, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.poly1d(*args, **kwargs)


class PowLaw1D(RegriddableModel1D):
    """One-dimensional power-law function.

    It is assumed that the independent axis is positive at all points.

    Attributes
    ----------
    gamma
        The photon index of the power law.
    ref
        As the reference point is degenerate with the amplitude, the
        ``alwaysfrozen`` attribute of the reference point is set so
        that it can not be varied during a fit.
    ampl
        The amplitude of the model.

    See Also
    --------
    Polynom1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * (x / ref)^(-gamma)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='powlaw1d'):
        self.gamma = Parameter(name, 'gamma', 1, -10, 10)
        self.ref = Parameter(name, 'ref', 1, alwaysfrozen=True)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.gamma, self.ref, self.ampl))

    def guess(self, dep, *args, **kwargs):
        ref = guess_reference(self.ref.min, self.ref.max, *args)
        param_apply_limits(ref, self.ref, **kwargs)
        norm = guess_amplitude_at_ref(self.ref.val, dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, pars, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        if kwargs['integrate']:
            # avoid numerical issues with C pow() function close to zero,
            # 0.0 +- ~1.e-14.  PowLaw1D integrated has multiple calls to
            # pow(X, 1.0 - gamma).  So gamma values close to 1.0 +- 1.e-10
            # should be be 1.0 to avoid errors propagating in the calculated
            # model.
            if sao_fcmp(pars[0], 1.0, 1.e-10) == 0:
                # pars { gamma, ref, ampl }
                pars[0] = 1.0

        return _modelfcts.powlaw(pars, *args, **kwargs)


class Scale1D(Const1D):
    """A constant model for one-dimensional data.

    This is a specialized form of the ``Const1D`` model where the
    ``integrate`` flag is turned off.

    Attributes
    ----------
    c0
        The amplitude of the model.

    See Also
    --------
    Box1D, Const1D, Delta1D, Scale2D, StepLo1D, StepHi1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl

    and for an integrated grid it is::

        f(xlo,xhi) = ampl * (xhi - xlo)

    This later case will only be used if the ``integrate``
    attribute is set to ``True``.
    """

    def __init__(self, name='scale1d'):
        Const1D.__init__(self, name)
        self.integrate = False


class Sin(RegriddableModel1D):
    """One-dimensional sine function.

    Attributes
    ----------
    period
        The period of the sine, in units of the independent axis.
    offset
        The offset (related to the phase) of the sine.
    ampl
        The amplitude of the sine.

    See Also
    --------
    Cos, Tan

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * sin (2 * pi * (x - offset) / period)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='sin'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.sin(*args, **kwargs)


class Sqrt(RegriddableModel1D):
    """One-dimensional square root function.

    Attributes
    ----------
    offset
        The offset of the model.
    ampl
        The amplitude of the model.

    See Also
    --------
    Exp, Exp10, Log, Log10, LogParabola

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * sqrt (x - offset)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='sqrt'):
        self.offset = Parameter(name, 'offset', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.offset, self.ampl))

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.sqrt(*args, **kwargs)


class StepHi1D(RegriddableModel1D):
    """One-dimensional step function.

    The model is flat above ``xcut``, where it is set to the ``ampl``
    parameter, and zero below this.

    Attributes
    ----------
    xcut
        The position of the step.
    ampl
        The amplitude of the step.

    See Also
    --------
    Box1D, Const1D, Delta1D, StepLo1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl if x >= xcut
             = 0       otherwise

    and for an integrated grid it is::

        f(xlo,xhi) = ampl * (xhi - xlo)  if xlo >= xcut
                   = ampl * (xhi - xcut) if xlo < xcut and xhi >= xcut
                   = 0                   if xhi < xcut

    """

    def __init__(self, name='stephi1d'):
        self.xcut = Parameter(name, 'xcut', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.xcut, self.ampl))

    def guess(self, dep, *args, **kwargs):
        cut = guess_bounds(args[0], False)
        norm = guess_amplitude(dep, *args)
        param_apply_limits(cut, self.xcut, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.stephi1d(*args, **kwargs)


class StepLo1D(RegriddableModel1D):
    """One-dimensional step function.

    The model is flat below ``xcut``, where it is set to the ``ampl``
    parameter, and zero above this.

    Attributes
    ----------
    xcut
        The position of the step.
    ampl
        The amplitude of the step.

    See Also
    --------
    Box1D, Const1D, Delta1D, StepHi1D

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl if x <= xcut
             = 0       otherwise

    and for an integrated grid it is::

        f(xlo,xhi) = ampl * (xhi - xlo)  if xhi <= xcut
                   = ampl * (xcut - xlo) if xlo <= xcut and xhi > xcut
                   = 0                   if xlo > xcut

    """

    def __init__(self, name='steplo1d'):
        self.xcut = Parameter(name, 'xcut', 0)
        self.ampl = Parameter(name, 'ampl', 1, 0)
        ArithmeticModel.__init__(self, name, (self.xcut, self.ampl))

    def guess(self, dep, *args, **kwargs):
        cut = guess_bounds(args[0], False)
        norm = guess_amplitude(dep, *args)
        param_apply_limits(cut, self.xcut, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.steplo1d(*args, **kwargs)


class Tan(RegriddableModel1D):
    """One-dimensional tan function.

    Attributes
    ----------
    period
        The period of the tangent, in units of the independent axis.
    offset
        The offset (related to the phase) of the tangent.
    ampl
        The amplitude of the tangent.

    See Also
    --------
    Cos, Sin

    Notes
    -----
    The functional form of the model for points is::

        f(x) = ampl * tan (2 * pi * (x - offset) / period)

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='tan'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))

    def guess(self, dep, *args, **kwargs):
        norm = guess_amplitude(dep, *args)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.tan(*args, **kwargs)


class Box2D(RegriddableModel2D):
    """Two-dimensional box function.

    The model is flat between the limits, where it is set to the
    ``ampl`` parameter. Outside this range the model is zero.

    Attributes
    ----------
    xlow
        The lower edge of the box (``x0`` axis).
    xhi
        The upper edge of the box (``x0`` axis).
    ylow
        The lower edge of the box (``x1`` axis).
    yhi
        The upper edge of the box (``x1`` axis).
    ampl
        The amplitude of the box.

    See Also
    --------
    Box1D, Const2D, Delta2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl if xlow <= x0 <= xhi
                           ylow <= x1 <= yhi
                 = 0       otherwise

    and for an integrated grid it is::

        f(x0lo,x0hi,x1lo,x1hi)
                 = 0            if x0hi <= xlow or x0lo >= xhi
                                or x1hi <= ylow or x1lo >= yhi
                 = ampl * g     where g is the fraction of the pixel
                                that falls within the region

    This behavior is different to how the amplitude is handled in
    other models, such as ``Const2D``.
    """

    def __init__(self, name='box2d'):
        self.xlow = Parameter(name, 'xlow', 0)
        self.xhi = Parameter(name, 'xhi', 0)
        self.ylow = Parameter(name, 'ylow', 0)
        self.yhi = Parameter(name, 'yhi', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.xlow, self.xhi, self.ylow, self.yhi,
                                  self.ampl))
        self.cache = 0

    def guess(self, dep, *args, **kwargs):
        xlo, xhi = guess_bounds(args[0])
        ylo, yhi = guess_bounds(args[1])
        norm = guess_amplitude2d(dep, *args)
        param_apply_limits(xlo, self.xlow, **kwargs)
        param_apply_limits(xhi, self.xhi, **kwargs)
        param_apply_limits(ylo, self.ylow, **kwargs)
        param_apply_limits(yhi, self.yhi, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.box2d(*args, **kwargs)


class Const2D(RegriddableModel2D, Const):
    """A constant model for two-dimensional data.

    Attributes
    ----------
    c0
        The amplitude of the model.

    See Also
    --------
    Box2D, Const1D, Delta2D, Scale2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0, x1) = ampl

    and for an integrated grid it is::

        f(x0lo, x1lo, x0hi, x1hi = ampl * (x0hi - x0lo)
                                        * (x1hi - x1lo)

    """

    def __init__(self, name='const2d'):
        Const.__init__(self, name)
        self.cache = 0

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.const2d(*args, **kwargs)


class Scale2D(Const2D):
    """A constant model for two-dimensional data.

    This is a specialized form of the ``Const2D`` model where the
    ``integrate`` flag is turned off.

    Attributes
    ----------
    c0
        The amplitude of the model.

    See Also
    --------
    Box2D, Const2D, Delta2D, Scale1D

    Notes
    -----
    The functional form of the model for points is::

        f(x0, x1) = ampl

    and for an integrated grid it is::

        f(x0lo, x1lo, x0hi, x1hi) = ampl * (x0hi - x0lo)
                                         * (x1hi - x1lo)

    This later case will only be used if the ``integrate``
    attribute is set to ``True``.
    """

    def __init__(self, name='scale2d'):
        Const2D.__init__(self, name)
        self.integrate = False
        self.cache = 0


class Delta2D(RegriddableModel2D):
    """Two-dimensional delta function.

    The delta function model is only non-zero at a single point
    (or pixel for integrated grids).

    Attributes
    ----------
    xpos
        The coordinate of the signal on the x0 axis.
    ypos
        The coordinate of the signal on the x1 axis.
    ampl
        The amplitude.

    See Also
    --------
    Box2D, Const2D, Delta1D

    Notes
    -----
    The functional form of the model for points is::

        f(x0, x1) = ampl if x0 == xpos and x1 == ypos
             = 0         otherwise

    and for an integrated grid it is::

        f(x0lo, x1lo, x0hi, x1hi) = ampl if x0lo <= xpos <= x0hi
                                            x1lo <= ypos <= x1hi
                                  = 0    otherwise

    This behavior is different to how the amplitude is handled in
    other models, such as ``Const2D``.
    """

    def __init__(self, name='delta2d'):
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.xpos, self.ypos, self.ampl))
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
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.delta2d(*args, **kwargs)


class Gauss2D(RegriddableModel2D):
    """Two-dimensional gaussian function.

    Attributes
    ----------
    fwhm
        The Full-Width Half Maximum of the gaussian along the major
        axis. It is related to the sigma value by:
        FWHM = sqrt(8 * log(2)) * sigma.
    xpos
        The center of the gaussian on the x0 axis.
    ypos
        The center of the gaussian on the x1 axis.
    ellip
        The ellipticity of the gaussian.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Gauss1D, NormGauss2D, SigmaGauss2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl * exp(-4 * log(2) * r(x0,x1)^2)

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                fwhm^2 * (1-ellip)^2

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

    def __init__(self, name='gauss2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2 * numpy.pi, 2 * numpy.pi,
                               -2 * numpy.pi, 4 * numpy.pi, 'radians',
                               frozen=True)
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
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.gauss2d(*args, **kwargs)


class SigmaGauss2D(Gauss2D):
    """Two-dimensional gaussian function (varying sigma).

    Attributes
    ----------
    sigma_a
        The sigma of the gaussian along the major axis.
    sigma_b
        The sigma of the gaussian along the minor axis.
    xpos
        The center of the gaussian on the x0 axis.
    ypos
        The center of the gaussian on the x1 axis.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the maximum peak of the model.

    See Also
    --------
    Gauss2D, NormGauss2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = ampl * exp(-r(x0,x1)^2 / 2)

        r(x0,x1)^2 = xoff(x0,x1)^2 + yoff(x0,x1)^2
                     -------------   -------------
                       sigma_a^2       sigma_b^2

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

    def __init__(self, name='sigmagauss2d'):
        self.sigma_a = Parameter(name, 'sigma_a', 10, tinyval, hard_min=tinyval)
        self.sigma_b = Parameter(name, 'sigma_b', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.theta = \
            Parameter(name, 'theta', 0, -2 * numpy.pi, 2 * numpy.pi,
                      -2 * numpy.pi, 4 * numpy.pi, 'radians', frozen=True)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name,
                                 (self.sigma_a, self.sigma_b, self.xpos,
                                  self.ypos, self.theta, self.ampl))
        self.cache = 0

    def guess(self, dep, *args, **kwargs):
        xpos, ypos = guess_position(dep, *args)
        norm = guess_amplitude2d(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(fwhm, self.sigma_a, **kwargs)
        param_apply_limits(fwhm, self.sigma_b, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.sigmagauss2d(*args, **kwargs)


class NormGauss2D(RegriddableModel2D):
    """Two-dimensional normalised gaussian function.

    Attributes
    ----------
    fwhm
        The Full-Width Half Maximum of the gaussian along the major
        axis. It is related to the sigma value by:
        FWHM = sqrt(8 * log(2)) * sigma.
    xpos
        The center of the gaussian on the x0 axis.
    ypos
        The center of the gaussian on the x1 axis.
    ellip
        The ellipticity of the gaussian.
    theta
        The angle of the major axis. It is in radians, measured
        counter-clockwise from the X0 axis (i.e. the line X1=0).
    ampl
        The amplitude refers to the integral of the model over the
        range -infinity to infinity for both axes.

    See Also
    --------
    Gauss2D, NormGauss1D, SigmaGauss2D

    Notes
    -----
    The functional form of the model for points is::

        f(x0,x1) = 4 * log(2) * ampl * exp(-4 * log(2) * r(x0,x1)^2)
                   -------------------------------------------------
                       pi * fwhm * fwhm * sqrt(1 - ellip * ellip)

        r(x0,x1)^2 = xoff(x0,x1)^2 * (1-ellip)^2 + yoff(x0,x1)^2
                     -------------------------------------------
                                 fwhm^2 * (1-ellip)^2

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

    def __init__(self, name='normgauss2d'):
        self.fwhm = Parameter(name, 'fwhm', 10, tinyval, hard_min=tinyval)
        self.xpos = Parameter(name, 'xpos', 0)
        self.ypos = Parameter(name, 'ypos', 0)
        self.ellip = Parameter(name, 'ellip', 0, 0, 0.999, 0, 0.9999,
                               frozen=True)
        self.theta = Parameter(name, 'theta', 0, -2 * numpy.pi, 2 * numpy.pi,
                               -2 * numpy.pi, 4 * numpy.pi, 'radians',
                               frozen=True)
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
        ampl = guess_amplitude2d(dep, *args)
        fwhm = guess_fwhm(dep, *args)
        param_apply_limits(xpos, self.xpos, **kwargs)
        param_apply_limits(ypos, self.ypos, **kwargs)
        param_apply_limits(fwhm, self.fwhm, **kwargs)

        # Apply normalization factor to guessed amplitude
        norm = (numpy.pi / _gfactor) * self.fwhm.val * self.fwhm.val * \
            numpy.sqrt(1.0 - (self.ellip.val * self.ellip.val))
        for key in ampl.keys():
            if ampl[key] is not None:
                ampl[key] *= norm
        param_apply_limits(ampl, self.ampl, **kwargs)

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.ngauss2d(*args, **kwargs)


class Polynom2D(RegriddableModel2D):
    """Two-dimensional polynomial function.

    The maximum order of the polynomial is 2.

    Attributes
    ----------
    c
        The constant term.
    cy1
        The coefficient for the x1 term.
    cy2
        The coefficient for the x1^2 term.
    cx1
        The coefficient for the x0 term.
    cx1y1
        The coefficient for the x0 x1 term.
    cx1y2
        The coefficient for the x0 x1^2 term.
    cx2
        The coefficient for the x0^2 term.
    cx2y1
        The coefficient for the x0^2 x1 term.
    cx2y2
        The coefficient for the x0^2 x1^2 term.

    See Also
    --------
    Polynom1D

    Notes
    -----
    The functional form of the model for points is::

        f(x,x1) = c + cx1 * x0 + cx2 * x0^2 +
                      cy1 * x1 + cy2 * x1^2 +
                      cx1y1 * x0 * x1 +
                      cx1y2 * x0 * x1^2 +
                      cx2y1 * x0^2 * x1 +
                      cx2y2 * x0^2 * x1^2

    and for an integrated grid it is the integral of this over
    the bin.
    """

    def __init__(self, name='polynom2d'):
        self.c = Parameter(name, 'c', 1)
        self.cx1 = Parameter(name, 'cx1', 0)
        self.cx2 = Parameter(name, 'cx2', 0)
        self.cy1 = Parameter(name, 'cy1', 0)
        self.cy2 = Parameter(name, 'cy2', 0)
        self.cx1y1 = Parameter(name, 'cx1y1', 0)
        self.cx1y2 = Parameter(name, 'cx1y2', 0)
        self.cx2y1 = Parameter(name, 'cx2y1', 0)
        self.cx2y2 = Parameter(name, 'cx2y2', 0)
        ArithmeticModel.__init__(self, name,
                                 (self.c, self.cy1, self.cy2, self.cx1,
                                  self.cx1y1, self.cx1y2, self.cx2,
                                  self.cx2y1, self.cx2y2))
        self.cache = 0

    def guess(self, dep, *args, **kwargs):
        x0min = args[0].min()
        x0max = args[0].max()
        x1min = args[1].min()
        x1max = args[1].max()
        ymin = dep.min()
        ymax = dep.max()

        ylo = 0
        if numpy.abs(ymin - 0) > DBL_EPSILON:
            ylo = -ymin
            if ymin < 0: ylo = ymin

        yhi = 0
        if numpy.abs(ymax - 0) > DBL_EPSILON:
            yhi = -ymax
            if ymax > 0: yhi = ymax

        dydx0 = (ymax - ymin) / (x0max - x0min)
        dydx1 = (ymax - ymin) / (x1max - x1min)
        dyd2x0 = (ymax - ymin) / ((x0max - x0min) * (x0max - x0min))
        dyd2x1 = (ymax - ymin) / ((x1max - x1min) * (x1max - x1min))
        dydx0dx1 = (ymax - ymin) / ((x0max - x0min) * (x1max - x1min))

        c     = {'val': (ymax + ymin) / 2., 'min': ylo, 'max': yhi}
        cx1   = {'val': 0., 'min': -100 * dydx0, 'max': 100 * dydx0}
        cy1   = {'val': 0., 'min': -100 * dydx1, 'max': 100 * dydx1}
        cx2   = {'val': 0., 'min': -100 * dyd2x0, 'max': 100 * dyd2x0}
        cy2   = {'val': 0., 'min': -100 * dyd2x1, 'max': 100 * dyd2x1}
        cx1y1 = {'val': 0., 'min': -100 * dydx0dx1, 'max': 100 * dydx0dx1}
        c22   = {'val': 0., 'min': ylo, 'max': yhi }

        param_apply_limits(c, self.c, **kwargs)
        param_apply_limits(cx1, self.cx1, **kwargs)
        param_apply_limits(cy1, self.cy1, **kwargs)
        param_apply_limits(cx2, self.cx2, **kwargs)
        param_apply_limits(cy2, self.cy2, **kwargs)
        param_apply_limits(cx1y1, self.cx1y1, **kwargs)
        param_apply_limits(c22, self.cx1y2, **kwargs)
        param_apply_limits(c22, self.cx2y1, **kwargs)
        param_apply_limits(c22, self.cx2y2, **kwargs)

    def calc(self, *args, **kwargs):
        kwargs['integrate'] = bool_cast(self.integrate)
        return _modelfcts.poly2d(*args, **kwargs)


class TableModel(ArithmeticModel):

    def __init__(self, name='tablemodel'):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        self.__x = None
        self.__y = None
        self.__filtered_y = None
        self.filename = None
        self.method = linear_interp  # interpolation method
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.ampl,))

    def __setstate__(self, state):
        self.__x = None
        self.__y = state.pop('_y', None)
        self.__filtered_y = state.pop('_filtered_y', None)
        self.filename = state.pop('_file', None)
        ArithmeticModel.__setstate__(self, state)

    def load(self, x, y):
        self.__y = y
        self.__x = x

        # Input grid is sorted!
        if x is not None:
            idx = numpy.asarray(x).argsort()
            self.__y = numpy.asarray(y)[idx]
            self.__x = numpy.asarray(x)[idx]

    def get_x(self):
        return self.__x

    def get_y(self):
        return self.__y

    def fold(self, data):
        mask = data.mask
        if self.__x is None and numpy.iterable(mask):
            if len(mask) != len(self.__y):
                raise ModelErr("filtermismatch", 'table model',
                               'data, (%s vs %s)' %
                               (len(self.__y), len(mask)))
            self.__filtered_y = self.__y[mask]

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):

        if self.__x is not None and self.__y is not None:
            return p[0] * interpolate(x0, self.__x, self.__y, function=self.method)

        elif (self.__filtered_y is not None and
              len(x0) == len(self.__filtered_y)):
            return p[0] * self.__filtered_y

        elif (self.__y is not None and
              len(x0) == len(self.__y)):
            return p[0] * self.__y

        raise ModelErr("filtermismatch", 'table model', 'data, (%s vs %s)' %
                       (len(self.__y), len(x0)))


class UserModel(ArithmeticModel):
    """Support for user-supplied models.

    The class is used by sherpa.ui.load_user_model and
    sherpa.astro.ui.load_user_model.
    """
    def __init__(self, name='usermodel', pars=None):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        self._y = []
        self._file = None
        if pars is None:
            self.ampl = Parameter(name, 'ampl', 1)
            pars = (self.ampl,)
        else:
            for par in pars:
                setattr(self, par.name, par)

        ArithmeticModel.__init__(self, name, pars)


class Integrator1D(CompositeModel, RegriddableModel1D):

    @staticmethod
    def wrapobj(obj):
        if isinstance(obj, ArithmeticModel):
            return obj
        return ArithmeticFunctionModel(obj)

    def __init__(self, model, *otherargs, **otherkwargs):
        self.model = self.wrapobj(model)
        self.otherargs = otherargs
        self.otherkwargs = otherkwargs
        self._errflag = 0
        CompositeModel.__init__(self,
                                ('integrate1d(%s)' % self.model.name),
                                (self.model,))

    def startup(self, cache=False):
        self.model.startup(cache)
        self._errflag = 1
        CompositeModel.startup(self, cache)

    def teardown(self):
        self.model.teardown()
        CompositeModel.teardown(self)

    def calc(self, p, xlo, xhi=None, **kwargs):
        if xhi is None:
            raise ModelErr('needsint')

        return _modelfcts.integrate1d(self.model.calc,
                                      p, xlo, xhi, **self.otherkwargs)


class Integrate1D(RegriddableModel1D):

    def __init__(self, name='integrate1d'):
        tol = numpy.finfo(float).eps
        self.epsabs = Parameter(name, 'epsabs', tol, alwaysfrozen=True)
        self.epsrel = Parameter(name, 'epsrel', 0, alwaysfrozen=True)
        self.maxeval = Parameter(name, 'maxeval', 10000, alwaysfrozen=True)
        ArithmeticModel.__init__(self, name, (self.epsabs, self.epsrel,
                                              self.maxeval))

    def __call__(self, model):
        return Integrator1D(model,
                            epsabs=self.epsabs.val,
                            epsrel=self.epsrel.val,
                            maxeval=int(self.maxeval.val),
                            logger=warning)
