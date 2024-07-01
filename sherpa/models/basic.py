#
#  Copyright (C) 2010, 2016, 2018 - 2024
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

import logging

import numpy

from sherpa.utils import bool_cast, interpolate, linear_interp, \
    sao_fcmp
from sherpa.utils.err import ModelErr
from sherpa.utils.guess import get_position, guess_amplitude, \
    guess_amplitude_at_ref, guess_amplitude2d, guess_bounds, \
    guess_fwhm, guess_position, guess_reference, param_apply_limits

from .parameter import Parameter, tinyval
from .model import ArithmeticModel, modelCacher1d, CompositeModel, \
    ArithmeticFunctionModel, RegriddableModel2D, RegriddableModel1D
from . import _modelfcts  # type: ignore

warning = logging.getLogger(__name__).warning


__all__ = ('Box1D', 'Const1D', 'Cos', 'Delta1D', 'Erf', 'Erfc', 'Exp', 'Exp10',
           'Gauss1D', 'Log', 'Log10', 'LogParabola', 'NormGauss1D', 'Poisson',
           'Polynom1D', 'PowLaw1D', 'Scale1D', 'Sin', 'Sqrt', 'StepHi1D',
           'StepLo1D', 'Tan',
           'Box2D', 'Const2D', 'Delta2D', 'Gauss2D', 'SigmaGauss2D',
           'NormGauss2D', 'Polynom2D', 'Scale2D', 'UserModel', 'TableModel',
           'Integrate1D')

DBL_EPSILON = numpy.finfo(float).eps

# This relies on the PyArg_ParseTypleAndKeywords calls in model_extension.hh
#
ALLOWED_KEYWORDS_1D = set(["pars", "xlo", "xhi", "integrate"])
ALLOWED_KEYWORDS_2D = set(["pars", "x0lo", "x1lo", "x0hi", "x1hi", "integrate"])


def clean_kwargs(allowed, model, kwargs):
    """Remove un-supported keywords for these models.

    Parameters
    ----------
    model : Model instance
        It must have an integrate field.
    kwargs : dict
        The input keyword arguments
    allowed : set of str
        The allowed keyword names.

    Returns
    -------
    allowed : dict
        The keyword arguments (names and values) that are allowed.

    """

    # TODO: remove the use of bool_cast here.
    #
    out = {"integrate": bool_cast(model.integrate)}
    for key in [k for k in kwargs.keys() if k in allowed]:
        out[key] = kwargs[key]

    return out


def clean_kwargs1d(model, kwargs):
    """Remove un-supported keywords for these models.

    The supported arguments are ALLOWED_KEYWORDS_1D.

    """

    return clean_kwargs(ALLOWED_KEYWORDS_1D, model, kwargs)


def clean_kwargs2d(model, kwargs):
    """Remove un-supported keywords for these models.

    The supported arguments are ALLOWED_KEYWORDS_2D.

    """

    return clean_kwargs(ALLOWED_KEYWORDS_2D, model, kwargs)


class Box1D(RegriddableModel1D):
    """One-dimensional box function.

    The model is flat between ``xlow`` and ``xhi`` (both limits are
    inclusive), where it is set to the ``ampl`` parameter. Outside
    this range the model is zero.

    .. versionchanged:: 4.16.0
       The default value for the xhi parameter has been changed from 0
       to 1, and the range of ampl now matches xlow and xhi rather
       than being set to -1 to 1.

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
    Box2D, Const1D, Delta1D, Scale1D, StepLo1D, StepHi1D

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
        self.xhi = Parameter(name, 'xhi', 1)
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.xlow, self.xhi, self.ampl))

    def guess(self, dep, *args, **kwargs):
        lo, hi = guess_bounds(args[0])
        norm = guess_amplitude(dep, *args)
        param_apply_limits(lo, self.xlow, **kwargs)
        param_apply_limits(hi, self.xhi, **kwargs)
        param_apply_limits(norm, self.ampl, **kwargs)

    @modelCacher1d
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.box1d(p, *args, **kwargs)


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
            if min < 0:
                ylo = 0
        yhi = 0
        if numpy.abs(max - 0) > DBL_EPSILON:
            yhi = -10 * max
            if max > 0:
                yhi = 100 * max
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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.const1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.cos(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.delta1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.erf(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.erfc(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.exp(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.exp10(p, *args, **kwargs)


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

    >>> from sherpa.models.basic import Gauss1D, NormGauss1D
    >>> m1 = Gauss1D()
    >>> m2 = NormGauss1D()
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

    >>> import numpy as np
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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.gauss1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.log(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.log10(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.logparabola(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.ngauss1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.poisson(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.poly1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        if kwargs['integrate']:
            # avoid numerical issues with C pow() function close to zero,
            # 0.0 +- ~1.e-14.  PowLaw1D integrated has multiple calls to
            # pow(X, 1.0 - gamma).  So gamma values close to 1.0 +- 1.e-10
            # should be be 1.0 to avoid errors propagating in the calculated
            # model.
            if sao_fcmp(p[0], 1.0, 1.e-10) == 0:
                # pars { gamma, ref, ampl }
                p[0] = 1.0

        return _modelfcts.powlaw(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.sin(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.sqrt(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.stephi1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.steplo1d(p, *args, **kwargs)


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
    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs1d(self, kwargs)
        return _modelfcts.tan(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.box2d(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.const2d(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.delta2d(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.gauss2d(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.sigmagauss2d(p, *args, **kwargs)


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

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.ngauss2d(p, *args, **kwargs)


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
            if ymin < 0:
                ylo = ymin

        yhi = 0
        if numpy.abs(ymax - 0) > DBL_EPSILON:
            yhi = -ymax
            if ymax > 0:
                yhi = ymax

        dydx0 = (ymax - ymin) / (x0max - x0min)
        dydx1 = (ymax - ymin) / (x1max - x1min)
        dyd2x0 = (ymax - ymin) / ((x0max - x0min) * (x0max - x0min))
        dyd2x1 = (ymax - ymin) / ((x1max - x1min) * (x1max - x1min))
        dydx0dx1 = (ymax - ymin) / ((x0max - x0min) * (x1max - x1min))

        c = {'val': (ymax + ymin) / 2., 'min': ylo, 'max': yhi}
        cx1 = {'val': 0., 'min': -100 * dydx0, 'max': 100 * dydx0}
        cy1 = {'val': 0., 'min': -100 * dydx1, 'max': 100 * dydx1}
        cx2 = {'val': 0., 'min': -100 * dyd2x0, 'max': 100 * dyd2x0}
        cy2 = {'val': 0., 'min': -100 * dyd2x1, 'max': 100 * dyd2x1}
        cx1y1 = {'val': 0., 'min': -100 * dydx0dx1, 'max': 100 * dydx0dx1}
        c22 = {'val': 0., 'min': ylo, 'max': yhi}

        param_apply_limits(c, self.c, **kwargs)
        param_apply_limits(cx1, self.cx1, **kwargs)
        param_apply_limits(cy1, self.cy1, **kwargs)
        param_apply_limits(cx2, self.cx2, **kwargs)
        param_apply_limits(cy2, self.cy2, **kwargs)
        param_apply_limits(cx1y1, self.cx1y1, **kwargs)
        param_apply_limits(c22, self.cx1y2, **kwargs)
        param_apply_limits(c22, self.cx2y1, **kwargs)
        param_apply_limits(c22, self.cx2y2, **kwargs)

    def calc(self, p, *args, **kwargs):
        kwargs = clean_kwargs2d(self, kwargs)
        return _modelfcts.poly2d(p, *args, **kwargs)


class TableModel(ArithmeticModel):
    """Tabulated values are linearly scaled and may be interpolated.

    The `load` method is used to read in the tabular data. There are
    two different modes:

      - both independent (x) and dependent (y) axes are set, in which
        case model evaluation will interpolate the requested grid onto
        the data, with the interpolation scheme controlled by the
        `method` attribute. The `fold` method does not need to be
        used.

      - the x axis is set to `None`, in which case the grid used to
        evaluate the model must have the same size as the y values, and
        the `fold` method must be sent the data object to be fit
        before a fit is made. In this case the `method` attribute is
        ignored.

    The model values - the y argument to `load` - are multiplied by
    the `ampl` parameter of the model.

    Attributes
    ----------
    ampl
        The linear scaling factor for the table values

    Notes
    -----
    The model has `ndim` set to `None` as it can be used for any
    dimension of data. However, this only makes sense for the second
    mode, that is when the x argument to `load` is `None`.  For the
    first mode, when x is not `None`, the model evaluation only uses
    the first component of the independent axes - so the low values
    for a Data1DInt object or the x0 values for a Data1D object.

    The model ignores the integrate setting.

    See Also
    --------
    Const1D, Scale1D

    Examples
    --------

    If the x array is given to the `load` call then the model can
    interpolate from the requested grid onto the load data (the
    default is linear interpolation):

    >>> tm = TableModel()
    >>> tm.load([10, 20, 25, 30], [14, 12, 17, 18])
    >>> tm.ampl = 10
    >>> tm([15, 20, 27])
    array([130., 120., 174.])

    The `method` attribute can be changed to select a different
    interpolation scheme: it requires a function that accepts (xout,
    xin, yin) and returns yout which are the interpolated values of
    xin,yin onto the grid xout:

    >>> from sherpa.utils import neville
    >>> tm.method = neville
    >>> tm([15, 20, 27])
    array([ 90.   , 120.   , 182.16])

    The model can be used for data when the independent axis is either
    not useful - such as for an image mask - or the data does not have
    a meaningful independent axis, as in the the independent variable
    being an index for a star and the dependent axis is a property of
    each star. In this case the x argument to `load` is set to `None`,
    which means that no interpolation is used and that the `fold`
    method must be used if the data has been filtered in any way, but
    it's safest to always use it:

    >>> from sherpa.models.basic import TableModel
    >>> from sherpa.data import Data1D
    >>> from sherpa.fit import Fit
    >>> d = Data1D('data', [1, 2, 3, 4, 5], [1.2, .4, 2.2, .3, 1.])
    >>> d.staterror = [.2, .2, .2, .2, .2]
    >>> tm1 = TableModel('tabmodel')
    >>> tm1.load(None, [.6, .2, 1.1, .2, .5])
    >>> tm1.ampl.val
    1.0
    >>> tm1.fold(d)
    >>> fit1 = Fit(d, tm1)
    >>> res1 = fit1.fit()
    >>> tm1.ampl.val
    1.9894736842102083

    In this case the `fold` method is necessary, to ensure that the
    fit only uses the valid data bins:

    >>> import numpy as np
    >>> y = np.ma.masked_invalid([np.nan, np.nan, 2.2, .3, 1.])
    >>> d = Data1D('data', [1, 2, 3, 4, 5], y)
    >>> d.staterror = [.2, .2, .2, .2, .2]
    >>> tm2 = TableModel('tabmodel')
    >>> tm2.load(None, [.6, .2, 1.1, .2, .5])
    >>> tm2.fold(d)
    >>> tm2.ampl.val
    1.0
    >>> fit2 = Fit(d, tm2)
    >>> res2 = fit2.fit()
    >>> tm2.ampl.val
    1.9866666666663104

    The masking also holds if the notice or ignore method has been
    used on the dataset:

    >>> d = Data1D('data', [1, 2, 3, 4, 5], [1.2, .4, 2.2, .3, 1])
    >>> d.staterror = [.2, .2, .2, .2, .2]
    >>> d.ignore(xhi=2)
    >>> tm3 = TableModel('tabmodel')
    >>> tm3.load(None, [.6, .2, 1.1, .2, .5])
    >>> tm3.fold(d)
    >>> tm3.ampl.val
    1.0
    >>> fit3 = Fit(d, tm3)
    >>> res = fit3.fit()
    >>> tm3.ampl.val
    1.9866666666663104

    """

    @property
    def method(self):
        """The interpolation method, used when x is not None in the load call.

        The method argument is a function that accepts arguments (xout,
        xin, yin) and returns the yout values from interpolating xout onto
        (xin, yin). The default is linear interpolation
        (sherpa.utils.linear_interp).
        """
        return self._method

    @method.setter
    def method(self, val):
        self._method = val

        # as the method affects the cache, clear it (we could skip
        # this if the method has not changed but is it worth it?)
        #
        self.cache_clear()

    def __init__(self, name='tablemodel'):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        #
        # NOTE added: actually, we can not have user-defined names
        # here thanks to the NoNewAttributesAfterInit base class.
        #
        self.__x = None
        self.__y = None
        self.__filtered_y = None
        self.filename = None
        self._method = linear_interp
        self.ampl = Parameter(name, 'ampl', 1)
        ArithmeticModel.__init__(self, name, (self.ampl,))

    def load(self, x, y):
        """Set the model values.

        Parameters
        ----------
        x, y : None or sequence
           The model values. It is expected that either both are
           given, and have the same number of elements, or that only y
           is set, although the model can be cleared by setting both
           to None. If x is given then the data is saved after being
           sorted into increasing order of x.

        See Also
        --------
        get_x, get_y

        """

        if x is not None and y is None:
            raise ModelErr("y must be set if x is set")

        # Clear the cache. We could avoid doing this if the
        # data has not changed but this is not worth the
        # complexity.
        #
        self.cache_clear()

        # A simplified version of sherpa.data._check, which is not
        # used here to avoid circular dependencies. Rather than raise
        # DataErr, we raise ModelErr with a similar message.
        #
        def _check(val):
            if val is None:
                return None

            val = numpy.asarray(val)
            if val.ndim != 1:
                raise ModelErr("Array must be 1D or None")

            # Check we can add 0 to the values. This is to try and
            # catch cases when a string column is passed to load. It
            # will not catch all possible problem cases. An
            # alternative would be to check
            #    isinstance(val[0], numbers.Number)
            # but it's not clear which is best.
            #
            # We could also check whether values are np.isfinite() but
            # users may want to read in bad values that we then always
            # mask out.
            #
            try:
                val + 0
            except Exception:
                raise ModelErr(f"Unable to treat array as numeric: {val}") from None

            return val

        self.__y = _check(y)
        self.__x = _check(x)

        # clear the filtered array
        self.__filtered_y = None

        if x is None:
            return

        nx = len(self.__x)
        ny = len(self.__y)
        if nx != ny:
            raise ModelErr(f"size mismatch between x and y: {nx} vs {ny}")

        # Ensure the data is sorted.
        #
        idx = self.__x.argsort()
        self.__y = self.__y[idx]
        self.__x = self.__x[idx]

    def get_x(self):
        """Return the independent axis or None.

        See Also
        --------
        get_y, load

        """
        return self.__x

    def get_y(self):
        """Return the dependent axis or None.

        See Also
        --------
        get_x, load

        """
        return self.__y

    def fold(self, data):
        """Ensure the model matches the data filter.

        This should be called after load, to ensure that any existing
        filter is applied correctly during a fit. It is only necessary
        for models where x is None in the load call.

        Parameters
        ----------
        data : sherpa.data.Data instance
           An object with a mask attribute.

        """

        if self.__y is None:
            raise ModelErr("The tablemodel's load method must be called first")

        # Clear out the setting. If needed it will get reset.
        #
        self.__filtered_y = None

        # If we are interpolating the data we do not care about the
        # data mask.
        #
        if self.__x is not None:
            return

        # What should we do with data.mask = {True, False}?
        #
        mask = data.mask
        if not numpy.iterable(mask):
            return

        # At this point we know mask is an iterable, so it should
        # match the y data of the model.
        #
        if len(mask) != len(self.__y):
            raise ModelErr("filtermismatch", 'table model',
                           f"data, ({len(self.__y)} vs {len(mask)})")

        self.__filtered_y = self.__y[mask]

    @modelCacher1d
    def calc(self, p, x0, x1=None, *args, **kwargs):
        """Evaluate the model.

        The load method must have been called first. If both x and y
        were given then the model is interpolated onto the x0 grid,
        otherwise x0 is only checked to see if it has the right size
        (fold should have been called if the data has been filtered).

        """
        if self.__y is None:
            raise ModelErr("The tablemodel's load method must be called first")

        if self.__x is not None:
            return p[0] * interpolate(x0, self.__x, self.__y, function=self.method)

        if (self.__filtered_y is not None and
              len(x0) == len(self.__filtered_y)):
            return p[0] * self.__filtered_y

        if len(x0) == len(self.__y):
            return p[0] * self.__y

        raise ModelErr("filtermismatch", 'table model',
                       f"data, ({len(self.__y)} vs {len(x0)})")


class UserModel(ArithmeticModel):
    """Support for user-supplied models.

    The class is used by sherpa.ui.load_user_model and
    sherpa.astro.ui.load_user_model.
    """
    def __init__(self, name='usermodel', pars=None):
        # these attributes should remain somewhat private
        # as not to conflict with user defined parameter names
        #
        # NOTE added: actually, we can not have user-defined names
        # here thanks to the NoNewAttributesAfterInit base class.
        #
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
    """Integrate a one-dimensional model.

    It is expected that instances of this class are created by the
    `Integrate1D` class and not directly.

    Attributes
    ----------
    model
        The model to be evaluated.
    otherargs
        Currently unused.
    otherkwargs
        Used to pass extra parameters to the integrator (currently
        `epsabs`, `epsrel`, `maxeval`, `errflag`, and `logger`).

    Raises
    ------
    sherpa.utils.err.ModelErr
        The model was evaluated on a "point" grid.

    See Also
    --------
    Integrate1D

    Examples
    --------

    Numericaly integrate the Gauss1D model across the bins defined by
    xlo and xhi, using an absolute tolerance of 1e-5. Note that you
    would not use this in practice, since the Gauss1D model already
    supports integrated bins, but it does allow us to compare the two
    approaches (the `y` and `ytrue` arrays).

    >>> import numpy as np
    >>> from sherpa.models.basic import Gauss1D, Integrator1D
    >>> gmdl = Gauss1D()
    >>> imdl = Integrator1D(gmdl, epsabs=1e-5)
    >>> xgrid = np.asarray([0, 0.5, 1, 2, 4])
    >>> xlo = xgrid[:-1]
    >>> xhi = xgrid[1:]
    >>> y = imdl(xlo, xhi)
    >>> ytrue = gmdl(xlo, xhi)
    >>> (y - ytrue) / ytrue
    array([-1.11278878e-16,  0.00000000e+00,  0.00000000e+00, -1.43150678e-16])

    """

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
                                f'integrate1d({self.model.name})',
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


# DOC note: we do not expose Integrator1D by default so it is not
# mentioned as a See Also link here.
#
class Integrate1D(RegriddableModel1D):
    """Integrate a model across each bin (one dimensional).

    Numerically integrate a one-dimensional model across each bin of a
    "histogram" dataset - that is one with low and high edges for each
    bin. The model to be integrated is supplied as an input to this
    model and any attribute changes must be made before this.

    Attributes
    ----------
    epsabs
        The maximum absolute difference allowed when integrating the
        model. This parameter is always frozen.
    epsrel
        The maximum relative difference allowed when integrating the
        model. This parameter is always frozen.
    maxeval
        The maximum number of iterations allowed per bin.  This
        parameter is always frozen.

    Notes
    -----
    If `imdl` is an instance of `Integrate1D` and `omdl` the model to
    integrate, then `imdl(omdl)` creates the integrated form. Note
    that changes to the `Integrate1D` parameters - such as `epsabs` -
    must be made before creating this integrated form. For example, to
    change the absolute tolerance to a value appropriate for 32-bit
    floats:

    >>> import numpy as np
    >>> omdl = Polynom1D()  # note: this class already supports integration
    >>> imdl = Integrate1D(name='imdl')
    >>> imdl.epsabs = np.finfo(np.float32).eps
    >>> mdl = imdl(omdl)

    Examples
    --------

    The Integrate1D model lets you use a one-dimensional model which
    can only be evaluated at a point in a case where the dataset has a
    low and high edge for the independent axis - such as a Data1DInt
    object. The Scale1D model is used as an example of a model which
    only evaluates at a point (as it always returns the scale value),
    but in actual code you should use the Const1D model rather than
    apply Integrate1D to Scale1D. The evaluation of the `smdl`
    instance returns the `c0` parameter value (set to 10) at each
    point, whereas the `mdl` instance integrates this across each bin,
    and so returning the bin width times `c0` for each bin. Note that
    the tolerance for the integration:

    - was changed from the default (tolerance for 64-bit float) to
      the 32-bit float tolerance (to avoid a warning message when
      evaluating mdl);

    - and must be changed before being applied to the model to
      integrate (smdl) in this case.

    >>> import numpy as np
    >>> from sherpa.models.basic import Scale1D, Integrate1D
    >>> xlo, xhi = [2, 5, 8], [4, 8, 12]
    >>> imdl = Integrate1D(name='imdl')
    >>> imdl.epsabs = np.finfo(np.float32).eps
    >>> smdl = Scale1D(name='smdl')
    >>> mdl = imdl(smdl)
    >>> smdl.c0 = 10
    >>> print(mdl)
    integrate1d(smdl)
       Param        Type          Value          Min          Max      Units
       -----        ----          -----          ---          ---      -----
       smdl.c0      thawed           10 -3.40282e+38  3.40282e+38
    >>> print(smdl(xlo, xhi))
    [10. 10. 10.]
    >>> print(mdl(xlo, xhi))
    [20. 30. 40.]

    """

    def __init__(self, name='integrate1d'):
        tol = DBL_EPSILON
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
