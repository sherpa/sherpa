#
#  Copyright (C) 2007, 2015, 2016, 2018 - 2024
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

"""Objects and utilities used by multiple Sherpa subpackages.

Code in this module should not be considered stable, as it may be
moved, changed, or removed.

"""

from collections.abc import Iterable
import inspect
import logging
import operator
import os
import pydoc
import string
import sys
from types import FunctionType, MethodType
from typing import Any, Callable, Generic, Optional, Sequence, \
    TypeVar
import warnings

import numpy as np

# Note: _utils.gsl_fcmp and _utils.ndtri are not exported from
#       this module; is this intentional?
from ._utils import hist1d, hist2d  # type: ignore
from . import _utils  # type: ignore
from . import _psf    # type: ignore

from .err import IOErr

# We re-export symbols from sherpa.utils modules but this will be
# removed at some point.
#
from .guess import _guess_ampl_scale, get_midpoint, get_peak, \
    get_position, get_valley, guess_amplitude, guess_amplitude2d, \
    guess_amplitude_at_ref, guess_bounds, guess_fwhm, get_fwhm, \
    guess_position, guess_radius, guess_reference, param_apply_limits
from .parallel import multi as _multi, ncpus as _ncpus, \
    parallel_map, parallel_map_funcs, run_tasks
from .random import poisson_noise


warning = logging.getLogger("sherpa").warning
debug = logging.getLogger("sherpa").debug

__all__ = ('NoNewAttributesAfterInit',
           '_guess_ampl_scale', 'apache_muller', 'bisection', 'bool_cast',
           'calc_ftest', 'calc_mlr', 'calc_total_error', 'create_expr',
           'create_expr_integrated',
           'dataspace1d', 'dataspace2d', 'demuller',
           'erf', 'export_method', 'extract_kernel',
           'filter_bins', 'gamma', 'get_fwhm',
           'get_keyword_defaults', 'get_keyword_names', 'get_midpoint',
           'get_num_args', 'get_peak', 'get_position', 'get_valley',
           'guess_amplitude', 'guess_amplitude2d', 'guess_amplitude_at_ref',
           'guess_bounds', 'guess_fwhm', 'guess_position', 'guess_radius',
           'guess_reference', 'histogram1d', 'histogram2d', 'igam', 'igamc',
           'incbet', 'interpolate', 'is_binary_file', 'Knuth_close',
           'lgam', 'linear_interp', 'nearest_interp',
           'neville', 'neville2d',
           'new_muller', 'normalize',
           'pad_bounding_box', 'parallel_map', 'parallel_map_funcs',
           'param_apply_limits', 'parse_expr', 'poisson_noise',
           'print_fields', 'rebin',
           'sao_arange', 'sao_fcmp', 'send_to_pager',
           'set_origin', 'sum_intervals', 'zeroin',
           'multinormal_pdf', 'multit_pdf', 'get_error_estimates', 'quantile',
           )


###############################################################################
#
# Types
#
###############################################################################


# This logic was found in several modules so centralize it. Note that
# this is not added to __all__.
#
def is_subclass(t1, t2):
    """Is t2 a subclass of t1 but not the same as t1?"""
    return inspect.isclass(t1) and issubclass(t1, t2) and (t1 is not t2)


###############################################################################


class NoNewAttributesAfterInit:
    """

    Prevents attribute deletion and setting of new attributes after
    __init__ has been called.  Derived classes must call
    NoNewAttributesAfterInit.__init__ after all other initialization.

    """

    __initialized = False  # Use name mangling

    def __init__(self) -> None:
        self.__initialized = True

    def __delattr__(self, name: str) -> None:
        if self.__initialized and hasattr(self, name):
            raise AttributeError(f"'{type(self).__name__}' object attribute '{name}' "
                                  "cannot be deleted")
        object.__delattr__(self, name)

    def __setattr__(self, name: str, val: Any) -> None:
        if self.__initialized and (not hasattr(self, name)):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        if self.__initialized and hasattr(self, name):
            cname = callable(getattr(self, name))
            cval = callable(val)
            if cname and not cval:
                raise AttributeError(f"'{type(self).__name__}' object attribute '{name}' "
                                     "cannot be replaced with a non-callable attribute")

            if not cname and cval:
                raise AttributeError(f"'{type(self).__name__}' object attribute '{name}' "
                                     "cannot be replaced with a callable attribute")

        object.__setattr__(self, name, val)


###############################################################################
#
# Compiled Utilities: _utils
#
###############################################################################


def calc_ftest(dof1, stat1, dof2, stat2):
    """Compare two models using the F test.

    The F-test is a model comparison test; that is, it is a test
    used to select from two competing models which best describes
    a particular data set. A model comparison test statistic, T,
    is created from the best-fit statistics of each fit; as with all
    statistics, it is sampled from a probability distribution p(T).
    The test significance is defined as the integral of p(T) from the
    observed value of T to infinity. The significance quantifies the
    probability that one would select the more complex model when in
    fact the null hypothesis is correct. See also `calc_mlr`.

    Parameters
    ----------
    dof1 : int or sequence of int
       degrees of freedom of the simple model
    stat1 : number or sequence of number
       best-fit chi-square statistic value of the simple model
    dof2 : int or sequence of int
       degrees of freedom of the complex model
    stat2 : number or sequence of number
       best-fit chi-square statistic value of the complex model

    Returns
    -------
    sig : number or ndarray
       The significance, or p-value. A standard threshold for
       selecting the more complex model is significance < 0.05 (the
       '95% criterion' of statistics).

    See Also
    --------
    calc_mlr, incbet

    Notes
    -----
    The F test uses the ratio of the reduced chi2, which follows
    the F-distribution, (stat1/dof1) / (stat2/dof2). The incomplete
    Beta function is used to calculate the integral of the tail of
    the F-distribution.

    The F test should only be used when:

     - the simpler of the two models is nested within the other;
       that is, one can obtain the simpler model by setting the extra
       parameters of the more complex model (often to zero or one);
     - the extra parameters have values sampled from normal
       distributions under the null hypothesis (i.e., if one samples
       many datasets given the null hypothesis and fits these data with
       the more complex model, the distributions of values for the
       extra parameters must be Gaussian);
     - those normal distributions are not truncated by parameter space
       boundaries;
     - the best-fit statistics are sampled from the chi-square
       distribution.

    See Protassov et al. 2002 [1]_ for more discussion.

    References
    ----------

    .. [1] Protassov et al., Statistics, Handle with Care: Detecting
           Multiple Model Components with the Likelihood Ratio Test,
           Astrophysical Journal, vol 571, pages 545-559, 2002,
           http://adsabs.harvard.edu/abs/2002ApJ...571..545P

    Examples
    --------
    >>> calc_ftest(11, 16.3, 10, 10.2)
    0.03452352914891555

    >>> calc_ftest([11, 11], [16.3, 16.3], [10, 9], [10.2, 10.5])
    array([0.03452353, 0.13819987])
    """

    return _utils.calc_ftest(dof1, stat1, dof2, stat2)


def calc_mlr(delta_dof, delta_stat):
    """Compare two models using the Maximum Likelihood Ratio test.

    The Maximum Likelihood Ratio (MLR) test is a model comparison
    test; that is, it is a test used to select from two competing
    models which best describes a particular data set. A model
    comparison test statistic, T, is created from the best-fit
    statistics of each fit; as with all statistics, it is sampled
    from a probability distribution p(T). The test significance is
    defined as the integral of p(T) from the observed value of T to
    infinity. The significance quantifies the probability that one
    would select the more complex model when in fact the null hypothesis
    is correct. See also `calc_ftest`.

    Parameters
    ----------
    delta_dof : int or sequence of int
       change in the number of degrees of freedom
    delta_stat : number or sequence of number
       change in the best-fit statistic value

    Returns
    -------
    sig : number or ndarray
       The significance, or p-value. A standard threshold for
       selecting the more complex model is significance < 0.05 (the
       '95% criterion' of statistics).

    See Also
    --------
    calc_ftest

    Notes
    -----
    The MLR test should only be used when:

     - the simpler of the two models is nested within the other;
       that is, one can obtain the simpler model by setting the extra
       parameters of the more complex model (often to zero or one);
     - the extra parameters have values sampled from normal
       distributions under the null hypothesis (i.e., if one samples
       many datasets given the null hypothesis and fits these data with
       the more complex model, the distributions of values for the
       extra parameters must be Gaussian);
     - those normal distributions are not truncated by parameter space
       boundaries;
     - the best-fit statistics for each fit are sampled from the
       chi-square distribution.

    See Protassov et al. 2002 [1]_ for more discussion.

    References
    ----------

    .. [1] Protassov et al., Statistics, Handle with Care: Detecting
           Multiple Model Components with the Likelihood Ratio Test,
           Astrophysical Journal, vol 571, pages 545-559, 2002,
           http://adsabs.harvard.edu/abs/2002ApJ...571..545P

    Examples
    --------

    In this example, the more-complex model has 2 extra degrees of
    freedom and a statistic value that is larger by 3.7. The MLR test
    does not provide any evidence that the complex model is a better
    fit to the data than the simple model since the result is much
    larger than 0.

    >>> calc_mlr(2, 3.7)
    0.15723716631362761
    """

    return _utils.calc_mlr(delta_dof, delta_stat)


def erf(x):
    """Calculate the error function.

    Parameters
    ----------
    x : scalar or array

    Returns
    -------
    val : scalar or array
       The error function of the input.

    See Also
    --------
    gamma

    Examples
    --------

    >>> erf(0)
    0.0

    >>> erf([1.0, 2.3])
    array([ 0.84270079,  0.99885682])
    """

    return _utils.erf(x)


def igamc(a, x):
    """Calculate the complement of the regularized incomplete Gamma function (upper).

    The function is defined using the regularized incomplete Gamma
    function - igam(a,x) - and the Gamma function - gamma(a) - as::

       igamc(a,x) = 1 - igam(a,x)
                  = 1 / gamma(a) Int_x^Inf e^(-t) t^(a-1) dt

    Parameters
    ----------
    a : scalar or array
       a > 0
    x : scalar or array
       x > 0

    Returns
    -------
    val : scalar or array
       The incomplete Gamma function of the input.

    See Also
    --------
    gamma, igam

    Notes
    -----
    In this implementation, which is provided by the Cephes Math
    Library [1]_, both arguments must be positive. The integral is
    evaluated by either a power series or continued fraction expansion,
    depending on the relative values of a and x. Using IEEE arithmetic,
    the relative errors are

    ========  ======  ========  =======  =======
     domain   domain  # trials   peak      rms
    ========  ======  ========  =======  =======
    0.5,100   0,100   200000    1.9e-14  1.7e-15
    0.01,0.5  0,100   200000    1.4e-13  1.6e-15
    ========  ======  ========  =======  =======

    References
    ----------

    .. [1] Cephes Math Library Release 2.0:  April, 1987.
           Copyright 1985, 1987 by Stephen L. Moshier.
           Direct inquiries to 30 Frost Street, Cambridge, MA 02140.

    Examples
    --------

    >>> igamc(1, 2)
    0.1353352832366127

    >>> igamc([1,1], [2,3])
    array([ 0.13533528,  0.04978707])
    """

    return _utils.igamc(a, x)


def igam(a, x):
    """Calculate the regularized incomplete Gamma function (lower).

    The function is defined using the complete Gamma function -
    gamma(a) - as::

       igam(a,x) = 1 / gamma(a) Int_0^x e^(-t) t^(a^-1) dt

    Parameters
    ----------
    a : scalar or array
       a > 0
    x : scalar or array
       x > 0

    Returns
    -------
    val : scalar or array
       The incomplete Gamma function of the input.

    See Also
    --------
    gamma, igamc

    Notes
    -----
    In this implementation, which is provided by the Cephes Math
    Library [1]_, both arguments must be positive. The integral is
    evaluated by either a power series or continued fraction expansion,
    depending on the relative values of a and x. Using IEEE arithmetic,
    the relative errors are

    ======  ========  =======  =======
    domain  # trials   peak      rms
    ======  ========  =======  =======
    0,30    200000    3.6e-14  2.9e-15
    0,100   300000    9.9e-14  1.5e-14
    ======  ========  =======  =======

    References
    ----------

    .. [1] Cephes Math Library Release 2.0:  April, 1987.
           Copyright 1985, 1987 by Stephen L. Moshier.
           Direct inquiries to 30 Frost Street, Cambridge, MA 02140.

    Examples
    --------

    >>> igam(1, 2)
    0.8646647167633873

    >>> igam([1,1], [2,3])
    array([ 0.86466472,  0.95021293])
    """

    return _utils.igam(a, x)


def incbet(a, b, x):
    """Calculate the incomplete Beta function.

    The function is defined as::

       sqrt(a+b)/(sqrt(a) sqrt(b)) Int_0^x t^(a-1) (1-t)^(b-1) dt

    and the integral from x to 1 can be obtained using the relation::

       1 - incbet(a, b, x) = incbet(b, a, 1-x)

    Parameters
    ----------
    a : scalar or array
       a > 0
    b : scalar or array
       b > 0
    x : scalar or array
       0 <= x <= 1

    Returns
    -------
    val : scalar or array
       The incomplete beta function calculated from the inputs.

    See Also
    --------
    calc_ftest

    Notes
    -----
    In this implementation, which is provided by the Cephes Math
    Library [1]_, the integral is evaluated by a continued fraction
    expansion or, when b*x is small, by a power series.

    Using IEEE arithmetic, the relative errors are (tested uniformly
    distributed random points (a,b,x) with a and b in 'domain' and
    x between 0 and 1):

    ========  ========  =======  =======
     domain   # trials   peak      rms
    ========  ========  =======  =======
    0,5       10000     6.9e-15  4.5e-16
    0,85      250000    2.2e-13  1.7e-14
    0,1000    30000     5.3e-12  6.3e-13
    0,1000    250000    9.3e-11  7.1e-12
    0,100000  10000     8.7e-10  4.8e-11
    ========  ========  =======  =======

    Outputs smaller than the IEEE gradual underflow threshold were
    excluded from these statistics.

    References
    ----------

    .. [1] Cephes Math Library Release 2.0:  April, 1987.
           Copyright 1985, 1987 by Stephen L. Moshier.
           Direct inquiries to 30 Frost Street, Cambridge, MA 02140.

    Examples
    --------

    >>> incbet(0.3, 0.6, 0.5)
    0.68786273145845922

    >>> incbet([0.3,0.3], [0.6,0.7], [0.5,0.4])
    array([ 0.68786273,  0.67356524])
    """

    return _utils.incbet(a, b, x)


def gamma(z):
    """Calculate the Gamma function.

    Parameters
    ----------
    z : scalar or array
       -171 <= z <- 171.6

    Returns
    -------
    val : scalar or array
       The gamma function of the input.

    See Also
    --------
    igam, lgam

    Notes
    -----
    This implementation is provided by the Cephes Math Library [1]_.
    Arguments ``|x| >= 34`` are reduced by recurrence and the function
    approximated by a rational function of degree 6/7 in the interval
    (2,3). Large arguments are handled by Stirling's formula. Large
    negative arguments are made positive using a reflection formula.

    Relative errors are

    ========  ========  =======  =======
     domain   # trials   peak      rms
    ========  ========  =======  =======
    -170,33   20000     2.3e-15  3.3e-16
    -33,33    20000     9.4e-16  2.2e-16
    33,171.6  20000     2.3e-15  3.2e-16
    ========  ========  =======  =======

    Errors for arguments outside the test range will be larger owing
    to amplification by the exponential function.

    References
    ----------

    .. [1] Cephes Math Library Release 2.0:  April, 1987.
           Copyright 1985, 1987 by Stephen L. Moshier.
           Direct inquiries to 30 Frost Street, Cambridge, MA 02140.

    Examples
    --------

    >>> gamma(2.3)
    1.1667119051981603

    >>> gamma([2.3,1.9])
    array([ 1.16671191,  0.96176583])
    """

    return _utils.gamma(z)


def lgam(z):
    """Calculate the log (base e) of the Gamma function.

    Parameters
    ----------
    z : scalar or array
       0 <= z <= 2.556348e305

    Returns
    -------
    val : scalar or array
       The log of the Gamma function of the input.

    See Also
    --------
    gamma, igam

    Notes
    -----
    This implementation is provided by the Cephes Math Library [1]_.
    For arguments greater than 13, the logarithm of the Gamma function
    is approximated by the logarithmic version of Stirling's formula
    using a polynomial approximation of degree 4. Arguments
    between -33 and +33 are reduced by recurrence to the interval [2,3]
    of a rational approximation. The cosecant reflection formula is
    employed for arguments less than -33.

    Relative errors are

    ===============  ========  =======  =======
        domain       # trials   peak      rms
    ===============  ========  =======  =======
    0,3              28000     5.4e-16  1.1e-16
    2.718,2.556e305  40000     3.5e-16  8.3e-17
    ===============  ========  =======  =======

    The error criterion was relative when the function magnitude was
    greater than one but absolute when it was less than one.

    The following test used the relative error criterion, though at
    certain points the relative error could be much higher than
    indicated.

    =======  ========  =======  =======
    domain   # trials   peak      rms
    =======  ========  =======  =======
    -200,-4  10000     4.8e-16  1.3e-16
    =======  ========  =======  =======

    References
    ----------

    .. [1] Cephes Math Library Release 2.0:  April, 1987.
           Copyright 1985, 1987 by Stephen L. Moshier.
           Direct inquiries to 30 Frost Street, Cambridge, MA 02140.

    Examples
    --------

    >>> lgam(104.56)
    380.21387239435785

    >>> lgam([104.56,2823.4])
    array([   380.21387239,  19607.42734396])
    """

    return _utils.lgam(z)


def sao_arange(start, stop, step=None):
    """Create a range of values between start and stop.

    See also `numpy.arange` and `numpy.linspace`.

    Parameters
    ----------
    start, stop : float
       The start and stop points.
    step : float or None, optional
       If not given the step size defaults to 1.0.

    Returns
    -------
    vals : NumPy array
       The values start, start + step, ... The last point
       is the first position where start + n * step >= stop,
       which means that it can include a point > stop.

    Examples
    --------

    >>> sao_arange(1, 3)
    array([ 1.,  2.,  3.])

    >>> sao_arange(1, 3, 0.6)
    array([ 1. ,  1.6,  2.2,  2.8,  3.4])

    """

    if step is None:
        return _utils.sao_arange(start, stop)

    return _utils.sao_arange(start, stop, step)


def sao_fcmp(x, y, tol):
    """Compare y to x, using an absolute tolerance.

    Parameters
    ----------
    x : number or array_like
       The expected value, or values.
    y : number or array_like
       The value, or values, to check. If x is an array, then
       y must be an array of the same size. If x is a scalar
       then y can be a scalar or an array.
    tol : number
       The absolute tolerance used for comparison.

    Returns
    -------
    flags : int or array_like
       0, 1, or -1 for each value in second. If the values match, then 0,
       otherwise -1 if the expected value (x) is less than the comparison
       value (y) or +1 if x is larger than y.

    See Also
    --------
    Knuth_close

    Examples
    --------

    >>> sao_fcmp(1, 1.01, 0.01)
    0

    >>> sao_fcmp(1, [0.9, 1, 1.1], 0.01)
    array([ 1,  0, -1], dtype=int32)

    >>> sao_fcmp([1.2, 2.3], [1.22, 2.29], 0.01)
    array([-1,  0], dtype=int32)
    """

    return _utils.sao_fcmp(x, y, tol)


def sum_intervals(src, indx0, indx1):
    """Sum up data within one or more pairs of indexes.

    Parameters
    ----------
    src : sequence of floats
       The data to be summed.
    indx0, indx1 : scalar or sequence of int
       The pair of indexes over which to sum the src array.
       The sizes of indx0 and indx1 must match, and each element of
       indx1 must be at least as large as the corresponding element
       in indx0.

    Returns
    -------
    val : scalar or array
       The sum of the src over the given interval ranges.

    Notes
    -----
    It is assumed that all indexes are valid. That is, they are in
    the range [0, length of src). This condition is not checked for.

    Examples
    --------

    >>> sum_intervals([1.1, 2.2, 3.3, 4.4], 1, 2)
    5.5

    >>> sum_intervals([1.1, -2.2, 3.3, 4.4], [1, 0], [3, 0])
    array([ 5.5,  1.1])

    """

    return _utils.sum_intervals(src, indx0, indx1)


def rebin(y0, x0lo, x0hi, x1lo, x1hi):
    """Rebin a histogram.

    Parameters
    ----------
    y0 : sequence of numbers
       The Y values of the histogram to rebin.
    x0lo, x0hi : sequence of numbers
       The lower and upper edges of the X values to rebin. They must match
       the size of `y0`.
    x1lo, x1hi : sequence of numbers
       The lower and upper edges of the X values of the output histogram.

    Returns
    -------
    yout : NumPy array of numbers
       The re-binned Y values (same size as `x1lo`).
    """

    return _utils.rebin(y0, x0lo, x0hi, x1lo, x1hi)


def neville(xout, xin, yin):
    """Polynomial one-dimensional interpolation using Neville's method.

    The scheme used for interpolation (Neville's method) is described
    at [1]_.

    Parameters
    ----------
    xout : array_like
       The positions at which to interpolate.
    xin : array_like
       The x values of the data to be interpolated. This must be
       sorted so that it is monotonically increasing.
    yin : array_like
       The y values of the data to interpolate (must be the same
       size as ``xin``).

    Returns
    -------
    yout : NumPy array of numbers
       The interpolated y values (same size as ``xout``).

    See Also
    --------
    interpolate, linear_interp, nearest_interp

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Neville%27s_algorithm

    Examples
    --------

    >>> import numpy as np
    >>> x = [1.2, 3.4, 4.5, 5.2]
    >>> y = [12.2, 14.4, 16.8, 15.5]
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = neville(xgrid, x, y)
    """

    return _utils.neville(xout, xin, yin)


###############################################################################
#
# Compiled Utilities: _psf
#
###############################################################################

def extract_kernel(kernel, dims_kern, dims_new, center, xlo, xhi, widths,
                   radial):
    """Extract the kernel.

    Parameters
    ----------
    kernel
    dims_kern
    dims_new
    center
    xlo
    xhi
    widths
    radial : int
       Set to 1 if using a radial profile, 0 otherwise.

    Returns
    -------
    out, dims, frac, lo, hi

    """

    return _psf.extract_kernel(kernel, dims_kern, dims_new, center,
                               xlo, xhi, widths, radial)


def normalize(xs):
    """Normalize an array.

    Parameters
    ----------
    xs : sequence
       The values to normalize. This must be a 1D array.

    Returns
    -------
    ns : ndarray
       The values of xs / sum of xs.

    """

    return _psf.normalize(xs)


def set_origin(dims, maxindex=None):
    """Return the position of the origin of the kernel.

    Parameters
    ----------
    dims : number or sequence
       The dimensions of the kernel. This should be a scalar or
       a one- or two-element sequence.
    maxindex : None or int, optional
       If given, then use this location - which is the index
       into the flattened array - as the center, otherwise
       use the center of the grid.

    Returns
    -------
    cs : ndarray or number
       The coordinates of the center, matching the input
       dims format.

    Examples
    --------

    >>> set_origin(12)
    5

    >>> set_origin([12])
    array([5])

    >>> set_origin([12], 4)
    array([4])

    >>> set_origin([12, 13])
    array([5, 6])

    >>> set_origin([12, 13], 42)
    array([6, 3])

    """

    if maxindex is None:
        return _psf.set_origin(dims)

    return _psf.set_origin(dims, maxindex)


def pad_bounding_box(kernel, mask):
    """Expand the kernel to match the mask.

    Parameters
    ----------
    kernel : numeric sequence
       The data to copy. The data is a 1D array.
    mask : int sequence
       The mask determines the size of the output and where to place
       the kernel values. It is expected that the number of non-zero
       mask elements matches the size of the `kernel` parameter.

    Returns
    -------
    nkernel : ndarray
       The output is the same size as the mask, and initialized to
       zero everywhere. Cells where the mask is non-zero are
       copied from the kernel.

    Examples
    --------

    >>> pad_bounding_box([1, 2, 3, 4], [1, 1, 0, 1, 1, 0, 0, 0, 0])
    array([ 1.,  2.,  0.,  3.,  4.,  0.,  0.,  0.,  0.])

    """

    return _psf.pad_bounding_box(kernel, mask)


###############################################################################
#
# Utilities
#
###############################################################################


# at what precisions do we assume equality in energy grids?
eps = np.finfo(np.float32).eps


def filter_bins(mins: Sequence[Optional[float]],
                maxes: Sequence[Optional[float]],
                axislist: Sequence[Sequence[float]],
                integrated: bool = False
                ) -> Optional[np.ndarray]:
    """What mask represents the given set of filters?

    The ranges are treated as inclusive at both ends if integrated is
    False, the default, otherwise the lower limit is inclusive but the
    upper limit is exclusive.

    Parameters
    ----------
    mins : sequence of values
       The minimum value of the valid range (elements may be None).
    maxes : sequence of values
       The maximum value of the valid range (elements may be None).
    axislist: sequence of arrays
       The axis to apply the range to. There must be the same
       number of elements in mins, maxes, and axislist.
       The number of elements of each element of axislist must
       also agree (the cell values do not need to match).
    integrated : bool, optional
       Is the data integrated (we have low and high bin edges)?  The
       default is False. When True it is expected that axislist
       contains a even number of rows, where the odd values are the
       low edges and the even values the upper edges, and that the
       mins and maxes only ever contain a single value, given in
       (None, hi) and (lo, None) ordering.

    Returns
    -------
    mask : ndarray or None
       A mask indicating whether the values are included (True) or
       excluded (False). If any of the input sequences are empty then
       None will be returned.

    Examples
    --------

    Calculate those points in xs which are in the range 1.5 <= x <= 4.

    >>> xs = [1, 2, 3, 4, 5]
    >>> filter_bins([1.5], [4], [xs])
    array([False,  True,  True,  True, False])

    Repeat the above calculation by combining filters for x >= 1.5
    and x <= 4 (note that the grid must be repeated for each
    filter):

    >>> filter_bins([1.5, None], [None, 4], [xs, xs])
    array([False,  True,  True,  True, False])

    For integrated data sets the lower and upper edges should be sent
    separately with the max and min limits, along with setting the
    integrated flag. The following selects the bins that cover the
    range 2 to 4 and 1.5 to 3.5:

    >>> xlo = [1, 2, 3, 4, 5]
    >>> xhi = [2, 3, 4, 5, 6]
    >>> filter_bins([None, 2], [4, None], [xlo, xhi], integrated=True)
    array([False,  True,  True,  False, False])
    >>> filter_bins([None, 1.5], [3.5, None], [xlo, xhi], integrated=True)
    array([ True,  True,  True, False, False])

    """

    mask = None

    def locheck(lo, axis):
        if integrated:
            return sao_fcmp(lo, axis, eps) < 0

        return sao_fcmp(lo, axis, eps) <= 0

    def hicheck(hi, axis):
        if integrated:
            return sao_fcmp(hi, axis, eps) > 0

        return sao_fcmp(hi, axis, eps) >= 0

    for lo, hi, axis in zip(mins, maxes, axislist):

        if (lo is None) and (hi is None):
            continue

        axis = np.asarray(axis)
        axismask = np.ones(axis.size, dtype=bool)

        if lo is not None:
            axismask &= locheck(lo, axis)

        if hi is not None:
            axismask &= hicheck(hi, axis)

        if mask is None:
            mask = axismask
        else:
            mask &= axismask

    return mask


def bool_cast(val):
    """Convert a string to a boolean.

    Parameters
    ----------
    val : bool, str or sequence
       The input value to decode.

    Returns
    -------
    flag : bool or ndarray
       True or False if val is considered to be a true or false term.
       If val is a sequence then the return value is an ndarray of
       the same size.

    Notes
    -----
    The string is compared in a case-insensitive manner to the
    following: 'true', 'on', 'yes', '1', 't', and 'y' for
    `True` and 'false', 'off', 'no', '0', 'f', and 'n' for `False`.

    If there is no match to the above then the default conversion
    provided by the `bool` routine is used.

    """

    if type(val) in (tuple, list, np.ndarray):
        return np.asarray([bool_cast(item) for item in val], bool)

    if type(val) == str:
        # since built in bool() only returns false for empty strings
        vlo = val.lower()
        if vlo in ('false', 'off', 'no', '0', 'f', 'n'):
            return False

        if vlo in ('true', 'on', 'yes', '1', 't', 'y'):
            return True

        raise TypeError(f"unknown boolean value: '{val}'")

    # use built in bool cast
    return bool(val)


# Telling the type system that the signature of meth is "the same" as
# the signature of the return value probably needs Python 3.10 for
# ParamSpec and then maybe Python 3.12 for the generic support.
#
def export_method(meth: Callable,
                  name: Optional[str] = None,
                  modname: Optional[str] = None
                  ) -> Callable:
    """
    Given a bound instance method, return a simple function that wraps
    it.  The only difference between the interface of the original
    method and the generated function is that the latter doesn't
    include 'self' in its argument list.  This means that when the
    wrapper function is called with an incorrect number of arguments,
    the error message does not include 'self' in argument counts.  The
    only reason to generate such a wrapper is to hide from a user the
    fact that they're using an instance method rather than a simple
    function.  If meth is not an instance method, it is returned
    unchanged.

    If name is None, the generated function will have the same name as
    the input method.  Otherwise, name must be a string containing the
    desired name of the new method.  If modname is not None, it must
    be a string and will be used as the module name for the generated
    function.  Note that the caller is responsible for assigning the
    returned function to an appropriate name in the calling scope.

    """

    if not isinstance(meth, MethodType):
        return meth

    # Most of the functionality here can be provided by
    # functools.wraps (and in fact would add more functionality, such
    # as the ability to handle keyword-only or positional-only
    # arguments). It would also likely require less maintenance.
    # However, it does not handle two important issues:
    #
    # a) when an error is raised it is reported as from Session.<name>
    #    rather than <name>
    #
    #    Attempts to change the __qualname__ field have not been
    #    successful. The wrapped function can include a check for an
    #    exception, manually removing any leading "Session."  text
    #    from the message, but it is not particularly good code.
    #
    # b) Error messages related to the number of arguments include the
    #    self argument (i.e. are 1 more than the user is told to
    #    expect), which is one of the main reasons for this routine.
    #

    # The only time name is not None appears to be in the tests, so
    # can this feature be removed?
    #
    if name is None:
        name = meth.__name__

    if name == meth.__name__:
        old_name = f'_old_{name}'
    else:
        old_name = meth.__name__

    defaults = meth.__defaults__
    doc = meth.__doc__

    def tostr(p):
        if p.kind == p.VAR_KEYWORD:
            return f"**{p.name}"

        if p.kind == p.VAR_POSITIONAL:
            return f"*{p.name}"

        return p.name

    # Ideally this would also identify when to add "/" or "*"
    # to indicate positional-only or keyword-only arguments.
    #
    sig = inspect.signature(meth)
    argspec = ",".join([tostr(p) for p in sig.parameters.values()])

    # Create a wrapper function with no default arguments
    g: dict[str, Any] = {old_name: meth}

    # The only time modname is None appears to be the test so can we
    # make it a required argument?
    #
    if modname is not None:
        g['__name__'] = modname

    fdef = f'def {name}({argspec}):  return {old_name}({argspec})'
    exec(fdef, g)

    # Create another new function from the one we just made, this time
    # adding the default arguments, doc string, and any annotations
    # from the original method.
    #
    # Why does this not change the __defaults__ field of new_meth
    # rather than creating a copy of it?
    #
    new_meth = g[name]
    new_meth = FunctionType(new_meth.__code__, new_meth.__globals__,
                            new_meth.__name__, defaults,
                            new_meth.__closure__)
    new_meth.__doc__ = doc
    new_meth.__annotations__ = meth.__annotations__

    return new_meth


def get_keyword_names(func, skip=0):
    """Return the names of the keyword arguments.

    Parameters
    ----------
    func
        The function to query.
    skip : int, optional
        The number of keyword arguments to skip.

    Returns
    -------
    names : list of str
        The names of the keyword arguments. It can be empty.

    See Also
    --------
    get_keyword_defaults, get_num_args

    """

    # This used to use getargspec but was changed to use inspect
    # since the former was removed briefly (circa Python 3.6).
    #
    sig = inspect.signature(func)
    kwargs = [p.name
              for p in sig.parameters.values()
              if p.kind == p.POSITIONAL_OR_KEYWORD and
              p.default != p.empty]

    return kwargs[skip:]


def get_keyword_defaults(func, skip=0):
    """Return the keyword arguments and their default values.

    Note that a similar function `sherpa.plot.backend_utils.get_keyword_defaults`
    exits, which differs from this one in that it deals with keyword-only
    arguments, not all arguments.

    Parameters
    ----------
    func
        The function to query.
    skip : int, optional
        The number of keyword arguments to skip.

    Returns
    -------
    vals : dict
        The keys are names of the keyword arguments, the values are
        the default value for that parameter. It can be empty.

    See Also
    --------
    get_keyword_names, get_num_args,
    `sherpa.plot.backend_utils.get_keyword_defaults`

    """

    # This used to use getargspec but was changed to use inspect
    # since the former was removed briefly (circa Python 3.6).
    #
    sig = inspect.signature(func)
    kwargs = [(p.name, p.default)
              for p in sig.parameters.values()
              if p.kind == p.POSITIONAL_OR_KEYWORD and
              p.default != p.empty]

    return dict(kwargs[skip:])


def get_num_args(func):
    """Return the number of arguments for a function.

    Parameters
    ----------
    func
        The function to query.

    Returns
    -------
    ntotal, npos, nkeyword : int, int, int
        The total number of arguments, the number of positional
        arguments, and the number of keyword arguments.

    See Also
    --------
    get_keyword_defaults, get_keyword_names

    """

    # This used to use getargspec but was changed to use inspect
    # since the former was removed briefly (circa Python 3.6).
    #
    sig = inspect.signature(func)
    posargs = [True
               for p in sig.parameters.values()
               if p.kind == p.POSITIONAL_OR_KEYWORD and
               p.default == p.empty]
    kwargs = [True
              for p in sig.parameters.values()
              if p.kind == p.POSITIONAL_OR_KEYWORD and
              p.default != p.empty]

    npos = len(posargs)
    nkw = len(kwargs)
    return (npos + nkw, npos, nkw)


def print_fields(names, vals, converters=None):
    """

    Given a list of strings names and mapping vals, where names is a
    subset of vals.keys(), return a listing of name/value pairs
    printed one per line in the format '<name> = <value>'.  If a value
    is a NumPy array, print it in the format
    '<data type name>[<array size>]'.  Otherwise, use str(value).

    """

    # This is the part of the deprecated typeNA dictionary Sherpa
    # would use up to v4.11.0. We included the dictionaty verbatim,
    # excluding the complex mapping which where wrong in typeNA.
    # Note only the class -> string mappings have been copied over.
    if converters is None:
        converters = {np.bool_: 'Bool',
                      np.bytes_: 'Bytes0',
                      np.complex128: 'Complex128',
                      np.complex64: 'Complex64',
                      np.datetime64: 'Datetime64',
                      np.float16: 'Float16',
                      np.float32: 'Float32',
                      np.float64: 'Float64',
                      np.int16: 'Int16',
                      np.int32: 'Int32',
                      np.int64: 'Int64',
                      np.int8: 'Int8',
                      np.object_: 'Object0',
                      np.str_: 'Str0',
                      np.timedelta64: 'Timedelta64',
                      np.uint16: 'UInt16',
                      np.uint32: 'UInt32',
                      np.uint64: 'UInt64',
                      np.uint8: 'UInt8',
                      np.void: 'Void0'
                      }
        try:
            converters[np.complex256] = 'Complex256'
        except AttributeError:
            pass
        try:
            converters[np.float128] = 'Float128'
        except AttributeError:
            pass

    width = max(len(n) for n in names)
    fmt = f'%-{width}s = %s'
    lines = []
    for n in names:
        v = vals[n]

        if isinstance(v, np.ndarray):
            v = f'{converters[v.dtype.type]}[{v.size}]'
        else:
            v = str(v)
        lines.append(fmt % (n, v))
    return '\n'.join(lines)


def create_expr(vals, mask=None, format='%s', delim='-'):
    """Create a string representation of a filter.

    Use the mask to convert the input values into a set of
    comma-separated filters - low value and high value, separated
    by the delimiter - that represent the data. If the mask is
    not given then the values must be "channel" values (that is,
    two values are consecutive if there difference is 1).

    Parameters
    ----------
    vals : sequence
        The values that represent the sequence if mask is not None,
        otherwise the selected channel numbers (in this case integer
        values).
    mask : sequence of bool or None, optional
        The mask setting for the full dataset, without any filtering
        applied. A value of True indicates the element is included
        and False means it is excluded.
    format : str, optional
        The format used to display each value.
    delim : str, optional
        The separator for a range.

    Raises
    ------
    ValueError
        If the ``vals`` and ``mask`` sequences do not match: the
        length of ``vals`` must equal the number of True values in
        ``mask``.

    See Also
    --------
    create_expr_integrated, parse_expr

    Examples
    --------

    >>> create_expr([1, 2, 3, 4])
    '1-4'

    >>> create_expr([1, 2, 4, 5, 7])
    '1-2,4-5,7'

    >>> create_expr([1, 2, 3, 4], [True, True, True, True])
    '1-4'

    >>> create_expr([0.1, 0.2, 0.4, 0.8], [True, True, True, True])
    '0.1-0.8'

    >>> create_expr([0.1, 0.2, 0.4, 0.8], [True, True, True, False, False, True])
    '0.1-0.4,0.8'

    """

    if len(vals) == 0:
        return ''

    if len(vals) == 1:
        return format % vals[0]

    if mask is None:
        seq = vals

    else:
        # Ensure we have a boolean array to make indexing behave sensibly
        # (NumPy 1.17 or so changed behavior related to this).
        #
        mask = np.asarray(mask, dtype=bool)

        # Ensure that the vals and mask array match: the number of
        # mask=True elements should equal the number of input values.
        #
        if sum(mask) != len(vals):
            raise ValueError("mask array mismatch with vals")

        # We only care about the difference between two consecutive
        # values, so it doesn't matter if index starts at 0 or 1.
        #
        index = np.arange(len(mask))
        seq = index[mask]

    exprs = []
    start = vals[0]

    # We follow create_expr_integrated but instead of having separate lo/hi
    # always we use the same array
    #
    startbins = vals[1:]
    endbins = vals[:-1]

    diffs = np.diff(seq)
    idxs, = np.where(diffs != 1)
    for idx in idxs:
        exprs.append((start, endbins[idx]))
        start = startbins[idx]

    exprs.append((start, vals[-1]))

    def filt(lo, hi):
        vstr = format % lo
        if lo == hi:
            return vstr

        return vstr + f"{delim}{format % hi}"

    return ",".join([filt(*expr) for expr in exprs])


def create_expr_integrated(lovals, hivals, mask=None,
                           format='%s', delim='-',
                           eps=np.finfo(np.float32).eps):
    """Create a string representation of a filter (integrated).

    Use the mask to convert the input values into a set of
    comma-separated filters - low value and high value, separated by
    the delimiter - that represent the data. Unlike `create_expr` this
    routine uses the lovals values for the start of the bin and
    hivals for the end of each bin, and assumes that contiguous bins
    should be combined.

    Parameters
    ----------
    lovals, hivals : sequence
        The lower and upper values of each bin. It is required that
        they are in ascending order and ``lovals`` < ``hivals``.
    mask : sequence of bool or None, optional
        The mask setting for the full dataset, without any filtering
        applied. A value of True indicates the element is included
        and False means it is excluded. Note that this is opposite to the
        numpy convention in numpy masked arrays.
    format : str, optional
        The format used to display each value.
    delim : str, optional
        The separator for a range.
    eps : number, optional
        This value is unused.

    Raises
    ------
    ValueError
        If the ``lovals`` and ``hivals`` sequences do not match.

    See Also
    --------
    create_expr, parse_expr

    Examples
    --------

    When there is no mask, or all mask values are True, we just show
    the full range:

    >>> create_expr_integrated([1, 2, 3, 4], [2, 3, 4, 5])
    '1-5'
    >>> create_expr_integrated([1, 2, 4, 5, 7], [2, 3, 5, 6, 8])
    '1-8'
    >>> create_expr_integrated([0.1, 0.2, 0.4, 0.8], [0.2, 0.4, 0.8, 1.0])
    '0.1-1.0'
    >>> create_expr_integrated([0.1, 0.2, 0.4, 0.8], [0.2, 0.4, 0.6, 1.0], [True, True, True, True])
    '0.1-1.0'

    If a mask is given then it defines the bins that are grouped
    together, even if the bins are not contiguous:

    >>> create_expr_integrated([1, 2, 4], [2, 3, 5], [True, True, False, True])
    '1-3,4-5'
    >>> create_expr_integrated([1, 3, 5], [2, 4, 6], [True, True, False])
    '1-4,5-6'

    More examples of the mask controlling the grouping:

    >>> create_expr_integrated([0.1, 0.2, 0.6, 0.8], [0.2, 0.4, 0.8, 1.0], [True, True, False, True, True])
    '0.1-0.4,0.6-1.0'
    >>> create_expr_integrated([0.1, 0.2, 0.4, 0.8], [0.2, 0.3, 0.5, 1.0], [True, True, False, True, False, True])
    '0.1-0.3,0.4-0.5,0.8-1.0'
    >>> create_expr_integrated([0.1, 0.2, 0.4, 0.8], [0.2, 0.3, 0.5, 1.0], [False, True, True, False, True, False, True, False])
    '0.1-0.3,0.4-0.5,0.8-1.0'

    An interesting case is that you can add a "break" between
    contiguous bins (this behavior may be changed):

    >>> create_expr_integrated([1, 2, 3, 4], [2, 3, 4, 5], [True, False, True, True, True])
    '1-2,2-5'

    """

    # Follow create_expr.
    #
    if len(lovals) != len(hivals):
        raise ValueError("hivals array mismatch with lovals")

    if len(lovals) == 0:
        return ''

    # To identify where there's a break we use an array of consecutive
    # integers that have missing data masked out.
    #
    if mask is None:
        seq = np.arange(len(lovals))
    else:
        mask = np.asarray(mask, dtype=bool)

        if sum(mask) != len(lovals):
            raise ValueError("mask array mismatch with lovals")

        seq = np.arange(len(mask))
        seq = seq[mask]

    out = format % lovals[0]

    startbins = lovals[1:]
    endbins = hivals[:-1]

    diffs = np.diff(seq)
    idxs, = np.where(diffs != 1)
    for idx in idxs:
        out += f"{delim}{format % endbins[idx]},{format % startbins[idx]}"

    out += f"{delim}{format % hivals[-1]}"
    return out


def parse_expr(expr):
    """Convert a filter expression into its parts.

    This is intended for parsing a notice or ignore expression
    given as a string.

    Parameters
    ----------
    expr : str
        The filter expression, of the form 'a:b' or a single number,
        separated by commas, and white space is ignored. The
        upper or lower limit of a pair may be ignored (e.g. 'a:' or
        ':b').

    Returns
    -------
    filters : list of pairs
        Each pair gives the lower- and upper-edge of the filter,
        using ``None`` to represent no limit.

    See Also
    --------
    create_expr, create_expr_int

    Notes
    -----
    There is no attempt to validate that the expression contains
    strictly ordered pairs, or that the pairs do not overlap, or
    that the lower- and upper-limits are in increasing numerical
    order. That is, the expression '5:7,:2,4:6,5:3' is allowed.

    Examples
    --------

    >>> parse_expr('0.5:7')
    [(0.5, 7.0)]

    >>> parse_expr('0.5:')
    [(0.5, None)]

    >>> parse_expr(':7')
    [(None, 7.0)]

    >>> parse_expr(':2, 4 : 5 ,7:8,10:')
    [(None, 2.0), (4.0, 5.0), (7.0, 8.0), (10.0, None)]

    >>> parse_expr('4')
    [(4.0, 4.0)]

    >>> parse_expr(' ')
    [(None, None)]

    """

    if expr is None or str(expr).strip() == '':
        return [(None, None)]

    res = []
    vals = str(expr).strip().split(',')
    for val in vals:
        lo, hi = None, None

        interval = val.strip().split(':')
        ninterval = len(interval)
        if ninterval == 1:
            lo = interval[0]
            if lo == '':
                lo = None
            hi = lo
        elif ninterval == 2:
            lo = interval[0]
            hi = interval[1]
            if lo == '':
                lo = None
            if hi == '':
                hi = None
        else:
            # This check exited but was never hit due to the way the
            # code was written. It now errors out if a user gives
            # a:b:c, whereas the old version would have just ignored
            # the ':c' part. Perhaps we should just keep dropping
            # it, in case there's existing code that assumes this?
            #
            raise TypeError("interval syntax requires a tuple, 'lo:hi'")

        if lo is not None:
            try:
                lo = float(lo)
            except ValueError:
                raise TypeError(f"Invalid lower bound '{lo}'") from None
        if hi is not None:
            try:
                hi = float(hi)
            except ValueError:
                raise TypeError(f"Invalid upper bound '{hi}'") from None

        res.append((lo, hi))

    return res


def calc_total_error(staterror=None, syserror=None):
    """Add statistical and systematic errors in quadrature.

    Parameters
    ----------
    staterror : array, optional
       The statistical error, or ``None``.
    syserror : array, optional
       The systematic error, or ``None``.

    Returns
    -------
    error : array or ``None``
       The errors, added in quadrature. If both ``staterror`` and
       ``syserror`` are ``None`` then the return value is ``None``.

    """

    if (staterror is None) and (syserror is None):
        error = None
    elif (staterror is not None) and (syserror is None):
        error = staterror
    elif (staterror is None) and (syserror is not None):
        error = syserror
    else:
        error = np.sqrt(staterror * staterror + syserror * syserror)
    return error


def quantile(sorted_array, f):
    """Return the quantile element from sorted_array, where f is [0,1]
    using linear interpolation.

    Based on the description of the GSL routine
    gsl_stats_quantile_from_sorted_data - e.g.
    http://www.gnu.org/software/gsl/manual/html_node/Median-and-Percentiles.html
    but all errors are my own.

    sorted_array is assumed to be 1D and sorted.
    """
    sorted_array = np.asarray(sorted_array)

    if len(sorted_array.shape) != 1:
        raise RuntimeError("Error: input array is not 1D")
    n = sorted_array.size

    q = (n - 1) * f
    i = int(np.floor(q))
    delta = q - i

    return (1.0 - delta) * sorted_array[i] + delta * sorted_array[i + 1]


def get_error_estimates(x, sorted=False):
    """Compute the median and (-1,+1) sigma values for the data.

    Parameters
    ----------
    x : array of numbers
       The input values.
    sorted : bool, optional
       If ``False``, the default, then ``x`` is assumed to not be sorted.

    Returns
    -------
    (median, lsig, usig)
       The median, value that corresponds to -1 sigma, and value that
       is +1 sigma, for the input distribution.

    Examples
    --------
    >>> (m, l, h) = get_error_estimates(x)

    """
    xs = np.asarray(x)
    if not sorted:
        xs.sort()
        xs = np.array(xs)

    sigfrac = 0.682689
    median = quantile(xs, 0.5)
    lval = quantile(xs, (1 - sigfrac) / 2.0)
    hval = quantile(xs, (1 + sigfrac) / 2.0)

    return (median, lval, hval)


def multinormal_pdf(x, mu, sigma):
    """The PDF of a multivariate-normal distribution.

    Returns the probability density function (PDF) of a
    multivariate normal [1]_ distribution.

    Parameters
    ----------
    x : array
      An array of length k.
    mu : array
      An array of length k.
    sigma : array
      A matrix of size (k,k). It must be symmetric and positive-definite.

    See Also
    --------
    multit_pdf

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Multivariate_normal_distribution

    """
    x = np.asarray(x)
    mu = np.asarray(mu)
    sigma = np.asarray(sigma)
    if x.size != mu.size:
        raise TypeError("x and mu sizes do not match")
    if mu.size != sigma.diagonal().size:
        raise TypeError("sigma shape does not match x")
    if np.min(np.linalg.eigvalsh(sigma)) <= 0:
        raise ValueError("sigma is not positive definite")
    if np.max(np.abs(sigma - sigma.T)) >= 1.e-9:
        raise ValueError("sigma is not symmetric")
    rank = mu.size
    coeff = 1.0 / (np.power(2.0 * np.pi, rank / 2.0) *
                   np.sqrt(np.abs(np.linalg.det(sigma))))

    xmu = np.asarray(x - mu)
    invsigma = np.asarray(np.linalg.inv(sigma))

    # The matrix multiplication looks backwards, but mu and x
    # are passed in already transposed.
    #
    #  mu = [[a,b,c]]
    #   x = [[d,e,f]]
    #
    out = coeff * np.exp(-0.5 * ((xmu @ invsigma) @ xmu.T))
    return float(out)


def multit_pdf(x, mu, sigma, dof):
    """The PDF of a multivariate student-t distribution.

    Returns the probability density function (PDF) of a
    multivariate student-t [1]_ distribution.

    Parameters
    ----------
    x : array
      An array of length k.
    mu : array
      An array of length k.
    sigma : array
      A matrix of size (k,k). It must be symmetric and positive-definite.
    dof : int

    See Also
    --------
    multinormal_pdf

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Multivariate_Student_distribution

    """
    n = float(dof)
    x = np.asarray(x)
    mu = np.asarray(mu)
    sigma = np.asarray(sigma)

    if x.size != mu.size:
        raise TypeError("x and mu sizes do not match")
    if mu.size != sigma.diagonal().size:
        raise TypeError("sigma shape does not match x")
    if np.min(np.linalg.eigvalsh(sigma)) <= 0:
        raise ValueError("sigma is not positive definite")
    if np.max(np.abs(sigma - sigma.T)) >= 1.e-9:
        raise ValueError("sigma is not symmetric")

    rank = mu.size
    npr = float(n + rank)
    coeff = (gamma(npr / 2.0) /
             (gamma(n / 2.0) * np.power(n, rank / 2.0) *
                 np.power(np.pi, rank / 2.0) *
                 np.sqrt(np.abs(np.linalg.det(sigma)))))

    xmu = np.asarray(x - mu)
    invsigma = np.asarray(np.linalg.inv(sigma))

    # The matrix multiplication looks backwards, but mu and x
    # are passed in already transposed.
    #
    #  mu = [[a,b,c]]
    #   x = [[d,e,f]]
    #
    term = 1.0 + 1.0 / n * ((xmu @ invsigma) @ xmu.T)
    out = coeff * np.power(term, -npr / 2.0)
    return float(out)


def dataspace1d(start, stop, step=1, numbins=None):
    """
    Populates an integrated grid

    if numbins is None (default) -> numpy.arange(start,stop,step)

    if numbins is not None -> numpy.linspace(start, stop, numbins)

    """
    if start >= stop:
        raise TypeError(f"input should be start < stop, found start={start} stop={stop}")

    if numbins is None:
        if step <= 0:
            raise TypeError(f"input should be step > 0, found step={step}")

        if step >= (stop - start):
            raise TypeError(
                f"input has produced less than 2 bins, found start={start} stop={stop} step={step}")

    # xx = np.arange(start, stop, step, dtype=float)
    # xx = sao_arange(start, stop, step)
    xx = None
    if numbins is not None:
        if numbins <= 1:
            raise TypeError(
                f"input should be numbins > 1, found numbins={numbins}")

        xx = np.linspace(start, stop, numbins + 1)
    else:
        xx = sao_arange(start, stop, step)

    xlo = np.array(xx[:-1])
    xhi = np.array(xx[1:])
    y = np.zeros(len(xlo), dtype=float)

    return xlo, xhi, y


def dataspace2d(dim):
    """
    Populates a blank image dataset
    """
    if not np.iterable(dim):
        raise TypeError("dim must be an array of dimensions")

    if len(dim) < 2:
        raise TypeError("dimensions for dataspace2d must be > 1")

    if dim[0] < 1 or dim[1] < 1:
        raise TypeError(f"dimensions should be > 0, found dim0 {dim[0]} dim1 {dim[1]}")

    x0 = np.arange(dim[0], dtype=float) + 1.0
    x1 = np.arange(dim[1], dtype=float) + 1.0

    x0, x1 = np.meshgrid(x0, x1)
    shape = tuple(x0.shape)
    x0 = x0.ravel()
    x1 = x1.ravel()
    y = np.zeros(np.prod(dim))

    return x0, x1, y, shape


def histogram1d(x, x_lo, x_hi):
    """Create a 1D histogram from a sequence of samples.

    See the `numpy.histogram` routine for a version with more options.

    .. versionchanged:: 4.15.0
       The x_lo and x_hi arguments are no-longer changed (sorted) by
       this routine.

    Parameters
    ----------
    x : sequence of numbers
       The array of samples
    x_lo : sequence of numbers
       The lower-edges of each bin.
    x_hi : sequence of numbers
       The upper-edges of each bin, which must be the same size
       as ``x_lo``.

    Returns
    -------
    y : NumPy array
       The number of samples in each histogram bin defined by
       the ``x_lo`` and ``x_hi`` arrays.

    Examples
    --------

    A simple example, calculating the histogram of 1000 values
    randomly distributed over [0, 1).

    >>> import numpy as np
    >>> rng = np.random.default_rng()
    >>> x = rng.random(1000)
    >>> edges = np.linspace(0, 1, 11)
    >>> xlo = edges[:-1]
    >>> xhi = edges[1:]
    >>> y = histogram1d(x, xlo, xhi)

    Given a list of samples, bin them up so that they can be used as
    the dependent axis (the value to be fitted) in a Sherpa data set:

    >>> dataspace1d(1, 10, 1)
    >>> lo, hi = get_indep()
    >>> n = histogram1d([2, 3, 2, 8, 5, 2], lo, hi)
    >>> set_dep(n)

    """

    x_lo = np.asarray(x_lo).copy()
    x_hi = np.asarray(x_hi).copy()

    x_lo.sort()
    x_hi.sort()

    return hist1d(np.asarray(x), x_lo, x_hi)


def histogram2d(x, y, x_grid, y_grid):
    """Create 2D histogram from a sequence of samples.

    See the `numpy.histogram2d` routine for a version with more options.

    .. versionchanged:: 4.15.0
       The x_grid and y_grid arguments are no-longer changed (sorted)
       by this routine.

    Parameters
    ----------
    x : sequence of numbers
       The array of samples (X coordinate)
    y : sequence of numbers
       The array of samples (Y coordinate), which must have the same
       size as the ``x`` sequence.
    x_grid : sequence of numbers
       The X bin edges.
    y_grid : sequence of numbers
       The Y bin edges.

    Returns
    -------
    y : NumPy array
       The number of samples in each histogram bin defined by
       the ``x_grid`` and ``y_grid`` arrays.

    Examples
    --------

    Given a list of coordinates (``xvals``, ``yvals``), bin
    them up so that they match the 5 by 10 pixel image
    data space. In this case the X grid is [1, 2, ..., 5]
    and the Y grid is [1, 2, .., 10].

    >>> dataspace2d([5, 10])
    >>> (xgrid, ygrid) = get_axes()
    >>> n = histogram2d(xvals, yvals, xgrid, ygrid)
    >>> set_dep(n)

    """
    x_grid = np.asarray(x_grid).copy()
    y_grid = np.asarray(y_grid).copy()

    x_grid.sort()
    y_grid.sort()

    vals = hist2d(np.asarray(x), np.asarray(y), x_grid, y_grid)
    return vals.reshape((len(x_grid), len(y_grid)))


def interp_util(xout, xin, yin):
    lenxin = len(xin)

    i1 = np.searchsorted(xin, xout)

    i1[i1 == 0] = 1
    i1[i1 == lenxin] = lenxin - 1

#     if 0 == i1:
#         i1 = 1
#     if lenxin == i1:
#         i1 = lenxin - 1

    x0 = xin[i1 - 1]
    x1 = xin[i1]
    y0 = yin[i1 - 1]
    y1 = yin[i1]
    return x0, x1, y0, y1


def linear_interp(xout, xin, yin):
    """Linear one-dimensional interpolation.

    Parameters
    ----------
    xout : array_like
       The positions at which to interpolate.
    xin : array_like
       The x values of the data to interpolate. This must be
       sorted so that it is monotonically increasing.
    yin : array_like
       The y values of the data to interpolate (must be the same
       size as ``xin``).

    Returns
    -------
    yout : NumPy array of numbers
       The interpolated y values (same size as ``xout``).

    See Also
    --------
    interpolate, nearest_interp, neville

    Examples
    --------
    >>> import numpy as np
    >>> x = np.asarray([1.2, 3.4, 4.5, 5.2])
    >>> y = np.asarray([12.2, 14.4, 16.8, 15.5])
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = linear_interp(xgrid, x, y)
    """

    x0, x1, y0, y1 = interp_util(xout, xin, yin)
    val = (xout - x0) / (x1 - x0) * (y1 - y0) + y0
    if np.isnan(val).any():
        # to handle the case where two adjacent elements of xout are equal
        return nearest_interp(xout, xin, yin)
    return val


def nearest_interp(xout, xin, yin):
    """Nearest-neighbor one-dimensional interpolation.

    Parameters
    ----------
    xout : array_like
       The positions at which to interpolate.
    xin : array_like
       The x values of the data to interpolate. This must be
       sorted so that it is monotonically increasing.
    yin : array_like
       The y values of the data to interpolate (must be the same
       size as ``xin``).

    Returns
    -------
    yout : NumPy array of numbers
       The interpolated y values (same size as ``xout``).

    See Also
    --------
    interpolate, linear_interp, neville

    Examples
    --------
    >>> import numpy as np
    >>> x = np.asarray([1.2, 3.4, 4.5, 5.2])
    >>> y = np.asarray([12.2, 14.4, 16.8, 15.5])
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = nearest_interp(xgrid, x, y)
    """

    x0, x1, y0, y1 = interp_util(xout, xin, yin)
    return np.where((np.abs(xout - x0) < np.abs(xout - x1)), y0, y1)


def interpolate(xout, xin, yin, function=linear_interp):
    """One-dimensional interpolation.

    Parameters
    ----------
    xout : array_like
       The positions at which to interpolate.
    xin : array_like
       The x values of the data to interpolate. This must be
       sorted so that it is monotonically increasing.
    yin : array_like
       The y values of the data to interpolate (must be the same
       size as ``xin``).
    function : func, optional
       The function to perform the interpolation. It accepts
       the arguments (xout, xin, yin) and returns the interpolated
       values. The default is to use linear interpolation.

    Returns
    -------
    yout : array_like
       The interpolated y values (same size as ``xout``).

    See Also
    --------
    linear_interp, nearest_interp, neville

    Examples
    --------

    Use linear interpolation to calculate the Y values for the
    ``xgrid`` array:

    >>> import numpy as np
    >>> x = np.asarray([1.2, 3.4, 4.5, 5.2])
    >>> y = np.asarray([12.2, 14.4, 16.8, 15.5])
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = interpolate(xgrid, x, y)

    Use Neville's algorithm for the interpolation:

    >>> ygrid = interpolate(xgrid, x, y, neville)
    """

    if not callable(function):
        raise TypeError(f"input function '{repr(function)}' is not callable")

    return function(xout, xin, yin)


def is_binary_file(filename):
    """Estimate if a file is a binary file.

    Parameters
    ----------
    filename : str
       The name of the file.

    Returns
    -------
    flag : bool
       Returns True if a non-printable character is found in the first
       1024 bytes of the file.

    Notes
    -----
    For this function, "binary" means the file contains a non-ASCII character.
    """

    # Originally "binary" was defined as a character not being in
    # string.printable. With Python 3, we can also use UnicodeDecodeError
    # as an indicator of a "binary" file, but the check against
    # string.printable is kept in, since this is more restrictive
    # than UnicodeDecodeError.
    #
    try:
        with open(filename, 'r') as fd:
            try:
                lines = fd.readlines(1024)
            except UnicodeDecodeError:
                return True

        if len(lines) == 0:
            return False

        # Are there any non-printable characters in the buffer?
        for line in lines:
            for char in line:
                if char not in string.printable:
                    return True

    except OSError as oe:
       raise IOErr('openfailed', f"unable to open {filename}: {oe}") from oe

    return False


################################# Neville2d ###################################


def neville2d(xinterp, yinterp, x, y, fval):
    """Polynomial two-dimensional interpolation using Neville's method.

    The scheme used for interpolation (Neville's method) is described
    at [1]_, where the interpolation is done first over the Y axis
    and then the X axis.

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Neville%27s_algorithm

    """

    nrow = fval.shape[0]
    # ncol = fval.shape[1]
    tmp = np.zeros(nrow)
    for row in range(nrow):
        tmp[row] = neville(yinterp, y, fval[row])
    return neville(xinterp, x, tmp)

################################## Hessian ####################################


class NumDeriv:

    def __init__(self, func, fval0):
        self.nfev, self.func = func_counter(func)
        self.fval_0 = fval0


class NumDerivCentralOrdinary(NumDeriv):
    """
    Subtract the following Taylor series expansion::

                                             2
                                  '         h  ''            3
    f( x +/- h ) = f( x ) +/-  h f  ( x ) + - f  ( x ) + O( h  )
                                            2

    gives::

                                              '            3
               f( x + h ) - f( x - h ) = 2 h f ( x ) + O( h  )

                 '
    solving for f ( x )::

                     '        f( x + h ) - f( x - h )       2
                    f ( x ) = ----------------------- + O( h  )
                                        2 h

    In addition to the truncation error of order h^2, there is a round off
    error due to the finite numerical precision ~ r f( x )::

             '        f( x + h ) - f( x - h )    r f( x )         2
            f ( x ) = ----------------------- + ---------  +  O( h  )
                                2 h                h

                            r      2
                 Error  ~=  -  + h
                            h
    minimizing the error by differentiating wrt h, the solve for h::

        h ~ r^1/3

    """

    def __init__(self, func, fval0=None):
        NumDeriv.__init__(self, func, fval0)

    def __call__(self, x, h):
        if 0.0 == h:
            return np.inf
        return (self.func(x + h) - self.func(x - h)) / (2.0 * h)


class NumDerivFowardPartial(NumDeriv):

    def __init__(self, func, fval0):
        NumDeriv.__init__(self, func, fval0)

    def __call__(self, x, h, *args):

        if 0.0 == h:
            h = pow(np.finfo(np.float32).eps, 1.0 / 3.0)

        ith = args[0]
        jth = args[1]

        ei = np.zeros(len(x), float)
        ej = np.zeros(len(x), float)

        deltai = h * abs(x[ith])
        if 0.0 == deltai:
            deltai = h
        ei[ith] = deltai

        deltaj = h * abs(x[jth])
        if 0.0 == deltaj:
            deltaj = h
        ej[jth] = deltaj

        fval = self.fval_0
        fval += self.func(x + ei + ej)
        fval -= self.func(x + ei)
        fval -= self.func(x + ej)
        fval /= deltai * deltaj
        return fval


class NumDerivCentralPartial(NumDeriv):
    """Add the following Taylor series expansion::

                                             2
                                  '         h  ''            3
    f( x +/- h ) = f( x ) +/-  h f  ( x ) + - f  ( x ) + O( h  )
                                            2
                   ''
    and solve for f  ( x ), gives::

              ''         f( x + h ) + f( x - h ) - 2 f( x )        2
             f  ( x ) = ------------------------------------ + O( h  )
                                         2
                                        h

    In addition to the truncation error of order h^2, there is a round off
    error due to the finite numerical precision ~ r f( x )::

     ''         f( x + h ) + f( x - h ) - 2 f( x )    r f( x )       2
    f  ( x ) = ------------------------------------ + -------- + O( h  )
                                2                        2
                               h                        h

                            r      2
                 Error  ~=  -  + h
                             2
                            h

    minimizing the error by differentiating wrt h, the solve for h::

        h ~ r^1/4
    """

    def __init__(self, func, fval0):
        NumDeriv.__init__(self, func, fval0)

    def __call__(self, x, h, *args):

        if 0.0 == h:
            h = pow(np.finfo(np.float32).eps, 1.0 / 3.0)

        ith = args[0]
        jth = args[1]

        ei = np.zeros(len(x), float)

        if ith == jth:

            delta = h * abs(x[ith])
            if 0.0 == delta:
                delta = h
            ei[ith] = delta

            fval = - 2.0 * self.fval_0
            fval += self.func(x + ei) + self.func(x - ei)
            fval /= delta * delta
            return fval

        ej = np.zeros(len(x), float)

        deltai = h * abs(x[ith])
        if 0.0 == deltai:
            deltai = h
        ei[ith] = deltai

        deltaj = h * abs(x[jth])
        if 0.0 == deltaj:
            deltaj = h
        ej[jth] = deltaj

        fval = self.func(x + ei + ej)
        fval -= self.func(x + ei - ej)
        fval -= self.func(x - ei + ej)
        fval += self.func(x - ei - ej)
        fval /= (4.0 * deltai * deltaj)
        return fval


class NoRichardsonExtrapolation:

    def __init__(self, sequence, verbose=False):
        self.sequence = sequence
        self.verbose = verbose

    def __call__(self, x, t, tol, maxiter, h, *args):
        self.sequence(x, h, *args)


class RichardsonExtrapolation(NoRichardsonExtrapolation):
    """From Wikipedia, the free encyclopedia
    In numerical analysis, Richardson extrapolation is a sequence acceleration
    method, used to improve the rate of convergence of a sequence. It is named
    after Lewis Fry Richardson, who introduced the technique in the early 20th
    century.[1][2] In the words of Birkhoff and Rota, '... its usefulness for
    practical computations can hardly be overestimated.'
    1. Richardson, L. F. (1911). \"The approximate arithmetical solution by
    finite differences of physical problems including differential equations,
    with an application to the stresses in a masonry dam \". Philosophical
    Transactions of the Royal Society of London, Series A 210.
    2. Richardson, L. F. (1927). \" The deferred approach to the limit \".
    Philosophical Transactions of the Royal Society of London, Series A 226:"""

    def __call__(self, x, t, tol, maxiter, h, *args):

        richardson = np.zeros((maxiter, maxiter), dtype=np.float64)
        richardson[0, 0] = self.sequence(x, h, *args)

        t_sqr = t * t
        for ii in range(1, maxiter):
            h /= t
            richardson[ii, 0] = self.sequence(x, h, *args)
            ii_1 = ii - 1
            for jj in range(1, ii + 1):
                # jjp1 = jj + 1  -- this variable is not used
                jj_1 = jj - 1
                factor = pow(t_sqr, jj)
                factor_1 = factor - 1
                richardson[ii, jj] = (factor * richardson[ii, jj_1] -
                                      richardson[ii_1, jj_1]) / factor_1
                arg_jj = richardson[ii, jj]
                arg_jj -= richardson[ii, jj_1]
                arg_ii = richardson[ii, jj]
                arg_ii -= richardson[ii_1, jj_1]
                if Knuth_close(richardson[ii, ii],
                               richardson[ii_1, ii_1], tol):
                    if self.verbose:
                        print_low_triangle(richardson, jj)
                    return richardson[ii, ii]

        if self.verbose:
            print_low_triangle(richardson, maxiter - 1)
        return richardson[maxiter - 1, maxiter - 1]


def hessian(func, par, extrapolation, algorithm, maxiter, h, tol, t):

    num_dif = algorithm(func, func(par))
    deriv = extrapolation(num_dif)
    npar = len(par)
    Hessian = np.zeros((npar, npar), dtype=np.float64)
    for ii in range(npar):
        for jj in range(ii + 1):
            answer = deriv(par, t, tol, maxiter, h, ii, jj)
            Hessian[ii, jj] = answer / 2.0
            Hessian[jj, ii] = Hessian[ii, jj]
    return Hessian, num_dif.nfev[0]


def print_low_triangle(matrix, num):
    # print matrix
    for ii in range(num):
        print(matrix[ii, 0], end=' ')
        for jj in range(1, ii + 1):
            print(matrix[ii, jj], end=' ')
        print()


def symmetric_to_low_triangle(matrix, num):
    low_triangle = []
    for ii in range(num):
        for jj in range(ii + 1):
            low_triangle.append(matrix[ii, jj])
    # print_low_triangle( matrix, num )
    # print low_triangle
    return low_triangle


############################### Root of all evil ##############################


def printf(format, *args):
    """Format args with the first argument as format string, and write.
    Return the last arg, or format itself if there are no args."""
    sys.stdout.write(str(format) % args)
    # WARNING: where is if_ meant to be defined?
    return if_(args, args[-1], format)


# With ParamSpec, added in Python 3.10, we might be able to annotate
# this so that we can match the arguments that `func` uses are the
# same as are sent to the __call__ method, although it might be easier
# to do after the generics changes added in Python 3.12.
#
T = TypeVar("T")


class FuncCounter(Generic[T]):
    """Store the number of times the function is called.

    .. versionadded:: 4.17.0

    """

    __slots__ = ("nfev", "func")

    def __init__(self, func: Callable[..., T]) -> None:
        self.nfev = 0
        self.func = func

    def __call__(self, *args) -> T:
        self.nfev += 1
        return self.func(*args)


def func_counter(func):
    """DEPRECATED.

    .. deprecated:: 4.17.0
       Use the `FuncCounter` class instead.

    """

    # This is FutureWarning rather than DeprecationWarning to make
    # sure users see the message.
    #
    warnings.warn("func_counter is deprecated in 4.17.0: use FuncCounter instead",
                  FutureWarning)

    def func_counter_wrapper(x, *args):
        nfev[0] += 1
        return func(x, *args)

    return nfev, func_counter_wrapper


def is_in(arg, seq):
    """DEPRECATED.

    .. deprecated:: 4.17.0
       Use the Python `in` operator instead.

    """
    # This is FutureWarning rather than DeprecationWarning to make
    # sure users see the message.
    #
    warnings.warn("is_in is deprecated in 4.17.0: use Python's in instead",
                  FutureWarning)
    for x in seq:
        if arg == x:
            return True

    return False


def is_iterable(arg):
    return isinstance(arg, (list, tuple, np.ndarray)) or np.iterable(arg)


# Can this return TypeGuard[Sequence]?
def is_iterable_not_str(arg: Any) -> bool:
    """It is iterable but not a string."""

    return not isinstance(arg, str) and isinstance(arg, Iterable)


def is_sequence(start, mid, end):
    return start < mid < end


def Knuth_close(x, y, tol, myop=operator.__or__):
    """Check whether two floating-point numbers are close together.

    See Also
    --------
    sao_fcmp

    Notes
    -----
    The following text was taken verbatim from [1]_:

    In most cases it is unreasonable to use an operator==(...)
    for a floating-point values equality check. The simple solution
    like ``abs(f1-f2) <= e`` does not work for very small or very big values.
    This floating-point comparison algorithm is based on the more
    confident solution presented by D. E. Knuth in 'The art of computer
    programming (vol II)'. For a given floating point values u and v and
    a tolerance e::

       | u - v | <= e * |u| and | u - v | <= e * |v|                    (1)

    defines a "very close with tolerance e" relationship between u and v::

       | u - v | <= e * |u| or   | u - v | <= e * |v|                   (2)

    defines a "close enough with tolerance e" relationship between
    u and v. Both relationships are commutative but are not transitive.
    The relationship defined by inequations (1) is stronger that the
    relationship defined by inequations (2) (i.e. (1) => (2) ).

    References
    ----------

    .. [1] http://www.boost.org/doc/libs/1_35_0/libs/test/doc/components/test_tools/floating_point_comparison.html#Introduction


    """

    diff = abs(x - y)
    if 0.0 == x or 0.0 == y:
        return diff <= tol
    return myop(diff <= tol * abs(x), diff <= tol * abs(y))


def safe_div(num, denom):

    dbl_max = sys.float_info.max
    dbl_min = sys.float_info.min

    # avoid overflow
    if denom < 1 and num > denom * dbl_max:
        return dbl_max

    # avoid underflow
    if 0.0 == num or denom > 1 and num < denom * dbl_min:
        return 0

    return num / denom


def Knuth_boost_close(x, y, tol, myop=operator.__or__):
    """ The following text was taken verbatim from:

    http://www.boost.org/doc/libs/1_35_0/libs/test/doc/components/test_tools/floating_point_comparison.html#Introduction

    In most cases it is unreasonable to use an operator==(...)
    for a floating-point values equality check. The simple solution
    like abs(f1-f2) <= e does not work for very small or very big values.
    This floating-point comparison algorithm is based on the more
    confident solution presented by D. E. Knuth in 'The art of computer
    programming (vol II)'. For a given floating point values u and v and
    a tolerance e:

    | u - v | <= e * |u| and | u - v | <= e * |v|                    (1)
    defines a "very close with tolerance e" relationship between u and v

    | u - v | <= e * |u| or   | u - v | <= e * |v|                   (2)
    defines a "close enough with tolerance e" relationship between
    u and v. Both relationships are commutative but are not transitive.
    The relationship defined by inequations (1) is stronger that the
    relationship defined by inequations (2) (i.e. (1) => (2) ).
    Because of the multiplication in the right side of inequations,
    that could cause an unwanted underflow condition, the implementation
    is using modified version of the inequations (1) and (2) where all
    underflow, overflow conditions could be guarded safely:

    | u - v | / |u| <= e and | u - v | / |v| <= e          	     (1`)
    | u - v | / |u| <= e or  | u - v | / |v| <= e                    (2`)"""

    diff = abs(x - y)
    if 0.0 == x or 0.0 == y:
        return diff <= tol
    diff_x = safe_div(diff, x)
    diff_y = safe_div(diff, y)
    return myop(diff_x <= tol, diff_y <= tol)


def list_to_open_interval(arg):
    if not np.iterable(arg):
        return arg

    return f'({arg[0]:e}, {arg[1]:e})'


class OutOfBoundErr(Exception):
    """Indicate an out-of-bounds exception in the error analysis"""
    # Should this just move to sherpa.estmethods?
    pass


class QuadEquaRealRoot:
    """ solve for the real roots of the quadratic equation:
    a * x^2 + b * x + c = 0"""

    def __call__(self, a, b, c):

        if 0.0 == a:
            #
            # 0 * x^2 + b * x + c = 0
            #

            if 0.0 != b:
                #
                # 0 * x^2 + b * x + c = 0
                # the folowing still works even if c == 0
                #
                answer = - c / b
                return [answer, answer]

            #
            # 0 * x^2 + 0 * x + c = 0
            #
            # a == 0, b == 0, so if c == 0 then all numbers work so
            # returning nan is not right. However if c != 0 then no
            # roots exist.
            #
            return [None, None]

        if 0.0 == b:

            #
            # a * x^2 + 0 * x + c = 0
            #
            if 0.0 == c:

                # a * x^2 + 0 * x + 0 = 0
                return [0.0, 0.0]

            # a * x^2 + 0 * x + c = 0
            if np.sign(a) == np.sign(c):
                return [None, None]

            answer = np.sqrt(c / a)
            return [-answer, answer]

        if 0.0 == c:

            #
            # a * x^2 + b * x + 0 = 0
            #
            return [0.0, - b / a]

        discriminant = b * b - 4.0 * a * c
        debug("disc=%s", discriminant)
        sqrt_disc = np.sqrt(discriminant)
        t = - (b + np.sign(b) * sqrt_disc) / 2.0
        return [c / t, t / a]


def bisection(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=48, tol=1.0e-6):
    """A basic root finding algorithm that uses standard bisection

    Bisection is a relatively slow method for root finding, but it guaranteed to
    work for a continuous function with a root in a bracketed interval; in other
    words the function must undergo a sign change between the bracketing values.

    See https://en.wikipedia.org/wiki/Bisection_method for a description of the
    bisection method.

    Parameters
    ----------
    fcn : callable
        The function with a root. The function signature is ``fcn(x, *args)``.
    xa : float
        Lower limit of the bracketing interval
    xb : float
        Upper limit of the bracketing interval
    fa : float or None
        Function value at ``xa``. This parameter is optional and can be passed
        to save time in cases where ``fcn(xa, *args)`` is already known and
        function evaluation takes a long time. If `None`, it will be
        calculated.
    fb : float or None
        Function value at ``xb``. This parameter is optional and can be passed
        to save time in cases where ``fcn(xb, *args)`` is already known and
        function evaluation takes a long time. If `None`, it will be
        calculated.
    args : tuple
        Additional parameters that will be passed through to ``fcn``.
    maxfev : int
        Maximal number of function evaluations
    tol : float
        The root finding algorithm stops if a value x with
        ``abs(fcn(x)) < tol`` is found.

    Returns
    -------
    out : list
        The output has the form of a list:
        ``[[x, fcn(x)], [x1, fcn(x1)], [x2, fcn(x2)], nfev]`` where ``x`` is
        the location of the root, and ``x1`` and ``x2`` are the previous
        steps. The function value for those steps is returned as well.
        ``nfev`` is the total number of function evaluations.
        If any of those values is not available, ``None`` will be returned
        instead.
    """
    myfcn = FuncCounter(fcn)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], myfcn.nfev]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], myfcn.nfev]

        if np.sign(fa) == np.sign(fb):
            # TODO: is this a useful message for the user?
            warning('%s: %s fa * fb < 0 is not met', __name__, fcn.__name__)
            return [[None, None], [[None, None], [None, None]],
                    myfcn.nfev]

        while myfcn.nfev < maxfev:

            if abs(fa) > tol and abs(fb) > tol:

                xc = (xa + xb) / 2.0
                fc = myfcn(xc, *args)

                if abs(xa - xb) < min(tol * abs(xb), tol / 10.0):
                    return [[xc, fc], [[xa, fa], [xb, fb]], myfcn.nfev]

                if np.sign(fa) != np.sign(fc):
                    xb, fb = xc, fc
                else:
                    xa, fa = xc, fc

            else:
                if abs(fa) <= tol:
                    return [[xa, fa], [[xa, fa], [xb, fb]], myfcn.nfev]

                return [[xb, fb], [[xa, fa], [xb, fb]], myfcn.nfev]

        xc = (xa + xb) / 2.0
        fc = myfcn(xc, *args)
        return [[xc, fc], [[xa, fa], [xb, fb]], myfcn.nfev]

    except OutOfBoundErr:
        return [[None, None], [[xa, fa], [xb, fb]], myfcn.nfev]


# Is this used at all?
def quad_coef(x, f):
    """
    p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
           = f( xc ) + A ( x - xc ) + B ( ( x - xc ) ( x - xc ) +
                                            ( x - xc ) ( xc - xb ) )
           = f( xc ) + ( A + B ( xc - xb ) ) ( x - xc ) + B ( x - xc )^2

           = f( xc ) + C ( x - xc ) + B ( x - xc )^2 ; C = A + B ( xc - xb )

           = f( xc ) + C x - C xc + B ( x^2  - 2 x xc + xc^2 )

           = B x^2 + ( C - 2 * B xc ) x + f( xc ) - C xc  + B xc^2

           = B x^2 + ( C - 2 * B x[2] ) x + f[ 2 ] + x[2] * ( B x[ 2 ] - C )
    """

    [B, C] = transformed_quad_coef(x, f)
    B_x2 = B * x[2]
    return [B, C - 2 * B_x2, f[2] + x[2] * (B_x2 - C)]


def transformed_quad_coef(x, f):
    """
    p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )

       where A and B are the divided differences::

                                   f( xc ) - f( xb )
                               A = -----------------
                                        xc - xb


                               1     ( f( xc ) - f( xb )   f( xb ) - f( xa ) )
                        B = -------  ( ----------------- - ----------------- )
                            xc - xa  (    xc - xb               xb - xa      )

        p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )
               = f( xc ) + A ( x - xc ) + B ( ( x - xc ) ( x - xc ) +
                                            ( x - xc ) ( xc - xb ) )
               = f( xc ) + ( A + B ( xc - xb ) ) ( x - xc ) + B ( x - xc )^2

               = f( xc ) + C ( x - xc ) + B ( x - xc )^2

        where  C = A + B ( xc - xb )

        The root of p( x ), using the quadratic formula::

                            1  (                                   )
                  x - xc = --- ( - C +/- sqrt( C^2 - 4 f( xc ) B ) )
                           2 B (                                   )

        Rationalize the numerator to avoid subtractive cancellation::

                                     2 f( xc )
                  x - xc = -------------------------------
                           C +/- sqrt( C^2 - 4 f( xc ) B )

        The sign should be chosen to maximize the denominator.  Therefore,
        the next point in the iteration is::

                                       2 f( xc )
                  x = xc - --------------------------------------
                           C + sgn( C ) sqrt( C^2 - 4 f( xc ) B )

                       {    -1,  x < 0
        where sgn(x) = {
                       {     1,  x >= 0
    """

    xa, xb, xc = x[0], x[1], x[2]
    fa, fb, fc = f[0], f[1], f[2]

    # What happens if xb_xa or xc_xa are 0? That is, either
    #     xa == xb
    #     xc == xa
    # Is the assumption that this just never happen?
    #
    xc_xb = xc - xb
    fc_fb = fc - fb
    A = fc_fb / xc_xb
    fb_fa = fb - fa
    xb_xa = xb - xa
    xc_xa = xc - xa
    B = (A - fb_fa / xb_xa) / xc_xa
    C = A + B * xc_xb
    return [B, C]


def _get_discriminant(xa, xb, xc, fa, fb, fc):
    """Wrap up code to transformed_quad_coef.

    This is common code that could be added to transformed_quad_coef
    but is left out at the moment, to make it easier to look back
    at code changes. There is no description of the parameters as
    the existing code has none.

    """

    [B, C] = transformed_quad_coef([xa, xb, xc], [fa, fb, fc])
    discriminant = max(C * C - 4.0 * fc * B, 0.0)
    return B, C, discriminant


def demuller(fcn, xa, xb, xc, fa=None, fb=None, fc=None, args=(),
             maxfev=32, tol=1.0e-6):
    """A root-finding algorithm using Muller's method.

    The algorithm is described at https://en.wikipedia.org/wiki/Muller%27s_method.

    ::

        p( x ) = f( xc ) + A ( x - xc ) + B ( x - xc ) ( x - xb )

    Notes
    -----
    The general case::

                                     2 f( x )
                                           n
               x   = x  -  ----------------------------------------
                n+1   n    C  + sgn( C  ) sqrt( C^2 - 4 f( x  ) B )
                            n         n          n          n    n

                           1     ( f( x  ) - f( x   )   f( x   ) - f( x   )  )
                                 (     n         n-1        n-1        n-2   )
                   B  = -------  ( ------------------ - -------------------  )
                    n   x - x    (    x - x                 x   - x          )
                         n   n-2 (     n   n-1               n-1   n-2       )


                        f( x  ) - f( x   )
                            n         n-1
                   A  = -----------------
                    n       x  - x
                             n    n-1

                   C  = A  + B ( x - x   )
                    n    n    n   n   n-1


    The convergence rate for Muller's method can be shown to be
    the real root of the cubic x - x^3, that is::

       p = (a + 4 / a + 1) / 3
       a = (19 + 3 sqrt(33))^1/3

    In other words: O(h^p) where p is approximately 1.839286755.

    Parameters
    ----------
    fcn : callable
        The function with a root. The function signature is ``fcn(x, *args)``.
    xa, xb, xc : float
        Muller's method requires three initial values.
    fa, fb, fc : float or None
        Function values at ``xa``, ``xb``, and ``xc``. These parameters are
        optional and can be passed
        to save time in cases where ``fcn(xa, *args)`` is already known and
        function evaluation takes a long time. If `None`, they will be
        calculated.
    args : tuple
        Additional parameters that will be passed through to ``fcn``.
    maxfev : int
        Maximal number of function evaluations
    tol : float
        The root finding algorithm stops if the function value a value x with
        ``abs(fcn(x)) < tol`` is found.

    Returns
    -------
    out : list
        The output has the form of a list:
        ``[[x, fcn(x)], [x1, fcn(x1)], [x2, fcn(x2)], nfev]`` where ``x`` is
        the location of the root, and ``x1`` and ``x2`` are the previous
        steps. The function value for those steps is returned as well.
        ``nfev`` is the total number of function evaluations.
        If any of those values is not available, ``None`` will be returned
        instead.
    """

    def is_nan(arg):
        if arg != arg:
            return True
        if arg is np.nan:
            return True
        return np.isnan(arg)

    myfcn = FuncCounter(fcn)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], myfcn.nfev]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], myfcn.nfev]

        if fc is None:
            fc = myfcn(xc, *args)
        if abs(fc) <= tol:
            return [[xc, fc], [[xc, fc], [xc, fc]], myfcn.nfev]

        while myfcn.nfev < maxfev:

            B, C, discriminant = _get_discriminant(xa, xb, xc, fa, fb, fc)
            if is_nan(B) or is_nan(C) or \
                    0.0 == C + np.sign(C) * np.sqrt(discriminant):
                return [[None, None], [[None, None], [None, None]],
                        myfcn.nfev]

            xd = xc - 2.0 * fc / (C + np.sign(C) * np.sqrt(discriminant))

            fd = myfcn(xd, *args)

            if abs(fd) <= tol:
                return [[xd, fd], [[None, None], [None, None]],
                        myfcn.nfev]

            xa = xb
            fa = fb
            xb = xc
            fb = fc
            xc = xd
            fc = fd

        return [[xd, fd], [[None, None], [None, None]],
                myfcn.nfev]

    except ZeroDivisionError:

        return [[xd, fd], [[None, None], [None, None]],
                myfcn.nfev]


def new_muller(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32, tol=1.e-6):
    '''Alternative implementation of Mueller's method for root finding

    Parameters
    ----------
    fcn : callable
        The function with a root. The function signature is ``fcn(x, *args)``.
    xa, xb: float
        Muller's method requires three initial values.
    fa, fb: float or None
        Function values at ``xa`` and ``xb``. These parameters are
        optional and can be passed
        to save time in cases where ``fcn(xa, *args)`` is already known and
        function evaluation takes a long time. If `None`, they will be
        calculated.
    args : tuple
        Additional parameters that will be passed through to ``fcn``.
    maxfev : int
        Maximal number of function evaluations
    tol : float
        The root finding algorithm stops if the function value a value x with
        ``abs(fcn(x)) < tol`` is found.

    Returns
    -------
    out : list
        The output has the form of a list:
        ``[[x, fcn(x)], [x1, fcn(x1)], [x2, fcn(x2)], nfev]`` where ``x`` is
        the location of the root, and ``x1`` and ``x2`` are the previous
        steps. The function value for those steps is returned as well.
        ``nfev`` is the total number of function evaluations.
        If any of those values is not available, ``None`` will be returned
        instead.
    '''

    myfcn = FuncCounter(fcn)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], myfcn.nfev]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], myfcn.nfev]

        if np.sign(fa) == np.sign(fb):
            warning('%s: %s fa * fb < 0 is not met', __name__, fcn.__name__)
            return [[None, None], [[None, None], [None, None]],
                    myfcn.nfev]

        while myfcn.nfev < maxfev:

            xc = (xa + xb) / 2.0
            fc = myfcn(xc, *args)
            if abs(fc) <= tol:
                return [[xc, fc], [[xa, fa], [xb, fb]], myfcn.nfev]

            B, C, discriminant = _get_discriminant(xa, xb, xc, fa, fb, fc)

            xd = xc - 2.0 * fc / (C + np.sign(C) * np.sqrt(discriminant))

            fd = myfcn(xd, *args)

            if abs(fd) <= tol:
                return [[xd, fd], [[xa, fa], [xb, fb]], myfcn.nfev]

            if np.sign(fa) != np.sign(fc):
                xb, fb = xc, fc
                continue

            if np.sign(fd) != np.sign(fc) and xc < xd:
                xa, fa = xc, fc
                xb, fb = xd, fd
                continue

            if np.sign(fb) != np.sign(fd):
                xa, fa = xd, fd
                continue

            if np.sign(fa) != np.sign(fd):
                xb, fb = xd, fd
                continue

            if np.sign(fc) != np.sign(fd) and xd < xc:
                xa, fa = xd, fd
                xb, fb = xc, fc
                continue

            if np.sign(fc) != np.sign(fd):
                xa, fa = xc, fc
                continue

        return [[xd, fd], [[xa, fa], [xb, fb]], myfcn.nfev]

    except (ZeroDivisionError, OutOfBoundErr):

        return [[xd, fd], [[xa, fa], [xb, fb]], myfcn.nfev]

#
# /*
#  * Licensed to the Apache Software Foundation (ASF) under one or more
#  * contributor license agreements.  See the NOTICE file distributed with
#  * this work for additional information regarding copyright ownership.
#  * The ASF licenses this file to You under the Apache License, Version 2.0
#  * (the "License"); you may not use this file except in compliance with
#  * the License.  You may obtain a copy of the License at
#  *
#  *      http://www.apache.org/licenses/LICENSE-2.0
#  *
#  * Unless required by applicable law or agreed to in writing, software
#  * distributed under the License is distributed on an "AS IS" BASIS,
#  * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  * See the License for the specific language governing permissions and
#  * limitations under the License.
#  */
#


def apache_muller(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32,
                  tol=1.0e-6):
    '''An alternative implementation of Muller's method for root finding.

    Unlike the rest of Sherpa, this method is available
    the Apache Software Foundation (ASF) licence - see code for this method
    for details.

    Parameters
    ----------
    fcn : callable
        The function with a root. The function signature is ``fcn(x, *args)``.
    xa, xb: float
        Muller's method requires three initial values.
    fa, fb: float or None
        Function values at ``xa`` and ``xb``. These parameters are
        optional and can be passed
        to save time in cases where ``fcn(xa, *args)`` is already known and
        function evaluation takes a long time. If `None`, they will be
        calculated.
    args : tuple
        Additional parameters that will be passed through to ``fcn``.
    maxfev : int
        Maximal number of function evaluations
    tol : float
        The root finding algorithm stops if the function value a value x with
        ``abs(fcn(x)) < tol`` is found.

    Returns
    -------
    out : list
        The output has the form of a list:
        ``[[x, fcn(x)], [x1, fcn(x1)], [x2, fcn(x2)], nfev]`` where ``x`` is
        the location of the root, and ``x1`` and ``x2`` are the previous
        steps. The function value for those steps is returned as well.
        ``nfev`` is the total number of function evaluations.
        If any of those values is not available, ``None`` will be returned
        instead.
    '''

    myfcn = FuncCounter(fcn)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], myfcn.nfev]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], myfcn.nfev]

        if np.sign(fa) == np.sign(fb):
            warning('%s: %s fa * fb < 0 is not met', __name__, fcn.__name__)
            return [[None, None], [[None, None], [None, None]],
                    myfcn.nfev]

        xc = (xa + xb) / 2.0
        fc = myfcn(xc, *args)
        if abs(fc) <= tol:
            return [[xc, fc], [[xc, fc], [xc, fc]], myfcn.nfev]

        xbest, fbest = xa, fa
        if abs(fb) < abs(fa):
            xbest, fbest = xb, fb
        if abs(fc) < abs(fbest):
            xbest, fbest = xc, fc

        oldx = 1.0e128
        while myfcn.nfev < maxfev:


            B, C, discriminant = _get_discriminant(xa, xb, xc, fa, fb, fc)
            den = np.sign(C) * np.sqrt(discriminant)
            xplus = xc - 2.0 * fc / (C + den)
            if C != den:
                xminus = xc - 2.0 * fc / (C - den)
            else:
                xminus = 1.0e128

            if is_sequence(xa, xplus, xb):
                x = xplus
            else:
                x = xminus

            # print 'xa=', xa, '\tx=', x, '\txb=', xb, '\txc=', xc

            # fubar = quad_coef( [xa,xb,xc], [fa,fb,fc] )
            # quad = QuadEquaRealRoot( )
            # print quad( fubar[0], fubar[1], fubar[2] )
            # print

            # sanity check
            if not is_sequence(xa, x, xb):
                x = (xa + xb) / 2.0

            y = myfcn(x, *args)

            if abs(y) < abs(fbest):
                xbest, fbest = x, y
            tolerance = min(tol * abs(x), tol)
            if abs(y) <= tol or abs(x - oldx) <= tolerance:
                return [[x, y], [[xa, fa], [xb, fb]], myfcn.nfev]

            mybisect = (x < xc and (xc - xa) > 0.95 * (xb - xa)) or \
                       (x > xc and (xb - xc) > 0.95 * (xb - xa)) or \
                       (x == xc)

            if not mybisect:
                if x > xc:
                    xa = xc
                    fa = fc
                if x < xc:
                    xb = xc
                    fb = fc
                xc, fc = x, y
                oldx = x

            else:
                xmid = (xa + xb) / 2.0
                fmid = myfcn(xmid, *args)
                if abs(fmid) < abs(fbest):
                    xbest, fbest = xmid, fmid

                if abs(fmid) <= tol:
                    return [[xmid, fmid], [[xa, fa], [xb, fb]], myfcn.nfev]
                if np.sign(fa) + np.sign(fmid) == 0:
                    xb = xmid
                    fb = fmid
                else:
                    xa = xmid
                    fa = fmid

                xc = (xa + xb) / 2.0
                fc = myfcn(xc, *args)
                if abs(fc) < abs(fbest):
                    xbest, fbest = xc, fc

                if abs(fc) <= tol:
                    return [[xc, fc], [[xa, fa], [xb, fb]], myfcn.nfev]
                oldx = 1.0e128

        #
        # maxfev has exceeded, return the minimum so far
        #
        return [[xbest, fbest], [[xa, fa], [xb, fb]], myfcn.nfev]

    #
    # Something drastic has happened
    #
    except (ZeroDivisionError, OutOfBoundErr):

        return [[xbest, fbest], [[xa, fa], [xb, fb]], myfcn.nfev]


def zeroin(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32, tol=1.0e-2):
    """Obtain a zero of a function of one variable using Brent's root finder.

    Return an approximate location for the root with accuracy::

       4*DBL_EPSILON*abs(x) + tol

    using the algorithm from [1]_.

    References
    ----------

    .. [1] G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
           computations. M., Mir, 1980, p.180 of the Russian edition

    Notes
    -----
    The function makes use of a bisection procedure combined with
    a linear or quadratic inverse interpolation.

    At each step the code operates three abscissae - a, b, and c:

      - b - the last and the best approximation to the root
      - a - the last but one approximation
      - c - the last but one or even an earlier approximation such that:

        1) ``|f(b)| <= |f(c)|``
        2) f(b) and f(c) have opposite signs, i.e. b and c encompass
           the root

    Given these abscissae, the code computes two new approximations,
    one by the  bisection procedure and the other one from interpolation
    (if a,b, and c are all different the quadratic interpolation is used,
    linear otherwise). If the approximation obtained by the interpolation
    looks reasonable (i.e. falls within the current interval [b,c], not
    too close to the end points of the interval), the point is accepted
    as a new approximation to the root. Otherwise, the result of the
    bissection is used.

    Parameters
    ----------
    fcn : callable
        The function with a root. The function signature is ``fcn(x, *args)``.
    xa : float
        Lower limit of the bracketing interval
    xb : float
        Upper limit of the bracketing interval
    fa : float or None
        Function value at ``xa``. This parameter is optional and can be passed
        to save time in cases where ``fcn(xa, *args)`` is already known and
        function evaluation takes a long time. If `None`, it will be
        calculated.
    fb : float or None
        Function value at ``xb``. This parameter is optional and can be passed
        to save time in cases where ``fcn(xb, *args)`` is already known and
        function evaluation takes a long time. If `None`, it will be
        calculated.
    args : tuple
        Additional parameters that will be passed through to ``fcn``.
    maxfev : int
        Maximal number of function evaluations
    tol : float
        The root finding algorithm stops if a value x with
        ``abs(fcn(x)) < tol`` is found.

    Returns
    -------
    out : list
        The output has the form of a list:
        ``[[x, fcn(x)], [x1, fcn(x1)], [x2, fcn(x2)], nfev]`` where ``x`` is
        the location of the root, and ``x1`` and ``x2`` are the previous
        steps. The function value for those steps is returned as well.
        ``nfev`` is the total number of function evaluations.
        If any of those values is not available, ``None`` will be returned
        instead.

    """

    myfcn = FuncCounter(fcn)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
            if abs(fa) <= tol:
                return [[xa, fa], [[xa, fa], [xb, fb]], myfcn.nfev]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xa, fa], [xb, fb]], myfcn.nfev]

        if np.sign(fa) == np.sign(fb):
            warning('%s: %s fa * fb < 0 is not met', __name__, fcn.__name__)
            return [[None, None], [[None, None], [None, None]], myfcn.nfev]

        # With NumPy 2.0 the casting rules changed, leading to some
        # behavioural changes in this code. The simplest fix was to
        # make sure DBL_EPSILON did not remain a np.float32 value.
        #
        xc = xa
        fc = fa
        DBL_EPSILON = float(np.finfo(np.float32).eps)
        while myfcn.nfev < maxfev:

            prev_step = xb - xa

            if abs(fc) < abs(fb):
                xa, fa = xb, fb
                xb, fb = xc, fc
                xc, fc = xa, fa

            tol_act = 2.0 * DBL_EPSILON * abs(xb) + tol / 2.0
            new_step = (xc - xb) / 2.0

            if abs(fb) <= tol:
                return [[xb, fb], [[xa, fa], [xb, fb]], myfcn.nfev]

            if abs(new_step) <= tol_act:
                if np.sign(fb) != np.sign(fa):
                    tmp = apache_muller(fcn, xa, xb, fa, fb, args=args,
                                        maxfev=maxfev - myfcn.nfev,
                                        tol=tol)
                    tmp[-1] += myfcn.nfev
                    return tmp

                if np.sign(fb) != np.sign(fc):
                    tmp = apache_muller(fcn, xb, xc, fb, fc, args=args,
                                        maxfev=maxfev - myfcn.nfev,
                                        tol=tol)
                    tmp[-1] += myfcn.nfev
                    return tmp

                return [[xb, fb], [[xa, fa], [xb, fb]], myfcn.nfev]

            if abs(prev_step) >= tol_act and abs(fa) > abs(fb):

                cb = xc - xb
                if xa == xc:
                    t1 = fb / fa
                    p = cb * t1
                    q = 1.0 - t1
                else:
                    t1 = fb / fc
                    t2 = fb / fa
                    q = fa / fc
                    p = t2 * (cb * q * (q - t1) - (xb - xa) * (t1 - 1.0))
                    q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0)

                if p > 0:
                    q = -q
                else:
                    p = -p

                if 2 * p < (1.5 * cb * q - abs(tol_act * q)) and \
                   2 * p < abs(prev_step * q):
                    new_step = p / q

            if abs(new_step) < tol_act:
                if new_step > 0:
                    new_step = tol_act
                else:
                    new_step = -tol_act
            xa = xb
            fa = fb
            xb += new_step
            fb = myfcn(xb, *args)

            if fb > 0 and fc > 0 or fb < 0 and fc < 0:
                xc = xa
                fc = fa

        return [[xb, fb], [[xa, fa], [xc, fc]], myfcn.nfev]

    except (ZeroDivisionError, OutOfBoundErr):
        return [[xb, fb], [[xa, fa], [xc, fc]], myfcn.nfev]


def public(f):
    """Use a decorator to avoid retyping function/class names.

    * Based on an idea by Duncan Booth:
    http://groups.google.com/group/comp.lang.python/msg/11cbb03e09611b8a
    * Improved via a suggestion by Dave Angel:
    http://groups.google.com/group/comp.lang.python/msg/3d400fb22d8a42e1

    See also https://bugs.python.org/issue26632

    """
    _all = sys.modules[f.__module__].__dict__.setdefault('__all__', [])
    if f.__name__ not in _all:  # Prevent duplicates if run from an IDE.
        _all.append(f.__name__)
    return f


def send_to_pager(txt: str, filename=None, clobber: bool = False) -> None:
    """Write out the given string, using pagination if supported.

    This used to call out to using less/more but now is handled
    by pydoc.pager

    Parameters
    ----------
    txt : str
        The text to display
    filename : str or StringIO or None, optional
        If not None, write the output to the given file or filelike
        object.
    clobber : bool, optional
        If filename is a string, then - when clobber is set - refuse
        to overwrite the file if it already exists.

    """

    if filename is None:
        pydoc.pager(txt)
        return

    # Have we been sent a StringIO-like object?
    #
    if hasattr(filename, 'write'):
        print(txt, file=filename)
        return

    # Assume a filename
    clobber = bool_cast(clobber)
    if os.path.isfile(filename) and not clobber:
        raise IOErr('filefound', filename)

    with open(filename, 'w', encoding="UTF-8") as fh:
        print(txt, file=fh)
