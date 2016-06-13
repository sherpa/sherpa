from __future__ import print_function
from __future__ import absolute_import
#
#  Copyright (C) 2007, 2015, 2016  Smithsonian Astrophysical Observatory
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

"""
Objects and utilities used by multiple Sherpa subpackages
"""

import operator
import inspect
from six.moves import zip as izip
from types import FunctionType as function
from types import MethodType as instancemethod
import string
import sys
import os
import importlib
import unittest
import numpy
import numpy.random
import numpy.fft
# Note: _utils.gsl_fcmp is not exported from this module; is this intentional?
from unittest import skipIf
from sherpa.utils._utils import calc_ftest, calc_mlr, igamc, igam, \
    incbet, gamma, lgam, erf, ndtri, sao_fcmp, rebin, \
    hist1d, hist2d, sum_intervals, neville, sao_arange
from sherpa.utils._psf import extract_kernel, normalize, set_origin, \
    pad_bounding_box

from sherpa import get_config
from six.moves.configparser import ConfigParser, NoSectionError
from six.moves import xrange

import logging
warning = logging.getLogger("sherpa").warning
debug = logging.getLogger("sherpa").debug

config = ConfigParser()
config.read(get_config())

_ncpu_val = "NONE"
try:
    _ncpu_val = config.get('parallel', 'numcores').strip().upper()
except NoSectionError:
    pass

_ncpus = None
if not _ncpu_val.startswith('NONE'):
    _ncpus = int(_ncpu_val)

_multi = False

try:
    import multiprocessing
    _multi = True

    if _ncpus is None:
        _ncpus = multiprocessing.cpu_count()
except Exception as e:
    warning("parallel processing is unavailable,\n" +
            "multiprocessing module failed with \n'%s'" % str(e))
    _ncpus = 1
    _multi = False

del _ncpu_val, config, get_config, ConfigParser, NoSectionError


__all__ = ('NoNewAttributesAfterInit', 'SherpaTestCase',
           '_guess_ampl_scale', 'apache_muller', 'bisection', 'bool_cast',
           'calc_ftest', 'calc_mlr', 'calc_total_error', 'create_expr',
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
           'neville', 'neville2d', 'requires_data', 'requires_fits', 'requires_package',
           'new_muller', 'normalize', 'numpy_convolve',
           'pad_bounding_box', 'parallel_map', 'param_apply_limits',
           'parse_expr', 'poisson_noise', 'print_fields', 'rebin',
           'sao_arange', 'sao_fcmp', 'set_origin', 'sum_intervals', 'zeroin',
           'multinormal_pdf', 'multit_pdf', 'get_error_estimates', 'quantile',
           'get_valid_args')

_guess_ampl_scale = 1.e+3

###############################################################################
#
# Types
#
###############################################################################


# Default numeric types (these match the typedefs in extension.hh)
SherpaInt = numpy.intp
SherpaUInt = numpy.uintp
SherpaFloat = numpy.float_

###############################################################################


class NoNewAttributesAfterInit(object):
    """

    Prevents attribute deletion and setting of new attributes after
    __init__ has been called.  Derived classes must call
    NoNewAttributesAfterInit.__init__ after all other initialization.

    """

    __initialized = False  # Use name mangling

    def __init__(self):
        self.__initialized = True

    def __delattr__(self, name):
        if self.__initialized and hasattr(self, name):
            raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                  "deleted") % (type(self).__name__, name))
        object.__delattr__(self, name)

    def __setattr__(self, name, val):
        if self.__initialized and (not hasattr(self, name)):
            raise AttributeError("'%s' object has no attribute '%s'" %
                                 (type(self).__name__, name))

        if self.__initialized and hasattr(self, name):
            if callable(getattr(self, name)) and not callable(val):
                raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                      "replaced with a non-callable attribute")
                                     % (type(self).__name__, name))
            elif not callable(getattr(self, name)) and callable(val):
                raise AttributeError(("'%s' object attribute '%s' cannot be " +
                                      "replaced with a callable attribute") %
                                     (type(self).__name__, name))

        object.__setattr__(self, name, val)


###############################################################################
#
# Testing stuff
#
###############################################################################

def _get_datadir():
    import os
    try:
        import sherpatest
        datadir = os.path.dirname(sherpatest.__file__)
    except ImportError:
        try:
            import sherpa
            datadir = os.path.join(os.path.dirname(sherpa.__file__), os.pardir,
                                   'sherpa-test-data', 'sherpatest')
            if not os.path.exists(datadir) or not os.listdir(datadir):
                # The dir is empty, maybe the submodule was not initialized
                datadir = None
        except ImportError:
            # neither sherpatest nor sherpa can be found, falling back to None
            datadir = None
    return datadir


class SherpaTestCase(unittest.TestCase):
    "Base class for Sherpa unit tests"

    # The location of the Sherpa test data (it is optional)
    datadir = _get_datadir()

    def make_path(self, *segments):
        """Add the segments onto the test data location.

        Parameters
        ----------
        *segments
           Path segments to combine together with the location of the
           test data.

        Returns
        -------
        fullpath : None or string
           The full path to the repository, or None if the
           data directory is not set.

        """
        if self.datadir is None:
            return None
        return os.path.join(self.datadir, *segments)

    # What is the benefit of this over numpy.testing.assert_allclose(),
    # which was added in version 1.5 of NumPy?
    def assertEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """Check that the values are equal within an absolute tolerance.

        Parameters
        ----------
        first : number or array_like
           The expected value, or values.
        second : number or array_like
           The value, or values, to check. If first is an array, then
           second must be an array of the same size. If first is
           a scalar then second can be a scalar or an array.
        tol : number
           The absolute tolerance used for comparison.
        msg : string
           The message to display if the check fails.

        """

        self.assertFalse(numpy.any(sao_fcmp(first, second, tol)), msg)

    def assertNotEqualWithinTol(self, first, second, tol=1e-7, msg=None):
        """Check that the values are not equal within an absolute tolerance.

        Parameters
        ----------
        first : number or array_like
           The expected value, or values.
        second : number or array_like
           The value, or values, to check. If first is an array, then
           second must be an array of the same size. If first is
           a scalar then second can be a scalar or an array.
        tol : number
           The absolute tolerance used for comparison.
        msg : string
           The message to display if the check fails.

        """

        self.assertTrue(numpy.all(sao_fcmp(first, second, tol)), msg)

    # for running regression tests from sherpa-test-data
    def run_thread(self, name, scriptname='fit.py'):
        """Run a regression test from the sherpa-test-data submodule.

        Parameters
        ----------
        name : string
           The name of the science thread to run (e.g., pha_read,
           radpro). The name should match the corresponding thread
           name in the sherpa-test-data submodule. See examples below.
        scriptname : string
           The suffix of the test script file name, usually "fit.py."

        Examples
        --------
        Regression test script file names have the structure
        "name-scriptname.py." By default, scriptname is set to "fit.py."
        For example, if one wants to run the regression test
        "pha_read-fit.py," they would write

        >>> run_thread("pha_read")

        If the regression test name is "lev3fft-bar.py," they would do

        >>> run_thread("lev3fft", scriptname="bar.py")

        """

        self.locals = {}
        cwd = os.getcwd()
        os.chdir(self.datadir)
        scriptname = name + "-" + scriptname
        try:
            execfile(scriptname, {}, self.locals)
        finally:
            os.chdir(cwd)


def requires_data(test_function):
    """
    Decorator for functions requiring external data (i.e. data not distributed with Sherpa
    itself) is missing.  This is used to skip tests that require such data.
    """
    condition = SherpaTestCase.datadir is None
    msg = "required test data missing"
    return skipIf(condition, msg)(test_function)


def _has_package_from_list(*packages):
    """
    Returns True if at least one of the ``packages`` args is importable.
    """
    for package in packages:
        try:
            importlib.import_module(package)
            return True
        except:
            pass
    return False


def requires_package(msg=None, *packages):
    """
    Decorator for test functions requiring specific packages.
    """
    condition = _has_package_from_list(*packages)
    msg = msg or "required module missing among {}.".format(
        ", ".join(packages))

    def decorator(test_function):
        return skipIf(not condition, msg)(test_function)
    return decorator


def requires_xspec(test_function):
    """Decorator for test functions requiring xspec"""
    return requires_package('xspec required', 'sherpa.astro.xspec')(test_function)


def requires_ds9(test_function):
    """Decorator for test functions requiring ds9"""
    return requires_package('ds9 required', 'sherpa.image.ds9_backend')(test_function)


def requires_fits(test_function):
    """
    Returns True if there is an importable backend for FITS I/O.
    Used to skip tests requiring fits_io
    """
    packages = ('pyfits',
                'astropy.io.fits',
                'pycrates',
                )
    msg = "FITS backend required"
    return requires_package(msg, *packages)(test_function)


def requires_pylab(test_function):
    """
    Returns True if the pylab module is available (pylab).
    Used to skip tests requiring matplotlib
    """
    packages = ('pylab',
                )
    msg = "matplotlib backend required"
    return requires_package(msg, *packages)(test_function)


###############################################################################
#
# Utilities
#
###############################################################################


# at what precisions do we assume equality in energy grids?
eps = numpy.finfo(numpy.float32).eps


def filter_bins(mins, maxes, axislist):
    mask = None

    for lo, hi, axis in izip(mins, maxes, axislist):

        if (lo is None) and (hi is None):
            continue

        if lo is None:
            # axismask = axis <= hi
            axismask = (sao_fcmp(hi, axis, eps) >= 0)

        elif hi is None:
            # axismask = axis >= lo
            axismask = (sao_fcmp(lo, axis, eps) <= 0)

        else:
            # axismask = (lo <= axis) & (axis <= hi)
            axismask = (sao_fcmp(lo, axis, eps) <= 0) & \
                       (sao_fcmp(hi, axis, eps) >= 0)

        if mask is None:
            mask = axismask
        else:
            mask &= axismask

    return mask


def bool_cast(val):
    "Converts a string (true|False|on|OFF|etc...)  to a boolean value"

    if type(val) in (tuple, list, numpy.ndarray):
        return numpy.asarray([bool_cast(item) for item in val], bool)

    elif type(val) == str:
        # since built in bool() only returns false for empty strings
        vlo = val.lower()
        if vlo in ('false', 'off', 'no', '0', 'f', 'n'):
            return False

        elif vlo in ('true', 'on', 'yes', '1', 't', 'y'):
            return True

        raise TypeError("unknown boolean value: '%s'" % str(val))

    else:
        # use built in bool cast
        return bool(val)


def export_method(meth, name=None, modname=None):
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

    if type(meth) is not instancemethod:
        return meth

    if name is None:
        name = meth.__name__

    if name == meth.__name__:
        old_name = '_old_' + name
    else:
        old_name = meth.__name__

    # Make an argument list string, removing 'self'
    argspec = inspect.getargspec(meth)
    defaults = meth.__defaults__
    doc = meth.__doc__
    argspec[0].pop(0)
    argspec = inspect.formatargspec(argspec[0], argspec[1], argspec[2])

    # Create a wrapper function with no default arguments
    g = {old_name: meth}
    if modname is not None:
        g['__name__'] = modname
    fdef = 'def %s%s:  return %s%s' % (name, argspec, old_name, argspec)
    exec(fdef, g)

    # Create another new function from the one we just made, this time
    # adding the default arguments and doc string from the original method
    new_meth = g[name]

    new_meth = function(new_meth.__code__, new_meth.__globals__,
                        new_meth.__name__, defaults,
                        new_meth.__closure__)
    new_meth.__doc__ = doc

    return new_meth


def get_keyword_names(func, skip=0):
    """

    Return a list containing the names of func's keyword arguments,
    skipping the first skip keywords.

    """

    argspec = inspect.getargspec(func)
    if argspec[3] is None:
        return []
    first = len(argspec[0]) - len(argspec[3])
    return argspec[0][first + skip:]


def get_keyword_defaults(func, skip=0):
    """

    Return a dictionary containing the default values of func's keyword
    arguments, skipping the first skip keywords.

    """

    argspec = inspect.getargspec(func)
    if argspec[3] is None:
        return {}
    first = len(argspec[0]) - len(argspec[3])
    return dict(izip(argspec[0][first + skip:], argspec[3][skip:]))


def get_num_args(func):
    """

    Return a tuple of the number of arguments.
    ( total number of args, number of non-keyword args, number of keyword args)

    """

    argspec = inspect.getargspec(func)
    num_args = 0
    num_kargs = 0

    if len(argspec[0]) != 0:
        num_args = len(argspec[0])

    if argspec[3] is not None:
        num_kargs = len(argspec[3])

    return (num_args, (num_args - num_kargs), num_kargs)


def print_fields(names, vals, converters={}):
    """

    Given a list of strings names and mapping vals, where names is a
    subset of vals.keys(), return a listing of name/value pairs
    printed one per line in the format '<name> = <value>'.  If a value
    is a NumPy array, print it in the format
    '<data type name>[<array size>]'.  Otherwise, use str(value).

    """

    width = max(len(n) for n in names)
    fmt = '%%-%ds = %%s' % width
    lines = []
    for n in names:
        v = vals[n]

        conv = None
        if converters:
            for t in converters:
                if isinstance(v, t):
                    conv = converters[t]
                    break

        if conv is not None:
            v = conv(v)
        elif isinstance(v, numpy.ndarray):
            v = '%s[%d]' % (numpy.typeNA[v.dtype.type], v.size)
        else:
            v = str(v)
        lines.append(fmt % (n, v))
    return '\n'.join(lines)


def create_expr(vals, mask=None, format='%s', delim='-'):
    """
    collapse a list of channels into an expression using hyphens
    and commas to indicate filtered intervals.
    """
    expr = []

    if len(vals) == 0:
        return ''
    elif len(vals) == 1:
        return format % vals[0]

    diffs = numpy.apply_along_axis(numpy.diff, 0, vals)
    if mask is not None:
        index = numpy.arange(len(mask))
        diffs = numpy.apply_along_axis(numpy.diff, 0, index[mask])

    for ii, delta in enumerate(diffs):
        if ii == 0:
            expr.append(format % vals[ii])
            if delta != 1 or len(diffs) == 1:
                expr.append(',')
            continue
        if delta == 1:
            if expr[-1] == ',':
                expr.append(format % vals[ii])
            if expr[-1] != delim:
                expr.append(delim)
        else:
            if not expr[-1] in (',', delim):
                expr.append(',')
            expr.append(format % vals[ii])
            expr.append(',')
    if len(expr) and expr[-1] in (',', delim):
        expr.append(format % vals[-1])

    return ''.join(expr)


def parse_expr(expr):
    """
    parse a filter expression into numerical components for notice/ignore
    e.g. ':2,4:5,7:8,10:'
    """
    res = []
    if expr is None or str(expr).strip() == '':
        res.append((None, None))
        return res
    vals = str(expr).strip().split(',')
    for val in vals:
        lo, hi = None, None
        interval = val.strip().split(':')
        if len(interval) == 1:
            lo = interval[0]
            if lo == '':
                lo = None
            hi = lo
        elif len(interval) > 1:
            lo = interval[0]
            hi = interval[1]
            if lo == '':
                lo = None
            if hi == '':
                hi = None
        else:
            raise TypeError("interval syntax requires a tuple, 'lo:hi'")
        if lo is not None:
            try:
                lo = float(lo)
            except ValueError:
                raise TypeError("Invalid lower bound '%s'" % str(lo))
        if hi is not None:
            try:
                hi = float(hi)
            except ValueError:
                raise TypeError("Invalid upper bound '%s'" % str(hi))
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
        error = numpy.sqrt(staterror * staterror + syserror * syserror)
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
    sorted_array = numpy.asarray(sorted_array)

    if len(sorted_array.shape) != 1:
        raise RuntimeError("Error: input array is not 1D")
    n = sorted_array.size

    q = (n - 1) * f
    i = numpy.int(numpy.floor(q))
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
    xs = numpy.asarray(x)
    if not sorted:
        xs.sort()
        xs = numpy.array(xs)

    sigfrac = 0.682689
    median = quantile(xs, 0.5)
    lval = quantile(xs, (1 - sigfrac) / 2.0)
    hval = quantile(xs, (1 + sigfrac) / 2.0)

    return (median, lval, hval)


def poisson_noise(x):
    """Draw samples from a Poisson distribution.

    Parameters
    ----------
    x : scalar or array
       The expectation value for the distribution.

    Returns
    -------
    out : scalar or array
       A random realisation of the input array, drawn from
       the Poisson distribution, as a `SherpaFloat`.

    Notes
    -----
    The distribution is calculated by `numpy.poisson.poisson`.

    Examples
    --------
    >>> poisson_noise([10, 20, 5])
    array([ 13.,  21.,   6.])

    """

    x = numpy.asarray(x)

    # Using numpy.where() and indexing doesn't work with 0-d arrays, so
    # handle them separately
    if x.shape == ():
        x = SherpaFloat(x)
        if x <= 0.0:
            x = 0.0
        else:
            x = numpy.random.poisson(x)
        return SherpaFloat(x)

    x_out = numpy.zeros(x.shape, SherpaFloat)
    good = numpy.where(x > 0.0)
    x_out[good] = numpy.random.poisson(x[good])

    return x_out


def multinormal_pdf(x, mu, sigma):
    """The PDF of a multivariate-normal.

    Returns the probability density function (PDF) of a
    multivariate normal [1]_.

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

    .. [1] http://en.wikipedia.org/wiki/Multivariate_normal_distribution

    """
    x = numpy.asarray(x)
    mu = numpy.asarray(mu)
    sigma = numpy.asarray(sigma)
    if x.size != mu.size:
        raise TypeError("x and mu sizes do not match")
    if mu.size != sigma.diagonal().size:
        raise TypeError("sigma shape does not match x")
    if numpy.min(numpy.linalg.eigvalsh(sigma)) <= 0:
        raise ValueError("sigma is not positive definite")
    if numpy.max(numpy.abs(sigma - sigma.T)) >= 1.e-9:
        raise ValueError("sigma is not symmetric")
    rank = mu.size
    coeff = 1.0 / (numpy.power(2.0 * numpy.pi, rank / 2.0) *
                   numpy.sqrt(numpy.abs(numpy.linalg.det(sigma))))
    xmu = numpy.mat(x - mu)
    invsigma = numpy.mat(numpy.linalg.inv(sigma))

    # The matrix multiplication looks backwards, but mu and x
    # are passed in already transposed.
    #
    #  mu = [[a,b,c]]
    #   x = [[d,e,f]]
    #
    return float(coeff * numpy.exp(-0.5 * ((xmu * invsigma) * xmu.T)))


def multit_pdf(x, mu, sigma, dof):
    """The PDF of a multivariate student-t.

    Returns the probability density function (PDF) of a
    multivariate student-t distribution [1]_.

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

    .. [1] http://en.wikipedia.org/wiki/Multivariate_Student_distribution

    """
    n = float(dof)
    x = numpy.asarray(x)
    mu = numpy.asarray(mu)
    sigma = numpy.asarray(sigma)

    if x.size != mu.size:
        raise TypeError("x and mu sizes do not match")
    if mu.size != sigma.diagonal().size:
        raise TypeError("sigma shape does not match x")
    if numpy.min(numpy.linalg.eigvalsh(sigma)) <= 0:
        raise ValueError("sigma is not positive definite")
    if numpy.max(numpy.abs(sigma - sigma.T)) >= 1.e-9:
        raise ValueError("sigma is not symmetric")

    rank = mu.size
    np = float(n + rank)
    coeff = (gamma(np / 2.0) /
             (gamma(n / 2.0) * numpy.power(n, rank / 2.0) *
                 numpy.power(numpy.pi, rank / 2.0) *
                 numpy.sqrt(numpy.abs(numpy.linalg.det(sigma)))))
    xmu = numpy.mat(x - mu)
    invsigma = numpy.mat(numpy.linalg.inv(sigma))

    # The matrix multiplication looks backwards, but mu and x
    # are passed in already transposed.
    #
    #  mu = [[a,b,c]]
    #   x = [[d,e,f]]
    #
    term = 1.0 + 1.0 / n * ((xmu * invsigma) * xmu.T)
    return float(coeff * numpy.power(term, -np / 2.0))


def _convolve(a, b):
    if len(a) != len(b):
        raise TypeError("Input arrays are not equal in length, a: %s b: %s" %
                        (len(a), len(b)))

    imag = numpy.fft.fft(a) * numpy.fft.fft(b)
    return numpy.asarray(numpy.fft.ifft(imag), dtype=SherpaFloat)


def numpy_convolve(a, b):

    if a.ndim > 1 or b.ndim > 1:
        raise TypeError("numpy_convolution is 1D only")

    c = numpy.concatenate((a, numpy.zeros(len(b))))
    d = numpy.concatenate((b, numpy.zeros(len(a))))

    if len(a) > len(b):
        return _convolve(c, d)[:len(a)]

    return _convolve(c, d)[:len(b)]


def dataspace1d(start, stop, step=1, numbins=None):
    """
    Populates an integrated grid

    if numbins is None (default) -> numpy.arange(start,stop,step)

    if numbins is not None -> numpy.linspace(start, stop, numbins)

    """
    if start >= stop:
        raise TypeError("input should be start < stop, found start=%s stop=%s" %
                        (start, stop))

    if numbins is None:
        if step <= 0:
            raise TypeError("input should be step > 0, found step=%s" % step)

        if step >= (stop - start):
            raise TypeError(
                "input has produced less than 2 bins, found start=%s stop=%s step=%s" % (start, stop, step))

    # xx = numpy.arange(start, stop, step, dtype=float)
    # xx = sao_arange(start, stop, step)
    xx = None
    if numbins is not None:
        if numbins <= 1:
            raise TypeError(
                "input should be numbins > 1, found numbins=%s" % numbins)

        xx = numpy.linspace(start, stop, numbins + 1)
    else:
        xx = sao_arange(start, stop, step)

    xlo = numpy.array(xx[:-1])
    xhi = numpy.array(xx[1:])
    y = numpy.zeros(len(xlo), dtype=float)

    return xlo, xhi, y


def dataspace2d(dim):
    """
    Populates a blank image dataset
    """
    if not numpy.iterable(dim):
        raise TypeError("dim must be an array of dimensions")

    if len(dim) < 2:
        raise TypeError("dimensions for dataspace2d must be > 1")

    if dim[0] < 1 or dim[1] < 1:
        raise TypeError("dimensions should be > 0, found dim0 %s dim1 %s"
                        % (dim[0], dim[1]))

    x0 = numpy.arange(dim[0], dtype=numpy.float) + 1.0
    x1 = numpy.arange(dim[1], dtype=numpy.float) + 1.0

    x0, x1 = numpy.meshgrid(x0, x1)
    shape = tuple(x0.shape)
    x0 = x0.ravel()
    x1 = x1.ravel()
    y = numpy.zeros(numpy.prod(dim))

    return x0, x1, y, shape


def histogram1d(x, x_lo, x_hi):
    """Create a 1D histogram from a binned grid (``x_lo``, ``xhi``)
    and array of samples (``x``).

    See the `numpy.histogram` routine for a version with more options.

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
    randomly distributed over [0,1).

    >>> x = np.random.random(1000)
    >>> edges = np.arange(0, 1.1, 0.1)
    >>> xlo = edges[:-1]
    >>> xhi = edges[1:]
    >>> y = histogram1d(x, xlo, xhi)

    Given a list of samples (``vals``), bin them up so that
    they can be used as the dependent axis (the value to
    be fitted) in a Sherpa data set:

    >>> dataspace1d(0.1, 10, 0.1)
    >>> (lo, hi) = get_indep()
    >>> n = histogram1d(vals, lo, hi)
    >>> set_dep(n)

    """

    x_lo = numpy.asarray(x_lo)
    x_hi = numpy.asarray(x_hi)

    x_lo.sort()
    x_hi.sort()

    return hist1d(numpy.asarray(x), x_lo, x_hi)


def histogram2d(x, y, x_grid, y_grid):
    """Create 21D histogram from a binned grid (``x_grid``, ``y_grid``)
    and array of samples (``x``, and ``y``).

    See the `numpy.histogram2d` routine for a version with more options.

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
    data space. In this case the X grid is [1,2,...,5]
    and the Y grid is [1,2,..,10].

    >>> dataspace2d([5, 10])
    >>> (xgrid, ygrid) = get_axes()
    >>> n = histogram2d(xvals, yvals, xgrid, ygrid)
    >>> set_dep(n)

    """
    x_grid = numpy.asarray(x_grid)
    y_grid = numpy.asarray(y_grid)

    x_grid.sort()
    y_grid.sort()

    vals = hist2d(numpy.asarray(x), numpy.asarray(y), x_grid, y_grid)
    return vals.reshape((len(x_grid), len(y_grid)))


def interp_util(xout, xin, yin):
    lenxin = len(xin)

    i1 = numpy.searchsorted(xin, xout)

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
    >>> x = [1.2, 3.4, 4.5, 5.2]
    >>> y = [12.2, 14.4, 16.8, 15.5]
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = linear_interp(xgrid, x, y)
    """

    x0, x1, y0, y1 = interp_util(xout, xin, yin)
    val = (xout - x0) / (x1 - x0) * (y1 - y0) + y0
    if numpy.isnan(val).any():
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
    >>> x = [1.2, 3.4, 4.5, 5.2]
    >>> y = [12.2, 14.4, 16.8, 15.5]
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = nearest_interp(xgrid, x, y)
    """

    x0, x1, y0, y1 = interp_util(xout, xin, yin)
    return numpy.where((numpy.abs(xout - x0) < numpy.abs(xout - x1)), y0, y1)


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

    >>> x = [1.2, 3.4, 4.5, 5.2]
    >>> y = [12.2, 14.4, 16.8, 15.5]
    >>> xgrid = np.linspace(2, 5, 5)
    >>> ygrid = interpolate(xgrid, x, y)

    Use Neville's algorithm for the interpolation:

    >>> ygrid = interpolate(xgrid, x, y, neville)
    """

    if not callable(function):
        raise TypeError("input function '%s' is not callable" %
                        repr(function))

    return function(xout, xin, yin)


def is_binary_file(filename):
    """Estimate if a file is a binary file.

    Returns True if a non-printable character is found in the first
    1024 bytes of the file.
    """
    fd = open(filename, 'r')
    lines = fd.readlines(1024)
    fd.close()

    if len(lines) == 0:
        return False

    # If a non-printable character is found in first 1024 --> binary
    for line in lines:
        for char in line:
            if char not in string.printable:
                return True

    return False


def get_midpoint(a):
    # return numpy.abs(a.max() - a.min())/2. + a.min()
    return numpy.abs(a.max() + a.min()) / 2.0


def get_peak(y, x, xhi=None):
    return x[y.argmax()]


def get_valley(y, x, xhi=None):
    return x[y.argmin()]


def get_fwhm(y, x, xhi=None):
    half_max_val = y.max() / 2.0
    x_max = x[y.argmax()]
    for ii, val in enumerate(y[:y.argmax()]):
        if val >= half_max_val:
            return 2.0 * numpy.abs(x[ii] - x_max)
    return x_max


def guess_fwhm(y, x, xhi=None, scale=1000):
    fwhm = get_fwhm(y, x, xhi)
    return {'val': fwhm, 'min': fwhm / scale, 'max': fwhm * scale}


def param_apply_limits(param_limits, par, limits=True, values=True):
    """

    apply the dictionary of guess values to parameter, also, save the
    defaults for rollback.

    """
    # only guess thawed parameters!
    if par.frozen:
        return

    if limits and values:
        default_val = par.val
        par.set(param_limits['val'], param_limits['min'], param_limits['max'],
                default_min=par.min, default_max=par.max)

        # set original value as default outside the property interface
        par._default_val = default_val

        # set guessed flag to enable user-defined limit reset
        par._guessed = True

    elif values:
        default_val = par.val
        par.set(param_limits['val'])

        # set original value as default outside the property interface
        par._default_val = default_val

    else:
        par.set(min=param_limits['min'], max=param_limits['max'],
                default_min=par.min, default_max=par.max)

        # set guessed flag to enable user-defined limit reset
        par._guessed = True


def get_amplitude_position(arr, mean=False):
    """
    Guess model amplitude and position of array
    returns (xval, xmin, xmax, xpos)

    """

    xpos = xmin = xmax = xval = 0
    amax = arr.max()
    amin = arr.min()
    if((amax > 0.0 and amin >= 0.0) or
       (amax > 0.0 and amin < 0.0 and (abs(amin) <= amax))):
        xpos = arr.argmax()
        if mean:
            xpos = numpy.where(arr == amax)

        xmax = amax * _guess_ampl_scale
        xmin = amax / _guess_ampl_scale
        xval = amax

    elif((amax > 0.0 and amin < 0.0 and abs(amin) > amax) or
         (amax == 0.0 and amin < 0.0) or (amax < 0.0)):
        xpos = arr.argmin()
        if mean:
            xpos = numpy.where(arr == amin)

        xmax = amin / _guess_ampl_scale
        xmin = amin * _guess_ampl_scale
        xval = amin
    elif (amax == 0.0 and amin == 0.0):
        xpos = arr.argmax()
        if mean:
            xpos = numpy.where(arr == amax)

        xmax = 100.0 / _guess_ampl_scale
        xmin = 0.0
        xval = 0.0

    return (xval, xmin, xmax, xpos)


def guess_amplitude(y, x, xhi=None):
    """
    Guess model parameter amplitude (val, min, max)

    """

    val, ymin, ymax, pos = get_amplitude_position(y)

    amin, amax = None, None
    if ymin != 0.0 or ymax != 0.0:
        amin = ymin
        amax = ymax

    if xhi is not None:
        binsize = numpy.abs(xhi[pos] - x[pos])
        if amin is not None:
            amin /= binsize
        if amax is not None:
            amax /= binsize
        val /= binsize

    return {'val': val, 'min': amin, 'max': amax}


def guess_amplitude_at_ref(r, y, x, xhi=None):
    """
    Guess model parameter amplitude (val, min, max)

    only valid for beta1d, bpl1d, powlaw1d

    """

    t = 1.0
    if x[1] > x[0] and r < x[0]:
        t = numpy.abs(y[0] + y[1]) / 2.0
    elif x[1] > x[0] and r > x[-1]:
        t = numpy.abs(y[-1] + y[-2]) / 2.0
    elif x[1] < x[0] and r > x[0]:
        t = numpy.abs(y[0] + y[1]) / 2.0
    elif x[1] < x[0] and r < x[-1]:
        t = numpy.abs(y[-1] + y[-2]) / 2.0
    else:
        for i in xrange(len(x) - 1):
            if ((r >= x[i] and r < x[i + 1]) or (r >= x[i + 1] and r < x[i])):
                t = numpy.abs(y[i] + y[i + 1]) / 2.0
                break

    if t == 0.0:
        totband = 0.0
        dv = 0.0
        i = 1
        for j in xrange(len(x) - 1):
            dv = x[i] - x[i - 1]
            t += y[i] * dv
            totband += dv
            i += 1
        t /= totband

    return {'val': t,
            'min': t / _guess_ampl_scale, 'max': t * _guess_ampl_scale}


def guess_amplitude2d(y, x0lo, x1lo, x0hi=None, x1hi=None):
    """
    Guess 2D model parameter amplitude (val, min, max)

    """

    limits = guess_amplitude(y, x0lo)

    # if (x0hi is not None and x1hi is not None):
    #     binsize = numpy.abs((x0hi[0]-x0lo[0])*(x1hi[0]-x1lo[0]))
    #     if limits['min'] is not None:
    #         limits['min'] /= binsize
    #     if limits['max'] is not None:
    #         limits['max'] /= binsize
    #     limits['val'] /= binsize

    return limits


def guess_reference(pmin, pmax, x, xhi=None):
    """
    Guess model parameter reference (val, min, max)

    """

    xmin = x.min()
    xmax = x.max()

    if xmin >= 1:
        pmin = 1
    if xmax <= 1:
        pmax = 1

    val = 0.0
    if xmin < 1.0 and xmax > 1.0:
        val = 1.0
    else:
        refval = numpy.floor((xmin + xmax) / 2.0)
        if refval < pmin or refval > pmax:
            refval = (xmin + xmax) / 2.0
        val = refval

    return {'val': val, 'min': None, 'max': None}


def get_position(y, x, xhi=None):
    """
    Get 1D model parameter positions pos (val, min, max)

    """
    pos = get_amplitude_position(y, mean=True)
    xpos = pos[3]

    val = numpy.mean(x[xpos])
    xmin = x.min()
    xmax = x.max()
    if xhi is not None:
        xmax = xhi.max()

    return {'val': val, 'min': xmin, 'max': xmax}


def guess_position(y, x0lo, x1lo, x0hi=None, x1hi=None):
    """
    Guess 2D model parameter positions xpos, ypos ({val0, min0, max0},
                                                   {val1, min1, max1})

    """

    # pos = int(y.argmax())
    # return the average of location of brightest pixels
    pos = numpy.where(y == y.max())

    x0, x1 = x0lo, x1lo
    r1 = {'val': numpy.mean(x0[pos]), 'min': x0.min()}
    r2 = {'val': numpy.mean(x1[pos]), 'min': x1.min()}
    if x0hi is None and x1hi is None:
        r1['max'] = x0.max()
        r2['max'] = x1.max()
    else:
        r1['max'] = x0hi.max()
        r2['max'] = x1hi.max()

    return (r1, r2)


def guess_bounds(x, xhi=True):
    """
    Guess model parameters xlo, xhi (val, min, max)

    """
    xmin = x.min()
    xmax = x.max()
    lo = xmin + (xmax - xmin) / 2.0
    if xhi:
        lo = xmin + (xmax - xmin) / 3.0
        hi = xmin + 2.0 * (xmax - xmin) / 3.0
        return ({'val': lo, 'min': xmin, 'max': xmax},
                {'val': hi, 'min': xmin, 'max': xmax})

    return {'val': lo, 'min': xmin, 'max': xmax}


def guess_radius(x0lo, x1lo, x0hi=None, x1hi=None):
    """
    Guess 2D model parameter radius (val, min, max)

    """
    # TODO: the following was the original code, but
    #   a) x1 isn't used
    #   b) there's no difference between the two branches
    # So, x0hi/x1hi are curently unused.
    #
    # if x0hi is None and x1hi is None:
    #     x0, x1 = x0lo, x1lo
    # else:
    #     x0, x1 = x0lo, x1lo
    x0 = x0lo

    delta = numpy.apply_along_axis(numpy.diff, 0, x0)[0]
    rad = numpy.abs(10 * delta)

    return {'val': rad,
            'min': rad / _guess_ampl_scale, 'max': rad * _guess_ampl_scale}


def split_array(arr, m):
    """Split array ``arr`` into ``m`` roughly equal chunks
    >>> split_array(range(27), 6)
    [[0, 1, 2, 3, 4],
     [5, 6, 7, 8],
     [9, 10, 11, 12, 13],
     [14, 15, 16, 17],
     [18, 19, 20, 21, 22],
     [23, 24, 25, 26]]

    >>> import numpy
    >>> split_array(numpy.arange(25), 6)
    [array([0, 1, 2, 3]),
     array([4, 5, 6, 7]),
     array([ 8,  9, 10, 11, 12]),
     array([13, 14, 15, 16]),
     array([17, 18, 19, 20]),
     array([21, 22, 23, 24])]

    >>> split_array(numpy.arange(30).reshape(5,-1), 3)
    [array([[ 0,  1,  2,  3,  4,  5],
           [ 6,  7,  8,  9, 10, 11]]),
    array([[12, 13, 14, 15, 16, 17]]),
    array([[18, 19, 20, 21, 22, 23],
          [24, 25, 26, 27, 28, 29]])]

    Author: Tom Aldcroft
      split_array() - originated from Python users working group
    """
    n = len(arr)
    idx = [int(round(i * n / float(m))) for i in range(m + 1)]
    return [arr[idx[i]:idx[i + 1]] for i in range(m)]


def worker(f, ii, chunk, out_q, err_q, lock):

    try:
        vals = map(f, chunk)
    except Exception as e:
        err_q.put(e)
        return

    # output the result and task ID to output queue
    out_q.put((ii, vals))


def run_tasks(procs, err_q, out_q, num):

    die = (lambda vals: [val.terminate() for val in vals
                         if val.exitcode is None])

    try:
        for proc in procs:
            proc.start()

        for proc in procs:
            proc.join()

    except KeyboardInterrupt as e:
        # kill all slave processes on ctrl-C
        die(procs)
        raise e

    if not err_q.empty():
        die(procs)
        raise err_q.get()

    results = [None] * num
    while not out_q.empty():
        idx, result = out_q.get()
        results[idx] = result

    # return list(numpy.concatenate(results))
    # Remove extra dimension added by split
    vals = []
    [vals.extend(result) for result in results]
    return vals


def parallel_map(function, sequence, numcores=None):
    """
    A parallelized version of the native Python map function that
    utilizes the Python multiprocessing module to divide and
    conquer sequence.

    parallel_map does not yet support multiple argument sequences.

    :param function: callable function that accepts argument from iterable
    :param sequence: iterable sequence
    :param numcores: number of cores to use
    """
    if not callable(function):
        raise TypeError("input function '%s' is not callable" %
                        repr(function))

    if not numpy.iterable(sequence):
        raise TypeError("input '%s' is not iterable" %
                        repr(sequence))

    size = len(sequence)

    if not _multi or size == 1 or (numcores is not None and numcores < 2):
        return map(function, sequence)

    if numcores is None:
        numcores = _ncpus

    # Returns a started SyncManager object which can be used for sharing
    # objects between processes. The returned manager object corresponds
    # to a spawned child process and has methods which will create shared
    # objects and return corresponding proxies.
    manager = multiprocessing.Manager()

    # Create FIFO queue and lock shared objects and return proxies to them.
    # The managers handles a server process that manages shared objects that
    # each slave process has access to.  Bottom line -- thread-safe.
    out_q = manager.Queue()
    err_q = manager.Queue()
    lock = manager.Lock()

    # if sequence is less than numcores, only use len sequence number of
    # processes
    if size < numcores:
        numcores = size

    # group sequence into numcores-worth of chunks
    sequence = split_array(sequence, numcores)

    procs = [multiprocessing.Process(target=worker,
                                     args=(function, ii, chunk, out_q, err_q, lock))
             for ii, chunk in enumerate(sequence)]

    return run_tasks(procs, err_q, out_q, numcores)


################################# Neville2d ###################################


def neville2d(xinterp, yinterp, x, y, fval):
    nrow = fval.shape[0]
    # ncol = fval.shape[1]
    tmp = numpy.zeros(nrow)
    for row in xrange(nrow):
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
            return numpy.Inf
        return (self.func(x + h) - self.func(x - h)) / (2.0 * h)


class NumDerivFowardPartial(NumDeriv):

    def __init__(self, func, fval0):
        NumDeriv.__init__(self, func, fval0)

    def __call__(self, x, h, *args):

        if 0.0 == h:
            h = pow(numpy.float_(numpy.finfo(numpy.float32)).eps, 1.0 / 3.0)

        ith = args[0]
        jth = args[1]

        ei = numpy.zeros(len(x), float)
        ej = numpy.zeros(len(x), float)

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
            h = pow(numpy.float_(numpy.finfo(numpy.float32)).eps, 1.0 / 3.0)

        ith = args[0]
        jth = args[1]

        ei = numpy.zeros(len(x), float)

        if ith == jth:

            delta = h * abs(x[ith])
            if 0.0 == delta:
                delta = h
            ei[ith] = delta

            fval = - 2.0 * self.fval_0
            fval += self.func(x + ei) + self.func(x - ei)
            fval /= delta * delta
            return fval

        else:

            ej = numpy.zeros(len(x), float)

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

    def __init__(self, sequence, verbose=False):
        self.sequence = sequence
        self.verbose = verbose

    def __call__(self, x, t, tol, maxiter, h, *args):

        richardson = numpy.zeros((maxiter, maxiter), dtype=numpy.float_)
        richardson[0, 0] = self.sequence(x, h, *args)

        t_sqr = t * t
        for ii in xrange(1, maxiter):
            h /= t
            richardson[ii, 0] = self.sequence(x, h, *args)
            ii_1 = ii - 1
            for jj in xrange(1, ii + 1):
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
    Hessian = numpy.zeros((npar, npar), dtype=numpy.float_)
    for ii in xrange(npar):
        for jj in xrange(ii + 1):
            answer = deriv(par, t, tol, maxiter, h, ii, jj)
            Hessian[ii, jj] = answer / 2.0
            Hessian[jj, ii] = Hessian[ii, jj]
    return Hessian, num_dif.nfev[0]


def print_low_triangle(matrix, num):
    # print matrix
    for ii in xrange(num):
        print(matrix[ii, 0], end=' ')
        for jj in xrange(1, ii + 1):
            print(matrix[ii, jj], end=' ')
        print()


def symmetric_to_low_triangle(matrix, num):
    low_triangle = []
    for ii in xrange(num):
        for jj in xrange(ii + 1):
            low_triangle.append(matrix[ii, jj])
    # print_low_triangle( matrix, num )
    # print low_triangle
    return low_triangle


############################### Root of all evil ##############################


def printf(format, *args):
    """Format args with the first argument as format string, and write.
    Return the last arg, or format itself if there are no args."""
    sys.stdout.write(str(format) % args)
    return if_(args, args[-1], format)


def func_counter(func):
    '''A function wrapper to count the number of times being called'''
    nfev = [0]

    def func_counter_wrapper(x, *args):
        nfev[0] += 1
        return func(x, *args)

    return nfev, func_counter_wrapper


def func_counter_history(func, history):
    '''A function wrapper to count the number of times being called'''
    nfev = [0]

    def func_counter_history_wrapper(x, *args):
        nfev[0] += 1
        y = func(x, *args)
        history[0].append(x)
        history[1].append(y)
        return y

    return nfev, func_counter_history_wrapper


def is_in(arg, seq):
    for x in seq:
        if arg == x:
            return True
    return False


def is_iterable(arg):
    return isinstance(arg, list) or isinstance(arg, tuple) \
        or isinstance(arg, numpy.ndarray) or numpy.iterable(arg)


def is_sequence(start, mid, end):
    return (start < mid) and (mid < end)


def Knuth_close(x, y, tol, myop=operator.__or__):
    """Check whether two floating-point numbers are close together.

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
    import sys
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
    if not numpy.iterable(arg):
        return arg
    str = '(%e, %e)' % (arg[0], arg[1])
    return str


def mysgn(arg):
    if arg == 0.0:
        return 0
    elif arg < 0.0:
        return -1
    else:
        return 1


# Is this ever used? It is checked for as an exception, so really should be
# derived from an exception, but is never thrown, as far as I can see.
class OutOfBoundErr:
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

            else:

                #
                # 0 * x^2 + 0 * x + c = 0
                #
                # a == 0, b == 0, so if c == 0 then all numbers work so
                # returning nan is not right. However if c != 0 then no
                # roots exist.
                #
                return [None, None]

        elif 0.0 == b:

            #
            # a * x^2 + 0 * x + c = 0
            #
            if 0.0 == c:

                # a * x^2 + 0 * x + 0 = 0
                return [0.0, 0.0]
            else:

                # a * x^2 + 0 * x + c = 0
                if mysgn(a) == mysgn(c):
                    return [None, None]
                answer = numpy.sqrt(c / a)
                return [-answer, answer]

        elif 0.0 == c:

            #
            # a * x^2 + b * x + 0 = 0
            #
            return [0.0, - b / a]

        else:

            discriminant = b * b - 4.0 * a * c
            # TODO: is this needed?
            debug("disc={}".format(discriminant))
            sqrt_disc = numpy.sqrt(discriminant)
            t = - (b + mysgn(b) * sqrt_disc) / 2.0
            return [c / t, t / a]


def bisection(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=48, tol=1.0e-6):

    history = [[], []]
    nfev, myfcn = func_counter_history(fcn, history)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], nfev[0]]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], nfev[0]]

        if mysgn(fa) == mysgn(fb):
            # TODO: is this a useful message for the user?
            warning(__name__ + ': ' + fcn.__name__ + ' fa * fb < 0 is not met')
            return [[None, None], [[None, None], [None, None]], nfev[0]]

        while nfev[0] < maxfev:

            if abs(fa) > tol and abs(fb) > tol:

                xc = (xa + xb) / 2.0
                fc = myfcn(xc, *args)

                if abs(xa - xb) < min(tol * abs(xb), tol / 10.0):
                    return [[xc, fc], [[xa, fa], [xb, fb]], nfev[0]]

                if mysgn(fa) != mysgn(fc):
                    xb, fb = xc, fc
                else:
                    xa, fa = xc, fc

            else:
                if abs(fa) <= tol:
                    return [[xa, fa], [[xa, fa], [xb, fb]], nfev[0]]
                else:
                    return [[xb, fb], [[xa, fa], [xb, fb]], nfev[0]]

        xc = (xa + xb) / 2.0
        fc = myfcn(xc, *args)
        return [[xc, fc], [[xa, fa], [xb, fb]], nfev[0]]

    except OutOfBoundErr:
        return [[None, None], [[xa, fa], [xb, fb]], nfev[0]]


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

    xc_xb = xc - xb
    fc_fb = fc - fb
    A = fc_fb / xc_xb
    fb_fa = fb - fa
    xb_xa = xb - xa
    xc_xa = xc - xa
    B = (A - fb_fa / xb_xa) / xc_xa
    C = A + B * xc_xb
    return [B, C]


def demuller(fcn, xa, xb, xc, fa=None, fb=None, fc=None, args=(),
             maxfev=32, tol=1.0e-6):
    """A root-finding algorithm using Muller's method [1]_.

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

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Muller%27s_method

    """

    def is_nan(arg):
        if arg != arg:
            return True
        if arg is numpy.nan:
            return True
        return numpy.isnan(arg)

    history = [[], []]
    nfev, myfcn = func_counter_history(fcn, history)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], nfev[0]]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], nfev[0]]

        if fc is None:
            fc = myfcn(xc, *args)
        if abs(fc) <= tol:
            return [[xc, fc], [[xc, fc], [xc, fc]], nfev[0]]

        while nfev[0] < maxfev:

            [B, C] = transformed_quad_coef([xa, xb, xc], [fa, fb, fc])

            discriminant = max(C * C - 4.0 * fc * B, 0.0)

            if is_nan(B) or is_nan(C) or \
                    0.0 == C + mysgn(C) * numpy.sqrt(discriminant):
                return [[None, None], [[None, None], [None, None]], nfev[0]]

            xd = xc - 2.0 * fc / (C + mysgn(C) * numpy.sqrt(discriminant))

            fd = myfcn(xd, *args)
            # print 'fd(%e)=%e' % (xd, fd)
            if abs(fd) <= tol:
                return [[xd, fd], [[None, None], [None, None]], nfev[0]]

            xa = xb
            fa = fb
            xb = xc
            fb = fc
            xc = xd
            fc = fd

        # print 'demuller(): maxfev exceeded'
        return [[xd, fd], [[None, None], [None, None]], nfev[0]]

    except ZeroDivisionError:

        # print 'demuller(): fixme ZeroDivisionError'
        # for x, y in izip( history[0], history[1] ):
        #     print 'f(%e)=%e' % ( x, y )
        return [[xd, fd], [[None, None], [None, None]], nfev[0]]


def new_muller(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32, tol=1.e-6):

    # This function does not appear to be used
    def regula_falsi(x0, x1, f0, f1):
        if f0 < f1:
            xl, fl = x0, f0
            xh, fh = x1, f1
        else:
            xl, fl = x1, f1
            xh, fh = x0, f0
        x = xl + (xh - xl) * fl / (fl - fh)
        if is_sequence(x0, x, x1):
            return x
        else:
            return (x0 + x1) / 2.0

    history = [[], []]
    nfev, myfcn = func_counter_history(fcn, history)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], nfev[0]]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], nfev[0]]

        if mysgn(fa) == mysgn(fb):
            # TODO: is this a useful message for the user?
            warning(__name__ + ': ' + fcn.__name__ + ' fa * fb < 0 is not met')
            return [[None, None], [[None, None], [None, None]], nfev[0]]

        while nfev[0] < maxfev:

            xc = (xa + xb) / 2.0
            fc = myfcn(xc, *args)
            if abs(fc) <= tol:
                return [[xc, fc], [[xa, fa], [xb, fb]], nfev[0]]

            tran = transformed_quad_coef([xa, xb, xc], [fa, fb, fc])
            B = tran[0]
            C = tran[1]

            discriminant = max(C * C - 4.0 * fc * B, 0.0)

            xd = xc - 2.0 * fc / (C + mysgn(C) * numpy.sqrt(discriminant))

            fd = myfcn(xd, *args)
            # print 'fd(%e)=%e' % (xd, fd)
            if abs(fd) <= tol:
                return [[xd, fd], [[xa, fa], [xb, fb]], nfev[0]]

            if mysgn(fa) != mysgn(fc):
                xb, fb = xc, fc
                continue

            if mysgn(fd) != mysgn(fc) and xc < xd:
                xa, fa = xc, fc
                xb, fb = xd, fd
                continue

            if mysgn(fb) != mysgn(fd):
                xa, fa = xd, fd
                continue

            if mysgn(fa) != mysgn(fd):
                xb, fb = xd, fd
                continue

            if mysgn(fc) != mysgn(fd) and xd < xc:
                xa, fa = xd, fd
                xb, fb = xc, fc
                continue

            if mysgn(fc) != mysgn(fd):
                xa, fa = xc, fc
                continue

        # print 'new_muller(): maxfev exceeded'
        return [[xd, fd], [[xa, fa], [xb, fb]], nfev[0]]

    except (ZeroDivisionError, OutOfBoundErr):

        # print 'new_muller(): fixme ZeroDivisionError'
        # for x, y in izip( history[0], history[1] ):
        #     print 'f(%e)=%e' % ( x, y )
        return [[xd, fd], [[xa, fa], [xb, fb]], nfev[0]]

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

    history = [[], []]
    nfev, myfcn = func_counter_history(fcn, history)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
        if abs(fa) <= tol:
            return [[xa, fa], [[xa, fa], [xa, fa]], nfev[0]]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xb, fb], [xb, fb]], nfev[0]]

        if mysgn(fa) == mysgn(fb):
            # TODO: is this a useful message for the user?
            warning(__name__ + ': ' + fcn.__name__ + ' fa * fb < 0 is not met')
            return [[None, None], [[None, None], [None, None]], nfev[0]]

        xc = (xa + xb) / 2.0
        fc = myfcn(xc, *args)
        # print 'MullerBound() fc(%.14e)=%.14e' % (xc,fc)
        if abs(fc) <= tol:
            return [[xc, fc], [[xc, fc], [xc, fc]], nfev[0]]

        xbest, fbest = xa, fa
        if abs(fb) < abs(fa):
            xbest, fbest = xb, fb
        if abs(fc) < abs(fbest):
            xbest, fbest = xc, fc

        oldx = 1.0e128
        while nfev[0] < maxfev:

            tran = transformed_quad_coef([xa, xb, xc], [fa, fb, fc])
            B = tran[0]
            C = tran[1]
            discriminant = max(C * C - 4.0 * fc * B, 0.0)
            den = mysgn(C) * numpy.sqrt(discriminant)
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
            # print 'MullerBound() y(%.14e)=%.14e' % (x,y)
            if abs(y) < abs(fbest):
                xbest, fbest = x, y
            tolerance = min(tol * abs(x), tol)
            if abs(y) <= tol or abs(x - oldx) <= tolerance:
                return [[x, y], [[xa, fa], [xb, fb]], nfev[0]]

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
                # print 'MullerBound() fmid(%.14e)=%.14e' % (xmid,fmid)
                if abs(fmid) <= tol:
                    return [[xmid, fmid], [[xa, fa], [xb, fb]], nfev[0]]
                if mysgn(fa) + mysgn(fmid) == 0:
                    xb = xmid
                    fb = fmid
                else:
                    xa = xmid
                    fa = fmid
                xc = (xa + xb) / 2.0
                fc = myfcn(xc, *args)
                if abs(fc) < abs(fbest):
                    xbest, fbest = xc, fc
                # print 'MullerBound() fc(%.14e)=%.14e' % (xc,fc)
                if abs(fc) <= tol:
                    return [[xc, fc], [[xa, fa], [xb, fb]], nfev[0]]
                oldx = 1.0e128

        #
        # maxfev has exceeded, return the minimum so far
        #
        return [[xbest, fbest], [[xa, fa], [xb, fb]], nfev[0]]

    #
    # Something drastic has happened
    #
    except (ZeroDivisionError, OutOfBoundErr):

        return [[xbest, fbest], [[xa, fa], [xb, fb]], nfev[0]]


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

    """

    history = [[], []]
    nfev, myfcn = func_counter_history(fcn, history)

    try:

        if fa is None:
            fa = myfcn(xa, *args)
            if abs(fa) <= tol:
                return [[xa, fa], [[xa, fa], [xb, fb]], nfev[0]]

        if fb is None:
            fb = myfcn(xb, *args)
        if abs(fb) <= tol:
            return [[xb, fb], [[xa, fa], [xb, fb]], nfev[0]]

        if mysgn(fa) == mysgn(fb):
            # TODO: is this a useful message for the user?
            warning(__name__ + ': ' + fcn.__name__ + ' fa * fb < 0 is not met')
            return [[None, None], [[None, None], [None, None]], nfev[0]]

        xc = xa
        fc = fa
        DBL_EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)
        while nfev[0] < maxfev:

            prev_step = xb - xa

            if abs(fc) < abs(fb):
                xa, fa = xb, fb
                xb, fb = xc, fc
                xc, fc = xa, fa

            tol_act = 2.0 * DBL_EPSILON * abs(xb) + tol / 2.0
            new_step = (xc - xb) / 2.0

            if abs(fb) <= tol:
                return [[xb, fb], [[xa, fa], [xb, fb]], nfev[0]]

            if abs(new_step) <= tol_act:
                if mysgn(fb) != mysgn(fa):
                    tmp = apache_muller(fcn, xa, xb, fa, fb, args=args,
                                        maxfev=maxfev - nfev[0], tol=tol)
                    tmp[-1] += nfev[0]
                    return tmp
                elif mysgn(fb) != mysgn(fc):
                    tmp = apache_muller(fcn, xb, xc, fb, fc, args=args,
                                        maxfev=maxfev - nfev[0], tol=tol)
                    tmp[-1] += nfev[0]
                    return tmp
                else:
                    return [[xb, fb], [[xa, fa], [xb, fb]], nfev[0]]

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
            # print 'fa(%f)=%f\tfb(%f)=%f\tfc(%f)=%f' % (xa,fa,xb,fb,xc,fc)
            if fb > 0 and fc > 0 or fb < 0 and fc < 0:
                xc = xa
                fc = fa

        return [[xb, fb], [[xa, fa], [xc, fc]], nfev[0]]

    except (ZeroDivisionError, OutOfBoundErr):
        return [[xb, fb], [[xa, fa], [xc, fc]], nfev[0]]


def get_valid_args(func):
    valid_args = func.__code__.co_varnames[:func.__code__.co_argcount]
    # number of keyword arguments
    kwargs_length = len(func.__defaults__) if func.__defaults__ else 0
    # because kwargs are last
    valid_kwargs = valid_args[-kwargs_length:] if kwargs_length else []
    return valid_kwargs


def public(f):
    """Use a decorator to avoid retyping function/class names.

    * Based on an idea by Duncan Booth:
    http://groups.google.com/group/comp.lang.python/msg/11cbb03e09611b8a
    * Improved via a suggestion by Dave Angel:
    http://groups.google.com/group/comp.lang.python/msg/3d400fb22d8a42e1
    """
    _all = sys.modules[f.__module__].__dict__.setdefault('__all__', [])
    if f.__name__ not in _all:  # Prevent duplicates if run from an IDE.
        _all.append(f.__name__)
    return f
