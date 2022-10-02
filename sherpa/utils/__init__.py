#
#  Copyright (C) 2007, 2015, 2016, 2018, 2019, 2020, 2021, 2022
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

"""
Objects and utilities used by multiple Sherpa subpackages
"""

from configparser import ConfigParser, NoSectionError
import inspect
import logging
import operator
import os
import pydoc
import string
import sys
from types import FunctionType as function
from types import MethodType as instancemethod

import numpy
import numpy.random
import numpy.fft

from sherpa import get_config
# Note: _utils.gsl_fcmp and _utils.ndtri are not exported from
#       this module; is this intentional?
from sherpa.utils._utils import hist2d
from sherpa.utils import _utils, _psf
from sherpa.utils.err import IOErr


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

    multiprocessing_start_method = config.get('multiprocessing', 'multiprocessing_start_method', fallback='fork')

    if multiprocessing_start_method not in ('fork', 'spawn', 'default'):
        raise ValueError('multiprocessing_start method must be one of "fork", "spawn", or "default"')

    if multiprocessing_start_method != 'default':
        multiprocessing.set_start_method(multiprocessing_start_method, force=True)

    _multi = True

    if _ncpus is None:
        _ncpus = multiprocessing.cpu_count()
except Exception as e:
    warning("parallel processing is unavailable,\n" +
            "multiprocessing module failed with \n'%s'" % str(e))
    _ncpus = 1
    _multi = False

del _ncpu_val, config, get_config, ConfigParser, NoSectionError


__all__ = ('NoNewAttributesAfterInit', 'SherpaFloat',
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
           'new_muller', 'normalize', 'numpy_convolve',
           'pad_bounding_box', 'parallel_map', 'parallel_map_funcs',
           'param_apply_limits', 'parse_expr', 'poisson_noise',
           'print_fields', 'rebin',
           'sao_arange', 'sao_fcmp', 'send_to_pager',
           'set_origin', 'sum_intervals', 'zeroin',
           'multinormal_pdf', 'multit_pdf', 'get_error_estimates', 'quantile')


_guess_ampl_scale = 1.e+3
"""The scaling applied to a value to create its range.

The minimum and maximum values for a range are calculated by
dividing and multiplying the value by ``_guess_ampl_scale``.
"""

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


class NoNewAttributesAfterInit():
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
