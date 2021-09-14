#
#  Copyright (C) 2007, 2016, 2018, 2019, 2020, 2021
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

"""Optimizing functions.

These functions take a callback, the current set of parameters, the
minimum and maximum parameter ranges, along with optional arguments,
and return a tuple containing

    status, parameters, statistic, message, dict

where ``status`` is a boolean indicating whether the optimisation
succeeded or not, parameters is the list of parameter values at the
best-fit location, the statistic value at this location, a string
message - when ``status`` is `False` this will give information on the
failure - and a dictionary which depends on the optimiser.

The callback should return the current statistic value and an array
of the statistic value per bin.

Notes
-----

Each optimizer has certain classes of problem where it is more, or
less, successful. For instance, the `neldermead` function should
only be used with chi-square based statistics.

Examples
--------

Fit a constant model to the array of values in ``y``, using a
least-square statistic:

>>> y = np.asarray([3, 2, 7])
>>> def cb(pars):
...     'Least-squares statistic value from fitting a constant model to y'
...     dy = y - pars[0]
...     dy *= dy
...     return (dy.sum(), dy)
...

This can be evaluated using the `neldermead` optimiser, starting at a
model value of 1 and bounded to the range 0 to 10000:

>>> res = neldermead(cb, [1], [0], [1e4])
>>> print(res)
(True, array([4.]), 14.0, 'Optimization terminated successfully', {'info': True, 'nfev': 98})
>>> print(f"Best-fit value: {res[1][0]}")
Best-fit value: 4.0

"""

import random
import numpy

from . import _saoopt
from sherpa.optmethods.ncoresde import ncoresDifEvo
from sherpa.optmethods.ncoresnm import ncoresNelderMead

from sherpa.utils import parallel_map, func_counter
from sherpa.utils._utils import sao_fcmp

__all__ = ('difevo', 'difevo_lm', 'difevo_nm', 'grid_search', 'lmdif',
           'minim', 'montecarlo', 'neldermead')


#
# Use FLT_EPSILON as default tolerance
#
EPSILON = numpy.float_(numpy.finfo(numpy.float32).eps)

#
# Maximum callback function value, used to indicate that the optimizer
# has exceeded parameter boundaries.  All the optimizers expect double
# precision arguments, so we use numpy.float_ instead of SherpaFloat.
#
FUNC_MAX = numpy.finfo(numpy.float_).max


def _check_args(x0, xmin, xmax):
    x = numpy.array(x0, numpy.float_)  # Make a copy
    xmin = numpy.asarray(xmin, numpy.float_)
    xmax = numpy.asarray(xmax, numpy.float_)

    if (x.shape != xmin.shape) or (x.shape != xmax.shape):
        raise TypeError('input array sizes do not match')

    _move_within_limits(x, xmin, xmax)

    return x, xmin, xmax


def _get_saofit_msg(maxfev, ierr):
    key = {
        0: (True, 'successful termination'),
        1: (False, 'improper input parameters'),
        2: (False, 'initial parameter value is out of bounds'),
        3: (False,
            ('number of function evaluations has exceeded maxfev=%d' %
             maxfev))
        }
    return key.get(ierr, (False, 'unknown status flag (%d)' % ierr))


def _move_within_limits(x, xmin, xmax):
    below = numpy.flatnonzero(x < xmin)
    if below.size > 0:
        x[below] = xmin[below]

    above = numpy.flatnonzero(x > xmax)
    if above.size > 0:
        x[above] = xmax[above]


def _my_is_nan(x):
    fubar = list(filter(lambda xx: xx != xx or xx is numpy.nan or numpy.isnan(xx) and numpy.isfinite(xx), x))
    return len(fubar) > 0


def _narrow_limits(myrange, xxx, debug):

    def double_check_limits(myx, myxmin, myxmax):
        for my_l, my_x, my_h in zip(myxmin, myx, myxmax):
            if my_x < my_l:
                print('x = ', my_x, ' is < lower limit = ', my_l)
            if my_x > my_h:
                print('x = ', my_x, ' is > upper limit = ', my_h)

    def raise_min_limit(range, xmin, x, debug=False):
        myxmin = numpy.asarray(list(map(lambda xx: xx - range * numpy.abs(xx), x)), numpy.float_)
        if debug:
            print()
            print('raise_min_limit: myxmin=%s' % myxmin)
            print('raise_min_limit: x=%s' % x)
        below = numpy.flatnonzero(myxmin < xmin)
        if below.size > 0:
            myxmin[below] = xmin[below]
        if debug:
            print('raise_min_limit: myxmin=%s' % myxmin)
            print('raise_min_limit: x=%s' % x)
            print()
        return myxmin

    def lower_max_limit(range, x, xmax, debug=False):
        myxmax = numpy.asarray(list(map(lambda xx: xx + range * numpy.abs(xx), x)), numpy.float_)
        if debug:
            print()
            print('lower_max_limit: x=%s' % x)
            print('lower_max_limit: myxmax=%s' % myxmax)
        above = numpy.flatnonzero(myxmax > xmax)
        if above.size > 0:
            myxmax[above] = xmax[above]
        if debug:
            print('lower_max_limit: x=%s' % x)
            print('lower_max_limit: myxmax=%s' % myxmax)
            print()
        return myxmax

    x = xxx[0]
    xmin = xxx[1]
    xmax = xxx[2]

    if debug:
        print('narrow_limits: xmin=%s' % xmin)
        print('narrow_limits: x=%s' % x)
        print('narrow_limits: xmax=%s' % xmax)
    myxmin = raise_min_limit(myrange, xmin, x, debug=False)
    myxmax = lower_max_limit(myrange, x, xmax, debug=False)

    if debug:
        print('range = %d' % myrange)
        print('narrow_limits: myxmin=%s' % myxmin)
        print('narrow_limits: x=%s' % x)
        print('narrow_limits: myxmax=%s\n' % myxmax)

    double_check_limits(x, myxmin, myxmax)

    return myxmin, myxmax


def _par_at_boundary(low, val, high, tol):
    for par_min, par_val, par_max in zip(low, val, high):
        if sao_fcmp(par_val, par_min, tol) == 0:
            return True
        if sao_fcmp(par_val, par_max, tol) == 0:
            return True
    return False


def _outside_limits(x, xmin, xmax):
    return (numpy.any(x < xmin) or numpy.any(x > xmax))


def _same_par(a, b):
    b = numpy.array(b, numpy.float_)
    same = numpy.flatnonzero(a < b)
    if same.size == 0:
        return 1
    return 0


def _set_limits(x, xmin, xmax):
    below = numpy.nonzero(x < xmin)
    if below.size > 0:
        return 1

    above = numpy.nonzero(x > xmax)
    if above.size > 0:
        return 1

    return 0


def difevo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
           seed=2005815, population_size=None, xprob=0.9,
           weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max(0.1, xprob)
    xprob = min(xprob, 1.0)

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max(0.1, weighting_factor)
    weighting_factor = min(weighting_factor, 1.0)

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * x.size

    de = _saoopt.difevo(verbose, maxfev, seed, population_size, ftol, xprob,
                        weighting_factor, xmin, xmax, x, fcn)
    fval = de[1]
    nfev = de[2]
    ierr = de[3]

    if verbose:
        print('difevo: f%s=%e in %d nfev' % (x, fval, nfev))

    status, msg = _get_saofit_msg(maxfev, ierr)
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})

    return rv


def difevo_lm(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
              seed=2005815, population_size=None, xprob=0.9,
              weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max(0.1, xprob)
    xprob = min(xprob, 1.0)

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max(0.1, weighting_factor)
    weighting_factor = min(weighting_factor, 1.0)

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * x.size

    de = _saoopt.lm_difevo(verbose, maxfev, seed, population_size, ftol,
                           xprob, weighting_factor, xmin, xmax,
                           x, fcn, numpy.asanyarray(fcn(x)).size)
    fval = de[1]
    nfev = de[2]
    ierr = de[3]

    status, msg = _get_saofit_msg(maxfev, ierr)
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})

    return rv


def difevo_nm(fcn, x0, xmin, xmax, ftol, maxfev, verbose, seed,
              population_size, xprob, weighting_factor):

    def stat_cb0(pars):
        return fcn(pars)[0]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max(0.1, xprob)
    xprob = min(xprob, 1.0)

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max(0.1, weighting_factor)
    weighting_factor = min(weighting_factor, 1.0)

    if population_size is None:
        population_size = max(population_size, 16 * x.size)

    if maxfev is None:
        maxfev = 1024 * population_size

    de = _saoopt.nm_difevo(verbose, maxfev, seed, population_size,
                           ftol, xprob, weighting_factor, xmin, xmax,
                           x, stat_cb0)
    fval = de[1]
    nfev = de[2]
    ierr = de[3]

    if verbose:
        print('difevo_nm: f%s=%e in %d nfev' % (x, fval, nfev))

    status, msg = _get_saofit_msg(maxfev, ierr)
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})

    return rv


def grid_search(fcn, x0, xmin, xmax, num=16, sequence=None, numcores=1,
                maxfev=None, ftol=EPSILON, method=None, verbose=0):
    """Grid Search optimization method.

    This method evaluates the fit statistic for each point in the
    parameter space grid; the best match is the grid point with the
    lowest value of the fit statistic. It is intended for use with
    template models as it is very inefficient for general models.

    Parameters
    ----------
    fcn : function reference
       Returns the current statistic and per-bin statistic value when
       given the model parameters.
    x0, xmin, xmax : sequence of number
       The starting point, minimum, and maximum values for each
       parameter.
    num : int
       The size of the grid for each parameter when `sequence` is
       `None`, so ``npar^num`` fits will be evaluated, where `npar` is
       the number of free parameters. The grid spacing is uniform.
    sequence : sequence of numbers or `None`
       The list through which to evaluate. Leave as `None` to use
       a uniform grid spacing as determined by the `num` attribute.
    numcores : int or `None`
       The number of CPU cores to use. The default is `1` and a
       value of `None` will use all the cores on the machine.
    maxfev : int or `None`
       The `maxfev` attribute if `method` is not `None`.
    ftol : number
       The `ftol` attribute if `method` is not `None`.
    method : str or `None`
       The optimization method to use to refine the best-fit
       location found using the grid search. If `None` then
       this step is not run.
    verbose: int
       The amount of information to print during the fit. The default
       is `0`, which means no output.

    Returns
    -------
    retval : tuple
       A boolean indicating whether the optimization succeeded, the
       best-fit parameter values, the best-fit statistic value, a
       string message indicating the status, and a dictionary
       returning information from the optimizer.

    """

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    npar = len(x)

    def func(pars):
        aaa = fcn(pars)[0]
        if verbose:
            print('f%s=%g' % (pars, aaa))
        return aaa

    def make_sequence(ranges, N):
        list_ranges = list(ranges)
        for ii in range(npar):
            list_ranges[ii] = tuple(list_ranges[ii]) + (complex(N),)
            list_ranges[ii] = slice(*list_ranges[ii])

        grid = numpy.mgrid[list_ranges]
        mynfev = pow(N, npar)
        grid = list(map(numpy.ravel, grid))
        sequence = []
        for index in range(mynfev):
            tmp = []
            for xx in range(npar):
                tmp.append(grid[xx][index])
            sequence.append(tmp)
        return sequence

    def eval_stat_func(xxx):
        return numpy.append(func(xxx), xxx)

    if sequence is None:
        ranges = []
        for index in range(npar):
            ranges.append([xmin[index], xmax[index]])
        sequence = make_sequence(ranges, num)
    else:
        if not numpy.iterable(sequence):
            raise TypeError("sequence option must be iterable")
        else:
            for seq in sequence:
                if npar != len(seq):
                    msg = "%s must be of length %d" % (seq, npar)
                    raise TypeError(msg)

    answer = eval_stat_func(x)
    sequence_results = list(parallel_map(eval_stat_func, sequence, numcores))
    for xresult in sequence_results[1:]:
        if xresult[0] < answer[0]:
            answer = xresult

    fval = answer[0]
    x = answer[1:]
    nfev = len(sequence_results) + 1
    ierr = 0
    status, msg = _get_saofit_msg(ierr, ierr)
    rv = (status, x, fval)
    rv += (msg, {'info': ierr, 'nfev': nfev})

    # TODO: should we just use case-insensitive comparison?
    if method in ['NelderMead', 'neldermead', 'Neldermead', 'nelderMead']:
        # re.search( '^[Nn]elder[Mm]ead', method ):
        nm_result = neldermead(fcn, x, xmin, xmax, ftol=ftol, maxfev=maxfev,
                               verbose=verbose)
        tmp_nm_result = list(nm_result)
        tmp_nm_result_4 = tmp_nm_result[4]
        tmp_nm_result_4['nfev'] += nfev
        rv = tuple(tmp_nm_result)

    if method in ['LevMar', 'levmar', 'Levmar', 'levMar']:
        # re.search( '^[Ll]ev[Mm]ar', method ):
        levmar_result = lmdif(fcn, x, xmin, xmax, ftol=ftol, xtol=ftol,
                              gtol=ftol, maxfev=maxfev, verbose=verbose)
        tmp_levmar_result = list(levmar_result)
        tmp_levmar_result_4 = tmp_levmar_result[4]
        tmp_levmar_result_4['nfev'] += nfev
        rv = tuple(tmp_levmar_result)

    return rv


#
# C-version of minim
#
def minim(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, step=None,
          nloop=1, iquad=1, simp=None, verbose=-1, reflect=True):

    # TODO: rework so do not have two stat_cb0 functions which
    #       are both used
    def stat_cb0(pars):
        return fcn(pars)[0]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if step is None:
        order = 'F' if numpy.isfortran(x) else 'C'
        step = 0.4*numpy.ones(x.shape, numpy.float_, order)
    if simp is None:
        simp = 1.0e-2 * ftol
    if maxfev is None:
        maxfev = 512 * len(x)

    orig_fcn = stat_cb0

    def stat_cb0(x_new):
        if _my_is_nan(x_new) or _outside_limits(x_new, xmin, xmax):
            return FUNC_MAX
        return orig_fcn(x_new)

    init = 0
    x, fval, neval, ifault = _saoopt.minim(reflect, verbose, maxfev, init, \
                                           iquad, simp, ftol, step, \
                                           xmin, xmax, x, stat_cb0)

    key = {
        0: (True, 'successful termination'),
        1: (False,
            'number of function evaluations has exceeded maxfev=%d' % maxfev),
        2: (False, 'information matrix is not +ve semi-definite'),
        3: (False, 'number of parameters is less than 1'),
        4: (False, 'nloop=%d is less than 1' % nloop)
        }
    status, msg = key.get(ifault, (False, 'unknown status flag (%d)' % ifault))

    rv = (status, x, fval, msg, {'info': ifault, 'nfev': neval})
    return rv


#
# Monte Carlo
#
def montecarlo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
               seed=74815, population_size=None, xprob=0.9,
               weighting_factor=0.8, numcores=1):
    """Monte Carlo optimization method.

    This is an implementation of the differential-evolution algorithm
    from Storn and Price (1997) [1]_. A population of fixed size -
    which contains n-dimensional vectors, where n is the number of
    free parameters - is randomly initialized.  At each iteration, a
    new n-dimensional vector is generated by combining vectors from
    the pool of population, the resulting trial vector is selected if
    it lowers the objective function.

    Parameters
    ----------
    fcn : function reference
       Returns the current statistic and per-bin statistic value when
       given the model parameters.
    x0, xmin, xmax : sequence of number
       The starting point, minimum, and maximum values for each
       parameter.
    ftol : number
       The function tolerance to terminate the search for the minimum;
       the default is sqrt(DBL_EPSILON) ~ 1.19209289551e-07, where
       DBL_EPSILON is the smallest number x such that ``1.0 != 1.0 +
       x``.
    maxfev : int or `None`
       The maximum number of function evaluations; the default value
       of `None` means to use ``8192 * n``, where `n` is the number of
       free parameters.
    verbose: int
       The amount of information to print during the fit. The default
       is `0`, which means no output.
    seed : int
       The seed for the random number generator.
    population_size : int or `None`
       The population of potential solutions is allowed to evolve to
       search for the minimum of the fit statistics. The trial
       solution is randomly chosen from a combination from the current
       population, and it is only accepted if it lowers the
       statistics.  A value of `None` means to use a value ``16 * n``,
       where `n` is the number of free parameters.
    xprob : num
       The crossover probability should be within the range [0.5,1.0];
       default value is 0.9. A high value for the crossover
       probability should result in a faster convergence rate;
       conversely, a lower value should make the differential
       evolution method more robust.
    weighting_factor: num
       The weighting factor should be within the range [0.5, 1.0];
       default is 0.8. Differential evolution is more sensitive to the
       weighting_factor then the xprob parameter. A lower value for
       the weighting_factor, coupled with an increase in the
       population_size, gives a more robust search at the cost of
       efficiency.
    numcores : int
       The number of CPU cores to use. The default is `1`.

    References
    ----------

    .. [1] Storn, R. and Price, K. "Differential Evolution: A Simple
           and Efficient Adaptive Scheme for Global Optimization over
           Continuous Spaces." J. Global Optimization 11, 341-359,
           1997.
           http://www.icsi.berkeley.edu/~storn/code.html

    """

    def stat_cb0(pars):
        return fcn(pars)[0]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = max(0.1, xprob)
    xprob = min(xprob, 1.0)

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = max(0.1, weighting_factor)
    weighting_factor = min(weighting_factor, 1.0)

    random.seed(seed)
    if seed is None:
        # pow(2,31) == 2147483648L
        seed = random.randint(0, 2147483648)
    if population_size is None:
        population_size = 12 * x.size

    if maxfev is None:
        maxfev = 8192 * population_size

    def myopt(myfcn, xxx, ftol, maxfev, seed, pop, xprob,
              weight, factor=4.0, debug=False):

        x = xxx[0]
        xmin = xxx[1]
        xmax = xxx[2]
        maxfev_per_iter = 512 * x.size

        def random_start(xmin, xmax):
            xx = []
            for ii in range(len(xmin)):
                xx.append(random.uniform(xmin[ii], xmax[ii]))
            return numpy.asarray(xx)

        ############################# NelderMead #############################
        mymaxfev = min(maxfev_per_iter, maxfev)
        if all(x == 0.0):
            mystep = list(map(lambda fubar: 1.2 + fubar, x))
        else:
            mystep = list(map(lambda fubar: 1.2 * fubar, x))
        if 1 == numcores:
            result = neldermead(myfcn, x, xmin, xmax, maxfev=mymaxfev,
                                ftol=ftol, finalsimplex=9, step=mystep)
            x = numpy.asarray(result[1], numpy.float_)
            nfval = result[2]
            nfev = result[4].get('nfev')
        else:
            ncores_nm = ncoresNelderMead()
            nfev, nfval, x = \
                ncores_nm(stat_cb0, x, xmin, xmax, ftol, mymaxfev, numcores)

        if verbose or debug:
            print('f_nm%s=%.14e in %d nfev' % (x, nfval, nfev))
        ############################# NelderMead #############################

        ############################## nmDifEvo #############################
        xmin, xmax = _narrow_limits(4 * factor, [x, xmin, xmax], debug=False)
        mymaxfev = min(maxfev_per_iter, maxfev - nfev)
        if 1 == numcores:
            result = difevo_nm(myfcn, x, xmin, xmax, ftol, mymaxfev, verbose,
                               seed, pop, xprob, weight)
            nfev += result[4].get('nfev')
            x = numpy.asarray(result[1], numpy.float_)
            nfval = result[2]
        else:
            ncores_de = ncoresDifEvo()
            mystep = None
            tmp_nfev, tmp_fmin, tmp_par = \
                ncores_de(stat_cb0, x, xmin, xmax, ftol, mymaxfev, mystep,
                          numcores, pop, seed, weight, xprob, verbose)
            nfev += tmp_nfev
            if tmp_fmin < nfval:
                nfval = tmp_fmin
                x = tmp_par

        if verbose or debug:
            print('f_de_nm%s=%.14e in %d nfev' % (x, nfval, nfev))
        ############################## nmDifEvo #############################

        ofval = FUNC_MAX
        while nfev < maxfev:

            xmin, xmax = _narrow_limits(factor, [x, xmin, xmax], debug=False)

            ############################ nmDifEvo #############################
            y = random_start(xmin, xmax)
            mymaxfev = min(maxfev_per_iter, maxfev - nfev)
            if numcores == 1:
                result = difevo_nm(myfcn, y, xmin, xmax, ftol, mymaxfev,
                                   verbose, seed, pop, xprob, weight)
                nfev += result[4].get('nfev')
                if result[2] < nfval:
                    nfval = result[2]
                    x = numpy.asarray(result[1], numpy.float_)
                if verbose or debug:
                    print('f_de_nm%s=%.14e in %d nfev' %
                          (x, result[2], result[4].get('nfev')))
            ############################ nmDifEvo #############################

            if debug:
                print('ofval=%.14e\tnfval=%.14e\n' % (ofval, nfval))

            if sao_fcmp(ofval, nfval, ftol) <= 0:
                return x, nfval, nfev
            ofval = nfval
            factor *= 2

        return x, nfval, nfev

    x, fval, nfev = myopt(fcn, [x, xmin, xmax], numpy.sqrt(ftol), maxfev,
                          seed, population_size, xprob, weighting_factor,
                          factor=2.0, debug=False)

    if nfev < maxfev:
        if all(x == 0.0):
            mystep = list(map(lambda fubar: 1.2 + fubar, x))
        else:
            mystep = list(map(lambda fubar: 1.2 * fubar, x))
        if 1 == numcores:
            result = neldermead(fcn, x, xmin, xmax,
                                maxfev=min(512*len(x), maxfev - nfev),
                                ftol=ftol, finalsimplex=9, step=mystep)

            x = numpy.asarray(result[1], numpy.float_)
            fval = result[2]
            nfev += result[4].get('nfev')
        else:
            ncores_nm = ncoresNelderMead()
            tmp_nfev, tmp_fmin, tmp_par = \
                ncores_nm(stat_cb0, x, xmin, xmax, ftol, maxfev - nfev,
                          numcores)
            nfev += tmp_nfev
            # There is a bug here somewhere using broyden_tridiagonal
            if tmp_fmin < fval:
                fval = tmp_fmin
                x = tmp_par
    ierr = 0
    if nfev >= maxfev:
        ierr = 3
    status, msg = _get_saofit_msg(maxfev, ierr)

    rv = (status, x, fval, msg, {'info': status, 'nfev': nfev})
    return rv


#
# Nelder Mead
#
def neldermead(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None,
               initsimplex=0, finalsimplex=9, step=None, iquad=1,
               verbose=0, reflect=True):
    """Nelder-Mead Simplex optimization method.

    The Nelder-Mead Simplex algorithm, devised by J.A. Nelder and
    R. Mead [1]_, is a direct search method of optimization for
    finding a local minimum of an objective function of several
    variables. The implementation of the Nelder-Mead Simplex algorithm is
    a variation of the algorithm outlined in [2]_ and [3]_. As noted,
    terminating the simplex is not a simple task:

    "For any non-derivative method, the issue of termination is
    problematical as well as highly sensitive to problem scaling.
    Since gradient information is unavailable, it is provably
    impossible to verify closeness to optimality simply by sampling f
    at a finite number of points.  Most implementations of direct
    search methods terminate based on two criteria intended to reflect
    the progress of the algorithm: either the function values at the
    vertices are close, or the simplex has become very small."

    "Either form of termination-close function values or a small
    simplex-can be misleading for badly scaled functions."

    Parameters
    ----------
    fcn : function reference
       Returns the current statistic and per-bin statistic value when
       given the model parameters.
    x0, xmin, xmax : sequence of number
       The starting point, minimum, and maximum values for each
       parameter.
    ftol : number
       The function tolerance to terminate the search for the minimum;
       the default is sqrt(DBL_EPSILON) ~ 1.19209289551e-07, where
       DBL_EPSILON is the smallest number x such that ``1.0 != 1.0 +
       x``.
    maxfev : int or `None`
       The maximum number of function evaluations; the default value
       of `None` means to use ``1024 * n``, where `n` is the number of
       free parameters.
    initsimplex : int
       Dictates how the non-degenerate initial simplex is to be
       constructed.  Default is `0`; see the "cases for initsimplex"
       section below for details.
    finalsimplex : int
       At each iteration, a combination of one of the following
       stopping criteria is tested to see if the simplex has converged
       or not.  Full details are in the "cases for finalsimplex"
       section below.
    step : array of number or `None`
       A list of length `n` (number of free parameters) to initialize
       the simplex; see the `initsimplex` for details. The default of
       `None` means to use a step of 0.4 for each free parameter.
    iquad : int
       A boolean flag which indicates whether a fit to a quadratic
       surface is done.  If iquad is set to `1` (the default) then a
       fit to a quadratic surface is done; if iquad is set to `0` then
       the quadratic surface fit is not done.  If the fit to the
       quadratic surface is not positive semi-definitive, then the
       search terminated prematurely.  The code to fit the quadratic
       surface was written by D. E. Shaw, CSIRO, Division of
       Mathematics & Statistics, with amendments by
       R. W. M. Wedderburn, Rothamsted Experimental Station, and Alan
       Miller, CSIRO, Division of Mathematics & Statistics.  See also
       [1]_.
    verbose : int
       The amount of information to print during the fit. The default
       is `0`, which means no output.
    reflect : bool
       When a parameter exceeds a limit should the parameter be
       reflected, so moved back within bounds (`True`, the default) or
       should the model evaluation return DBL_MAX, causing the current
       set of parameters to be excluded from the simplex.

    Notes
    -----

    The `initsimplex` option determines how the non-degenerate initial
    simplex is to be constructed:

    - when `initsimplex` is `0`:

      Then x_(user_supplied) is one of the vertices of the simplex.
      The other `n` vertices are::

        for ( int i = 0; i &lt; n; ++i ) {
          for ( int j = 0; j &lt; n; ++j )
            x[ i + 1 ][ j ] = x_[ j ];
            x[ i + 1 ][ i ] = x_[ i ] + step[ i ];
        }

      where step[i] is the ith element of the option step.

    - if `initsimplex` is `1`:

      Then x_(user_supplied) is one of the vertices of the simplex.
      The other `n` vertices are::

                    { x_[j] + pn,   if i - 1 != j
                    {
        x[i][j]  =  {
                    {
                    { x_[j] + qn,   otherwise

      for 1 <= i <= n, 0 <= j < n and::

        pn = ( sqrt( n + 1 ) - 1 + n ) / ( n * sqrt(2) )
        qn = ( sqrt( n + 1 ) - 1 ) / ( n * sqrt(2) )

    The `finalsimplex` option determines whether the simplex has
    converged:

    - case a (if the max length of the simplex is small enough)::

        max( | x_i - x_0 | ) <= ftol max( 1, | x_0 | )
        1 <= i <= n

    - case b (if the standard deviation the simplex is < `ftol`)::

         n           -   2
        ===   ( f  - f )
        \        i                    2
        /     -----------     <=  ftol
        ====   sqrt( n )
        i = 0

    - case c (if the function values are close enough)::

        f_0  < f_(n-1)     within ftol

    The combination of the above stopping criteria are:

    - case 0: same as case a

    - case 1: case a, case b and case c have to be met

    - case 2: case a and either case b or case c have to be met.

    The `finalsimplex` value controls which of these criteria need to
    hold:

    - if ``finalsimplex=0`` then convergence is assumed if case 1 is met.

    - if ``finalsimplex=1`` then convergence is assumed if case 2 is met.

    - if ``finalsimplex=2`` then convergence is assumed if case 0 is met
      at two consecutive iterations.

    - if ``finalsimplex=3`` then convergence is assumed if case 0 then
      case 1 are met on two consecutive iterations.

    - if ``finalsimplex=4`` then convergence is assumed if case 0 then
      case 1 then case 0 are met on three consecutive iterations.

    - if ``finalsimplex=5`` then convergence is assumed if case 0 then
      case 1 then case 0 are met on three consecutive iterations.

    - if ``finalsimplex=6`` then convergence is assumed if case 1 then
      case 1 then case 0 are met on three consecutive iterations.

    - if ``finalsimplex=7`` then convergence is assumed if case 2 then
      case 1 then case 0 are met on three consecutive iterations.

    - if ``finalsimplex=8`` then convergence is assumed if case 0 then
      case 2 then case 0 are met on three consecutive iterations.

    - if ``finalsimplex=9`` then convergence is assumed if case 0 then
      case 1 then case 1 are met on three consecutive iterations.

    - if ``finalsimplex=10`` then convergence is assumed if case 0 then
      case 2 then case 1 are met on three consecutive iterations.

    - if ``finalsimplex=11`` then convergence is assumed if case 1 is
      met on three consecutive iterations.

    - if ``finalsimplex=12`` then convergence is assumed if case 1 then
      case 2 then case 1 are met on three consecutive iterations.

    - if ``finalsimplex=13`` then convergence is assumed if case 2 then
      case 1 then case 1 are met on three consecutive iterations.

    - otherwise convergence is assumed if case 2 is met on three
      consecutive iterations.

    References
    ----------

    .. [1] "A simplex method for function minimization", J.A. Nelder
           and R. Mead (Computer Journal, 1965, vol 7, pp 308-313)
           https://doi.org/10.1093%2Fcomjnl%2F7.4.308

    .. [2] "Convergence Properties of the Nelder-Mead Simplex
           Algorithm in Low Dimensions", Jeffrey C. Lagarias, James
           A. Reeds, Margaret H. Wright, Paul E. Wright , SIAM Journal
           on Optimization, Vol. 9, No. 1 (1998), pages 112-147.
           http://citeseer.ist.psu.edu/3996.html

    .. [3] "Direct Search Methods: Once Scorned, Now Respectable"
           Wright, M. H. (1996) in Numerical Analysis 1995
           (Proceedings of the 1995 Dundee Biennial Conference in
           Numerical Analysis, D.F. Griffiths and G.A. Watson, eds.),
           191-208, Addison Wesley Longman, Harlow, United Kingdom.
           http://citeseer.ist.psu.edu/155516.html

    """

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    order = 'F' if numpy.isfortran(x) else 'C'
    if step is None or (numpy.iterable(step) and len(step) != len(x)):
        step = 1.2 * numpy.ones(x.shape, numpy.float_, order)
    elif numpy.isscalar(step):
        step = step * numpy.ones(x.shape, numpy.float_, order)

    def stat_cb0(pars):
        return fcn(pars)[0]

    #
    # A safeguard just in case the initial simplex is outside the bounds
    #
    orig_fcn = stat_cb0

    def stat_cb0(x_new):
        if _my_is_nan(x_new) or _outside_limits(x_new, xmin, xmax):
            return FUNC_MAX
        return orig_fcn(x_new)

    # for internal use only
    debug = False

    if numpy.isscalar(finalsimplex) and numpy.iterable(finalsimplex) == 0:
        finalsimplex = int(finalsimplex)
        if 0 == finalsimplex:
            finalsimplex = [1]
        elif 1 == finalsimplex:
            finalsimplex = [2]
        elif 2 == finalsimplex:
            finalsimplex = [0, 0]
        elif 3 == finalsimplex:
            finalsimplex = [0, 1]
        elif 4 == finalsimplex:
            finalsimplex = [0, 1, 0]
        elif 5 == finalsimplex:
            finalsimplex = [0, 2, 0]
        elif 6 == finalsimplex:
            finalsimplex = [1, 1, 0]
        elif 7 == finalsimplex:
            finalsimplex = [2, 1, 0]
        elif 8 == finalsimplex:
            finalsimplex = [1, 2, 0]
        elif 9 == finalsimplex:
            finalsimplex = [0, 1, 1]
        elif 10 == finalsimplex:
            finalsimplex = [0, 2, 1]
        elif 11 == finalsimplex:
            finalsimplex = [1, 1, 1]
        elif 12 == finalsimplex:
            finalsimplex = [1, 2, 1]
        elif 13 == finalsimplex:
            finalsimplex = [2, 1, 1]
        else:
            finalsimplex = [2, 2, 2]
    elif (not numpy.isscalar(finalsimplex) and
          numpy.iterable(finalsimplex) == 1):
        pass
    else:
        finalsimplex = [2, 2, 2]

    finalsimplex = numpy.asarray(finalsimplex, numpy.int_)

    if maxfev is None:
        maxfev = 1024 * len(x)

    if debug:
        print('opfcts.py neldermead() finalsimplex=%s\tisscalar=%s\titerable=%d' % (finalsimplex, numpy.isscalar(finalsimplex), numpy.iterable(finalsimplex)))

    def simplex(verbose, maxfev, init, final, tol, step, xmin, xmax, x,
                myfcn, debug, ofval=FUNC_MAX):

        tmpfinal = final[:]
        if len(final) >= 3:
            # get rid of the last entry in the list
            tmpfinal = final[0:-1]

        xx, ff, nf, er = _saoopt.neldermead(verbose, maxfev, init, tmpfinal,
                                            tol, step, xmin, xmax, x, myfcn)

        if debug:
            print('finalsimplex=%s, nfev=%d:\tf%s=%.20e' % (tmpfinal, nf, xx, ff))

        if len(final) >= 3 and ff < 0.995 * ofval and nf < maxfev:
            myfinal = [final[-1]]
            x, fval, nfev, err = simplex(verbose, maxfev-nf, init, myfinal, tol,
                                         step, xmin, xmax, x, myfcn, debug,
                                         ofval=ff)
            return x, fval, nfev + nf, err
        else:
            return xx, ff, nf, er

    x, fval, nfev, ier = simplex(verbose, maxfev, initsimplex, finalsimplex,
                                 ftol, step, xmin, xmax, x, stat_cb0, debug)
    if debug:
        print('f%s=%e in %d nfev' % (x, fval, nfev))

    info = 1
    covarerr = None
    if len(finalsimplex) >= 3 and 0 != iquad:
        nelmea = minim(fcn, x, xmin, xmax, ftol=10.0*ftol,
                       maxfev=maxfev - nfev - 12, iquad=1, reflect=reflect)
        nelmea_x = numpy.asarray(nelmea[1], numpy.float_)
        nelmea_nfev = nelmea[4].get('nfev')
        info = nelmea[4].get('info')
        covarerr = nelmea[4].get('covarerr')
        nfev += nelmea_nfev
        minim_fval = nelmea[2]
        if minim_fval < fval:
            x = nelmea_x
            fval = minim_fval
        if debug:
            print('minim: f%s=%e %d nfev, info=%d' % (x, fval, nelmea_nfev, info))

    if nfev >= maxfev:
        ier = 3
    key = {
        0: (True, 'Optimization terminated successfully'),
        1: (False, 'improper input parameters'),
        2: (False, 'improper values for x, xmin or xmax'),
        3: (False,
            'number of function evaluations has exceeded %d' % maxfev)
        }
    status, msg = key.get(ier,
                          (False, 'unknown status flag (%d)' % ier))

    rv = (status, x, fval)
    print_covar_err = False
    if print_covar_err and covarerr is not None:
        rv += (msg, {'covarerr': covarerr, 'info': status, 'nfev': nfev})
    else:
        rv += (msg, {'info': status, 'nfev': nfev})
    return rv


def lmdif(fcn, x0, xmin, xmax, ftol=EPSILON, xtol=EPSILON, gtol=EPSILON,
          maxfev=None, epsfcn=EPSILON, factor=100.0, numcores=1, verbose=0):
    """Levenberg-Marquardt optimization method.

    The Levenberg-Marquardt method is an interface to the MINPACK
    subroutine lmdif to find the local minimum of nonlinear least
    squares functions of several variables by a modification of the
    Levenberg-Marquardt algorithm [1]_.

    Parameters
    ----------
    fcn : function reference
       Returns the current statistic and per-bin statistic value when
       given the model parameters.
    x0, xmin, xmax : sequence of number
       The starting point, minimum, and maximum values for each
       parameter.
    ftol : number
       The function tolerance to terminate the search for the minimum;
       the default is FLT_EPSILON ~ 1.19209289551e-07, where
       FLT_EPSILON is the smallest number x such that ``1.0 != 1.0 +
       x``. The conditions are satisfied when both the actual and
       predicted relative reductions in the sum of squares are, at
       most, ftol.
    xtol : number
       The relative error desired in the approximate solution; default
       is FLT_EPSILON ~ 1.19209289551e-07, where FLT_EPSILON
       is the smallest number x such that ``1.0 != 1.0 + x``. The
       conditions are satisfied when the relative error between two
       consecutive iterates is, at most, `xtol`.
    gtol : number
       The orthogonality desired between the function vector and the
       columns of the jacobian; default is FLT_EPSILON ~
       1.19209289551e-07, where FLT_EPSILON is the smallest number x
       such that ``1.0 != 1.0 + x``. The conditions are satisfied when
       the cosine of the angle between fvec and any column of the
       jacobian is, at most, `gtol` in absolute value.
    maxfev : int or `None`
       The maximum number of function evaluations; the default value
       of `None` means to use ``1024 * n``, where `n` is the number of
       free parameters.
    epsfcn : number
       This is used in determining a suitable step length for the
       forward-difference approximation; default is FLT_EPSILON
       ~ 1.19209289551e-07, where FLT_EPSILON is the smallest number
       x such that ``1.0 != 1.0 + x``. This approximation assumes that
       the relative errors in the functions are of the order of
       `epsfcn`. If `epsfcn` is less than the machine precision, it is
       assumed that the relative errors in the functions are of the
       order of the machine precision.
    factor : int
       Used in determining the initial step bound; default is 100. The
       initial step bound is set to the product of `factor` and the
       euclidean norm of diag*x if nonzero, or else to factor itself.
       In most cases, `factor` should be from the interval (.1,100.).
    numcores : int
       The number of CPU cores to use. The default is `1`.
    verbose: int
       The amount of information to print during the fit. The default
       is `0`, which means no output.

    References
    ----------

    .. [1] J.J. More, "The Levenberg Marquardt algorithm:
           implementation and theory," in Lecture Notes in Mathematics
           630: Numerical Analysis, G.A. Watson (Ed.),
           Springer-Verlag: Berlin, 1978, pp.105-116.

    """

    class fdJac:

        def __init__(self, func, fvec, pars):
            self.func = func
            self.fvec = fvec
            epsmch = numpy.finfo(float).eps
            self.eps = numpy.sqrt(max(epsmch, epsfcn))
            self.h = self.calc_h(pars)
            self.pars = numpy.copy(pars)
            return

        def __call__(self, param):
            wa = self.func(param[1:])
            return (wa - self.fvec) / self.h[int(param[0])]

        def calc_h(self, pars):
            nn = len(pars)
            h = numpy.empty((nn,))
            for ii in range(nn):
                h[ii] = self.eps * pars[ii]
                if h[ii] == 0.0:
                    h[ii] = self.eps
                if pars[ii] + h[ii] > xmax[ii]:
                    h[ii] = - h[ii]
            return h

        def calc_params(self):
            h = self.calc_h(self.pars)
            params = []
            for ii in range(len(h)):
                tmp_pars = numpy.copy(self.pars)
                tmp_pars[ii] += h[ii]
                tmp_pars = numpy.append(ii, tmp_pars)
                params.append(tmp_pars)
            return tuple(params)

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if maxfev is None:
        maxfev = 256 * len(x)

    def stat_cb0(pars):
        return fcn(pars)[0]

    def stat_cb1(pars):
        return fcn(pars)[1]

    def fcn_parallel(pars, fvec):
        fd_jac = fdJac(stat_cb1, fvec, pars)
        params = fd_jac.calc_params()
        fjac = parallel_map(fd_jac, params, numcores)
        return numpy.concatenate(fjac)

    num_parallel_map, fcn_parallel_counter = func_counter(fcn_parallel)

    # TO DO: reduce 1 model eval by passing the resulting 'fvec' to cpp_lmdif
    m = numpy.asanyarray(stat_cb1(x)).size

    error = []

    n = len(x)
    fjac = numpy.empty((m*n,))

    x, fval, nfev, info, fjac = \
        _saoopt.cpp_lmdif(stat_cb1, fcn_parallel_counter, numcores, m, x, ftol,
                          xtol, gtol, maxfev, epsfcn, factor, verbose, xmin,
                          xmax, fjac)

    if info > 0:
        fjac = numpy.reshape(numpy.ravel(fjac, order='F'), (m, n), order='F')

        if m != n:
            covar = fjac[:n, :n]
        else:
            covar = fjac

        if _par_at_boundary(xmin, x, xmax, xtol):
            nm_result = neldermead(fcn, x, xmin, xmax, ftol=numpy.sqrt(ftol),
                                   maxfev=maxfev-nfev, finalsimplex=2, iquad=0,
                                   verbose=0)
            nfev += nm_result[4]['nfev']
            x = nm_result[1]
            fval = nm_result[2]

    if error:
        raise error.pop()

    if 0 == info:
        info = 1
    elif info >= 1 or info <= 4:
        info = 0
    else:
        info = 3
    status, msg = _get_saofit_msg(maxfev, info)

    if info == 0:
        rv = (status, x, fval, msg, {'info': info, 'nfev': nfev,
                                     'covar': covar,
                                     'num_parallel_map': num_parallel_map[0]})
    else:
        rv = (status, x, fval, msg, {'info': info, 'nfev': nfev,
                                     'num_parallel_map': num_parallel_map[0]})
    return rv
