#
#  Copyright (C) 2007, 2016, 2018 - 2024
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

>>> import numpy as np
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

import numpy as np

from sherpa.optmethods.ncoresde import ncoresDifEvo
from sherpa.optmethods.ncoresnm import ncoresNelderMead

from sherpa.utils import FuncCounter
from sherpa.utils.parallel import parallel_map
from sherpa.utils._utils import sao_fcmp  # type: ignore
from sherpa.utils import random
from sherpa.utils.types import ArrayType

from . import _saoopt  # type: ignore

__all__ = ('difevo', 'difevo_lm', 'difevo_nm', 'grid_search', 'lmdif',
           'minim', 'montecarlo', 'neldermead')


#
# Use FLT_EPSILON as default tolerance
#
EPSILON = np.float64(np.finfo(np.float32).eps)

#
# Maximum callback function value, used to indicate that the optimizer
# has exceeded parameter boundaries.  All the optimizers expect double
# precision arguments, so we use np.float64 instead of SherpaFloat.
#
FUNC_MAX = np.finfo(np.float64).max


def _check_args(x0: ArrayType,
                xmin: ArrayType,
                xmax: ArrayType
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert to ndarray, check shape, and ensure x is within (xmin,xmax).

    Thr x0 array is copied (so that changes to it do not affect the
    input x0 argument).

    """

    x = np.array(x0, np.float64)  # Make a copy
    xmin = np.asarray(xmin, np.float64)
    xmax = np.asarray(xmax, np.float64)

    if (x.shape != xmin.shape) or (x.shape != xmax.shape):
        raise TypeError('input array sizes do not match')

    xclip = np.clip(x, xmin, xmax)
    return xclip, xmin, xmax


def _get_saofit_msg(maxfev: int,
                    ierr: int
                    ) -> tuple[bool, str]:
    key = {
        0: (True, 'successful termination'),
        1: (False, 'improper input parameters'),
        2: (False, 'initial parameter value is out of bounds'),
        3: (False,
            f'number of function evaluations has exceeded maxfev={maxfev}')
        }
    return key.get(ierr, (False, f'unknown status flag ({ierr})'))


def _raise_min_limit(factor: float,
                     xmin: np.ndarray,
                     x: np.ndarray
                     ) -> np.ndarray:
    """Calculate the new minimum limits."""

    myxmin = x - factor * np.abs(x)
    return np.clip(myxmin, xmin, None)


def _lower_max_limit(factor: float,
                     x: np.ndarray,
                     xmax: np.ndarray
                     ) -> np.ndarray:
    """Calculate the new maximum limits."""

    myxmax = x + factor * np.abs(x)
    return np.clip(myxmax, None, xmax)


def _narrow_limits(factor: float,
                   x: np.ndarray,
                   xmin: np.ndarray,
                   xmax: np.ndarray
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Do we need to change the limits?"""

    myxmin = _raise_min_limit(factor, xmin, x)
    myxmax = _lower_max_limit(factor, x, xmax)
    return myxmin, myxmax


def _par_at_boundary(low, val, high, tol):
    for par_min, par_val, par_max in zip(low, val, high):
        if sao_fcmp(par_val, par_min, tol) == 0:
            return True
        if sao_fcmp(par_val, par_max, tol) == 0:
            return True
    return False


def _outside_limits(x, xmin, xmax):
    return (np.any(x < xmin) or np.any(x > xmax))


def difevo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
           seed=2005815, population_size=None, xprob=0.9,
           weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = float(np.clip(xprob, 0.1, 1.0))

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = float(np.clip(weighting_factor, 0.1, 1.0))

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
        print(f'difevo: f{x}={fval:e} in {nfev} nfev')

    status, msg = _get_saofit_msg(maxfev, ierr)
    return (status, x, fval, msg, {'info': ierr, 'nfev': nfev})


def difevo_lm(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
              seed=2005815, population_size=None, xprob=0.9,
              weighting_factor=0.8):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = float(np.clip(xprob, 0.1, 1.0))

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = float(np.clip(weighting_factor, 0.1, 1.0))

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * x.size

    # TODO: can we not just call x.size rather than
    #       np.asanyarray(fcn(x)).size for the last argument?
    #
    de = _saoopt.lm_difevo(verbose, maxfev, seed, population_size, ftol,
                           xprob, weighting_factor, xmin, xmax,
                           x, fcn, np.asanyarray(fcn(x)).size)
    fval = de[1]
    nfev = de[2]
    ierr = de[3]

    status, msg = _get_saofit_msg(maxfev, ierr)
    return (status, x, fval, msg, {'info': ierr, 'nfev': nfev})


def difevo_nm(fcn, x0, xmin, xmax, ftol, maxfev, verbose, seed,
              population_size, xprob, weighting_factor):

    def stat_cb0(pars):
        return fcn(pars)[0]

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    # make sure that the cross over prob is within [0.1,1.0]
    xprob = float(np.clip(xprob, 0.1, 1.0))

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = float(np.clip(weighting_factor, 0.1, 1.0))

    if population_size is None:
        population_size = 16 * x.size

    if maxfev is None:
        maxfev = 1024 * population_size

    de = _saoopt.nm_difevo(verbose, maxfev, seed, population_size,
                           ftol, xprob, weighting_factor, xmin, xmax,
                           x, stat_cb0)
    fval = de[1]
    nfev = de[2]
    ierr = de[3]

    if verbose:
        print('difevo_nm: f{x}={fval:e} in {nfev} nfev')

    status, msg = _get_saofit_msg(maxfev, ierr)
    return (status, x, fval, msg, {'info': ierr, 'nfev': nfev})


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
            print(f'f{pars}={aaa:g}')
        return aaa

    def make_sequence(ranges, N):
        list_ranges = list(ranges)
        for ii in range(npar):
            list_ranges[ii] = tuple(list_ranges[ii]) + (complex(N),)
            list_ranges[ii] = slice(*list_ranges[ii])

        grid = np.mgrid[list_ranges]
        mynfev = pow(N, npar)
        grid = list(map(np.ravel, grid))
        sequence = []
        for index in range(mynfev):
            tmp = []
            for xx in range(npar):
                tmp.append(grid[xx][index])
            sequence.append(tmp)
        return sequence

    def eval_stat_func(xxx):
        return np.append(func(xxx), xxx)

    if sequence is None:
        ranges = []
        for index in range(npar):
            ranges.append([xmin[index], xmax[index]])
        sequence = make_sequence(ranges, num)
    else:
        if not np.iterable(sequence):
            raise TypeError("sequence option must be iterable")

        for seq in sequence:
            if npar != len(seq):
                raise TypeError(f"{seq} must be of length {npar}")

    answer = eval_stat_func(x)
    sequence_results = parallel_map(eval_stat_func, sequence, numcores)
    for xresult in sequence_results[1:]:
        if xresult[0] < answer[0]:
            answer = xresult

    x = answer[1:]
    nfev = len(sequence_results) + 1

    # TODO: should we just use case-insensitive comparison?
    if method in ['NelderMead', 'neldermead', 'Neldermead', 'nelderMead']:
        # re.search( '^[Nn]elder[Mm]ead', method ):
        nm_result = neldermead(fcn, x, xmin, xmax, ftol=ftol, maxfev=maxfev,
                               verbose=verbose)
        (status, x, fval, msg, imap) = nm_result
        imap['nfev'] += nfev
        return (status, x, fval, msg, imap)

    if method in ['LevMar', 'levmar', 'Levmar', 'levMar']:
        # re.search( '^[Ll]ev[Mm]ar', method ):
        levmar_result = lmdif(fcn, x, xmin, xmax, ftol=ftol, xtol=ftol,
                              gtol=ftol, maxfev=maxfev, verbose=verbose)
        (status, x, fval, msg, imap) = levmar_result
        imap['nfev'] += nfev
        return (status, x, fval, msg, imap)

    fval = answer[0]
    ierr = 0
    status, msg = _get_saofit_msg(ierr, ierr)
    return (status, x, fval, msg, {'info': ierr, 'nfev': nfev})


#
# C-version of minim
#
def minim(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, step=None,
          nloop=1, iquad=1, simp=None, verbose=-1, reflect=True):

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if step is None:
        step = np.full(x.shape, 0.4, dtype=np.float64)

    if simp is None:
        simp = 1.0e-2 * ftol

    if maxfev is None:
        maxfev = 512 * len(x)

    def stat_cb0(x_new):
        if np.isnan(x_new).any() or _outside_limits(x_new, xmin, xmax):
            return FUNC_MAX
        return fcn(x_new)[0]

    init = 0
    x, fval, neval, ifault = _saoopt.minim(reflect, verbose, maxfev, init, \
                                           iquad, simp, ftol, step, \
                                           xmin, xmax, x, stat_cb0)

    key = {
        0: (True, 'successful termination'),
        1: (False,
            f'number of function evaluations has exceeded maxfev={maxfev}'),
        2: (False, 'information matrix is not +ve semi-definite'),
        3: (False, 'number of parameters is less than 1'),
        4: (False, f'nloop={nloop} is less than 1')
        }
    status, msg = key.get(ifault, (False, f'unknown status flag ({ifault})'))

    return (status, x, fval, msg, {'info': ifault, 'nfev': neval})


#
# Monte Carlo
#
def montecarlo(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None, verbose=0,
               seed=74815, population_size=None, xprob=0.9,
               weighting_factor=0.8, numcores=1, rng=None):
    """Monte Carlo optimization method.

    This is an implementation of the differential-evolution algorithm
    from Storn and Price (1997) [1]_. A population of fixed size -
    which contains n-dimensional vectors, where n is the number of
    free parameters - is randomly initialized.  At each iteration, a
    new n-dimensional vector is generated by combining vectors from
    the pool of population, the resulting trial vector is selected if
    it lowers the objective function.

    .. versionchanged:: 4.16.0
       The rng parameter was added.

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
    seed : int or None
       The seed for the random number generator. If not set then the
       rng parameter is used to create a seed value.
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
    rng : np.random.Generator, np.random.RandomState, or None, optional
       Determines how the random numbers are created. If set to None
       then the routines from `numpy.random` are used, and so can be
       controlled by calling `numpy.random.seed`.

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
    xprob = float(np.clip(xprob, 0.1, 1.0))

    # make sure that weighting_factor is within [0.1,1.0]
    weighting_factor = float(np.clip(weighting_factor, 0.1, 1.0))

    # Do we need to create a seed?
    #
    # We use a seed up to
    #   pow(2, 31) == 2147483648L
    #
    if seed is None:
        seed = random.integers(rng, 2147483648)

    if population_size is None:
        population_size = 12 * x.size

    if maxfev is None:
        maxfev = 8192 * population_size

    def myopt(myfcn, xxx, ftol, maxfev, seed, pop, xprob,
              weight, factor=4.0):

        x = xxx[0]
        xmin = xxx[1]
        xmax = xxx[2]
        maxfev_per_iter = 512 * x.size

        def random_start(xmin, xmax):
            xx = []
            for ii in range(len(xmin)):
                xx.append(random.uniform(rng, xmin[ii], xmax[ii]))
            return np.asarray(xx)

        ############################# NelderMead #############################
        mymaxfev = min(maxfev_per_iter, maxfev)
        if all(x == 0.0):
            mystep = 1.2 + x
        else:
            mystep = 1.2 * x

        if 1 == numcores:
            result = neldermead(myfcn, x, xmin, xmax, maxfev=mymaxfev,
                                ftol=ftol, finalsimplex=9, step=mystep)
            x = np.asarray(result[1], np.float64)
            nfval = result[2]
            nfev = result[4].get('nfev')
        else:
            ncores_nm = ncoresNelderMead()
            nfev, nfval, x = \
                ncores_nm(stat_cb0, x, xmin, xmax, ftol, mymaxfev, numcores)

        if verbose:
            print(f'f_nm{x}={nfval:.14e} in {nfev} nfev')

        ############################# NelderMead #############################

        ############################## nmDifEvo #############################
        xmin, xmax = _narrow_limits(4 * factor, x, xmin, xmax)
        mymaxfev = min(maxfev_per_iter, maxfev - nfev)
        if 1 == numcores:
            result = difevo_nm(myfcn, x, xmin, xmax, ftol, mymaxfev, verbose,
                               seed, pop, xprob, weight)
            nfev += result[4].get('nfev')
            x = np.asarray(result[1], np.float64)
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

        if verbose:
            print(f'f_de_nm{x}={nfval:.14e} in {nfev} nfev')

        ############################## nmDifEvo #############################

        ofval = FUNC_MAX
        while nfev < maxfev:

            xmin, xmax = _narrow_limits(factor, x, xmin, xmax)

            ############################ nmDifEvo #############################
            y = random_start(xmin, xmax)
            mymaxfev = min(maxfev_per_iter, maxfev - nfev)

            if numcores == 1:
                # TODO: should this update the seed somehow?
                result = difevo_nm(myfcn, y, xmin, xmax, ftol, mymaxfev,
                                   verbose, seed, pop, xprob, weight)
                nfev += result[4].get('nfev')
                if result[2] < nfval:
                    nfval = result[2]
                    x = np.asarray(result[1], np.float64)
                if verbose:
                    print(f'f_de_nm{x}={result[2]:.14e} in {result[4].get("nfev")} nfev')

            ############################ nmDifEvo #############################

            if sao_fcmp(ofval, nfval, ftol) <= 0:
                return x, nfval, nfev

            ofval = nfval
            factor *= 2

        return x, nfval, nfev

    x, fval, nfev = myopt(fcn, [x, xmin, xmax], np.sqrt(ftol), maxfev,
                          seed, population_size, xprob, weighting_factor,
                          factor=2.0)

    if nfev < maxfev:
        if all(x == 0.0):
            mystep = 1.2 + x
        else:
            mystep = 1.2 * x

        if 1 == numcores:
            result = neldermead(fcn, x, xmin, xmax,
                                maxfev=min(512*len(x), maxfev - nfev),
                                ftol=ftol, finalsimplex=9, step=mystep)

            x = np.asarray(result[1], np.float64)
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
    return (status, x, fval, msg, {'info': status, 'nfev': nfev})


#
# Nelder Mead
#
def neldermead(fcn, x0, xmin, xmax, ftol=EPSILON, maxfev=None,
               initsimplex=0, finalsimplex=9, step=None, iquad=1,
               verbose=0, reflect=True):
    r"""Nelder-Mead Simplex optimization method.

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

    if step is None or (np.iterable(step) and len(step) != len(x)):
        step = np.full(x.shape, 1.2, dtype=np.float64)
    elif np.isscalar(step):
        step = np.full(x.shape, step, dtype=np.float64)

    # A safeguard just in case the initial simplex is outside the bounds
    #
    def stat_cb0(x_new):
        if np.isnan(x_new).any() or _outside_limits(x_new, xmin, xmax):
            return FUNC_MAX
        return fcn(x_new)[0]

    if np.isscalar(finalsimplex) and not np.iterable(finalsimplex):
        farg = int(finalsimplex)
        if 0 == farg:
            finalsimplex_ary = [1]
        elif 1 == farg:
            finalsimplex_ary = [2]
        elif 2 == farg:
            finalsimplex_ary = [0, 0]
        elif 3 == farg:
            finalsimplex_ary = [0, 1]
        elif 4 == farg:
            finalsimplex_ary = [0, 1, 0]
        elif 5 == farg:
            finalsimplex_ary = [0, 2, 0]
        elif 6 == farg:
            finalsimplex_ary = [1, 1, 0]
        elif 7 == farg:
            finalsimplex_ary = [2, 1, 0]
        elif 8 == farg:
            finalsimplex_ary = [1, 2, 0]
        elif 9 == farg:
            finalsimplex_ary = [0, 1, 1]
        elif 10 == farg:
            finalsimplex_ary = [0, 2, 1]
        elif 11 == farg:
            finalsimplex_ary = [1, 1, 1]
        elif 12 == farg:
            finalsimplex_ary = [1, 2, 1]
        elif 13 == farg:
            finalsimplex_ary = [2, 1, 1]
        else:
            finalsimplex_ary = [2, 2, 2]
    elif (not np.isscalar(finalsimplex) and np.iterable(finalsimplex)):
        # support for finalsimplex being a sequence is not documented
        # and not tested
        finalsimplex_ary = finalsimplex
    else:
        finalsimplex_ary = [2, 2, 2]

    fsimplex = np.asarray(finalsimplex_ary, np.int_)

    if maxfev is None:
        maxfev = 1024 * len(x)

    def simplex(verbose, maxfev, init, final, tol, step, xmin, xmax, x,
                myfcn, ofval=FUNC_MAX):

        tmpfinal = final[:]
        if len(final) >= 3:
            # get rid of the last entry in the list
            tmpfinal = final[0:-1]

        xx, ff, nf, er = _saoopt.neldermead(verbose, maxfev, init, tmpfinal,
                                            tol, step, xmin, xmax, x, myfcn)

        if len(final) >= 3 and ff < 0.995 * ofval and nf < maxfev:
            myfinal = [final[-1]]
            x, fval, nfev, err = simplex(verbose, maxfev-nf, init, myfinal, tol,
                                         step, xmin, xmax, x, myfcn,
                                         ofval=ff)
            return x, fval, nfev + nf, err

        return xx, ff, nf, er

    x, fval, nfev, ier = simplex(verbose, maxfev, initsimplex, fsimplex,
                                 ftol, step, xmin, xmax, x, stat_cb0)

    covarerr = None
    if len(fsimplex) >= 3 and 0 != iquad:
        nelmea = minim(fcn, x, xmin, xmax, ftol=10.0*ftol,
                       maxfev=maxfev - nfev - 12, iquad=1, reflect=reflect)
        nelmea_x = np.asarray(nelmea[1], np.float64)
        nelmea_nfev = nelmea[4].get('nfev')
        covarerr = nelmea[4].get('covarerr')
        nfev += nelmea_nfev
        minim_fval = nelmea[2]

        # Have we found a better location?
        if minim_fval < fval:
            x = nelmea_x
            fval = minim_fval

    if nfev >= maxfev:
        ier = 3

    key = {
        0: (True, 'Optimization terminated successfully'),
        1: (False, 'improper input parameters'),
        2: (False, 'improper values for x, xmin or xmax'),
        3: (False,
            f'number of function evaluations has exceeded {maxfev}')
        }
    status, msg = key.get(ier,
                          (False, f'unknown status flag ({ier})'))

    imap = {'info': status, 'nfev': nfev}
    print_covar_err = False
    if print_covar_err and covarerr is not None:
        imap['covarerr'] = covarerr

    return (status, x, fval, msg, imap)


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
            epsmch = np.finfo(float).eps
            self.eps = np.sqrt(max(epsmch, epsfcn))
            self.h = self.calc_h(pars)
            self.pars = np.copy(pars)

        def __call__(self, param):
            wa = self.func(param[1:])
            return (wa - self.fvec) / self.h[int(param[0])]

        def calc_h(self, pars):
            nn = len(pars)
            h = np.empty((nn,))
            for ii in range(nn):
                h[ii] = self.eps * pars[ii]
                if h[ii] == 0.0:
                    h[ii] = self.eps
                if pars[ii] + h[ii] > xmax[ii]:
                    h[ii] = - h[ii]
            return h

        def calc_params(self):
            params = []
            for ii in range(len(self.h)):
                tmp_pars = np.copy(self.pars)
                tmp_pars[ii] += self.h[ii]
                tmp_pars = np.append(ii, tmp_pars)
                params.append(tmp_pars)
            return tuple(params)

    x, xmin, xmax = _check_args(x0, xmin, xmax)

    if maxfev is None:
        maxfev = 256 * len(x)

    def stat_cb1(pars):
        return fcn(pars)[1]

    def fcn_parallel(pars, fvec):
        fd_jac = fdJac(stat_cb1, fvec, pars)
        params = fd_jac.calc_params()
        fjac = parallel_map(fd_jac, params, numcores)
        return np.concatenate(fjac)

    fcn_parallel_counter = FuncCounter(fcn_parallel)

    # TO DO: reduce 1 model eval by passing the resulting 'fvec' to cpp_lmdif
    m = np.asanyarray(stat_cb1(x)).size

    n = len(x)
    fjac = np.empty((m*n,))

    x, fval, nfev, info, fjac = \
        _saoopt.cpp_lmdif(stat_cb1, fcn_parallel_counter, numcores, m, x, ftol,
                          xtol, gtol, maxfev, epsfcn, factor, verbose, xmin,
                          xmax, fjac)

    if info > 0:
        fjac = np.reshape(np.ravel(fjac, order='F'), (m, n), order='F')

        if m != n:
            covar = fjac[:n, :n]
        else:
            covar = fjac

        if _par_at_boundary(xmin, x, xmax, xtol):
            nm_result = neldermead(fcn, x, xmin, xmax, ftol=np.sqrt(ftol),
                                   maxfev=maxfev-nfev, finalsimplex=2, iquad=0,
                                   verbose=0)
            nfev += nm_result[4]['nfev']
            x = nm_result[1]
            fval = nm_result[2]

    if 0 == info:
        info = 1
    elif info >= 1 or info <= 4:
        info = 0
    else:
        info = 3
    status, msg = _get_saofit_msg(maxfev, info)

    imap = {'info': info, 'nfev': nfev,
            'num_parallel_map': fcn_parallel_counter.nfev}
    if info == 0:
        imap['covar'] = covar

    return (status, x, fval, msg, imap)
