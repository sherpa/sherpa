#
#  Copyright (C) 2007, 2015, 2016, 2019-2021, 2023-2025
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

from __future__ import annotations

from collections.abc import Callable, Sequence
import logging
from typing import Any, Protocol, SupportsFloat
import warnings

import numpy as np
from numpy.linalg import LinAlgError

from sherpa.stats import StatCallback
from sherpa.utils import NoNewAttributesAfterInit, \
    FuncCounter, OutOfBoundErr, Knuth_close, \
    print_fields, is_iterable, list_to_open_interval, quad_coef, \
    demuller, zeroin
from sherpa.utils.parallel import SupportsLock, SupportsProcess, \
    SupportsQueue, multi, ncpus, context, process_tasks
from sherpa.utils.types import ArrayType, FitFunc, StatFunc

from . import _est_funcs  # type: ignore


# TODO: this should not be set globally
_ = np.seterr(invalid='ignore')


__all__ = ('EstNewMin', 'Covariance', 'Confidence',
           'Projection', 'est_success', 'est_failure', 'est_hardmin',
           'est_hardmax', 'est_hardminmax', 'est_newmin', 'est_maxiter',
           'est_hitnan')


warning = logging.getLogger(__name__).warning


# The return type for the estimation routines.
#
# It looks like the limits can be numbers or None.
#
ArrayVal = Sequence[SupportsFloat | None] | np.ndarray
EstReturn = tuple[ArrayVal, ArrayVal,
                  Sequence[int | None] | np.ndarray,
                  int, np.ndarray | None]


est_success = 0
est_failure = 1
est_hardmin = 2
est_hardmax = 3
est_hardminmax = 4
est_newmin = 5
est_maxiter = 6
est_hitnan = 7


# This class is used to send around the "new parameter values" in
# sherpa.fit.Fit.est_errors.
#
class EstNewMin(Exception):
    "Reached a new minimum fit statistic"
    pass


class EstMethod(NoNewAttributesAfterInit):
    """Estimate errors on a set of parameters.

    .. versionchanged:: 4.17.1
       The estfunc argument is now unused and will be removed in a
       later release.

    """

    # defined pre-instantiation for pickling
    config = {'sigma': 1,
              'eps': 0.01,
              'maxiters': 200,
              'soft_limits': False}

    def __init__(self,
                 name: str,
                 estfunc: Any = None
                 ) -> None:
        self.name = name

        if estfunc is not None:
            warnings.warn("EstMethod: the estfunc argument is deprecated",
                          DeprecationWarning)

        # config should be defined pre-instantiation for pickling
        # however, for some unknown reason membership in self.__dict__
        # requires declaration in __init__()
        self.config = self.config.copy()

        super().__init__()

    def __getattr__(self, name):
        if name in self.__dict__.get('config', ()):
            return self.config[name]
        raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

    def __setattr__(self, name, val):
        if name in self.__dict__.get('config', ()):
            self.config[name] = val
        else:
            NoNewAttributesAfterInit.__setattr__(self, name, val)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} error-estimation method instance '{self.name}'>"

    def __str__(self) -> str:
        # Put name first always
        keylist = list(self.config.keys())
        keylist = ['name'] + keylist
        full_config = {'name': self.name}
        full_config.update(self.config)

        return print_fields(keylist, full_config)

    def __setstate__(self, state):
        self.__dict__.update(state)

        # obtain config values from object class
        # TODO: pylint points out there's no name set for
        #       the initialization call
        self.__dict__['config'] = getattr(self.__class__(), 'config', {})

        # update new config dict with user defined from old
        self.__dict__['config'].update(state.get('config', {}))

    def compute(self,
                statfunc: StatFunc,
                fitfunc: FitFunc,
                *,
                pars: np.ndarray,
                parmins: np.ndarray,
                parmaxes: np.ndarray,
                parhardmins: np.ndarray,
                parhardmaxes: np.ndarray,
                limit_parnums: np.ndarray,  # integers
                freeze_par: Callable,
                thaw_par: Callable,
                report_progress: Callable,
                get_par_name: Callable,
                statargs: Any = None,
                statkwargs: Any = None
                ) -> EstReturn:
        """Estimate the error range.

        .. versionchanged:: 4.17.1
           All arguments other than statfunc and fitfunc are now
           keyword-only.

        Notes
        -----
        The statargs and statkwargs arguments are currently unused.

        """

        # This does not use abc.abstractmethod as inheriting from
        # NoNewAttributesAfterInit.
        #
        raise NotImplementedError()


class Covariance(EstMethod):
    """The covariance method for estimating errors."""

    def __init__(self, name: str = 'covariance') -> None:
        super().__init__(name)

    def compute(self,
                statfunc: StatFunc,
                fitfunc: FitFunc,
                *,
                pars: np.ndarray,
                parmins: np.ndarray,
                parmaxes: np.ndarray,
                parhardmins: np.ndarray,
                parhardmaxes: np.ndarray,
                limit_parnums: np.ndarray,
                freeze_par: Callable,
                thaw_par: Callable,
                report_progress: Callable,
                get_par_name: Callable,
                statargs: Any = None,
                statkwargs: Any = None
                ) -> EstReturn:
        """Estimate the error range.

        .. versionchanged:: 4.17.1
           All arguments other than statfunc and fitfunc are now
           keyword-only.

        Notes
        -----
        The statargs and statkwargs arguments are currently unused.

        """

        if statargs is not None or statkwargs is not None:
            warning("statargs/kwargs set but values unused")

        stat_cb = StatCallback(statfunc)

        def fit_cb(scb, pars, parmins, parmaxes, i):
            # parameter i is a no-op usually
            return fitfunc(scb, pars, parmins, parmaxes)[2]

        # remin means reminimize -- *generally* not done (only
        # proj needs to reminimize).
        #
        # covar still needs to pass a reminimize flag to
        # get_one_sided_interval; a value less than zero is
        # interpreted as "never reminimize", so pass that value here.

        remin = -1.0
        tol = -1.0
        return covariance(pars,
                          parmins=parmins,
                          parmaxes=parmaxes,
                          parhardmins=parhardmins,
                          parhardmaxes=parhardmaxes,
                          sigma=self.sigma,
                          eps=self.eps,
                          tol=tol,
                          maxiters=self.maxiters,
                          remin=remin,
                          limit_parnums=limit_parnums,
                          stat_cb=stat_cb,
                          fit_cb=fit_cb,
                          report_progress=report_progress)



class Confidence(EstMethod):
    """The confidence method for estimating errors."""

    # defined pre-instantiation for pickling
    _added_config = {'remin': 0.01,
                     'fast': False,
                     'parallel': True,
                     'numcores': ncpus,
                     'maxfits': 5,
                     'max_rstat': 3,
                     'tol': 0.2,
                     'verbose': False,
                     'openinterval': False}

    def __init__(self, name: str = 'confidence') -> None:
        super().__init__(name)

        # Update EstMethod.config dict with Confidence specifics
        self.config.update(self._added_config)

    def compute(self,
                statfunc: StatFunc,
                fitfunc: FitFunc,
                *,
                pars: np.ndarray,
                parmins: np.ndarray,
                parmaxes: np.ndarray,
                parhardmins: np.ndarray,
                parhardmaxes: np.ndarray,
                limit_parnums: np.ndarray,  # integers
                freeze_par: Callable,
                thaw_par: Callable,
                report_progress: Callable,
                get_par_name: Callable,
                statargs: Any = None,
                statkwargs: Any = None
                ) -> EstReturn:

        """Estimate the error range.

        .. versionchanged:: 4.17.1
           All arguments other than statfunc and fitfunc are now
           keyword-only.
        Notes
        -----
        The statargs and statkwargs arguments are currently unused.


        """

        if statargs is not None or statkwargs is not None:
            warning("statargs/kwargs set but values unused")

        stat_cb = StatCallback(statfunc)

        def fit_cb(pars, parmins, parmaxes, i):
            # freeze model parameter i
            (current_pars,
             current_parmins,
             current_parmaxes) = freeze_par(pars, parmins, parmaxes, i)

            fit_pars = fitfunc(statfunc, current_pars,
                               current_parmins,
                               current_parmaxes)[1]
            # If stat is not chi-squared, and fit method is
            # lmdif, need to recalculate stat at end, just
            # like in sherpa/sherpa/fit.py:fit()
            stat = statfunc(fit_pars)[0]
            # stat = fitfunc(scb, pars, parmins, parmaxes)[2]
            # thaw model parameter i
            thaw_par(i)
            return stat

        #
        # convert stat call back to have the same signature as fit call back
        #
        def stat_cb_extra_args(fcn):
            def stat_cb_wrapper(x, *args):
                return fcn(x)
            return stat_cb_wrapper

        statcb = stat_cb_extra_args(stat_cb)
        if 1 == len(pars):
            fitcb = statcb
        else:
            fitcb = fit_cb

        return confidence(pars,
                          parmins=parmins,
                          parmaxes=parmaxes,
                          parhardmins=parhardmins,
                          parhardmaxes=parhardmaxes,
                          sigma=self.sigma,
                          eps=self.eps,
                          tol=self.tol,
                          maxiters=self.maxiters,
                          remin=self.remin,
                          verbose=self.verbose,
                          limit_parnums=limit_parnums,
                          stat_cb=statcb,
                          fit_cb=fitcb,
                          report_progress=report_progress,
                          get_par_name=get_par_name,
                          do_parallel=self.parallel,
                          numcores=self.numcores,
                          open_interval=self.openinterval)



class Projection(EstMethod):
    """The projection method for estimating errors."""

    # defined pre-instantiation for pickling
    _added_config = {'remin': 0.01,
                     'fast': False,
                     'parallel': True,
                     'numcores': ncpus,
                     'maxfits': 5,
                     'max_rstat': 3,
                     'tol': 0.2}

    def __init__(self, name: str = 'projection') -> None:
        super().__init__(name)

        # Update EstMethod.config dict with Projection specifics
        self.config.update(self._added_config)

    def compute(self,
                statfunc: StatFunc,
                fitfunc: FitFunc,
                *,
                pars: np.ndarray,
                parmins: np.ndarray,
                parmaxes: np.ndarray,
                parhardmins: np.ndarray,
                parhardmaxes: np.ndarray,
                limit_parnums: np.ndarray,  # integers
                freeze_par: Callable,
                thaw_par: Callable,
                report_progress: Callable,
                get_par_name: Callable,
                statargs: Any = None,
                statkwargs: Any = None
                ) -> EstReturn:
        """Estimate the error range.

        .. versionchanged:: 4.17.1
           All arguments other than statfunc and fitfunc are now
           keyword-only.

        Notes
        -----
        The statargs and statkwargs arguments are currently unused.

        """

        if statargs is not None or statkwargs is not None:
            warning("statargs/kwargs set but values unused")

        if fitfunc is None:
            raise TypeError("fitfunc should not be none")

        stat_cb = StatCallback(statfunc)

        def fit_cb(pars, parmins, parmaxes, i):
            # freeze model parameter i
            (current_pars,
             current_parmins,
             current_parmaxes) = freeze_par(pars, parmins, parmaxes, i)
            fit_pars = fitfunc(statfunc, current_pars,
                               current_parmins,
                               current_parmaxes)[1]
            # If stat is not chi-squared, and fit method is
            # lmdif, need to recalculate stat at end, just
            # like in sherpa/sherpa/fit.py:fit()
            stat = statfunc(fit_pars)[0]
            # stat = fitfunc(scb, pars, parmins, parmaxes)[2]
            # thaw model parameter i
            thaw_par(i)
            return stat

        return projection(pars,
                          parmins=parmins,
                          parmaxes=parmaxes,
                          parhardmins=parhardmins,
                          parhardmaxes=parhardmaxes,
                          sigma=self.sigma,
                          eps=self.eps,
                          tol=self.tol,
                          maxiters=self.maxiters,
                          remin=self.remin,
                          limit_parnums=limit_parnums,
                          stat_cb=stat_cb,
                          fit_cb=fit_cb,
                          report_progress=report_progress,
                          get_par_name=get_par_name,
                          do_parallel=self.parallel,
                          numcores=self.numcores)


def covariance(pars: np.ndarray,
               parmins: np.ndarray,
               parmaxes: np.ndarray,
               parhardmins: np.ndarray,
               parhardmaxes: np.ndarray,
               sigma: float,
               eps: float,
               tol: float,
               maxiters: int,
               remin: float,
               limit_parnums: np.ndarray,  # integers
               stat_cb: Callable,
               fit_cb: Callable,
               report_progress: Callable
               ) -> EstReturn:
    """Estimate errors using the covariance method."""

    # Do nothing with tol
    # Do nothing with report_progress (generally fast enough we don't
    # need to report back per-parameter progress)

    # Even though we only want limits on certain parameters, we have to
    # compute the matrix for *all* thawed parameters.  So we will do that,
    # and then pick the parameters of interest out of the result.

    try:
        info = _est_funcs.info_matrix(pars, parmins, parmaxes, parhardmins,
                                      parhardmaxes, sigma, eps, maxiters,
                                      remin, stat_cb)
    except EstNewMin as emin:
        # catch the EstNewMin exception and attach the modified
        # parameter values to the exception obj.  These modified
        # parvals determine the new lower statistic.
        raise EstNewMin(pars) from emin

    # Invert matrix, take its square root and multiply by sigma to get
    # parameter uncertainties; parameter uncertainties are the
    # diagonal elements of the matrix.

    # Use simpler matrix inversion function from numpy.  If that
    # doesn't work, assume it's an ill-conditioned or singular matrix,
    # and call pinv from numpy -- pinv will call the SVD function to
    # invert the matrix.  But call pinv *only* when inv is shown not
    # to work in a particular case -- use inv by default!

    # The reason for this is that pinv can give back very strange
    # results, when you don't *need* to use pinv, and it *also*
    # happens that the ratio between smallest and largest diagonal
    # elements approaches the machine precision for the data type.
    # The result is that errors that come from the largest diagonal
    # element are ludicrously small; you can't have a parameter value
    # of order 1.0, and an error of order 10^-30, for example.  The
    # simpler inv function for inverting matrices does not appear to
    # have the same issue.

    try:
        inv_info = np.linalg.inv(info)

    except LinAlgError:
        # catch the SVD exception and exit gracefully
        inv_info = np.zeros_like(info)
        inv_info[:] = np.nan

    except:
        try:
            inv_info = np.linalg.pinv(info)
        except LinAlgError:
            # catch the SVD exception and exit gracefully
            inv_info = np.zeros_like(info)
            inv_info[:] = np.nan

    diag = (sigma * np.sqrt(inv_info)).diagonal()

    # limit_parnums lists the indices of the array pars, that
    # correspond to the parameters of interest.  We will pick out
    # the diagonal elements corresponding to entries in limits_parnums,
    # and return only those bounds to the user.
    upper_bounds = []
    lower_bounds = []
    error_flags = []
    for num in limit_parnums:
        eflag = est_success
        ubound = diag[num]
        lbound = -diag[num]

        # What happens when lbound or ubound is NaN? This is
        # presumably why the code is written as it is below (e.g. a
        # pass if the values can be added to pars[num]).
        #
        if pars[num] + ubound < parhardmaxes[num]:
            pass
        else:
            ubound = np.nan
            eflag = est_hardmax
        if pars[num] + lbound > parhardmins[num]:
            pass
        else:
            lbound = np.nan
            if eflag == est_hardmax:
                eflag = est_hardminmax
            else:
                eflag = est_hardmin
        upper_bounds.append(ubound)
        lower_bounds.append(lbound)
        error_flags.append(eflag)

    return (np.array(lower_bounds), np.array(upper_bounds),
            np.array(error_flags), 0, inv_info)


def projection(pars: np.ndarray,
               parmins: np.ndarray,
               parmaxes: np.ndarray,
               parhardmins: np.ndarray,
               parhardmaxes: np.ndarray,
               sigma: float,
               eps: float,
               tol: float,
               maxiters: int,
               remin: float,
               limit_parnums: np.ndarray,  # integers
               stat_cb: Callable,
               fit_cb: Callable,
               report_progress: Callable,
               get_par_name: Callable,
               do_parallel: bool,
               numcores: int
               ) -> EstReturn:
    """Estimate errors using the projection method."""

    # Number of parameters to be searched on (*not* number of thawed
    # parameters, just number we are searching on)
    numsearched = len(limit_parnums)

    # _est_funcs.projection can be called on any subset of the thawed
    # parameters.  So we made a change here to call _est_funcs.projection
    # once per parameter number listed in limit_parnums, instead of
    # calling _est_funcs.projection once, with the original limit_parnums
    # array.  This way, we can report pack progress after the confidence
    # limit search is completed for each parameter, without descending
    # into the C++ code.
    #
    # It does mean we have to take apart the tuple returned by each call
    # to _est_funcs.projection; take the data we've pulled out, and
    # upon exiting the while loop, constructing a new tuple to return.
    # SMD 03/17/2009

    proj_func = _est_funcs.projection

    # LocalEstFunction
    def func(counter: int,  # unused
             singleparnum: int,
             lock: SupportsLock | None = None
             ) -> tuple[SupportsFloat, SupportsFloat, int, int, None]:
        try:
            singlebounds = proj_func(pars, parmins, parmaxes,
                                     parhardmins, parhardmaxes,
                                     sigma, eps, tol, maxiters,
                                     remin, [singleparnum], stat_cb,
                                     fit_cb)
        except EstNewMin as emin:
            # catch the EstNewMin exception and attach the modified
            # parameter values to the exception obj.  These modified
            # parvals determine the new lower statistic.
            raise EstNewMin(pars) from emin

        if lock is not None:
            lock.acquire()

        report_progress(singleparnum, singlebounds[0], singlebounds[1])

        if lock is not None:
            lock.release()

        return (singlebounds[0][0], singlebounds[1][0], singlebounds[2][0],
                singlebounds[3], None)

    if numsearched < 2 or not multi or numcores < 2:
        do_parallel = False

    if not do_parallel:
        lower_limits = np.zeros(numsearched)
        upper_limits = np.zeros(numsearched)
        eflags = np.zeros(numsearched, dtype=int)
        nfits = 0
        for i, pnum in enumerate(limit_parnums):
            singlebounds = func(i, pnum)
            lower_limits[i] = singlebounds[0]
            upper_limits[i] = singlebounds[1]
            eflags[i] = singlebounds[2]
            nfits = nfits + singlebounds[3]

        return (lower_limits, upper_limits, eflags, nfits, None)

    return parallel_est(func, limit_parnums, pars, numcores)

#################################confidence###################################


class ConfArgs:
    """The class ConfArgs is responsible for the arguments to the fit
    call back function."""

    def __init__(self,
                 xpars: ArrayType,
                 smin: ArrayType,
                 smax: ArrayType,
                 hmin: ArrayType,
                 hmax: ArrayType,
                 target_stat: SupportsFloat
                 ) -> None:
        self.ith_par = 0
        self.xpars = np.array(xpars, copy=True)
        self.slimit = (np.array(smin, copy=True),
                       np.array(smax, copy=True))
        self.hlimit = (np.array(hmin, copy=True),
                       np.array(hmax, copy=True))
        self.target_stat = target_stat

    def __call__(self
                 ) -> tuple[int, np.ndarray, np.ndarray, np.ndarray, SupportsFloat]:
        return (self.ith_par, self.xpars, self.slimit, self.hlimit,
                self.target_stat)

    def __str__(self) -> str:
        a2s = np.array2string
        msg = ''
        msg += '# smin = ' + a2s(self.slimit[0], precision=6) + '\n'
        msg += '# smax = ' + a2s(self.slimit[1], precision=6) + '\n'
        msg += '# hmin = ' + a2s(self.hlimit[0], precision=6) + '\n'
        msg += '# hmax = ' + a2s(self.hlimit[1], precision=6) + '\n#\n'
        msg += '# Note: for the intermediate steps, the notation:\n'
        msg += '         par.name -/+: f( x ) = stat\n'
        msg += '# ==> `stat` is the statistic when parameter `par.name` is frozen at `x`\n'
        msg += '# while searching for the `lower/upper` confidence level, respectively.\n#'
        return msg

    def get_par(self) -> SupportsFloat:
        """return the current (worked on) par"""
        return self.xpars[self.ith_par]

    def get_hlimit(self, dir: int) -> SupportsFloat:
        """ return the current (worked on) hard limit"""
        return self.hlimit[dir][self.ith_par]

    def get_slimit(self, dir: int) -> SupportsFloat:
        """ return the current (worked on) soft limit"""
        return self.slimit[dir][self.ith_par]


class ConfBlog:

    def __init__(self, blogger, prefix, verbose, lock) -> None:
        self.blogger = blogger
        self.prefix = prefix
        self.verbose = verbose
        self.lock = lock


class Limit:
    """Represent a limit.

    This should not be used directly: use `LowerLimit` or `UpperLimit`
    instead.

    """

    def __init__(self, limit: SupportsFloat) -> None:
        self.limit = limit


class LowerLimit(Limit):
    "A lower limit"

    def __str__(self) -> str:
        return f'LowerLimit: limit={self.limit:e}'

    def is_beyond_limit(self, x) -> bool:
        return x < self.limit


class UpperLimit(Limit):
    "An upper limit"

    def __str__(self) -> str:
        return f'UpperLimit: limit={self.limit:e}'

    def is_beyond_limit(self, x: SupportsFloat) -> bool:
        # The float calls are for typing checks.
        return float(x) > float(self.limit)


class ConfBracket:
    """The class ConfBracket is responsible for bracketing the root within
    the interval (a,b) where f(a)*f(b) < 0.0"""

    neg_pos = (-1, 1)

    def __init__(self, myargs, trial_points) -> None:
        self.myargs = myargs
        self.trial_points = trial_points
        self.fcn: Callable | None = None

    # TODO: rename the dir and iter arguments
    def __call__(self,
                 dir: int,
                 iter,
                 step_size,
                 open_interval,
                 maxiters,
                 tol,
                 bloginfo: ConfBlog
                 ) -> ConfRootNone:

        #
        # Either 1) a root has been found (ConfRootZero), 2) an interval
        # where the root has been confined (ConfRootBracket) or 3) No
        # possible chance for a root (ConfRootNone), ie by trying points
        # upto/beyond the hard limit and no chance for a root has been found.
        #
        return self.find(dir, iter, step_size, open_interval,
                         maxiters, tol, bloginfo)

    # TODO: rename the dir and iter arguments
    # TODO: the iter argument gets over-ridden below, so do we need it? (and rename it)
    def find(self,
             dir: int,
             iter,
             step_size,
             open_interval,
             maxiters, tol,
             bloginfo: ConfBlog,
             base=2.0
             ) -> ConfRootNone:

        assert self.fcn is not None, 'callback func has not been set'

        hlimit = [LowerLimit(self.myargs.get_hlimit(dir)),
                  UpperLimit(self.myargs.get_hlimit(dir))]
        slimit = [LowerLimit(self.myargs.get_slimit(dir)),
                  UpperLimit(self.myargs.get_slimit(dir))]

        xxx = self.trial_points[0]
        fff = self.trial_points[1]

        conf_step = ConfStep(xxx, fff)

        mymaxiters = min(maxiters, 16)

        plateau = 0
        max_plateau_iter = 5

        try:

            # TODO: rename iter argument
            for iter in range(mymaxiters):

                if 0 == iter:
                    x = conf_step.covar(dir, iter, step_size, base)
                elif 1 == iter:
                    x = conf_step.secant(dir, iter, step_size, base)
                    # x = conf_step.covar( dir, iter, step_size, base )
                else:
                    x = conf_step.quad(dir, iter, step_size, base, bloginfo)

                if x is None or np.isnan(x):
                    return ConfRootNone()

                # Make sure x is not beyond the **hard** limit
                if hlimit[dir].is_beyond_limit(x):

                    x = hlimit[dir].limit
                    f = self.fcn(x, self.myargs())
                    # print 'find(): beyond hard limit: f(%.14e)=%.14e' % (x,f)
                    if abs(f) <= tol:
                        return ConfRootZero(x)

                    if f >= 0.0:
                        return ConfRootBracket(self.fcn, self.trial_points,
                                               open_interval)

                    return ConfRootNone()

                if slimit[dir].is_beyond_limit(x):

                    f = self.fcn(x, self.myargs())
                    # print 'find(): beyond soft limit: f(%.14e)=%.14e' % (x,f)
                    if abs(f) <= tol:
                        return ConfRootZero(x)

                    if f >= 0.0:
                        return ConfRootBracket(self.fcn, self.trial_points,
                                               open_interval)

                    if f < fff[-2]:
                        # if the fit beyond the soft limit is a better fit
                        # then the confidence for the parameter does not exist
                        return ConfRootNone()

                else:

                    f = self.fcn(x, self.myargs())
                    # print 'find(): f(%.14e)=%.14e' % (x,f)
                    if abs(f) <= tol:
                        return ConfRootZero(x)

                    if f >= 0.0:
                        return ConfRootBracket(self.fcn, self.trial_points,
                                               open_interval)

                if Knuth_close(fff[-2], fff[-1], 1.0e-6):
                    plateau += 1
                    if plateau > max_plateau_iter:
                        # print 'find( %d ): plateau = %d', (iter,plateau)
                        return ConfRootNone(None)
                # else:
                #    if plateau > 0:
                #        plateau -= 1

            return ConfRootNone(None)

        except OutOfBoundErr:
            return ConfRootNone()


class ConfRootNone:
    """The base class for the root of the confidence interval"""

    def __init__(self,
                 root: SupportsFloat | None = None
                 ) -> None:
        """If self.root == None, then
        1) points up to the hard limits were tried and it was not possible
        to bracketed the solution.
        2) a parameter beyond the soft limit has been tried and the new stat
        was found to be **less** then the initial minimum"""
        self.root = root

    def __call__(self,
                 tol,
                 bloginfo: ConfBlog
                 ) -> SupportsFloat | None:
        return self.root

    def __str__(self) -> str:
        return 'No possible root exist'


class ConfRootBracket(ConfRootNone):
    """The class contains the bracket where the confidence root has
    been bracketed, ie where f(a)*f(b) < 0"""

    def __init__(self,
                 fcn,
                 trial_points,
                 open_interval
                 ) -> None:
        super().__init__(root=None)
        self.fcn = fcn
        self.trial_points = trial_points
        self.open_interval = open_interval

    def __call__(self,
                 tol,
                 bloginfo: ConfBlog
                 ) -> SupportsFloat | tuple[SupportsFloat, SupportsFloat] | None:

        def warn_user_about_open_interval(listval):

            if bloginfo.lock is not None:
                bloginfo.lock.acquire()

            if 0 == bloginfo.verbose:
                prefix = '%s ' % bloginfo.prefix.lstrip()
            else:
                prefix = '%s ' % bloginfo.prefix

            interval = list_to_open_interval(listval)
            bloginfo.blogger.info(prefix +
                                  'WARNING: The confidence level lies within '
                                  + interval)

            if bloginfo.lock is not None:
                bloginfo.lock.release()

        xxx = self.trial_points[0]
        fff = self.trial_points[1]

        if np.sign(fff[-2]) == np.sign(fff[-1]):
            self.root = None
            return None

        answer = zeroin(self.fcn, xxx[-2], xxx[-1], fa=fff[-2],
                        fb=fff[-1], maxfev=32, tol=tol)
        if abs(answer[0][1]) > tol:
            xafa = answer[1][0]
            xa = xafa[0]
            fa = xafa[1]
            xbfb = answer[1][1]
            xb = xbfb[0]
            fb = xbfb[1]
            if np.sign(fa) != np.sign(fb):
                if not self.open_interval:
                    warn_user_about_open_interval([xa, xb])
                    return (xa + xb) / 2.0

                return min(xa, xb), max(xa, xb)

            return None

        self.root = answer[0][0]
        return self.root

    def __str__(self) -> str:
        msg = 'root is within the interval ( f(%e)=%e, f(%e)=%e )' \
            % (self.trial_points[0][-2], self.trial_points[1][-2],
               self.trial_points[0][-1], self.trial_points[1][-1], )
        return msg


class ConfRootZero(ConfRootNone):
    """The class with the root/zero of the confidence interval"""

    def __str__(self) -> str:
        return f'root = {self.root:e}'


class ConfStep:

    def __init__(self, xtrial, ftrial) -> None:
        self.xtrial = xtrial
        self.ftrial = ftrial

    def covar(self, dir, iter, stepsize, base):
        return self.xtrial[-1] + \
            ConfBracket.neg_pos[dir] * pow(base, iter) * stepsize
        # return self.xtrial[ 0 ] + \
        #       ConfBracket.neg_pos[ dir ] * pow( base, iter ) * stepsize

    def Halley(self, coeffs, x, maxfev=8, tol=1.0e-3):

        for nfev in range(maxfev):

            ax = coeffs[0] * x
            fval = (ax + coeffs[1]) * x + coeffs[2]
            if abs(fval) <= tol:
                return [x, fval]

            fdif = 2.0 * ax + coeffs[1]
            fdif2 = 2.0

            numer = 2.0 * fval * fdif
            denom = 2.0 * fdif * fdif - fval * fdif2
            x -= numer / denom
            nfev += 1

        return [x, fval]

    def is_same_dir(self, dir, current_pos, proposed_pos) -> bool:
        delta = proposed_pos - current_pos
        return np.sign(delta) == ConfBracket.neg_pos[dir]

    def quad(self, dir, iter, step_size, base, bloginfo):

        coeffs = quad_coef(self.xtrial[-3:], self.ftrial[-3:])
        delta = ConfBracket.neg_pos[dir]
        delta *= abs(self.xtrial[-1] - self.xtrial[-2])
        lastx = self.xtrial[-1]

        mroot = demuller(np.poly1d(coeffs), lastx + delta,
                         lastx + 2 * delta, lastx + 3 * delta,
                         tol=1.0e-2)

        xroot = mroot[0][0]
        if xroot is None or np.isnan(xroot):
            return self.covar(dir, iter, step_size, base)

        try:
            [xroot, froot] = self.Halley(coeffs, xroot, tol=1.0e-3)
        except ZeroDivisionError:
            xroot = None

        if xroot is not None and not np.isnan(xroot) and \
           self.is_same_dir(dir, self.xtrial[-1], xroot):
            return xroot

        return self.covar(dir, iter, step_size, base)

    def secant(self, dir, iter, step_size, base):

        xb = self.xtrial[-2]
        fb = self.ftrial[-2]
        xa = self.xtrial[-1]
        fa = self.ftrial[-1]

        if abs(fb) > abs(fa) or 0.0 == fa:
            return self.covar(dir, iter, step_size, base)

        s = fb / fa
        p = (xa - xb) * s
        if 1.0 == s:
            return self.covar(dir, iter, step_size, base)
        q = 1.0 - s
        x = xb - p / q

        if self.is_same_dir(dir, xa, x):
            return x

        return self.covar(dir, iter, step_size, base)


def confidence(pars: np.ndarray,
               parmins: np.ndarray,
               parmaxes: np.ndarray,
               parhardmins: np.ndarray,
               parhardmaxes: np.ndarray,
               sigma: float,
               eps: float,
               tol: float,
               maxiters: int,
               remin: float,
               verbose: bool,      # different to covariance/projection
               limit_parnums: np.ndarray,  # integers
               stat_cb: Callable,
               fit_cb: Callable,
               report_progress: Callable,
               get_par_name: Callable,
               do_parallel: bool,
               numcores: int,
               open_interval
               ) -> EstReturn:

    def get_prefix(index, name, minus_plus):
        '''To print the prefix/indent when verbose is on'''
        blank = 3 * index * ' '
        return [f"{blank}{name} {mtext}:" for mtext in minus_plus]

    def get_delta_root(arg, dir, par_at_min):

        my_neg_pos = ConfBracket.neg_pos[dir]

        if is_iterable(arg):
            # return map( lambda x: my_neg_pos * abs( x - par_at_min ), arg )
            return arg

        if arg is not None:
            arg -= par_at_min
            return my_neg_pos * abs(arg)

        return arg

    def get_step_size(error_scales, upper_scales, index, par):

        if 0 != error_scales[index]:
            # if covar error is NaN then set it to fraction of the par value.
            ith_covar_err = 0.0625 * abs(par)
        else:
            ith_covar_err = abs(upper_scales[index])
        if 0.0 == ith_covar_err:
            # just in case covar and/or par is 0
            ith_covar_err = 1.0e-6

        return ith_covar_err

    def monitor_func(fcn, history):
        def myfunc(x, *args):
            fval = fcn(x, *args)
            history[0].append(x)
            history[1].append(fval)
            return fval

        return myfunc

    def print_status(myblog, verbose, prefix, answer, lock):

        if lock is not None:
            lock.acquire()

        if 0 == verbose:
            msg = '%s\t' % prefix.lstrip()
        else:
            msg = '%s\t' % prefix

        if is_iterable(answer):
            msg += list_to_open_interval(answer)
        elif answer is None:
            msg += '-----'
        else:
            msg += '%g' % answer
        myblog(msg)

        if lock is not None:
            lock.release()

    #
    # Work in the translated coordinate. Hence the 'errors/confidence'
    # are the zeros/roots in the translated coordinate system.
    #
    def translated_fit_cb(fcn, myargs):
        def translated_fit_cb_wrapper(x, *args):
            hlimit = myargs.hlimit
            slimit = myargs.slimit
            hmin = hlimit[0]
            hmax = hlimit[1]
            xpars = myargs.xpars
            ith_par = myargs.ith_par
            # The parameter must be within the hard limits
            if x < hmin[ith_par] or x > hmax[ith_par]:
                raise OutOfBoundErr

            smin = slimit[0]
            smax = slimit[1]
            orig_ith_xpar = xpars[ith_par]
            xpars[ith_par] = x
            translated_stat = fcn(
                xpars, smin, smax, ith_par) - myargs.target_stat
            xpars[ith_par] = orig_ith_xpar
            return translated_stat
        return translated_fit_cb_wrapper

    def verbose_fitcb(fcn, bloginfo):
        if 0 == bloginfo.verbose:
            return fcn

        def verbose_fcn(x, *args):
            fval = fcn(x, *args)
            msg = '%s f( %e ) =' % (bloginfo.prefix, x)
            if fval is None:
                msg = '%s None' % msg
            else:
                msg = '%s %e' % (msg, fval)
            bloginfo.blogger.info(msg)
            return fval
        return verbose_fcn

    sherpablog = logging.getLogger('sherpa')  # where to print progress report

    # Get minimum fit statistic, and calculate target statistic value
    orig_min_stat = stat_cb(pars)
    delta_stat = sigma * sigma
    target_stat = orig_min_stat + delta_stat

    lower_scales = None
    upper_scales = None
    error_scales = None
    nfits = 0
    results = None

    try:
        (lower_scales, upper_scales, error_scales, nfits,
         results) = covariance(pars, parmins, parmaxes, parhardmins,
                               parhardmaxes, 1.0, eps, tol, maxiters,
                               remin, limit_parnums, stat_cb,
                               fit_cb, report_progress)
    except EstNewMin as e:
        raise e
    except:
        error_scales = np.full(len(pars), est_hardminmax)

    myargs = ConfArgs(pars, parmins, parmaxes, parhardmins, parhardmaxes,
                      target_stat)

    if 0 != verbose:
        msg = '#\n# f' + np.array2string(np.asarray(pars), precision=6)
        msg += ' = %e\n' % orig_min_stat
        msg += '# sigma = %e\n' % sigma
        msg += '# target_stat = %e\n' % target_stat
        msg += '# tol = %e\n' % eps
        msg += '%s' % myargs
        sherpablog.info(msg)

    # LocalEstFunc
    def func(counter: int,
             singleparnum: int,
             lock: SupportsLock | None = None
             ) -> tuple[SupportsFloat, SupportsFloat, int, int, None]:

        counter_cb = FuncCounter(fit_cb)

        #
        # These are the bounds to be returned by this method
        #
        conf_int = [[], []]
        error_flags = []

        #
        # If the user has requested a specific parameter to be
        # calculated then 'ith_par' represents the index of the
        # free parameter to deal with.
        #
        myargs.ith_par = singleparnum

        fitcb = translated_fit_cb(counter_cb, myargs)

        par_name = get_par_name(myargs.ith_par)

        ith_covar_err = get_step_size(error_scales, upper_scales, counter,
                                      pars[myargs.ith_par])

        trial_points = [[], []]
        fitcb = monitor_func(fitcb, trial_points)

        bracket = ConfBracket(myargs, trial_points)

        # the parameter name is set, may as well get the prefix
        prefix = get_prefix(counter, par_name, ['-', '+'])

        myfitcb = [verbose_fitcb(fitcb,
                                 ConfBlog(sherpablog, prefix[0], verbose, lock)),
                   verbose_fitcb(fitcb,
                                 ConfBlog(sherpablog, prefix[1], verbose, lock))]

        for dirn in range(2):

            # trial_points stores the history of the points for the
            # parameter which has been evaluated in order to locate
            # the root. Note the first point is 'given' since the info
            # of the minimum is crucial to the search.
            #
            bracket.trial_points[0].append(pars[myargs.ith_par])
            bracket.trial_points[1].append(- delta_stat)

            myblog = ConfBlog(sherpablog, prefix[dirn], verbose, lock)

            # have to set the callback func otherwise disaster.
            bracket.fcn = myfitcb[dirn]
            root = bracket(dirn, iter, ith_covar_err, open_interval, maxiters,
                           eps, myblog)

            myzero = root(eps, myblog)

            delta_zero = get_delta_root(myzero, dirn, pars[myargs.ith_par])

            conf_int[dirn].append(delta_zero)

            status_prefix = get_prefix(counter, par_name, ['lower bound',
                                                           'upper bound'])
            print_status(myblog.blogger.info, verbose, status_prefix[dirn],
                         delta_zero, lock)

        # This should really set the error flag appropriately.
        error_flags.append(est_success)

        return (conf_int[0][0], conf_int[1][0], error_flags[0],
                counter_cb.nfev, None)

    if len(limit_parnums) < 2 or not multi or numcores < 2:
        do_parallel = False

    if not do_parallel:
        lower_limits = []
        upper_limits = []
        eflags = []
        nfits = 0
        for i, lpar in enumerate(limit_parnums):
            lower_limit, upper_limit, flags, nfit, extra = func(
                i, lpar)
            lower_limits.append(lower_limit)
            upper_limits.append(upper_limit)
            eflags.append(flags)
            nfits += nfit

        return (lower_limits, upper_limits, eflags, nfits, None)

    return parallel_est(func, limit_parnums, pars, numcores)

#################################confidence###################################


class LocalEstFunc(Protocol):
    """Process a single parameter."""

    def __call__(self,
                 counter: int,
                 singleparnum: int,
                 lock: SupportsLock | None = None
                 ) -> tuple[SupportsFloat, SupportsFloat, int, int, None]:
        ...


def parallel_worker(estfunc: LocalEstFunc,
                    pars: np.ndarray,
                    out_q: SupportsQueue,
                    err_q: SupportsQueue,
                    lock: SupportsLock,
                    parids: np.ndarray,
                    parnums: np.ndarray
                    ) -> None:
    """Estimate errors for a set of parameters.

    The return value is sent back via the out_q parameter.

    """

    results = []
    for parid, singleparnum in zip(parids, parnums):
        try:
            result = estfunc(parid, singleparnum, lock)
            results.append((parid, result))
        except EstNewMin:
            # catch the EstNewMin exception and include the exception
            # class and the modified parameter values to the error queue.
            # These modified parvals determine the new lower statistic.
            # The exception class will be re-raised with the
            # parameter values attached.  C++ Python exceptions are not
            # picklable for use in the queue.
            err_q.put(EstNewMin(pars))
            return
        except Exception as e:
            err_q.put(e)
            return

    out_q.put(results)


def parallel_est(estfunc: LocalEstFunc,
                 limit_parnums: np.ndarray,  # integers
                 pars: np.ndarray,
                 numcores: int = ncpus
                 ) -> EstReturn:
    """Run a function on a sequence of inputs in parallel.

    A specialized version of sherpa.utils.parallel.parallel_map.

    Parameters
    ----------
    estfunc : function
       This function accepts three arguments and returns a value.
    limit_parnums : sequence
    pars : sequence
       The current parameter values
    numcores : int, optional
       The number of calls to ``function`` to run in parallel. When
       set to ``None``, all the available CPUs on the machine - as
       set either by the 'numcores' setting of the 'parallel' section
       of Sherpa's preferences or by `multiprocessing.cpu_count` - are
       used.

    Returns
    -------
    ans : tuple

    """

    # Help the type checkers, as this should not be called when
    # multiprocessing (hence context) is not available.
    #
    assert context is not None

    # See sherpa.utils.parallel for a discussion of how multiprocessing is
    # being used to run code in parallel.
    #
    manager = context.Manager()
    out_q = manager.Queue()
    err_q = manager.Queue()

    # Unlike sherpa.utils.parallel.parallel_map, these routines do use a lock
    # for serializing screen output.
    #
    lock = manager.Lock()

    size = len(limit_parnums)
    parids = np.arange(size)

    # if len(limit_parnums) is less than numcores, only use length number of
    # processes
    if size < numcores:
        numcores = size

    # group limit_parnums into numcores-worth of chunks
    limit_chunk = np.array_split(limit_parnums, numcores)
    parids_chunk = np.array_split(parids, numcores)

    tasks = [context.Process(target=parallel_worker,
                             args=(estfunc, pars, out_q, err_q, lock,
                                   parid, parnum))
             for parid, parnum in zip(parids_chunk, limit_chunk)]

    return run_tasks(tasks, out_q, err_q, size)


def run_tasks(tasks: Sequence[SupportsProcess],
              out_q: SupportsQueue,
              err_q: SupportsQueue,
              size: int
              ) -> EstReturn:
    """Run the processes, exiting early if necessary, and return the results.

    A specialized version of sherpa.utils.parallel.run_tasks (note the
    different order of the queues).

    Parameters
    ----------
    procs : list of multiprocessing.Process tasks
        The processes to run.
    out_q, err_q : manager.Queue
        The success and error and channels used by the processes.
    size : int
        The number of results being returned (can be larger than len(procs)).

    Returns
    -------
    result : tuple

    """

    process_tasks(tasks, err_q)

    lower_limits = size * [None]
    upper_limits = size * [None]
    eflags = size * [None]
    nfits = 0

    while not out_q.empty():
        for parid, singlebounds in out_q.get():
            # Have to guarantee that the tuple returned by projection
            # is always (array, array, array, int) for this to work.
            lower_limits[parid] = singlebounds[0]
            upper_limits[parid] = singlebounds[1]
            eflags[parid] = singlebounds[2]
            nfits += singlebounds[3]

    return (lower_limits, upper_limits, eflags, nfits, None)
