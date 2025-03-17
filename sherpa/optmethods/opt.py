#
#  Copyright (C) 2019 - 2021, 2023 - 2025
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

from collections.abc import Callable, Sequence
from typing import Protocol

import numpy as np

from sherpa.utils import Knuth_close, FuncCounter
from sherpa.utils.parallel import SupportsQueue, \
    multi, context, run_tasks, create_seeds
from sherpa.utils.random import RandomType, uniform
from sherpa.utils.types import ArrayType


__all__ = ('Opt', 'MyNcores', 'SimplexRandom', 'SimplexNoStep',
           'SimplexStep')


FUNC_MAX = float(np.finfo(np.float64).max)

# Some code assumes the first argument is an ndarray, but is this always
# the case? For now assume that the data is always sent as an ndarray.
#
# sherpa.stats.StatCallback converts a StatFunc into an OptimizerFunc.
#
OptimizerFunc = Callable[[np.ndarray], float]

MyOptOutput = tuple[int, float, np.ndarray]

class WorkerFunc(Protocol):

    def __call__(self,
                 func: OptimizerFunc,
                 x: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float,
                 maxnfev: int | None,
                 rng: RandomType | None = None
                 ) -> MyOptOutput:
        ...


class MyNcores:
    """Support distributed processing."""

    def __init__(self) -> None:
        if multi is False:
            raise TypeError("multicores not available")

    def calc(self,
             funcs: Sequence[WorkerFunc],
             numcores: int,  # TODO: this is currently unused
             fcn: OptimizerFunc,
             x: np.ndarray,
             xmin: np.ndarray,
             xmax: np.ndarray,
             tol: float,
             maxnfev: int | None,
             rng: RandomType | None = None
             ) -> list[MyOptOutput]:
        """Apply each function to the arguments, running in parallel."""

        nfuncs = len(funcs)
        if nfuncs == 0:
            raise TypeError("funcs can not be empty")

        for func in funcs:
            if not callable(func):
                raise TypeError(f"input func '{repr(func)}' is not callable")

        # See sherpa.utils.parallel for the logic used here.
        # At this point we can assume that context is not None,
        # since multi is True.
        #
        # Should this shard by numcores?
        #
        assert context is not None
        manager = context.Manager()
        out_q = manager.Queue()
        err_q = manager.Queue()

        seeds = create_seeds(rng, nfuncs)

        procs = [context.Process(target=self.my_worker,
                                 args=(func, ii, out_q, err_q,
                                       fcn, x, xmin, xmax, tol, maxnfev),
                                 kwargs={"rng": np.random.default_rng(seed)}
                                 )
                 for ii, (func, seed) in enumerate(zip(funcs, seeds))]

        return run_tasks(procs, err_q, out_q)

    def my_worker(self,
                  opt: WorkerFunc,
                  idval: int,
                  out_q: SupportsQueue[tuple[int, list[MyOptOutput]]],
                  err_q: SupportsQueue[Exception],
                  fcn: OptimizerFunc,
                  x: np.ndarray,
                  xmin: np.ndarray,
                  xmax: np.ndarray,
                  tol: float,
                  maxnfev: int | None,
                  rng: RandomType | None = None
                  ) -> None:
        raise NotImplementedError("my_worker has not been implemented")


class Opt:
    """Base optimisation class.

    .. versionchanged:: 4.17.0
       The class structure has been changed (e.g. `nfev` is now a
       scalar and not a single-element list).

    """

    def __init__(self,
                 func: OptimizerFunc,
                 xmin: ArrayType,
                 xmax: ArrayType
                 ) -> None:
        self.xmin = np.asarray(xmin)
        self.xmax = np.asarray(xmax)
        self.npar = len(self.xmin)
        if len(self.xmax) != self.npar:
            raise ValueError("xmin and xmax must be the same size")

        # The function count should only be done on "valid" ranges,
        # hence the ordering here.
        #
        self.func_count = FuncCounter(func)
        self.func = self.func_bounds(self.func_count, self.npar,
                                     self.xmin, self.xmax)

    # The sub-classes have different arguments, but do have maxnfev
    # and ftol as common values.
    #
    def __call__(self,
                 maxnfev: int,
                 ftol: float,
                 *args,
                 **kwargs) -> MyOptOutput:
        raise NotImplementedError

    @property
    def nfev(self) -> int:
        """How many evaluations of the function have been made?

        This does not count proposed values which were outside the
        parameter limits.

        """
        return self.func_count.nfev

    def _outside_limits(self,
                        x: np.ndarray,
                        xmin: np.ndarray,
                        xmax: np.ndarray
                        ) -> bool:
        return bool(np.any(x < xmin) or np.any(x > xmax))

    # We should be able to take these parameters from the class, or
    # re-write this logic.
    #
    def func_bounds(self,
                    func: OptimizerFunc,
                    npar: int,
                    xmin: ArrayType | None = None,
                    xmax: ArrayType | None = None
                    ) -> OptimizerFunc:
        """In order to keep the current number of function evaluations:
        func_counter should be called before func_bounds. For example,
        the following code
          x0 = [-1.2, 1.0]
          xmin = [-10.0, -10.0]
          xmax = [10.0, 10.0]
          nfev, rosen = func_counter(Rosenbrock)
          rosenbrock = func_bounds(rosen, xmin, xmax)
          print rosenbrock([-15.0, 1.0]), nfev[0]
        should output:
        inf 0"""
        if xmin is not None and xmax is not None:
            xmin = np.asarray(xmin)
            xmax = np.asarray(xmax)
        elif xmin is not None:
            xmin = np.asarray(xmin)
            xmax = np.asarray([np.inf for ii in xmin])
        elif xmax is not None:
            xmax = np.asarray(xmax)
            xmin = np.asarray([- np.inf for ii in xmax])
        else:
            xmin = np.asarray([- np.inf for ii in range(npar)])
            xmax = np.asarray([np.inf for ii in range(npar)])

        def func_bounds_wrapper(x):
            if self._outside_limits(x, xmin, xmax):
                return FUNC_MAX

            return func(x)

        return func_bounds_wrapper


# The random state is, as of 4.17.1, now set when calling __init__:
# either rng is set or it is created based on the seed argument.  The
# self.rng field is then used to access random numbers, and the seed
# argument in other routines, such as init, is ignored.
#
# The simplex field is a npop by npar array, where each row
# contains the parameter values for a given index/pop.
#
# The __get/setitem__ calls allow the object to index into
# the pars array by pop number. The statistic needs to be
# set/get manually via the fctvals field.
#
class SimplexBase:
    """

    .. versionchanged:: 4.17.1
       The init routine has been reworked and is now sent the starting
       simplex and the statistic value is now stored separately from
       the simplex. The seed argument is now ignored other than in the
       initialization code, where it is only used if the rng argument
       is not set.  Some of the arguments must now be set by name.

    """

    def __init__(self,
                 func: OptimizerFunc,
                 npop: int,
                 xpar: ArrayType,
                 xmin: ArrayType,
                 xmax: ArrayType,
                 *,
                 step: np.ndarray | None,
                 seed: int | None,
                 factor: float | None,
                 rng: RandomType | None = None
                 ) -> None:

        assert npop > 0  # safety check

        self.func = func
        self.xmin = np.asarray(xmin)
        self.xmax = np.asarray(xmax)
        self.npar = len(xpar)
        if rng is None:
            self.rng = np.random.default_rng(seed)
        else:
            self.rng = rng

        xpar_np = np.asarray(xpar)

        simplex = np.empty((npop, self.npar))
        simplex[0] = xpar_np
        simplex = self.init(npop=npop, xpar=xpar_np, simplex=simplex,
                            step=step, factor=factor)
        self.simplex, self.fctvals = self.eval_simplex(npop, simplex)

    def __getitem__(self, index):
        return self.simplex[index]

    def __setitem__(self, index, val) -> None:
        self.simplex[index] = val

    def calc_centroid(self) -> np.ndarray:
        return np.mean(self.simplex[:-1, :], 0)

    def check_convergence(self,
                          ftol: float,
                          method: int
                          ) -> bool:

        def are_func_vals_close_enough() -> bool:
            return Knuth_close(self.fctvals[0], self.fctvals[-1],
                               ftol)

        def is_fct_stddev_small_enough() -> bool:
            fval_std = np.std(self.fctvals)
            # Force float comparison to avoid type checker warnings, both for
            # the arguments to the comparison and the return value.
            return float(fval_std) < float(ftol)

        def is_max_length_small_enough() -> bool:
            """

               max  || x  - x  || <= tol max(1.0, || x ||)
                        i    0                         0

            where 1 <= i <= n"""

            max_xi_x0 = -1.0
            x0 = self.simplex[0]

            npar_plus_1 = len(x0) + 1
            for ii in range(1, npar_plus_1):
                xi = self.simplex[ii]
                xi_x0 = xi - x0
                max_xi_x0 = max(max_xi_x0, np.dot(xi_x0, xi_x0))

            # force float conversion to avoid type checker warnings
            return max_xi_x0 <= float(ftol) * max(1.0, np.dot(x0, x0))

        if 0 == method:
            if is_max_length_small_enough():
                return True
        elif 2 == method:
            if not is_max_length_small_enough():
                return False
            stddev = is_fct_stddev_small_enough()
            fctval = are_func_vals_close_enough()
            return stddev and fctval
        else:
            if not is_max_length_small_enough():
                return False
            stddev = is_fct_stddev_small_enough()
            fctval = are_func_vals_close_enough()
            return stddev or fctval

        return False

        # TODO: what is this code meant to be doing as it is unreachable?
        #
        # It has therefore been commented out.
        #
        # num = 2.0 * abs(self.simplex[0, -1] - self.simplex[-1, -1])
        # denom = abs(self.simplex[0, -1]) + abs(self.simplex[-1, -1]) + 1.0
        # if num / denom > ftol:
        #     return False
        #
        # func_vals = [col[-1] for col in self.simplex]
        # if np.std(func_vals) > ftol:
        #     return False
        #
        # return True

    def eval_simplex(self,
                     npop: int,
                     simplex: np.ndarray
                     ) -> tuple[np.ndarray, np.ndarray]:

        fctvals = np.zeros(npop)
        for ii in range(npop):
            fctvals[ii] = self.func(simplex[ii])

        return self.sort_me(simplex, fctvals)

    def init(self,
             *,
             npop: int,
             xpar: np.ndarray,
             simplex: np.ndarray,
             step: np.ndarray | None,
             factor: float | None
             ) -> np.ndarray:
        """Initialize the class.

        .. versionchanged:: 4.17.1
           The arguments must now all be given by name and the seed
           argument has been removed.

        """
        raise NotImplementedError("init has not been implemented")

    # This mutates the input simplex argument (and returns it).
    #
    def init_random_simplex(self,
                            xpar: np.ndarray,
                            simplex: np.ndarray,
                            *,
                            start: int,
                            npop: int,
                            factor: float | None
                            ) -> np.ndarray:
        """Initialize the simplex.

        .. versionchanged:: 4.17.1
           The seed argument has been removed. Most of the arguments
           must now be set by name.

        """

        # This could be done before changing the seed, but code may
        # require the current behavior.
        #
        if start >= npop:
            return simplex

        # At this point it is assumed the caller has set factor.
        #
        assert factor is not None, "caller did not set factor"
        deltas = factor * np.abs(np.asarray(xpar))
        xmins = np.maximum(self.xmin, xpar - deltas)
        xmaxs = np.minimum(self.xmax, xpar + deltas)
        for ii in range(start, npop):
            simplex[ii] = uniform(self.rng, xmins, xmaxs)

        return simplex

    def move_vertex(self,
                    centroid: np.ndarray,
                    coef: float
                    ) -> tuple[np.ndarray, float]:
        vertex = (1.0 + coef) * centroid - coef * self.simplex[self.npar]
        return vertex, self.func(vertex)

    def shrink(self, shrink_coef: float) -> None:

        self.simplex[1:] = self.simplex[0] + \
            shrink_coef * (self.simplex[1:] - self.simplex[0])

        for idx, simplex in enumerate(self.simplex[1:], 1):
            self.fctvals[idx] = self.func(simplex)

    @staticmethod
    def sort_me(simp: np.ndarray,
                fctvals: np.ndarray
                ) -> tuple[np.ndarray, np.ndarray]:
        """Reorder the simplex by the statistic value (low to high)"""

        idx = np.argsort(fctvals)
        return simp[idx], fctvals[idx]

    def sort(self) -> None:
        self.simplex, self.fctvals = self.sort_me(self.simplex,
                                                  self.fctvals)


class SimplexNoStep(SimplexBase):

    def init(self,
             *,
             npop: int,
             xpar: np.ndarray,
             simplex: np.ndarray,
             step: np.ndarray | None,
             factor: float | None
             ) -> np.ndarray:
        for ii in range(self.npar):
            tmp = np.copy(xpar)
            if 0.0 == tmp[ii]:
                tmp[ii] = 2.5e-4
            else:
                tmp[ii] *= 1.05
            simplex[ii + 1] = tmp

        npar1 = self.npar + 1
        return self.init_random_simplex(xpar, simplex, start=npar1,
                                        npop=npop, factor=factor)


class SimplexStep(SimplexBase):

    def init(self,
             *,
             npop: int,
             xpar: np.ndarray,
             simplex: np.ndarray,
             step: np.ndarray | None,
             factor: float | None
             ) -> np.ndarray:

        # This should raise an error, but as it is low-level code, just
        # use an assert.
        assert step is not None

        for ii in range(self.npar):
            tmp = xpar[ii] + step[ii]
            simplex[ii + 1] = tmp

        npar1 = self.npar + 1
        return self.init_random_simplex(xpar, simplex, start=npar1,
                                        npop=npop, factor=factor)


class SimplexRandom(SimplexBase):

    def init(self,
             *,
             npop: int,
             xpar: np.ndarray,
             simplex: np.ndarray,
             step: np.ndarray | None,
             factor: float | None
             ) -> np.ndarray:

        return self.init_random_simplex(xpar, simplex, start=1,
                                        npop=npop, factor=factor)
