#
#  Copyright (C) 2019-2021, 2023-2025
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

from abc import ABCMeta, abstractmethod
from collections.abc import Sequence

import numpy as np

from sherpa.utils.parallel import SupportsQueue, ncpus

from . import _saoopt  # type: ignore
from .opt import MyNcores, Opt, OptimizerFunc, WorkerFunc, OptOutput, \
    SimplexBase, SimplexStep

__all__ = ('ncoresNelderMead', )


# In the following it is unclear when the following hold
#
#   finalsimplex: int
#   finalsimplex: int | None
#   finalsimplex: Sequence[int] | int | None
#
# (and other variants). The current typing ignores the list option
# even though some code ends up setting a list.
#


EPSILON = float(np.finfo(np.float32).eps)


class MyNelderMead(Opt):
    """

    .. versionchanged:: 4.18.0
       The calling convention has changed to match its superclass.

    """

    def __init__(self,
                 fcn: OptimizerFunc,
                 xmin: np.ndarray,
                 xmax: np.ndarray
                 ) -> None:
        super().__init__(fcn, xmin, xmax)

        self.expansion_coef = 2.0          # chi
        self.contraction_coef = 0.5          # gamma
        self.reflection_coef = 1.0          # rho
        self.shrink_coef = 0.5          # sigma

    def __call__(self,
                 maxnfev: int,  # TODO: can this be None?
                 ftol: float,
                 xpar: np.ndarray,
                 step: np.ndarray,
                 finalsimplex: int,
                 verbose: int
                 ) -> OptOutput:

        # SimplexStep does not use factor when npop=npar + 1.
        #
        npar = len(xpar)
        simplex = SimplexStep(func=self.func, npop=npar + 1,
                              xpar=xpar, xmin=self.xmin,
                              xmax=self.xmax, step=step, seed=None,
                              factor=None)
        return self.optimize(xpar, simplex, maxnfev, ftol,
                             finalsimplex, verbose)

    def contract_in_out(self,
                        simplex: SimplexBase,
                        centroid: np.ndarray,
                        reflection_pt: np.ndarray,
                        rho_gamma: float,
                        contraction_coef: float,
                        badindex: int,
                        maxnfev: int,  # this is not used
                        verbose: int
                        ) -> bool:

        if simplex[badindex - 1, -1] <= reflection_pt[-1] and \
               reflection_pt[-1] < simplex[badindex, -1]:

            outside_contraction_pt = \
                simplex.move_vertex(centroid, rho_gamma)

            if outside_contraction_pt[-1] <= reflection_pt[-1]:
                simplex[badindex] = outside_contraction_pt
                if verbose > 2:
                    print('\taccept outside contraction point')
                    # and terminate the iteration
                return False

            # otherwise, go to step 5 (perform a shrink).
            return True

        if reflection_pt[-1] >= simplex[badindex, -1]:

            inside_contraction_pt = simplex.move_vertex(centroid,
                                                        - contraction_coef)

            if inside_contraction_pt[-1] < simplex[badindex, -1]:
                simplex[badindex] = inside_contraction_pt
                if verbose > 2:
                    print('\taccept inside contraction point')
                return False

            # otherwise, go to step 5 (perform a shrink).
            return True

        print('something is wrong with contract_in_out')
        return True

    def optimize(self,
                 xpar: np.ndarray,
                 simplex: SimplexBase,
                 maxnfev: int,
                 tol: float,
                 finalsimplex: int,
                 verbose: int
                 ) -> OptOutput:

        rho_chi = self.reflection_coef * self.expansion_coef
        rho_gamma = self.reflection_coef * self.contraction_coef

        # This is presumably because SimplexStep was called with
        # npop=npar + 1, so bad_index corresponds to the last pop
        # element.
        #
        bad_index = len(xpar)

        while self.nfev < maxnfev:

            simplex.sort()
            if verbose > 2:
                print(f'f{simplex[0, :-1]}={simplex[0, -1]:e}')

            if simplex.check_convergence(tol, finalsimplex):
                break

            centroid = simplex.calc_centroid()
            reflection_pt = simplex.move_vertex(centroid, self.reflection_coef)

            if simplex[0, -1] <= reflection_pt[-1] and \
                    reflection_pt[-1] < simplex[bad_index - 1, -1]:
                simplex[bad_index] = reflection_pt
                if verbose > 2:
                    print('\taccept reflection point')

            elif reflection_pt[-1] < simplex[0, -1]:
                expansion_pt = simplex.move_vertex(centroid, rho_chi)

                if expansion_pt[-1] < reflection_pt[-1]:
                    simplex[bad_index] = expansion_pt
                    if verbose > 2:
                        print('\taccept expansion point')
                else:
                    simplex[bad_index] = reflection_pt
                    if verbose > 2:
                        print('\taccept reflection point')

            elif self.contract_in_out(simplex, centroid, reflection_pt,
                                      rho_gamma, self.contraction_coef,
                                      bad_index, maxnfev, verbose):
                simplex.shrink(self.shrink_coef)

        best_vertex = simplex[0]
        best_par = best_vertex[:-1]
        best_val = best_vertex[-1]
        return OptOutput(nfev=self.nfev,
                         statval=best_val,
                         pars=best_par)


class NelderMeadBase(metaclass=ABCMeta):

    def __init__(self) -> None:
        self.nfev = 0
        self.fmin = np.inf
        self.par = np.nan

        # TODO: this should either be configurable or something that is
        # only changed when needed. See issue #1238
        #
        np.seterr(over='ignore', divide='ignore', under='ignore',
                  invalid='ignore')

    @abstractmethod
    def __call__(self,
                 fcn: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = 1.0e-6,
                 maxnfev: int | None = None,
                 step: np.ndarray | None = None,
                 finalsimplex: int | None = 1,
                 verbose: int = 0
                 ) -> OptOutput:
        raise NotImplementedError()

    def get_maxnfev(self, maxnfev: int | None, npar: int) -> int:
        if maxnfev is None:
            return 512 * npar

        return maxnfev


class NelderMead0(NelderMeadBase):

    def __call__(self,
                 fcn: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = 1.0e-6,
                 maxnfev: int | None = None,  # TODO: can we drop the None?
                 step: np.ndarray | None = None,
                 finalsimplex: int | None = 1,
                 verbose: int = 0
                 ) -> OptOutput:
        return self.neldermead0(fcn, xpar, xmin, xmax, step=step,
                                finalsimplex=finalsimplex,
                                maxnfev=maxnfev, tol=tol,
                                verbose=verbose)

    def calc_step(self, x: np.ndarray) -> np.ndarray:
        return 1.2 * x

    def neldermead0(self,
                    fcn: OptimizerFunc,
                    xpar: np.ndarray,
                    xmin: np.ndarray,
                    xmax: np.ndarray,
                    *,
                    step: np.ndarray | None = None,
                    finalsimplex: int | None = 1,
                    maxnfev: int | None = None,
                    tol: float = 1.0e-6,
                    verbose: int = 0
                    ) -> OptOutput:
        """

        .. versionchanged:: 4.18.0
           Most of the arguments must now be given by name.

        """
        x0 = np.asarray(xpar)
        maxnfev = self.get_maxnfev(maxnfev, len(x0))

        my_nm = MyNelderMead(fcn, xmin, xmax)
        if step is None:
            step = self.calc_step(x0)

        return my_nm(xpar=x0, maxnfev=maxnfev, ftol=tol, step=step,
                     finalsimplex=finalsimplex, verbose=verbose)


class NelderMead1(NelderMead0):

    def calc_step(self, x: np.ndarray) -> np.ndarray:
        return x + 1.2


class NelderMead2(NelderMead0):

    def calc_step(self, x: np.ndarray) -> np.ndarray:
        return abs(x)


class NelderMead3(NelderMead0):

    def __call__(self,
                 fcn: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = EPSILON,
                 maxnfev: int | None = None,
                 step: np.ndarray | None = None,
                 finalsimplex: int | None = None,
                 verbose: int = 0
                 ) -> OptOutput:

        # Avoid having a mutable argument
        if finalsimplex is None:
            finalsimplex = [0, 1, 1]

        x0 = np.asarray(xpar)
        n = len(x0)
        if step is None:
            step = np.full(n, 1.2)

        maxnfev = self.get_maxnfev(maxnfev, n)
        init = 0
        par, fmin, nfev, _ = \
            _saoopt.neldermead(verbose, maxnfev, init, finalsimplex,
                               tol, step, xmin, xmax, x0, fcn)

        return OptOutput(nfev=nfev,
                         statval=fmin,
                         pars=par)


class NelderMead4(NelderMead0):

    def __call__(self,
                 fcn: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = EPSILON,
                 maxnfev: int | None = None,
                 step: np.ndarray | None = None,
                 finalsimplex: int | None = None,
                 verbose: int = 0,
                 reflect: bool = True
                 ) -> OptOutput:

        # Avoid having a mutable argument
        if finalsimplex is None:
            finalsimplex = [0, 1, 1]

        x0 = np.asarray(xpar)
        n = len(x0)
        if step is None:
            step = abs(x0) + 1.2

        maxnfev = self.get_maxnfev(maxnfev, n)
        init = 0
        x0, _, nfev, _ = \
            _saoopt.neldermead(verbose, maxnfev, init, finalsimplex, tol, step,
                               xmin, xmax, x0, fcn)
        iquad = 1
        simp = 1.0e-2 * tol
        step = np.full(n, 0.4)
        par, fmin, tmpnfev, _ = \
            _saoopt.minim(reflect, verbose, maxnfev - nfev, init, iquad, simp,
                          tol*10, step, xmin, xmax, x0, fcn)
        nfev += tmpnfev
        return OptOutput(nfev=nfev,
                         statval=fmin,
                         pars=par)


class NelderMead5(NelderMead0):

    def __call__(self,
                 fcn: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = 1.0e-6,
                 maxnfev: int | None = None,
                 step: np.ndarray | None = None,
                 finalsimplex: int | None = 1,
                 verbose: int = 0,
                 reflect: bool = True
                 ) -> OptOutput:
        init = 0
        iquad = 1
        simp = 1.0e-2 * tol
        x0 = np.asarray(xpar)
        n = len(x0)
        if step is None:
            step = np.full(n, 0.4)

        maxnfev = self.get_maxnfev(maxnfev, n)
        par, fmin, nfev, _ = \
            _saoopt.minim(reflect, verbose, maxnfev, init, iquad, simp, tol*10,
                          step, xmin, xmax, x0, fcn)
        return OptOutput(nfev=nfev,
                         statval=fmin,
                         pars=par)


# This is only used by tests/test_opt_original.py when run directly,
# not via pytest. Left commented out to make it easier to find.
#
# class NelderMead6(NelderMeadBase):
#
#     # TODO: do we really need this internal class?
#     class MyNelderMead6(MyNelderMead):
#
#         def __call__(self, x, maxnfev, tol, step, finalsimplex, verbose):
#             npar = len(x)
#             simplex = SimplexNoStep(func=self.func, npop=npar + 1,
#                                     xpar=x, xmin=self.xmin,
#                                     xmax=self.xmax, step=None,
#                                     seed=None, factor=None)
#             return self.optimize(x, simplex, maxnfev, tol,
#                                  finalsimplex, verbose)
#
#     def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6,  maxnfev=None,
#                  step=None, finalsimplex=1, verbose=0):
#         my_nm_6 = NelderMead6.MyNelderMead6(fcn, xmin, xmax)
#         if maxnfev is None:
#             maxnfev = 512 * len(x)
#         return my_nm_6(x, maxnfev, tol, step, finalsimplex, verbose)


# This is only used by tests/test_opt_original.py when run directly,
# not via pytest. Left commented out to make it easier to find.
#
# class NelderMead7(NelderMeadBase):
#
#     # TODO: do we really need this internal class?
#     class MyNelderMead7(MyNelderMead):
#
#         def __call__(self, x, maxnfev, tol, step, finalsimplex, verbose):
#             npar = len(x)
#             factor = 2
#             simplex = SimplexRandom(func=self.func, npop=npar + 1, xpar=x,
#                                     xmin=self.xmin, xmax=self.xmax,
#                                     step=None, seed=None, factor=factor)
#             return self.optimize(x, simplex, maxnfev, tol,
#                                  finalsimplex, verbose)
#
#     def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6,  maxnfev=None,
#                  step=None, finalsimplex=1, verbose=0):
#         my_nm_7 = NelderMead7.MyNelderMead7(fcn, xmin, xmax)
#         if maxnfev is None:
#             maxnfev = 512 * len(x)
#         return my_nm_7(x, maxnfev, tol, step, finalsimplex, verbose)


class nmNcores(MyNcores):

    def my_worker(self,
                  opt: WorkerFunc,
                  idval: int,
                  out_q: SupportsQueue[tuple[int, list[OptOutput]]],
                  err_q: SupportsQueue[Exception],
                  fcn: OptimizerFunc,
                  x: np.ndarray,
                  xmin: np.ndarray,
                  xmax: np.ndarray,
                  tol: float,
                  maxnfev: int | None
                  ) -> None:
        try:
            vals = opt(fcn, x, xmin, xmax, tol, maxnfev)
        except Exception as e:
            err_q.put(e)
            return

        # The output queue is sent the positon and a list of
        # results. In this case only one result has been generated so
        # convert it to a list.
        #
        out_q.put((idval, [vals]))


class ncoresNelderMead:
    """

    .. versionchanged:: 4.18.0
       The default for the algo parameter is now None.

    """

    def __init__(self,
                 algo: Sequence[NelderMeadBase] | None = None
                 ) -> None:
        if algo is None:
            self.algo = [NelderMead0(), NelderMead1(), NelderMead2(),
                         NelderMead3(), NelderMead4(), NelderMead5()]
            # NelderMead6(), NelderMead7()]):

        else:
            self.algo = algo

    def __call__(self,
                 fcn: OptimizerFunc,
                 x: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = EPSILON,
                 maxnfev: int | None = None,
                 numcores=ncpus
                 ) -> OptOutput:

        nm_ncores = nmNcores()
        results = nm_ncores.calc(self.algo, numcores, fcn, x, xmin,
                                 xmax, tol, maxnfev)
        return self.unpack_results(results)

    def unpack_results(self,
                       results: list[OptOutput]
                       ) -> OptOutput:
        """

        .. versionchanged:: 4.18.0
           The num argument has been removed.

        """

        # As self.algo is not empty then results is not empty, but
        # this is hard to assert with types.
        #
        assert len(results) > 0

        nfev = results[0].nfev
        fmin = results[0].statval
        par = results[0].pars
        for res in results[1:]:
            nfev += res.nfev
            if res.statval < fmin:
                fmin = res.statval
                par = res.pars

        return OptOutput(nfev=nfev,
                         statval=fmin,
                         pars=par)


# This code is currently unused. Left commented out to make it easier
# to find.
#
# class ncoresNelderMeadRecursive(ncoresNelderMead):
#     """ As noted in the paper, terminating the simplex is not a simple task:
#     For any non-derivative method, the issue of termination is problematical as
#     well as highly sensitive to problem scaling. Since gradient information is
#     unavailable, it is provably impossible to verify closeness to optimality
#     simply by sampling f at a finite number of points. Most implementations
#     of direct search methods terminate based on two criteria intended to
#     reflect the progress of the algorithm: either the function values at the
#     vertices are close, or the simplex has become very small. """
#
#     def __call__(self,
#                  fcn: OptimizerFunc,
#                  x: np.ndarray,
#                  xmin: np.ndarray,
#                  xmax: np.ndarray,
#                  tol: float = EPSILON,
#                  maxnfev: int | None = None,
#                  numcores=ncpus
#                  ) -> OptOutput:
#
#         return self.calc(fcn, x, xmin, xmax, tol, maxnfev, numcores)
#
#     def calc(self,
#              fcn: OptimizerFunc,
#              x: np.ndarray,
#              xmin: np.ndarray,
#              xmax: np.ndarray,
#              tol: float = EPSILON,
#              maxnfev: int | None = None,
#              numcores=ncpus,
#              fval: float = np.inf,
#              nfev: int = 0
#              ) -> OptOutput:
#
#         nm_ncores = nmNcores()
#         results = nm_ncores.calc(self.algo, numcores, fcn, x, xmin,
#                                  xmax, tol, maxnfev)
#         tmp_nfev, fmin, par = self.unpack_results(results)
#         nfev += tmp_nfev
#         if fmin < fval:
#             return self.calc(fcn, par, xmin, xmax, tol, maxnfev,
#                              numcores, fval=fmin, nfev=nfev)
#
#         # TODO: shouldn't this return fmin rather than fval?
#         return nfev, fval, par
