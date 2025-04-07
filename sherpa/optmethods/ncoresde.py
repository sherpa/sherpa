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

from typing import Any

import numpy as np

from sherpa.utils.parallel import parallel_map, ncpus
from sherpa.utils import random

from .ncoresnm import ncoresNelderMead
from .opt import FUNC_MAX, Opt, OptimizerFunc, MyOptOutput, SimplexRandom


class Key2:

    def __init__(self, n: int = 12) -> None:
        self.nbit = n
        self.max_arg2 = 2**n - 1

    def calc(self, arg1: int, arg2: int) -> int:
        if arg2 > self.max_arg2:
            raise ValueError(f"arg2 ({arg2}) must be < {self.max_arg2}")

        key = arg1 + 1
        key <<= self.nbit
        key += arg2
        return key

    def parse(self, key: int) -> tuple[int, int]:
        arg1 = key
        arg1 >>= self.nbit
        arg2 = arg1
        arg2 <<= self.nbit
        arg2 = key - arg2
        return arg1 - 1, arg2


class Strategy:
    """Create a trial set of parameters.

    The RNG is sent when the strategy is created, rather than when
    called, as it may be called in parallel, and so needs a unique
    generator.

    """

    def __init__(self,
                 func: OptimizerFunc,
                 npar: int,
                 npop: int,
                 sfactor: float,
                 xprob: float,
                 rng: random.RandomType | None = None
                 ) -> None:
        self.func = func
        self.npar = npar
        self.npop = npop
        self.sfactor = sfactor
        self.xprob = xprob
        self.rng = rng

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        raise NotImplementedError

    # The arg argument contains just the parameter values.
    #
    def calc(self,
             # do we need to send in arg like this as it looks to
             # be more than just the pars?
             arg: np.ndarray,
             pop: Any  # unused
             ) -> MyOptOutput:
        funcval = self.func(arg)
        nfev = 1 if funcval != FUNC_MAX else 0
        return (nfev, funcval, arg)

    def init(self, num: int) -> np.ndarray:
        return random.choice(self.rng, range(self.npop), num)


class Strategy0(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        # Although only two numbers are needed, leave as is since the
        # code has been tested using this call.
        _, r2, r3 = self.init(3)
        trial = pop[icurrent].copy()
        n = random.integers(self.rng, self.npar)

        base = pop[0]
        # Pull out the terms which do not change (so not 'trial').
        delta = pop[r2] - pop[r3]

        for _ in range(self.npar):
            trial[n] = base[n] + self.sfactor * delta[n]
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy1(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        # Although only two numbers are needed, leave as is since the
        # code has been tested using this call.
        _, r2, r3 = self.init(3)
        trial = pop[icurrent].copy()

        base = trial
        delta = pop[r2] - pop[r3]

        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = base[n] + self.sfactor * delta[n]
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy2(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2 = self.init(2)
        trial = pop[icurrent].copy()

        base = trial
        delta = pop[0] + pop[r1] - pop[r2]

        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = base[n] + self.sfactor * (delta[n] - trial[n])
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy3(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2, r3, r4 = self.init(4)
        trial = pop[icurrent].copy()

        base = pop[0]
        delta = pop[r1] + pop[r2] - pop[r3] - pop[r4]

        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = base[n] + self.sfactor * delta[n]
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy4(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2, r3, r4, r5 = self.init(5)
        trial = pop[icurrent].copy()

        base = pop[r5]
        delta = pop[r1] + pop[r2] - pop[r3] - pop[r4]

        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = base[n] + self.sfactor * delta[n]
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy5(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        # Although only two numbers are needed, leave as is since the
        # code has been tested using this call.
        _, r2, r3 = self.init(3)
        trial = pop[icurrent].copy()

        base = pop[0]
        delta = pop[r2] - pop[r3]

        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = base[n] + self.sfactor * delta[n]
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy6(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2, r3 = self.init(3)
        trial = pop[icurrent].copy()

        base = pop[r1]
        delta = pop[r2] - pop[r3]

        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = base[n] + self.sfactor * delta[n]
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy7(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2 = self.init(2)
        trial = pop[icurrent].copy()

        base = trial
        delta = pop[0] + pop[r1] - pop[r2]

        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = base[n] + self.sfactor * (delta[n] - trial[n])
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy8(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        # Although only three numbers are needed, leave as is since the
        # code has been tested using this call.
        _, r2, r3, r4 = self.init(4)
        trial = pop[icurrent].copy()

        base = pop[0]
        delta = pop[r2] - pop[r3] - pop[r4]

        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = base[n] + self.sfactor * delta[n]
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy9(Strategy):

    def __call__(self,
                 pop: SimplexRandom,
                 icurrent: int
                 ) -> MyOptOutput:
        r1, r2, r3, r4, r5 = self.init(5)
        trial = pop[icurrent].copy()

        base = pop[r5]
        delta = pop[r1] + pop[r2] - pop[r3] - pop[r4]

        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = base[n] + self.sfactor * delta[n]
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class MyDifEvo(Opt):
    """

    .. versionchanged:: 4.17.1
       Some of the methods now use more-structured return types.

    """

    def __init__(self,
                 func: OptimizerFunc,
                 xpar: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 npop: int,
                 sfactor: float,
                 xprob: float,
                 step: np.ndarray | None,
                 seed: int,
                 rng: random.RandomType | None = None
                 ) -> None:
        super().__init__(func, xmin, xmax)

        self.ncores_nm = ncoresNelderMead()
        self.key2 = Key2()
        self.npop = min(npop, 4096)

        self.seed = seed
        if rng is None:
            # Create the RNG if not set. For now use the seed, but
            # perhaps it should be called with no argument.
            #
            self.rng = np.random.default_rng(seed)

        else:
            self.rng = rng

        # Create separate RNGs for the strategy elements, following
        # https://numpy.org/doc/stable/reference/random/parallel.html,
        # to allow this code to run with parallel_map.
        #
        # This approach is used even when rng is None, and it always
        # uses default_rng to create the RNG.
        #
        strats = [Strategy0, Strategy1, Strategy2, Strategy3, Strategy4,
                  Strategy5, Strategy6, Strategy7, Strategy8, Strategy9]

        # See also sherpa.utils.parallel.create_seeds.
        #
        # Should the seed for SeedSequence be created from the RNG
        # (so, use create_seeds) rather than using a hard-coded value?
        #
        sseeds = np.random.SeedSequence(seed).spawn(len(strats))
        self.strategies = [strat(self.func, self.npar, npop, sfactor, xprob,
                                 rng=np.random.default_rng(sseed))
                           for strat, sseed in zip(strats, sseeds)]

        xpar_np = np.asarray(xpar)
        xmin_np = np.asarray(xmin)
        xmax_np = np.asarray(xmax)
        if step is None:
            step = xpar_np * 1.2 + 1.2

        factor = 10

        # TODO: should this create a separate rng from self.rng?
        self.polytope = SimplexRandom(func=func, npop=npop,
                                      xpar=xpar_np, xmin=xmin_np,
                                      xmax=xmax_np, step=step,
                                      seed=seed, factor=factor,
                                      rng=self.rng)
        self.local_opt = self.ncores_nm.algo

    # Only used by DifEvo which is currently unused.
    #
    # def __call__(self, maxnfev, ftol):
    #
    #     # Set the seed if RNG is not sent in. This used to change
    #     # random.seed but now changes the NumPy version.
    #     #
    #     if self.rng is None:
    #         numpy.random.seed(self.seed)
    #
    #     mypop = self.polytope
    #     npop_1 = self.npop - 1
    #     while self.nfev < maxnfev:
    #         for pop_index in range(self.npop):
    #             key = self.calc_key([pop_index])
    #             # trial = self.all_strategies(pop_index)
    #             trial = self.all_strategies(key[0])
    #             if trial[-1] < mypop[npop_1][-1]:
    #                 mypop[npop_1] = trial[1:]
    #                 self.polytope.sort()
    #
    #         if self.check_convergence(mypop, ftol, 0):
    #             break
    #
    #     best_vertex = mypop[0]
    #     best_par = best_vertex[:-1]
    #     best_val = best_vertex[-1]
    #     return self.nfev, best_val, best_par

    def all_strategies(self, key: int) -> MyOptOutput:
        _, index = self.key2.parse(key)

        mypop = self.polytope
        best_trial = self.strategies[0](mypop, index)
        for ii in range(1, len(self.strategies)):
            trial = self.strategies[ii](mypop, index)
            if trial[1] < best_trial[1]:
                best_trial = trial

        if best_trial[1] < mypop.fctvals[0]:
            return self.apply_local_opt(best_trial, index)

        return best_trial

    def apply_local_opt(self,
                        arg: MyOptOutput,
                        index: int
                        ) -> MyOptOutput:
        local_opt = self.local_opt[index % len(self.local_opt)]
        return local_opt(self.func, arg[2], self.xmin, self.xmax,
                         rng=self.rng)

    def calc_key(self,
                 indices,
                 start: int = 0,
                 end: int = 65536
                 ) -> np.ndarray:
        result = np.empty(len(indices), dtype=np.int64)
        for ii, index in enumerate(indices):
            # want to generate [start, end)
            rand = random.integers(self.rng, end - start) + start
            result[ii] = self.key2.calc(rand, index)
        return result

    # Only used by DifEvo which is currently unused.
    #
    # def check_convergence(self, mypop, ftol, npar):
    #     fval_std = numpy.std([col[-1] for col in mypop])
    #     if fval_std < ftol:
    #         return True
    #     return False


class ncoresMyDifEvo(MyDifEvo):
    """

    .. versionchanged:: 4.17.1
       The calling convention has been changed to match its superclass.

    """

    def __call__(self,
                 maxnfev: int,
                 ftol: float,
                 numcores: int | None = ncpus
                 ) -> MyOptOutput:

        nfev = 0

        # Set the seed if RNG is not sent in. This used to change
        # random.seed but now changes the NumPy version.
        #
        if self.rng is None:
            np.random.seed(self.seed)

        mypop = self.polytope
        old_fval = np.inf
        while nfev < maxnfev:

            # all_strategies has been set up so that each strategy has
            # its own generator, so it's okay to use parallel_map here,
            # and not have to think about parallel_map_rng.
            #
            keys = self.calc_key(range(self.npop))
            results: list[MyOptOutput] = \
                parallel_map(self.all_strategies, keys, numcores)

            for index, (nfev_r, fval_r, pars_r) in enumerate(results):
                nfev += int(nfev_r)
                if fval_r < mypop.fctvals[index]:
                    mypop[index] = pars_r
                    mypop.fctvals[index] = fval_r

            self.polytope.sort()
            if self.polytope.check_convergence(ftol, 0):
                break

            best_fval = mypop.fctvals[0]
            if best_fval < old_fval:
                best_par = mypop[0]
                tmp_nfev, tmp_fval, tmp_par = \
                    self.ncores_nm(self.func, best_par, self.xmin,
                                   self.xmax, tol=ftol)
                nfev += tmp_nfev
                if tmp_fval < best_fval:
                    mypop[1] = tmp_par
                    mypop.fctvals[1] = tmp_fval
                    self.polytope.sort()
                    old_fval = tmp_fval
                else:
                    old_fval = best_fval

        return nfev, self.polytope.fctvals[0], self.polytope[0]


# This is only used by tests/test_opt_original.py when run directly,
# not via pytest.
#
#  class DifEvo:
#
#      # The classes tend to take rng as an argument when constructing the
#      # object, so follow that approach here.
#      #
#      def __init__(self, rng=None):
#          self.rng = rng
#
#      def __call__(self, fcn, x, xmin, xmax, step=None, maxnfev=None, tol=1.0e-6,
#                   npop=None, seed=45, sfactor=0.85, xprob=0.7, verbose=0):
#
#          npar = len(x)
#          if npop is None:
#              npop = 10 * npar
#          npop = max(npop, npar * 4)
#          if maxnfev is None:
#              maxnfev = 8192 * npar
#
#          # TODO: this over-writes the seed argument, which looks wrong
#          seed = 123
#          mydifevo = MyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob, step,
#                              seed, rng=self.rng)
#          return mydifevo(maxnfev, tol)


class ncoresDifEvo:
    """

    .. versionchanged:: 4.17.1
       The rng argument is now set when calling the class, not when
       creating it.

    """

    def __call__(self,
                 fcn: OptimizerFunc,
                 x: np.ndarray,
                 xmin: np.ndarray,
                 xmax: np.ndarray,
                 tol: float = 1.0e-6,
                 maxnfev: int | None = None,
                 step: np.ndarray | None = None,
                 numcores: int | None = None,
                 npop: int | None = None,
                 seed: int = 23,
                 sfactor: float = 0.85,
                 xprob: float = 0.7,
                 verbose: Any = 0,  # unused
                 rng: random.RandomType | None = None
                 ) -> MyOptOutput:

        npar = len(x)
        if npop is None:
            npop = 10 * npar
        npop = min(npop, npar * 24)
        if maxnfev is None:
            maxnfev = 8192 * npar

        mydifevo = ncoresMyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob,
                                  step, seed, rng=rng)
        return mydifevo(ftol=tol, maxnfev=maxnfev, numcores=numcores)


# This is only used by tests/test_opt_original.py when run directly,
# not via pytest. The code fails so has been commented out.
#
# class ncoresDifEvoNelderMead:
#
#     def __init__(self, rng=None):
#         self.ncores_nm = ncoresNelderMead()
#         self.rng = rng
#
#     def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6, maxnfev=None, step=None# ,
#                  numcores=None, npop=None, seed=23, sfactor=0.85, xprob=0.7,
#                  verbose=0):
#
#         nfev, nm_fmin, nm_par = \
#             self.ncores_nm(fcn, x, xmin, xmax, tol, maxnfev, numcores)
#
#         npar = len(x)
#         if npop is None:
#             npop = 12 * npar
#         npop = max(npop, npar * 32)
#         if maxnfev is None:
#             maxnfev = 8192 * npar
#
#         # TODO: the seed argument is not sent in
#         mydifevo = \
#             ncoresMyDifEvo(fcn, nm_par, xmin, xmax, npop, sfactor, xprob,
#                            step, rng=self.rng)
#         de_nfev, de_fmin, de_par = \
#             mydifevo(tol, maxnfev - nfev, step, seed, numcores)
#         nfev += de_nfev
#
#         if nm_fmin < de_fmin:
#             my_fmin = nm_fmin
#             my_par = nm_par
#         else:
#             my_fmin = de_fmin
#             my_par = de_par
#         nm_nfev, nm_fmin, nm_par = self.ncores_nm(fcn, my_par, xmin, xmax, tol# ,
#                                                   maxnfev - nfev, numcores)
#         nfev += nm_nfev
#
#         if nm_fmin < my_fmin:
#             my_fmin = nm_fmin
#             my_par = nm_par
#
#         return nfev, my_fmin, my_par
