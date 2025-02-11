#
#  Copyright (C) 2019, 2020, 2021, 2023
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

import numpy

from sherpa.optmethods.ncoresnm import ncoresNelderMead
from sherpa.optmethods.opt import Opt, SimplexRandom
from sherpa.utils.parallel import parallel_map, ncpus
from sherpa.utils import random


class Key2:

    def __init__(self, n=12):
        self.nbit = n
        self.max_arg2 = 2**n - 1

    def calc(self, arg1, arg2):
        if arg2 > self.max_arg2:
            raise ValueError(f"arg2 ({arg2}) must be < {self.max_arg2}")

        key = arg1 + 1
        key <<= self.nbit
        key += arg2
        return key

    def parse(self, key):
        arg1 = key
        arg1 >>= self.nbit
        arg2 = arg1
        arg2 <<= self.nbit
        arg2 = key - arg2
        return arg1 - 1, arg2


class Strategy:

    def __init__(self, func, npar, npop, sfactor, xprob, rng=None):
        self.func = func
        self.npar = npar
        self.npop = npop
        self.sfactor = sfactor
        self.xprob = xprob
        self.rng = rng

    def calc(self, arg, pop):
        arg[-1] = self.func(arg[:-1])
        tmp = numpy.empty(self.npar + 2)
        tmp[1:] = arg[:]
        if numpy.finfo(numpy.float64).max == arg[-1]:
            tmp[0] = 0
        else:
            tmp[0] = 1
        return tmp

    def init(self, num):
        return random.choice(self.rng, range(self.npop), num)


class Strategy0(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = pop[0][n] + self.sfactor * (pop[r2][n] - pop[r3][n])
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy1(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = trial[n] + self.sfactor * (pop[r2][n] - pop[r3][n])
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy2(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2 = self.init(2)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = trial[n] + self.sfactor * (pop[0][n] - trial[n]) + \
                self.sfactor * (pop[r1][n] - pop[r2][n])
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy3(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4 = self.init(4)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = pop[0][n] + \
                (pop[r1][n] + pop[r2][n] - pop[r3][n] - pop[r4][n]) * \
                self.sfactor
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy4(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4, r5 = self.init(5)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for _ in range(self.npar):
            trial[n] = pop[r5][n] + \
                (pop[r1][n] + pop[r2][n] - pop[r3][n] - pop[r4][n]) * \
                self.sfactor
            n = (n + 1) % self.npar
            if random.random(self.rng) > self.xprob:
                break

        return self.calc(trial, pop)


class Strategy5(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[0][n] + \
                    self.sfactor * (pop[r2][n] - pop[r3][n])
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy6(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[r1][n] + self.sfactor * \
                    (pop[r2][n] - pop[r3][n])
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy7(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2 = self.init(2)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] += self.sfactor * ((pop[0][n] - trial[n]) +
                                            (pop[r1][n] - pop[r2][n]))
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy8(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4 = self.init(4)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[0][n] + \
                    self.sfactor * (pop[r2][n] - pop[r3][n] - pop[r4][n])
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class Strategy9(Strategy):

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4, r5 = self.init(5)
        trial = numpy.array(pop[icurrent][:])
        n = random.integers(self.rng, self.npar)
        for counter in range(self.npar):
            if random.random(self.rng) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[r5][n] + \
                    self.sfactor * (pop[r1][n] + pop[r2][n] - pop[r3][n] -
                                    pop[r4][n])
                n = (n + 1) % self.npar

        return self.calc(trial, pop)


class MyDifEvo(Opt):

    def __init__(self, func, xpar, xmin, xmax, npop, sfactor, xprob, step,
                 seed, rng=None):
        Opt.__init__(self, func, xmin, xmax)
        self.ncores_nm = ncoresNelderMead()
        self.key2 = Key2()
        self.npop = min(npop, 4096)
        self.seed = seed

        # Create separate RNGs for the strategy elements, following
        # https://numpy.org/doc/stable/reference/random/parallel.html,
        # to allow this code to run with parallel_map.
        #
        # This approach is used even when rng is None, and it always
        # uses default_rng to create the RNG.
        #
        strats = [Strategy0, Strategy1, Strategy2, Strategy3, Strategy4,
                  Strategy5, Strategy6, Strategy7, Strategy8, Strategy9]
        sseeds = numpy.random.SeedSequence(seed).spawn(len(strats))
        self.strategies = [strat(self.func, self.npar, npop, sfactor, xprob,
                                 rng=numpy.random.default_rng(sseed))
                           for strat, sseed in zip(strats, sseeds)]

        xpar = numpy.asarray(xpar)
        if step is None:
            step = xpar * 1.2 + 1.2
        factor = 10
        self.rng = rng
        self.polytope = SimplexRandom(func=func, npop=npop, xpar=xpar,
                                      xmin=xmin, xmax=xmax, step=step,
                                      seed=seed, factor=factor, rng=rng)
        self.local_opt = self.ncores_nm.algo

    def __call__(self, maxnfev, ftol):

        # Set the seed if RNG is not sent in. This used to change
        # random.seed but now changes the NumPy version.
        #
        if self.rng is None:
            numpy.random.seed(self.seed)

        mypop = self.polytope
        npop_1 = self.npop - 1
        while self.nfev < maxnfev:
            for pop_index in range(self.npop):
                key = self.calc_key([pop_index])
                # trial = self.all_strategies(pop_index)
                trial = self.all_strategies(key[0])
                if trial[-1] < mypop[npop_1][-1]:
                    mypop[npop_1] = trial[1:]
                    self.polytope.sort()

            if self.check_convergence(mypop, ftol, 0):
                break

        best_vertex = mypop[0]
        best_par = best_vertex[:-1]
        best_val = best_vertex[-1]
        return self.nfev, best_val, best_par

    def all_strategies(self, key):
        rand, index = self.key2.parse(key)

        # Set the seed if RNG is not sent in. This used to change
        # random.seed but now changes the NumPy version. It is not
        # clear if the int is needed but leave in for now (Python 3.11
        # needed it when it was random.seed).
        #
        if self.rng is None:
            numpy.random.seed(int(self.seed))

        mypop = self.polytope
        best_trial = self.strategies[0](mypop, index)
        for ii in range(1, len(self.strategies)):
            trial = self.strategies[ii](mypop, index)
            if trial[-1] < best_trial[-1]:
                best_trial = trial

        if best_trial[-1] < mypop[0][-1]:
            best_trial = self.apply_local_opt(best_trial, index)
        return best_trial

    def apply_local_opt(self, arg, index):
        local_opt = self.local_opt[index % len(self.local_opt)]
        result = local_opt(self.func, arg[1:-1], self.xmin, self.xmax)
        tmp = numpy.append(result[0], result[2])
        result = numpy.append(tmp, result[1])
        return result

    def calc_key(self, indices, start=0, end=65536):
        result = numpy.empty(len(indices), dtype=numpy.int64)
        for ii, index in enumerate(indices):
            # want to generate [start, end)
            rand = random.integers(self.rng, end - start) + start
            result[ii] = self.key2.calc(rand, index)
        return result

    def check_convergence(self, mypop, ftol, npar):
        fval_std = numpy.std([col[-1] for col in mypop])
        if fval_std < ftol:
            return True
        return False


class ncoresMyDifEvo(MyDifEvo):

    def __call__(self, tol, maxnfev, numcores=ncpus):
        nfev = 0

        # Set the seed if RNG is not sent in. This used to change
        # random.seed but now changes the NumPy version.
        #
        if self.rng is None:
            numpy.random.seed(self.seed)

        mypop = self.polytope
        old_fval = numpy.inf
        while nfev < maxnfev:

            # all_strategies has been set up so that each strategy has
            # its own generator, so it's okay to use parallel_map here,
            # and not have to think about parallel_map_rng.
            #
            keys = self.calc_key(range(self.npop))
            results = \
                parallel_map(self.all_strategies, keys, numcores)

            for index, result in enumerate(results):
                nfev += int(result[0])
                if result[-1] < mypop[index][-1]:
                    mypop[index] = result[1:]

            self.polytope.sort()
            if self.polytope.check_convergence(tol, 0):
                break

            best = mypop[0]
            best_fval = best[-1]
            if best_fval < old_fval:
                best_par = best[:-1]
                tmp_nfev, tmp_fval, tmp_par = \
                    self.ncores_nm(self.func, best_par, self.xmin, self.xmax,
                                   tol)
                nfev += tmp_nfev
                if tmp_fval < best_fval:
                    best_par = numpy.append(tmp_par, tmp_fval)
                    mypop[1] = best_par[:]
                    self.polytope.sort()
                    old_fval = tmp_fval
                else:
                    old_fval = best_fval

        best_vertex = self.polytope[0]
        best_par = best_vertex[:-1]
        best_fval = best_vertex[-1]
        return nfev, best_fval, best_par


class DifEvo:

    # The classes tend to take rng as an argument when constructing the
    # object, so follow that approach here.
    #
    def __init__(self, rng=None):
        self.rng = rng

    def __call__(self, fcn, x, xmin, xmax, step=None, maxnfev=None, tol=1.0e-6,
                 npop=None, seed=45, sfactor=0.85, xprob=0.7, verbose=0):

        npar = len(x)
        if npop is None:
            npop = 10 * npar
        npop = max(npop, npar * 4)
        if maxnfev is None:
            maxnfev = 8192 * npar

        # TODO: this over-writes the seed argument, which looks wrong
        seed = 123
        mydifevo = MyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob, step,
                            seed, rng=self.rng)
        return mydifevo(maxnfev, tol)


class ncoresDifEvo:

    # The classes tend to take rng as an argument when constructing the
    # object, so follow that approach here.
    #
    def __init__(self, rng=None):
        self.rng = rng

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6, maxnfev=None, step=None,
                 numcores=None, npop=None, seed=23, sfactor=0.85, xprob=0.7,
                 verbose=0):

        npar = len(x)
        if npop is None:
            npop = 10 * npar
        npop = min(npop, npar * 24)
        if maxnfev is None:
            maxnfev = 8192 * npar

        mydifevo = ncoresMyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob,
                                  step, seed, rng=self.rng)
        return mydifevo(tol, maxnfev, numcores)


class ncoresDifEvoNelderMead:

    def __init__(self, rng=None):
        self.ncores_nm = ncoresNelderMead()
        self.rng = rng

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6, maxnfev=None, step=None,
                 numcores=None, npop=None, seed=23, sfactor=0.85, xprob=0.7,
                 verbose=0):

        nfev, nm_fmin, nm_par = \
            self.ncores_nm(fcn, x, xmin, xmax, tol, maxnfev, numcores)

        npar = len(x)
        if npop is None:
            npop = 12 * npar
        npop = max(npop, npar * 32)
        if maxnfev is None:
            maxnfev = 8192 * npar

        # TODO: the seed argument is not sent in
        mydifevo = \
            ncoresMyDifEvo(fcn, nm_par, xmin, xmax, npop, sfactor, xprob,
                           step, rng=self.rng)
        de_nfev, de_fmin, de_par = \
            mydifevo(tol, maxnfev - nfev, step, seed, numcores)
        nfev += de_nfev

        if nm_fmin < de_fmin:
            my_fmin = nm_fmin
            my_par = nm_par
        else:
            my_fmin = de_fmin
            my_par = de_par
        nm_nfev, nm_fmin, nm_par = self.ncores_nm(fcn, my_par, xmin, xmax, tol,
                                                  maxnfev - nfev, numcores)
        nfev += nm_nfev

        if nm_fmin < my_fmin:
            my_fmin = nm_fmin
            my_par = nm_par

        return nfev, my_fmin, my_par
