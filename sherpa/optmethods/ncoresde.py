#!/usr/bin/env python

#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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
import random

from sherpa.optmethods.ncoresnm import ncoresNelderMead
from sherpa.optmethods.opt import Opt, Polytope, tst_opt
from sherpa.utils import parallel_map, _ncpus, Knuth_close, sao_fcmp


class Key2:

    def __init__(self, n=12):
        self.nbit = n
        self.max_arg2 = 2**n - 1
        return

    def calc(self, arg1, arg2):
        if arg2 > self.max_arg2:
            str = "arg2 ({}) must be < {}".format(arg2, self.max_arg2)
            raise ValueError(str)
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

    def __init__(self, func, npar, npop, sfactor, xprob):
        self.func = func
        self.npar  = npar
        self.npop = npop
        self.sfactor = sfactor
        self.xprob = xprob
        return

    def calc(self, arg, pop):
        arg[-1] = self.func(arg[:-1])
        tmp = numpy.empty(self.npar + 2)
        tmp[1:] = arg[:]
        if numpy.float_(numpy.finfo(numpy.float_).max) == arg[-1]:
            tmp[0] = 0
        else:
            tmp[0] = 1
        return tmp

    def init(self, num):
        return random.sample(range(self.npop), num)


class Strategy0(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for _ in range(self.npar):
            trial[n] = pop[0][n] + self.sfactor * (pop[r2][n] - pop[r3][n])
            n = (n + 1) % self.npar
            if random.uniform(0, 1) > self.xprob:
                break
        return self.calc(trial, pop)


class Strategy1(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for _ in range(self.npar):
            trial[n] = trial[n] + self.sfactor * (pop[r2][n] - pop[r3][n])
            n = (n + 1) % self.npar
            if random.uniform(0, 1) > self.xprob:
                break
        return self.calc(trial, pop)


class Strategy2(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2 = self.init(2)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for _ in range(self.npar):
            trial[n] = trial[n] + self.sfactor* (pop[0][n] - trial[n]) + \
                self.sfactor* (pop[r1][n] - pop[r2][n])
            n = (n + 1) % self.npar
            if random.uniform(0, 1) > self.xprob:
                break
        return self.calc(trial, pop)


class Strategy3(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4 = self.init(4)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for _ in range(self.npar):
            trial[n] = pop[0][n]  + \
                (pop[r1][n] + pop[r2][n] - pop[r3][n] - pop[r4][n]) * \
                self.sfactor
            n = (n + 1) % self.npar
            if random.uniform(0, 1) > self.xprob:
                break
        return self.calc(trial, pop)


class Strategy4(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4, r5 = self.init(5)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for _ in range(self.npar):
            trial[n] = pop[r5][n]  + \
                (pop[r1][n] + pop[r2][n] - pop[r3][n] - pop[r4][n]) * \
                self.sfactor
            n = (n + 1) % self.npar
            if random.uniform(0, 1) > self.xprob:
                break
        return self.calc(trial, pop)


class Strategy5(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for counter in range(self.npar):
            if random.uniform(0, 1) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[0][n]  + \
                    self.sfactor *(pop[r2][n] - pop[r3][n])
                n = (n + 1) % self.npar
        return self.calc(trial, pop)


class Strategy6(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3 = self.init(3)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for counter in range(self.npar):
            if random.uniform(0, 1) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[r1][n]  + self.sfactor * \
                    (pop[r2][n] - pop[r3][n])
                n = (n + 1) % self.npar
        return self.calc(trial, pop)


class Strategy7(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2 = self.init(2)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for counter in range(self.npar):
            if random.uniform(0, 1) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] += self.sfactor * ((pop[0][n] - trial[n]) + \
                                             (pop[r1][n] - pop[r2][n]))
                n = (n + 1) % self.npar
        return self.calc(trial, pop)


class Strategy8(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4 = self.init(4)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for counter in range(self.npar):
            if random.uniform(0, 1) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[0][n]  + \
                    self.sfactor * (pop[r2][n] - pop[r3][n] - pop[r4][n])
                n = (n + 1) % self.npar
        return self.calc(trial, pop)


class Strategy9(Strategy):

    def __init__(self, func, npar, npop, sfactor, xprob):
        Strategy.__init__(self, func, npar, npop, sfactor, xprob)
        return

    def __call__(self, pop, icurrent):
        r1, r2, r3, r4, r5 = self.init(5)
        trial = numpy.array(pop[icurrent][:])
        n = random.randint(0, self.npar - 1)
        for counter in range(self.npar):
            if random.uniform(0, 1) < self.xprob or \
                    counter == self.npar - 1:
                trial[n] = pop[r5][n] + \
                    self.sfactor * (pop[r1][n] +pop[r2][n] - pop[r3][n] - \
                                         pop[r4][n])
                n = (n + 1) % self.npar
        return self.calc(trial, pop)


class MyDifEvo(Opt):

    def __init__(self, func, xpar, xmin, xmax, npop, sfactor, xprob, step,
                 seed):
        Opt.__init__(self, func, xmin, xmax)
        self.ncores_nm = ncoresNelderMead()
        self.key2 = Key2()
        self.npop = min(npop, 4096)
        self.seed = seed
        self.strategies = \
            (Strategy0(self.func, self.npar, npop, sfactor, xprob),
             Strategy1(self.func, self.npar, npop, sfactor, xprob),
             Strategy2(self.func, self.npar, npop, sfactor, xprob),
             Strategy3(self.func, self.npar, npop, sfactor, xprob),
             Strategy4(self.func, self.npar, npop, sfactor, xprob),
             Strategy5(self.func, self.npar, npop, sfactor, xprob),
             Strategy6(self.func, self.npar, npop, sfactor, xprob),
             Strategy7(self.func, self.npar, npop, sfactor, xprob),
             Strategy8(self.func, self.npar, npop, sfactor, xprob),
             Strategy9(self.func, self.npar, npop, sfactor, xprob))

        xpar = numpy.asarray(xpar)
        if step is None:
            step = xpar * 1.2 + 1.2
        self.polytope = Polytope(func, npop, xpar, xmin, xmax, step, seed)
        self.local_opt = self.ncores_nm.algo
        return

    def __call__(self, maxnfev, ftol):
        
        random.seed(self.seed)
        mypop = self.polytope
        npop_1 = self.npop - 1
        while self.nfev[0] < maxnfev:
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
        return self.nfev[0], best_val, best_par

    def all_strategies(self, key):
        rand, index = self.key2.parse(key)
        random.seed(rand)
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
        if numpy.isinf(result[1]) is False:
            arg[0] += result[0]
            arg[1:-1] = result[2][:]
            arg[-1] = result[1]
        return arg

    def calc_key(self, indices, start=0, end=65536):
        result = numpy.empty(len(indices), dtype=numpy.int64)
        for ii, index in enumerate(indices):
            rand = random.randint(start, end)
            result[ii] = self.key2.calc(rand, index)
        return result
            
    def check_convergence(self, mypop, ftol, npar):
        fval_std = numpy.std([col[-1] for col in mypop])
        if fval_std < ftol:
            return True
        return False

    
class ncoresMyDifEvo(MyDifEvo):

    def __init__(self, func, xpar, xmin, xmax, npop, sfactor, xprob, step,
                 seed):
        MyDifEvo.__init__(self, func, xpar,xmin, xmax, npop, sfactor, xprob,
                          step, seed)
        return

    def __call__(self, tol, maxnfev, numcores=_ncpus):
        nfev = 0
        random.seed(self.seed)
        mypop = self.polytope
        old_fval = numpy.inf
        while nfev < maxnfev:

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

    def __init__(self):
        pass
    
    def __call__(self, fcn, x, xmin, xmax, step=None, maxnfev=None, tol=1.0e-6,
                 npop=None, seed=45, sfactor=0.85, xprob=0.7, verbose=0):

        npar = len(x)
        if npop is None:
            npop = 10 * npar
        npop = max(npop, npar * 4)
        if maxnfev is None:
            maxnfev = 8192 * npar

        mydifevo = MyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob, step)
        return mydifevo(maxnfev, tol, seed, step)


class ncoresDifEvo:

    def __init__(self):
        pass

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6, maxnfev=None, step=None,
                 numcores=None, npop=None, seed=23, sfactor=0.85, xprob=0.7,
                 verbose=0):

        npar = len(x)
        if npop is None:
            npop = 10 * npar
        npop = min(npop, npar * 24)
        if maxnfev is None:
            maxnfev = 8192 * npar

        mydifevo = ncoresMyDifEvo(fcn, x, xmin, xmax, npop, sfactor, xprob, \
                                  step, seed)
        return mydifevo(tol, maxnfev, numcores)


class ncoresDifEvoNelderMead:

    def __init__(self):
        self.ncores_nm = ncoresNelderMead()        
        return
    
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
        mydifevo = \
            ncoresMyDifEvo(fcn, nm_par, xmin, xmax, npop, sfactor, xprob, step)
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


if '__main__' == __name__:

    from sherpa.optmethods.opt import tst_opt, Rosenbrock
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--difevo', action="store_true",
                        default=False, help='run simple difevo', dest="difevo")
    parser.add_argument('-c', '--combine', action="store_true",
                        default=False, help='run nm & difevo', dest="combine")
    parser.add_argument('-N', action="store", dest="num", default=4, type=int)

    options = parser.parse_args()
    print('options =', options)

    npar = options.num
    if options.difevo:
        tst_opt([DifEvo()], npar)
    elif options.combine:
        tst_opt([ncoresDifEvoNelderMead()], npar)
    else:
        tst_opt([ncoresDifEvo()], npar)
