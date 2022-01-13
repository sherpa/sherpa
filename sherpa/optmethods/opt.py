#
#  Copyright (C) 2019, 2020, 2021, 2022
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

import multiprocessing
import random

import numpy as np

from sherpa.utils import Knuth_close, _multi, _ncpus, run_tasks, func_counter


__all__ = ('Opt', 'MyNcores', 'SimplexRandom', 'SimplexNoStep',
           'SimplexStep')


# import autograd.numpy as np
# from autograd import hessian
# from numpy.linalg.linalg import LinAlgError
# def autograd_covar(func, par):
#     try:
#         hess = hessian(func)(par)
#         result = 2.0 * np.linalg.inv(hess)
#         return result
#     except LinAlgError:
#         return None

class MyNcores:

    def __init__(self):
        if _multi is False:
            raise TypeError("multicores not available")

    def calc(self, funcs, numcores, *args, **kwargs):

        for func in funcs:
            if not callable(func):
                raise TypeError("input func '%s' is not callable" % repr(func))

        if numcores is None:
            numcores = _ncpus
        num_funcs = len(funcs)
        numcores = min(numcores, num_funcs)

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
        procs = []

        for id, func in enumerate(funcs):
            myargs = (func, id, out_q, err_q, lock) + args
            try:
                procs.append(multiprocessing.Process(target=self.my_worker,
                                                     args=myargs))
            except NotImplementedError as nie:
                raise nie
        return run_tasks(procs, err_q, out_q, num_funcs)

    def my_worker(self, *args):
        raise NotImplementedError("my_worker has not been implemented")


class Opt:

    def __init__(self, func, xmin, xmax):
        self.npar = len(xmin)
        self.nfev, self.func = \
            self.func_counter_bounds_wrappers(func, self.npar, xmin, xmax)
        self.xmin = np.asarray(xmin)
        self.xmax = np.asarray(xmax)

    def _outside_limits(self, x, xmin, xmax):
        return (np.any(x < xmin) or np.any(x > xmax))

    def func_bounds(self, func, npar, xmin=None, xmax=None):
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

        def func_bounds_wrapper(x, *args):
            if self._outside_limits(x, xmin, xmax):
                return np.finfo(np.float_).max
            return func(x, *args)

        return func_bounds_wrapper

    def func_counter_bounds_wrappers(self, func, npar, xmin, xmax):
        """Wraps the calls to func_counter then func_bounds"""
        nfev, afunc = func_counter(func)
        myfunc = self.func_bounds(afunc, npar, xmin, xmax)
        return nfev, myfunc


class SimplexBase:

    def __init__(self, func, npop, xpar, xmin, xmax, step, seed, factor):
        self.func = func
        self.xmin = xmin
        self.xmax = xmax
        self.npar = len(xpar)
        self.simplex = self.init(npop, xpar, step, seed, factor)

    def __getitem__(self, index):
        return self.simplex[index]

    def __setitem__(self, index, val):
        self.simplex[index] = val

    def calc_centroid(self):
        return np.mean(self.simplex[:-1, :], 0)

    def check_convergence(self, ftol, method):

        def are_func_vals_close_enough(tolerance):
            smallest_fct_val = self.simplex[0, -1]
            largest_fct_val = self.simplex[-1, -1]
            return Knuth_close(smallest_fct_val, largest_fct_val,
                               tolerance)

        def is_fct_stddev_small_enough(tolerance):
            fval_std = np.std([col[-1] for col in self.simplex])
            if fval_std < ftol:
                return True
            return False

        def is_max_length_small_enough(tolerance):
            """

               max  || x  - x  || <= tol max(1.0, || x ||)
                        i    0                         0

            where 1 <= i <= n"""

            max_xi_x0 = -1.0
            x0 = self.simplex[0, :-1]

            npar_plus_1 = len(x0) + 1
            for ii in range(1, npar_plus_1):
                xi = self.simplex[ii, :-1]
                xi_x0 = xi - x0
                max_xi_x0 = max(max_xi_x0, np.dot(xi_x0, xi_x0))
            return max_xi_x0 <= tolerance * max(1.0, np.dot(x0, x0))

        if 0 == method:
            if is_max_length_small_enough(ftol):
                return True
        elif 2 == method:
            if not is_max_length_small_enough(ftol):
                return False
            stddev = is_fct_stddev_small_enough(ftol)
            fctval = are_func_vals_close_enough(ftol)
            return stddev and fctval
        else:
            if not is_max_length_small_enough(ftol):
                return False
            stddev = is_fct_stddev_small_enough(ftol)
            fctval = are_func_vals_close_enough(ftol)
            return stddev or fctval
        return False

        # NOTE: this code is never used
        num = 2.0 * abs(self.simplex[0, -1] - self.simplex[-1, -1])
        denom = abs(self.simplex[0, -1]) + abs(self.simplex[-1, -1]) + 1.0
        if (num / denom > ftol):
            return False
        func_vals = [col[-1] for col in self.simplex]
        if np.std(func_vals) > ftol:
            return False
        return True

    def eval_simplex(self, npop, simplex):
        for ii in range(npop):
            simplex[ii][-1] = self.func(simplex[ii][:-1])
        return self.sort_me(simplex)

    def init(self, npop, xpar, step, seed, factor):
        raise NotImplementedError("init has not been implemented")

    def init_random_simplex(self, xpar, simplex, start, npop, seed, factor):
        random.seed(seed)
        for ii in range(start, npop):
            simplex[ii][:-1] = \
                np.array([random.uniform(max(self.xmin[jj],
                                             xpar[jj]-factor*abs(xpar[jj])),
                                         min(self.xmax[jj],
                                             xpar[jj]+factor*abs(xpar[jj])))
                          for jj in range(self.npar)])
        return simplex

    def move_vertex(self, centroid, coef):
        vertex = (1.0 + coef) * centroid - coef * self.simplex[self.npar]
        vertex[-1] = self.func(vertex[:-1])
        return vertex

    def shrink(self, shrink_coef):
        npars_plus_1 = self.npar + 1
        for ii in range(1, npars_plus_1):
            self.simplex[ii] = \
                self.simplex[0] + shrink_coef * \
                (self.simplex[ii] - self.simplex[0])
            self.simplex[ii, -1] = self.func(self.simplex[ii, :-1])

    def sort_me(self, simp):
        myshape = simp.shape
        tmp = np.array(sorted(simp, key=lambda arg: arg[-1]))
        tmp.reshape(myshape)
        return tmp

    def sort(self):
        self.simplex = self.sort_me(self.simplex)


class SimplexNoStep(SimplexBase):

    def init(self, npop, xpar, step, seed, factor):
        npar1 = self.npar + 1
        simplex = np.empty((npop, npar1))
        simplex[0][:-1] = np.copy(xpar)
        for ii in range(self.npar):
            tmp = np.copy(xpar)
            if 0.0 == tmp[ii]:
                tmp[ii] = 2.5e-4
            else:
                tmp[ii] *= 1.05
            simplex[ii+1][:-1] = tmp[:]
        simplex = \
            self.init_random_simplex(xpar, simplex, npar1, npop, seed, factor)
        return self.eval_simplex(npop, simplex)


class SimplexStep(SimplexBase):

    def init(self, npop, xpar, step, seed, factor):
        npar1 = self.npar + 1
        simplex = np.empty((npop, npar1))
        simplex[0][:-1] = np.copy(xpar)
        for ii in range(self.npar):
            tmp = xpar[ii] + step[ii]
            simplex[ii + 1][:-1] = tmp
        simplex = \
            self.init_random_simplex(xpar, simplex, npar1, npop, seed, factor)
        return self.eval_simplex(npop, simplex)


class SimplexRandom(SimplexBase):

    def init(self, npop, xpar, step, seed, factor):
        npar1 = self.npar + 1
        simplex = np.empty((npop, npar1))
        simplex[0][:-1] = np.copy(xpar)
        simplex = self.init_random_simplex(xpar, simplex, 1, npop, seed,
                                           factor)
        return self.eval_simplex(npop, simplex)
