#!/usr/bin/env python

#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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
           'SimplexStep', 'tst_opt', 'tst_unc_opt')


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
            myargs=(func, id, out_q, err_q, lock) + args
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
        return

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
                FUNC_MAX = np.float_(np.finfo(np.float_).max)
                return FUNC_MAX
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
        return

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

    def init(self, npop, xpar, step, seed):
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
        tmp = np.array(sorted(simp, key=lambda arg:arg[-1]))
        tmp.reshape(myshape)
        return tmp

    def sort(self):
        self.simplex = self.sort_me(self.simplex)

class SimplexNoStep(SimplexBase):

    def __init__(self, func, npop, xpar, xmin, xmax, step, seed, factor):
        SimplexBase.__init__(self, func, npop, xpar, xmin, xmax, step, seed, factor)
        return

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

    def __init__(self, func, npop, xpar, xmin, xmax, step, seed, factor):
        SimplexBase.__init__(self, func, npop, xpar, xmin, xmax, step, seed, factor)
        return

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

    def __init__(self, func, npop, xpar, xmin, xmax, step, seed, factor):
        SimplexBase.__init__(self, func, npop, xpar, xmin, xmax, step, seed, factor)
        return

    def init(self, npop, xpar, step, seed, factor):
        npar1 = self.npar + 1
        simplex = np.empty((npop, npar1))
        simplex[0][:-1] = np.copy(xpar)
        simplex = self.init_random_simplex(xpar, simplex, 1, npop, seed, factor)
        return self.eval_simplex(npop, simplex)


def Ackley(x):
    """Ackley(0, ..., 0) = 0"""
    n = x.shape[0]
    a = 20
    b = 0.2
    c = 2 * np.pi
    s1 = 0
    s2 = 0
    for ii in range(n):
        s1 += x[ii] * x[ii]
        s2 += np.cos(c * x[ii])
    tmp = -a*np.exp(-b*np.sqrt(s1/n))-np.exp(s2/n)+a+np.exp(1)
    return tmp

def Beale(x):
    """Beale(3, 0.5) = 0."""
    a = pow(1.5-x[0]*(1-x[1]), 2.0)
    b = pow(2.25-x[0]*(1-x[1]*x[1]), 2.0)
    c = pow(2.625-x[0]*(1-x[1]*x[1]*x[1]), 2.0)
    return a + b + c

def Bohachevsky1(x):
    """Bohchevsky1(0, 0) = 0"""
    return x[0]*x[0] + 2*x[1]*x[1] - 0.3*np.cos(3*np.pi*x[0]) - \
        0.4 * np.cos(4*np.pi*x[1]) + 0.7

def Bohachevsky2(x):
    """Bohachevsky2(0, 0) = 0"""
    return x[0]*x[0] + 2*x[1]*x[1] - \
        0.3*np.cos(3*np.pi*x[0])*np.cos(4*np.pi*x[1]) + 0.3

def Bohachevsky3(x):
    """Bohachevsky3(0, 0) = 0"""
    tmp = x[0]*x[0] + 2*x[1]*x[1] - \
        0.3*np.cos(3*np.pi*x[0]+4*np.pi*x[1]) + 0.3
    return tmp

def Booth(x):
    """Booth(1, 3) = 0"""
    return pow(x[0]+2*x[1]-7, 2.0) + pow(2*x[0]+x[1]-5, 2.0)

def BoxBetts(x):
    """BoxBetts(1, 10, 1) = 0"""
    fval = 0.0
    for ii in range(10):
        e0 = np.exp( -0.1 * ii * x[0] )
        e1 = np.exp( -0.1 * ii * x[1] )
        e2 = np.exp( -0.1 * ii ) - np.exp( - ii )
        tmp = e0 - e1 - e2 * x[2]
        fval += tmp * tmp
    return fval

def Branin(x):
    """Branin(-pi, 12.275) = 0.397887
    Branin(pi, 2.275) = 0.397887
    Branin(9.42478, 2.475)  = 0.397887"""
    return \
        pow(x[1]-(5.1/(4*np.pi*np.pi))*x[0]*x[0]+5*x[0]/np.pi-6,2.0) \
        + 10*(1-1/(8*np.pi))*np.cos(x[0])+10

def BrownBadlyScaled(x):
    """BrownBadlyScaled(1.0e6, 2.0e-6, ...,1.0e6, 2.0e-6) = 0"""
    n = len(x)
    fval = 0.0
    for ii in range(0, n, 2):
        fvec0 = x[ii] - 1.0e6
        fvec1 = x[ii + 1] - 2.0e-6
        fvec2 = x[ii] * x[ii + 1] - 2.0
        fval += fvec0 * fvec0 + fvec1 * fvec1 + fvec2 * fvec2
    return fval

def Colville(x):
    """Colville(1, 1, 1, 1) = 0"""
    if 0.0 == x[1]:
        return 1.0e35
    tmp = 100 * pow(x[0]*x[0]-x[1],2.0) + pow(x[0]-1,2.0) + \
        pow(x[2]-1, 2.0) + 90 * pow(x[2]*x[2]-x[3], 2.0) + \
        10.1 * (pow(x[1]-1, 2.0) + pow(x[3]-1,2.0)) + 19.8*(1./x[1])*(x[3]-1)
    return tmp

def DixonPrice(x):
    """DixonPrice(2^((2^i-2) / 2^i)) = 0"""
    npar = len(x)
    jj = range(2, npar + 1)
    x2 = 2 * x**2
    return sum(jj * (x2[1:] - x[:-1])**2) + (x[0] - 1)**2
    s = 0

def Easom(x):
    """Easom(pi, pi) = -1"""
    return -np.cos(x[0])*np.cos(x[1]) * \
        np.exp(-pow(x[0]-np.pi,2) - pow(x[1]-np.pi,2.0))

def FreudensteinRoth(x):
    """FreudensteinRoth(5, 4, 5, 4, ...., 5, 4) = 0"""
    n = len(x)
    tmp = 0
    for ii in range(0, n, 2):
        tmp += pow(-13 + x[ii] + ((5 - x[ii+1]) * x[ii+1] - 2) *
                   x[ii+1], 2.0) + \
                   pow(-29 + x[ii] + ((x[ii+1] + 1) *
                                      x[ii+1] - 14) * x[ii+1], 2.0)
    return tmp

def GoldsteinPrice(x):
    """GoldsteinPrice(0, 1) = 3"""
    a = 1 + pow(x[0]+x[1]+1,2) * \
        (19-14*x[0]+3*x[0]*x[0]-14*x[1]+6*x[0]*x[1]+3*x[1]*x[1])
    b = 30 + \
        pow(2*x[0]-3*x[1],2.0) * \
        (18-32*x[0]+12*x[0]*x[0]+48*x[1]-36*x[0]*x[1]+27*x[1]*x[1])
    return a * b

def Griewank(x):
    """Griewank(0, 0, ..., 0) = 0"""
    n = len(x)
    s = 0
    p = 1
    for ii in range(n):
        s += x[ii] * x[ii]
    for ii in range(n):
        p *= np.cos(x[ii]) / np.sqrt(ii + 1)
    return s / 4000.0 - p + 1

def Hump(x):
    """Hump((0.0898, -0.7126) = Hump(-0.0898, 0.7126) = 0"""
    return 1.0316285 + 4 * x[0] * x[0] - 2.1 * pow(x[0], 4) + \
        pow(x[0], 6) / 3 + x[0] * x[1] - 4 * x[1]* x[1] + \
        4 * pow(x[1], 4)

def Levy(x):
    """Levy(1, 1, ..., 1) = 0"""
    n = len(x)
    z = np.empty(n)
    for ii in range(n):
        z[ii] = 1+(x[ii]-1)/4
    s = np.sin(np.pi * z[0]) * np.sin(np.pi * z[0])
    for ii in range(n - 1):
        s = s + \
            pow(z[ii]-1, 2) * \
            (1 + 10 * pow(np.sin(np.pi * z[ii] + 1), 2.))
    return s + pow(z[n - 1]-1, 2) * \
        (1 + pow(np.sin(2 * np.pi * z[n - 1]) , 2))

def Matyas(x):
    """Matyas(0, 0) = 0"""
    return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1]

def McCormick(x):
    """McCormick( -0.547197553, -1.54719756 ) = -1.91"""
    a = x[0] + x[1]
    b = x[0] - x[1]
    fval = np.sin(a) + b * b - 1.5 * x[0] + 2.5 * x[1] + 1.0
    return fval

def Paviani(x):
    """Paviani( 9.35026583, 9.35026583, ..., 9.35026583, 9.35026583 ) = -45.7"""
    mul = 1.0
    fval = 0.0
    n = len(x)
    for ii in range(n):
        a = np.log(x[ii] - 2.0)
        b = np.log(10.0 - x[ii ])
        fval += a * a +  b * b
        mul *= x[ii]
    fval -= pow(mul, 0.2)
    return fval

def Rastrigin(x):
    """Rastrigin(0, 0, ..., 0) = 0"""
    if not np.isscalar(x[0]):
        N = len(x[0])
        return [10*N + sum(xi**2 - 10*np.cos(2*np.pi*xi)) for xi in x]
    N = len(x)
    return 10*N + sum(x**2 - 10*np.cos(2*np.pi*x))

def Rosenbrock(x):
    """Rosenbrock(1, 1, ..., 1) = 0"""
    x = np.asarray(x)
    val = np.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)
    return val

def Schwefel(x):
    """Schwefel(1, , 1, ...., 1) = 0"""
    n = len(x)
    s = sum(- x * np.sin(np.sqrt(abs(x))))
    return 418.9829 * n + s

def Shubert(x):
    """Shubert(x) = -186.7309"""
    s1 = 0
    s2 = 0
    for ii in range(5):
        s1 += (ii + 1) * np.cos((ii + 2) * x[0] + ii + 1)
        s2 += (ii + 1) * np.cos((ii + 2) * x[1] + ii + 1)
    return s1 * s2

def Sphere(x):
    """Sphere(0, 0, ..., 0) = 0"""
    n = len(x)
    s = 0
    for ii in range(n):
        s += x[ii] * x[ii]
    return s

def SumSquares(x):
    """SumSquares(0, 0, ..., 0) = 0"""
    n = len(x)
    s = 0
    for ii in range(n):
        s += (ii + 1) * x[ii] * x[ii]
    return s

def Trid(x):
    """Trid(n == 6) = -50, Trid(n = 10) = -200"""
    n = len(x)
    s1 = 0
    s2 = 0
    for ii in range(n):
        s1 += pow(x[ii] - 1, 2.0)

    for ii in range(1, n):
        s2 += x[ii] * x[ii - 1]
    return s1 - s2

def Zakharov(x):
    """Zakharov(0, 0, ..., 0) = 0"""
    n = len(x)
    jj = np.arange(1.0, n + 1)
    s2 = sum( jj * x ) / 2
    return sum(x**2) + s2**2 + s2**4

def pack_result(arg):
    par = arg[2:]
    result = np.asarray(arg[:2])
    result = np.append(result, par)
    return result

def myprint(name, func, result):
    try:
        if len(result) == 3:
            result = pack_result(result)
        else:
            result = pack_result(result[1:])
        nfev = int(result[0])
        fmin = result[1]
        par = result[2:]
        print(name, par, ' = ', fmin, 'in', nfev, 'nfevs')
        if func is None:
            return
        # print('autograd covariance:\n', autograd_covar(func, par))
    except NotImplementedError as nie:
        raise nie

def tst_unc_opt(algorithms, npar):
    """
    More, J.J., Garbow, B.S. and Hillstrom, K.E., Testing Unconstrained Optimization Software, ACM Trans. Math. Software 7 (1981), 17-41.
    """

    from sherpa.optmethods import _tstoptfct

    def tst_algo(opt, fcn, name, num):
        def func_wrapper(arg):
            return fcn(arg)[0]
        x0, xmin, xmax, fmin = _tstoptfct.init(name, num)
        result = opt(func_wrapper, x0, xmin, xmax)
        opt_name = opt.__class__.__name__
        print(opt_name, result[2], '=', result[1], 'in', result[0], 'nfevs')

    def rosenbrock(name, opt):
        tst_algo(opt, _tstoptfct.rosenbrock, name, npar)
    name = 'rosenbrock'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        rosenbrock(name, algo)

    def freudenstein_roth(name, opt):
        tst_algo(opt, _tstoptfct.freudenstein_roth, name, npar)
    name = 'freudenstein_roth'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        freudenstein_roth(name, algo)

    def powell_badly_scaled(name, opt):
        tst_algo(opt, _tstoptfct.powell_badly_scaled, name, npar)
    name = 'powell_badly_scaled'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        powell_badly_scaled(name, algo)

    def brown_badly_scaled(name, opt):
        tst_algo(opt, _tstoptfct.brown_badly_scaled, name, npar)
    name = 'brown_badly_scaled'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        brown_badly_scaled(name, algo)

    def beale(name, opt):
        tst_algo(opt, _tstoptfct.beale, name, npar)
    name = 'beale'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        beale(name, algo)

    def jennrich_sampson(name, opt):
        tst_algo(opt, _tstoptfct.jennrich_sampson, name, npar)
    name = 'jennrich_sampson'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        jennrich_sampson(name, algo)

    def helical_valley(name, opt):
        tst_algo(opt, _tstoptfct.helical_valley, name, 3)
    name = 'helical_valley'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        helical_valley(name, algo)

    def bard(name, opt):
        tst_algo(opt, _tstoptfct.bard, name, 3 * npar)
    name = 'bard'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3 * npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        bard(name, algo)

    def gaussian(name, opt):
        tst_algo(opt, _tstoptfct.gaussian, name, 3)
    name = 'gaussian'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        gaussian(name, algo)

    def meyer(name, opt):
        tst_algo(opt, _tstoptfct.meyer, name, 3)
    name = 'meyer'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        meyer(name, algo)

    def gulf_research_development(name, opt):
        tst_algo(opt, _tstoptfct.gulf_research_development, name, 3)
    name = 'gulf_research_development'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        gulf_research_development(name, algo)

    def box3d(name, opt):
        tst_algo(opt, _tstoptfct.box3d, name, 3)
    name = 'box3d'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 3)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        box3d(name, algo)

    def powell_singular(name, opt):
        tst_algo(opt, _tstoptfct.powell_singular, name, 4 * npar)
    name = 'powell_singular'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4 * npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        powell_singular(name, algo)

    def wood(name, opt):
        tst_algo(opt, _tstoptfct.wood, name, 4 * npar)
    name = 'wood'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4 * npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        wood(name, algo)

    def kowalik_osborne(name, opt):
        tst_algo(opt, _tstoptfct.kowalik_osborne, name, 4)
    name = 'kowalik_osborne'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        kowalik_osborne(name, algo)

    def brown_dennis(name, opt):
        tst_algo(opt, _tstoptfct.brown_dennis, name, 4)
    name = 'brown_dennis'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        brown_dennis(name, algo)

    def osborne1(name, opt):
        tst_algo(opt, _tstoptfct.osborne1, name, 5)
    name = 'osborne1'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 5)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        osborne1(name, algo)

    def biggs(name, opt):
        tst_algo(opt, _tstoptfct.biggs, name, 6)
    name = 'biggs'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 6)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        if algo.__class__.__name__ == 'Midnight':
            print('Minuit aborts skip test')
            continue
        biggs(name, algo)

    def osborne2(name, opt):
        tst_algo(opt, _tstoptfct.osborne2, name, 11)
    name = 'osborne2'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 11)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        osborne2(name, algo)

    def watson(name, opt):
        tst_algo(opt, _tstoptfct.watson, name, 6)
    name = 'watson'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 6)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        watson(name, algo)

    def penaltyI(name, opt):
        tst_algo(opt, _tstoptfct.penaltyI, name, 4)
    name = 'penaltyI'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        penaltyI(name, algo)

    def penaltyII(name, opt):
        tst_algo(opt, _tstoptfct.penaltyII, name, 4)
    name = 'penaltyII'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, 4)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        penaltyII(name, algo)

    def variably_dimensioned(name, opt):
        tst_algo(opt, _tstoptfct.variably_dimensioned, name, npar)
    name = 'variably_dimensioned'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        variably_dimensioned(name, algo)

    def trigonometric(name, opt):
        tst_algo(opt, _tstoptfct.trigonometric, name, npar)
    name = 'trigonometric'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        trigonometric(name, algo)

    def brown_almost_linear(name, opt):
        tst_algo(opt, _tstoptfct.brown_almost_linear, name, npar)
    name = 'brown_almost_linear'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        brown_almost_linear(name, algo)

    def discrete_boundary(name, opt):
        tst_algo(opt, _tstoptfct.discrete_boundary, name, npar)
    name = 'discrete_boundary'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        discrete_boundary(name, algo)

    def discrete_integral(name, opt):
        tst_algo(opt, _tstoptfct.discrete_integral, name, npar)
    name = 'discrete_integral'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        discrete_integral(name, algo)

    def broyden_tridiagonal(name, opt):
        tst_algo(opt, _tstoptfct.broyden_tridiagonal, name, npar)
    name = 'broyden_tridiagonal'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        broyden_tridiagonal(name, algo)

    def broyden_banded(name, opt):
        tst_algo(opt, _tstoptfct.broyden_banded, name, npar)
    name = 'broyden_banded'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        broyden_banded(name, algo)

    def linear_fullrank(name, opt):
        tst_algo(opt, _tstoptfct.linear_fullrank, name, npar)
    name = 'linear_fullrank'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        linear_fullrank(name, algo)

    def linear_fullrank1(name, opt):
        tst_algo(opt, _tstoptfct.linear_fullrank, name, npar)
    name = 'linear_fullrank1'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        linear_fullrank1(name, algo)

    def linear_fullrank0cols0rows(name, opt):
        tst_algo(opt, _tstoptfct.linear_fullrank0cols0rows, name, npar)
    name = 'linear_fullrank0cols0rows'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        if algo.__class__.__name__ == 'Midnight':
            print('Minuit aborts skip test')
            continue
        linear_fullrank0cols0rows(name, algo)

    def chebyquad(name, opt):
        tst_algo(opt, _tstoptfct.chebyquad, name, npar)
    name = 'chebyquad'
    x0, xmin, xmax, fmin = _tstoptfct.init(name, npar)
    print('\n', name, ' fmin =', fmin)
    for algo in algorithms:
        chebyquad(name, algo)


def tst_opt(algorithms, npar):

    def tst_algos(func, x0, xmin, xmax):
        print('\n', func.__doc__)
        for algo in algorithms:
            result = algo(func, x0, xmin, xmax)
            myprint(algo.__class__.__name__, func, result)
        return

    xmin = npar * [-32.768]
    xmax = npar * [32.768]
    x0 = npar * [12.3]
    tst_algos(Ackley, x0, xmin, xmax)

    xmin = [-4.5, -4.5]
    xmax = [4.5, 4.5]
    x0 = [-1.0, 2.0]
    tst_algos(Beale, x0, xmin, xmax)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [-12, 10]
    tst_algos(Bohachevsky1, x0, xmin, xmax)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [12, 10]
    tst_algos(Bohachevsky2, x0, xmin, xmax)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [-61.2, 51.0]
    tst_algos(Bohachevsky3, x0, xmin, xmax)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-6.2, 5.0]
    tst_algos(Booth, x0, xmin, xmax)

    xmin = [0.9, 9.0, 0.9]
    xmax = [1.2, 11.2, 1.2]
    x0      = [(xmin[0] + xmax[0]) * 0.5, (xmin[1] + xmax[1]) * 0.5,
               (xmin[2] + xmax[2]) * 0.5]
    tst_algos(BoxBetts, x0, xmin, xmax)

    xmin = [-5, 0]
    xmax = [10, 15]
    x0	 = [-3.2, 5.0]
    tst_algos(Branin, x0, xmin, xmax)

    xmin = npar * [-1.0e2]
    xmax = npar * [1.0e9]
    x0 = npar * [1]
    tst_algos(BrownBadlyScaled, x0, xmin, xmax)

    xmin = [-10, -10, -10, -10]
    xmax = [10, 10, 10, 10.]
    x0	 = [-3.2, -5.0, -6.0, -1.0]
    tst_algos(Colville, x0, xmin, xmax)

    xmin = npar * [-10]
    xmax = npar * [10]
    x0 = npar * [-1.]
    tst_algos(DixonPrice, x0, xmin, xmax)

    xmin = [-100, -100]
    xmax = [100, 100]
    # x0	 = [5.5, 5.5]
    x0	 = [25., 25.]
    tst_algos(Easom, x0, xmin, xmax)

    xmin = npar * [-1000, -1000]
    xmax = npar * [1000, 1000]
    x0       = npar * [0.5, -2]
    tst_algos(FreudensteinRoth, x0, xmin, xmax)

    xmin = [-2, -2]
    xmax = [2, 2]
    x0	 = [-1, 1]
    tst_algos(GoldsteinPrice, x0, xmin, xmax)

    xmin = npar * [-600]
    xmax = npar * [600]
    x0 = npar * [-100.]
    tst_algos(Griewank, x0, xmin, xmax)

    xmin = [-5, -5]
    xmax = [5, 5]
    x0	 = [-3.2, 5.0]
    tst_algos(Hump, x0, xmin, xmax)

    xmin = npar * [-10]
    xmax = npar * [10]
    x0 = npar * [-5.]
    tst_algos(Levy, x0, xmin, xmax)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-3.2, 5.0]
    tst_algos(Matyas, x0, xmin, xmax)

    xmin = [-1.5, -3.0]
    xmax = [4.0, 4.0]
    x0     = [0.0, 0.0]
    tst_algos(McCormick, x0, xmin, xmax)

    xmin = npar * [2.001]
    xmax = npar * [9.999]
    x0	 = npar * [5.0]
    tst_algos(Paviani, x0, xmin, xmax)

    xmin = npar * [-5.12]
    xmax = npar * [5.12]
    x0	 = npar * [-2.0]
    tst_algos(Rastrigin, x0, xmin, xmax)

    xmin = npar * [-1000, -1000]
    xmax = npar * [1000, 1000]
    x0	 = npar * [-1.2, 1.0]
    tst_algos(Rosenbrock, x0, xmin, xmax)

    xmin = npar * [-500]
    xmax = npar * [500]
    x0	 = npar * [-200]
    tst_algos(Schwefel, x0, xmin, xmax)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-2.0, 5.0]
    tst_algos(Shubert, x0, xmin, xmax)

    xmin = npar * [-5.12]
    xmax = npar * [5.12]
    x0	 = npar * [-2.0]
    tst_algos(Sphere, x0, xmin, xmax)

    xmin = npar * [-10]
    xmax = npar * [10]
    x0	 = npar * [-2.0]
    tst_algos(SumSquares, x0, xmin, xmax)

    if npar == 6 or npar == 10:
        xmin = npar * [- npar * npar]
        xmax = npar * [npar * npar]
        x0 = npar * [10]
        tst_algos(Trid, x0, xmin, xmax)

    xmin = npar * [-5, -5]
    xmax = npar * [10, 10]
    x0   = npar * [0.5, -2]
    tst_algos(Zakharov, x0, xmin, xmax)

if '__main__' == __name__:

    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-N", "--npar", dest="npar", default=10,
                      type=int, help="set npar")
    (options, args) = parser.parse_args()
    npar = options.npar

    x0 = np.array(npar * [-1.2, 1.0])
    xmin = npar * [-1000, -1000]
    xmax = npar * [1000, 1000]
    factor = 10
    seed = 234
    simp = SimplexNoStep(Rosenbrock, len(x0) + 1, x0, xmin, xmax, None,
                         seed, factor)
    print('simp =\n', simp.simplex)

    simp = SimplexStep(Rosenbrock, len(x0) + 2, x0, xmin, xmax, x0 + 1.2,
                       seed, factor)
    print('simp =\n', simp.simplex)

    simp = SimplexRandom(Rosenbrock, len(x0) + 5, x0, xmin, xmax, x0 + 1.2,
                         seed, factor)
    print('simp =\n', simp.simplex)
