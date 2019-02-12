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

import multiprocessing
import numpy
import random

from sherpa.utils import Knuth_close, _multi, _ncpus, run_tasks

__all__ = ('TraceCalls', 'Opt', 'Polytope', 'Simplex', 'MyNcores', 'tst_opt')


#####################################################
import sys
from functools import wraps
class TraceCalls(object):
    """ Use as a decorator on functions that should be traced. Several
        functions can be decorated - they will all be indented according
        to their call depth.
    """
    def __init__(self, stream=sys.stdout, indent_step=2, show_ret=False):
        self.stream = stream
        self.indent_step = indent_step
        self.show_ret = show_ret

        # This is a class attribute since we want to share the indentation
        # level between different traced functions, in case they call
        # each other.
        TraceCalls.cur_indent = 0

    def __call__(self, fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            indent = ' ' * TraceCalls.cur_indent
            argstr = ', '.join(
                [repr(a) for a in args] +
                ["%s=%s" % (a, repr(b)) for a, b in kwargs.items()])
            self.stream.write('%s%s(%s)\n' % (indent, fn.__name__, argstr))

            TraceCalls.cur_indent += self.indent_step
            ret = fn(*args, **kwargs)
            TraceCalls.cur_indent -= self.indent_step

            if self.show_ret:
                self.stream.write('%s--> %s\n' % (indent, ret))
            return ret
        return wrapper
#####################################################


class MyNcores:

    def __init__(self):
        if _multi is False:
            raise TypeError("multicores not available")

    def calc(self, funcs, numcores, *args, **kwargs):

        for func in funcs:
            if not callable(func):
                raise TypeError("input func '%s' is not callable" % repr(func))

        # fix me
        if numcores is None:
            numcores = _ncpus
        numcores = min(numcores, len(funcs))
        
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
            myargs=(func, id, *args, out_q, err_q, lock)
            try:
                procs.append(multiprocessing.Process(target=self.my_worker,
                                                     args=myargs))
            except NotImplementedError as nie:
                raise nie
        return run_tasks(procs, err_q, out_q, numcores)

    def my_worker(self, *args):
        raise NotImplementedError("my_worker has not been implemented")

    
class Opt:

    def __init__(self, func, xmin, xmax):
        self.npar = len(xmin)
        self.nfev, self.func = \
            self.func_counter_bounds_wrappers(func, self.npar, xmin, xmax)
        self.xmin = numpy.asarray(xmin)
        self.xmax = numpy.asarray(xmax)
        return

    def _outside_limits(self, x, xmin, xmax):
        return (numpy.any(x < xmin) or numpy.any(x > xmax))

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
            xmin = numpy.asarray(xmin)
            xmax = numpy.asarray(xmax)
        elif xmin is not None:
            xmin = numpy.asarray(xmin)
            xmax = numpy.asarray([numpy.inf for ii in xmin])
        elif xmax is not None:
            xmax = numpy.asarray(xmax)
            xmin = numpy.asarray([- numpy.inf for ii in xmax])
        else:
            xmin = numpy.asarray([- numpy.inf for ii in range(npar)])
            xmax = numpy.asarray([numpy.inf for ii in range(npar)])

        def func_bounds_wrapper(x, *args):
            if self._outside_limits(x, xmin, xmax):
                FUNC_MAX = numpy.float_(numpy.finfo(numpy.float_).max)
                return FUNC_MAX
            return func(x, *args)

        return func_bounds_wrapper

    def func_counter(self, func):
        """A function wrapper to count the number of times being called"""
        nfev = [0]
        def func_counter_wrapper(x, *args):
            nfev[0] += 1
            tmp = func(x, *args)
            # print(func.__name__, x, '=', tmp, 'at', nfev[0])
            return tmp
        return nfev, func_counter_wrapper

    def func_counter_bounds_wrappers(self, func, npar, xmin, xmax):
        """Wraps the calls to func_counter then func_bounds"""
        nfev, afunc = self.func_counter(func)
        myfunc = self.func_bounds(afunc, npar, xmin, xmax)
        return nfev, myfunc

    def sort(self, pop):
        myshape = pop.shape
        tmp = numpy.array(sorted(pop, key=lambda arg:arg[-1]))
        tmp.reshape(myshape)
        return tmp


class Polytope:

    def __init__(self, func, num, xpar, xmin, xmax, step, seed=123):
        self.func = func
        self.xmin = xmin
        self.xmax = xmax
        self.npar = len(xpar)
        self.polytope = self.init(num, xpar, xmin, xmax, step, seed)
        return

    def __getitem__(self, index):
        return self.polytope[index]

    def __setitem__(self, index, val):
        self.polytope[index] = val

    def calc_centroid(self):
        return numpy.mean(self.polytope[:-1, :], 0)

    def init(self, num, xpar, xmin, xmax, step, seed):
        npar_plus_1 = self.npar + 1
        simplex = numpy.zeros((num, npar_plus_1))

        simplex[0][:-1] = xpar[:]
        if step is None:
            for ii in range(self.npar):
                tmp = xpar[:]
                if 0.0 == ii:
                    tmp[ii] = 2.5e-4
                else:
                    tmp[ii] *= 1.05
                simplex[ii + 1][:-1] = tmp[:]
        else:
            for ii in range(self.npar):
                tmp = xpar[ii] + step[ii]
                simplex[ii + 1][:-1] = tmp

        random.seed(seed)
        tmp = 10
        for ii in range(npar_plus_1, num):
            simplex[ii][:-1] = \
                numpy.array([random.uniform(max(xmin[jj],
                                                xpar[jj]-tmp*abs(xpar[jj])),
                                            min(xmax[jj],
                                                xpar[jj]+tmp*abs(xpar[jj]))) \
                                   for jj in range(self.npar)])

        for ii in range(num):
            simplex[ii][-1] = self.func(simplex[ii][:-1])
        return self.mysort(simplex)

    # @TraceCalls(indent_step=4, show_ret=True)
    def check_convergence(self, ftol, method):

        def are_func_vals_close_enough(tolerance):
            smallest_fct_val = self.polytope[0, -1]
            largest_fct_val = self.polytope[-1, -1]
            return Knuth_close(smallest_fct_val, largest_fct_val,
                                tolerance)

        def is_fct_stddev_small_enough(tolerance):
            fval_std = numpy.std([col[-1] for col in self.polytope])
            if fval_std < ftol:
                return True
            return False

        def is_max_length_small_enough(tolerance):
            """

               max  || x  - x  || <= tol max(1.0, || x ||)
                        i    0                         0

            where 1 <= i <= n"""

            max_xi_x0 = -1.0
            x0 = self.polytope[0, :-1]

            npar_plus_1 = len(x0) + 1
            for ii in range(1, npar_plus_1):
                xi = self.polytope[ii, :-1]
                xi_x0 = xi - x0
                max_xi_x0 = max(max_xi_x0, numpy.dot(xi_x0, xi_x0))
            return max_xi_x0 <= tolerance * max(1.0, numpy.dot(x0, x0))

        if 0 == method:
            if is_max_length_small_enough(ftol):
                return True
        elif 2 == method:
            if False == is_max_length_small_enough(ftol):
                return False
            stddev = is_fct_stddev_small_enough(ftol)
            fctval = are_func_vals_close_enough(ftol)
            return stddev and fctval
        else:
            if False == is_max_length_small_enough(ftol):
                return False
            stddev = is_fct_stddev_small_enough(ftol)
            fctval = are_func_vals_close_enough(ftol)
            return stddev or fctval
        return False

        num = 2.0 * abs(self.polytope[0, -1] - self.polytope[-1, -1])
        denom = abs(self.polytope[0, -1]) + abs(self.polytope[-1, -1]) + 1.0
        if (num / denom > ftol):
            return False
        func_vals = [col[-1] for col in self.polytope]
        if numpy.std(func_vals) > ftol:
            return False
        return True

    def mysort(self, pop):
        myshape = pop.shape
        tmp = numpy.array(sorted(pop, key=lambda arg:arg[-1]))
        tmp.reshape(myshape)
        return tmp

    def sort(self):
        self.polytope = self.mysort(self.polytope)

        
class Simplex(Polytope):

    def __init__(self, func, num, xpar, xmin, xmax, step):
        Polytope.__init__(self, func, num, xpar, xmin, xmax, step)
        return

    ##@TraceCalls(indent_step=4, show_ret=True)
    def calc_centroid(self):
        return numpy.mean(self.polytope[:self.npar, :], 0)


    ##@TraceCalls(indent_step=4, show_ret=True)
    def move_vertex(self, centroid, coef):
        vertex = (1.0 + coef) * centroid - coef * self.polytope[self.npar]
        vertex[-1] = self.func(vertex[:-1])
        return vertex

    ##@TraceCalls(indent_step=4, show_ret=True)
    def shrink(self, shrink_coef):
        npars_plus_1 = self.npar + 1
        for ii in range(1, npars_plus_1):
            self.polytope[ii] = \
                self.polytope[0] + shrink_coef * \
                (self.polytope[ii] - self.polytope[0])
            self.polytope[ii, -1] = self.func(self.polytope[ii, :-1])

def Ackley(x):
    """Ackley(0, ..., 0) = 0"""
    n = len(x); a = 20; b = 0.2; c = 2*numpy.pi; s1 = 0; s2 = 0;
    for ii in range(n):
        s1 += x[ii] * x[ii]
        s2 += numpy.cos(c * x[ii])
    tmp = -a*numpy.exp(-b*numpy.sqrt(s1/n))-numpy.exp(s2/n)+a+numpy.exp(1)
    return tmp

def Beale(x):
    """Beale(3, 0.5) = 0."""
    a = pow(1.5-x[0]*(1-x[1]), 2.0)
    b = pow(2.25-x[0]*(1-x[1]*x[1]), 2.0)
    c = pow(2.625-x[0]*(1-x[1]*x[1]*x[1]), 2.0)
    return a + b + c

def Bohachevsky1(x):
    """Bohchevsky1(0, 0) = 0"""
    return x[0]*x[0] + 2*x[1]*x[1] - 0.3*numpy.cos(3*numpy.pi*x[0]) - \
        0.4 * numpy.cos(4*numpy.pi*x[1]) + 0.7

def Bohachevsky2(x):
    """Bohachevsky2(0, 0) = 0"""
    return x[0]*x[0] + 2*x[1]*x[1] - \
        0.3*numpy.cos(3*numpy.pi*x[0])*numpy.cos(4*numpy.pi*x[1]) + 0.3

def Bohachevsky3(x):
    """Bohachevsky3(0, 0) = 0"""
    tmp = x[0]*x[0] + 2*x[1]*x[1] - \
        0.3*numpy.cos(3*numpy.pi*x[0]+4*numpy.pi*x[1]) + 0.3
    return tmp

def Booth(x):
    """Booth(1, 3) = 0"""
    return pow(x[0]+2*x[1]-7, 2.0) + pow(2*x[0]+x[1]-5, 2.0)

def BoxBetts(x):
    """BoxBetts(1, 10, 1) = 0"""
    fval = 0.0
    for ii in range(10):
        e0 = numpy.exp( -0.1 * ii * x[0] )
        e1 = numpy.exp( -0.1 * ii * x[1] )
        e2 = numpy.exp( -0.1 * ii ) - numpy.exp( - ii )
        tmp = e0 - e1 - e2 * x[2]
        fval += tmp * tmp
    return fval

def Branin(x):
    """Branin(-pi, 12.275) = 0.397887
    Branin(pi, 2.275) = 0.397887
    Branin(9.42478, 2.475)  = 0.397887"""
    return \
        pow(x[1]-(5.1/(4*numpy.pi*numpy.pi))*x[0]*x[0]+5*x[0]/numpy.pi-6,2.0) \
        + 10*(1-1/(8*numpy.pi))*numpy.cos(x[0])+10

def BrownBadlyScaled(x):
    """BrownBadlyScaled(1.0e6, 2.0e-6, ...,1.0e6, 2.0e-6) = 0"""
    n = len(x)
    fval = 0.0
    for ii in range(0, n, 2):
        fvec0 = x[ii] - 1.0e6;
        fvec1 = x[ii + 1] - 2.0e-6;
        fvec2 = x[ii] * x[ii + 1] - 2.0;
        fval += fvec0 * fvec0 +fvec1 * fvec1 + fvec2 * fvec2
    return fval

#@TraceCalls(indent_step=4, show_ret=True)
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
    return -numpy.cos(x[0])*numpy.cos(x[1]) * \
        numpy.exp(-pow(x[0]-numpy.pi,2) - pow(x[1]-numpy.pi,2.0))

def FreudensteinRoth(x):
    """FreudensteinRoth(5, 4, 5, 4, ...., 5, 4) = 0"""
    n = len(x)
    tmp = 0
    for ii in range(0, n, 2):
        tmp += pow(-13 + x[ii] + ((5 - x[ii+1]) * x[ii+1] - 2) * \
                        x[ii+1], 2.0) + \
                        pow(-29 + x[ii] + ((x[ii+1] + 1) * \
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
        p *= numpy.cos(x[ii]) / numpy.sqrt(ii + 1)
    return s / 4000.0 - p + 1

def Hump(x):
    """Hump((0.0898, -0.7126) = Hump(-0.0898, 0.7126) = 0"""
    return 1.0316285 + 4 * x[0] * x[0] - 2.1 * pow(x[0], 4) + \
        pow(x[0], 6) / 3 + x[0] * x[1] - 4 * x[1]* x[1] + \
        4 * pow(x[1], 4)

def Levy(x):
    """Levy(1, 1, ..., 1) = 0"""
    n = len(x)
    z = numpy.empty(n)
    for ii in range(n):
        z[ii] = 1+(x[ii]-1)/4
    s = numpy.sin(numpy.pi * z[0]) * numpy.sin(numpy.pi * z[0])
    for ii in range(n - 1):
        s = s + \
            pow(z[ii]-1, 2) * \
            (1 + 10 * pow(numpy.sin(numpy.pi * z[ii] + 1), 2.))
    return s + pow(z[n - 1]-1, 2) * \
        (1 + pow(numpy.sin(2 * numpy.pi * z[n - 1]) , 2))

def Matyas(x):
    """Matyas(0, 0) = 0"""
    return 0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1]

def McCormick(x):
    """McCormick( -0.547197553, -1.54719756 ) = -1.91"""
    a = x[0] + x[1]
    b = x[0] - x[1]
    fval = numpy.sin(a) + b * b - 1.5 * x[0] + 2.5 * x[1] + 1.0
    return fval

def Paviani(x):
    """Paviani( 9.35026583, 9.35026583, ..., 9.35026583, 9.35026583 ) = -45.7"""
    mul = 1.0
    fval = 0.0
    n = len(x)
    for ii in range(n):
        a = numpy.log(x[ii] - 2.0)
        b = numpy.log(10.0 - x[ii ])
        fval += a * a +  b * b
        mul *= x[ii]
    fval -= pow(mul, 0.2)
    return fval

def Rastrigin(x):
    """Rastrigin(0, 0, ..., 0) = 0"""
    if not numpy.isscalar(x[0]):
        N = len(x[0])
        return [10*N + sum(xi**2 - 10*numpy.cos(2*numpy.pi*xi)) for xi in x]
    N = len(x)
    return 10*N + sum(x**2 - 10*numpy.cos(2*numpy.pi*x))


    n = len(x)
    s = 0
    for ii in range(n):
        s += x[ii] * x[ii] - 10.0 * numpy.cos(2.0 * numpy.pi * x[ii])
    return 10 * n + s

def Rosenbrock(x):
    """Rosenbrock(1, 1, ..., 1) = 0"""
    x = numpy.asarray(x)
    val = numpy.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0,axis=0)
    return val

def Schwefel(x):
    """Schwefel(1, , 1, ...., 1) = 0"""
    n = len(x)
    s = sum(- x * numpy.sin(numpy.sqrt(abs(x))))
    return 418.9829 * n + s

def Shubert(x):
    """Shubert(x) = -186.7309"""
    s1 = 0
    s2 = 0
    for ii in range(5):
        s1 += (ii + 1) * numpy.cos((ii + 2) * x[0] + ii + 1)
        s2 += (ii + 1) * numpy.cos((ii + 2) * x[1] + ii + 1)
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
    jj = numpy.arange(1.0, n + 1)
    s2 = sum( jj * x ) / 2
    return sum(x**2) + s2**2 + s2**4

def pack_result(arg):
    par = arg[2:]
    result = numpy.asarray(arg[:2])
    result = numpy.append(result, par)
    return result

def myprint(name, result):
    try:
        if len(result) == 3:
            result = pack_result(result)
        else:
            result = pack_result(result[1:])
        print(name, result[2:], ' = ', result[1], 'in', int(result[0]), 'nfevs')
    except NotImplementedError as nie:
        raise nie
    
def tst_opt(algorithms, num):

    from sherpa.optmethods import _tstoptfct
    from sherpa.optmethods import optfcts

    xmin = num * [-32.768]
    xmax = num * [32.768]
    x0 = num * [12.3]
    print('\n', Ackley.__doc__)
    for algorithm in algorithms:
        result = algorithm(Ackley, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)
        
    xmin = [-4.5, -4.5]
    xmax = [4.5, 4.5]
    x0 = [-1.0, 2.0]
    print('\n', Beale.__doc__)
    for algorithm in algorithms:
        result = algorithm(Beale, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [-12, 10]
    print('\n', Bohachevsky1.__doc__)
    for algorithm in algorithms:
        result = algorithm(Bohachevsky1, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [12, 10]
    print('\n', Bohachevsky2.__doc__)
    for algorithm in algorithms:
        result = algorithm(Bohachevsky2, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [-61.2, 51.0]
    print('\n', Bohachevsky3.__doc__)
    for algorithm in algorithms:
        result = algorithm(Bohachevsky3, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-6.2, 5.0]
    print('\n', Booth.__doc__)
    for algorithm in algorithms:
        result = algorithm(Booth, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [0.9, 9.0, 0.9]
    xmax = [1.2, 11.2, 1.2]
    x0      = [(xmin[0] + xmax[0]) * 0.5, (xmin[1] + xmax[1]) * 0.5,
               (xmin[2] + xmax[2]) * 0.5]
    print('\n', BoxBetts.__doc__)
    for algorithm in algorithms:
        result = algorithm(BoxBetts, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-5, 0]
    xmax = [10, 15]
    x0	 = [-3.2, 5.0]
    print('\n', Branin.__doc__)
    for algorithm in algorithms:
        result = algorithm(Branin, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-1.0e2]
    xmax = num * [1.0e9]
    x0 = num * [1]
    print('\n', BrownBadlyScaled.__doc__)
    for algorithm in algorithms:
        result = algorithm(BrownBadlyScaled, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)
        
    xmin = [-10, -10, -10, -10]
    xmax = [10, 10, 10, 10.]
    x0	 = [-3.2, -5.0, -6.0, -1.0]
    print('\n', Colville.__doc__)
    for algorithm in algorithms:
        result = algorithm(Colville, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-10]
    xmax = num * [10]
    x0 = num * [-1.]
    print('\n', DixonPrice.__doc__)
    for algorithm in algorithms:
        result = algorithm(DixonPrice, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-100, -100]
    xmax = [100, 100]
    x0	 = [25.0, 25.0]
    print('\n', Easom.__doc__)
    for algorithm in algorithms:
        result = algorithm(Easom, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-1000, -1000]
    xmax = num * [1000, 1000]
    x0       = num * [0.5, -2]
    print('\n', FreudensteinRoth.__doc__)
    for algorithm in algorithms:
        result = algorithm(FreudensteinRoth, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-2, -2]
    xmax = [2, 2]
    x0	 = [-1, 1]
    print('\n', GoldsteinPrice.__doc__)
    for algorithm in algorithms:
        result = algorithm(GoldsteinPrice, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-600]
    xmax = num * [600]
    x0 = num * [-100.]
    print('\n', Griewank.__doc__)
    for algorithm in algorithms:
        result = algorithm(Griewank, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-5, -5]
    xmax = [5, 5]
    x0	 = [-3.2, 5.0]
    print('\n', Hump.__doc__)
    for algorithm in algorithms:
        result = algorithm(Hump, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-10]
    xmax = num * [10]
    x0 = num * [-5.]
    print('\n', Levy.__doc__)
    for algorithm in algorithms:
        result = algorithm(Levy, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-3.2, 5.0]
    print('\n', Matyas.__doc__)
    for algorithm in algorithms:
        result = algorithm(Matyas, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = [-1.5, -3.0]
    xmax = [4.0, 4.0]
    x0     = [0.0, 0.0]
    print('\n', McCormick.__doc__)
    for algorithm in algorithms:
        result = algorithm(McCormick, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [2.001]
    xmax = num * [9.999]
    x0	 = num * [5.0]
    print('\n', Paviani.__doc__)
    for algorithm in algorithms:
        result = algorithm(Paviani, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-5.12]
    xmax = num * [5.12]
    x0	 = num * [-2.0]
    print('\n', Rastrigin.__doc__)
    for algorithm in algorithms:
        result = algorithm(Rastrigin, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-1000, -1000]
    xmax = num * [1000, 1000]
    x0	 = num * [-1.2, 1.0]
    print('\n', Rosenbrock.__doc__)
    for algorithm in algorithms:
        result = algorithm(Rosenbrock, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    # xmin = num * [-500]
    # xmax = num * [500]
    # x0	 = num * [-200]
    # print('\n', Schwefel.__doc__)
    # for algorithm in algorithms:
    #     result = algorithm(Schwefel, x0, xmin, xmax)
    #     myprint(algorithm.__class__.__name__, result)

    xmin = [-10, -10]
    xmax = [10, 10]
    x0	 = [-2.0, 5.0]
    print('\n', Shubert.__doc__)
    for algorithm in algorithms:
        result = algorithm(Shubert, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-5.12]
    xmax = num * [5.12]
    x0	 = num * [-2.0]
    print('\n', Sphere.__doc__)
    for algorithm in algorithms:
        result = algorithm(Sphere, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    xmin = num * [-10]
    xmax = num * [10]
    x0	 = num * [-2.0]
    print('\n', SumSquares.__doc__)
    for algorithm in algorithms:
        result = algorithm(SumSquares, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    if num == 6 or num == 10:
        xmin = num * [- num * num]
        xmax = num * [num * num]
        x0 = num * [10]
        print('\n', Trid.__doc__)
        for algorithm in algorithms:
            result = algorithm(Trid, x0, xmin, xmax)
            myprint(algorithm.__class__.__name__, result)

    xmin = num * [-5, -5]
    xmax = num * [10, 10]
    x0   = num * [0.5, -2]
    print('\n', Zakharov.__doc__)
    for algorithm in algorithms:
        result = algorithm(Zakharov, x0, xmin, xmax)
        myprint(algorithm.__class__.__name__, result)

    npar=18
    x0, xmin, xmax, fmin = _tstoptfct.init('broyden_banded', npar)
    print('\nfmin =', fmin)
    def func(arg):
        return _tstoptfct.broyden_banded(arg)[0]
    for algo in algorithms:
        result = algo(func, x0, xmin, xmax)
        print('f', result[2], '=', result[1], 'at nfev =', result[0])

    npar=18
    x0, xmin, xmax, fmin = _tstoptfct.init('broyden_tridiagonal', npar)
    print('\nfmin =', fmin)
    def func(arg):
        return _tstoptfct.broyden_tridiagonal(arg)[0]
    for algo in algorithms:
        result = algo(func, x0, xmin, xmax)
        print('f', result[2], '=', result[1], 'at nfev =', result[0])

    npar = 11
    x0, xmin, xmax, fmin = _tstoptfct.init('osborne2', npar)
    print('\nfmin =', fmin)
    def func(arg):
        return _tstoptfct.osborne2(arg)[0]
    for algo in algorithms:
        result = algo(func, x0, xmin, xmax)
        print('f', result[2], '=', result[1], 'at nfev =', result[0])        
