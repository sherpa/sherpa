#
#  Copyright (C) 2019- 2021, 2023, 2024
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

import numpy as np

from sherpa.utils.parallel import ncpus
from sherpa.optmethods import _saoopt
from sherpa.optmethods.opt import MyNcores, Opt, SimplexNoStep, SimplexStep, \
    SimplexRandom

__all__ = ('ncoresNelderMead', )

EPSILON = np.float64(np.finfo(np.float32).eps)


class MyNelderMead(Opt):

    def __init__(self, fcn, xmin, xmax):
        Opt.__init__(self, fcn, xmin, xmax)
        self.expansion_coef = 2.0          # chi
        self.contraction_coef = 0.5          # gamma
        self.reflection_coef = 1.0          # rho
        self.shrink_coef = 0.5          # sigma
        self.simplex = None

    def __call__(self, xpar, maxnfev, tol, step, finalsimplex, verbose):

        npar = len(xpar)
        simplex = SimplexStep(func=self.func, npop=npar + 1,
                              xpar=xpar, xmin=self.xmin,
                              xmax=self.xmax, step=step, seed=None,
                              factor=None)
        return self.optimize(xpar, simplex, maxnfev, tol,
                             finalsimplex, verbose)

    def contract_in_out(self, simplex, centroid, reflection_pt, rho_gamma,
                        contraction_coef, badindex, maxnfev, verbose):

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

    def optimize(self, xpar, simplex, maxnfev, tol, finalsimplex, verbose):

        rho_chi = self.reflection_coef * self.expansion_coef
        rho_gamma = self.reflection_coef * self.contraction_coef
        bad_index = len(xpar)

        while self.nfev < maxnfev:

            simplex.sort()
            if verbose > 2:
                print(f'f{simplex[0, :-1]}={simplex[0, -1]:e}')

            centroid = simplex.calc_centroid()

            if simplex.check_convergence(tol, finalsimplex):
                break

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

            else:
                if self.contract_in_out(simplex, centroid, reflection_pt,
                                        rho_gamma, self.contraction_coef,
                                        bad_index, maxnfev, verbose):
                    simplex.shrink(self.shrink_coef)

        best_vertex = simplex[0]
        best_par = best_vertex[:-1]
        best_val = best_vertex[-1]
        return self.nfev, best_val, best_par


class NelderMeadBase:

    def __init__(self):
        self.nfev = 0
        self.fmin = np.inf
        self.par = np.nan

        # TODO: this should either be configurable or something that is
        # only changed when needed. See issue #1238
        #
        np.seterr(over='ignore', divide='ignore', under='ignore',
                  invalid='ignore')

    def __call__(self, fcn, xpar, xmin, xmax, tol=1.0e-6,  maxnfev=None,
                 step=None, finalsimplex=1, verbose=0):
        num = len(xpar)
        return self.nfev, self.fmin, num * [self.par]

    def get_maxnfev(self, maxnfev, npar):
        if maxnfev is None:
            return 512 * npar

        return maxnfev


class NelderMead0(NelderMeadBase):

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6,  maxnfev=None, step=None,
                 finalsimplex=1, verbose=0):
        return self.neldermead0(fcn, x, xmin, xmax, step, finalsimplex, maxnfev,
                                tol, verbose)

    def calc_step(self, x):
        return 1.2 * x

    def neldermead0(self, fcn, x0, xmin, xmax, step=None, finalsimplex=1,
                    maxnfev=None, tol=1.0e-6, verbose=0):
        x0 = np.asarray(x0)
        maxnfev = self.get_maxnfev(maxnfev, len(x0))

        my_nm = MyNelderMead(fcn, xmin, xmax)
        if step is None:
            step = self.calc_step(x0)

        return my_nm(x0, maxnfev, tol, step, finalsimplex, verbose)


class NelderMead1(NelderMead0):

    def calc_step(self, x):
        return x + 1.2


class NelderMead2(NelderMead0):

    def calc_step(self, x):
        return abs(x)


class NelderMead3(NelderMead0):

    def __call__(self, fcn, x0, xmin, xmax, tol=EPSILON,  maxnfev=None,
                 step=None, finalsimplex=None, verbose=0):

        # Avoid having a mutable argument
        if finalsimplex is None:
            finalsimplex = [0, 1, 1]

        x0 = np.asarray(x0)
        n = len(x0)
        if step is None:
            step = n * [1.2]
        maxnfev = self.get_maxnfev(maxnfev, n)
        init = 0
        par, fmin, nfev, err = \
            _saoopt.neldermead(verbose, maxnfev, init, finalsimplex, tol, step,
                               xmin, xmax, x0, fcn)
        return nfev, fmin, par


class NelderMead4(NelderMead0):

    def __call__(self, fcn, x0, xmin, xmax, tol=EPSILON,  maxnfev=None,
                 step=None, finalsimplex=None, verbose=0, reflect=True):

        # Avoid having a mutable argument
        if finalsimplex is None:
            finalsimplex = [0, 1, 1]

        x0 = np.asarray(x0)
        n = len(x0)
        if step is None:
            step = abs(x0) + 1.2
        maxnfev = self.get_maxnfev(maxnfev, n)
        init = 0
        x0, fval, nfev, err = \
            _saoopt.neldermead(verbose, maxnfev, init, finalsimplex, tol, step,
                               xmin, xmax, x0, fcn)
        iquad = 1
        simp = 1.0e-2 * tol
        step = n * [0.4]
        self.par, self.fmin, tmpnfev, ifault = \
            _saoopt.minim(reflect, verbose, maxnfev - nfev, init, iquad, simp,
                          tol*10, step, xmin, xmax, x0, fcn)
        self.nfev = nfev + tmpnfev
        return self.nfev, self.fmin, self.par


class NelderMead5(NelderMead0):

    def __call__(self, fcn, x0, xmin, xmax, tol=1.0e-6,  maxnfev=None,
                 step=None, finalsimplex=1, verbose=0, reflect=True):
        init = 0
        iquad = 1
        simp = 1.0e-2 * tol
        x0 = np.asarray(x0)
        n = len(x0)
        if step is None:
            step = n * [0.4]
        maxnfev = self.get_maxnfev(maxnfev, n)
        par, fmin, nfev, ifault = \
            _saoopt.minim(reflect, verbose, maxnfev, init, iquad, simp, tol*10,
                          step, xmin, xmax, x0, fcn)
        return nfev, fmin, par


class NelderMead6(NelderMeadBase):

    # TODO: do we really need this internal class?
    class MyNelderMead6(MyNelderMead):

        def __call__(self, x, maxnfev, tol, step, finalsimplex, verbose):
            npar = len(x)
            simplex = SimplexNoStep(func=self.func, npop=npar + 1,
                                    xpar=x, xmin=self.xmin,
                                    xmax=self.xmax, step=None,
                                    seed=None, factor=None)
            return self.optimize(x, simplex, maxnfev, tol,
                                 finalsimplex, verbose)

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6,  maxnfev=None,
                 step=None, finalsimplex=1, verbose=0):
        my_nm_6 = NelderMead6.MyNelderMead6(fcn, xmin, xmax)
        if maxnfev is None:
            maxnfev = 512 * len(x)
        return my_nm_6(x, maxnfev, tol, step, finalsimplex, verbose)


class NelderMead7(NelderMeadBase):

    # TODO: do we really need this internal class?
    class MyNelderMead7(MyNelderMead):

        def __call__(self, x, maxnfev, tol, step, finalsimplex, verbose):
            npar = len(x)
            factor = 2
            simplex = SimplexRandom(func=self.func, npop=npar + 1, xpar=x,
                                    xmin=self.xmin, xmax=self.xmax,
                                    step=None, seed=None, factor=factor)
            return self.optimize(x, simplex, maxnfev, tol,
                                 finalsimplex, verbose)

    def __call__(self, fcn, x, xmin, xmax, tol=1.0e-6,  maxnfev=None,
                 step=None, finalsimplex=1, verbose=0):
        my_nm_7 = NelderMead7.MyNelderMead7(fcn, xmin, xmax)
        if maxnfev is None:
            maxnfev = 512 * len(x)
        return my_nm_7(x, maxnfev, tol, step, finalsimplex, verbose)


class nmNcores(MyNcores):

    def my_worker(self, opt, idval, out_q, err_q,
                  fcn, x, xmin, xmax, tol, maxnfev):
        try:
            vals = opt(fcn, x, xmin, xmax, tol, maxnfev)
        except Exception as e:
            err_q.put(e)
            return
        # output the result and task ID to output queue
        out_q.put((idval, vals))


class ncoresNelderMead:

    # TODO: using a list as an argument triggers pylint dangerous-default-value check
    def __init__(self, algo=[NelderMead0(), NelderMead1(), NelderMead2(),
                             NelderMead3(), NelderMead4(), NelderMead5()]):
        # NelderMead6(), NelderMead7()]):
        self.algo = algo

    def __call__(self, fcn, x, xmin, xmax, tol=EPSILON, maxnfev=None,
                 numcores=ncpus):

        num_algo = len(self.algo)
        nm_ncores = nmNcores()
        results = nm_ncores.calc(self.algo, numcores, fcn, x, xmin, xmax, tol, maxnfev)
        return self.unpack_results(num_algo, results)

    def unpack_results(self, num, results):
        nfev = results[0]
        fmin = results[1]
        par = results[2]
        solution_at = 0
        for ii in range(1, num):
            index = ii * 3
            nfev += results[index]
            # print(ii, 'unpack_results: f', par, '=', fmin, '@', nfev, 'nfevs')
            if results[index + 1] < fmin:
                fmin = results[index + 1]
                par = results[index + 2]
                solution_at = ii
        # print('unpack_results: solution_@ =', solution_at)
        return nfev, fmin, par


class ncoresNelderMeadRecursive(ncoresNelderMead):
    """ As noted in the paper, terminating the simplex is not a simple task:
    For any non-derivative method, the issue of termination is problematical as
    well as highly sensitive to problem scaling. Since gradient information is
    unavailable, it is provably impossible to verify closeness to optimality
    simply by sampling f at a finite number of points. Most implementations
    of direct search methods terminate based on two criteria intended to
    reflect the progress of the algorithm: either the function values at the
    vertices are close, or the simplex has become very small. """

    # TODO: using a list as an argument triggers pylint dangerous-default-value check
    #
    # algo is the same as used in ncoresNelderMead but leave as is in case
    # there is a need to have a different set of classes.
    #
    def __init__(self, algo=[NelderMead0(), NelderMead1(), NelderMead2(),
                             NelderMead3(), NelderMead4(), NelderMead5()]):
        ncoresNelderMead.__init__(self, algo)

    def __call__(self, fcn, x, xmin, xmax, tol=EPSILON, maxnfev=None,
                 numcores=ncpus):

        return self.calc(fcn, x, xmin, xmax, tol, maxnfev, numcores)

    def calc(self, fcn, x, xmin, xmax, tol=EPSILON, maxnfev=None,
             numcores=ncpus, fval=np.inf, nfev=0):

        num_algo = len(self.algo)
        nm_ncores = nmNcores()
        results = nm_ncores.calc(self.algo, numcores, fcn, x, xmin, xmax, tol, maxnfev)
        tmp_nfev, fmin, par = self.unpack_results(num_algo, results)
        nfev += tmp_nfev
        # print('ncoresNelderMead::calc f', par, ' = ', fmin, '@', nfev)
        if fmin < fval:
            return self.calc(fcn, par, xmin, xmax, tol, maxnfev, numcores, fmin, nfev)

        return nfev, fval, par
