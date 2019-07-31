from __future__ import print_function
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


from sherpa.utils import _ncpus
from sherpa.optmethods.ncoresnm import ncoresNelderMead
from sherpa.optmethods.ncoresde import ncoresDifEvo
import sherpa.optmethods.opt as tstopt

import pytest

NUMPAR = 10
NCORES_NM = ncoresNelderMead()
NCORES_DE = ncoresDifEvo()


def print_result(name, f, x, nfev):
    print('%s(%s) = %g in %d nfev' % (name, x, f, nfev))


# Mapping from SherpaTestCase.assertEqualWithinTol to
#              pytest.approx
# and noting that the tolerance is an absolute tolerance,
# not relative in SherpaTestCase.
#
def tst_opt(opt, fcn, x0, xmin, xmax, fmin, tol=1e-2):
    def func(arg):
        return fcn(arg)[0]
    nfev, fval, par = opt(func, x0, xmin, xmax)
    assert fval == pytest.approx(fmin, abs=tol)


def tst_func(opt, func, x0, xmin, xmax, fmin, tol=1e-2):
    nfev, fval, par = opt(func, x0, xmin, xmax)
    assert fval == pytest.approx(fmin, abs=tol)


def test_Ackley():
    xmin = NUMPAR * [-32.768]
    xmax = NUMPAR * [32.768]
    x0 = NUMPAR * [12.3]
    if _ncpus != 1:
        tst_func(NCORES_NM, tstopt.Ackley, x0, xmin, xmax, 0.0)

    tst_func(NCORES_DE, tstopt.Ackley, x0, xmin, xmax, 0.0)


def test_Bohachevsky1():
    xmin = [-100, -100]
    xmax = [100, 100]
    x0 = [-12, 10]

    fmin = 0.0

    tst_func(NCORES_NM, tstopt.Bohachevsky1, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Bohachevsky1, x0, xmin, xmax, fmin)


def test_Bohachevsky2():
    xmin = [-100, -100]
    xmax = [100, 100]
    x0 = [12, 10]

    fmin = 0.0

    tst_func(NCORES_NM, tstopt.Bohachevsky2, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Bohachevsky2, x0, xmin, xmax, fmin)


def test_Bohachevsky3():
    xmin = [-100, -100]
    xmax = [100, 100]
    x0 = [-61.2, 51.0]

    fmin = 0.0
    tst_func(NCORES_NM, tstopt.Bohachevsky3, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Bohachevsky3, x0, xmin, xmax, fmin)


def test_Booth():
    xmin = [-10, -10]
    xmax = [10, 10]
    x0 = [-6.2, 5.0]

    fmin = 0.0
    tst_func(NCORES_NM, tstopt.Booth, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Booth, x0, xmin, xmax, fmin)


def test_BoxBetts():
    xmin = [0.9, 9.0, 0.9]
    xmax = [1.2, 11.2, 1.2]
    x0 = [(xmin[0] + xmax[0]) * 0.5, (xmin[1] + xmax[1]) * 0.5,
          (xmin[2] + xmax[2]) * 0.5]

    fmin = 0.0
    tst_func(NCORES_NM, tstopt.BoxBetts, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.BoxBetts, x0, xmin, xmax, fmin)


def test_Branin():
    xmin = [-5, 0]
    xmax = [10, 15]
    x0 = [-3.2, 5.0]

    fmin = 0.397887
    tst_func(NCORES_NM, tstopt.Branin, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Branin, x0, xmin, xmax, fmin)


def test_DixonPrixe():
    xmin = NUMPAR * [-10]
    xmax = NUMPAR * [10]
    x0 = NUMPAR * [-1.]

    fmin = 0.0

    tst_func(NCORES_DE, tstopt.DixonPrice, x0, xmin, xmax, fmin)


def test_Easom():
    xmin = [-100, -100]
    xmax = [100, 100]
    x0 = [25.0, 25.0]

    fmin = -1.0

    tst_func(NCORES_DE, tstopt.Easom, x0, xmin, xmax, fmin)


def test_FreudensteinRoth():
    xmin = NUMPAR * [-1000, -1000]
    xmax = NUMPAR * [1000, 1000]
    x0 = NUMPAR * [0.5, -2]

    fmin = 0.0

    tst_func(NCORES_NM, tstopt.FreudensteinRoth, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.FreudensteinRoth, x0, xmin, xmax, fmin)


def test_GoldsteinPrice():
    xmin = [-2, -2]
    xmax = [2, 2]
    x0 = [-1, 1]

    fmin = 3.0

    tst_func(NCORES_NM, tstopt.GoldsteinPrice, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.GoldsteinPrice, x0, xmin, xmax, fmin)


def test_Hump():
    xmin = [-5, -5]
    xmax = [5, 5]
    x0 = [-3.2, 5.0]

    fmin = 0.0

    tst_func(NCORES_NM, tstopt.Hump, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Hump, x0, xmin, xmax, fmin)


def test_Matyas():
    xmin = [-10, -10]
    xmax = [10, 10]
    x0 = [-3.2, 5.0]

    fmin = 0.0

    tst_func(NCORES_NM, tstopt.Matyas, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Matyas, x0, xmin, xmax, fmin)


# TODO: mark this as only-run on multi-core systems
def test_McCormick():
    xmin = [-1.5, -3.0]
    xmax = [4.0, 4.0]
    x0 = [0.0, 0.0]

    fmin = -1.91

    if _ncpus != 1:
        tst_func(NCORES_NM, tstopt.McCormick, x0, xmin, xmax, fmin)

    # tst_func(NCORES_DE, tstopt.McCormick, x0, xmin, xmax, fmin)


def test_Paviani():
    xmin = NUMPAR * [2.001]
    xmax = NUMPAR * [9.999]
    x0 = NUMPAR * [5.0]

    # ans = -45.7

    # From
    # http://infinity77.net/global_optimization/test_functions_nd_P.html#go_benchmark.Paviani
    # which gives -45.7784684040686
    #
    fmin = -45.778

    tst_func(NCORES_NM, tstopt.Paviani, x0, xmin, xmax, fmin)
    tst_func(NCORES_DE, tstopt.Paviani, x0, xmin, xmax, fmin)


def test_Rastrigin():
    xmin = NUMPAR * [-5.12]
    xmax = NUMPAR * [5.12]
    x0 = NUMPAR * [-2.0]
    if _ncpus != 1:
        tst_func(NCORES_NM, tstopt.Rastrigin, x0, xmin, xmax, 0.0)

    tst_func(NCORES_DE, tstopt.Rastrigin, x0, xmin, xmax, 0.0)


def test_Shubert():
    xmin = [-10, -10]
    xmax = [10, 10]
    x0 = [-2.0, 5.0]

    fmin = -186.7309

    if _ncpus != 1:
        tst_func(NCORES_NM, tstopt.Shubert, x0, xmin, xmax, fmin)

    tst_func(NCORES_DE, tstopt.Shubert, x0, xmin, xmax, fmin)


#
# comment out a few time consuming tests, anyone
# modifying the ncores opt should uncomment the tests
#
# def test_Levy():
#     xmin = NUMPAR * [-10]
#     xmax = NUMPAR * [10]
#     x0 = NUMPAR * [-5.]
#     tst_func(NCORES_NM, tstopt.Levy, x0, xmin, xmax, 0.0)
#     tst_func(NCORES_DE, tstopt.Levy, x0, xmin, xmax, 0.0)
# def test_Sphere():
#     xmin = NUMPAR * [-5.12]
#     xmax = NUMPAR * [5.12]
#     x0 = NUMPAR * [-2.0]
#     tst_func(NCORES_NM, tstopt.Sphere, x0, xmin, xmax, 0.0)
#     tst_func(NCORES_DE, tstopt.Sphere, x0, xmin, xmax, 0.0)

# def test_SumSquares():
#     xmin = NUMPAR * [-10]
#     xmax = NUMPAR * [10]
#     x0 = NUMPAR * [-2.0]
#     tst_func(NCORES_NM, tstopt.SumSquares, x0, xmin, xmax, 0.0)
#     tst_func(NCORES_DE, tstopt.SumSquares, x0, xmin, xmax, 0.0)

# def test_Zakharov():
#     xmin = NUMPAR * [-5, -5]
#     xmax = NUMPAR * [10, 10]
#     x0   = NUMPAR * [0.5, -2]
#     tst_func(NCORES_NM, tstopt.Zakharov, x0, xmin, xmax, 0.0)
#     tst_func(NCORES_DE, tstopt.Zakharov, x0, xmin, xmax, 0.0)
