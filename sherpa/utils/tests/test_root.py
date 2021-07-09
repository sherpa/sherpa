#
#  Copyright (C) 2007, 2016, 2018, 2020, 2021  Smithsonian Astrophysical Observatory
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

import sys

import numpy

import pytest

from sherpa.utils import demuller, bisection, new_muller, apache_muller, \
    zeroin


def sqr(x, *args):
    return x * x


def prob1(x, *args):
    return numpy.cos(x) - pow(x, 3.0)


def prob2(x, *args):
    return pow(x, 4.0) - x - 1.0


def prob3(x, *args):
    return (x + 3.0) * (x - 1.0) * (x - 1.0)


def prob4(x, *args):
    return numpy.exp(x) + x - 2.0


def prob5(x, *args):
    return x + numpy.exp(x)


def prob6(x, *args):
    return numpy.sqrt(x) + x * x - 5.0


def prob7(x, *args):
    return x + numpy.sin(x) - 4.0


def prob8(x, *args):
    return x * x - numpy.exp(x)


def prob9(x, *args):
    return numpy.sin(x) + numpy.cos(x)


def prob10(x, *args):
    return x * x - 2.0


def prob11(x, *args):
    x2 = sqr(x)
    x4 = sqr(x2)
    x8 = sqr(x4)
    return x8 * x


def prob12(x, *args):
    if x == 0.0:
        return 0.0
    return x - numpy.exp(-1.0 / (x * x))


def prob13(x, *args):
    x2 = sqr(x - 1.0)
    x4 = x2 * x2
    return x * x4 * (x - 1.0)


def prob14(x, *args):
    return (x - 5.0) * (x - 2.0)


def prob15(x, *args):
    return (x - 50.0) * (x + 20.0)


def prob16(x, *args):
    return (sqr(x) - 2.0) * x - 5.0


def prob17(x, *args):
    return sqr(x - 1.0e-5) - 3.0


def prob18(x, *args):
    return numpy.cos(x) - x


def prob19(x, *args):
    return numpy.sin(x) - x


def prob20(x, *args):
    return x * x * x - 2.0 * x - 5.0


def prob21(x, *args):
    return numpy.sin(x) - 0.5 * x


def prob22(x, *args):
    return 2 * x - numpy.exp(-x)


def prob23(x, *args):
    return x * numpy.exp(-x)


def prob24(x, *args):
    return numpy.exp(x) - 1. / (100 * x * x)


def prob25(x, *args):
    tmp = x - 1
    return (x + 3) * tmp * tmp


def prob26(x, *args):
    return numpy.exp(x) - 2 - 1 / (10 * x)**2 - 2 / (100 * x)**3


def prob27(x, *args):
    return x**3


def prob28(x, *args):
    return numpy.cos(x) - x


# root is 1.83 the convergence rate for Muller's method
def prob29(x, *args):
    return pow(x, 3.0) - x * x - x - 1


# muller gets it wrong but muller_bound got it right
def prob30(x, *args):
    return numpy.exp(x) - 1


def prob31(x, *args):
    return x * numpy.exp(x) - 1


def prob32(x, *args):
    return 11.0 * pow(x, 11.0) - 1.0


def prob33(x, *args):
    return numpy.exp((1.0 * x + 7) * x - 30.0) - 1.0


def prob34(x, *args):
    return 1.0 / x - numpy.sin(x) + 1.0


def prob35(x, *args):
    return (x*x - 2.0) * x - 5.0


def prob36(x, *args):
    return 1.0 / x - 1.0


def prob37(x, *args):
    return numpy.log(x)


def prob38(x, *args):
    return x - numpy.exp(numpy.sin(x))


def repeller(x, *args):
    return 20.0 * x / (100.0 * x * x + 1.0)


def pinhead(x, *args):
    epsilon = 0.00001
    if epsilon == 0.0:
        return (16.0 - x**4) / (16.0 * x**4)

    return (16.0 - x**4) / (16.0 * x**4 + epsilon)


#    Sample results:
#
#    E        M      X
#    -----  ---  ----------
#    0.100    5    5.554589
#    0.200    5    6.246908
#    0.300    5    7.134960
#    0.400    5    8.313903
#    0.500    5    9.950063
#    0.600    5   12.356653
#    0.700    5   16.167990
#    0.800    5   22.656579
#    0.900    5   33.344447
#    0.990    5   45.361023
#    0.990    1   24.725822
#    0.990   33   89.722155
#    0.750   70  110.302
#    0.990    2   32.361007


def kepler(x, *args):
    e = 0.8
    m = 5.0
    pi = 3.141592653589793
    return pi * (x - m) / 180.0 - e * numpy.sin(pi * x / 180.0)


def wallis(x, *args):
    return x**3 - 2*x - 5


# NOTE: unused
def thinpole(x, *args):
    pi = 3.141592653589793
    return 3.0 * x**2 + 1.0 + (numpy.log((x - pi)**2)) / pi**4


def muller_convergence_rate(x):
    return x - pow(x, 3.0) / 3.0


def demuller2(fcn, xa, xb, fa=None, fb=None, args=(), maxfev=32,
              tol=1.0e-6):
    return demuller(fcn, xa, xb, (xa+xb)/2.0, args=args, maxfev=maxfev,
                    tol=tol)


@pytest.mark.parametrize('method', [bisection, demuller2, new_muller, apache_muller, zeroin])
@pytest.mark.parametrize('problem,a,b',
                         [(prob1, 0.0, 1.0),
                          (prob2, -5.0, 1.2),
                          (prob3, -10.0, 10.0),
                          (prob4, -1.0, 1.0),
                          (prob5, -1.0, 1.0),
                          (prob6, 0.0, 3.0),
                          (prob7, 0.0, 5.0),
                          (prob8, -5.0, 0.0),
                          (prob9, 0.0, 3.14),
                          (prob10, 0.0, 2.0),
                          (prob11, -1.0, 1.5),
                          (prob12, -1.0, 1.0),
                          (prob13, -0.5, 0.99),
                          (prob14, 1.0, 3.0),
                          (prob15, 15.0, 75.0),
                          (prob16, 2.0, 3.0),
                          (prob17, -1.0, 3.0),
                          (prob18, -1.0, 3.0),
                          (prob19, -1.0, 3.0),
                          (prob20, 2.0, 3.0),
                          (prob21, -1000.0, 1000.0),
                          (prob22, -10.0, 100.0),
                          (prob23, -10.0, 100.0),
                          (prob24, 1.0e-4, 20.0),
                          (prob25, -1000.0, 1000.0),
                          (prob26, 1.0e-4, 20.0),
                          (prob27, -1000.0, 1000.0),
                          (prob28, -10.0, 10.0),
                          (prob29, 0, 10),
                          (prob30, -50, 100.0),
                          (prob31, -1.0, 1.0),
                          (prob32, 0.1, 0.9),
                          (prob33, 2.8, 3.1),
                          (prob34, -1.3, -0.5),
                          (prob35, 2.0, 3.0),
                          (prob36, 0.5, 1.5),
                          (prob37, 0.5, 5.0),
                          (prob38, 1.0, 4.0),
                          (repeller, -10.0, 10.0),
                          (pinhead, 1.0e-5, 10.0),
                          (kepler, -175.0, 185.0),
                          (wallis, 2.0, 3.0),
                          (muller_convergence_rate, 1.7, 10.0)
                         ])
def test_case(method, problem, a, b):
    """Test a given problem"""

    verbose = False
    tol = 1.0e-6
    if verbose:
        sys.stdout.write('\n')

    # demuller2 cannot solve prob30 & pinhead so may as well skip them
    if method.__name__ == "demuller2" and problem.__name__ in ['prob30', 'pinhead']:
        pytest.skip('Known failure case')

    result = method(problem, a, b, tol=tol)
    [root, froot] = result[0]
    if verbose:
        msg = '%s:\t%s(%e) = %e in %d nfevs\n' % (method.__name__,
                                                  problem.__name__,
                                                  root, froot, result[-1])
        sys.stdout.write(msg)

    assert abs(froot) <= tol
