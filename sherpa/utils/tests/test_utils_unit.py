#
#  copyright (c) 2016  smithsonian astrophysical observatory
#
#
#  this program is free software; you can redistribute it and/or modify
#  it under the terms of the gnu general public license as published by
#  the free software foundation; either version 3 of the license, or
#  (at your option) any later version.
#
#  this program is distributed in the hope that it will be useful,
#  but without any warranty; without even the implied warranty of
#  merchantability or fitness for a particular purpose.  see the
#  gnu general public license for more details.
#
#  you should have received a copy of the gnu general public license along
#  with this program; if not, write to the free software foundation, inc.,
#  51 franklin street, fifth floor, boston, ma 02110-1301 usa.
#
import pytest
import numpy

from sherpa.utils import _utils
from numpy.testing import assert_almost_equal, assert_array_equal


def test_calc_ftest():
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    assert_almost_equal(0.475500468, _utils.calc_ftest(10, 1.0, 20, 2.0))


def test_calc_mlr():
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    assert_almost_equal(0.415880186, _utils.calc_mlr(5, 5.0))


def test_calc_erf():
    """
    used math.erf as oracle

    py3-todo: should we use the built-in python function instead?
    """
    assert_almost_equal(0.520499877, _utils.erf(0.5))


def test_calc_gamma():
    """
    used math.gamma as oracle

    py3-todo: should we use the built-in python function instead?
    """
    assert_almost_equal(1.7724538509, _utils.gamma(0.5))


@pytest.mark.parametrize("tolerance, expected", [
    (0.01, [0, 0]),
    (0.00001, [-1, 0]),
    (0.000001, [-1, 1]),
])
def test_calc_gsl_fcmp(tolerance, expected):
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    assert_array_equal(expected, _utils.gsl_fcmp([0.12345, 0.54321], [0.12346, 0.54320], tolerance))


@pytest.mark.parametrize("tolerance, expected", [
    (0.01, [0, 0]),
    (0.00001, [-1, 0]),
    (0.000001, [-1, 1]),
])
def test_calc_sao_fcmp(tolerance, expected):
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration

    py3-todo: what is the difference between sao_fcmp and gsl_fcmp?
    """
    assert_array_equal(expected, _utils.sao_fcmp([0.12345, 0.54321], [0.12346, 0.54320], tolerance))


def test_calc_hist1d():
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    x, xlo, xhi = numpy.arange(0, 21)/2, range(0, 20, 1), range(1, 21, 1)
    expected = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert_array_equal(expected, _utils.hist1d(x, xlo, xhi))


def test_calc_hist2d():
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    x = numpy.arange(10)//2
    xx, yy = numpy.meshgrid(x, x)
    xin = xx.ravel()
    yin = yy.ravel()

    expected = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 0, 0,
       0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 4, 0, 4, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,
       0, 4, 0, 4, 0, 4, 0, 4]

    assert_array_equal(expected, _utils.hist2d(xin, yin, x, x))


@pytest.mark.parametrize("a, x, expected", [
    (2, 1.5, 0.4421745996),
    ([1, 2, 3], [1.0, 1.5, 2.0], [0.6321205588, 0.4421745996, 0.323323583]),
])
def test_igam(a, x, expected):
    """
    used scipy.special.gamminc as oracle
    """
    assert_almost_equal(expected, _utils.igam(a, x))


@pytest.mark.parametrize("a, x, expected", [
    (2, 1.5, 0.557825400),
    ([1, 2, 3], [1.0, 1.5, 2.0], [0.36787944, 0.5578254, 0.67667642]),
])
def test_igamc(a, x, expected):
    """
    used scipy.special.gammincc as oracle
    """
    assert_almost_equal(expected, _utils.igamc(a, x))


@pytest.mark.parametrize("a, b, x, expected", [
    (2, 3.5, 0.5, 0.756932043),
    ([1, 2, 3], [1.0, 1.5, 2.0], [0.3, 0.5, 0.75], [0.3, 0.38128157, 0.73828125]),
])
def test_incbet(a, b, x, expected):
    """
    used scipy.special.betainc as oracle
    """
    assert_almost_equal(expected, _utils.incbet(a, b, x))


@pytest.mark.parametrize("z, expected", [
    (30, 71.25703896),
    ([1.0, 1.5, 25.7, 132.3], [0, -0.1207822376, 57.0337539076, 512.4720673873])
])
def test_lgam(z, expected):
    """
    used scipy.special.gammaln as oracle
    """
    assert_almost_equal(expected, _utils.lgam(z))


@pytest.mark.parametrize("x, expected", [
    (0.1, -1.2815515655446004),
    ([0.1, 0.4, 0.95], [-1.28155157, -0.2533471, 1.64485363])
])
def test_ndtri(x, expected):
    """
    used scipy.stats.norm.ppf as oracle.

    py3-todo: this function is missing docs in python
    """
    assert_almost_equal(expected, _utils.ndtri(x))


def test_neville():
    """
    used sherpa (python2) as oracle.
    """
    x = [1.2, 3.4, 4.5, 5.2]
    y = [12.2, 14.4, 16.8, 15.5]
    xgrid = numpy.linspace(2, 5, 5)
    expected = [10.7775023, 12.24227718, 14.7319838, 16.60004786, 16.19989505]
    assert_almost_equal(expected, _utils.neville(xgrid, x, y))


def test_rebin():
    """
    used sherpa (python2) as oracle
    """
    y0 = [1, 3, 4, 10, 0]
    x0lo = [0, 2, 4, 6,  8]
    x0hi = [2, 4, 6, 8, 10]
    x1lo = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    x1hi = [1,  2,  3,  4,  5,  6,  7,  8,  9, 10]

    expected = [0.5,  0.5,  1.5,  1.5,  2.,  2.,  5.,  5.,  0.,  0.]

    assert_array_equal(expected, _utils.rebin(y0, x0lo, x0hi, x1lo, x1hi))


def test_sao_arange():
    """
    used sherpa (python2) as oracle

    py3-todo: this function has no python docs
    """
    assert_array_equal([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.], _utils.sao_arange(0, 10, 1))


# py3-todo: skipping test for sum_intervals, not enough info to make a quick test.
