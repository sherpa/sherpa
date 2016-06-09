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
    arrays = ([0.12345, 0.54321], [0.12346, 0.54320])
    assert_array_equal(expected, _utils.gsl_fcmp(*arrays, tolerance))


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
    arrays = ([0.12345, 0.54321], [0.12346, 0.54320])
    assert_array_equal(expected, _utils.sao_fcmp(*arrays, tolerance))


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


def test_igam():
    """
    used scipy.special.gamminc as oracle
    """
    assert_almost_equal(0.4421745996, _utils.igam(2, 1.5))
