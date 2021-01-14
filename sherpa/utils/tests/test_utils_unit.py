#
#  Copyright (C) 2016, 2018, 2019, 2020  Smithsonian Astrophysical Observatory
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
import pytest
import numpy

from sherpa.utils import _utils, is_binary_file, pad_bounding_box
from sherpa.utils.testing import requires_data
from numpy.testing import assert_almost_equal, assert_array_equal, \
    assert_array_almost_equal


def test_calc_ftest():
    """
    this test was created using sherpa 4.8.1 (2507414) as an oracle for the python 3 migration
    """
    assert_almost_equal(0.03452352914891555,
                        _utils.calc_ftest(11, 16.3, 10, 10.2))


def test_calc_ftest2():
    assert_array_almost_equal(2 * [0.03452352914891555],
                              _utils.calc_ftest(2 * [11], 2 * [16.3],
                                                2 * [10], 2 * [10.2]))


def test_calc_ftest_chi2_equal_0():
    with pytest.raises(TypeError) as excinfo:
        _utils.calc_ftest(11, 16.3, 10, 0.0)
    assert 'chisq_2[0] cannot be equal to 0:' in str(excinfo.value)


def test_calc_ftest_dof2_equal_0():
    with pytest.raises(TypeError) as excinfo:
        _utils.calc_ftest(11, 16.3, 0, 10.0)
    assert 'dof_2[0] cannot be equal to 0:' in str(excinfo.value)


def test_calc_ftest_dof1_equal_dof2():
    with pytest.raises(TypeError) as excinfo:
        _utils.calc_ftest(11, 16.3, 11, 10.0)
    assert 'dof_1[0] cannot be equal to dof_2[0]:' in str(excinfo.value)


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


def test_rebin_partial_overlap():
    """Does rebin handle the case when the output grid is > input?

    This only tests the easy case (when the model evaluates to 0
    for the missing range), but is added whilst investigating #722
    (the handling of regrid for integrated 1D models). It is intended
    to show the problem in #722 comes from the regrid code, not the
    rebin routine.
    """

    # Truth is Box1D evaluated with xlow=3 xhi=5 ampl=4 evaluated
    # on the integrated bin defined by arange(1, 7, 0.5).
    #
    # However, the model is evaluated on the grid arange(2.1, 6, 0.5)
    # and then rebinned onto the requested grid. So, the edge bins
    # are partially filled.
    #
    xdataspace = numpy.arange(1, 7, 0.5)
    xdlo, xdhi = xdataspace[:-1], xdataspace[1:]

    xgrid = numpy.arange(2.1, 6, 0.5)
    xglo, xghi = xgrid[:-1], xgrid[1:]

    # Manually calculated
    ygrid = numpy.asarray([0, 0.4, 2, 2, 2, 1.6, 0])

    yans = _utils.rebin(ygrid, xglo, xghi, xdlo, xdhi)

    # Manually calculated
    yexp = numpy.asarray([0, 0, 0, 0.32, 1.68, 2, 2, 1.68, 0.32, 0, 0])
    assert yans == pytest.approx(yexp)


def test_sao_arange():
    """
    used sherpa (python2) as oracle

    py3-todo: this function has no python docs
    """
    assert_array_equal([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.], _utils.sao_arange(0, 10, 1))


# py3-todo: skipping test for sum_intervals, not enough info to make a quick test.

@requires_data
def test_is_binary_file(make_data_path):
    ascii = make_data_path("gauss2d.dat")
    pha = make_data_path("3c273.pi")

    assert is_binary_file(pha)
    assert not is_binary_file(ascii)


def test_pad_bounding_box_fail():
    """The mask size must be >= kernel."""

    kernel = numpy.arange(12)
    mask = numpy.ones(10)
    with pytest.raises(TypeError) as excinfo:

        pad_bounding_box(kernel, mask)

    emsg = 'kernel size: 12 is > than mask size: 10'
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("mask, expected",
                         [([1, 1, 1, 1, 1], [1, 2, 3, 4, 5]),
                          ([1, 1, 1, 1, 0], [1, 2, 3, 4, 0]),
                          ([0, 1, 1, 1, 1], [0, 1, 2, 3, 4]),
                          ([0, 1, 0, 1, 1], [0, 1, 0, 2, 3]),
                          ([0, 0, 0, 0, 1], [0, 0, 0, 0, 1]),
                          ([0, 0, 0, 0, 0], [0, 0, 0, 0, 0]),
                          ([0, 0, 0, 1, 1, 1, 1, 1, 0, 0],
                           [0, 0, 0, 1, 2, 3, 4, 5, 0, 0]),
                          ([1, 0, 0, 1, 0, 1, 1, 0, 0, 1],
                           [1, 0, 0, 2, 0, 3, 4, 0, 0, 5])
                          ])
def test_pad_bounding_box(mask, expected):
    """Basic tests"""

    kernel = numpy.asarray([1, 2, 3, 4, 5])

    # The tests use equality rather than approximate equality
    # since the values are just being copied around by the
    # code, not manipulated.
    #
    exp = numpy.asarray(expected).astype(numpy.float64)

    ans = pad_bounding_box(kernel, mask)
    assert_array_equal(ans, exp)


def test_pad_bounding_box_mask_too_large():
    """What happens when the mask has more valid elements
    than the kernel? At present the code will read out-of-bounds
    elements rather than use 0 elements, which means that the
    test below should fail (in all but the most unlikely case
    of reading into memory which happens to be zeros).

    Perhaps the code should catch this condition and error out
    instead?
    """

    kernel = numpy.arange(5)
    mask = numpy.ones(10)

    ans = pad_bounding_box(kernel, mask)

    # Assume that the "extra" mask elements get mapped to 0 (or
    # ignored).
    exp = numpy.asarray([0, 1, 2, 3, 4, 0, 0, 0, 0, 0]).astype(numpy.float64)
    assert_array_equal(ans, exp)
