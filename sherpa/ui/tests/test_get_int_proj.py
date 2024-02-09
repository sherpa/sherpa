#
#  Copyright (C) 2019 - 2021, 2024
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

import pytest

import numpy as np

from sherpa import ui
from sherpa.utils.err import ConfidenceErr


def cmp_args(result):
    assert result.min == -5
    assert result.max == 0
    assert result.nloop == 21
    assert result.delv is None
    assert result.fac == 1
    assert result.log is False


def cmp_x(result):

    expected = np.array([-5., -4.75, -4.5, -4.25, -4., -3.75, -3.5, -3.25, -3.,
                         -2.75, -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75,
                         -0.5, -0.25, 0.])
    assert result.x == pytest.approx(expected)


def cmp_args_log(result):
    assert result.min == 0.1
    assert result.max == 0.9
    assert result.nloop == 10
    assert result.delv is None
    assert result.fac == 1
    assert result.log is True


def cmp_x_log(result):

    expected = np.logspace(-1, np.log10(0.9), 10)
    assert result.x == pytest.approx(expected)


# Also used in the log tests, hence a global.
#
X = [-13, -5, -3, 2, 7, 12]

@pytest.fixture
def setUp(clean_ui, hide_logging):

    y = [102.3, 16.7, -0.6, -6.7, -9.9, 33.2]
    dy = np.ones(6) * 5
    ui.load_arrays(1, X, y, dy)
    ui.set_source(ui.polynom1d.poly)
    ui.thaw(poly.c1, poly.c2)
    ui.fit()
    # Just check we are fitting 6 data points with 3 parameters
    assert ui.get_fit_results().dof == 3


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
@pytest.mark.parametrize("val", [-1, 0, 1])
def test_get_int_xxx_positive_nloop(setUp, method, val):
    """nloop must be > 1"""

    with pytest.raises(ConfidenceErr,
                       match="Nloop parameter must be > 1"):
        method(poly.c1, recalc=True, min=-5, max=0, nloop=val)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_limits_sensible(setUp, method):
    """max must be > min"""

    with pytest.raises(ConfidenceErr,
                       match="Bad parameter limits"):
        method(poly.c1, recalc=True, min=0, max=0, nloop=10)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_limits_scalar(setUp, method):
    """max and min must be scalars"""

    with pytest.raises(ConfidenceErr,
                       match="Parameter limits must be scalars"):
        method(poly.c1, recalc=True, min=(-5, 2), max=(0, 4), nloop=10)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
@pytest.mark.parametrize("minval,maxval", [(0, 1), (-2, -1), (-2, 0)])
def test_get_int_xxx_limits_positive_when_log(setUp, method, minval, maxval):
    """max and min must be positive when log scale"""

    with pytest.raises(ConfidenceErr,
                       match="Log scale must be on positive boundaries"):
        method(poly.c1, recalc=True, min=minval, max=maxval, nloop=10,
               log=True)


# FAIL: this currently fails with ZeroDivisionError.
@pytest.mark.xfail
@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_delv_positive(setUp, method):
    """delv must be positive"""

    with pytest.raises(ConfidenceErr,
                       match="delv parameter must be > 0"):
        method(poly.c1, recalc=True, min=-5, max=-1, delv=-0.1)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_not_leastsq(setUp, method):
    """can not use the least-square statistic"""

    ui.set_stat("leastsq")
    with pytest.raises(ConfidenceErr,
                       match="leastsq is inappropriate for confidence limit estimation"):
        method(poly.c1, recalc=True, min=-5, max=0, nloop=10)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_not_thawed(setUp, method):
    """can not use on a frozen parameter"""

    with pytest.raises(ConfidenceErr,
                       match="Frozen parameter poly.c4 cannot be used for interval "):
        method(poly.c4, recalc=True, min=1, max=3, nloop=10)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_get_int_xxx_no_recalc(setUp, method, check_str):
    """what is the output like when recalc=False (the default)?"""

    res = method(poly.c1, min=-2, max=-1, nloop=5, recalc=False)
    assert res.x is None
    assert res.y is None
    assert res.min is None
    assert res.max is None
    assert res.nloop == 20
    assert res.delv is None
    assert res.fac == 1
    assert not res.log

    check_str(str(res),
              ["x     = None",
               "y     = None",
               "min   = None",
               "max   = None",
               "nloop = 20",
               "delv  = None",
               "fac   = 1",
               "log   = False"
               ])


@pytest.mark.parametrize("method,expected",
                         [(ui.get_int_unc, "[ 6.7979,17.4686,36.1392]"),
                          (ui.get_int_proj, "[ 6.7782,17.3734,35.9119]")
                          ])
def test_get_int_xxx_show(setUp, method, expected, check_str):
    """Check string output"""

    res = method(poly.c1, min=-2, max=-1, nloop=3, recalc=True)
    check_str(str(res),
              ["x     = [-2. ,-1.5,-1. ]",
               f"y     = {expected}  # doctest: +FLOAT_CMP",
               "min   = -2",
               "max   = -1",
               "nloop = 3",
               "delv  = None",
               "fac   = 1",
               "log   = False"
               ])


def test_get_int_proj(setUp):
    # check we get the poly.c1 results with recalc=True
    ui.int_proj(poly.c0)
    result = ui.get_int_proj(poly.c1, recalc=True, min=-5, max=0, nloop=21)
    cmp_args(result)
    cmp_x(result)

    expected = np.array([110.0185, 90.493, 72.9534, 57.3995, 43.8316,
                         32.2494, 22.6531, 15.0427, 9.4181, 5.7794, 4.1265,
                         4.4594, 6.7782, 11.0829, 17.3734, 25.6497, 35.9119,
                         48.1599, 62.3938, 78.6135, 96.8191])
    assert result.y == pytest.approx(expected, abs=1.0e-4)


def test_get_int_unc(setUp):
    # check we get the poly.c1 results with recalc=True
    ui.int_unc(poly.c0)
    result = ui.get_int_unc(poly.c1, recalc=True, min=-5, max=0, nloop=21)
    cmp_args(result)
    cmp_x(result)

    expected = np.array([110.774, 91.1094, 73.4447, 57.78, 44.1153, 32.4507,
                         22.786, 15.1213, 9.4566, 5.7919, 4.1273, 4.4626,
                         6.7979, 11.1332, 17.4686, 25.8039, 36.1392, 48.4745,
                         62.8099, 79.1452, 97.4805])
    assert result.y == pytest.approx(expected, abs=1.0e-4)


@pytest.mark.xfail
def test_get_int_unc_log(setUp):
    """Use log scaling to check. See issue #1561

    Note that c2 is the only positive parameter we can test with.
    """

    result = ui.get_int_unc(poly.c2, recalc=True, log=True,
                            min=0.1, max=0.9, nloop=10)

    cmp_args_log(result)
    cmp_x_log(result)

    # We can work out the expected vaues as we know all the parameter
    # values.
    #
    expected = []
    for x in result.x:
        poly.c2 = x
        expected.append(ui.calc_stat(1))

    assert result.y == pytest.approx(expected, abs=1.0e-4)


@pytest.mark.xfail
def test_get_int_proj_log(setUp):
    """Use log scaling to check. See issue #1561

    Note that c2 is the only positive parameter we can test with.
    """

    result = ui.get_int_proj(poly.c2, recalc=True, log=True,
                            min=0.1, max=0.9, nloop=10)

    cmp_args_log(result)
    cmp_x_log(result)

    # Unlike the int_unc command it's harder to a-priori know the
    # expected values, so treat as a regression test. We can compare
    # to the linear case to check that they make sense.
    #
    expected = np.asarray([150.37621233, 129.76054302, 105.71703208,
                           78.72851931, 50.31198137, 23.87142974,
                           6.14353848, 9.62456249, 56.61597458,
                           185.932959])
    assert result.y == pytest.approx(expected, abs=1.0e-4)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
@pytest.mark.parametrize("fac", [1, 3])
def test_int_xxx_fac(setUp, method, fac):
    """What happens if we use fac to calculate the limits"""

    c = 0.47827334
    dc = 0.0312677
    expected = [c - fac * dc, c, c + fac * dc]

    result = method(poly.c2, recalc=True, fac=fac, nloop=3)
    assert result.x == pytest.approx(expected)


@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_int_xxx_delv(setUp, method):
    """What happens if we set delv"""

    result = method(poly.c2, recalc=True, min=0.1, max=0.9, delv=0.3)
    assert result.x == pytest.approx([0.1, 0.4, 0.7, 1])


@pytest.mark.xfail
@pytest.mark.parametrize("method", [ui.get_int_unc, ui.get_int_proj])
def test_int_xxx_delv_log(setUp, method):
    """What happens if we set delv and log?"""

    result = method(poly.c2, recalc=True, log=True, min=0.1, max=0.9,
                    delv=0.3)

    expected = 10**np.asarray([-1, -0.7, -0.4, -0.1])
    assert result.x == pytest.approx(expected)
