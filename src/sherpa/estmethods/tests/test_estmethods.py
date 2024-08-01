#
#  Copyright (C) 2007, 2018, 2021, 2023
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

import re

import numpy

import pytest

from sherpa.estmethods import Confidence, Covariance, Projection


# Test data arrays -- together this makes a line best fit with a
# 1D Gaussian function.
x = numpy.array(
    [758,  759,  760,  761,  762,  763,  764,  765,  766,  767,  768,  769,
     770,  771,  772,  773,  774,  775,  776,  777,  778,  779,  780,  781,
     782,  783,  784,  785,  786,  787,  788,  789,  790,  791,  792,  793,
     794,  795,  796,  797,  798,  799,  800,  801,  802,  803,  804,  805,
     806,  807,  808,  809,  810,  811,  812,  813,  814,  815,  816,  817,
     818,  819,  820,  821,  822,  823,  824,  825,  826,  827,  828,  829,
     830,  831,  832,  833,  834,  835,  836,  837,  838,  839,  840,  841,
     842,  843,  844,  845,  846,  847,  848,  849,  850,  851,  852,  853,
     854,  855,  856,  857,  858,  859,  860,  861,  862,  863,  864,  865,
     866,  867,  868,  869,  870,  871,  872,  873,  874,  875,  876,  877,
     878,  879,  880,  881])

y = numpy.array(
    [1,    0,    1,    0,    1,    2,    0,    3,    0,    1,    0,    5,
     0,    5,    3,    5,    6,    5,   11,   14,   11,   13,   12,   21,
     15,   24,   20,   29,   32,   43,   47,   49,   50,   64,   60,   72,
     61,   73,   83,   98,   99,  100,   94,   92,  121,  107,  126,  107,
     112,  123,  114,  126,  113,   86,  111,  126,   95,  119,   93,  119,
     93,   89,   75,   80,   71,   68,   59,   54,   61,   37,   21,   33,
     37,   32,   31,   22,   19,   25,   14,   13,   12,   10,    7,   10,
     5,    4,    8,    1,    5,    2,    1,    3,    1,    5,    0,    3,
     1,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
     0,    0,    0,    1])

# These parameter values are known to be the best-fit parameter values
# of a Gaussian to the x- and y-arrays above.
fittedpars = numpy.array([32.3536, 807.863, 117.826])
limit_parnums = numpy.array([0, 1, 2])
maxpars = numpy.array([200, 900, 200])
minpars = numpy.array([1, 0, 0])
hardmaxpars = numpy.array([1.0e+120, 1.0e+120, 1.0e+120])
hardminpars = numpy.array([1.0e-120, -1.0e+120, -1.0e+120])

gfactor = 4.0 * 0.6931471805599453094172321214581766


def freeze_par(pars, parmins, parmaxes, i):
    return (pars, parmins, parmaxes)


def thaw_par(i):
    pass


def report_progress(i, lower, upper):
    pass


def get_par_name(ii):
    pass


# Here's the 1D Gaussian function to use to generate predicted
# data, which will be compared to the y-array above
def gauss_func(p):
    return p[2] * numpy.exp(-1.0 * gfactor *
                            (x - p[1]) * (x - p[1]) / p[0] / p[0])


# Compare the y-array above to values calculated with gauss_func,
# and calculate chi-squared, using Gehrel's approximation for
# estimating errors.
def stat(p):
    errors = 1.0 + numpy.sqrt(y + 0.75)
    fvec = (y - gauss_func(p)) / errors
    return ((fvec * fvec).sum(), )


# Easiest "fit" function for unit tests is actually just to
# return current parameter values.  In this test, that will
# have the effect of making projection act like the old
# uncertainty method, but that is OK for the unit test.
def fitter(scb, pars, parmins, parmaxs):
    return (1, pars, scb(pars)[0])


def test_covar_failures1():
    with pytest.raises(TypeError):
        Covariance().compute(None, fitter, fittedpars,
                             minpars, maxpars,
                             hardminpars, hardmaxpars,
                             limit_parnums, freeze_par, thaw_par,
                             report_progress)


def test_covar_failures2():
    with pytest.raises(RuntimeError):
        Covariance().compute(stat, fitter, None, minpars, maxpars,
                             hardminpars, hardmaxpars, limit_parnums, freeze_par,
                             thaw_par, report_progress, get_par_name)


def test_covar_failures3():
    with pytest.raises(RuntimeError):
        Covariance().compute(stat, fitter, numpy.array([1, 2]),
                             minpars, maxpars, hardminpars, hardmaxpars,
                             limit_parnums, freeze_par, thaw_par,
                             report_progress, get_par_name)


def test_covar_failures4():
    with pytest.raises(RuntimeError):
        Covariance().compute(stat, fitter, fittedpars, numpy.array([1, 2]),
                             maxpars, hardminpars, hardmaxpars, limit_parnums,
                             freeze_par, thaw_par, report_progress, get_par_name)


def test_projection_failures1():
    with pytest.raises(TypeError):
        Projection().compute(stat, None, fittedpars, minpars, maxpars,
                             hardminpars, hardmaxpars, limit_parnums, freeze_par,
                             thaw_par, report_progress, get_par_name)


def test_projection_failures2():
    with pytest.raises(RuntimeError):
        Projection().compute(stat, fitter, None, minpars, maxpars,
                             hardminpars, hardmaxpars, limit_parnums,
                             freeze_par, thaw_par, report_progress, get_par_name)


def test_projection_failures3():
    with pytest.raises(RuntimeError):
        Projection().compute(stat, fitter, numpy.array([1, 2]),
                             minpars, maxpars, hardminpars, hardmaxpars,
                             limit_parnums, freeze_par, thaw_par,
                             report_progress, get_par_name)


def test_projection_failures4():
    with pytest.raises(RuntimeError):
        Projection().compute(stat, fitter, fittedpars, numpy.array([1, 2]),
                             maxpars, hardminpars, hardmaxpars, limit_parnums,
                             freeze_par, thaw_par, report_progress, get_par_name)


# Unlike projection, covariance does not have a "parallel" option.
#
def test_covar():
    standard = numpy.array([[0.4935702,  0.06857833, numpy.nan],
                            [0.06857833, 0.26405554, numpy.nan],
                            [numpy.nan,  numpy.nan,  2.58857314]])

    cov = Covariance()
    results = cov.compute(stat, None, fittedpars,
                          minpars, maxpars,
                          hardminpars, hardmaxpars,
                          limit_parnums, freeze_par, thaw_par,
                          report_progress, get_par_name)

    # These tests used to have a tolerance of 1e-4 but it appears to
    # be able to use a more-restrictive tolerance.
    #
    expected = standard.diagonal()
    assert results[1] == pytest.approx(expected)


# There is no guarantee we can run with parallel=True but try to do so.
#
@pytest.mark.parametrize("parallel", [True, False])
def test_projection(parallel):
    standard_elo = numpy.array([-0.39973743, -0.26390339, -2.08784716])
    standard_ehi = numpy.array([0.39580942,  0.26363223,  2.08789851])

    proj = Projection()
    proj.parallel = parallel
    assert proj.parallel == parallel

    results = proj.compute(stat, fitter, fittedpars,
                           minpars, maxpars,
                           hardminpars, hardmaxpars,
                           limit_parnums, freeze_par, thaw_par,
                           report_progress, get_par_name)

    # These tests used to have a tolerance of 1e-4 but it appears to
    # be able to use a more-restrictive tolerance (given that the
    # "fitter" doesn't do a fit here).
    #
    assert results[0] == pytest.approx(standard_elo)
    assert results[1] == pytest.approx(standard_ehi)


@pytest.mark.parametrize("cls,name",
                         [(Covariance, "Covariance"),
                          (Confidence, "Confidence"),
                          (Projection, "Projection")])
def test_estmethod_repr(cls, name):
    """Simple check"""
    m = cls()
    assert repr(m) == f"<{name} error-estimation method instance '{name.lower()}'>"


def test_estmethod_str_covariance(check_str):
    """Simple check."""
    m = Covariance()
    check_str(str(m),
              ["name        = covariance",
               "sigma       = 1",
               "eps         = 0.01",
               "maxiters    = 200",
               "soft_limits = False"])


def test_estmethod_str_confidence(check_str):
    """Simple check."""
    m = Confidence()
    check_str(str(m),
              ["name         = confidence",
               "sigma        = 1",
               "eps          = 0.01",
               "maxiters     = 200",
               "soft_limits  = False",
               "remin        = 0.01",
               "fast         = False",
               "parallel     = True",
               re.compile(r"^numcores     ?= \d+$"),
               "maxfits      = 5",
               "max_rstat    = 3",
               "tol          = 0.2",
               "verbose      = False",
               "openinterval = False"
               ])


def test_estmethod_str_projection(check_str):
    """Simple check."""
    m = Projection()
    check_str(str(m),
              ["name        = projection",
               "sigma       = 1",
               "eps         = 0.01",
               "maxiters    = 200",
               "soft_limits = False",
               "remin       = 0.01",
               "fast        = False",
               "parallel    = True",
               re.compile(r"^numcores     ?= \d+$"),
               "maxfits     = 5",
               "max_rstat   = 3",
               "tol         = 0.2"
               ])
