#
#  Copyright (C) 2019 - 2021, 2023, 2024
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

import logging
import warnings

import numpy as np

import pytest

from sherpa.fit import Fit
from sherpa.data import Data1D, Data1DAsymmetricErrs, Data2D
from sherpa.optmethods import LevMar
from sherpa.utils.testing import requires_data, requires_fits
from sherpa.models import Const1D, Const2D, PowLaw1D
from sherpa.stats import Chi2Gehrels
from sherpa.sim import ReSampleData
from sherpa.utils.logging import SherpaVerbosity
from sherpa.astro import ui


RESULTS_BENCH_AVG = {
    'rstat': 1.4361549916463103,
    'qval': 0.015715747140941,
    'numpoints': 61,
    'succeeded': 1,
    'dof': 59,
    'istatval': 248263792.37785548,
    'statval': 84.73314450713231,
    'parnames': ('p1.gamma', 'p1.ampl'),
    'parvals': np.array([-0.5983957573984262, 332.5326662832957]),
    'parerrs': np.array([2.07753267e-02, 3.82647498e+01])
}

RESULTS_BENCH_RMS = {
    'rstat': 0.6654892642692887,
    'qval': 0.9776806145696532,
    'numpoints': 61,
    'succeeded': 1,
    'dof': 59,
    'istatval': 114123785.83302237,
    'statval': 39.26386659188803,
    'parnames': ('p1.gamma', 'p1.ampl'),
    'parvals': np.array([-0.5984085135559107, 333.537172525073]),
    'parerrs': np.array([3.01965392e-02, 5.57390091e+01])
}

# The RESMPLE_BENCH* values were calculated before the code was
# changed so that each fit starts at the best-fit location,
# rather than starting at the best-fit location from the
# previous fit. The relative tolerances have been changed
# to 1e-4 to keep these values as the reference.
#
RESAMPLE_BENCH = np.array([-4.27921009e-01, 1.54966801e-01,
                           1.67211651e+02, 9.01622079e+01])

RESAMPLE_BENCH_10 = \
    {'p1.gamma': np.array([-0.5328294503498137, -0.5001748057127161,
                           -0.3215221918987224, -0.574564441010874,
                           -0.43821524305367615, -0.426883933868221,
                           -0.3708382753872252, -0.5092469458936105,
                           -0.4460289613152441, -0.4350294472721625]),
     'p1.ampl': np.array([234.30388600765593, 192.1359747458908,
                          70.35141613845244, 291.1140203785552,
                          138.67404346944733, 128.917712169535,
                          93.98512569969482, 207.7714117742187,
                          142.61963597065656, 136.00140089942673])}


def compare(bench, results):
    for key in ["numpoints", "succeeded", "dof"]:
        assert getattr(results, key) == bench[key]

    for key in ["rstat", "qval", "istatval", "statval", "parvals"]:
        assert getattr(results, key) == pytest.approx(bench[key])

    parerrs = np.sqrt(results.extra_output['covar'].diagonal())
    assert parerrs == pytest.approx(bench['parerrs'])


def fit_asymmetric_err(bench, data):
    model = PowLaw1D('p1')
    fit = Fit(data, model, Chi2Gehrels(), LevMar())
    results = fit.fit()
    compare(bench, results)
    return model


def resample_data(data, bench, results_bench):
    model = fit_asymmetric_err(results_bench, data)
    rdata = ReSampleData(data, model)
    result = rdata(niter=100, seed=123,
                   rng=np.random.RandomState(123))

    gamma = result['p1.gamma']
    ampl = result['p1.ampl']
    got = np.array([np.average(gamma),
                    np.std(gamma),
                    np.average(ampl),
                    np.std(ampl)])
    assert got == pytest.approx(bench, rel=1e-4)

    assert model.thawedpars == pytest.approx(results_bench['parvals'])


def rms(a, b):
    """combine errors"""
    return np.sqrt(a * a + b * b)


# Basic data values from the gro*txt file. Note the data is unordered.
#
GRO_X = np.asarray(
    [0.00310, 0.00284, 0.00413, 0.00469, 0.00451, 0.00465, 0.00598,
     0.00604, 0.00605, 0.00559, 0.00555, 0.00567, 0.00589, 0.00455,
     0.00467, 0.00409, 0.00448, 0.00436, 0.00371, 0.00303, 0.00311,
     0.00295, 0.00319, 0.00265, 0.00199, 0.00172, 0.00170, 0.00155,
     0.00155, 0.00143, 0.00157, 0.00111, 0.00082, 0.00122, 0.00090,
     0.00097, 0.00095, 0.00068, 0.00044, 0.00059, 0.00036, 0.00093,
     0.00099, 0.00023, 0.00046, 0.00044, 0.00032, 0.00046, 0.00040,
     0.00042, 0.00041, 0.00036, 0.00041, 0.00037, 0.00051, 0.00032,
     0.00019, 0.00043, 0.00042, 0.00041, 0.00029])

GRO_Y = np.asarray(
    [13.383, 13.352, 14.717, 14.275, 13.169, 12.591, 14.218, 16.317,
     15.154, 13.961, 17.997, 15.528, 15.946, 11.378, 9.605, 11.602,
     10.888, 11.471, 10.348, 10.349, 10.009, 10.035, 10.231, 10.611,
     7.987, 7.618, 7.189, 10.245, 8.336, 7.157, 7.382, 5.667,
     4.278, 6.179, 5.213, 6.174, 5.570, 5.129, 6.937, 7.539,
     3.220, 5.956, 14.412, 5.162, 3.093, 3.214, 2.827, 2.242,
     2.489, 2.812, 2.952, 3.436, 2.167, 2.748, 3.073, 2.177,
     1.146, 3.575, 0.510, 4.194, 1.393])

# These are from gro_delta.txt
GRO_YLO = np.asarray(
    [1.004, 0.915, 1.725, 0.757, 0.848, 0.423, 0.97, 0.465, 0.771,
     1.119, 1.039, 1.054, 0.85, 1.391, 0.701, 0.785, 0.533, 1.013,
     2.133, 2.881, 0.998, 0.953, 1.454, 1.684, 2.087, 0.848, 0.787,
     1.8, 2.546, 0.733, 1.182, 0.589, 0.835, 0.827, 0.488, 2.667,
     0.79, 0.73, 6.937, 1.248, 0.435, 4.843, 14.402, 1.605, 0.505,
     1.854, 0.502, 2.197, 0.504, 1.854, 0.564, 3.436, 2.167, 0.49,
     1.353, 0.488, 0.352, 3.575, 0.204, 1.704, 1.393])

GRO_YHI = np.asarray(
    [1.619, 0.971, 0.976, 0.818, 1.15, 1.247, 0.47, 0.558, 0.866,
     1.111, 1.17, 0.531, 0.629, 0.669, 2.437, 1.147, 1.117, 0.749,
     1.078, 2.417, 0.541, 0.748, 0.786, 0.721, 1.258, 0.525, 1.072,
     1.015, 1.448, 1.711, 1.147, 1.543, 2.396, 1.862, 1.127, 0.989,
     0.911, 0.761, 12.461, 1.609, 1.083, 3.88, 8.217, 79.687, 0.991,
     1.013, 1.25, 2.128, 0.468, 2.062, 0.803, 3.207, 1.507, 0.782,
     3.887, 0.426, 0.858, 2.476, 3.019, 2.878, 1.025])


@requires_data
@requires_fits
@pytest.mark.parametrize("filename, delta", [('gro.txt', False),
                                             ('gro_delta.txt', True)])
def test_load_ascii(filename, delta, make_data_path):
    """Check the arrays get read in."""
    infile = make_data_path(filename)
    ui.load_ascii_with_errors(1, infile, delta=delta)
    data = ui.get_data(1)

    assert isinstance(data, Data1DAsymmetricErrs)
    N = 61
    assert len(data.x) == N
    assert len(data.y) == N
    assert len(data.staterror) == N
    assert data.syserror is None
    assert len(data.elo) == N
    assert len(data.ehi) == N

    assert data.x == pytest.approx(GRO_X)
    assert data.y == pytest.approx(GRO_Y)
    assert data.elo == pytest.approx(GRO_YLO)
    assert data.ehi == pytest.approx(GRO_YHI)

    ymid = (GRO_YLO + GRO_YHI) / 2
    assert data.staterror == pytest.approx(ymid)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename, delta", [('gro.txt', False),
                                             ('gro_delta.txt', True)])
def test_load_ascii(filename, delta, make_data_path):
    infile = make_data_path(filename)
    ui.load_ascii_with_errors(1, infile, delta=delta)
    data = ui.get_data(1)
    fit_asymmetric_err(RESULTS_BENCH_AVG, data)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename, delta", [('gro.txt', False),
                                             ('gro_delta.txt', True)])
def test_load_ascii_defaultid(filename, delta, make_data_path):
    """Use the default id"""
    infile = make_data_path(filename)
    ui.load_ascii_with_errors(infile, delta=delta)
    data = ui.get_data()
    fit_asymmetric_err(RESULTS_BENCH_AVG, data)


@requires_data
@requires_fits
@pytest.mark.parametrize("filename, delta", [('gro.txt', False),
                                             ('gro_delta.txt', True)])
def test_load_ascii_rms(filename, delta, make_data_path):
    infile = make_data_path(filename)
    ui.load_ascii_with_errors(1, infile, delta=delta, func=rms)
    data = ui.get_data(1)
    fit_asymmetric_err(RESULTS_BENCH_RMS, data)


@requires_data
@requires_fits
def test_constructor_avg(make_data_path):
    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True)
    base = ui.get_data(1)
    data = Data1DAsymmetricErrs(2, base.x, base.y, base.elo, base.ehi,
                                base.staterror, base.syserror)
    fit_asymmetric_err(RESULTS_BENCH_AVG, data)


@requires_data
@requires_fits
def test_constructor_rms(make_data_path):
    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True, func=rms)
    base = ui.get_data(1)
    data = Data1DAsymmetricErrs(2, base.x, base.y, base.elo, base.ehi,
                                base.staterror, base.syserror)
    fit_asymmetric_err(RESULTS_BENCH_RMS, data)


@requires_data
@requires_fits
def test_resample_avg(make_data_path):
    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True)
    base = ui.get_data(1)
    data = Data1DAsymmetricErrs(2, base.x, base.y, base.elo, base.ehi,
                                base.staterror, base.syserror)
    resample_data(data, RESAMPLE_BENCH, RESULTS_BENCH_AVG)


@requires_data
@requires_fits
def test_resample_rmd(make_data_path):
    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True, func=rms)
    base = ui.get_data(1)
    data = Data1DAsymmetricErrs(2, base.x, base.y, base.elo, base.ehi,
                                base.staterror, base.syserror)
    resample_data(data, RESAMPLE_BENCH, RESULTS_BENCH_RMS)


@requires_data
@requires_fits
def test_warning(make_data_path):

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    powlaw1d = PowLaw1D('p1')
    ui.set_model(powlaw1d)
    ui.fit()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ui.resample_data(1, 3)
        assert len(w) == 0


@requires_data
@requires_fits
def test_ui(make_data_path, clean_astro_ui):

    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True)

    ui.set_stat('leastsq')
    ui.set_model('powlaw1d.p1')
    ui.fit()

    ui.set_rng(np.random.RandomState(123))
    sample = ui.resample_data(1, 10, seed=123)
    for p in ['p1.gamma', 'p1.ampl']:
        assert sample[p] == pytest.approx(RESAMPLE_BENCH_10[p], rel=1e-4)


@requires_data
@requires_fits
def test_ui_plot_data(make_data_path, clean_astro_ui, all_plot_backends):
    """Just check we can call plot_data."""

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    ui.plot_data()


@requires_data
@requires_fits
def test_to_plot(make_data_path, clean_astro_ui):
    """What does to_plot return?

    We use the ui code to make it easy to create the class.
    """

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)

    data = ui.get_data()
    (x, y, yerr, xerr, xlabel, ylabel) = data.to_plot()

    N = 61
    assert x.shape == (N,)
    assert y.shape == (N,)
    assert len(yerr) == 2
    assert yerr[0].shape == (N,)
    assert yerr[1].shape == (N,)

    assert x == pytest.approx(GRO_X)
    assert y == pytest.approx(GRO_Y)
    assert yerr[0] == pytest.approx(GRO_YLO)
    assert yerr[1] == pytest.approx(GRO_YHI)


@requires_data
@requires_fits
def test_ui_get_data_plot(make_data_path, clean_astro_ui):
    """Just check we can call get_data_plot"""

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    plot = ui.get_data_plot()

    N = 61
    assert plot.x.shape == (N,)
    assert plot.y.shape == (N,)
    assert len(plot.yerr) == 2
    assert plot.yerr[0].shape == (N,)
    assert plot.yerr[1].shape == (N,)

    assert plot.x == pytest.approx(GRO_X)
    assert plot.y == pytest.approx(GRO_Y)
    assert plot.yerr[0] == pytest.approx(GRO_YLO)
    assert plot.yerr[1] == pytest.approx(GRO_YHI)


@requires_data
@requires_fits
def test_ui_plot_filtered_data(make_data_path, clean_astro_ui, all_plot_backends):
    """Just check we can call plot_data after a filter

    See issue #2090
    """

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    ui.ignore(lo=0.005)
    ui.plot_data()


@requires_data
@requires_fits
def test_to_plot_filtered(make_data_path, clean_astro_ui):
    """What does to_plot return?

    We use the ui code to make it easy to create the class.

    See issue #2090
    """

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    ui.ignore(lo=0.005)

    data = ui.get_data()
    (x, y, yerr, xerr, xlabel, ylabel) = data.to_plot()

    N = 54
    assert x.shape == (N,)
    assert y.shape == (N,)
    assert len(yerr) == 2
    assert yerr[0].shape == (N,)
    assert yerr[1].shape == (N,)

    idx = GRO_X < 0.005
    assert x == pytest.approx(GRO_X[idx])
    assert y == pytest.approx(GRO_Y[idx])
    assert yerr[0] == pytest.approx(GRO_YLO[idx])
    assert yerr[1] == pytest.approx(GRO_YHI[idx])


@requires_data
@requires_fits
def test_ui_get_data_plot_filtered(make_data_path, clean_astro_ui):
    """Just check we can call get_data_plot after a filter

    See issue #2090
    """

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    ui.ignore(lo=0.005)
    plot = ui.get_data_plot()

    N = 54
    assert plot.x.shape == (N,)
    assert plot.y.shape == (N,)
    assert len(plot.yerr) == 2
    assert plot.yerr[0].shape == (N,)
    assert plot.yerr[1].shape == (N,)

    idx = GRO_X < 0.005
    assert plot.x == pytest.approx(GRO_X[idx])
    assert plot.y == pytest.approx(GRO_Y[idx])
    assert plot.yerr[0] == pytest.approx(GRO_YLO[idx])
    assert plot.yerr[1] == pytest.approx(GRO_YHI[idx])


@requires_data
@requires_fits
def test_ui_get_staterror(make_data_path, clean_astro_ui):
    """What does get_staterror return?"""

    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    err = ui.get_staterror()

    mid = (GRO_YLO + GRO_YHI) / 2
    assert err == pytest.approx(mid)


@requires_data
@requires_fits
def test_ui_get_staterror_filtered(make_data_path, clean_astro_ui):
    """What does get_staterror return after a filter?

    See issue #2090 (although is not a part of this issue).
    """

    LOVAL = 0.0023
    HIVAL = 0.0033
    infile = make_data_path('gro.txt')
    ui.load_ascii_with_errors(1, infile)
    ui.notice(LOVAL, HIVAL)
    err = ui.get_staterror(filter=True)

    assert len(err) == 7
    mid = (GRO_YLO + GRO_YHI) / 2
    idx = (GRO_X >= LOVAL) & (GRO_X < HIVAL)
    assert err == pytest.approx(mid[idx])


def test_zero_case():
    """Check what happens when values can be near -1. See #740"""

    xs = np.arange(1, 6)
    ys = np.asarray([0.5, -0.5, 0.3, 0.2, -0.1])
    dyl = np.asarray([2, 2, 2, 2, 2])
    dyh = np.asarray([2, 2, 2, 2, 2])

    data = Data1DAsymmetricErrs('zero', xs, ys, dyl, dyh)
    mdl = Const1D('flat')

    bestfit = Fit(data, mdl).fit()

    rd = ReSampleData(data, mdl)

    # Both approaches should give the same results (as setting
    # niter explicitly).
    #
    # res = rd.call(niter=10, seed=47)
    res = rd(niter=10, seed=47,
             rng=np.random.RandomState(47))

    # Technically this depends on random chance, but it is rather
    # unlikely we'd get 10 -1 values here. Relax the -1 value and
    # just ensure we have more than 1 unique value here (can probably
    # assert nvs == 10 but technically we allow repeated values).
    #
    vs = np.unique(res['flat.c0'])
    nvs = len(vs)
    assert nvs > 1

    # Minimal testing of the other return values. We do assume that
    # we have not found a better fit than the fit.
    #
    samples = res['samples']
    stats = res['statistic']
    assert samples.shape == (10, 5)
    assert stats.shape == (10, )
    assert (stats >= bestfit.statval).all()


def test_resample_supports_data1d(caplog, check_str):

    orig = Data1D("orig", [1, 2, 3], [4, 2, 5], [0.1, 0.2, 0.5])
    model = Const1D("mdl")
    resampler = ReSampleData(orig, model)

    with SherpaVerbosity("INFO"):
        # seed is unused of rng is set
        res = resampler(niter=10, seed=None,
                        rng=np.random.RandomState(123))

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa"
    assert lvl == logging.INFO
    check_str(msg,
              ["mdl.c0 : avg = 1.0343987745935705 , std = 0.5208696279243179  # doctest: +FLOAT_CMP"])


def test_resample_fails_unsupported_data():

    orig = Data2D("orig", [1, 2, 3], [1, 1, 2], [4, 3, 2])
    model = Const2D("mdl")

    with pytest.raises(NotImplementedError,
                       match="^ReSampleData is only implemented for 1D data, got <class 'sherpa.data.Data2D'> instead.$"):
        ReSampleData(orig, model)
