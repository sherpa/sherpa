#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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

import warnings

import numpy as np

import pytest

from sherpa.fit import Fit
from sherpa.data import Data1DAsymmetricErrs
from sherpa.optmethods import LevMar
from sherpa.utils.testing import requires_data, requires_fits
from sherpa.models import Const1D, PowLaw1D
from sherpa.stats import Chi2Gehrels
from sherpa.sim import ReSampleData
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
    result = rdata(niter=100, seed=123)

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
    data = ui.get_data(1)
    powlaw1d = PowLaw1D('p1')
    ui.set_model(powlaw1d)
    ui.fit()
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        ui.resample_data(1, 3)
        assert len(w) == 0


@requires_data
@requires_fits
def test_ui(make_data_path):

    infile = make_data_path('gro_delta.txt')
    ui.load_ascii_with_errors(1, infile, delta=True)

    ui.set_stat('leastsq')
    ui.set_model('powlaw1d.p1')
    ui.fit()
    sample = ui.resample_data(1, 10, seed=123)
    for p in ['p1.gamma', 'p1.ampl']:
        assert sample[p] == pytest.approx(RESAMPLE_BENCH_10[p], rel=1e-4)


def test_zero_case():
    """Check what happens when values can be near -1. See #740"""

    xs = np.arange(1, 6)
    ys = np.asarray([0.5, -0.5, 0.3, 0.2, -0.1])
    dyl = np.asarray([2, 2, 2, 2, 2])
    dyh = np.asarray([2, 2, 2, 2, 2])

    data = Data1DAsymmetricErrs('zero', xs, ys, dyl, dyh)
    mdl = Const1D('flat')

    Fit(data, mdl).fit()

    rd = ReSampleData(data, mdl)

    # Both approaches should give the same results (as setting
    # niter explicitly).
    #
    # res = rd.call(niter=10, seed=47)
    res = rd(niter=10, seed=47)

    # Technically this depends on random chance, but it is rather
    # unlikely we'd get 10 -1 values here. Relax the -1 value and
    # just ensure we have more than 1 unique value here (can probably
    # assert nvs == 10 but technically we allow repeated values).
    #
    vs = np.unique(res['flat.c0'])
    nvs = len(vs)
    assert nvs > 1

    # minimal testing of the samples return value
    samples = res['samples']
    assert samples.shape == (10, 5)
