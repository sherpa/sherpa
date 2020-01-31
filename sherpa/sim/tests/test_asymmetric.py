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

import numpy as np
import warnings

import pytest

from sherpa.fit import Fit
from sherpa.data import Data1DAsymmetricErrs
from sherpa.optmethods import LevMar
from sherpa.utils.testing import SherpaTestCase, requires_data, requires_fits
from sherpa.models import Const1D, PowLaw1D
from sherpa.stats import Chi2Gehrels
from sherpa.estmethods import Covariance
from sherpa.sim import ReSampleData
from sherpa.astro import ui

@requires_data
@requires_fits
class test_sim(SherpaTestCase):

    _results_bench_avg = {
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

    _results_bench_rms = {
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

    _resample_bench = np.array([-0.4697257926643954, 0.075012829992575,
                                   177.2066436604025, 60.50264184911246])

    _resample_bench_10 = \
        {'p1.gamma': np.array([-0.2944668147469516, -0.4802566080292259,
                               -0.4529002661840189, -0.5238444562856526,
                               -0.3549169211972059, -0.2911998274429656,
                               -0.5136099861403167, -0.4989971477937067,
                               -0.530802555677153, -0.5645228923616015]),
         'p1.ampl': np.array([60.75374932719593, 172.73752299705234,
                              149.67141746787001, 222.60620753449422,
                              87.33459341923499, 59.25767249100198,
                              212.86125707290233, 194.83475286415504,
                              233.3323867603084, 275.7531516242366])}
    def setUp(self):
        self.method = LevMar()
        self.stat = Chi2Gehrels()
        self.est = Covariance()
        self.gro_fname = self.make_path('gro.txt')
        self.gro_delta_fname = self.make_path('gro_delta.txt')
        return

    def cmp(self, bench, results, tol=1.0e-3):
        for key in ["numpoints", "succeeded", "dof"]:
            assert bench[key] == int(getattr(results, key))
        for key in ["rstat", "qval", "istatval", "statval"]:
            assert np.allclose(float(bench[key]),
                                        float(getattr(results, key)), tol, tol)
        self.assertEqualWithinTol(bench['parvals'], results.parvals, tol)
        parerrs = np.sqrt(results.extra_output['covar'].diagonal())
        self.assertEqualWithinTol(bench['parerrs'], parerrs, tol)

    def fit_asymmetric_err(self, bench, data):
        model = PowLaw1D('p1')
        fit = Fit(data, model, self.stat, self.method, self.est)
        results = fit.fit()
        self.cmp(bench, results)
        return model

    def test_gro_ascii(self):
        ui.load_ascii_with_errors(1, self.gro_fname, delta=False)
        data = ui.get_data(1)
        self.fit_asymmetric_err(self._results_bench_avg, data)

    def test_gro_delta(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        data = ui.get_data(1)
        self.fit_asymmetric_err(self._results_bench_avg, data)

    def test_AsymmetricErrs_avg(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        tmp = ui.get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.fit_asymmetric_err(self._results_bench_avg, data)

    def rms(self, a, b):
        return np.sqrt(a * a + b * b)

    def test_gro_ascii_rms(self):
        ui.load_ascii_with_errors(1, self.gro_fname, func=self.rms,
                               delta=False)
        data = ui.get_data(1)
        self.fit_asymmetric_err(self._results_bench_rms, data)

    def test_gro_delta_rms(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, func=self.rms,
                               delta=True)
        data = ui.get_data(1)
        self.fit_asymmetric_err(self._results_bench_rms, data)

    def test_AsymmetricErrs_rms(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, func=self.rms,
                                     delta=True)
        tmp = ui.get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.fit_asymmetric_err(self._results_bench_rms, data)

    def cmp_resample_data(self, bench, result, tol=1.0e-3):
        gamma = result['p1.gamma']
        ampl = result['p1.ampl']
        self.assertEqualWithinTol(bench, np.array([np.average(gamma),
                                                   np.std(gamma),
                                                   np.average(ampl),
                                                   np.std(ampl)]), tol)

    def resample_data(self, data, bench, results_bench, tol=1.0e-3):
        model = self.fit_asymmetric_err(results_bench, data)
        rd = ReSampleData(data, model)
        result = rd(niter=100, seed=123)
        self.cmp_resample_data(bench, result)
        self.assertEqualWithinTol(results_bench['parvals'],
                                  model.thawedpars, tol)

    def test_AsymmetricErrors_resample_avg(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        tmp = ui.get_data(1)
        data = Data1DAsymmetricErrs(1, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.resample_data(data, self._resample_bench,
                           self._results_bench_avg)

    def test_AsymmetricErrors_resample_rms(self):
        ui.load_ascii_with_errors(1, self.gro_delta_fname, delta=True,
                               func=self.rms)
        tmp = ui.get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.resample_data(data, self._resample_bench,
                           self._results_bench_rms)

    def test_ui(self, tol=1.0e-3):
        # from shepa.astro.ui import *
        ui.load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        ui.set_stat('leastsq')
        ui.set_model('powlaw1d.p1')
        ui.fit()
        sample = ui.resample_data(1, 10, seed=123)
        self.assertEqualWithinTol(self._resample_bench_10['p1.gamma'],
                                  sample['p1.gamma'])
        self.assertEqualWithinTol(self._resample_bench_10['p1.ampl'],
                                  sample['p1.ampl'])

    def test_warning(self):
        ui.load_ascii_with_errors(1, self.gro_fname)
        data = ui.get_data(1)
        powlaw1d = PowLaw1D('p1')
        ui.set_model(powlaw1d)
        fit = Fit(data, powlaw1d)
        results = fit.fit()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            ui.resample_data(1, 3)
            assert len(w) == 0


def test_zero_case():
    """Check what happens when values can be near -1."""

    xs = np.arange(1, 6)
    ys = np.asarray([0.5, -0.5, 0.3, 0.2, -0.1])
    dyl = np.asarray([2, 2, 2, 2, 2])
    dyh = np.asarray([2, 2, 2, 2, 2])

    data = Data1DAsymmetricErrs('zero', xs, ys, dyl, dyh, dyh)
    mdl = Const1D('flat')

    f = Fit(data, mdl)

    rd = ReSampleData(data, mdl)
    res = rd.call(niter=10, seed=47)

    # Technically this depends on random chance, but it is rather
    # unlikely we'd get 10 -1 values here. Relax the -1 value and
    # just ensure we have more than 1 unique value here (can probably
    # assert nvs == 10 but technically we allow repeated values).
    #
    vs = np.unique(res['flat.c0'])
    nvs = len(vs)
    assert nvs > 1
