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

import numpy

from sherpa.fit import Fit
from sherpa.data import Data1DAsymmetricErrs
from sherpa.optmethods import LevMar
from sherpa.utils.testing import SherpaTestCase
from sherpa.models import PowLaw1D
from sherpa.stats import Chi2Gehrels
from sherpa.estmethods import Covariance
from sherpa.astro.ui import load_ascii_with_errors, get_data
from sherpa.sim import ReSampleData

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
        'parvals': numpy.array([-0.5983957573984262, 332.5326662832957]),
        'parerrs': numpy.array([2.07753267e-02, 3.82647498e+01])
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
        'parvals': numpy.array([-0.5984085135559107, 333.537172525073]),
        'parerrs': numpy.array([3.01965392e-02, 5.57390091e+01])
    }

    _resample_bench = numpy.array([-0.4697257926643954, 0.075012829992575,
                                   177.2066436604025, 60.50264184911246])
    
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
            assert numpy.allclose(float(bench[key]),
                                        float(getattr(results, key)), tol, tol)
        for index, val in enumerate(bench['parvals']):
            assert numpy.allclose(val, results.parvals[index], tol, tol)

        parerrs = numpy.sqrt(results.extra_output['covar'].diagonal())
        for index, val in enumerate(bench['parerrs']):
            assert numpy.allclose(val, parerrs[index], tol, tol)
        
    def fit_asymmetric_err(self, bench, data):
        model = PowLaw1D('p1')
        fit = Fit(data, model, self.stat, self.method, self.est)
        results = fit.fit()
        self.cmp(bench, results)
        
    def test_gro_ascii(self):
        load_ascii_with_errors(1, self.gro_fname, delta=False)
        data = get_data(1)
        self.fit_asymmetric_err(self._results_bench_avg, data)

    def test_gro_delta(self):
        load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        data = get_data(1)        
        self.fit_asymmetric_err(self._results_bench_avg, data)

    def test_AsymmetricErrs_avg(self):
        load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        tmp = get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.fit_asymmetric_err(self._results_bench_avg, data)        

    def rms(self, a, b):
        return numpy.sqrt(a * a + b * b)
        
    def test_gro_ascii_rms(self):
        load_ascii_with_errors(1, self.gro_fname, func=self.rms,
                               delta=False)
        data = get_data(1)
        self.fit_asymmetric_err(self._results_bench_rms, data)

    def test_gro_delta_rms(self):
        load_ascii_with_errors(1, self.gro_delta_fname, func=self.rms,
                               delta=True)
        data = get_data(1)
        self.fit_asymmetric_err(self._results_bench_rms, data)

    def test_AsymmetricErrs_rms(self):
        load_ascii_with_errors(1, self.gro_delta_fname, func=self.rms,
                                     delta=True)
        tmp = get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.fit_asymmetric_err(self._results_bench_rms, data)        

    def cmp_resample_data(self, bench, rd, tol=1.0e-3):
        for a, b in zip(bench, rd):
            assert numpy.allclose(bench, rd, tol, tol)

    def resample_data(self, data, bench):
        model = PowLaw1D('p1')        
        rd = ReSampleData(data, model)
        result = rd(niter=100)
        self.cmp_resample_data(bench, result)
        
    def test_AsymmetricErros_resample_avg(self):
        load_ascii_with_errors(1, self.gro_delta_fname, delta=True)
        tmp = get_data(1)
        data = Data1DAsymmetricErrs(1, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.resample_data(data, self._resample_bench)

    def test_AsymmetricErros_resample_rms(self):
        load_ascii_with_errors(1, self.gro_delta_fname, delta=True,
                               func=self.rms)
        tmp = get_data(1)
        data = Data1DAsymmetricErrs(2, tmp.x, tmp.y, tmp.elo,
                                    tmp.ehi, tmp.staterror, tmp.syserror)
        self.resample_data(data, self._resample_bench)
        
