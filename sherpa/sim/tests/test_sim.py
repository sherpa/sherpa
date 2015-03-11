# 
#  Copyright (C) 2011  Smithsonian Astrophysical Observatory
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
import logging

from sherpa.data import Data1D
from sherpa.models import Gauss1D, PowLaw1D
from sherpa.fit import Fit
from sherpa.stats import Cash, Chi2DataVar, CStat
from sherpa.optmethods import NelderMead, LevMar
from sherpa.estmethods import Covariance
from sherpa.sim import *

from sherpa.utils import SherpaTestCase, SherpaTest


_max  = numpy.finfo(numpy.float32).max
_tiny = numpy.finfo(numpy.float32).tiny
_eps  = numpy.finfo(numpy.float32).eps


class test_sim(SherpaTestCase):

    _fit_results_bench = {
        'rstat': 89.29503933428586,
        'qval': 0.0,
        'succeeded': 1,
        'numpoints': 100,
        'dof': 95,
        'nfev': 93,
        'statval': 8483.0287367571564,
        'parnames': ['p1.gamma', 'p1.ampl', 'g1.fwhm',
                     'g1.pos', 'g1.ampl'], 
        'parvals': numpy.array(
            [1.0701938169914813,
             9.1826254677279469,
             2.5862083052721028,
             2.601619746022207,
             47.262657692418749])
        }

    _x = numpy.arange(0.1, 10.1, 0.1)
    _y = numpy.array(
        [ 114, 47, 35, 30, 40, 27, 30, 26, 24, 20, 26, 35,
          29, 28, 34, 36, 43, 39, 33, 47, 44, 46, 53, 56,
          52, 53, 49, 57, 49, 36, 33, 42, 49, 45, 42, 32,
          31, 34, 18, 24, 25, 11, 17, 17, 11,  9,  8,  5,
           4, 10,  3,  4,  6,  3,  0,  2,  4,  4,  0,  1,
           2,  0,  3,  3,  0,  2,  1,  2,  3,  0,  1,  0,
           1,  0,  0,  1,  3,  3,  0,  2,  0,  0,  1,  2,
           0,  1,  0,  1,  1,  0,  1,  1,  1,  1,  1,  1,
           1,  0,  1,  0
          ]
        )
    _err = numpy.ones(100)*0.4


    def setUp(self):
        data = Data1D('fake', self._x, self._y, self._err)

        g1 = Gauss1D('g1')
        g1.fwhm.set(1.0, _tiny, _max, frozen=False)
        g1.pos.set(1.0, -_max, _max, frozen=False)
        g1.ampl.set(1.0, -_max, _max, frozen=False)
        p1 = PowLaw1D('p1')
        p1.gamma.set(1.0, -10, 10, frozen=False)
        p1.ampl.set(1.0, 0.0, _max, frozen=False)
        p1.ref.set(1.0, -_max, _max, frozen=True)      
        model = p1 + g1

        method = LevMar()
        method.config['maxfev'] = 10000
        method.config['ftol'] = float(_eps)
        method.config['epsfcn'] = float(_eps)
        method.config['gtol'] = float(_eps)
        method.config['xtol'] = float(_eps)
        method.config['factor'] = float(100)

        self.fit = Fit(data, model, Chi2DataVar(), method, Covariance())
        results = self.fit.fit()

        for key in ["succeeded", "numpoints", "nfev"]:
            assert self._fit_results_bench[key] == int(getattr(results, key))

        for key in ["rstat", "qval", "statval", "dof"]:
            assert numpy.allclose(float(self._fit_results_bench[key]),
                                  float(getattr(results, key)),
                                  1.e-7, 1.e-7)

        for key in ["parvals"]:
            try:
                assert numpy.allclose(self._fit_results_bench[key],
                                      getattr(results, key),
                                      1.e-4, 1.e-4)
            except AssertionError:
                print 'parvals bench: ', self._fit_results_bench[key]
                print 'parvals fit:   ', getattr(results, key)
                print 'results', results
                raise

        covresults = self.fit.est_errors()
        self.dof = results.dof
        self.mu = numpy.array(results.parvals)
        self.cov = numpy.array(covresults.extra_output)
        self.num = 10


    def test_student_t(self):
        multivariate_t(self.mu, self.cov, self.dof, self.num)
        
    def test_cauchy(self):
        multivariate_cauchy(self.mu, self.cov, self.num)       

    def test_parameter_scale_vector(self):
        ps = ParameterScaleVector()
        ps.get_scales(self.fit)

    def test_parameter_scale_matrix(self):
        ps = ParameterScaleMatrix()
        ps.get_scales(self.fit)

    def test_uniform_parameter_sample(self):
        up = UniformParameterSampleFromScaleVector()
        up.get_sample(self.fit, num=self.num)
        
    def test_normal_parameter_sample_vector(self):
        np = NormalParameterSampleFromScaleVector()
        np.get_sample(self.fit, num=self.num)

    def test_normal_parameter_sample_matrix(self):
        np = NormalParameterSampleFromScaleMatrix()
        np.get_sample(self.fit, num=self.num)

    def test_t_parameter_sample_matrix(self):
        np = StudentTParameterSampleFromScaleMatrix()
        np.get_sample(self.fit, self.dof, num=self.num)

    def test_uniform_sample(self):
        up = UniformSampleFromScaleVector()
        up.get_sample(self.fit, num=self.num)
        
    def test_normal_sample_vector(self):
        np = NormalSampleFromScaleVector()
        np.get_sample(self.fit, num=self.num)

    def test_normal_sample_matrix(self):
        np = NormalSampleFromScaleMatrix()
        np.get_sample(self.fit, num=self.num)

    def test_t_sample_matrix(self):
        np = StudentTSampleFromScaleMatrix()
        np.get_sample(self.fit, self.num, self.dof)

    def test_uniform_sample(self):
        uniform_sample(self.fit, num=self.num)

    def test_normal_sample(self):
        normal_sample(self.fit, num=self.num, correlate=False)

    def test_normal_sample_correlated(self):
        normal_sample(self.fit, num=self.num, correlate=True)

    def test_t_sample(self):
        t_sample(self.fit, self.num, self.dof)

    def test_lrt(self):
        results = LikelihoodRatioTest.run(self.fit, self.fit.model.lhs,
                                          self.fit.model, niter=25)


    def test_mh(self):

        self.fit.method = NelderMead()
        self.fit.stat = Cash()
        results = self.fit.fit()
        results = self.fit.est_errors()
        cov = results.extra_output

        mcmc = MCMC()
        
        samplers = mcmc.list_samplers()
        priors = mcmc.list_priors()
        for par in self.fit.model.pars:
            mcmc.set_prior(par, flat)
            prior = mcmc.get_prior(par)

        sampler = mcmc.get_sampler()
        name = mcmc.get_sampler_name()

        mcmc.set_sampler('MH')

        opt = mcmc.get_sampler_opt('defaultprior')
        mcmc.set_sampler_opt('defaultprior', opt)
        #mcmc.set_sampler_opt('verbose', True)

        log = logging.getLogger("sherpa")
        level = log.level
        log.setLevel(logging.ERROR)
        stats, accept, params = mcmc.get_draws(self.fit, cov, niter=1e2)
        log.setLevel(level)


    def test_metropolisMH(self):

        self.fit.method = NelderMead()
        self.fit.stat = CStat()
        results = self.fit.fit()
        results = self.fit.est_errors()
        cov = results.extra_output

        mcmc = MCMC()
        mcmc.set_sampler('MetropolisMH')
        #mcmc.set_sampler_opt('verbose', True)
        
        log = logging.getLogger("sherpa")
        level = log.level
        log.setLevel(logging.ERROR)
        stats, accept, params = mcmc.get_draws(self.fit, cov, niter=1e2)
        log.setLevel(level)
 

    def tearDown(self):
        pass


if __name__ == '__main__':

    import sherpa.sim as sim
    SherpaTest(sim).test()
