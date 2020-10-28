#
#  Copyright (C) 2011, 2016, 2018, 2020  Smithsonian Astrophysical Observatory
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

from collections import namedtuple
import logging

import numpy
import pytest

from sherpa.data import Data1D
from sherpa.models import Gauss1D, PowLaw1D
from sherpa.fit import Fit
from sherpa.stats import Cash, Chi2DataVar, CStat
from sherpa.optmethods import NelderMead, LevMar
from sherpa.estmethods import Covariance
from sherpa import sim


_max  = numpy.finfo(numpy.float32).max
_tiny = numpy.finfo(numpy.float32).tiny
_eps  = numpy.finfo(numpy.float32).eps


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
_err = numpy.ones(100) * 0.4


@pytest.fixture
def setup():
    data = Data1D('fake', _x, _y, _err)

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

    fit = Fit(data, model, Chi2DataVar(), method, Covariance())
    results = fit.fit()

    for key in ["succeeded", "numpoints", "nfev"]:
        assert _fit_results_bench[key] == int(getattr(results, key))

    for key in ["rstat", "qval", "statval", "dof"]:
        # used rel and abs tol of 1e-7 with numpy allclose
        assert float(getattr(results, key)) == pytest.approx(_fit_results_bench[key])

    for key in ["parvals"]:
        try:
            # used rel and abs tol of 1e-4 with numpy allclose
            assert getattr(results, key) == pytest.approx(_fit_results_bench[key])
        except AssertionError:
            print('parvals bench: ', _fit_results_bench[key])
            print('parvals fit:   ', getattr(results, key))
            print('results', results)
            raise

    fields = ['data', 'model', 'method', 'fit', 'results',
              'covresults', 'dof', 'mu', 'num']
    out = namedtuple('Results', fields)

    out.data = data
    out.model = model
    out.method = method
    out.fit = fit
    out.results = results
    out.covresults = fit.est_errors()
    out.dof = results.dof
    out.mu = numpy.array(results.parvals)
    out.cov = numpy.array(out.covresults.extra_output)
    out.num = 10
    return out


def test_student_t(setup):
    sim.multivariate_t(setup.mu, setup.cov, setup.dof, setup.num)


def test_cauchy(setup):
    sim.multivariate_cauchy(setup.mu, setup.cov, setup.num)


def test_parameter_scale_vector(setup):
    ps = sim.ParameterScaleVector()
    ps.get_scales(setup.fit)


def test_parameter_scale_matrix(setup):
    ps = sim.ParameterScaleMatrix()
    ps.get_scales(setup.fit)


def test_uniform_parameter_sample(setup):
    up = sim.UniformParameterSampleFromScaleVector()
    up.get_sample(setup.fit, num=setup.num)


def test_normal_parameter_sample_vector(setup):
    np = sim.NormalParameterSampleFromScaleVector()
    np.get_sample(setup.fit, num=setup.num)


def test_normal_parameter_sample_matrix(setup):
    np = sim.NormalParameterSampleFromScaleMatrix()
    np.get_sample(setup.fit, num=setup.num)


def test_t_parameter_sample_matrix(setup):
    np = sim.StudentTParameterSampleFromScaleMatrix()
    np.get_sample(setup.fit, setup.dof, num=setup.num)


def test_uniform_sample(setup):
    up = sim.UniformSampleFromScaleVector()
    up.get_sample(setup.fit, num=setup.num)


# This used to have the same name as test_uniform_sample,
# so it has been renamed. Neither apparent to do much
# testing.
#
def test_uniform_sample2(setup):
    sim.uniform_sample(setup.fit, num=setup.num)


def test_normal_sample_vector(setup):
    np = sim.NormalSampleFromScaleVector()
    np.get_sample(setup.fit, num=setup.num)


def test_normal_sample_matrix(setup):
    np = sim.NormalSampleFromScaleMatrix()
    np.get_sample(setup.fit, num=setup.num)


def test_t_sample_matrix(setup):
    np = sim.StudentTSampleFromScaleMatrix()
    np.get_sample(setup.fit, setup.num, setup.dof)


def test_normal_sample(setup):
    sim.normal_sample(setup.fit, num=setup.num, correlate=False)


def test_normal_sample_correlated(setup):
    sim.normal_sample(setup.fit, num=setup.num, correlate=True)

def test_t_sample(setup):
    sim.t_sample(setup.fit, setup.num, setup.dof)

def test_lrt(setup):
    results = sim.LikelihoodRatioTest.run(setup.fit, setup.fit.model.lhs,
                                          setup.fit.model, niter=25)


def test_mh(setup):

    setup.fit.method = NelderMead()
    setup.fit.stat = Cash()
    results = setup.fit.fit()
    results = setup.fit.est_errors()
    cov = results.extra_output

    mcmc = sim.MCMC()

    samplers = mcmc.list_samplers()
    priors = mcmc.list_priors()
    for par in setup.fit.model.pars:
        mcmc.set_prior(par, sim.flat)
        prior = mcmc.get_prior(par)

    sampler = mcmc.get_sampler()
    name = mcmc.get_sampler_name()

    mcmc.set_sampler('MH')

    opt = mcmc.get_sampler_opt('defaultprior')
    mcmc.set_sampler_opt('defaultprior', opt)
    # mcmc.set_sampler_opt('verbose', True)

    log = logging.getLogger("sherpa")
    level = log.level
    log.setLevel(logging.ERROR)
    try:
        stats, accept, params = mcmc.get_draws(setup.fit, cov, niter=1e2)
    finally:
        log.setLevel(level)

def test_metropolisMH(setup):

    setup.fit.method = NelderMead()
    setup.fit.stat = CStat()
    results = setup.fit.fit()
    results = setup.fit.est_errors()
    cov = results.extra_output

    mcmc = sim.MCMC()
    mcmc.set_sampler('MetropolisMH')
    # mcmc.set_sampler_opt('verbose', True)

    log = logging.getLogger("sherpa")
    level = log.level
    log.setLevel(logging.ERROR)
    try:
        stats, accept, params = mcmc.get_draws(setup.fit, cov, niter=1e2)
    finally:
        log.setLevel(level)
