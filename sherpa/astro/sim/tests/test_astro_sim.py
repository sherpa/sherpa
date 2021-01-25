#
#  Copyright (C) 2011, 2015, 2018, 2020, 2021
#                Smithsonian Astrophysical Observatory
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

import numpy as np

import pytest

from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.astro import sim

from sherpa.astro.instrument import Response1D
from sherpa.fit import Fit
from sherpa.models.parameter import Parameter
from sherpa.sim import inverse2
from sherpa.stats import CStat
from sherpa.optmethods import NelderMead
from sherpa.estmethods import Covariance

logger = logging.getLogger('sherpa')


@pytest.fixture
def setup(make_data_path):
    from sherpa.astro.io import read_pha
    from sherpa.astro.xspec import XSwabs, XSpowerlaw

    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.CRITICAL)

    pha = make_data_path("refake_0934_1_21_1e4.fak")

    simarf = make_data_path("aref_sample.fits")
    pcaarf = make_data_path("aref_Cedge.fits")

    data = read_pha(pha)
    data.ignore(None, 0.3)
    data.ignore(7.0, None)

    rsp = Response1D(data)
    abs1 = XSwabs('abs1')
    p1 = XSpowerlaw('p1')
    model = rsp(abs1 * p1)

    abs1.nh = 0.092886
    p1.phoindex = 0.994544
    p1.norm = 9.26369

    fit = Fit(data, model, CStat(), NelderMead(), Covariance())

    yield {'simarf': simarf,
           'pcaarf': pcaarf,
           'niter': 10,
           'fit': fit}

    # Reset the logger
    logger.setLevel(old_level)


@requires_xspec
@requires_data
@requires_fits
def test_pragbayes_simarf(setup):
    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler('PragBayes')
    mcmc.set_sampler_opt("simarf", setup['simarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 7)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    niter = setup['niter']
    stats, accept, params = mcmc.get_draws(fit, cov, niter=niter)

    assert params.shape == (3, niter + 1)

    # try:
    #     assert (covar_results.parmaxes < params.std(1)).all()
    # except AssertionError:
    #     print 'covar: ', str(covar_results.parmaxes)
    #     print 'param: ', str(params.std(1))
    #     raise


@requires_xspec
@requires_data
@requires_fits
def test_fullbayes_simarf_fails(setup):
    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler('FullBAYes')
    mcmc.set_sampler_opt("simarf", setup['simarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 7)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    niter = setup['niter']
    with pytest.raises(TypeError) as exc:
        mcmc.get_draws(fit, cov, niter=niter)

    assert str(exc.value) == 'Simulation ARF must be PCA for FullBayes not SIM1DAdd'


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize("sampler", ["pragBayes", "fullbayes"])
def test_pragbayes_pcaarf(sampler, setup):
    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler(sampler)
    mcmc.set_sampler_opt("simarf", setup['pcaarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 5)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    niter = setup['niter']
    stats, accept, params = mcmc.get_draws(fit, cov, niter=niter)

    assert params.shape == (3, niter + 1)

    # try:
    #     assert (covar_results.parmaxes < params.std(1)).all()
    # except AssertionError:
    #     print 'covar: ', str(covar_results.parmaxes)
    #     print 'param: ', str(params.std(1))
    #     raise


@requires_xspec
@requires_data
@requires_fits
@pytest.mark.parametrize("sampler", ["pragBayes", "fullbayes"])
def test_pragbayes_pcaarf_limits(sampler, setup, caplog, reset_seed):
    """Try and trigger limit issues.

    """

    from sherpa.astro.xspec import XSAdditiveModel, XSMultiplicativeModel, \
        XSwabs, XSpowerlaw

    # Set the seed for the RNG. The seed was adjusted to try and make
    # sure the coverage was "good" (i.e. hits parts of
    # sherpa/astro/sim/*bayes.py) whilst still passing the test and
    # reducing the runtime.  This is not a guarantee that this is the
    # "fastest" seed, just that it's one of the better ones I've seen.
    #
    np.random.seed(0x723c)

    class HackAbs(XSwabs):
        """Restrict hard limits"""

        def __init__(self, name='wabs'):
            self.nH = Parameter(name, 'nH', 0.1, 0, 1, 0, 1, '10^22 atoms / cm^2')
            XSMultiplicativeModel.__init__(self, name, (self.nH, ))

    class HackPowerLaw(XSpowerlaw):
        """Restrict hard limits"""
        def __init__(self, name='powerlaw'):
            self.PhoIndex = Parameter(name, 'PhoIndex', 1., 0.95, 1.05, 0.95, 1.05)
            self.norm = Parameter(name, 'norm', 9.2, 8.8, 9.7, 8.8, 9.7)
            XSAdditiveModel.__init__(self, name, (self.PhoIndex, self.norm))

    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler(sampler)
    mcmc.set_sampler_opt("simarf", setup['pcaarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 5)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    # Restrict the parameter values to try and trigger some
    # invalid proposal steps. It's not obvious how the soft,
    # hard, and prior function values all interact.
    #
    myabs = HackAbs()
    mypl = HackPowerLaw()

    pvals = np.asarray(covar_results.parvals)
    pmins = np.asarray(covar_results.parmins)
    pmaxs = np.asarray(covar_results.parmaxes)

    fit.model = myabs * mypl

    fit.model.thawedpars = pvals
    fit.model.thawedparmins = pvals + 2 * pmins  # pmins are < 0
    fit.model.thawedparmaxes = pvals + 2 * pmaxs

    # weight values away from the best-fit (does this actually
    # help?)
    #
    for par in fit.model.pars:
        mcmc.set_prior(par, inverse2)

    niter = setup['niter']
    with caplog.at_level(logging.INFO, logger='sherpa'):
        # Do nothing with the warning at the moment, which could be
        # a RuntimeWarning about the covariance matrix not being
        # positive-semidefinite. This is just needed to make sure
        # we don't trigger the default warning check.
        #
        with pytest.warns(Warning):
            stats, accept, params = mcmc.get_draws(fit, cov, niter=niter)

    # This is a lower bound, in case there's any messages from
    # the sampling (either before or after displaying the
    # 'Using Priors' message).
    #
    nrecords = len(caplog.record_tuples)
    assert nrecords > 3

    i = 0
    while caplog.record_tuples[i][2] != 'Using Priors:':
        i += 1
        assert i < nrecords

    assert i < (nrecords - 3)

    assert caplog.record_tuples[i + 1][2].startswith('wabs.nH: <function inverse2 at ')
    assert caplog.record_tuples[i + 2][2].startswith('powerlaw.PhoIndex: <function inverse2 at ')
    assert caplog.record_tuples[i + 3][2].startswith('powerlaw.norm: <function inverse2 at ')

    # It is not guaranteed what limits/checks we hit
    #
    have_hard_limit = False
    have_reject = False
    for loc, lvl, msg in caplog.record_tuples[i + 4:]:
        if msg.startswith('Draw rejected: parameter boundary exception'):
            have_reject = True
            assert lvl == logging.INFO

        elif msg.startswith('hard '):
            have_hard_limit = True
            assert lvl == logging.WARNING

    assert have_hard_limit
    assert have_reject
