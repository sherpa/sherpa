#
#  Copyright (C) 2011, 2015, 2018, 2020
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

import pytest

from sherpa.utils.testing import requires_data, requires_fits, requires_xspec
from sherpa.astro import sim

from sherpa.astro.instrument import Response1D
from sherpa.fit import Fit
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
           'fit': fit}

    # Reset the logger
    logger.setLevel(old_level)


@requires_xspec
@requires_data
@requires_fits
def test_pragbayes_simarf(setup):
    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler("PragBayes")
    mcmc.set_sampler_opt("simarf", setup['simarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 7)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    niter = 10
    stats, accept, params = mcmc.get_draws(fit, cov, niter=niter)
    # try:
    #     assert (covar_results.parmaxes < params.std(1)).all()
    # except AssertionError:
    #     print 'covar: ', str(covar_results.parmaxes)
    #     print 'param: ', str(params.std(1))
    #     raise


@requires_xspec
@requires_data
@requires_fits
def test_pragbayes_pcaarf(setup):
    fit = setup['fit']

    mcmc = sim.MCMC()
    mcmc.set_sampler("pragBayes")
    mcmc.set_sampler_opt("simarf", setup['pcaarf'])
    mcmc.set_sampler_opt("p_M", 0.5)
    mcmc.set_sampler_opt("nsubiter", 5)

    covar_results = fit.est_errors()
    cov = covar_results.extra_output

    niter = 10
    stats, accept, params = mcmc.get_draws(fit, cov, niter=niter)
    # try:
    #     assert (covar_results.parmaxes < params.std(1)).all()
    # except AssertionError:
    #     print 'covar: ', str(covar_results.parmaxes)
    #     print 'param: ', str(params.std(1))
    #     raise
