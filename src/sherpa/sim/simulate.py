#
#  Copyright (C) 2010, 2016, 2019 - 2021, 2023, 2024
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


"""
Classes for PPP simulations
"""

from copy import deepcopy
import logging
from typing import Sequence

import numpy as np

from sherpa.stats import Cash, CStat, WStat
from sherpa.optmethods import NelderMead
from sherpa.estmethods import Covariance
from sherpa.fit import Fit
from sherpa.plot import arr2str
from sherpa.sim.sample import NormalParameterSampleFromScaleMatrix
from sherpa.utils import NoNewAttributesAfterInit
from sherpa.utils.parallel import parallel_map_rng
from sherpa.utils.random import poisson_noise
from sherpa.utils.types import ArrayType


# TODO: this should probably use __name__ rather than "sherpa"
logger = logging.getLogger("sherpa")
debug = logger.debug

_tol = np.finfo(float).eps

__all__ = ('LikelihoodRatioTest', 'LikelihoodRatioResults')


class LikelihoodRatioResults(NoNewAttributesAfterInit):
    """The results of a likelihood ratio comparison simulation.

    .. versionchanged:: 4.15.1
       The parnames and parvals attributes have been added. They are
       intended to debug problem cases and so are not displayed by
       default.

    Attributes
    ----------
    ratios : numpy array
       The likelihood ratio for each simulation.
    stats : numpy array
       The fit statistic for the null and alternative models for each
       simulation. The shape is (nsim, 2).
    samples : numpy array
       The parameter samples array for each simulation, with shape
       (nsim, npar).
    lr : number
       The likelihood ratio of the observed data for the null and
       alternate models.
    ppp : number
       The p value of the observed data for the null and alternate
       models.
    null : number
       The fit statistic of the null model on the observed data.
    alt : number
       The fit statistic of the alternate model on the observed data.
    parnames : list of str
       The thawed parameters in the alternate model.
    parvals : ndarray
       The parameter values for each iteration for the alternate
       model, matching the order of parnames. The shape is (nsim,
       len(parnames)).

    Notes
    -----

    The parvals field is useful to check that the simulations have not
    got stuck with certain parameter sets, for instance if the ratio
    value drops to ~ 0 and stays there. If this is the case then the
    analysis can be re-run after adjusting the range (the min or max
    range) of the parameters in question.

    """

    def __init__(self,
                 ratios: ArrayType,
                 stats: ArrayType,
                 samples: ArrayType,
                 lr: float,
                 ppp: float,
                 null: float,
                 alt: float,
                 parnames: Sequence[str],
                 parvals: np.ndarray
                 ) -> None:

        self.ratios = np.asarray(ratios)
        self.stats = np.asarray(stats)
        self.samples = np.asarray(samples)
        self.lr = float(lr)
        self.ppp = float(ppp)
        self.null = float(null)
        self.alt = float(alt)
        self.parnames = parnames
        self.parvals = np.asarray(parvals)
        if len(self.parnames) != self.parvals.shape[1]:
            # This is an unlikely error so do not make it a Sherpa error case
            raise ValueError(f"len(parnames) = {len(self.parnames)}  "
                             f"parvals.shape = {self.parvals.shape}")

        super().__init__()

    def __repr__(self) -> str:
        return '<Likelihood ratio results instance>'

    def __str__(self) -> str:
        samples = arr2str(self.samples)
        stats = arr2str(self.stats)
        ratios = arr2str(self.ratios)

        output = '\n'.join([
                f'samples = {samples}',
                f'stats   = {stats}',
                f'ratios  = {ratios}',
                f'null    = {repr(self.null)}',
                f'alt     = {repr(self.alt)}',
                f'lr      = {repr(self.lr)}',
                f'ppp     = {repr(self.ppp)}'
                ])

        return output

    def format(self) -> str:
        """Convert the object to a string representation for display purposes.

        Returns
        -------
        txt : str
           A string representation of the data stored in the object.

        """
        s = 'Likelihood Ratio Test\n'
        s += f'null statistic   =  {self.null}\n'
        s += f'alt statistic    =  {self.alt}\n'
        s += f'likelihood ratio =  {self.lr}\n'
        if self.ppp == 0.0:
            s += f'p-value          <  {1./len(self.samples)}'
        else:
            s += f'p-value          =  {self.ppp}'
        return s


class LikelihoodRatioTestWorker:
    """
    Worker class for LikelihoodRatioTest
    """
    def __init__(self, null_fit, alt_fit, null_vals, alt_vals):
        self.null_fit = null_fit
        self.alt_fit = alt_fit
        self.null_vals = null_vals
        self.alt_vals = alt_vals

        # Store the original values
        self.null_thawedpars = self.null_fit.model.thawedpars
        self.alt_thawedpars = self.alt_fit.model.thawedpars

    def __call__(self, proposal, rng=None):
        try:
            return LikelihoodRatioTest.calculate(self.null_fit, self.alt_fit,
                                                 proposal, self.null_vals,
                                                 self.alt_vals, rng=rng)
        finally:
            # Ensure the parameters are reset
            self.alt_fit.model.thawedpars = self.alt_thawedpars
            self.null_fit.model.thawedpars = self.null_thawedpars


class LikelihoodRatioTest(NoNewAttributesAfterInit):
    """Likelihood Ratio Test.

    The likelihood ratio [1]_, D, is defined as::

                  (   likelihood for null model	)
        D = -2 ln -----------------------------------
                  ( likelihood for alternative model)

          = -2 ln(likelihood for null model) +
             2 ln(likelihood for alternative model)

    Since the Cash and C fit statistics are defined as -2
    ln(likelihood), the equation reduces to::

        D = statistic for null model -
            statistic for alternative model

    .. versionchanged:: 4.17.0
       The run method can now be called when using the WStat
       statistic.

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Likelihood-ratio_test

    """

    @staticmethod
    def calculate(nullfit, altfit, proposal, null_vals, alt_vals, rng=None):

        # FIXME: only null perturbed?
        nullfit.model.thawedpars = proposal

        # Fake using poisson_noise with null
        fake = poisson_noise(nullfit.data.eval_model(nullfit.model), rng=rng)

        # Set faked data for both nullfit and altfit
        nullfit.data.set_dep(fake)

        # Start the faked fit at initial null best-fit values
        # nullfit.model.thawedpars = null_vals

        # Fit with null model
        nullfr = nullfit.fit()
        debug(nullfr.format())

        null_stat = nullfr.statval
        debug("statistic null = %s", repr(null_stat))

        # nullfit and altfit BOTH point to same faked dataset
        assert id(nullfit.data) == id(altfit.data)
        assert (nullfit.data.get_dep() == altfit.data.get_dep()).all()

        # Start the faked fit at the initial alt best-fit values
        # altfit.model.thawedpars = alt_vals

        debug("proposal: %s", repr(proposal))
        debug("alt model")
        debug(str(altfit.model))

        # Set alt model and fit
        altfr = altfit.fit()
        debug(altfr.format())

        debug(str(altfit.model))

        alt_stat = altfr.statval
        debug("statistic alt = %s", repr(alt_stat))

        LR = -(alt_stat - null_stat)
        debug("LR = %s", repr(LR))

        return [null_stat, alt_stat, LR, altfit.model.thawedpars]

    @staticmethod
    def run(fit, null_comp, alt_comp, conv_mdl=None,
            stat=None, method=None,
            niter=500, numcores=None, rng=None):
        if stat is None:
            stat = CStat()
        if method is None:
            method = NelderMead()

        # There should be a better way to determine this (by querying the
        # stat object in some way). We could use
        #
        #   isinstance(stat, sherpa.stats.Likelihood)
        #
        # which would seem to be a more-appropriate check.
        #
        if not isinstance(stat, (Cash, CStat, WStat)):
            raise TypeError("Sherpa fit statistic must be Cash, CStat, or "
                            "WStat for likelihood ratio test")

        niter = int(niter)

        alt = alt_comp
        null = null_comp

        oldaltvals = alt.thawedpars
        oldnullvals = null.thawedpars

        data = fit.data

        if conv_mdl is not None:
            # Copy the PSF
            null_conv_mdl = deepcopy(conv_mdl)

            alt = conv_mdl(alt_comp)
            if hasattr(conv_mdl, 'fold'):
                conv_mdl.fold(data)

            # Convolve the null model
            null = null_conv_mdl(null_comp)
            if hasattr(null_conv_mdl, 'fold'):
                null_conv_mdl.fold(data)

        nullfit = Fit(data, null, stat, method, Covariance())

        # Fit with null model
        nullfit_results = nullfit.fit()
        debug(nullfit_results.format())

        null_stat = nullfit_results.statval
        null_vals = nullfit_results.parvals

        # Calculate niter samples using null best-fit and covariance
        sampler = NormalParameterSampleFromScaleMatrix()
        samples = sampler.get_sample(nullfit, mycov=None, num=niter, rng=rng)

        # Fit with alt model, null component starts at null's best fit params.
        altfit = Fit(data, alt, stat, method, Covariance())
        altfit_results = altfit.fit()
        debug(altfit_results.format())

        alt_stat = altfit_results.statval
        alt_vals = altfit_results.parvals

        LR = -(alt_stat - null_stat)

        olddep = data.get_dep(filter=False)
        try:
            worker = LikelihoodRatioTestWorker(nullfit, altfit,
                                               null_vals, alt_vals)
            statistics = parallel_map_rng(worker, samples, rng=rng,
                                          numcores=numcores)
        finally:
            data.set_dep(olddep)
            alt.thawedpars = oldaltvals
            null.thawedpars = oldnullvals

        debug("statistic null = %s", repr(null_stat))
        debug("statistic alt = %s", repr(alt_stat))
        debug("LR = %s", repr(LR))

        lrs = []
        stats = []
        thawedpars = []
        for statrow in statistics:
            stats.append(statrow[0:2])
            lrs.append(statrow[2])
            thawedpars.append(statrow[3])

        stats = np.asarray(stats)
        lrs = np.asarray(lrs)
        thawedpars = np.asarray(thawedpars)

        pppvalue = np.sum(lrs > LR) / (1.0 * niter)
        debug('ppp value = %s', str(pppvalue))

        pars = [p.fullname for p in altfit.model.get_thawed_pars()]
        return LikelihoodRatioResults(lrs, stats, samples, LR,
                                      pppvalue, null_stat, alt_stat,
                                      pars, thawedpars)
