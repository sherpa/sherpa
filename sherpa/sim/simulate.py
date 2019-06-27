#
#  Copyright (C) 2010, 2016  Smithsonian Astrophysical Observatory
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

from sherpa.stats import Cash, CStat
from sherpa.optmethods import NelderMead
from sherpa.estmethods import Covariance
from sherpa.utils import parallel_map, poisson_noise, NoNewAttributesAfterInit
from sherpa.fit import Fit
from sherpa.sim.sample import NormalParameterSampleFromScaleMatrix

import time
import numpy

from copy import deepcopy

import logging
logger = logging.getLogger("sherpa")
debug = logger.debug
info = logger.info

_tol = numpy.finfo(numpy.float).eps

__all__ = ('LikelihoodRatioTest', 'LikelihoodRatioResults')


class LikelihoodRatioResults(NoNewAttributesAfterInit):
    """The results of a likelihood ratio comparison simulation.

    Attributes
    ----------
    samples : numpy array
       The parameter samples array for each simulation.
    stats : numpy array
       The fit statistic for the null and alternative models
       for each simulation. The shape is (nsim, 2).
    ratios : numpy array
       The likelihood ratio for each simulation.
    null : number
       The fit statistic of the null model on the observed data.
    alt : number
       The fit statistic of the alternate model on the observed data.
    lr : number
       The likelihood ratio of the observed data for the null and
       alternate models.
    ppp : number
       The p value of the observed data for the null and alternate
       models.

    """

    _fields = ('ratios', 'stats', 'samples', 'lr', 'ppp', 'null', 'alt')

    def __init__(self, ratios, stats, samples, lr, ppp, null, alt):
        self.ratios = numpy.asarray(ratios)
        self.stats = numpy.asarray(stats)
        self.samples = numpy.asarray(samples)
        self.lr = float(lr)
        self.ppp = float(ppp)
        self.null = float(null)
        self.alt = float(alt)
        NoNewAttributesAfterInit.__init__(self)


    def __repr__(self):
        return '<Likelihood ratio results instance>'


    def __str__(self):
        samples = self.samples
        if self.samples is not None:
            samples = numpy.array2string(self.samples, separator=',', precision=4, suppress_small=False)

        stats = self.stats
        if self.stats is not None:
            stats = numpy.array2string(self.stats, separator=',', precision=4, suppress_small=False)

        ratios = self.ratios
        if self.ratios is not None:
            ratios = numpy.array2string(self.ratios, separator=',', precision=4, suppress_small=False)

        output = '\n'.join([
                'samples = %s' % samples,
                'stats   = %s' % stats,
                'ratios  = %s' % ratios,
                'null    = %s' % repr(self.null),
                'alt     = %s' % repr(self.alt),
                'lr      = %s' % repr(self.lr),
                'ppp     = %s' % repr(self.ppp)
                ])

        return output


    def format(self):
        """Convert the object to a string representation for display purposes.

        Returns
        -------
        txt : str
           A string representation of the data stored in the object.

        """
        s =  'Likelihood Ratio Test\n'
        s += 'null statistic   =  %s\n' % str(self.null)
        s += 'alt statistic    =  %s\n' % str(self.alt)
        s += 'likelihood ratio =  %s\n' % str(self.lr)
        if self.ppp == 0.0:
            s += 'p-value          <  %s' % str(1./len(self.samples))
        else:
            s += 'p-value          =  %s' % str(self.ppp)
        return s


class LikelihoodRatioTestWorker(object):
    """
    Worker class for LikelihoodRatioTest
    """
    def __init__(self, null_fit, alt_fit, null_vals, alt_vals):
        self.null_fit = null_fit
        self.alt_fit = alt_fit
        self.null_vals = null_vals
        self.alt_vals = alt_vals

    def __call__(self, proposal):
        return LikelihoodRatioTest.calculate(self.null_fit, self.alt_fit, proposal, self.null_vals, self.alt_vals)


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

    References
    ----------

    .. [1] http://en.wikipedia.org/wiki/Likelihood-ratio_test

    """

    @staticmethod
    def calculate(nullfit, altfit, proposal, null_vals, alt_vals):

        # FIXME: only null perturbed?
        nullfit.model.thawedpars = proposal

        # Fake using poisson_noise with null
        fake = poisson_noise(nullfit.data.eval_model(nullfit.model))

                # Set faked data for both nullfit and altfit
        nullfit.data.set_dep(fake)

        # Start the faked fit at initial null best-fit values
        #nullfit.model.thawedpars = null_vals

        # Fit with null model
        nullfr = nullfit.fit()
        debug(nullfr.format())

        null_stat = nullfr.statval
        debug("statistic null = " + repr(null_stat))

        # nullfit and altfit BOTH point to same faked dataset
        assert ( id(nullfit.data) == id(altfit.data) )
        assert (nullfit.data.get_dep() == altfit.data.get_dep()).all()

        # Start the faked fit at the initial alt best-fit values
        #altfit.model.thawedpars = alt_vals

        debug("proposal: " + repr(proposal))
        debug("alt model")
        debug(str(altfit.model))

        # Set alt model and fit
        altfr = altfit.fit()
        debug(altfr.format())

        debug(str(altfit.model))

        alt_stat = altfr.statval
        debug("statistic alt = " + repr(alt_stat))

        LR = -(alt_stat - null_stat)
        debug("LR = " + repr(LR))

        return [null_stat, alt_stat, LR]


    @staticmethod
    def run(fit, null_comp, alt_comp, conv_mdl=None,
            stat=None, method=None,
            niter=500, numcores=None):
        if stat is None:   stat = CStat()
        if method is None: method = NelderMead()

        if not isinstance(stat, (Cash, CStat)):
            raise TypeError("Sherpa fit statistic must be Cash or CStat" +
                            " for likelihood ratio test")

        niter = int(niter)

        alt  = alt_comp
        null = null_comp

        oldaltvals = numpy.array(alt.thawedpars)
        oldnullvals = numpy.array(null.thawedpars)

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
        samples = sampler.get_sample(nullfit, None, niter)

        # Fit with alt model, null component starts at null's best fit params.
        altfit = Fit(data, alt, stat, method, Covariance())
        altfit_results = altfit.fit()
        debug(altfit_results.format())

        alt_stat = altfit_results.statval
        alt_vals = altfit_results.parvals

        LR = -(alt_stat - null_stat)

        olddep = data.get_dep(filter=False)
        try:
            statistics = parallel_map(
                LikelihoodRatioTestWorker(nullfit, altfit, null_vals, alt_vals),
                samples,
                numcores
            )
        finally:
            data.set_dep(olddep)
            alt.thawedpars = list(oldaltvals)
            null.thawedpars = list(oldnullvals)

        debug("statistic null = " + repr(null_stat))
        debug("statistic alt = " + repr(alt_stat))
        debug("LR = " + repr(LR))

        statistics = numpy.asarray(statistics)

        pppvalue = numpy.sum( statistics[:,2] > LR ) / (1.0*niter)

        debug('ppp value = '+str(pppvalue))

        return LikelihoodRatioResults(statistics[:,2], statistics[:,0:2],
                                      samples, LR, pppvalue, null_stat,
                                      alt_stat)
