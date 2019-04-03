#
#  Copyright (C) 2011, 2015, 2016, 2018  Smithsonian Astrophysical Observatory
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

"""Monte-Carlo Markov Chain support for low-count data (Poisson statistics).

The ``sherpa.sim`` module provides support for exploring the posterior
probability density of parameters in a fit to low-count data, for
which Poisson statistics hold, using a Bayesian algorithm and a
Monte-Carlo Markov Chain (MCMC). It was originally known as the
pyBLoCXS (python Bayesian Low-Count X-ray Spectral) package [1]_, but
has since been incorporated into Sherpa.

The Sherpa UI modules - e.g. `sherpa.ui` and `sherpa.astro.ui` - provide
many of the routines described below (e.g. ``list_samplers``).

Acknowledgements
----------------

The original version of the code was developed by the CHASC
Astro-Statistics collaboration http://hea-www.harvard.edu/AstroStat/,
and was called pyBLoCXS. It has since been developed by the
Chandra X-ray Center and weas added to Sherpa in version 4.5.1.

Overview
--------

The algorithm explores parameter space at a suspected minimum -
i.e. after a standard Sherpa fit. It supports a flexible definition of priors
and allows for variations in the calibration information. It can be used to
compute posterior predictive p-values for the likelihood ratio test
[2]_. Future versions will allow for the incorporation of calibration
uncertainty [3]_.

MCMC is a complex computational technique that requires some sophistication on
the part of its users to ensure that it both converges and explores the
posterior distribution properly. The pyBLoCXS code has been tested with a
number of simple single-component spectral models. It should be used with
great care in more complex settings. The code is based on the methods in
[4]_ but employs a different MCMC sampler than is described in that article.
A general description of the techniques employed along with their
convergence diagnostics can be found in the Appendices of [4]_
and in [5]_.

Jumping Rules
-------------

The jumping rule determines how each step in the Monte-Carlo Markov
Chain is calculated. The setting can be changed using ``set_sampler``.
The ``sherpa.sim`` module provides the following rules, which may
be augmented by other modules:

- ``MH`` uses a Metropolis-Hastings jumping rule that is a multivariate
  t-distribution with user-specified degrees of freedom centered on the
  best-fit parameters, and with multivariate scale determined by the
  ``covar`` function applied at the best-fit location.

- ``MetropolisMH`` mixes this Metropolis Hastings jumping rule with a
  Metropolis jumping rule centered at the current draw, in both cases
  drawing from the same t-distribution as used with ``MH``. The
  probability of using the best-fit location as the start of the jump
  is given by the ``p_M`` parameter of the rule (use ``get_sampler`` or
  ``get_sampler_opt`` to view and ``set_sampler_opt`` to set this value),
  otherwise the jump is from the previous location in the chain.

Options for the sampler are retrieved and set by ``get_sampler`` or
``get_sampler_opt``, and ``set_sampler_opt`` respectively. The list of
available samplers is given by ``list_samplers``.

Choosing the parameter values
-----------------------------

By default, the prior on each parameter is taken to be flat, varying
from the parameter minima to maxima values. This prior can be changed
using the ``set_prior`` function, which can set the prior for a
parameter to a function or Sherpa model. The list of currently set
prior-parameter pairs is returned by the ``list_priors`` function, and the
prior function associated with a particular Sherpa model parameter may be
accessed with ``get_prior``.

Running the chain
-----------------

The ``get_draws`` function runs a pyBLoCXS chain using fit information
associated with the specified data set(s), and the currently set sampler and
parameter priors, for a specified number of iterations. It returns an array of
statistic values, an array of acceptance Booleans, and a 2-D array of
associated parameter values.

Analyzing the results
---------------------

The module contains several routines to visualize the results of the chain,
including ``plot_trace``, ``plot_cdf``, and ``plot_pdf``, along with
``sherpa.utils.get_error_estimates`` for calculating the limits from a
parameter chain.

References
----------

.. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/

.. [2] "Statistics, Handle with Care: Detecting Multiple Model Components
       with the Likelihood Ratio Test", Protassov et al., 2002, ApJ, 571, 545
       http://adsabs.harvard.edu/abs/2002ApJ...571..545P

.. [3] "Accounting for Calibration Uncertainties in X-ray Analysis:
       Effective Areas in Spectral Fitting", Lee et al., 2011, ApJ, 731, 126
       http://adsabs.harvard.edu/abs/2011ApJ...731..126L

.. [4] "Analysis of Energy Spectra with Low Photon Counts via Bayesian
       Posterior Simulation", van Dyk et al. 2001, ApJ, 548, 224
       http://adsabs.harvard.edu/abs/2001ApJ...548..224V

.. [5] Chapter 11 of Gelman, Carlin, Stern, and Rubin
       (Bayesian Data Analysis, 2nd Edition, 2004, Chapman & Hall/CRC).

Examples
--------

Analysis proceeds as normal, up to the point that a good fit has
been determined, as shown below (note that a Poisson likelihood,
such as the ``cash`` statistic, must be used)::

    >>> from sherpa.astro import ui
    >>> ui.load_pha('src.pi')
    >>> ui.notice(0.5, 7)
    >>> ui.set_source(ui.xsphabs.gal * ui.xspowerlaw.pl)
    >>> ui.set_stat('cash')
    >>> ui.set_method('simplex')
    >>> ui.fit()
    >>> ui.covar()

Once the best-fit location has been determined (which may require
multiple calls to ``fit``), the chain can be run. In this example
the default sampler (``MetropolisMH``) and default parameter priors
(flat, varying between the minimum and maximum values) are used,
as well as the default number of iterations (1000)::

    >>> stats, accept, params = ui.get_draws()

The ``stats`` array contains the fit statistic for each iteration
(the first element of these arrays is the starting point of the chain,
so there will be 1001 elements in this case). The "trace" - i.e.
statistic versus iteration - can be plotted using::

    >>> ui.plot_trace(stats, name='stat')

The ``accept`` array indicates whether, at each iteration, the proposed
jump was accepted, (``True``) or if the previous iterations parameter
values are used. This can be used to look at the acceptance rate for
the chain (dropping the last element and a burn-in period, which
here is arbitrarily taken to be 100)::

    >>> nburn = 100
    >>> arate = accept[nburn:-1].sum() * 1.0 / (len(accept) - nburn - 1)
    >>> print("acceptance rate = {}".format(arate))

The trace of the parameter values can also be displayed; in this
example a burn-in period has not been removed)::

    >>> par1 = params[:, 0]
    >>> par2 = params[:, 1]
    >>> ui.plot_trace(par1, name='par1')
    >>> ui.plot_trace(par2, name='par2')

The cumulative distribution can also be viewed::

    >>> ui.plot_cdf(par1[nburn:], name='par1')

as well as the probability density::

    >>> ui.plot_[pdf(par2[nburn:], name='par2')

The traces can be used to estimate the credible interval for a
parameter::

    >>> from sherpa.utils import get_error_estimates
    >>> pval, plo, phi = get_error_estimates(par1[nburn:])

"""

from six import string_types

# Although all this module needs is the following import
#   from sherpa.sim.mh import LimitError, MetropolisMH, MH, Sampler, Walk
# it looks like the following modules are being re-exported by this
# one, so the 'from blah import *' lines can not easily be removed.
#
from sherpa.sim.simulate import *
from sherpa.sim.sample import *
from sherpa.sim.mh import *
from sherpa.utils import NoNewAttributesAfterInit, get_keyword_defaults, \
    sao_fcmp
from sherpa.stats import Cash, CStat, WStat, LeastSq

from sherpa.fit import Fit
from sherpa.data import Data1D, Data1DAsymmetricErrs
from sherpa.optmethods import LevMar

import numpy
import logging
info = logging.getLogger("sherpa").info
_log = logging.getLogger("sherpa")

_tol = numpy.finfo(numpy.float).eps


def flat(x):
    """The flat prior (returns 1 everywhere)."""

    return 1.0


def inverse(x):
    """Returns the inverse of x."""

    prior = 1.0 / x
    return prior


def inverse2(x):
    """Returns the invers of x^2."""

    prior = 1.0 / (x * x)
    return prior


_samplers = dict(metropolismh=MetropolisMH, mh=MH)
_walkers = dict(metropolismh=Walk, mh=Walk)


class MCMC(NoNewAttributesAfterInit):
    """

    High-level UI to pyBLoCXS that joins the loop in 'Walk' with the jumping
    rule in 'Sampler'.  Implements a user interface for configuration.  This
    class implements a calc_stat() function using the Sherpa interface to
    'Fit'.

    """
    __samplers = _samplers.copy()

    __walkers = _walkers.copy()

    def _get_sampler(self):
        return self._sampler

    def _set_sampler(self, sampler):
        self._sampler = sampler
        self._sampler_opt = get_keyword_defaults(sampler.init)

    sampler = property(_get_sampler, _set_sampler)

    def _get_walker(self):
        return self._walker

    def _set_walker(self, walker):
        self._walker = walker

    walker = property(_get_walker, _set_walker)

    def _get_sampler_opt(self, opt):
        return self._sampler_opt[opt]

    def _set_sampler_opt(self, opt, val):
        self._sampler_opt[opt] = val

    def __init__(self):
        self.priors = {}
        self._walker = Walk
        self._sampler = MetropolisMH
        self._sampler_opt = get_keyword_defaults(MetropolisMH.init)
        self.sample = None
        self.walk = lambda: None
        NoNewAttributesAfterInit.__init__(self)

    def __getstate__(self):
        state = self.__dict__.copy()
        del state['walk']
        del state['sample']
        return state

    def __setstate__(self, state):
        self.walk = lambda: None
        self.sample = None
        self.__dict__.update(state)

    # ## DOC-TODO: include examples
    def list_priors(self):
        """Return the priors set for model parameters, if any.

        Returns
        -------
        priors : dict
           The dictionary of mappings between parameters (keys)
           and prior functions (values) created by `set_prior`.

        See Also
        --------
        get_prior : Return the prior function for a parameter.
        set_prior : Set the prior function to use with a parameter.

        """
        return self.priors

    def get_prior(self, par):
        """Return the prior function for a parameter.

        Parameters
        ----------
        par : sherpa.models.parameter.Parameter
           A parameter of a model instance.

        Returns
        -------
        prior :
           The function or parameter instance set by
           a previous call to `set_prior`.

        Raises
        ------
        ValueError
           If a prior has not been set for the parameter.

        See Also
        --------
        set_prior : Set the prior function to use with a parameter.

        Examples
        --------

        >>> pfunc = get_prior(bgnd.c0)

        """
        prior = self.priors.get(par.fullname, None)
        if prior is None:
            raise ValueError("prior function has not been set for '%s'" %
                             par.fullname)
        return prior

    # ## DOC-TODO: should set_sampler_opt be mentioned here?
    def set_prior(self, par, prior):
        """Set the prior function to use with a parameter.

        The default prior used by ``get_draws`` for each parameter
        is flat, varying between the hard minimum and maximum
        values of the parameter (as given by the ``hard_min`` and
        ``hard_max`` attributes of the parameter object). The ``set_prior``
        function is used to change the form of the prior for a
        parameter.

        Parameters
        ----------
        par : sherpa.models.parameter.Parameter instance
           A parameter of a model instance.
        prior : function or sherpa.models.model.Model instance
           The function to use for a prior. It must accept a
           single argument and return a value of the same size
           as the input.

        See Also
        --------
        get_draws : Run the MCMC algorithm.
        get_prior : Set the prior function to use with a parameter.
        set_sampler : Set the MCMC sampler.

        Examples
        --------

        Set the prior for the ``kT`` parameter of the ``therm`` instance
        to be a gaussian, centered on 1.7 keV and with a FWHM of 0.35
        keV::

        >>> create_model_component('xsapec', 'therm')
        >>> create_model_component('gauss1d', 'p_temp')
        >>> p_temp.pos = 1.7
        >>> p_temp.fwhm = 0.35
        >>> set_prior(therm.kT, p_temp)

        Create a function (``lognorm``) and use it as the prior of the
        ``nH`` parameter of the ``abs1`` instance::

            >>> def lognorm(x):
            ...     nh = 20
            ...     sigma = 0.5  # use a sigma of 0.5
            ...     # nH is in units of 10^-22 so convert
            ...     dx = np.log10(x) + 22 - nh
            ...     norm = sigma / np.sqrt(2 * np.pi)
            ...     return norm * np.exp(-0.5 * dx * dx / (sigma * sigma))
            ...
            >>> create_model_component('xsphabs', 'abs1')
            >>> set_prior(abs1.nH, lognorm)

        """

        # NOTE: the second piece of code is indented in the example
        #       above because otherwise sphinx seems to think that the
        #       colon at the end of the "def lognorm" line ends the
        #       code block.

        self.priors[par.fullname] = prior

    def list_samplers(self):
        """List the samplers available for MCMC analysis with ``get_draws``.

        Returns
        -------
        samplers : list of str
           A list of the names (in lower case) that can be used with
           `set_sampler`.

        See Also
        --------
        get_sampler_name : Return the name of the current MCMC sampler.
        set_sampler : Set the MCMC sampler.

        Notes
        -----
        The available samplers depends on what modules have been
        loaded, and may be more extensive than the example output
        below.

        Examples
        --------

        >>> list_samplers()
        ['metropolismh', 'mh']

        """
        return list(self.__samplers.keys())

    def set_sampler(self, sampler):
        """Set the MCMC sampler.

        The sampler determines the type of jumping rule to
        be used by ``get_draws``.

        Parameters
        ----------
        sampler : str or sherpa.sim.Sampler instance
           When a string, the name of the sampler to use (case
           insensitive). The supported options are given by the
           `list_samplers` function.

        See Also
        --------
        get_draws : Run the MCMC algorithm.
        list_samplers : List the MCMC samplers.
        set_sampler : Set the MCMC sampler.
        set_sampler_opt : Set an option for the current MCMC sampler.

        Notes
        -----
        The jumping rules are:

        MH
           The Metropolis-Hastings rule, which always jumps from the
           best-fit location, even if the previous iteration had moved
           away from it.

        MetropolisMH
           This is the Metropolis with Metropolis-Hastings algorithm,
           that jumps from the best-fit with probability ``p_M``,
           otherwise it jumps from the last accepted jump. The
           value of ``p_M`` can be changed using `set_sampler_opt`.

        PragBayes
           This is used when the effective area calibration
           uncertainty is to be included in the calculation. At each
           nominal MCMC iteration, a new calibration product is
           generated, and a series of N (the ``nsubiters`` option) MCMC
           sub-iteration steps are carried out, choosing between
           Metropolis and Metropolis-Hastings types of samplers with
           probability ``p_M``.  Only the last of these sub-iterations
           are kept in the chain.  The ``nsubiters`` and ``p_M`` values
           can be changed using `set_sampler_opt`.

        FullBayes
           Another sampler for use when including uncertainties due
           to the effective area.

        Examples
        --------

        >>> set_sampler('metropolismh')

        """
        if isinstance(sampler, string_types):

            # case insensitive
            sampler = str(sampler).lower()

            if sampler not in self.__samplers:
                raise TypeError("Unknown sampler '%s'" % sampler)

            self.sampler = self.__samplers.get(sampler)
            self.walker = self.__walkers.get(sampler, Walk)

        elif issubclass(sampler, Sampler):
            self.sampler = sampler
            self.walker = self.__walkers.get(sampler, Walk)

        else:
            raise TypeError("Unknown sampler '%s'" % sampler)

    def get_sampler(self):
        """Return the current MCMC sampler options.

        Returns
        -------
        options : dict
           A copy of the options for the chosen sampler.  Use
           `set_sampler_opt` to change these values. The fields depend
           on the current sampler.

        See Also
        --------
        get_sampler_name : Return the name of the current MCMC sampler.
        get_sampler_opt : Return an option of the current MCMC sampler.
        set_sampler : Set the MCMC sampler.
        set_sampler_opt : Set an option for the current MCMC sampler.

        """
        return self._sampler_opt.copy()

    def get_sampler_name(self):
        """Return the name of the current MCMC sampler.

        Returns
        -------
        name : str

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_sampler : Set the MCMC sampler.

        Examples
        --------

        >>> get_sampler_name()
        'MetropolisMH'

        """
        return self.sampler.__name__

    def get_sampler_opt(self, opt):
        """Return an option of the current MCMC sampler.

        Returns
        -------
        opt : str
           The name of the option. The fields depend on the current
           sampler.

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_sampler_opt : Set an option for the current MCMC sampler.

        Examples
        --------

        >>> get_sampler_opt('log')
        False

        """
        return self._get_sampler_opt(opt)

    def set_sampler_opt(self, opt, value):
        """Set an option for the current MCMC sampler.

        Parameters
        ----------
        opt : str
           The option to change. Use `get_sampler` to view the
           available options for the current sampler.
        value
           The value for the option.

        See Also
        --------
        get_sampler : Return the current MCMC sampler options.
        set_prior: Set the prior function to use with a parameter.
        set_sampler : Set the MCMC sampler.

        Notes
        -----
        The options depend on the sampler. The options include:

        defaultprior
           Set to ``False`` when the default prior (flat, between the
           parameter's soft limits) should not be used. Use
           `set_prior` to set the form of the prior for each
           parameter.

        inv
           A bool, or array of bools, to indicate which parameter is
           on the inverse scale.

        log
           A bool, or array of bools, to indicate which parameter is
           on the logarithm (natural log) scale.

        original
           A bool, or array of bools, to indicate which parameter is
           on the original scale.

        p_M
           The proportion of jumps generatd by the Metropolis
           jumping rule.

        priorshape
           An array of bools indicating which parameters have a
           user-defined prior functions set with `set_prior`.

        scale
           Multiply the output of `covar` by this factor and
           use the result as the scale of the t-distribution.

        Examples
        --------

        >>> set_sampler_opt('scale', 3)


        """
        self._set_sampler_opt(opt, value)

    def get_draws(self, fit, sigma, niter=1000):
        """Run the pyBLoCXS MCMC algorithm.

        The function runs a Markov Chain Monte Carlo (MCMC) algorithm
        designed to carry out Bayesian Low-Count X-ray Spectral
        (BLoCXS) analysis. Unlike many MCMC algorithms, it is
        designed to explore the parameter space at the suspected
        statistic minimum (i.e.  after using `fit`). The return
        values include the statistic value, parameter values, and a
        flag indicating whether the row represents a jump from the
        current location or not.

        Parameters
        ----------
        fit
           The Sherpa fit object to use.
        sigma
           The covariance matrix, defined at the best-fit parameter
           values.
        niter : int, optional
           The number of draws to use. The default is ``1000``.

        Returns
        -------
        stats, accept, params
           The results of the MCMC chain. The stats and accept arrays
           contain niter+1 elements, with the first row being the
           starting values. The params array has (nparams,niter+1)
           elements, where nparams is the number of free parameters in
           the model expression, and the first column contains the
           values that the chain starts at. The accept array contains
           boolean values, indicating whether the jump, or step, was
           accepted (``True``), so the parameter values and statistic
           change, or it wasn't, in which case there is no change to
           the previous row.

        """
        if not isinstance(fit.stat, (Cash, CStat, WStat)):
            raise ValueError("Fit statistic must be cash, cstat or " +
                             "wstat, not %s" % fit.stat.name)

        _level = _log.getEffectiveLevel()
        mu = fit.model.thawedpars
        dof = len(mu)

        info('Using Priors:')
        priors = []
        for par in fit.model.pars:
            if not par.frozen:
                name = par.fullname
                # assume all parameters have flat priors
                func = flat
                if name in self.priors:
                    # update the parameter priors with user defined values
                    func = self.priors[name]
                priors.append(func)
                info(": ".join([name, str(func)]))

        sampler = self._sampler
        walker = self._walker
        sampler_kwargs = self._sampler_opt.copy()
        sampler_kwargs['priors'] = priors

        oldthawedpars  = numpy.array(mu)
        thawedparmins  = fit.model.thawedparhardmins
        thawedparmaxes = fit.model.thawedparhardmaxes

        def calc_stat(proposed_params):

            # automatic rejection outside hard limits
            mins  = sao_fcmp(proposed_params, thawedparmins, _tol)
            maxes = sao_fcmp(thawedparmaxes, proposed_params, _tol)
            if -1 in mins or -1 in maxes:
                raise LimitError('Sherpa parameter hard limit exception')

            level = _log.getEffectiveLevel()

            try:
                # ignore warning from Sherpa about hard limits
                _log.setLevel(logging.CRITICAL)

                # soft limits are ignored, hard limits rejected.
                # proposed values beyond hard limit default to limit.
                fit.model.thawedpars = proposed_params

                # Calculate statistic on proposal, use likelihood
                proposed_stat = -0.5 * fit.calc_stat()

                # _log.setLevel(level)

            except:
                # set the model back to original state on exception
                fit.model.thawedpars = oldthawedpars
                raise
            finally:
                # set the logger back to previous level
                _log.setLevel(level)

            return proposed_stat

        try:
            fit.model.startup()
            self.sample = sampler(calc_stat, sigma, mu, dof, fit)
            self.walk = walker(self.sample, niter)
            stats, accept, params = self.walk(**sampler_kwargs)
        finally:
            fit.model.teardown()

            # set the model back to original state
            fit.model.thawedpars = oldthawedpars

            # set the logger back to previous level
            _log.setLevel(_level)

        # Change to Sherpa statistic convention
        stats = -2.0 * stats

        return (stats, accept, params)


class ReSampleData(NoNewAttributesAfterInit):

    def __init__(self, data, model):
        self.data = data
        self.model = model
        NoNewAttributesAfterInit.__init__(self)
        return

    def __call__(self, niter=1000, seed=123):

        pars = {}
        pars_index = {}
        index = 0
        for par in self.model.pars:
            if par.frozen is False:
                pars[par.name] = []
                pars_index[index] = par.name
                index += 1

        data = self.data
        y = data.y
        x = data.x
        if type(data) == Data1DAsymmetricErrs:
            y_l = y - data.elo
            y_h = y + data.ehi
        elif isinstance(data, (Data1D,)):
            y_l = data.staterror
            y_h = data.staterror
        else:
            msg ="{0} {1}".format(ReSampleData.__name__, type(data))
            raise NotImplementedError(msg)
        
        numpy.random.seed(seed)
        for j in range(niter):
            ry = []
            for i in range(len(y_l)): 
                a = y_l[i]
                b = y_h[i]
                r = -1
                while r < a or r > b:
                    sigma = b - y[i]
                    u = numpy.random.random_sample()
                    if u < 0.5:
                        sigma=y[i]-a
                    r = numpy.random.normal(loc=y[i],scale=sigma,size=None)
                    if u < 0.5 and r > y[i]:
                        r = -1
                    if u > 0.5 and r < y[i]:
                        r = -1
                ry.append(r)

            # fit is performed for each simulated data point
            fit = Fit(Data1D('tmp', x, ry), self.model, LeastSq( ), LevMar())
            fit_result = fit.fit()

            for index, val in enumerate(fit_result.parvals):
                name = pars_index[index]
                pars[name].append(val)
            
        result = []
        for index, name in pars_index.items():
            result.append(numpy.average(pars[name]))
            result.append(numpy.std(pars[name]))

        print(result)
        return result
