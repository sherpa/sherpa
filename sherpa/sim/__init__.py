# 
#  Copyright (C) 2011, 2015  Smithsonian Astrophysical Observatory
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


from sherpa.sim.simulate import *
from sherpa.sim.sample import *
from sherpa.sim.mh import *

from sherpa.utils import NoNewAttributesAfterInit, get_keyword_defaults, \
    sao_fcmp
from sherpa.stats import Cash, CStat

import numpy
import logging
info = logging.getLogger("sherpa").info
_log = logging.getLogger("sherpa")
_level = _log.getEffectiveLevel()

_tol = numpy.finfo(numpy.float).eps


def flat(x):
    return 1.0

def inverse(x):
    prior = 1.0/x
    return prior

def inverse2(x):
    prior = 1.0/(x*x)
    return prior

_samplers = dict(metropolismh=MetropolisMH, mh=MH)
_walkers = dict(metropolismh=Walk, mh=Walk)


class MCMC(NoNewAttributesAfterInit):
    """

    High-level UI to pyBLoCXS that joins the loop in 'Walk' with the jumping
    rule in 'Sampler'.  Implements a user interface for configuration.  This
    class implements a calc_stat() function using the Sherpa interface to 'Fit'.

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
        self.walk = lambda : None
        NoNewAttributesAfterInit.__init__(self)


    def __getstate__(self):
        state = self.__dict__.copy()
        del state['walk']
        del state['sample']
        return state


    def __setstate__(self, state):
        self.walk = lambda : None
        self.sample = None
        self.__dict__.update(state)


    ### DOC-TODO: include examples
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


    ### DOC-TODO: should set_sampler_opt be mentioned here?
    def set_prior(self, par, prior):
        """Set the prior function to use with a parameter.

        The pyBLoCXS Markov Chain Monte Carlo (MCMC) algorithm [1]_
        supports Bayesian Low-Count X-ray Spectral analysis. By
        default, a flat prior is used for each parameter in the fit,
        varying between its soft minimum and maximum values.  The
        `set_prior` function is used to change the form of the prior
        for a parameter.

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
        get_draws : Run the pyBLoCXS MCMC algorithm.
        get_prior : Set the prior function to use with a parameter.
        set_sampler : Set the pyBLoCXS sampler.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        Set the prior for the `kT` parameter of the `therm` instance
        to be a gaussian, centered on 1.7 keV and with a FWHM of 0.35
        keV:

        >>> create_model_component('xsapec', 'therm')
        >>> create_model_component('gauss1d', 'p_temp')
        >>> p_temp.pos = 1.7
        >>> p.temo_fwhm = 0.35
        >>> set_prior(therm.kT, p_temp)

        Create a function (`lognorm`) and use it as the prior the the
        `nH` parameter of the `abs1` instance:

        >>> create_model_component('xsphabs', 'abs1')
        >>> def lognorm(x):
           # center on 10^20 cm^2 with a sigma of 0.5
           sigma = 0.5
           x0 = 20
           # nH is in units of 10^-22 so convert
           dx = np.log10(x) + 22 - x0
           norm = sigma / np.sqrt(2 * np.pi)
           return norm * np.exp(-0.5*dx*dx/(sigma*sigma))

        >>> set_prior(abs1.nH, lognorm)

        """
        self.priors[par.fullname] = prior


    def list_samplers(self):
        """List the pyBLoCXS samplers.

        Returns
        -------
        samplers : list of str
           A list of the names (in lower case) that can be used with
           `set_sampler`.

        See Also
        --------
        get_sampler_name : Return the name of the current pyBLoCXS sampler.
        set_sampler : Set the pyBLoCXS sampler.

        Examples
        --------

        >>> list_samplers()
        ['metropolismh', 'fullbayes', 'mh', 'pragbayes']

        """
        return self.__samplers.keys()


    def set_sampler(self, sampler):
        """Set the pyBLoCXS sampler.

        The sampler determines the type of jumping rule to
        be used when running the MCMC analysis.

        Parameters
        ----------
        sampler : str or sherpa.sim.Sampler instance
           When a string, the name of the sampler to use (case
           insensitive). The supported options are given by the
           `list_samplers` function.

        See Also
        --------
        get_draws : Run the pyBLoCXS MCMC algorithm.
        list_samplers : List the pyBLoCXS samplers.
        set_sampler : Set the pyBLoCXS sampler.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        Notes
        -----
        The jumping rules are [1]_:

        MH
           The 'MH' option refers to the Metropolis-Hastings rule,
           which always jumps from the best-fit location.

        MetropolisMH
           This is the Metropolis with Metropolis-Hastings algorithm,
           that jumps from the best-fit with probability 'p_M',
           otherwise it jumps from the last accepted jump. The
           value of `p_M` can be changed using `set_sampler_opt`.

        PragBayes
           This is used when the effective area calibration
           uncertainty is to be included in the calculation. At each
           nominal MCMC iteration, a new calibration product is
           generated, and a series of N (the `nsubiters` option) MCMC
           sub-iteration steps are carried out, choosing between
           Metropolis and Metropolis-Hastings types of samplers with
           probability `p_M`.  Only the last of these sub-iterations
           are kept in the chain.  The `nsubiters` and `p_M` values
           can be changed using `set_sampler_opt`.

        FullBayes
           Another sampler for use when including uncertainties due
           to the effective area.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        >>> set_sampler('metropolismh')

        """
        if isinstance(sampler, basestring):

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
        """Return the current pyBLoCXS sampler options.

        Returns
        -------
        options : dict
           A copy of the options for the chosen sampler.  Use
           `set_sampler_opt` to change these values. The fields depend
           on the current sampler.

        See Also
        --------
        get_sampler_name : Return the name of the current pyBLoCXS sampler.
        get_sampler_opt : Return an option of the current pyBLoCXS sampler.
        set_sampler : Set the pyBLoCXS sampler.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        """
        #return self.sampler
        return self._sampler_opt.copy()


    def get_sampler_name(self):
        """Return the name of the current pyBLoCXS sampler.

        Returns
        -------
        name : str

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_sampler : Set the pyBLoCXS sampler.

        Examples
        --------

        >>> get_sampler_name()
        'MetropolisMH'

        """
        return self.sampler.__name__


    def get_sampler_opt(self, opt):
        """Return an option of the current pyBLoCXS sampler.

        Returns
        -------
        opt : str
           The name of the option. The fields depend on the current
           sampler.

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_sampler_opt : Set an option for the current pyBLoCXS sampler.

        Examples
        --------

        >>> get_sampler_opt('log')
        False

        """
        return self._get_sampler_opt(opt)


    def set_sampler_opt(self, opt, value):
        """Set an option for the current pyBLoCXS sampler.

        Parameters
        ----------
        opt : str
           The option to change. Use `get_sampler` to view the
           available options for the current sampler.
        value :
           The value for the option.

        See Also
        --------
        get_sampler : Return the current pyBLoCXS sampler options.
        set_prior: Set the prior function to use with a parameter.
        set_sampler : Set the pyBLoCXS sampler.

        Notes
        -----
        The options depend on the sampler [1]_. The options include:

        `defaultprior`
           Set to `False` when the default prior (flat, between the
           parameter's soft limits) should not be used. Use
           `set_prior` to set the form of the prior for each
           parameter.

        `inv`
           A bool, or array of bools, to indicate which parameter is
           on the inverse scale.

        `log`
           A bool, or array of bools, to indicate which parameter is
           on the logarithm (natural log) scale.

        `original`
           A bool, or array of bools, to indicate which parameter is
           on the original scale.

        `p_M`
           The proportion of jumps generatd by the Metropolis
           jumping rule.

        `priorshape`
           An array of bools indicating which parameters have a
           user-defined prior functions set with `set_prior`.

        `scale`
           Multiply the output of `covar` by this factor and
           use the result as the scale of the t-distribution.

        References
        ----------

        .. [1] http://hea-www.harvard.edu/AstroStat/pyBLoCXS/#high-level-user-interface-functions

        Examples
        --------

        >>> set_sampler_opt('scale', 3)


        """
        self._set_sampler_opt(opt, value)




    def get_draws(self, fit, sigma, niter=1000):
        """ 
        Run pyblocxs using current sampler and current sampler configuration
        options for *niter* number of iterations.  The results are returned as a
        3-tuple of Numpy ndarrays.  The tuple specifys an array of statistic
        values, an array of acceptance flags, and a 2-D array of associated
        parameter values.


        `fit`    Sherpa fit object
        `sigma`  Covariance matrix, centered on the best-fit parameter values
        `niter`  Number of iterations, default = 1000

        returns a tuple of ndarrays e.g. (stats, accept, params)

        Examples
        --------

        >>> stats, accept, params = get_draws(fit, niter=1e4)

        """
        if not isinstance(fit.stat, (Cash, CStat)):
            raise ValueError("Fit statistic must be cash or cstat, not %s" %
                             fit.stat.name)

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
                #print'hard limit exception'
                raise LimitError('Sherpa parameter hard limit exception')

            level = _log.getEffectiveLevel()

            try:
                # ignore warning from Sherpa about hard limits
                _log.setLevel(logging.CRITICAL)

                # soft limits are ignored, hard limits rejected.
                # proposed values beyond hard limit default to limit.
                fit.model.thawedpars = proposed_params

                # Calculate statistic on proposal, use likelihood
                proposed_stat = -0.5*fit.calc_stat()

                #_log.setLevel(level)

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
        stats = -2.0*stats

        return (stats, accept, params)
