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


    def list_priors(self):
        """
        List the dictionary of currently set prior functions for the set
        of thawed Sherpa model parameters

        """
        return str(self.priors)


    def get_prior(self, par):
        """
        Get the prior function set for the input Sherpa parameter

        `par`    Sherpa model parameter

        returns associated prior function

        """
        prior = self.priors.get(par.fullname, None)
        if prior is None:
            raise ValueError("prior function has not been set for '%s'" %
                             par.fullname)
        return prior


    def set_prior(self, par, prior):
        """
        Set the prior function for an associated input Sherpa parameter

        `par`    Sherpa model parameter
        `prior`  function pointer to run as prior on associated parameter

        returns None

        Example:

        set_prior(abs1.nh, inverse2)

        """
        self.priors[par.fullname] = prior


    def list_samplers(self):
        """
        List the dictionary of available MCMC Sampler classes
        """       
        return self.__samplers.keys()


    def set_sampler(self, sampler):
        """
        Set the sampler type for use within pyblocxs

        `sampler`   String indicating sampler type, default="MetropolisMH" or
                    class of type Sampler.

        returns None

        Example:

        set_sampler("MH")

        set_sampler(MH)

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
        """
        Return the current sampler's dictionary of configuration options

        returns a dictionary of configuration options

        """
        #return self.sampler
        return self._sampler_opt.copy()


    def get_sampler_name(self):
        """
        Return the current sampler type 

        returns a string indicating the current sampler type
        """
        return self.sampler.__name__


    def get_sampler_opt(self, opt):
        """
        Return the value of an input sampler configuration option

        `opt`     String indicating option

        returns option value

        Example:

        get_sampler_opt("log")
        True

        """
        return self._get_sampler_opt(opt)


    def set_sampler_opt(self, opt, value):
        """
        Set the sampler configuration option for use within pyblocxs

        `opt`     String indicating option
        `value`   Option value

        returns None

        Example:

        set_sampler_opt("log", True)

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

        Example:

        stats, accept, params = get_draws(fit, niter=1e4)

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
