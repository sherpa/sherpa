#
#  Copyright (C) 2011, 2016, 2019, 2020  Smithsonian Astrophysical Observatory
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
pyBLoCXS is a sophisticated Markov chain Monte Carlo (MCMC) based algorithm
designed to carry out Bayesian Low-Count X-ray Spectral (BLoCXS) analysis in the
Sherpa environment. The code is a Python extension to Sherpa that explores
parameter space at a suspected minimum using a predefined Sherpa model to
high-energy X-ray spectral data. pyBLoCXS includes a flexible definition of
priors and allows for variations in the calibration information. It can be used
to compute posterior predictive p-values for the likelihood ratio test (see
Protassov et al., 2002, ApJ, 571, 545). Future versions will allow for the
incorporation of calibration uncertainty (Lee et al., 2011, ApJ, 731, 126).

MCMC is a complex computational technique that requires some sophistication on
the part of its users to ensure that it both converges and explores the
posterior distribution properly. The pyBLoCXS code has been tested with a number
of simple single-component spectral models. It should be used with great care in
more complex settings. Readers interested in Bayesian low-count spectral
analysis should consult van Dyk et al. (2001, ApJ, 548, 224). pyBLoCXS is based
on the methods in van Dyk et al. (2001) but employs a different MCMC sampler
than is described in that article. In particular, pyBLoCXS has two sampling
modules. The first uses a Metropolis-Hastings jumping rule that is a
multivariate t-distribution with user specified degrees of freedom centered on
the best spectral fit and with multivariate scale determined by the Sherpa
function, covar(), applied to the best fit. The second module mixes this
Metropolis Hastings jumping rule with a Metropolis jumping rule centered at the
current draw, also sampling according to a t-distribution with user specified
degrees of freedom and multivariate scale determined by a user specified scalar
multiple of covar() applied to the best fit.

A general description of the MCMC techniques we employ along with their
convergence diagnostics can be found in Appendices A.2 - A.4 of van Dyk et
al. (2001) and in more detail in Chapter 11 of Gelman, Carlin, Stern, and Rubin
(Bayesian Data Analysis, 2nd Edition, 2004, Chapman & Hall/CRC).

http://hea-www.harvard.edu/AstroStat/pyBLoCXS/
"""

# The pyBLoCXS code base is cleanly separable from Sherpa!

import numpy as np
import logging
import math
import inspect

logger = logging.getLogger("sherpa")
info = logger.info
debug = logger.debug
error = logger.error

__all__ = ('LimitError', 'MetropolisMH', 'MH', 'Sampler',
           'Walk', 'dmvt', 'dmvnorm')


class LimitError(Exception):
    pass

class CovarError(Exception):
    pass


def rmvt(mu, sigma, dof):
    """
    Sampling the non-central multivariate Student's t distribution
    using deviates from multivariate normal and chi-squared distributions
    Source: Kshirsagar method taken from function `rmvt` in R package `mvtnorm`.
    http://cran.r-project.org/web/packages/mvtnorm/index.html

    `mu`     the current sample
    `sigma`  covariance matrix
    `dof`    degrees of freedom

    returns a sample from the multivariate t distribution in the shape of `mu`

    """

    if dof < 1:
        raise ValueError("The degrees of freedom must be > 0")

    zero_vec = np.zeros_like(mu)
    q = np.random.chisquare(dof, 1)[0]
    nsample  = np.random.multivariate_normal(zero_vec, sigma)
    proposal = mu + nsample / np.sqrt(q / dof)
    return proposal


def dmvt(x, mu, sigma, dof, log=True, norm=False):
    """

    Probability Density of a multi-variate Student's t distribution
    """

    #if np.min( np.linalg.eigvalsh(sigma))<=0 :
    #    raise ValueError("Error: sigma is not positive definite")
    if np.max( np.abs(sigma-sigma.T))>=1e-9 :
        raise ValueError("Error: sigma is not symmetric")

    p = mu.size

    # log density unnormalized
    val = (-0.5*np.log(np.linalg.det(sigma)) - (dof+p)/2.0*
            np.log( dof + np.dot( x-mu, np.dot(
                    np.linalg.inv(sigma), x-mu ) ) ) )

    # log density normalized
    if norm:
        lgam = math.lgamma
        val += (lgam((dof+p)/2.) - lgam(dof/2.) - (p/2.) *
                np.log(np.pi) + (dof/2.) * np.log(dof))

    # density
    if not log:
        val = np.exp(val)

    return val


def dmvnorm(x, mu, sigma, log=True):
    """

    Probability Density of a multi-variate Normal distribution
    """

    #if np.min( np.linalg.eigvalsh(sigma))<=0 :
    #    raise ValueError("Error: sigma is not positive definite")
    if np.max( np.abs(sigma-sigma.T))>=1e-9 :
        raise ValueError("Error: sigma is not symmetric")

    # log density
    logdens = (-mu.size/2.0*np.log(2*np.pi)-
                1/2.0*np.log( np.linalg.det(sigma) )-1/2.0 *
                np.dot( x-mu, np.dot(np.linalg.inv(sigma), x-mu ) ) )

    if log:
        return logdens

    # density
    dens = np.exp( logdens )
    return dens



# def progress_bar(current, total, tstart, name=None):
#     """simple progress in percent"""

#     if not sys.stdout.isatty():
#         return

#     percent = 0.0
#     if current != 0 and total != 0:
#         percent = current*100./float(total)

#     output = '\r%.1f %%' % percent
#     if name is not None:
#         output = '\r%s: %.1f %%' % (name, percent)

#     sys.stdout.write(output)
#     if current == total:
#         sys.stdout.write(' finished in %s secs \n' % (time.time()-tstart))

#     sys.stdout.flush()



class Walk():

    def __init__(self, sampler=None, niter=1000):
        self._sampler = sampler
        self.niter = int(niter)

    def set_sampler(self, sampler):
        self._sampler = sampler

    def __call__(self, **kwargs):

        if self._sampler is None:
            raise AttributeError("sampler object has not been set, "+
                                 "please use set_sampler()")

        pars, stat = self._sampler.init(**kwargs)

        # setup proposal variables
        npars = len(pars)
        niter = self.niter
        nelem = niter+1

        proposals = np.zeros((nelem,npars), dtype=np.float)
        proposals[0] = pars.copy()

        stats = np.zeros(nelem, dtype=np.float)
        stats[0] = stat

        acceptflag = np.zeros(nelem, dtype=np.bool)

        # Iterations
        # - no burn in at present
        # - the 0th element of the params array is the input value
        # - we loop until all parameters are within the allowable
        #   range; should there be some check to ensure we are not
        #   rejecting a huge number of proposals, which would indicate
        #   that the limits need increasing or very low s/n data?
        #

        #tstart = time.time()

        try:
            for ii in range(niter):

                #progress_bar(ii, niter, tstart, self._sampler.__class__.__name__)

                jump = ii+1

                current_params = proposals[ii]
                current_stat   = stats[ii]

                # Assume proposal is rejected by default
                proposals[jump] = current_params
                stats[jump]  = current_stat
                #acceptflag[jump] = False

                # Draw a proposal

                try:
                    proposed_params = self._sampler.draw(current_params)
                except CovarError:
                    error("Covariance matrix failed! " + str(proposed_params))
                    # automatically reject if the covar is malformed
                    self._sampler.reject()
                    continue

                proposed_params = np.asarray(proposed_params)
                try:
                    proposed_stat = self._sampler.calc_stat(proposed_params)
                except LimitError:
                    # automatically reject the proposal if outside hard limits
                    self._sampler.reject()
                    continue

                # Accept this proposal?
                if self._sampler.accept(current_params, current_stat,
                                         proposed_params, proposed_stat):
                    proposals[jump] = proposed_params
                    stats[jump] = proposed_stat
                    acceptflag[jump] = True

                else:
                    self._sampler.reject()
        finally:
            self._sampler.tear_down()
            #progress_bar(niter, niter, tstart, self._sampler.__class__.__name__)

        params = proposals.transpose()
        return (stats, acceptflag, params)


class Sampler():

    def __init__(self):

        # get the initial keyword argument defaults;
        # it looks like inspect.getargspec is not being removed
        # in Python 3.6, but use the signature function if it is
        # available to avoid possible warnings.
        try:
            sig = inspect.signature(self.init)
            opts = [(p.name, p.default)
                    for p in sig.parameters.values()
                    if p.kind == p.POSITIONAL_OR_KEYWORD and
                    p.default != p.empty]

        except AttributeError:
            argspec = inspect.getargspec(self.init)
            first = len(argspec[0]) - len(argspec[3])
            opts = zip(argspec[0][first:], argspec[3][0:])

        self._opts = dict(opts)
        self.walk = None

    def init(self):
        raise NotImplementedError

    def draw(self, current, **kwargs):
        raise NotImplementedError

    def accept(self, current, current_stat, proposal, proposal_stat, **kwargs):
        raise NotImplementedError

    def reject(self):
        raise NotImplementedError

    def calc_stat(self, proposed_params):
        raise NotImplementedError

    def tear_down(self):
        raise NotImplementedError


class MH(Sampler):
    """ The Metropolis Hastings Sampler """

    def __init__(self, fcn, sigma, mu, dof, *args):
        self.fcn = fcn
        self._dof = dof
        self._mu = np.array(mu)
        self._sigma = np.array(sigma)

        self.accept_func = None
        self.currently_metropolis = False
        self.prior = None

        # MH tunable parameters
        self.log = False
        self.inv = False
        self.defaultprior = True
        self.priorshape = False
        self.originalscale = True
        self.scale = 1
        self.prior_funcs = ()
        self.sigma_m = False
        Sampler.__init__(self)


    def calc_fit_stat(self, proposed_params):
        return self.fcn(proposed_params)


    def init(self, log=False, inv=False, defaultprior=True, priorshape=False,
             priors=(), originalscale=True, scale=1, sigma_m=False):

        if self._sigma is None or self._mu is None:
            raise AttributeError('sigma or mu is None, initialization failed')

        self.prior = np.ones(self._mu.size)
        self.defaultprior = defaultprior
        self.priorshape = np.array(priorshape)
        self.originalscale = np.array(originalscale)

        self.scale = scale
        self.prior_funcs = priors

        debug(str(self.prior_funcs))

        # if not default prior, prior calculated at each iteration
        if not defaultprior:
            if self.priorshape.size != self._mu.size:
                raise ValueError(
                    "If not using default prior, must specify a " +
                    "function for the prior on each parameter")
            if self.originalscale.size != self._mu.size:
                raise ValueError(
                    "If not using default prior, must specify the " +
                    "scale on which the prior is defined for each parameter")

        self.jacobian = np.zeros(self._mu.size, dtype=bool)
        # jacobian needed if transforming parameter but prior for parameter
        # on original scale
        if not defaultprior:
            # if log transformed but prior on original scale, jacobian
            # for those parameters is needed
            if np.sum( log*self.originalscale ) > 0:
                self.jacobian[ log*self.originalscale ] = True
            if np.sum( inv*self.originalscale ) > 0:
                self.jacobian[ inv*self.originalscale ] = True

        self.log = np.array(log)
        if self.log.size == 1:
            self.log = np.tile(self.log, self._mu.size)

        self.inv = np.array(inv)
        if self.inv.size == 1:
            self.inv = np.tile(self.inv, self._mu.size)

        if np.sum(log*inv) > 0:
            raise TypeError(
                "Cannot specify both log and inv transformation for the same " +
                "parameter")

        debug("Running Metropolis-Hastings")

        current = self._mu.copy()
        stat = self.calc_fit_stat(current)

        # include prior
        stat = self.update(stat, self._mu)

        self.initial_stat = stat

        # using delta method to create proposal distribution on log scale for
        # selected parameters
        if np.sum(self.log) > 0:
            logcovar = self._sigma.copy()
            logcovar[:,self.log]= logcovar[:,self.log]/self._mu[self.log]
            logcovar[self.log]= (logcovar[self.log].T/self._mu[self.log]).T
            self._sigma = np.copy(logcovar)
            self._mu[self.log]=np.log(self._mu[self.log])
            current[self.log]=np.log( current[self.log])

        # using delta method to create proposal distribution on inverse scale
        # for selected parameters
        if np.sum(self.inv) > 0:
            invcovar = self._sigma.copy()
            invcovar[:,self.inv] = invcovar[:,self.inv]/(
                                   -1.0*np.power(self._mu[self.inv],2))
            invcovar[self.inv] = (invcovar[self.inv].T/(
                                  -1.0*np.power(self._mu[self.inv],2))).T
            self._sigma = np.copy(invcovar)
            self._mu[self.inv]=1.0/(self._mu[self.inv])
            current[self.inv]=1.0/( current[self.inv])

        self.rejections=0

        self.sigma_m = sigma_m
        if np.mean(sigma_m) == False:
            self.sigma_m = self._sigma.copy()

        return (current, stat)


    def update(self, stat, mu, init=True):
        """ include prior """
        if not self.defaultprior:
            x = mu.copy()
            if np.sum(self.originalscale) < mu.size:
                for j in range(mu.size):
                    if self.log[j]*(1-self.originalscale[j])>0:
                        x[j] = np.log(x[j])
                    if self.inv[j]*(1-self.originalscale[j])>0:
                        x[j] = 1.0 / x[j]

            for ii, func in enumerate(self.prior_funcs):
                if self.priorshape[ii]:
                    self.prior[ii] = func(x[ii])

        # If no prior then
        # 0.0 == np.sum(np.log(np.ones(mu.size)))
        stat += np.sum(np.log(self.prior))

        if np.sum(self.log*self.jacobian) > 0:
            stat += np.sum( np.log( mu[self.log*self.jacobian] ) )
        if np.sum(self.inv*self.jacobian) > 0:
            stat_temp = np.sum(2.0*np.log(np.abs(mu[self.inv*self.jacobian])))
            if init:
                stat += stat_temp
            else:
                stat -= stat_temp
        return stat


    def draw(self, current):
        """Create a new set of parameter values using the t distribution.

        Given the best-guess (mu) and current (current) set of
        parameters, along with the covariance matrix (sigma),
        return a new set of parameters.
        """
        proposal = self.mh(current)
        self.accept_func = self.accept_mh
        return proposal


    def mh(self, current):
        """ MH jumping rule """

        # The current proposal is ignored here.
        # MH jumps from the best-fit parameter values at each iteration
        proposal = rmvt(self._mu, self._sigma, self._dof)
        return proposal


    def dmvt(self, x, log=True, norm=False):
        return dmvt(x, self._mu, self._sigma, self._dof, log, norm)


    def accept_mh(self, current, current_stat, proposal, proposal_stat):
        alpha = np.exp(proposal_stat + self.dmvt(current) -
                       current_stat - self.dmvt(proposal))
        return alpha


    def accept(self, current, current_stat, proposal, proposal_stat, **kwargs):
        """
        Should the proposal be accepted (using the Cash statistic and the
        t distribution)?
        """
        alpha = self.accept_func(current, current_stat, proposal, proposal_stat)
        u = np.random.uniform(0,1,1)
        return u <= alpha


    def reject(self):
        ### added for test
        self.rejections += 1


    def calc_stat(self, proposed_params):

        if np.sum(self.log)>0:
            proposed_params[self.log]=np.exp(proposed_params[self.log])
        if np.sum(self.inv)>0:
            proposed_params[self.inv]=1.0/proposed_params[self.inv]

        proposed_stat = self.calc_fit_stat(proposed_params)

        #putting parameters back on log scale
        if np.sum(self.log)>0:
            proposed_params[self.log] = np.log(proposed_params[self.log])
        #putting parameters back on inverse scale
        if np.sum(self.inv)>0:
            proposed_params[self.inv] = 1.0/proposed_params[self.inv]

        # include prior
        proposed_stat = self.update(proposed_stat, proposed_params, False)

        return proposed_stat


    def tear_down(self):
        pass


class MetropolisMH(MH):
    """ The Metropolis Metropolis-Hastings Sampler """

    def __init__(self, fcn, sigma, mu, dof, *args):
        MH.__init__(self, fcn, sigma, mu, dof, *args)

        # count the p_M
        self.num_mh = 0
        self.num_metropolis = 0


    def init(self, log=False, inv=False, defaultprior=True, priorshape=False,
             priors=(), originalscale=True, scale=1, sigma_m=False, p_M=.5):

        debug("Running Metropolis with Metropolis-Hastings")

        self.p_M = p_M

        debug("X ~ uniform(0,1) <= %.2f --> Metropolis" % float(p_M))
        debug("X ~ uniform(0,1) >  %.2f --> Metropolis-Hastings" % float(p_M))

        return MH.init(self, log, inv, defaultprior, priorshape, priors,
                       originalscale, scale, sigma_m)


    def draw(self, current):
        """Create a new set of parameter values using the t distribution.

        Given the best-guess (mu) and current (current) set of
        parameters, along with the covariance matrix (sigma),
        return a new set of parameters.
        """
        u = np.random.uniform(0,1,1)
        proposal = None
        if u <= self.p_M:
            proposal = self.metropolis(current)
            self.accept_func = self.accept_metropolis
            self.num_metropolis += 1
        else:
            proposal = self.mh(current)
            self.accept_func = self.accept_mh
            self.num_mh += 1

        return proposal


    def metropolis(self, current):
        """ Metropolis Jumping Rule """

        # Metropolis with MH jumps from the current accepted parameter
        # proposal at each iteration
        proposal = rmvt(current, self.sigma_m*self.scale, self._dof)
        return proposal


    def accept_metropolis(self, current, current_stat, proposal, proposal_stat):
        alpha = np.exp( proposal_stat - current_stat)
        return alpha


    def tear_down(self):
        num = float(self.num_metropolis + self.num_mh)
        if num > 0:
            debug("p_M: %g, Metropolis: %g%%" % (self.p_M, 100 * self.num_metropolis/num))
            debug("p_M: %g, Metropolis-Hastings: %g%%" % (self.p_M, 100 * self.num_mh/num))
