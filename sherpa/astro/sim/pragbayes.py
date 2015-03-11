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


import logging
import numpy as np
import time
from sherpa.fit import Fit
from sherpa.optmethods import NelderMead, LevMar
from sherpa.estmethods import Covariance
from sherpa.sim.mh import MetropolisMH, rmvt, CovarError, Walk, LimitError
read_table_blocks = None
try:
    from sherpa.astro.io import read_table_blocks
except:
    read_table_blocks = None

logger = logging.getLogger("sherpa")
info = logger.info


__all__=['PragBayes', 'PCA1DAdd', 'SIM1DAdd', 'ARFSIMFactory',
         'WalkWithSubIters']


class ARFSIMFactory(object):
    
    def __call__(self, filename):
        return self.read(filename)

    def read(self, filename):
        filename, cols, hdr = read_table_blocks(filename)
        
        emethod = None
        for key in hdr.keys():
            if hdr[key].has_key('EMETHOD'):
                emethod = hdr[key]['EMETHOD'].strip().upper()

        if emethod is not None and emethod.startswith('PCA1DADD'):
            bias      = cols[2]['BIAS']
            component = cols[3]['COMPONENT']
            fvariance = cols[3]['FVARIANCE']
            eigenval  = cols[3]['EIGENVAL']
            eigenvec  = cols[3]['EIGENVEC']
            return PCA1DAdd(bias, component, fvariance, eigenval, eigenvec)

        elif emethod is not None and emethod.startswith('SIM1DADD'):
            bias      = cols[2]['BIAS']
            component = cols[3]['COMPONENT']
            simcomp   = cols[3]['SIMCOMP'] 
            return SIM1DAdd(bias, component, simcomp)

        raise TypeError("Unknown simulation ARF '%s'" % filename)



class PCA1DAdd(object):

    def __init__(self, bias, component, fvariance, eigenval, eigenvec):
        self.bias      = bias
        self.component = component
        self.fvariance = fvariance
        self.eigenval  = eigenval
        self.eigenvec  = eigenvec
        self.ncomp     = len(self.component)
        self.rrout     = None

    def add_deviations(self, specresp, rrin=None, rrsig=None):
        # copy the old ARF (use new memory for deviations)
        new_arf = np.add(specresp, self.bias)

        rrout = np.random.standard_normal(self.ncomp)
        if rrin is not None and rrsig is not None:
            rrout = rrin + rrsig * rrout
        self.rrout = rrout

        tmp = self.eigenvec * self.eigenval[:,np.newaxis] * rrout[:,np.newaxis]
        return np.add(new_arf, tmp.sum(axis=0), new_arf)


class SIM1DAdd(object):

    def __init__(self, bias, component, simcomp):
        self.bias      = bias
        self.component = component
        self.simcomp   = simcomp
        self.ncomp     = len(self.component)


    def add_deviations(self, specresp):
        # copy the old ARF (use new memory for deviations)
        new_arf = np.add(specresp, self.bias)
        # Include the perturbed effective area in each iteration.
        rr = np.random.randint(self.ncomp)
        return np.add(new_arf, self.simcomp[rr], new_arf)


def search_arfs(fit):
    datasets = [fit.data]
    if hasattr(fit.data, 'datasets'):
        datasets = fit.data.datasets

    srcarfs = {}
    bkgarfs = {}

    for ii, data in enumerate(datasets):
        if not hasattr(data, 'response_ids'):
            #raise TypeError("dataset does not contain an ARF, dataset must be PHA")
            continue

        srcarfs[ii] = {}
        for resp_id in data.response_ids:
            arf, rmf = data.get_response(resp_id)
            srcarfs[ii][resp_id] = arf

        bkgarfs[ii] = {}

        # Update only the source ARF for now

        # FIXME: do we include the background ARFs
        # for bkg_id in data.background_ids:
        #     bkg = data.get_background(bkg_id)
        #     bkgarfs[ii][bkg_id] = {}
        #     for bkg_resp_id in bkg.response_ids:
        #         barf, brmf = bkg.get_response(bkg_resp_id)
        #         bkgarfs[ii][bkg_id][bkg_resp_id] = barf

    return srcarfs, bkgarfs


def flatten_arfs(src, bkg):
    arfs = []
    for srcid in src.keys():
        for respid in src[srcid].keys():
            arfs.append(src[srcid][respid])

    # Update only the source ARF for now

    # FIXME: do we include the background ARFs
    # for srcid in bkg.keys():
    #     for bkgid in bkg[srcid].keys():
    #         for respid in bkg[srcid][bkgid].keys():
    #             arfs.append(bkg[srcid][bkgid][respid])

    return arfs



class WalkWithSubIters(Walk):

    def __init__(self, sampler=None, niter=1000):
        self._sampler = sampler
        self.niter = int(niter)
        self.nsubiter = 1

    def set_sampler(self, sampler):
        self._sampler = sampler

    def __call__(self, **kwargs):

        subiters = kwargs.pop('nsubiter', None)
        if subiters is not None:
            self.nsubiter = int(subiters)

        if self._sampler is None:
            raise AttributeError("sampler object has not been set, "+
                                 "please use set_sampler()")

        pars, stat = self._sampler.init(**kwargs)

        # setup proposal variables
        npars = len(pars)
        nsubiter = int(self.nsubiter)
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
            for ii in xrange(niter):
                jump = ii+1

                current_params = proposals[ii]
                current_stat   = stats[ii]

                # Perturb the associated PHA ARF with deviations
                self._sampler.perturb_arf(current_params, current_stat)

                # Initialize the likelihood using perturbed ARF and current
                # set of parameter values
                stats[ii] = self._sampler.calc_stat(proposals[ii])

                # Assume proposal is rejected by default
                proposals[jump] = current_params
                stats[jump]  = current_stat            

                for jj in xrange(nsubiter):

                    #progress_bar(ii*nsubiter+jj, niter*nsubiter, tstart, self._sampler.__class__.__name__)
                    
                    # Draw a proposal
                    try:
                        proposed_params = self._sampler.draw(current_params)
                    except CovarError:
                        info("Draw rejected: covariance matrix failed. " + str(proposed_params))
                        # automatically reject if the covar is malformed
                        self._sampler.reject()
                        continue

                    proposed_params = np.asarray(proposed_params)
                    try:
                        proposed_stat = self._sampler.calc_stat(proposed_params)
                    except LimitError:
                        info("Draw rejected: parameter boundary exception: " + str(proposed_params))
                        # automatically reject the proposal if outside hard limits
                        self._sampler.reject()
                        continue

                    # Accept this proposal?
                    if self._sampler.accept(current_params, current_stat,
                                             proposed_params, proposed_stat):
                        proposals[jump] = np.array(proposed_params)
                        stats[jump] = proposed_stat
                        acceptflag[jump] = True

                        # Mini convergence step
                        current_params = proposed_params
                        current_stat   = proposed_stat

                    else:
                        self._sampler.reject()
                        acceptflag[jump] = False

                        # If acceptance fails do we start back at the last accepted 
                        # iteration or the last accepted subiteration?
                        current_params = proposals[ii]
                        current_stat   = stats[ii]

                        # Proposal is rejected
                        proposals[jump] = np.array(current_params)
                        stats[jump]  = current_stat


        finally:
            self._sampler.tear_down()
            #progress_bar(niter, niter, tstart, self._sampler.__class__.__name__)

        params = proposals.transpose()
        return (stats, acceptflag, params)



class PragBayes(MetropolisMH):

    def __init__(self, fcn, sigma, mu, dof, fit, *args):
        MetropolisMH.__init__(self, fcn, sigma, mu, dof, *args)

        self._fit = fit
        if hasattr(fit.model, 'teardown'):
            fit.model.teardown()
        self.srcarfs, self.bkgarfs = search_arfs(fit)

        # Save a copy of original ARF
        self.arfs = flatten_arfs(self.srcarfs, self.bkgarfs)
        self.backup_arfs = [arf.specresp.copy() for arf in self.arfs]

        self.simarf = None
        

    def init(self, log=False, inv=False, defaultprior=True, priorshape=False,
             priors=(), originalscale=True, scale=1, sigma_m=False, p_M=.5,
             simarf=None, nsubiters=10):

        # Note that nsubiters is used as a dummy parameter to indicate the default
        # value.  See the function WalkWithSubIters.__call__() 

        if isinstance(simarf, (PCA1DAdd, SIM1DAdd)):
            self.simarf = simarf
        else:
            self.simarf = ARFSIMFactory()(simarf)

        return MetropolisMH.init(self, log, inv, defaultprior, priorshape,
                                 priors, originalscale, scale, sigma_m, p_M)


    def fit(self, current):
        self._fit.model.thawedpars = current

        # nm = NelderMead()
        # nm.config['iquad'] = 0
        # nm.config['finalsimplex'] = 1

        #lm = LevMar()
        #lm.config['maxfev'] = 5

        cv = Covariance()

        # Use the fit method defined before called get_draws().  This way the user
        # does not have to pass in the fitting method and method options.
        fit = Fit(self._fit.data, self._fit.model, self._fit.stat, self._fit.method,
                  cv)

        fit_result = fit.fit()
        covar_result = fit.est_errors()
        sigma = np.array(covar_result.extra_output)

        if np.isnan(sigma).any():
            raise CovarError("NaNs found in covariance matrix")

        # cache the fitting scales
        self._sigma = sigma


    def mh(self, current):
        """ MH jumping rule """
        
        # The current proposal is ignored here.
        # MH jumps from the current best-fit parameter values using the 
        # covariance scale from the sub-iteration fit
        
        # If the ARF is updated then sigma will be None, then
        # refit and run covariance.  fit() will set the _sigma.
        if self._sigma is None:
            self.fit(current)

        # Use _sigma from fit in the first subiteration.
        # Reset for each iteration.
        proposal = rmvt(current, self._sigma, self._dof)
        return proposal


    def perturb_arf(self, current_params, current_stat):
        if self.simarf is not None:
            # add deviations starting with original ARF for each iter
            for specresp, arf in zip(self.backup_arfs, self.arfs):
                arf.specresp = self.simarf.add_deviations(specresp)

            # When ARF is updated, set scale to None
            self._sigma = None


    def tear_down(self):
        MetropolisMH.tear_down(self)
        fit = self._fit

        # Restore ARF to original state
        for specresp, arf in zip(self.backup_arfs, self.arfs):
            arf.specresp = specresp
