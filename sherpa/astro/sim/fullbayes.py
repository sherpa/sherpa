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


import numpy as np
from sherpa.sim.mh import MetropolisMH, dmvnorm
from sherpa.astro.sim.pragbayes import PragBayes, PCA1DAdd, ARFSIMFactory


__all__ = ('FullBayes',)


class FullBayes(PragBayes):

    def __init__(self, fcn, sigma, mu, dof, fit, *args):
        PragBayes.__init__(self, fcn, sigma, mu, dof, fit, *args)
        self.arf_dicts = [{'current' : arf.specresp.copy(), 'old_rr' : None }
                          for arf in self.arfs]


    def init(self, log=False, inv=False, defaultprior=True, priorshape=False,
             priors=(), originalscale=True, scale=1, sigma_m=False, p_M=.5,
             simarf=None, p_M_arf=.5, sigma_arf=0.1):
        # nsubiters is missing from init() to indicate that nubiters=1 for 
        # full bayes

        if isinstance(simarf, (PCA1DAdd,)):
            self.simarf = simarf
        else:
            self.simarf = ARFSIMFactory()(simarf)

        if not isinstance(self.simarf, (PCA1DAdd,)):
            raise TypeError("Simulation ARF must be PCA for FullBayes" + 
                            " not %s" % type(self.simarf).__name__)

        self.accept_arfs = [0]
        self.p_M_arf = p_M_arf
        self.rrsig = sigma_arf
        return MetropolisMH.init(self, log, inv, defaultprior, priorshape,
                                 priors, originalscale, scale, sigma_m, p_M)


    def _update_arf(self, arf, specresp, current_params, current_stat, arf_dict):

        u = np.random.uniform(0,1,1)
        if u > self.p_M_arf:

            ncomp = self.simarf.ncomp
            old_rr = arf_dict['old_rr']
            if old_rr is None:
                old_rr = np.random.standard_normal(ncomp)

            # Assume the ARF is accepted by default
            # Update the ARFs with new deviates

            arf.specresp = self.simarf.add_deviations(specresp, old_rr, self.rrsig)
            new_rr = self.simarf.rrout

            stat_temp = self.calc_fit_stat(self._mu)

            mu0  = np.repeat(0,ncomp)
            sig0 = np.diag(np.repeat(1,ncomp))
            accept_pr = dmvnorm(new_rr,mu0,sig0)-dmvnorm(old_rr,mu0,sig0)
            accept_pr += stat_temp - current_stat
            accept_pr = np.exp( accept_pr )

            uu = np.random.uniform(0,1,1)
            if accept_pr > uu:
                arf_dict['old_rr'] = new_rr
                arf_dict['current'] = arf.specresp.copy()
                self.accept_arf()

                # When ARF is updated, set scale to None to signal a fit
                self._sigma = None

            else:
                # Restore to previous ARF
                arf.specresp = arf_dict['current']
                self.reject_arf()

        else:
            # Assume the ARF is accepted by default
            # Update the ARFs with new deviates

            arf.specresp = self.simarf.add_deviations(specresp)
            
            stat_temp = self.calc_fit_stat(self._mu)
            accept_pr=0
            accept_pr += stat_temp - current_stat
            accept_pr = np.exp( accept_pr )

            uu = np.random.uniform(0,1,1)
            if accept_pr > uu:
                arf_dict['current'] = arf.specresp.copy()
                self.accept_arf()

                # When ARF is updated, set scale to None to signal a fit
                self._sigma = None

            else:
                # Restore to previous ARF
                arf.specresp = arf_dict['current']
                self.reject_arf()


    def accept_arf(self):
        # Accept the updated ARF deviates
        self.accept_arfs.append(self.accept_arfs[-1] + 1)


    def reject_arf(self):
        # Reject
        self.accept_arfs.append(self.accept_arfs[-1])


    def perturb_arf(self, current_params, current_stat):
        if self.simarf is not None:         

            # add deviations starting with original ARF for each iter
            for specresp, arf, arf_dict in zip(self.backup_arfs, self.arfs, self.arf_dicts):
                self._update_arf(arf, specresp, current_params, current_stat, arf_dict)
