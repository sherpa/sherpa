#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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
import numpy
import logging
logger = logging.getLogger("sherpa")

from sherpa.utils.testing import SherpaTestCase, requires_data, \
    requires_fits, requires_xspec
from sherpa.astro import ui

@requires_data
@requires_xspec
@requires_fits
class test_mhtest(SherpaTestCase):

    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.pha = 'p1_113+3380+3381+3382+3399+15293+14418_HRC-LEG_1_xspec.pha'
        self.rmf = 'hrcf00113_repro_leg_p1.rmf'
        self.bkg = 'p1_113+3380+3381+3382+3399+15293+14418_HRC-LEG_1_bkg.pha'
        self.PCAfname = 'aref_hrcsletg_8.fits'        
        self.fixarfname = 'p1_113+3380+3381+3382+3399+15293+14418_HRC-LEG_1.arf'
        ui.clean()
        
    def tearDown(self):
        if hasattr(self, '_old_logger_level'):
            logger.setLevel(self._old_logger_level)
        ui.clean()
        
    def test_me(self, tol=1.0e-3):
        ui.load_data(self.make_path(self.pha))        
        ui.load_rmf(self.make_path(self.rmf))
        ui.load_bkg(self.make_path(self.bkg))
        ui.set_analysis("wave")
        ui.ignore()
        ui.notice(25., 59.5)
        ui.notice(68., 80.)
        bkg_arf = ui.get_bkg_arf()
        bkg_arf.specresp = bkg_arf.specresp * 0 + 1.0
        bkg_scale = ui.get_bkg_scale()
        arf = ui.get_arf()
        ui.set_full_model(arf("xsphabs.abs1*xsbbody.bb1") + bkg_scale*bkg_arf("const1d.bkg_c1 * polynom1d.pp"))
        bkg_c1.c0 = 78.195
        abs1.nH = 0.00791702
        ui.set_par(bb1.kt, frozen=False, val=0.0620321, min=0.001, max=0.2)
        ui.set_par(bb1.norm, frozen=False, val=0.000285548)

        # define background model parameters
        ui.set_bkg_full_model(bkg_arf("polynom1d.pp"))
        ui.freeze(pp)
        pp.c0 = -2.08129
        pp.c1 = 0.349143
        pp.c2 = -0.0246803
        pp.c3 = 0.000968759
        pp.c4 = -2.312e-5
        pp.c5 = 3.4386e-7
        pp.c6 = -3.11659e-9
        pp.c7 = 1.57691e-11
        pp.c8 = -3.41824e-14
        ui.set_stat("cash")
        ui.fit()
        ui.covariance()
        fit_result = ui.get_fit_results()
        parnames = fit_result.parnames
        parvals = numpy.asarray(fit_result.parvals)
        covar_result = ui.get_covar_results()
        covar = covar_result.extra_output
        PCAfname = self.make_path(self.PCAfname)
        fixARFname = self.make_path(self.fixarfname)
        result = \
            ui.mh_sampling_newdata(parnames, parvals, 4, PCAfname, \
                                   fixARFname, niter=100, covar_matrix=covar, \
                                   improvedBayes=True, p_M_arf=0.5, comp=8, \
                                   sd_arf=.1, scale=1)
        assert numpy.allclose(187622.3458427461, numpy.average(result[0:10]),
                              tol, tol)
        assert numpy.allclose(573174.5150765756, numpy.std(result[0:10]),
                              tol, tol)
