#
#  Copyright (C) 2007, 2015, 2016, 2017, 2018, 2019
#     Smithsonian Astrophysical Observatory
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
import warnings
from numpy import sqrt
from pytest import approx
import pytest

from sherpa.utils.testing import SherpaTestCase, requires_data, \
    requires_fits, requires_xspec, requires_group

import sherpa.astro.ui as ui
from sherpa.astro.data import DataPHA

logger = logging.getLogger('sherpa')

try:
    from sherpa.astro import xspec
    has_xspec = True
except ImportError:
    has_xspec = False


# This value should match the [ogip] minimum_energy setting in the
# Sherpa configuration file.
#
EMIN = 1.0e-10

# has_xspec = has_package_from_list("sherpa.astro.xspec")


@requires_data
class test_threads(SherpaTestCase):

    def setUp(self):
        self.is_crates_io = False
        try:
            import sherpa.astro.io
            if "sherpa.astro.io.crates_backend" == sherpa.astro.io.backend.__name__:
                self.is_crates_io = True
        except ImportError:
            self.is_crates_io = False

        self.old_state = ui._session.__dict__.copy()
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.CRITICAL)

        # Store XSPEC settings, if applicable
        if has_xspec:
            self.old_xspec = xspec.get_xsstate()
            # As of XSPEC 12.10.1, it is safest to explicitly set
            # the state to a known value. For now just pick the
            # "easy" settings.
            #
            xspec.set_xsabund('angr')
            xspec.set_xsxsect('bcmc')

    def tearDown(self):
        ui._session.__dict__.update(self.old_state)
        logger.setLevel(self.old_level)

        if has_xspec:
            xspec.set_xsstate(self.old_xspec)

    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        super(test_threads, self).run_thread(name, scriptname=scriptname)

    @requires_fits
    @requires_xspec
    def test_pha_intro(self):
        self.run_thread('pha_intro')
        # astro.ui imported as ui, instead of
        # being in global namespace
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.0790393, rel=1e-4)
        assert covarerr[1] == approx(1.4564e-05, rel=1e-4)
        assert fit_results.statval == approx(37.9079, rel=1e-4)
        assert fit_results.rstat == approx(0.902569, rel=1e-4)
        assert fit_results.qval == approx(0.651155, rel=1e-4)
        assert self.locals['p1'].gamma.val == approx(2.15852, rel=1e-4)
        assert self.locals['p1'].ampl.val == approx(0.00022484, rel=1e-4)

        assert ui.calc_photon_flux() == approx(0.000469964, rel=1e-4)
        assert ui.calc_energy_flux() == approx(9.614847e-13, rel=1e-4)
        assert ui.calc_data_sum() == approx(706.85714092, rel=1e-4)
        assert ui.calc_model_sum() == approx(638.45693377, rel=1e-4)
        assert ui.calc_source_sum() == approx(0.046996409, rel=1e-4)

        calc = ui.eqwidth(self.locals['p1'], ui.get_source())
        assert calc == approx(-0.57731725, rel=1e-4)

        calc = ui.calc_kcorr([1, 1.2, 1.4, 1.6, 1.8, 2], 0.5, 2)
        expected = [0.93341286, 0.93752836, 0.94325233,
                    0.94990140, 0.95678054, 0.96393515]
        assert calc == approx(expected, rel=1e-4)

        self.assertEqual(ui.get_fit_results().nfev, 22)
        self.assertEqual(ui.get_fit_results().numpoints, 44)
        self.assertEqual(ui.get_fit_results().dof, 42)

    @requires_fits
    def test_pha_read(self):
        self.run_thread('pha_read')
        self.assertEqual(type(ui.get_data()), DataPHA)

    @requires_fits
    def test_basic(self):
        # In data1.dat for this test, there is a comment with one
        # word at the beginning -- deliberately would break when reading
        # with DM ASCII kernel, but passes because we have Sherpa code
        # to bypass that.
        self.run_thread('basic')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.0192539, rel=1e-4)
        assert covarerr[1] == approx(0.00392255, rel=1e-4)
        assert fit_results.statval == approx(151.827, rel=1e-4)
        assert fit_results.rstat == approx(16.8697, rel=1e-4)
        assert fit_results.qval == approx(3.68798e-28, rel=1e-4)
        assert self.locals['m1'].c0.val == approx(1.49843, rel=1e-4)
        assert self.locals['m1'].c1.val == approx(0.1447, rel=1e-4)
        assert self.locals['m1'].c2.val == approx(0.0322936, rel=1e-4)
        assert self.locals['m1'].c3.val == approx(-0.00277729, rel=1e-4)
        assert self.locals['m2'].c0.val == approx(1.75548, rel=1e-4)
        assert self.locals['m2'].c1.val == approx(0.198455, rel=1e-4)
        self.assertEqual(fit_results.nfev, 9)
        self.assertEqual(fit_results.numpoints, 11)
        self.assertEqual(fit_results.dof, 9)

    @requires_fits
    @requires_xspec
    def test_simultaneous(self):
        self.run_thread('simultaneous')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.397769, rel=1e-3)
        assert covarerr[1] == approx(0.486058, rel=1e-3)
        assert covarerr[2] == approx(1.48213e-05, rel=1e-3)
        assert covarerr[3] == approx(1.54245e-05, rel=1e-3)
        assert fit_results.statval == approx(7.4429, rel=1e-4)
        assert fit_results.rstat == approx(0.531636, rel=1e-4)
        assert fit_results.qval == approx(0.916288, rel=1e-4)
        assert self.locals['abs1'].nh.val == approx(0.898162, rel=1e-2)
        assert self.locals['pl1'].gamma.val == approx(1.645, rel=1e-4)
        self.assertEqualWithinTol(self.locals['pl1'].ampl.val,
                                  2.28323e-05, 1e-3)
        self.assertEqualWithinTol(self.locals['pl2'].ampl.val,
                                  2.44585e-05, 1e-3)
        self.assertEqual(fit_results.numpoints, 18)
        self.assertEqual(fit_results.dof, 14)

    @requires_fits
    @requires_xspec
    def test_sourceandbg(self):
        self.run_thread('sourceandbg')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.012097, rel=1e-3)
        assert covarerr[1] == approx(0, rel=1e-3)
        assert covarerr[2] == approx(0.000280678, rel=1e-3)
        assert covarerr[3] == approx(0.00990783, rel=1e-3)
        assert covarerr[4] == approx(2.25746e-07, rel=1e-3)
        assert fit_results.statval == approx(947.5, rel=1e-4)
        assert fit_results.rstat == approx(0.715094, rel=1e-4)
        assert fit_results.qval == approx(1, rel=1e-4)
        assert self.locals['a1'].nh.val == approx(0.0342266, rel=1e-2)
        assert self.locals['b1'].kt.val == approx(20, rel=1e-2)
        self.assertEqualWithinTol(self.locals['b1'].norm.val,
                                  0.00953809, 1e-2)
        self.assertEqualWithinTol(self.locals['b2'].kt.val,
                                  0.563109, 1e-2)
        self.assertEqualWithinTol(self.locals['b2'].norm.val,
                                  1.16118e-05, 1e-2)
        self.assertEqual(fit_results.numpoints, 1330)
        self.assertEqual(fit_results.dof, 1325)

    @requires_fits
    def test_spatial(self):
        self.run_thread('spatial')
        self.assertEqualWithinTol(ui.get_fit_results().statval,
                                  -59229.749441, 1e-4)
        assert self.locals['g1'].fwhm.val == approx(61.5615, rel=1e-4)
        assert self.locals['g1'].xpos.val == approx(4070.45, rel=1e-4)
        assert self.locals['g1'].ypos.val == approx(4251.35, rel=1e-4)
        assert self.locals['g1'].ampl.val == approx(22.1269, rel=1e-4)
        assert self.locals['g2'].fwhm.val == approx(6.20409, rel=1e-4)
        assert self.locals['g2'].xpos.val == approx(4070.78, rel=1e-4)
        assert self.locals['g2'].ypos.val == approx(4249.33, rel=1e-4)
        assert self.locals['g2'].ampl.val == approx(226.563, rel=1e-4)
        # self.assertEqual(ui.get_fit_results().nfev, 371)
        self.assertEqual(ui.get_fit_results().numpoints, 4881)
        self.assertEqual(ui.get_fit_results().dof, 4877)

    @requires_fits
    def test_radpro(self):
        self.run_thread('radpro')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(9.37345, rel=1e-4)
        assert covarerr[1] == approx(0.512596, rel=1e-4)
        assert covarerr[2] == approx(0.0691102, rel=1e-4)
        assert fit_results.statval == approx(217.450, rel=1e-4)
        assert fit_results.rstat == approx(6.21287, rel=1e-4)
        assert fit_results.qval == approx(0.0, rel=1e-4)
        assert self.locals['src'].r0.val == approx(125.829, rel=1e-4)
        assert self.locals['src'].beta.val == approx(4.1633, rel=1e-4)
        assert self.locals['src'].xpos.val == approx(0.0, rel=1e-4)
        assert self.locals['src'].ampl.val == approx(4.42821, rel=1e-4)
        self.assertEqual(fit_results.nfev, 92)
        self.assertEqual(fit_results.numpoints, 38)
        self.assertEqual(fit_results.dof, 35)

    def test_radpro_dm(self):
        # This test is completely redundant to test_radpro above.
        # The only difference is that here I test if using DM syntax
        # in the file name returns the correct arrays from CRATES.
        # Since any other I/O backend would not understand DM syntax,
        # I make this test a no-op unless CRATES is used as the I/O
        # backend.  SMD 05/23/13
        if self.is_crates_io:
            self.run_thread('radpro_dm')

            fres = ui.get_fit_results()
            covarerr = sqrt(fres.extra_output['covar'].diagonal())
            assert covarerr[0] == approx(9.37344776, rel=1e-4)
            assert covarerr[1] == approx(0.51259645, rel=1e-4)
            assert covarerr[2] == approx(0.06911017, rel=1e-4)
            assert fres.statval == approx(217.450, rel=1e-4)
            assert fres.rstat == approx(6.21287, rel=1e-4)
            assert fres.qval == approx(0.0, rel=1e-4)

            srcmdl = self.locals['src']
            assert srcmdl.r0.val == approx(125.829, rel=1e-4)
            assert srcmdl.beta.val == approx(4.1633, rel=1e-4)
            assert srcmdl.xpos.val == approx(0.0, rel=1e-4)
            assert srcmdl.ampl.val == approx(4.42821, rel=1e-4)

            self.assertEqual(fres.nfev, 92)
            self.assertEqual(fres.numpoints, 38)
            self.assertEqual(fres.dof, 35)

    @requires_fits
    def test_psf2d(self):
        self.run_thread('psf')
        assert ui.get_fit_results().statval == approx(4066.78, rel=1e-4)
        assert self.locals['g1'].fwhm.val == approx(2.80117, rel=1e-2)
        assert self.locals['g1'].ypos.val == approx(77.2271, rel=1e-2)
        assert self.locals['g1'].xpos.val == approx(88.661, rel=1e-2)
        assert self.locals['g1'].ampl.val == approx(166.649, rel=1e-2)
        # self.assertEqual(ui.get_fit_results().nfev, 342)
        self.assertEqual(ui.get_fit_results().numpoints, 4899)
        self.assertEqual(ui.get_fit_results().dof, 4895)

    @requires_fits
    def test_fpsf2d(self):
        self.run_thread('fpsf')

        fres = ui.get_fit_results()
        assert fres.statval == approx(-4053.6635, rel=1e-4)

        # assert self.locals['b1'].xlow.val == approx(-4.70832, rel=1e-4)
        # assert self.locals['b1'].xhi.val == approx(164.687, rel=1e-4)
        # assert self.locals['b1'].ylow.val == approx(0.83626, rel=1e-4)
        # assert self.locals['b1'].yhi.val == approx(142.603, rel=1e-4)
        # assert self.locals['b1'].ampl.val == approx(0.956766, rel=1e-4)

        g1mdl = self.locals['g1']
        assert g1mdl.fwhm.val == approx(6.420237, rel=1e-4)
        assert g1mdl.xpos.val == approx(88.940712, rel=1e-4)
        assert g1mdl.ypos.val == approx(76.577265, rel=1e-4)
        assert g1mdl.ampl.val == approx(36344.48324, rel=1e-4)

        # self.assertEqual(fres.nfev, 978)
        self.assertEqual(fres.numpoints, 4899)
        self.assertEqual(fres.dof, 4895)

    @requires_fits
    def test_radpro_psf(self):
        self.run_thread('radpro_psf')
        assert ui.get_fit_results().statval == approx(200.949, rel=1e-4)
        assert self.locals['src'].r0.val == approx(83.0997, rel=1e-4)
        assert self.locals['src'].beta.val == approx(2.97737, rel=1e-4)
        assert self.locals['src'].ampl.val == approx(5.27604, rel=1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 48)
        self.assertEqual(ui.get_fit_results().numpoints, 38)
        self.assertEqual(ui.get_fit_results().dof, 35)

    @requires_fits
    def test_linepro(self):
        self.run_thread('linepro')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.176282, rel=1e-4)
        assert covarerr[1] == approx(0.0019578, rel=1e-4)
        assert covarerr[2] == approx(0.495889, rel=1e-4)
        assert fit_results.statval == approx(203.34, rel=1e-4)
        assert self.locals['b1'].r0.val == approx(4.25557, rel=1e-4)
        assert self.locals['b1'].beta.val == approx(0.492232, rel=1e-4)
        assert self.locals['b1'].ampl.val == approx(11.8129, rel=1e-4)
        self.assertEqual(fit_results.nfev, 17)
        self.assertEqual(fit_results.numpoints, 75)
        self.assertEqual(fit_results.dof, 72)

    @requires_fits
    def test_kernel(self):
        self.run_thread('kernel')
        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.210895, rel=1e-4)
        assert covarerr[1] == approx(0.00154839, rel=1e-4)
        assert covarerr[2] == approx(0.0223859, rel=1e-4)
        assert fit_results.statval == approx(98.5793, rel=1e-4)
        assert self.locals['b1'].r0.val == approx(19.2278, rel=1e-4)
        assert self.locals['b1'].beta.val == approx(0.555464, rel=1e-4)
        assert self.locals['b1'].ampl.val == approx(1.93706, rel=1e-4)
        self.assertEqual(fit_results.nfev, 21)
        self.assertEqual(fit_results.numpoints, 75)
        self.assertEqual(fit_results.dof, 72)

    @requires_fits
    @requires_xspec
    def test_spectrum(self):
        self.run_thread('spectrum')

        fres = ui.get_fit_results()
        covarerr = sqrt(fres.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.00148391, rel=1e-3)
        assert covarerr[1] == approx(0.0011518, rel=1e-3)
        assert covarerr[2] == approx(0.00377755, rel=1e-3)
        assert covarerr[3] == approx(0.00370543, rel=1e-3)
        assert covarerr[4] == approx(0.0016608, rel=1e-3)
        assert fres.statval == approx(0.0496819, rel=1e-4)
        assert self.locals['abs2'].nh.val == approx(1.1015, rel=1e-4)
        assert self.locals['mek1'].kt.val == approx(0.841025, rel=1e-4)
        assert self.locals['mek1'].norm.val == approx(0.699761, rel=1e-4)
        assert self.locals['mek2'].kt.val == approx(2.35845, rel=1e-4)
        assert self.locals['mek2'].norm.val == approx(1.03724, rel=1e-4)
        self.assertEqual(fres.numpoints, 446)
        self.assertEqual(fres.dof, 441)

    @requires_fits
    def test_histo(self):
        self.run_thread('histo')
        assert ui.get_fit_results().statval == approx(14.7264, rel=1e-4)
        assert self.locals['g1'].fwhm.val == approx(0.0232473, rel=1e-4)
        assert self.locals['g1'].pos.val == approx(1.26713, rel=1e-4)
        assert self.locals['g1'].ampl.val == approx(40.4503, rel=1e-4)
        # self.assertEqual(ui.get_fit_results().nfev, 19)
        self.assertEqual(ui.get_fit_results().numpoints, 50)
        self.assertEqual(ui.get_fit_results().dof, 47)

    @requires_fits
    @requires_xspec
    def test_xmm(self):
        self.run_thread('xmm')

        fres = ui.get_fit_results()
        covarerr = sqrt(fres.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.954993, rel=1e-3)
        assert covarerr[1] == approx(0.142357, rel=1e-3)
        assert covarerr[2] == approx(0.00038775, rel=1e-3)
        assert fres.statval == approx(118.085, rel=1e-4)
        self.assertEqualWithinTol(self.locals['intrin'].nh.val,
                                  11.0769, 1e-2)
        self.assertEqualWithinTol(self.locals['phard'].phoindex.val,
                                  1.49055, 1e-2)
        self.assertEqualWithinTol(self.locals['phard'].norm.val,
                                  0.00140301, 1e-2)

        self.assertEqual(fres.nfev, 95)
        self.assertEqual(fres.numpoints, 162)
        self.assertEqual(fres.dof, 159)

    @requires_fits
    # As of CIAO 4.5, can filter on channel number, even when
    # data are grouped! Test results should exactly match CIAO 4.4
    # fit results in grouped/fit.py
    def test_grouped_ciao4_5(self):
        self.run_thread('grouped_ciao4.5')

        fres = ui.get_fit_results()
        covarerr = sqrt(fres.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0.104838, rel=1e-4)
        assert covarerr[1] == approx(2.43937e-05, rel=1e-4)
        assert fres.statval == approx(18.8316, rel=1e-4)
        self.assertEqual(fres.numpoints, 46)

        aamdl = self.locals['aa']
        assert aamdl.gamma.val == approx(1.83906, rel=1e-4)
        assert aamdl.ampl.val == approx(0.000301258, rel=1e-4)

    @requires_fits
    @requires_xspec
    def test_proj(self):
        self.run_thread('proj')
        # Fit 1
        self.assertEqualWithinTol(self.locals['fit_res1'].statval,
                                  64.3803, 1e-2)
        assert self.locals['g'].fwhm.val == approx(32.3536, rel=1e-2)
        assert self.locals['g'].pos.val == approx(807.863, rel=1e-2)
        assert self.locals['g'].ampl.val == approx(117.826, rel=1e-2)

        # Covar 1
        self.assertEqualWithinTol(self.locals['covar_res1'].parmins[0],
                                  -0.49357, 1e-2)
        self.assertEqualWithinTol(self.locals['covar_res1'].parmins[1],
                                  -0.264056, 1e-2)
        self.assertEqualWithinTol(self.locals['covar_res1'].parmins[2],
                                  -2.58857, 1e-2)
        self.assertEqualWithinTol(self.locals['covar_res1'].parmaxes[0],
                                  0.49357, 1e-2)
        self.assertEqualWithinTol(self.locals['covar_res1'].parmaxes[1],
                                  0.264056, 1e-2)
        self.assertEqualWithinTol(self.locals['covar_res1'].parmaxes[2],
                                  2.58857, 1e-2)

        # Projection 1
        self.assertEqualWithinTol(self.locals['proj_res1'].parmins[0],
                                  -0.492815, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res1'].parmins[1],
                                  -0.263945, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res1'].parmins[2],
                                  -2.57957, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res1'].parmaxes[0],
                                  0.49421, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res1'].parmaxes[1],
                                  0.26389, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res1'].parmaxes[2],
                                  2.59695, 1e-2)

        # Fit 2
        self.assertEqualWithinTol(self.locals['fit_res2'].statval,
                                  108.048, 1e-2)
        assert self.locals['xs1'].nh.val == approx(0.0532586, rel=1e-2)
        assert self.locals['p'].gamma.val == approx(1.4816, rel=1e-2)
        assert self.locals['p'].ampl.val == approx(0.302343, rel=1e-2)

        # Covar 2
        self.assertEqualWithinTol(self.locals['covar_res2'].parmins[0],
                                  -0.0379975, 5e-1)
        self.assertEqualWithinTol(self.locals['covar_res2'].parmins[1],
                                  -0.0877478, 5e-1)
        self.assertEqualWithinTol(self.locals['covar_res2'].parmins[2],
                                  -0.0831945, 5e-1)
        self.assertEqualWithinTol(self.locals['covar_res2'].parmaxes[0],
                                  0.0379975, 5e-1)
        self.assertEqualWithinTol(self.locals['covar_res2'].parmaxes[1],
                                  0.0877478, 5e-1)
        self.assertEqualWithinTol(self.locals['covar_res2'].parmaxes[2],
                                  0.0831945, 5e-1)

        # Projection 2
        self.assertEqualWithinTol(self.locals['proj_res2'].parmins[0],
                                  -0.0385636, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res2'].parmins[1],
                                  -0.0891261, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res2'].parmins[2],
                                  -0.0737413, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res2'].parmaxes[0],
                                  0.0388651, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res2'].parmaxes[1],
                                  0.0896566, 1e-2)
        self.assertEqualWithinTol(self.locals['proj_res2'].parmaxes[2],
                                  0.0981627, 1e-2)

    @requires_fits
    @requires_xspec
    def test_proj_bubble(self):
        xspec.set_xsxsect('bcmc')

        self.run_thread('proj_bubble')

        fit_results = ui.get_fit_results()
        covarerr = sqrt(fit_results.extra_output['covar'].diagonal())
        assert covarerr[0] == approx(0, rel=1e-4)
        assert covarerr[1] == approx(8.74608e-07, rel=1e-2)

        # Fit -- Results from reminimize
        assert self.locals['mek1'].kt.val == approx(17.8849, rel=1e-2)
        assert self.locals['mek1'].norm.val == approx(4.15418e-06, rel=1e-2)

        # Fit -- Results from reminimize

        # The fit results change in XSPEC 12.10.0 since the mekal model
        # was changed (FORTRAN to C++). A 1% difference is used for the
        # parameter ranges from covar and proj (matches the tolerance for
        # the fit results).

        # Covar
        #
        # TODO: should this check that parmaxes is -1 * parmins instead?
        covar = ui.get_covar_results()
        assert covar.parmins[0] == approx(-0.328832, rel=0.1)
        assert covar.parmins[1] == approx(-8.847916e-7, rel=0.01)
        assert covar.parmaxes[0] == approx(0.328832, rel=0.1)
        assert covar.parmaxes[1] == approx(8.847916e-7, rel=0.01)

        # Proj -- Upper bound of kT can't be found
        #
        proj = ui.get_proj_results()
        assert proj.parmins[0] == approx(-12.048069, rel=0.01)
        assert proj.parmins[1] == approx(-9.510913e-07, rel=0.01)
        assert proj.parmaxes[1] == approx(2.403640e-06, rel=0.01)

        assert proj.parmaxes[0] is None

    # New tests based on SDS threads -- we should catch these errors
    # (if any occur) so SDS doesn't waste time tripping over them.
    @requires_fits
    @requires_xspec
    def test_counts(self):
        self.run_thread('counts')
        assert self.locals['counts_data1'] == approx(52701.0, rel=1e-4)
        assert self.locals['counts_data2'] == approx(25032.0, rel=1e-4)
        self.assertEqualWithinTol(self.locals['counts_model1'],
                                  73226263.55355, 1e-4)
        self.assertEqualWithinTol(self.locals['eflux1'],
                                  1.495482e-08, 1e-4)
        self.assertEqualWithinTol(self.locals['counts_model2'],
                                  46082543.21529, 1e-4)
        self.assertEqualWithinTol(self.locals['eflux2'],
                                  1.39662954483e-08, 1e-3)
        self.assertEqualWithinTol(self.locals['pflux1'],
                                  1.6178938637, 1e-2)

    @requires_fits
    @requires_xspec
    def test_stats_all(self):
        self.run_thread('stats_all')
        self.assertEqualWithinTol(self.locals['stat_lsqr'],
                                  213746.236464, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_chi2m'],
                                  1232.0330242, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_cash'],
                                  -1163566.90688, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_chi2c'],
                                  1411.60090961, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_chi2g'],
                                  972.388468358, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_chi2d'],
                                  1204.69363458, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_chi2x'],
                                  1204.69363458, 1e-4)
        self.assertEqualWithinTol(self.locals['stat_cstat'],
                                  1210.56896183, 1e-4)

    @requires_fits
    def test_lev3fft(self):
        self.run_thread('lev3fft', scriptname='bar.py')
        self.assertEqualWithinTol(self.locals['src'].fwhm.val,
                                  0.044178, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].xpos.val,
                                  150.016, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ypos.val,
                                  2.66493839, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ampl.val,
                                  1.56090546, 1e-4)
        self.assertEqualWithinTol(self.locals['bkg'].c0.val,
                                  -1.513700715, 1e-4)

        fres = ui.get_fit_results()
        assert fres.istatval == approx(19496.3, rel=1e-4)
        assert fres.statval == approx(592.32647, rel=1e-4)
        self.assertEqual(fres.numpoints, 3307)
        self.assertEqual(fres.dof, 3302)

    @requires_fits
    @requires_xspec
    @requires_group
    def test_setfullmodel(self):
        self.run_thread('setfullmodel')

    @requires_fits
    def test_bug13537(self):
        self.run_thread('bug13537')

    @requires_fits
    @requires_xspec
    def test_xmm2(self):

        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            self.run_thread('xmm2')

        assert len(ws) == 2
        cats = set([w.category for w in ws])
        assert cats == set([UserWarning])

        # The order of reading the ARF and RMF is not guaranteed,
        # so do not force it here when testing the two warning
        # messages.
        #
        arffile = 'MNLup_2138_0670580101_EMOS1_S001_spec.arf'
        rmffile = 'MNLup_2138_0670580101_EMOS1_S001_spec.rmf'
        emsg_arf = "The minimum ENERG_LO in the ARF " + \
                   "'{}' ".format(arffile) + \
                   "was 0 and has been replaced by {}".format(EMIN)
        emsg_rmf = "The minimum ENERG_LO in the RMF " + \
                   "'{}' ".format(rmffile) + \
                   "was 0 and has been replaced by {}".format(EMIN)

        emsgs = set([emsg_arf, emsg_rmf])
        wmsgs = set([str(w.message) for w in ws])
        assert wmsgs == emsgs

        assert ui.get_data().channel[0] == approx(1.0, rel=1e-4)
        rmf = ui.get_rmf()
        arf = ui.get_arf()
        self.assertEqual(rmf.detchans, 800)
        self.assertEqual(len(rmf.energ_lo), 2400)
        self.assertEqual(len(rmf.energ_hi), 2400)
        self.assertEqual(len(rmf.n_grp), 2400)
        self.assertEqual(len(rmf.f_chan), 2394)
        self.assertEqual(len(rmf.n_chan), 2394)
        self.assertEqual(len(rmf.matrix), 1281216)
        self.assertEqual(rmf.offset, 0)
        self.assertEqual(len(rmf.e_min), 800)
        self.assertEqual(len(rmf.e_max), 800)
        self.assertEqual(len(arf.energ_lo), 2400)
        self.assertEqual(len(arf.energ_hi), 2400)
        self.assertEqual(len(arf.specresp), 2400)

        etol = EMIN / 100.0
        assert rmf.energ_lo[0] == approx(EMIN, rel=etol)
        assert arf.energ_lo[0] == approx(EMIN, rel=etol)


@requires_data
@requires_fits
@requires_xspec
def test_missmatch_arf(make_data_path):
    ui.load_pha(1, make_data_path("source1.pi"))
    ui.load_bkg(1, make_data_path("back1.pi"))
    ui.load_arf(1, make_data_path("arf_1024.fits"))
    ui.load_rmf(1, make_data_path("rmf_1024.fits"))
    ui.set_method('levmar')
    ui.set_model(ui.powlaw1d.p1 * ui.xswabs.abs1)
    ui.set_par('p1.ampl', 0.0001)
    ui.set_stat('cash')
    ui.fit()
    parvals = ui.get_fit_results().parvals
    assert parvals[0] == approx(1.47969, rel=1.0e-3)
    assert parvals[1] == approx(0.0019491, rel=1.0e-3)
    assert parvals[2] == approx(2.35452, rel=1.0e-3)


# This was test_threads.test_pileup
#
@requires_data
@requires_fits
@requires_xspec
@pytest.mark.usefixtures("clean_astro_ui")
def test_thread_pileup(run_thread):

    models = run_thread('pileup')

    fr = ui.get_fit_results()
    covarerr = sqrt(fr.extra_output['covar'].diagonal())

    assert fr.statval == approx(53.6112, rel=1e-4)
    assert fr.rstat == approx(1.44895, rel=1e-4)
    assert fr.qval == approx(0.0379417, rel=1e-4)
    assert fr.numpoints == 42
    assert fr.dof == 37

    jdp = models['jdp']
    assert jdp.alpha.val == approx(0.522593, rel=1e-1)
    assert jdp.f.val == approx(0.913458, rel=1e-2)

    abs1 = models['abs1']
    assert abs1.nh.val == approx(6.12101, rel=1e-2)

    power = models['power']
    assert power.gamma.val == approx(1.41887, rel=1e-2)
    assert power.ampl.val == approx(0.00199457, rel=1e-2)

    # Move covariance checks to the end of the parameters
    # to check whether the fit is the same, even if the
    # covariance is different.
    #
    assert covarerr[0] == approx(684.056, rel=1e-4)
    assert covarerr[1] == approx(191.055, rel=1e-3)
    assert covarerr[2] == approx(0.632061, rel=1e-3)
    assert covarerr[3] == approx(0.290159, rel=1e-3)
    assert covarerr[4] == approx(1.62529, rel=1e-3)

    # Issue #294 was a problem with serializing the pileup model
    # after a fit in Python 3 (but not Python 2). Add some basic
    # validation that the conversion to a string works. For the
    # pileup model we expect the standard model layout - e.g.
    #
    #   jdp
    #   paramter headers
    #   ---- ----- ...
    #   jdp.alpha ...
    #   ...
    #   jdp.nterms ...
    #   <blank line>
    #   1: ...
    #   ...
    #   7: ...
    #   *** pileup fraction: value
    #
    lines = str(jdp).split('\n')
    assert len(lines) == 19
    assert lines[10].strip() == ''
    assert lines[11].startswith('   1: ')
    assert lines[17].startswith('   7: ')
    assert lines[18].startswith('   *** pileup fraction: ')
