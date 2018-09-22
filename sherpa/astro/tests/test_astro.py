#
#  Copyright (C) 2007, 2015, 2016, 2017, 2018  Smithsonian Astrophysical Observatory
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

from numpy.testing import assert_allclose

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
            if ("sherpa.astro.io.crates_backend" ==
                sherpa.astro.io.backend.__name__):
                self.is_crates_io = True
        except:
            self.is_crates_io = False
        self.old_state = ui._session.__dict__.copy()
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.CRITICAL)

        # Store XSPEC settings, if applicable
        if has_xspec:
            self.old_xspec = xspec.get_xsstate()

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
        self.assertEqualWithinTol(ui.get_fit_results().statval, 37.9079, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().rstat, 0.902569, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().qval, 0.651155, 1e-4)
        self.assertEqualWithinTol(self.locals['p1'].gamma.val, 2.15852, 1e-4)
        self.assertEqualWithinTol(self.locals['p1'].ampl.val, 0.00022484, 1e-4)

        self.assertEqualWithinTol(ui.calc_photon_flux(), 0.000469964, 1e-4)
        self.assertEqualWithinTol(ui.calc_energy_flux(), 9.614847e-13, 1e-4)
        self.assertEqualWithinTol(ui.calc_data_sum(), 706.85714092, 1e-4)
        self.assertEqualWithinTol(ui.calc_model_sum(), 638.45693377, 1e-4)
        self.assertEqualWithinTol(ui.calc_source_sum(), 0.046996409, 1e-4)

        calc = ui.eqwidth(self.locals['p1'], ui.get_source())
        self.assertEqualWithinTol(calc, -0.57731725, 1e-4)

        calc = ui.calc_kcorr([1, 1.2, 1.4, 1.6, 1.8, 2], 0.5, 2)
        expected = [0.93341286, 0.93752836, 0.94325233,
                    0.94990140, 0.95678054, 0.96393515]
        self.assertEqualWithinTol(calc, expected, 1e-4)

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
        self.assertEqualWithinTol(ui.get_fit_results().statval, 151.827, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().rstat, 16.8697, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().qval, 3.68798e-28, 1e-4)
        self.assertEqualWithinTol(self.locals['m1'].c0.val, 1.49843, 1e-4)
        self.assertEqualWithinTol(self.locals['m1'].c1.val, 0.1447, 1e-4)
        self.assertEqualWithinTol(self.locals['m1'].c2.val, 0.0322936, 1e-4)
        self.assertEqualWithinTol(self.locals['m1'].c3.val, -0.00277729, 1e-4)
        self.assertEqualWithinTol(self.locals['m2'].c0.val, 1.75548, 1e-4)
        self.assertEqualWithinTol(self.locals['m2'].c1.val, 0.198455, 1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 9)
        self.assertEqual(ui.get_fit_results().numpoints, 11)
        self.assertEqual(ui.get_fit_results().dof, 9)

    @requires_fits
    @requires_xspec
    def test_simultaneous(self):
        self.run_thread('simultaneous')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 7.4429, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().rstat, 0.531636, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().qval, 0.916288, 1e-4)
        self.assertEqualWithinTol(self.locals['abs1'].nh.val, 0.898162, 1e-2)
        self.assertEqualWithinTol(self.locals['pl1'].gamma.val, 1.645, 1e-4)
        self.assertEqualWithinTol(self.locals['pl1'].ampl.val,
                                  2.28323e-05, 1e-3)
        self.assertEqualWithinTol(self.locals['pl2'].ampl.val,
                                  2.44585e-05, 1e-3)
        self.assertEqual(ui.get_fit_results().numpoints, 18)
        self.assertEqual(ui.get_fit_results().dof, 14)

    @requires_fits
    @requires_xspec
    def test_sourceandbg(self):
        self.run_thread('sourceandbg')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 947.5, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().rstat, 0.715094, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().qval, 1, 1e-4)
        self.assertEqualWithinTol(self.locals['a1'].nh.val, 0.0342266, 1e-2)
        self.assertEqualWithinTol(self.locals['b1'].kt.val, 20, 1e-2)
        self.assertEqualWithinTol(self.locals['b1'].norm.val,
                                  0.00953809, 1e-2)
        self.assertEqualWithinTol(self.locals['b2'].kt.val,
                                  0.563109, 1e-2)
        self.assertEqualWithinTol(self.locals['b2'].norm.val,
                                  1.16118e-05, 1e-2)
        self.assertEqual(ui.get_fit_results().numpoints, 1330)
        self.assertEqual(ui.get_fit_results().dof, 1325)

    @requires_fits
    def test_spatial(self):
        self.run_thread('spatial')
        self.assertEqualWithinTol(ui.get_fit_results().statval,
                                  -59229.749441, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].fwhm.val, 61.5615, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].xpos.val, 4070.45, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].ypos.val, 4251.35, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].ampl.val, 22.1269, 1e-4)
        self.assertEqualWithinTol(self.locals['g2'].fwhm.val, 6.20409, 1e-4)
        self.assertEqualWithinTol(self.locals['g2'].xpos.val, 4070.78, 1e-4)
        self.assertEqualWithinTol(self.locals['g2'].ypos.val, 4249.33, 1e-4)
        self.assertEqualWithinTol(self.locals['g2'].ampl.val, 226.563, 1e-4)
        # self.assertEqual(ui.get_fit_results().nfev, 371)
        self.assertEqual(ui.get_fit_results().numpoints, 4881)
        self.assertEqual(ui.get_fit_results().dof, 4877)

    @requires_fits
    @requires_xspec
    def test_pileup(self):
        self.run_thread('pileup')

        fr = ui.get_fit_results()
        self.assertEqualWithinTol(fr.statval, 53.6112, 1e-4)
        self.assertEqualWithinTol(fr.rstat, 1.44895, 1e-4)
        self.assertEqualWithinTol(fr.qval, 0.0379417, 1e-4)
        self.assertEqual(fr.numpoints, 42)
        self.assertEqual(fr.dof, 37)

        jdp = self.locals['jdp']
        self.assertEqualWithinTol(jdp.alpha.val, 0.522593, 1e-1)
        self.assertEqualWithinTol(jdp.f.val, 0.913458, 1e-2)

        abs1 = self.locals['abs1']
        self.assertEqualWithinTol(abs1.nh.val, 6.12101, 1e-2)

        power = self.locals['power']
        self.assertEqualWithinTol(power.gamma.val, 1.41887, 1e-2)
        self.assertEqualWithinTol(power.ampl.val, 0.00199457, 1e-2)

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
        self.assertEqual(len(lines), 19)
        self.assertEqual(lines[10].strip(), '')
        self.assertTrue(lines[11].startswith('   1: '))
        self.assertTrue(lines[17].startswith('   7: '))
        self.assertTrue(lines[18].startswith('   *** pileup fraction: '))

    @requires_fits
    def test_radpro(self):
        self.run_thread('radpro')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 217.450, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().rstat, 6.21287, 1e-4)
        self.assertEqualWithinTol(ui.get_fit_results().qval, 0.0, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].r0.val, 125.829, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].beta.val, 4.1633, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].xpos.val, 0.0, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ampl.val, 4.42821, 1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 92)
        self.assertEqual(ui.get_fit_results().numpoints, 38)
        self.assertEqual(ui.get_fit_results().dof, 35)

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
            self.assertEqualWithinTol(fres.statval, 217.450, 1e-4)
            self.assertEqualWithinTol(fres.rstat, 6.21287, 1e-4)
            self.assertEqualWithinTol(fres.qval, 0.0, 1e-4)

            srcmdl = self.locals['src']
            self.assertEqualWithinTol(srcmdl.r0.val, 125.829, 1e-4)
            self.assertEqualWithinTol(srcmdl.beta.val, 4.1633, 1e-4)
            self.assertEqualWithinTol(srcmdl.xpos.val, 0.0, 1e-4)
            self.assertEqualWithinTol(srcmdl.ampl.val, 4.42821, 1e-4)

            self.assertEqual(fres.nfev, 92)
            self.assertEqual(fres.numpoints, 38)
            self.assertEqual(fres.dof, 35)

    @requires_fits
    def test_psf2d(self):
        self.run_thread('psf')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 4066.78, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].fwhm.val, 2.80117, 1e-2)
        self.assertEqualWithinTol(self.locals['g1'].ypos.val, 77.2271, 1e-2)
        self.assertEqualWithinTol(self.locals['g1'].xpos.val, 88.661, 1e-2)
        self.assertEqualWithinTol(self.locals['g1'].ampl.val, 166.649, 1e-2)
        # self.assertEqual(ui.get_fit_results().nfev, 342)
        self.assertEqual(ui.get_fit_results().numpoints, 4899)
        self.assertEqual(ui.get_fit_results().dof, 4895)

    @requires_fits
    def test_fpsf2d(self):
        self.run_thread('fpsf')

        fres = ui.get_fit_results()
        self.assertEqualWithinTol(fres.statval, -4053.6635, 1e-4)

        # self.assertEqualWithinTol(self.locals['b1'].xlow.val, -4.70832, 1e-4)
        # self.assertEqualWithinTol(self.locals['b1'].xhi.val, 164.687, 1e-4)
        # self.assertEqualWithinTol(self.locals['b1'].ylow.val, 0.83626, 1e-4)
        # self.assertEqualWithinTol(self.locals['b1'].yhi.val, 142.603, 1e-4)
        # self.assertEqualWithinTol(self.locals['b1'].ampl.val, 0.956766, 1e-4)

        g1mdl = self.locals['g1']
        self.assertEqualWithinTol(g1mdl.fwhm.val, 6.420237, 1e-4)
        self.assertEqualWithinTol(g1mdl.xpos.val, 88.940712, 1e-4)
        self.assertEqualWithinTol(g1mdl.ypos.val, 76.577265, 1e-4)
        self.assertEqualWithinTol(g1mdl.ampl.val, 36344.48324, 1e-4)

        # self.assertEqual(fres.nfev, 978)
        self.assertEqual(fres.numpoints, 4899)
        self.assertEqual(fres.dof, 4895)

    @requires_fits
    def test_radpro_psf(self):
        self.run_thread('radpro_psf')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 200.949, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].r0.val, 83.0997, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].beta.val, 2.97737, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ampl.val, 5.27604, 1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 48)
        self.assertEqual(ui.get_fit_results().numpoints, 38)
        self.assertEqual(ui.get_fit_results().dof, 35)

    @requires_fits
    def test_linepro(self):
        self.run_thread('linepro')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 203.34, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].r0.val, 4.25557, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].beta.val, 0.492232, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].ampl.val, 11.8129, 1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 17)
        self.assertEqual(ui.get_fit_results().numpoints, 75)
        self.assertEqual(ui.get_fit_results().dof, 72)

    @requires_fits
    def test_kernel(self):
        self.run_thread('kernel')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 98.5793, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].r0.val, 19.2278, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].beta.val, 0.555464, 1e-4)
        self.assertEqualWithinTol(self.locals['b1'].ampl.val, 1.93706, 1e-4)
        self.assertEqual(ui.get_fit_results().nfev, 21)
        self.assertEqual(ui.get_fit_results().numpoints, 75)
        self.assertEqual(ui.get_fit_results().dof, 72)

    @requires_fits
    @requires_xspec
    def test_spectrum(self):
        self.run_thread('spectrum')

        fres = ui.get_fit_results()
        self.assertEqualWithinTol(fres.statval, 0.0496819, 1e-4)
        self.assertEqualWithinTol(self.locals['abs2'].nh.val, 1.1015, 1e-4)
        self.assertEqualWithinTol(self.locals['mek1'].kt.val, 0.841025, 1e-4)
        self.assertEqualWithinTol(self.locals['mek1'].norm.val, 0.699761, 1e-4)
        self.assertEqualWithinTol(self.locals['mek2'].kt.val, 2.35845, 1e-4)
        self.assertEqualWithinTol(self.locals['mek2'].norm.val, 1.03724, 1e-4)
        self.assertEqual(fres.numpoints, 446)
        self.assertEqual(fres.dof, 441)

    @requires_fits
    def test_histo(self):
        self.run_thread('histo')
        self.assertEqualWithinTol(ui.get_fit_results().statval, 14.7264, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].fwhm.val, 0.0232473, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].pos.val, 1.26713, 1e-4)
        self.assertEqualWithinTol(self.locals['g1'].ampl.val, 40.4503, 1e-4)
        # self.assertEqual(ui.get_fit_results().nfev, 19)
        self.assertEqual(ui.get_fit_results().numpoints, 50)
        self.assertEqual(ui.get_fit_results().dof, 47)

    @requires_fits
    @requires_xspec
    def test_xmm(self):
        self.run_thread('xmm')

        fres = ui.get_fit_results()
        self.assertEqualWithinTol(fres.statval, 118.085, 1e-4)
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
        self.assertEqualWithinTol(fres.statval, 18.8316, 1e-4)
        self.assertEqual(fres.numpoints, 46)

        aamdl = self.locals['aa']
        self.assertEqualWithinTol(aamdl.gamma.val, 1.83906, 1e-4)
        self.assertEqualWithinTol(aamdl.ampl.val, 0.000301258, 1e-4)

    @requires_fits
    @requires_xspec
    def test_proj(self):
        self.run_thread('proj')
        # Fit 1
        self.assertEqualWithinTol(self.locals['fit_res1'].statval,
                                  64.3803, 1e-2)
        self.assertEqualWithinTol(self.locals['g'].fwhm.val, 32.3536, 1e-2)
        self.assertEqualWithinTol(self.locals['g'].pos.val, 807.863, 1e-2)
        self.assertEqualWithinTol(self.locals['g'].ampl.val, 117.826, 1e-2)

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
        self.assertEqualWithinTol(self.locals['xs1'].nh.val, 0.0532586, 1e-2)
        self.assertEqualWithinTol(self.locals['p'].gamma.val, 1.4816, 1e-2)
        self.assertEqualWithinTol(self.locals['p'].ampl.val, 0.302343, 1e-2)

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
        self.run_thread('proj_bubble')

        # Fit -- Results from reminimize

        # TODO: given that some valurs << 1e-2 using atol=1e-2 is not
        # likely to be very useful here.
        #
        # The fit results change in XSPEC 12.10.0 since the mekal model
        # was changed (FORTRAN to C++). The following are from prior
        # to XSPEC 12.10.0:
        #
        kt = 17.8849
        norm = 4.15418e-6
        kt_min_covar = -0.328832
        norm_min_covar = -8.847916e-7
        kt_max_covar = 0.328832
        norm_max_covar = 8.847916e-7
        kt_min_proj = -12.048069
        norm_min_proj = -9.510913e-07
        norm_max_proj = 2.403640e-06

        # Given that the tolerance is now quite large, and not much
        # use for the normalization, should the test be updated?
        #
        # atol = 1e-2 used prior to XSPEC 12.10.0
        atol = 0.057

        mek1 = self.locals['mek1']
        assert_allclose(mek1.kt.val, kt, atol=atol, rtol=0)
        assert_allclose(mek1.norm.val, norm, atol=atol, rtol=0)

        # Covar
        #
        # TODO: should this check that parmaxes is -1 * parmins instead?
        covar = ui.get_covar_results()
        assert_allclose(covar.parmins[0], kt_min_covar, atol=atol, rtol=0)
        assert_allclose(covar.parmins[1], norm_min_covar, atol=atol, rtol=0)
        assert_allclose(covar.parmaxes[0], kt_max_covar, atol=atol, rtol=0)
        assert_allclose(covar.parmaxes[1], norm_max_covar, atol=atol, rtol=0)

        # Proj -- Upper bound of kT can't be found
        #
        proj = ui.get_proj_results()
        assert_allclose(proj.parmins[0], kt_min_proj, atol)
        assert_allclose(proj.parmins[1], norm_min_proj, atol)
        self.assertTrue(proj.parmaxes[0] is None, msg=str(proj.parmaxes[0]))
        assert_allclose(proj.parmaxes[1], norm_max_proj, atol)

    # New tests based on SDS threads -- we should catch these errors
    # (if any occur) so SDS doesn't waste time tripping over them.
    @requires_fits
    @requires_xspec
    def test_counts(self):
        self.run_thread('counts')
        self.assertEqualWithinTol(self.locals['counts_data1'], 52701.0, 1e-4)
        self.assertEqualWithinTol(self.locals['counts_data2'], 25032.0, 1e-4)
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
                                  0.04418584, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].xpos.val,
                                  150.016, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ypos.val,
                                  2.66493839, 1e-4)
        self.assertEqualWithinTol(self.locals['src'].ampl.val,
                                  1.56090546, 1e-4)
        self.assertEqualWithinTol(self.locals['bkg'].c0.val,
                                  -1.513700715, 1e-4)

        fres = ui.get_fit_results()
        self.assertEqualWithinTol(fres.istatval, 19496.3, 1e-4)
        self.assertEqualWithinTol(fres.statval, 592.32647, 1e-4)
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

        self.assertEqualWithinTol(ui.get_data().channel[0], 1.0, 1e-4)
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
        self.assertEqualWithinTol(rmf.energ_lo[0], EMIN, etol)
        self.assertEqualWithinTol(arf.energ_lo[0], EMIN, etol)
