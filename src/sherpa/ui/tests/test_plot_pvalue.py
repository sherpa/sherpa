#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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

from numpy import testing

from sherpa.testing import SherpaTestCase, requires_data, \
    requires_xspec, requires_fits
import logging
logger = logging.getLogger("sherpa")
from sherpa.astro import ui

@requires_fits
@requires_xspec
@requires_data
class test_plot_pvalue(SherpaTestCase):
    def setUp(self):
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        if hasattr(self, "_old_logger_level"):
            logger.setLevel(self._old_logger_level)

    def compare_results(self, arg1, arg2, atol=1.0e-4, rtol=1.0e-4):
        for a, b in zip(arg1, arg2):
            testing.assert_allclose(a, b, rtol=rtol, atol=atol)

    def test_rsp(self):
        fname = self.make_path("qso.pi")
        ui.load_pha(fname)
        ui.set_stat("chi2xspecvar")
        ui.set_method("neldermead")
        ui.group_counts(10)
        ui.notice(0.3,8)
        ui.set_model("xsphabs.abs1*xspowerlaw.p1")
        ui.set_model("abs1*(p1+gauss1d.g1)")
        g1.pos = 3.
        ui.freeze(g1.pos)
        g1.fwhm=0.1
        ui.freeze(g1.fwhm)
        ui.set_stat('cstat')
        ui.fit()
        ui.plot_pvalue(p1,p1+g1, num=100)
        tmp = ui.get_pvalue_results()
        expected = [210.34566845619273, 207.66618095925094, 2.679487496941789]
        self.compare_results(expected, [tmp.null, tmp.alt, tmp.lr])
