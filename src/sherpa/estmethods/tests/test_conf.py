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

import numpy as np
import logging

from sherpa.fit import Fit
from sherpa.data import Data1D
from sherpa.models.basic import Polynom1D
from sherpa.estmethods import Confidence
from sherpa.utils.testing import SherpaTestCase
from sherpa import ui

logger = logging.getLogger("sherpa")

class test_estmethods(SherpaTestCase):

    _x = [-13, -5, -3, 2, 7, 12]
    _y = [102.3, 16.7, -0.6, -6.7, -9.9, 33.2]
    _e = np.ones(6) * 5

    _conf_bench = {
        'parnames' : ('mdl.c0', 'mdl.c1', 'mdl.c2'),
        'parvals' : np.array([-9.384889507344322, -2.4169154937347925,
                              0.47827334261018023]),
        'parmins' : np.array([-2.917507940156449, -0.250889317129555,
                              -0.03126766429871808]),
        'parmaxes' : np.array([2.91750794015645, 0.250889317129555,
                               0.03126766429871808])
        }
    
    def setUp(self):
        # defensive programming (one of the tests has been seen to fail
        # when the whole test suite is run without this)
        ui.clean()
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)
        self.data = Data1D('tst', self._x, self._y, self._e)
        self.mdl = Polynom1D('mdl')
        
    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def cmp_results(self, result, tol=1.0e-3):
        for ii in range(len(self._conf_bench['parnames'])):
            assert self._conf_bench['parnames'][ii] == result.parnames[ii]
        self.assertEqualWithinTol(self._conf_bench['parvals'],
                                  result.parvals, tol)
        self.assertEqualWithinTol(self._conf_bench['parmins'],
                                  result.parmins, tol)
        self.assertEqualWithinTol(self._conf_bench['parmaxes'],
                                  result.parmaxes, tol)
        
    def tst_low_level(self, thaw_c1):
        if thaw_c1:
            self.mdl.c1.thaw()
        self.mdl.c2.thaw()
        f = Fit(self.data, self.mdl, estmethod=Confidence())
        self.mdl.c2 = 1
        f.fit()
        if not thaw_c1:
            self.mdl.c1.thaw()
            f.fit()
        result = f.est_errors()
        self.cmp_results(result)
        
    def tst_ui(self, thaw_c1):
        ui.load_arrays(1, self._x, self._y, self._e)
        ui.set_source(1, ui.polynom1d.mdl)
        if  thaw_c1:
            ui.thaw(mdl.c1)
        ui.thaw(mdl.c2)
        mdl.c2 = 1
        ui.fit()
        if not thaw_c1:
            ui.thaw(mdl.c1)
            ui.fit()
        ui.conf()
        result = ui.get_conf_results()
        self.cmp_results(result)
        
    def test_low_level_True(self):
        self.tst_low_level(True)

    def test_low_level_False(self):
        self.tst_low_level(False)

    def test_ui_True(self):
        self.tst_ui(True)

    def test_ui_False(self):
        self.tst_ui(False)        
