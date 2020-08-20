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

import numpy as np
import logging
import pytest

from sherpa.testing import SherpaTestCase
from sherpa import ui

logger = logging.getLogger('sherpa')

class test_get_int_proj(SherpaTestCase):

    _get_int_uncproj_recalc_bench = {
        'x': np.array([-5., -4.75, -4.5, -4.25, -4., -3.75, -3.5, -3.25, -3.,
                      -2.75, -2.5, -2.25, -2., -1.75, -1.5, -1.25, -1., -0.75,
                       -0.5, -0.25, 0.]),
        'proj': np.array([110.0185, 90.493 , 72.9534, 57.3995, 43.8316,
                          32.2494, 22.6531, 15.0427, 9.4181, 5.7794, 4.1265,
                          4.4594, 6.7782, 11.0829, 17.3734, 25.6497, 35.9119,
                          48.1599, 62.3938, 78.6135, 96.8191]),
        'unc': np.array([110.774, 91.1094, 73.4447, 57.78, 44.1153, 32.4507,
                         22.786, 15.1213, 9.4566, 5.7919, 4.1273, 4.4626,
                         6.7979, 11.1332, 17.4686, 25.8039, 36.1392, 48.4745,
                         62.8099, 79.1452, 97.4805]),
        'min': -5,
        'max': 0,
        'nloop': 21,
        'delv': None,
        'fac': 1,
        'log': False
    }

    def setUp(self):
        # defensive programming (one of the tests has been seen to fail
        # when the whole test suite is run without this)
        ui.clean()
        self._old_logger_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        x = [-13, -5, -3, 2, 7, 12];
        y = [102.3, 16.7, -0.6, -6.7, -9.9, 33.2];
        dy = np.ones(6) * 5
        ui.load_arrays(1, x, y, dy)
        ui.set_source(ui.polynom1d.poly)
        poly.c1.thaw()
        poly.c2.thaw()
        ui.int_proj(poly.c0)
        ui.fit()

    def tearDown(self):
        ui.clean()

        try:
            logger.setLevel(self._old_logger_level)
        except AttributeError:
            pass

    def cmp_values(self, arg1, arg2):
        for a, b in zip(arg1, arg2):
            assert a == pytest.approx(b, abs=1.0e-3)

    def cmp_args(self, result):
        assert self._get_int_uncproj_recalc_bench["min"] == result.min
        assert self._get_int_uncproj_recalc_bench["max"] == result.max
        assert self._get_int_uncproj_recalc_bench["nloop"] == result.nloop
        assert self._get_int_uncproj_recalc_bench["delv"] == result.delv
        assert self._get_int_uncproj_recalc_bench["fac"] == result.fac
        assert self._get_int_uncproj_recalc_bench["log"] == result.log

    def test_get_int_proj(self):
        result = ui.get_int_proj(poly.c1, recalc=True, min=-5, max=0, nloop=21)
        self.cmp_args(result)
        self.cmp_values(self._get_int_uncproj_recalc_bench["x"], result.x)
        self.cmp_values(self._get_int_uncproj_recalc_bench["proj"], result.y)

    def test_get_int_unc(self):
        result = ui.get_int_unc(poly.c1, recalc=True, min=-5, max=0, nloop=21)
        self.cmp_args(result)
        self.cmp_values(self._get_int_uncproj_recalc_bench["x"], result.x)
        self.cmp_values(self._get_int_uncproj_recalc_bench["unc"], result.y)
