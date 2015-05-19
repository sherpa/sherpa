# 
#  Copyright (C) 2013, 2015  Smithsonian Astrophysical Observatory
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



import unittest
from sherpa.utils import SherpaTest, SherpaTestCase, test_data_missing
import sherpa.astro.ui as ui
import logging
import os
import numpy
logger = logging.getLogger("sherpa")

class test_more_ui(SherpaTestCase):
    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        self.locals = {}
        os.chdir(os.path.join(self.datadir, 'ciao4.3', name))
        execfile(scriptname, {}, self.locals)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.img = self.datadir + '/img.fits'
        self.pha = self.datadir + '/threads/simultaneous/pi2286.fits'
        self.rmf = self.datadir + '/threads/simultaneous/rmf2286.fits'
        self.nan = self.datadir + '/ciao4.3/filternan/with_nan.fits'
        logger.setLevel(logging.ERROR)

    #bug 12784
    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_filter_nan(self):
        self.run_thread('filternan')
        self.assertFalse(numpy.isnan(ui.get_fit_results().statval))

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        SherpaTest(ui).test(datadir=sys.argv[1])
