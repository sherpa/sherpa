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
from sherpa.utils import SherpaTest, SherpaTestCase
from sherpa.utils import requires_data, requires_fits
from sherpa.astro import ui
import logging
import os
import numpy
logger = logging.getLogger("sherpa")


@requires_data
@requires_fits
class test_more_ui(SherpaTestCase):
    def assign_model(self, name, obj):
        self.locals[name] = obj

    def run_thread(self, name, scriptname='fit.py'):
        ui.clean()
        ui.set_model_autoassign_func(self.assign_model)
        self.locals = {}
        cwd = os.getcwd()
        os.chdir(self.make_path('ciao4.3', name))
        try:
            execfile(scriptname, {}, self.locals)
        finally:
            os.chdir(cwd)

    def setUp(self):
        self.img = self.make_path('img.fits')
        self.pha = self.make_path('threads/simultaneous/pi2286.fits')
        self.rmf = self.make_path('threads/simultaneous/rmf2286.fits')
        self.nan = self.make_path('ciao4.3/filternan/with_nan.fits')
        self.loggingLevel = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

    def tearDown(self):
        if hasattr(self, 'loggingLevel'):
            logger.setLevel(self.loggingLevel)

    # bug 12784
    def test_filter_nan(self):
        self.run_thread('filternan')
        self.assertFalse(numpy.isnan(ui.get_fit_results().statval))

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
