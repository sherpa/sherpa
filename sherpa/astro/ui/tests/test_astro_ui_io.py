#
#  Copyright (C) 2015  Smithsonian Astrophysical Observatory
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

# break out some of the I/O tests in the UI layer into a separate set
# of tests, in part to reduce file size, but also because some of these
# may be better placed in tests of the sherpa.astro.io module, once that
# becomes possible

import unittest

from sherpa.utils import SherpaTest, SherpaTestCase
from sherpa.utils import requires_data, requires_fits
from sherpa.astro import ui
from sherpa.astro.data import DataPHA
from sherpa.astro.instrument import ARF1D, RMF1D

import logging
logger = logging.getLogger("sherpa")


@requires_fits
@requires_data
class test_load_pha3_gzip(SherpaTestCase):
    """Handle a .gz FITS PHA Type 3 file"""

    longMessage = True

    def setUp(self):
        # hide warning messages from file I/O
        self._old_logger_level = logger.level
        logger.setLevel(logging.ERROR)

        ui.clean()

        self.head = self.make_path('acisf01575_001N001_r0085')

    def tearDown(self):
        logger.setLevel(self._old_logger_level)

    def validate_pha(self, idval):
        """Check that the PHA dataset in id=idval is
        as expected.
        """

        self.assertEqual(ui.list_data_ids(), [idval])

        pha = ui.get_data(idval)
        self.assertIsInstance(pha, DataPHA)

        arf = ui.get_arf(idval)
        self.assertIsInstance(arf, ARF1D)

        rmf = ui.get_rmf(idval)
        self.assertIsInstance(rmf, RMF1D)

        bpha = ui.get_bkg(idval, bkg_id=1)
        self.assertIsInstance(bpha, DataPHA)

        barf = ui.get_arf(idval, bkg_id=1)
        self.assertIsInstance(barf, ARF1D)

        brmf = ui.get_rmf(idval, bkg_id=1)
        self.assertIsInstance(brmf, RMF1D)

        # normally the background data set would have a different name,
        # but this is a  PHA Type 3 file.
        # self.assertEqual(pha.name, bpha.name)
        self.assertEqual(arf.name, barf.name)
        self.assertEqual(rmf.name, brmf.name)

    def testReadExplicit(self):
        """Include .gz in the file name"""

        idval = 12
        fname = self.head + '_pha3.fits.gz'
        ui.load_pha(idval, fname)

        self.validate_pha(idval)

        # TODO: does this indicate that the file name, as read in,
        #       should have the .gz added to it to match the data
        #       read in, or left as is?
        pha = ui.get_data(idval)
        bpha = ui.get_bkg(idval, bkg_id=1)
        self.assertEqual(pha.name, bpha.name + '.gz')

    def testReadImplicit(self):
        """Exclude .gz from the file name"""

        idval = "13"
        fname = self.head + '_pha3.fits'
        ui.load_pha(idval, fname)

        self.validate_pha(idval)

        pha = ui.get_data(idval)
        bpha = ui.get_bkg(idval, bkg_id=1)
        self.assertEqual(pha.name, bpha.name)

if __name__ == '__main__':

    import sys
    if len(sys.argv) > 1:
        datadir = sys.argv[1]
    else:
        datadir = None

    SherpaTest(ui).test(datadir=datadir)
