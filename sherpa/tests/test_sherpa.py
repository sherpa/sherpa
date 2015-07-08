# 
#  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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
import os.path
import sherpa
from sherpa.utils import SherpaTest, SherpaTestCase, test_data_missing
from sherpa import ui

class test_sherpa(SherpaTestCase):

    def test_include_dir(self):
        incdir = os.path.join(sherpa.get_include(), 'sherpa')
        self.assert_(os.path.isdir(incdir))

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def setUp(self):
        self.agn2 = self.datadir + '/ciao4.3/faulty_load_data/agn2'
        self.agn2_fixed = self.datadir + '/ciao4.3/faulty_load_data/agn2_fixed'
        self.template_idx = self.datadir + '/ciao4.3/faulty_load_data/table.txt'

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_not_reading_header_without_comment(self):
        self.assertRaises(ValueError, ui.load_data, self.agn2)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_reading_floats(self):
        ui.load_data(self.agn2_fixed)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_reading_strings(self):
        ui.load_data(self.template_idx, require_floats=False)

    @unittest.skipIf(test_data_missing(), "required test data missing")
    def test_require_float(self):
        self.assertRaises(ValueError, ui.load_data, self.agn2)
