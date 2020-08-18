#
#  Copyright (C) 2007, 2015, 2016, 2018  Smithsonian Astrophysical Observatory
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

import os.path
import sherpa
from sherpa.utils.testing import SherpaTestCase, requires_data
from sherpa import ui


@requires_data
class test_sherpa(SherpaTestCase):

    def test_include_dir(self):
        incdir = os.path.join(sherpa.get_include(), 'sherpa')
        self.assertTrue(os.path.isdir(incdir))

    def setUp(self):
        self.agn2 = self.make_path('agn2')
        self.agn2_fixed = self.make_path('agn2_fixed')
        self.template_idx = self.make_path('table.txt')

    def test_not_reading_header_without_comment(self):
        self.assertRaises(ValueError, ui.load_data, self.agn2)

    def test_reading_floats(self):
        ui.load_data(self.agn2_fixed)

    def test_reading_strings(self):
        ui.load_data(self.template_idx, require_floats=False)

    def test_require_float(self):
        self.assertRaises(ValueError, ui.load_data, self.agn2)
