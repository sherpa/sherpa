#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

from sherpa.utils import SherpaTestCase
from sherpa.astro.utils._region import region_mask, Region
from sherpa.astro import ui
from tempfile import NamedTemporaryFile


class TestRegion(SherpaTestCase):
    def setUp(self):
        ui.clean()
        self.region_string = 'Circle(256 ,256 , 25)'
        self.normalized_string = 'Circle(256,256,25)'
        self.x = [1, 2, 3]
        self.y = [1, 2, 3]

    def test_region_string(self):
        r = Region(self.region_string)
        self.assertEqual(self.normalized_string, str(r))

    def test_region_mask_string(self):
        r, m = region_mask(None, self.region_string, self.x, self.y, False, False)
        self.assertEqual(self.normalized_string, str(r))
        assert all(val == 0 for val in m)

    def test_region_mask_file(self):
        with NamedTemporaryFile(bufsize=0) as f:
            f.write(self.normalized_string)
            r, m = region_mask(None, f.name, self.x, self.y, True, False)
        self.assertEqual(self.normalized_string, str(r))
        assert all(val == 0 for val in m)
