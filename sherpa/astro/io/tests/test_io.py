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

from sherpa.utils import SherpaTestCase, has_fits_support, has_package_from_list
from unittest import skipIf

from sherpa.astro import ui


class test_89_issues(SherpaTestCase):
    @skipIf(not has_fits_support() or not has_package_from_list("sherpa.astro.xspec"), "fits support required")
    def test_mod_fits(self):
        ui.clean()
        tablemodelfile = self.make_path("xspec", "tablemodel", "RCS.mod")
        ui.load_table_model("tmod", tablemodelfile)
        tmod = ui.get_model_component("tmod")
        self.assertEqual("xstablemodel.tmod", tmod.name)
