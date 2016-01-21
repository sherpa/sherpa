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

from sherpa.utils.test import SherpaTestCase, requires_data, \
    requires_fits, requires_xspec
from sherpa.astro import ui

from tempfile import NamedTemporaryFile
import warnings


class test_89_issues(SherpaTestCase):
    def setUp(self):
        ui.clean()

    @requires_data
    @requires_fits
    @requires_xspec
    def test_mod_fits(self):
        tablemodelfile = self.make_path("xspec-tablemodel-RCS.mod")
        ui.load_table_model("tmod", tablemodelfile)
        tmod = ui.get_model_component("tmod")
        self.assertEqual("xstablemodel.tmod", tmod.name)

    @requires_fits
    def test_warnings_are_gone_arrays(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with NamedTemporaryFile() as f:
                ui.load_arrays(1, [1,2,3], [4,5,6])
                ui.save_data(1, f.name, ascii=True, clobber=True)
            with NamedTemporaryFile() as f:
                ui.save_data(1, f.name, ascii=False, clobber=True)
            assert len(w) == 0

    @requires_fits
    @requires_data
    def test_warnings_are_gone_pha(self):
        pha = self.make_path("3c273.pi")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with NamedTemporaryFile() as f:
                ui.load_pha(pha)
                ui.save_data(1, f.name, ascii=False, clobber=True)
            assert len(w) == 0
