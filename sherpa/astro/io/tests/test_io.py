#
#  Copyright (C) 2016, 2017  Smithsonian Astrophysical Observatory
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

from sherpa.utils import SherpaTestCase, requires_data, requires_fits, \
    requires_xspec
from sherpa.astro import ui

from tempfile import NamedTemporaryFile
import warnings

try:
    from astropy.utils.exceptions import AstropyDeprecationWarning
    have_astropy = True
except ImportError:
    have_astropy = False


def strip_astropy_clobber_deprecation(ws):
    """Remove the AstroPy deprecation warning from the warning list.

    Parameters
    ----------
    ws : seq of warnings
        The warnings thrown by the code

    Returns
    -------
    ws : list of warnings
        A copy of the input sequence, excluding any that are
        an Astropy Deprecation warning for the clobber parameter
        introduced in AstroPy 1.3.0.
    """

    wmsg = '"clobber" was deprecated in version 1.3 and will be ' + \
           'removed in a future version. Use argument "overwrite" instead.'

    def keep(w):
        # it is not clear how stable _category_name is, but it is
        # simpler than checking against w.category
        if w._category_name != 'AstropyDeprecationWarning':
            return True
        return str(w.message) != wmsg

    return [w for w in ws if keep(w)]


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
        ui.load_arrays(1, [1, 2, 3], [4, 5, 6])
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with NamedTemporaryFile() as f:
                ui.save_data(1, f.name, ascii=True, clobber=True)
            with NamedTemporaryFile() as f:
                ui.save_data(1, f.name, ascii=False, clobber=True)
            w = strip_astropy_clobber_deprecation(w)
            assert len(w) == 0

    @requires_fits
    @requires_data
    def test_warnings_are_gone_pha(self):
        pha = self.make_path("3c273.pi")
        ui.load_pha(pha)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            with NamedTemporaryFile() as f:
                ui.save_data(1, f.name, ascii=False, clobber=True)
            w = strip_astropy_clobber_deprecation(w)
            assert len(w) == 0
