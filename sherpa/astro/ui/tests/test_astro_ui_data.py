#
#  Copyright (C) 2025
#  Smithsonian Astrophysical Observatory
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
import logging
import numpy as np
import pytest
from tempfile import NamedTemporaryFile

from sherpa.astro.io.wcs import WCS
from sherpa.astro import ui
from sherpa.utils.testing import requires_fits


@requires_fits
def test_DataIMG_filter(caplog, clean_astro_ui):
    """Test filtering of DataIMG

    This is a slightly modified version of https://github.com/sherpa/sherpa/issues/1880
    """
    y = np.arange(6).reshape(2, 3) * 10 + 10

    sky = WCS("physical", "LINEAR", [100, 200], [1, 1], [10, 10])
    eqpos = WCS("world", "WCS", [30, 50], [100, 200], [-0.1, 0.1])
    data = ui.DataIMG.from_2d_array_with_wcs("faked", y=y,
                    sky=sky, eqpos=eqpos)

    ui.set_data(data)
    ui.set_coord("physical")

    # this should, based on DS9, remove the 50 bin
    ui.ignore2d("circle(110, 210, 6)")

    assert ui.get_dep(filter=True) == pytest.approx([10, 40, 20, 30, 60])

    # switch back to logical for the independent axis
    ui.set_coord("logical")

    g0, g1 = ui.get_indep()
    assert g0 == pytest.approx([1 ,2, 1, 2, 1, 2])
    assert g1 == pytest.approx([1, 1, 2, 2, 3, 3])

    with NamedTemporaryFile(suffix=".fits") as tmpf:
        with caplog.at_level(logging.WARNING):
            ui.save_image(tmpf.name, ascii=False, clobber=True)
        assert "Region filter has been removed from 'faked'" in caplog.text
        ui.load_image("copy", tmpf.name)

    # Filter is not saved
    assert ui.get_dep("copy", filter=False) == pytest.approx([10, 40, 20, 50, 30, 60])
    assert ui.get_dep("copy", filter=True) == pytest.approx([10, 40, 20, 50, 30, 60])

    h0, h1 = ui.get_indep("copy")
    assert h0 == pytest.approx([1, 2, 1, 2, 1, 2])
    assert h1 == pytest.approx([1, 1, 2, 2, 3, 3])

    ui.set_coord("copy", "physical")
    ui.ignore2d_id("copy", "circle(110, 210, 6)")
    assert ui.get_dep("copy", filter=True) == pytest.approx([10, 40, 20, 30, 60])