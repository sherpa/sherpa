#
#  Copyright (C) 2024
#  MIT
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
"""This module collects tests specific to bokeh."""

import os
import pytest
from sherpa import plot
from sherpa.astro import ui
from sherpa.utils.testing import requires_data, requires_fits


@requires_fits
@requires_data
def test_bokeh_delchi(caplog, clean_ui):
    """Test a specific multi-panel plot with the bokeh backend.

    This also implicitly tests setting a line_color, because
    delchi put a black horizontal line on the plot.

    This is a regression test for https://github.com/sherpa/sherpa/issues/1975
    """
    path = os.path.join(
        os.path.dirname(__file__),
        "..",
        "..",
        "astro",
        "datastack",
        "tests",
        "data",
        "3c273.pi",
    )
    ui.load_pha(path)
    ui.set_source(1, ui.powlaw1d(name="powlaw"))

    bokehbackend = pytest.importorskip("sherpa.plot.bokeh_backend")
    newback = bokehbackend.BokehBackend()
    with plot.TemporaryPlottingBackend(newback):
        ui.plot_fit_delchi()
