#
#  Copyright (C) 2021  MIT
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
from contextlib import contextmanager
import json

import numpy as np

import pytest

from sherpa import data
from sherpa import plot
from sherpa.utils.testing import requires_bokeh



# This might go into general utils, but for now, I'm developing it here
@contextmanager
def plottingbackend(backend):
    oldbackend = plot.backend
    plot.backend = backend
    try:
        yield
    finally:
        plot.backend = backend

        
@pytest.fixture
def use_bokeh():
    from sherpa.plot import bokeh_backend
    oldbackend = plot.backend
    plot.backend = bokeh_backend
    try:
        yield
    finally:
        plot.backend = oldbackend

        
# short explanation why those are separate fixtures.
# or make them one?
@requires_bokeh
# Need decorator or with statement to set plotting backend here
def test_basic_bokeh_plot(use_bokeh):
    '''At this point, this test is just a template on how to write tests
    for bokeh output'''

    # Imports here so that the module imports correctly
    from bokeh.embed import json_item
    from sherpa.plot import backend
    
    x1 = np.asarray([100, 200, 600, 1200])
    y1 = np.asarray([2000, 2100, 1400, 3050])

    d1 = data.Data1D('somelongname', x1, y1)

    plot1 = plot.DataPlot()
    plot1.prepare(d1)
    plot1.plot()

    j = json_item(backend.current_fig['current_axes'])
    out =  json.dumps(json_item(backend.current_fig['current_axes']))

    assert '"axis_label": "y"' in out
    assert '"x": [100, 200, 600, 1200]' in out
    assert 'somelongname' in out
