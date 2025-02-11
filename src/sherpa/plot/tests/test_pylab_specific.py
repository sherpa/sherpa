#
#  Copyright (C) 2022, 2023
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
'''This module collects tests specific to matplotlib.

The purpose here is not to throw all tests that work with matplotlib into
a single file, but instead to provide a space for tests that do not so much
test sherpa plotting functionality, but that just look at specific
implementation details in pylab_backend.py.

That is not always a clear line, so as a rule of thumb, existing tests
remain in their old place, only newly written tests will be placed here.
'''

import numpy as np
import pytest
from sherpa import data
from sherpa import plot


def test_PylabErrorArea(caplog):
    '''Test that PylabErrorArea is functional.

    This does not test that its function is useful, only that it works at all.
    This test essentially duplicates an example from the documentation and can
    be removed if the code in the documentation is tested in CI in the future.
    '''
    x1 = np.asarray([100, 200, 600, 1200])
    y1 = np.asarray([2000, 2100, 1400, 3050])

    dy1 = np.asarray([100, 50, 200, 300])
    d2 = data.Data1D('errors', x1, y1, dy1)

    # When we have more tests, it might be better to check it matplotlib is
    # available at the module level, but for now this is easier.
    mpl = pytest.importorskip("matplotlib")
    from sherpa.plot import pylab_area_backend

    newback = pylab_area_backend.PylabErrorArea()
    with plot.TemporaryPlottingBackend(newback):
        plot2 = plot.DataPlot()
        plot2.prepare(d2)
        plot2.plot()

    # Check we don't see any warnings from changing the backend
    assert len(caplog.records) == 0
