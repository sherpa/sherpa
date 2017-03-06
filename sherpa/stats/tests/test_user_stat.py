#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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
import numpy as np
import pytest

from sherpa.astro import ui
from sherpa.utils.err import StatErr


def test_341():
    def calc_stat(data, _model):
        return 3.235, np.ones_like(data.get_y())

    xdata = [1, 2, 3]
    ydata = xdata

    ui.load_arrays(1, xdata, ydata)

    ui.set_model(1, 'polynom1d.p')

    ui.load_user_stat('customstat', calc_stat, lambda x: np.ones_like(x))
    ui.set_stat(customstat)

    try:
        ui.fit(1)
    except StatErr:
        pytest.fail("Call should not be throwing and exception (bug #341)")

    # Test the result is what we made the user stat return
    assert 3.235 == ui.get_fit_results().statval
