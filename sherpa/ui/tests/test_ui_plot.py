#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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

"""
Basic tests of the plot functionality in sherpa.ui. The idea is that
these - or most of them - can be run even when there is no plot backend.
"""

from sherpa import ui
from sherpa.plot import DataPlot, FitPlot, ModelPlot

import pytest


def example_data():
    """Create an example data set."""

    return ui.Data1D('example', [10, 20, 40, 90], [10, 40, 30, 50])

def example_model():
    """Create an example model."""

    ui.create_model_component('const1d', 'cpt')
    cpt = ui.get_model_component('cpt')
    cpt.c0 = 35
    return cpt


def setup_example(idval):
    """Set up a simple dataset for use in the tests.

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier.
    """

    d = example_data()
    m = example_model()
    if idval is None:
        ui.set_data(d)
        ui.set_source(m)

    else:
        ui.set_data(idval, d)
        ui.set_source(idval, m)


@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_fit_plot(idval):
    """Basic testing of get_fit_plot
    """

    setup_example(idval)
    if idval is None:
        f = ui.get_fit_plot()
    else:
        f = ui.get_fit_plot(idval)

    # Should we be checking exact class?
    #
    # No - because the mock_chips pytest routine redefines these
    # classes and so makes the test fail (or at least that's my
    # assumption).
    #
    # assert isinstance(f, FitPlot)
    # assert isinstance(f.dataplot, DataPlot)
    # assert isinstance(f.modelplot, ModelPlot)
    
    dp = f.dataplot
    mp = f.modelplot

    # Are we guaranteed that arrays are Numpy arrays?
    def are_equal(a, b):
        r = a == b
        try:
            return all(r)
        except TypeError:
            return r
        
    assert are_equal(dp.x, mp.x)
    assert are_equal(dp.y, [10, 40, 30, 50])
    assert are_equal(mp.y, [35, 35, 35, 35])

    assert dp.xerr is None
    assert mp.xerr is None
    assert mp.yerr is None

    assert dp.title == 'example'
    assert mp.title == 'Model'


@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_fit_plot(idval):
    """How bad do things fail on Travis with this?

    The idea is just to see what the failure modes on Travis are
    for the "no-plot-backend" cases, so there's no actual test
    to see if the plot did anything.
    """

    setup_example(idval)
    if idval is None:
        ui.plot_fit()
    else:
        ui.plot_fit(idval)
