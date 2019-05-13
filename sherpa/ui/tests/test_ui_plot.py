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
from sherpa.utils.testing import requires_plotting

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


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("get_prefs", [ui.get_data_plot_prefs,
                                       ui.get_model_plot_prefs])
def test_plot_prefs_xxx(get_prefs):
    """Can we change and reset a preference.

    Pick the 'xlog' field, since the assumption is that: a) this
    defaults to 'False'; b) each plot type has this setting; c) we
    do not need to check all settings.

    Some tests fail due to missing plot preferences when there's
    no plotting backend (e.g. missing 'xlog' settings), so skip
    these tests in this case.
    """

    prefs1 = get_prefs()
    assert not prefs1['xlog']
    prefs1['xlog'] = True

    prefs2 = get_prefs()
    assert prefs2['xlog']

    ui.clean()
    prefs3 = get_prefs()
    assert prefs1['xlog']
    assert prefs2['xlog']
    assert not prefs3['xlog']


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("pfunc", [ui.plot_data,
                                   ui.plot_model,
                                   ui.plot_source,
                                   ui.plot_resid,
                                   ui.plot_ratio,
                                   ui.plot_delchi,
                                   ui.plot_fit,
                                   ui.plot_fit_resid,
                                   ui.plot_fit_delchi])
def test_fit_plot_xxx(idval, pfunc):
    """Can we call a plot_xxx routine?

    Currently this just tests that the call succeeds. There is no
    test to see if the plot did anything.

    Some tests fail due to missing plot preferences when there's
    no plotting backend (e.g. missing 'xlog' settings), so skip
    these tests in this case.
    """

    setup_example(idval)
    if idval is None:
        pfunc()
    else:
        pfunc(idval)
