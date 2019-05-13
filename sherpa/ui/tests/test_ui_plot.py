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
from sherpa.stats import Chi2Gehrels
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

    assert dp.x == pytest.approx(mp.x)
    assert dp.y == pytest.approx([10, 40, 30, 50])
    assert mp.y == pytest.approx([35, 35, 35, 35])

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


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_plot_data_change(idval):
    """Can we plot, change data, plot and see the difference?

    This relies on checking the plot structure (returned by get_data_plot)
    to approximate what the plot would be like.
    """

    setup_example(idval)
    if idval is None:
        ui.plot_data()
        pvals1 = ui.get_data_plot()
        dset = ui.get_data()
    else:
        ui.plot_data(idval)
        pvals1 = ui.get_data_plot(idval)
        dset = ui.get_data(idval)

    # do not test the plot_prefs field
    yold = [10, 40, 30, 50]
    assert pvals1.xlabel == 'x'
    assert pvals1.ylabel == 'y'
    assert pvals1.title == 'example'
    assert pvals1.x == pytest.approx([10, 20, 40, 90])
    assert pvals1.y == pytest.approx(yold)
    assert pvals1.xerr is None

    assert pvals1.yerr == pytest.approx(Chi2Gehrels.calc_staterror(yold))

    # Modify the data values; rely on changing the dataset object
    # directly means that we do not need to call set_data here.
    #
    ynew = [12, 45, 33, 49]
    dset.y = ynew

    # Check the new value (could just check pvals1, but make the
    # explicit call to get_data_plot)
    #
    if idval is None:
        ui.plot_data()
        pvals2 = ui.get_data_plot(idval)
    else:
        ui.plot_data(idval)
        pvals2 = ui.get_data_plot(idval)

    assert pvals2.xlabel == 'x'
    assert pvals2.ylabel == 'y'
    assert pvals2.title == 'example'
    assert pvals2.x == pytest.approx([10, 20, 40, 90])
    assert pvals2.y == pytest.approx(ynew)
    assert pvals2.xerr is None

    # Should use approximate equality here
    assert pvals2.yerr == pytest.approx(Chi2Gehrels.calc_staterror(ynew))

    # just check that the previous value has been updated too
    assert pvals1.y == pytest.approx(pvals2.y)


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_plot_data_replot(idval):
    """Can we plot, change data, plot with reploat and see no difference?

    This relies on accessing the undelying session object directly,
    as ui.get_data_plot always recreates the plot objects, which
    isn't helpful here.

    """

    setup_example(idval)
    if idval is None:
        ui.plot_data()
        pvals1 = ui.get_data_plot()
        dset = ui.get_data()
    else:
        ui.plot_data(idval)
        pvals1 = ui.get_data_plot(idval)
        dset = ui.get_data(idval)

    # the fields returned by get_data_plot have already been tested
    # by test_plot_data_change, so no need to repeat this here.
    #

    # Modify the data values; rely on changing the dataset object
    # directly means that we do not need to call set_data here.
    #
    yold = [10, 40, 30, 50]
    ynew = [12, 45, 33, 49]
    dset.y = ynew

    # Check the new value
    #
    if idval is None:
        ui.plot_data(replot=True)
    else:
        ui.plot_data(idval, replot=True)

    # access the data-plot object directly, since get_data_plot
    # would cause this to be changed.
    pvals2 = ui._session._dataplot

    assert pvals2.xlabel == 'x'
    assert pvals2.ylabel == 'y'
    assert pvals2.title == 'example'
    assert pvals2.x == pytest.approx([10, 20, 40, 90])
    assert pvals2.y == pytest.approx(yold)
    assert pvals2.xerr is None

    # Should use approximate equality here
    assert pvals2.yerr == pytest.approx(Chi2Gehrels.calc_staterror(yold))
