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


_data_x = [10, 20, 40, 90]
_data_y = [10, 40, 30, 50]
      
def example_data():
    """Create an example data set."""

    return ui.Data1D('example', _data_x, _data_y)

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


def calc_errors(ys):
    """Return errors for ys using the default statistic.

    Consolidate this code to make it easier to change if the
    default statistic ever changes.
    """

    return Chi2Gehrels.calc_staterror(ys)


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

    assert pvals1.yerr == pytest.approx(calc_errors(yold))

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
    assert pvals2.yerr == pytest.approx(calc_errors(ynew))

    # just check that the previous value has been updated too
    assert pvals1.y == pytest.approx(pvals2.y)


def change_example(idval):
    """Change the example y values (created by setup_example)"""

    d = ui.get_data(idval)
    d.y = [12, 45, 33, 49]


def change_model(idval):
    """Change the example model values (created by setup_model)"""

    cpt = ui.get_model_component('cpt')
    cpt.c0 = 41


def change_fit(idval):
    """Change both the model and the data."""

    change_example(idval)
    change_model(idval)


def check_example(xlabel='x'):
    """Check that the data plot has not changed"""

    dplot = ui._session._dataplot

    assert dplot.xlabel == xlabel
    assert dplot.ylabel == 'y'
    assert dplot.title == 'example'
    assert dplot.x == pytest.approx(_data_x)
    assert dplot.y == pytest.approx(_data_y)
    assert dplot.xerr is None

    # Should use approximate equality here
    assert dplot.yerr == pytest.approx(calc_errors(_data_y))

    
def check_model_plot(plot, title='Model', xlabel='x'):
    """Helper for check_model/source"""

    assert plot.xlabel == xlabel
    assert plot.ylabel == 'y'
    assert plot.title == title
    assert plot.x == pytest.approx(_data_x)
    assert plot.y == pytest.approx([35 for x in _data_x])
    assert plot.xerr is None
    assert plot.yerr is None


def check_model(xlabel='x'):
    """Check that the model plot has not changed"""

    check_model_plot(ui._session._modelplot,
                     title='Model', xlabel=xlabel)


def check_source():
    """Check that the source plot has not changed"""

    check_model_plot(ui._session._sourceplot,
                     title='Source')

    
def check_resid(title='Residuals for example'):
    """Check that the resid plot has not changed"""

    rplot = ui._session._residplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 35 for y in _data_y])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y))

    
def check_ratio():
    """Check that the ratio plot has not changed"""

    rplot = ui._session._ratioplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == 'Ratio of Data to Model for example'
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 35 for y in _data_y])
    assert rplot.xerr is None
    dy = [dy / 35 for dy in calc_errors(_data_y)]
    assert rplot.yerr == pytest.approx(dy)

    
def check_delchi(title='Sigma Residuals for example'):
    """Check that the delchi plot has not changed"""

    rplot = ui._session._delchiplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y)
    assert rplot.y == pytest.approx([(y - 35) / dy
                                     for y,dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_y])

    
def check_chisqr():
    """Check that the chisqr plot has not changed"""

    rplot = ui._session._chisqrplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == '$\\chi^2$'
    assert rplot.title == '$\\chi^2$ for example'
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y)
    assert rplot.y == pytest.approx([((y - 35) / dy)**2
                                     for y,dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr is None


def check_fit():
    """Check that the fit plot has not changed"""

    check_example()
    check_model()


def check_fit_resid():
    """Check that the fit + resid plot has not changed"""

    check_example(xlabel='')
    check_model(xlabel='')
    check_resid(title='')


def check_fit_delchi():
    """Check that the fit + delchi plot has not changed"""

    check_example(xlabel='')
    check_model(xlabel='')
    check_delchi(title='')


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,changefunc,checkfunc",
                         [(ui.plot_data, change_example, check_example),
                          (ui.plot_model, change_model, check_model),
                          (ui.plot_source, change_model, check_source),
                          (ui.plot_resid, change_model, check_resid),
                          (ui.plot_ratio, change_example, check_ratio),
                          (ui.plot_delchi, change_example, check_delchi),
                          (ui.plot_chisqr, change_example, check_chisqr),
                          (ui.plot_fit, change_fit, check_fit),
                          (ui.plot_fit_resid, change_fit, check_fit_resid),
                          (ui.plot_fit_delchi, change_fit, check_fit_delchi),
                         ])
def test_plot_xxx_replot(idval, plotfunc, changefunc, checkfunc):
    """Can we plot, change data, plot with reploat and see a difference?

    This relies on accessing the undelying session object directly,
    as the user-supplied accessors either aren't suited (e.g.
    ui.get_data_plot always recreates the plot structure) or don't
    exist.

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier to use
    plotfunc
        The function to call to create the plot. If idval is None it
        is called with no argument, otherwise with idval.
    changefunc
        The function to call to change the setup (e.g. data or model).
        It is called with idval.
    checkfunc
        The function which performs the checks on the plot. It is called
        with no argument.
    """

    setup_example(idval)
    if idval is None:
        plotfunc()
    else:
        plotfunc(idval)

    changefunc(idval)
    
    # Recreate the plot
    #
    if idval is None:
        plotfunc(replot=True)
    else:
        plotfunc(idval, replot=True)

    checkfunc()
