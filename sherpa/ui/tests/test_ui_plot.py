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

# the chips plot tests mean that we can't test the plot instances
# from sherpa.plot import DataPlot, FitPlot, ModelPlot

from sherpa.stats import Chi2Gehrels
from sherpa.utils.testing import requires_plotting

import pytest


_data_x = [10, 20, 40, 90]
_data_y = [10, 40, 30, 50]
_data_y2 = [12, 45, 33, 49]


def example_data():
    """Create an example data set."""

    # Note: we copy the x and y arrays just so there's no accidental
    # aliasing.
    #
    x = [d for d in _data_x]
    y = [d for d in _data_y]
    return ui.Data1D('example', x, y)


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


def change_example(idval):
    """Change the example y values (created by setup_example)"""

    d = ui.get_data(idval)
    # copy the values to ensure _data_y2 isn't changed by accident
    d.y = [d for d in _data_y2]


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


def check_example_changed(xlabel='x'):
    """Check that the data plot has changed

    Assumes change_example has been called
    """

    dplot = ui._session._dataplot

    assert dplot.xlabel == xlabel
    assert dplot.ylabel == 'y'
    assert dplot.title == 'example'
    assert dplot.x == pytest.approx(_data_x)
    assert dplot.y == pytest.approx(_data_y2)
    assert dplot.xerr is None

    # Should use approximate equality here
    assert dplot.yerr == pytest.approx(calc_errors(_data_y2))


def check_model_plot(plot, title='Model', xlabel='x', modelval=35):
    """Helper for check_model/source"""

    assert plot.xlabel == xlabel
    assert plot.ylabel == 'y'
    assert plot.title == title
    assert plot.x == pytest.approx(_data_x)
    assert plot.y == pytest.approx([modelval for x in _data_x])
    assert plot.xerr is None
    assert plot.yerr is None


def check_model(xlabel='x'):
    """Check that the model plot has not changed"""

    check_model_plot(ui._session._modelplot,
                     title='Model', xlabel=xlabel)


def check_model_changed(xlabel='x'):
    """Check that the model plot has changed

    Assumes change_model has been called
    """

    check_model_plot(ui._session._modelplot,
                     title='Model', xlabel=xlabel,
                     modelval=41)


def check_source():
    """Check that the source plot has not changed"""

    check_model_plot(ui._session._sourceplot,
                     title='Source')


def check_source_changed():
    """Check that the source plot has changed

    Assumes change_model has been called
    """

    check_model_plot(ui._session._sourceplot,
                     title='Source', modelval=41)


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


def check_resid_changed(title='Residuals for example'):
    """Check that the resid plot has changed

    Assumes that change_model has been called
    """

    rplot = ui._session._residplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 41 for y in _data_y])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y))


def check_resid_changed2(title='Residuals for example'):
    """Check that the resid plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui._session._residplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 41 for y in _data_y2])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y2))


def check_ratio(title='Ratio of Data to Model for example'):
    """Check that the ratio plot has not changed"""

    rplot = ui._session._ratioplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 35 for y in _data_y])
    assert rplot.xerr is None
    dy = [dy / 35 for dy in calc_errors(_data_y)]
    assert rplot.yerr == pytest.approx(dy)


def check_ratio_changed():
    """Check that the ratio plot has changed

    Assumes that change_example has been called
    """

    rplot = ui._session._ratioplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == 'Ratio of Data to Model for example'
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 35 for y in _data_y2])
    assert rplot.xerr is None
    dy = [dy / 35 for dy in calc_errors(_data_y2)]
    assert rplot.yerr == pytest.approx(dy)


def check_ratio_changed2(title='Ratio of Data to Model for example'):
    """Check that the ratio plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui._session._ratioplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 41 for y in _data_y2])
    assert rplot.xerr is None
    dy = [dy / 41 for dy in calc_errors(_data_y2)]
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
                                     for y, dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_y])


def check_delchi_changed(title='Sigma Residuals for example'):
    """Check that the delchi plot has changed

    Assumes that change_example has been called
    """

    rplot = ui._session._delchiplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([(y - 35) / dy
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_x])


def check_delchi_changed2(title='Sigma Residuals for example'):
    """Check that the delchi plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui._session._delchiplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([(y - 41) / dy
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_x])


def check_chisqr():
    """Check that the chisqr plot has not changed"""

    rplot = ui._session._chisqrplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == '$\\chi^2$'
    assert rplot.title == '$\\chi^2$ for example'
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y)
    assert rplot.y == pytest.approx([((y - 35) / dy)**2
                                     for y, dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr is None


def check_chisqr_changed():
    """Check that the chisqr plot has changed

    Assumes that change_example has been called
    """

    rplot = ui._session._chisqrplot
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == '$\\chi^2$'
    assert rplot.title == '$\\chi^2$ for example'
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([((y - 35) / dy)**2
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr is None


def check_fit():
    """Check that the fit plot has not changed"""

    check_example()
    check_model()


def check_fit_changed():
    """Check that the fit plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed()
    check_model_changed()


def check_fit_resid():
    """Check that the fit + resid plot has not changed"""

    check_example(xlabel='')
    check_model(xlabel='')
    check_resid(title='')


def check_fit_resid_changed():
    """Check that the fit + resid plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(xlabel='')
    check_model_changed(xlabel='')
    check_resid_changed2(title='')


def check_fit_ratio():
    """Check that the fit + ratio plot has not changed"""

    check_example(xlabel='')
    check_model(xlabel='')
    check_ratio(title='')


def check_fit_ratio_changed():
    """Check that the fit + ratio plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(xlabel='')
    check_model_changed(xlabel='')
    check_ratio_changed2(title='')


def check_fit_delchi():
    """Check that the fit + delchi plot has not changed"""

    check_example(xlabel='')
    check_model(xlabel='')
    check_delchi(title='')


def check_fit_delchi_changed():
    """Check that the fit + delchi plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(xlabel='')
    check_model_changed(xlabel='')
    check_delchi_changed2(title='')


_plot_all = [
    {'plot': ui.plot_data, 'change': change_example,
     'check': check_example, 'check_changed': check_example_changed},
    {'plot': ui.plot_model, 'change': change_model,
     'check': check_model, 'check_changed': check_model_changed},
    {'plot': ui.plot_source, 'change': change_model,
     'check': check_source, 'check_changed': check_source_changed},
    {'plot': ui.plot_resid, 'change': change_model,
     'check': check_resid, 'check_changed': check_resid_changed},
    {'plot': ui.plot_ratio, 'change': change_example,
     'check': check_ratio, 'check_changed': check_ratio_changed},
    {'plot': ui.plot_delchi, 'change': change_example,
     'check': check_delchi, 'check_changed': check_delchi_changed},
    {'plot': ui.plot_chisqr, 'change': change_example,
     'check': check_chisqr, 'check_changed': check_chisqr_changed},
    {'plot': ui.plot_fit, 'change': change_fit,
     'check': check_fit, 'check_changed': check_fit_changed},
    {'plot': ui.plot_fit_resid, 'change': change_fit,
     'check': check_fit_resid, 'check_changed': check_fit_resid_changed},
    {'plot': ui.plot_fit_ratio, 'change': change_fit,
     'check': check_fit_ratio, 'check_changed': check_fit_ratio_changed},
    {'plot': ui.plot_fit_delchi, 'change': change_fit,
     'check': check_fit_delchi, 'check_changed': check_fit_delchi_changed}]

_plot_opts = [(p['plot'], p['check']) for p in _plot_all]
_plot_replot_opts = [(p['plot'], p['change'], p['check']) for p in _plot_all]
_plot_change_opts = [(p['plot'], p['change'], p['check_changed'])
                     for p in _plot_all]


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("pfunc, checkfunc", _plot_opts)
def test_plot_xxx(idval, pfunc, checkfunc):
    """Can we call a plot_xxx routine?

    There is limited testing that the plot call worked (this
    tests that the underlying data objects in the UI session
    were updated, not that the plot was created by the backend).

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier to use
    plotfunc
        The function to call to create the plot. If idval is None it
        is called with no argument, otherwise with idval.
    checkfunc
        The function which performs the checks on the plot. It is called
        with no argument.

    See Also
    --------
    test_plot_xxx_change, test_plot_xxx_replot

    """

    setup_example(idval)
    if idval is None:
        pfunc()
    else:
        pfunc(idval)

    checkfunc()


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,changefunc,checkfunc", _plot_replot_opts)
def test_plot_xxx_replot(idval, plotfunc, changefunc, checkfunc):
    """Can we plot, change data, plot with replot and see no difference?

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

    See Also
    --------
    test_plot_xxx, test_plot_xxx_change

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


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,changefunc,checkfunc", _plot_change_opts)
def test_plot_xxx_change(idval, plotfunc, changefunc, checkfunc):
    """Can we plot, change data, plot and see a difference?

    Unlike test_plot_xxx_replot, this does not set replot to True,
    so it should see the changed data in the plot structures.

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

    See Also
    --------
    test_plot_xxx, test_plot_xxx_replot

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
        plotfunc()
    else:
        plotfunc(idval)

    checkfunc()


_dplot = (ui.get_data_plot_prefs, "_dataplot", ui.plot_data)
_mplot = (ui.get_model_plot_prefs, "_modelplot", ui.plot_model)


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("getprefs,attr,plotfunc",
                         [_dplot, _mplot])
def test_prefs_change_session_objects(getprefs, attr, plotfunc):
    """Is a plot-preference change also reflected in the session object?

    This is intended to test an assumption that will be used in the
    plot_fit_xxx routines rather than of an explicit user-visible
    behavior. The test may be "obvious behavior" given how
    get_data_plot_prefs works, but DJB wanted to ensure this
    behavior/assumption was tested.
    """

    # This has to be retrieved here, rather than passed in in the
    # parametrize list, as the ui._session object is changed by
    # the clean_ui fixture.
    #
    session = getattr(ui._session, attr)

    # All but the last assert are just to check things are behaving
    # as expected (and stuck into one routine rather than have a
    # bunch of tests that repeat a subset of this test)
    #
    prefs = getprefs()
    assert not prefs['xlog']
    assert not session.plot_prefs['xlog']
    assert session.x is None

    prefs['xlog'] = True
    assert session.plot_prefs['xlog']

    setup_example(None)
    plotfunc()

    assert session.plot_prefs['xlog']
    assert session.x is not None

    prefs['xlog'] = False

    # The aim of the test is to check that the session plot object
    # has been updated with the new preference setting.
    #
    assert not session.plot_prefs['xlog']


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
def test_prefs_change_session_objects_fit():
    """Is plot-preference change reflected in the fitplot session object?

    This is test_prefs_change_session_objects but for the _fitplot
    attribute. This test encodes the current behavior - so we can see
    if things change in the future - rather than being a statement
    about what we expect/want to happen.
    """

    plotobj = ui._session._fitplot
    assert plotobj.dataplot is None
    assert plotobj.modelplot is None

    dprefs = ui.get_data_plot_prefs()
    mprefs = ui.get_model_plot_prefs()

    # Ensure we are actually changing a preference setting
    #
    assert not dprefs['xlog']
    assert not mprefs['ylog']

    dprefs['xlog'] = True
    mprefs['ylog'] = True

    setup_example(12)
    ui.plot_fit(12)

    # We have already checked this in previous tests, but
    # just in case
    #
    assert ui._session._dataplot.plot_prefs['xlog']
    assert ui._session._modelplot.plot_prefs['ylog']

    # Now check that the fit plot has picked up these changes;
    # the simplest way is to check that the data/model plots
    # are now referencing the underlying _data/_model plot
    # attributes. An alternative would be to check that
    # plotobj.dataplot.plot_prefs['xlog'] is True
    # which is less restrictive, but for not check the
    # equality
    #
    assert plotobj.dataplot is ui._session._dataplot
    assert plotobj.modelplot is ui._session._modelplot
