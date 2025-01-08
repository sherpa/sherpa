#
#  Copyright (C) 2019 - 2024
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

"""
Basic tests of the plot functionality in sherpa.ui. The idea is that
these - or most of them - can be run even when there is no plot backend.

There is a start to make these tests run on both Session classes,
but there's much to do.
"""

import logging

import numpy as np

import pytest

from sherpa import ui
from sherpa.ui.utils import Session as BaseSession
from sherpa.astro.ui.utils import Session as AstroSession

from sherpa.data import Data1D, Data1DInt, Data2D
from sherpa.models import basic
import sherpa.plot
from sherpa.stats import Chi2Gehrels
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, \
    IdentifierErr, PlotErr


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


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_fit_plot(idval, clean_ui):
    """Basic testing of get_fit_plot
    """

    setup_example(idval)
    if idval is None:
        f = ui.get_fit_plot()
    else:
        f = ui.get_fit_plot(idval)

    assert isinstance(f, sherpa.plot.FitPlot)
    assert isinstance(f.dataplot, sherpa.plot.DataPlot)
    assert isinstance(f.modelplot, sherpa.plot.ModelPlot)

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


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["data", "model"])
@pytest.mark.parametrize("arg", [None, 1, 'foo'])
def test_plot_prefs_xxx(session, ptype, arg):
    """Can we change and reset a preference.

    Pick the 'xlog' field, since the assumption is that: a) this
    defaults to 'False'; b) each plot type has this setting; c) we
    do not need to check all settings.

    """

    s = session()

    get_prefs = getattr(s, 'get_{}_plot_prefs'.format(ptype))

    prefs1 = get_prefs(arg)
    assert not prefs1['xlog']
    prefs1['xlog'] = True

    prefs2 = get_prefs(arg)
    assert prefs2['xlog']

    s.clean()
    prefs3 = get_prefs(arg)
    assert prefs1['xlog']
    assert prefs2['xlog']
    assert not prefs3['xlog']


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["data", "model"])
def test_plot_prefs_xxx_data1dint(session, ptype):
    """Data1DInt is different to Data1D.
    """

    s = session()

    get_prefs = getattr(s, 'get_{}_plot_prefs'.format(ptype))

    s.load_arrays(1, [1, 2, 4], [2, 3, 5], [4, 5, 10],
                  Data1DInt)

    s.load_arrays(2, [1, 2, 4], [4, 5, 10],
                  Data1D)

    # when there's no dataset it defaults to the plot, not histogram
    # prefs
    prefs = get_prefs('bob')
    assert 'xerrorbars' in prefs
    assert not prefs['xlog']
    prefs['xlog'] = True

    # It's not easy to check the difference between
    # point and histogram preferences.
    #
    # I also check xerrorbars as we want this for histograms.
    #
    prefs = get_prefs()
    assert 'xerrorbars' in prefs
    assert not prefs['xlog']

    prefs = get_prefs(2)
    assert 'xerrorbars' in prefs
    assert prefs['xlog']


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


def check_example(idval, xlabel='x'):
    """Check that the data plot has not changed"""

    dplot = ui.get_data_plot(id=idval, recalc=False)

    assert dplot.xlabel == xlabel
    assert dplot.ylabel == 'y'
    assert dplot.title == 'example'
    assert dplot.x == pytest.approx(_data_x)
    assert dplot.y == pytest.approx(_data_y)
    assert dplot.xerr is None

    # Should use approximate equality here
    assert dplot.yerr == pytest.approx(calc_errors(_data_y))


def check_example_changed(idval, xlabel='x'):
    """Check that the data plot has changed

    Assumes change_example has been called
    """

    dplot = ui.get_data_plot(id=idval, recalc=False)

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


def check_model(idval, xlabel='x'):
    """Check that the model plot has not changed"""

    mplot = ui.get_model_plot(id=idval, recalc=False)
    check_model_plot(mplot, title='Model', xlabel=xlabel)


def check_model_changed(idval, xlabel='x'):
    """Check that the model plot has changed

    Assumes change_model has been called
    """

    mplot = ui.get_model_plot(id=idval, recalc=False)
    check_model_plot(mplot, title='Model', xlabel=xlabel, modelval=41)


def check_source(idval):
    """Check that the source plot has not changed"""

    splot = ui.get_source_plot(id=idval, recalc=False)
    check_model_plot(splot, title='Source')


def check_source_changed(idval):
    """Check that the source plot has changed

    Assumes change_model has been called
    """

    splot = ui.get_source_plot(id=idval, recalc=False)
    check_model_plot(splot, title='Source', modelval=41)


def check_resid(idval, title='Residuals for example'):
    """Check that the resid plot has not changed"""

    rplot = ui.get_resid_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 35 for y in _data_y])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y))


def check_resid_changed(idval, title='Residuals for example'):
    """Check that the resid plot has changed

    Assumes that change_model has been called
    """

    rplot = ui.get_resid_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 41 for y in _data_y])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y))


def check_resid_changed2(idval, title='Residuals for example'):
    """Check that the resid plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui.get_resid_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data - Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y - 41 for y in _data_y2])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx(calc_errors(_data_y2))


def check_ratio(idval, title='Ratio of Data to Model for example'):
    """Check that the ratio plot has not changed"""

    rplot = ui.get_ratio_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 35 for y in _data_y])
    assert rplot.xerr is None
    dy = [dy / 35 for dy in calc_errors(_data_y)]
    assert rplot.yerr == pytest.approx(dy)


def check_ratio_changed(idval):
    """Check that the ratio plot has changed

    Assumes that change_example has been called
    """

    rplot = ui.get_ratio_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == 'Ratio of Data to Model for example'
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 35 for y in _data_y2])
    assert rplot.xerr is None
    dy = [dy / 35 for dy in calc_errors(_data_y2)]
    assert rplot.yerr == pytest.approx(dy)


def check_ratio_changed2(idval, title='Ratio of Data to Model for example'):
    """Check that the ratio plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui.get_ratio_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Data / Model'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)
    assert rplot.y == pytest.approx([y / 41 for y in _data_y2])
    assert rplot.xerr is None
    dy = [dy / 41 for dy in calc_errors(_data_y2)]
    assert rplot.yerr == pytest.approx(dy)


def check_delchi(idval, title='Sigma Residuals for example'):
    """Check that the delchi plot has not changed"""

    rplot = ui.get_delchi_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y)
    assert rplot.y == pytest.approx([(y - 35) / dy
                                     for y, dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_y])


def check_delchi_changed(idval, title='Sigma Residuals for example'):
    """Check that the delchi plot has changed

    Assumes that change_example has been called
    """

    rplot = ui.get_delchi_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([(y - 35) / dy
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_x])


def check_delchi_changed2(idval, title='Sigma Residuals for example'):
    """Check that the delchi plot has changed

    Assumes that change_example and change_model has been called
    """

    rplot = ui.get_delchi_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == 'Sigma'
    assert rplot.title == title
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([(y - 41) / dy
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr == pytest.approx([1.0 for y in _data_x])


def check_chisqr(idval):
    """Check that the chisqr plot has not changed"""

    rplot = ui.get_chisqr_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == '$\\chi^2$'
    assert rplot.title == '$\\chi^2$ for example'
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y)
    assert rplot.y == pytest.approx([((y - 35) / dy)**2
                                     for y, dy in zip(_data_y, dy)])
    assert rplot.xerr is None
    assert rplot.yerr is None


def check_chisqr_changed(idval):
    """Check that the chisqr plot has changed

    Assumes that change_example has been called
    """

    rplot = ui.get_chisqr_plot(id=idval, recalc=False)
    assert rplot.xlabel == 'x'
    assert rplot.ylabel == '$\\chi^2$'
    assert rplot.title == '$\\chi^2$ for example'
    assert rplot.x == pytest.approx(_data_x)

    dy = calc_errors(_data_y2)
    assert rplot.y == pytest.approx([((y - 35) / dy)**2
                                     for y, dy in zip(_data_y2, dy)])
    assert rplot.xerr is None
    assert rplot.yerr is None


def check_fit(idval):
    """Check that the fit plot has not changed"""

    check_example(idval)
    check_model(idval)


def check_fit_changed(idval):
    """Check that the fit plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(idval)
    check_model_changed(idval)


def check_fit_resid(idval):
    """Check that the fit + resid plot has not changed"""

    check_example(idval, xlabel='')
    check_model(idval, xlabel='')
    check_resid(idval, title='')


def check_fit_resid_changed(idval):
    """Check that the fit + resid plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(idval, xlabel='')
    check_model_changed(idval, xlabel='')
    check_resid_changed2(idval, title='')


def check_fit_ratio(idval):
    """Check that the fit + ratio plot has not changed"""

    check_example(idval, xlabel='')
    check_model(idval, xlabel='')
    check_ratio(idval, title='')


def check_fit_ratio_changed(idval):
    """Check that the fit + ratio plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(idval, xlabel='')
    check_model_changed(idval, xlabel='')
    check_ratio_changed2(idval, title='')


def check_fit_delchi(idval):
    """Check that the fit + delchi plot has not changed"""

    check_example(idval, xlabel='')
    check_model(idval, xlabel='')
    check_delchi(idval, title='')


def check_fit_delchi_changed(idval):
    """Check that the fit + delchi plot has changed

    Assumes that change_fit has been called
    """

    check_example_changed(idval, xlabel='')
    check_model_changed(idval, xlabel='')
    check_delchi_changed2(idval, title='')


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


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("pfunc, checkfunc", _plot_opts)
def test_plot_xxx(idval, pfunc, checkfunc, clean_ui):
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

    checkfunc(idval)


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,changefunc,checkfunc", _plot_replot_opts)
def test_plot_xxx_replot(idval, plotfunc, changefunc, checkfunc, clean_ui):
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

    checkfunc(idval)


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,changefunc,checkfunc", _plot_change_opts)
def test_plot_xxx_change(idval, plotfunc, changefunc, checkfunc, clean_ui):
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

    checkfunc(idval)


_dplot = (ui.get_data_plot_prefs, ui.get_data_plot, ui.plot_data)
_mplot = (ui.get_model_plot_prefs, ui.get_model_plot, ui.plot_model)


@pytest.mark.parametrize("getprefs,getplot,plotfunc",
                         [_dplot, _mplot])
def test_prefs_change_session_objects(getprefs, getplot, plotfunc, clean_ui):
    """Is a plot-preference change also reflected in the session object?

    This is intended to test an assumption that will be used in the
    plot_fit_xxx routines rather than of an explicit user-visible
    behavior. The test may be "obvious behavior" given how
    get_data_plot_prefs works, but DJB wanted to ensure this
    behavior/assumption was tested. The setup has changed since
    the test was written.
    """

    # This used to access the plot object directly, but is now
    # using an accessor function.
    #
    session = getplot(recalc=False)

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


def test_prefs_change_session_objects_fit(clean_ui):
    """Is plot-preference change reflected in the fitplot session object?

    This is test_prefs_change_session_objects but for the fit plot
    object. This test encodes the current behavior - so we can see
    if things change in the future - rather than being a statement
    about what we expect/want to happen.
    """

    plotobj = ui.get_fit_plot(recalc=False)
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
    assert ui.get_data_plot(id=12, recalc=False).plot_prefs['xlog']
    assert ui.get_model_plot(id=12, recalc=False).plot_prefs['ylog']

    # Now check that the fit plot has picked up these changes;
    # the simplest way is to check that the data/model plots
    # are now referencing the underlying _data/_model plot
    # attributes. An alternative would be to check that
    # plotobj.dataplot.plot_prefs['xlog'] is True
    # which is less restrictive, but for now check the
    # equality
    #
    assert plotobj.dataplot is ui.get_data_plot(id=12, recalc=False)
    assert plotobj.modelplot is ui.get_model_plot(id=12, recalc=False)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["data", "kernel", "model" ,"psf", "ratio", "resid", "source"])
def test_get_plot_prefs_returns_something(session, ptype):
    """Check this returns something.

    We do not require a plotting backend so the returned dictionary
    can be empty.

    """

    s = session()
    p = s.get_plot_prefs(ptype)
    assert isinstance(p, dict)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_plot_prefs_does_not_like_fit(session):
    """Check this errors out.

    We could support "fit" but the plot_prefs settings are currently
    unused so it is likely to be confusing to users, so we error out.

    """

    # This is a different error to test_get_plot_prefs_fails
    s = session()
    with pytest.raises(ArgumentErr,
                       match="^Use 'data' or 'model' instead of 'fit'$"):
        s.get_plot_prefs("fit")


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["", "not-a-plot", None, 1])
def test_get_plot_prefs_fails(session, ptype):
    """Check this call fails."""

    s = session()
    with pytest.raises(PlotErr,
                       match=r"^Plot type '.*' not found in \['data', "):
        s.get_plot_prefs(ptype)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_plot_prefs_recognizes_datatype(session):
    """Check that get_plot_prefs recognizes the data type."""

    s = session()
    s.dataspace1d(1, 5, id=1)
    s.dataspace1d(1, 5, id=2, dstype=Data1D)

    # The idea is that we add a key for id=1, which should be the
    # "histogram type" which will not be present for the "line type"
    # plot used for id=2. This is an attempt to be agnostic of the
    # plotting backend, including having no backend (when the default
    # preferences may be empty).
    #
    key = "NOT_AN_EXPECTED_VALID_KEY"
    p1 = s.get_plot_prefs("data", 1)
    assert key not in p1  # if this fails the key name needs to be changed
    p1[key] = 23

    # check it's remembered
    assert key in s.get_plot_prefs("data", 1)

    # check it isn't in p2
    p2 = s.get_plot_prefs("data", 2)
    assert key not in p2


@pytest.mark.parametrize("plotfunc", [ui.plot_cdf, ui.plot_pdf])
def test_plot_xdf(plotfunc):
    """Very basic check we can call plot_cdf/pdf

    This can be run even without a plotting backend available.
    """

    pvals = [0.1, 0.2, 0.1, 0.4, 0.3, 0.2, 0.1, 0.6]
    plotfunc(pvals)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_psf(session):
    """Very basic check we can call plot_psf/get_psf_plot

    This can be run even without a plotting backend available.
    """

    s = session()
    s._add_model_types(basic)

    x = np.arange(10, 40, 2)
    y = np.ones(x.size)
    s.load_arrays(1, x, y)

    psfmdl = s.create_model_component('gauss1d', 'psfmdl')
    psfmdl.fwhm = 5
    psfmdl.ampl = 2
    s.load_psf('psf', psfmdl)
    s.set_psf('psf')

    s.plot_psf()

    plotobj = s.get_psf_plot()
    assert isinstance(plotobj, sherpa.plot.PSFPlot)
    assert plotobj.title == 'gauss1d.psfmdl'
    assert plotobj.xlabel == 'x'
    assert plotobj.ylabel == 'y'
    assert plotobj.xerr is None
    assert plotobj.yerr is None
    assert plotobj.x == pytest.approx(x)

    yexp = psfmdl(x)
    assert plotobj.y == pytest.approx(yexp)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_psf_plot_recalc(session):
    """get_psf_plot with recalc=False

    This can be run even without a plotting backend available.
    """

    s = session()
    s._add_model_types(basic)

    x = np.arange(10, 40, 2)
    y = np.ones(x.size)
    s.load_arrays(1, x, y)

    psfmdl = s.create_model_component('gauss1d', 'psfmdl')
    psfmdl.fwhm = 5
    psfmdl.ampl = 2
    s.load_psf('psf', psfmdl)
    s.set_psf('psf')

    yexp = psfmdl(x)

    # ensure the plotobj is set
    s.get_psf_plot()

    # Change the PSF
    psfmdl.fwhm = 2
    psfmdl.ampl = 10

    yexp2 = psfmdl(x)

    # sanity check
    assert (yexp != yexp2).all()

    plotobj = s.get_psf_plot(recalc=False)
    assert plotobj.y == pytest.approx(yexp)

    plotobj = s.get_psf_plot()
    assert plotobj.y == pytest.approx(yexp2)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_kernel(session, caplog, plot_backends):
    """Very basic check we can call plot_kernel/get_kernel_plot

    This can be run even without a plotting backend available.
    """

    s = session()
    s._add_model_types(basic)

    x = np.arange(10, 40, 2)
    y = np.ones(x.size)
    s.load_arrays(1, x, y)

    psfmdl = s.create_model_component('gauss1d', 'psfmdl')
    psfmdl.fwhm = 5
    psfmdl.ampl = 2
    s.load_psf('psf', psfmdl)
    s.set_psf('psf')

    # TODO: check screen putput
    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.plot_kernel()

    assert len(caplog.records) == 1
    logname, loglvl, logmsg = caplog.record_tuples[0]
    assert logname == 'sherpa.instrument'
    assert loglvl == logging.INFO
    assert logmsg == 'PSF frac: 1.0'

    plotobj = s.get_kernel_plot()
    assert isinstance(plotobj, sherpa.plot.PSFKernelPlot)
    assert plotobj.title == 'PSF Kernel'
    assert plotobj.xlabel == 'x'
    assert plotobj.ylabel == 'y'
    assert plotobj.xerr is None
    assert plotobj.yerr is None
    assert plotobj.x == pytest.approx(x)

    yexp = psfmdl(x)
    yexp /= yexp.sum()
    assert plotobj.y == pytest.approx(yexp)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_kernel_plot_recalc(session):
    """get_kernel_plot with recalc=False

    This can be run even without a plotting backend available.
    """

    s = session()
    s._add_model_types(basic)

    x = np.arange(10, 40, 2)
    y = np.ones(x.size)
    s.load_arrays(1, x, y)

    psfmdl = s.create_model_component('gauss1d', 'psfmdl')
    psfmdl.fwhm = 5
    psfmdl.ampl = 2
    s.load_psf('psf', psfmdl)
    s.set_psf('psf')

    yexp = psfmdl(x)
    yexp /= yexp.sum()

    # ensure the plotobj is set
    s.get_kernel_plot()

    # Change the PSF
    psfmdl.fwhm = 2
    psfmdl.ampl = 10

    yexp2 = psfmdl(x)
    yexp2 /= yexp2.sum()

    # sanity check
    assert (yexp != yexp2).all()

    plotobj = s.get_kernel_plot(recalc=False)
    assert plotobj.y == pytest.approx(yexp)

    plotobj = s.get_kernel_plot()
    assert plotobj.y == pytest.approx(yexp2)


@pytest.mark.parametrize("getplot,plotfunc",
                         [(ui.get_resid_plot, ui.plot_resid),
                          (ui.get_delchi_plot, ui.plot_delchi)
                          ])
def test_plot_resid_ignores_ylog(getplot, plotfunc, clean_ui):
    """Do the plot_resid-family of routines ignore the ylog setting?

    Note that plot_chisqr is not included in support for ignoring
    ylog (since the data should be positive in this case).
    """

    prefs = getplot(recalc=False).plot_prefs

    setup_example(None)

    ui.set_ylog()
    assert prefs['ylog']

    plotfunc(ylog=True)

    # Note that the ylog setting has been removed (to reflect
    # what was displayed).
    #
    assert not prefs['ylog']


@pytest.mark.parametrize("getplot,plotfunc",
                         [(ui.get_resid_plot, ui.plot_fit_resid),
                          (ui.get_delchi_plot, ui.plot_fit_delchi)
                          ])
def test_plot_fit_resid_ignores_ylog(getplot, plotfunc, clean_ui):
    """Do the plot_resid-family of routines ignore the ylog setting?"""

    rprefs = getplot(recalc=False).plot_prefs
    dprefs = ui.get_data_plot(recalc=False).plot_prefs

    setup_example(None)

    ui.set_ylog()
    assert rprefs['ylog']
    assert dprefs['ylog']

    plotfunc(ylog=True)

    # Note that the ylog setting has been removed (to reflect
    # what was displayed), for the residual-style component only.
    #
    assert not rprefs['ylog']
    assert dprefs['ylog']


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("plottype", ["plot", "contour"])
def test_plot_contour_error_out_invalid(plottype, session):
    """plot()/contour() error out if argument name is invalid

    When it's not the first argument it will cause an error
    because of an invalid identifier.
    """

    s = session()
    func = getattr(s, plottype)
    with pytest.raises(PlotErr,
                       match="Plot type 'fooflan flim flam' not found in"):
        func("fooflan flim flam")


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_single(session, requires_pylab):
    """Can we call plot() with a single plot type?

    There's no real way to test this without a backend.
    """

    from matplotlib import pyplot as plt

    s = session()

    x = np.asarray([10, 20, 30, 45, 55, 70])
    y = np.asarray([5, 20, 15, 2, 17, 16])

    s.load_arrays(1, x, y)
    s.set_source(basic.Const1D())

    s.plot("data")

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_subplotspec().get_geometry() == (1, 1, 0, 0)
    assert ax.get_title() == ''
    assert ax.xaxis.get_label().get_text() == 'x'
    assert ax.yaxis.get_label().get_text() == 'y'

    plt.close()

    s.plot("model", 1)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_subplotspec().get_geometry() == (1, 1, 0, 0)
    assert ax.get_title() == 'Model'
    assert ax.xaxis.get_label().get_text() == 'x'
    assert ax.yaxis.get_label().get_text() == 'y'

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_multiple(session, requires_pylab):
    """Can we call plot() with multiple plot types?

    Also tests out sending in a kwarg.

    There's no real way to test this without a backend.
    """

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    xlo = np.asarray([10, 20, 30, 50, 55])
    xhi = np.asarray([20, 30, 45, 55, 70])
    y = np.asarray([20, 15, 2, 17, 16])

    # use Data1DInt as test_plot_single used Data1D
    #
    s.load_arrays('tst', xlo, xhi, y, Data1DInt)

    mdl = s.create_model_component('const1d', 'mdl')
    s.set_source('tst', mdl)

    # pick an odd number to plot
    s.plot("data", "tst", "model", "tst",
           "fit", "tst", "source", "tst",
           "ratio", "tst", alpha=0.8)

    # Note: I wanted to try source_component but it is not
    # clear it works when the id is not the default.

    fig = plt.gcf()
    assert len(fig.axes) == 5

    for idx, (ax, title, ylabel) in enumerate(zip(fig.axes,
                                                  ['', 'Model', '', 'Source',
                                                   'Ratio of Data to Model'],
                                                  ['y', 'y', 'y', 'y',
                                                   'Data / Model'])):

        assert ax.get_subplotspec().get_geometry() == (2, 3, idx, idx)
        assert ax.get_title() == title
        assert ax.xaxis.get_label().get_text() == 'x'
        assert ax.yaxis.get_label().get_text() == ylabel

        assert len(ax.lines) > 0
        assert ax.lines[0].get_alpha() == 0.8

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_contour_single(session, requires_pylab):
    """Can we call contour() with a single plot type?

    There's no real way to test this without a backend.
    """

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    mdl = s.create_model_component('gauss2d', 'mdl')
    mdl.xpos = 10
    mdl.ypos = 0
    mdl.fwhm = 3
    mdl.ampl = 100

    s.set_source(mdl)

    s.contour("data")

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_subplotspec().get_geometry() == (1, 1, 0, 0)
    assert ax.get_title() == ''
    assert ax.xaxis.get_label().get_text() == 'x0'
    assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()

    s.contour("model", 1)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_subplotspec().get_geometry() == (1, 1, 0, 0)
    assert ax.get_title() == 'Model'
    assert ax.xaxis.get_label().get_text() == 'x0'
    assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_contour_multiple(session, requires_pylab):
    """Can we call contour() with multiple plot types?

    There's no real way to test this without a backend.
    """

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    mdl = s.create_model_component('gauss2d', 'mdl')
    mdl.xpos = 10
    mdl.ypos = 0
    mdl.fwhm = 3
    mdl.ampl = 100

    s.set_source(mdl)

    s.contour("data", "model", "source", "fit", "ratio")

    fig = plt.gcf()
    assert len(fig.axes) == 5

    for idx, (ax, title) in enumerate(zip(fig.axes,
                                          ['', 'Model', 'Source', '',
                                           'Ratio of Data to Model'])):

        assert ax.get_subplotspec().get_geometry() == (2, 3, idx, idx)
        assert ax.get_title() == title
        assert ax.xaxis.get_label().get_text() == 'x0'
        assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("plotfunc,title,pcls",
                         [("data", "", sherpa.plot.DataContour),
                          ("model", "Model", sherpa.plot.ModelContour),
                          ("source", "Source", sherpa.plot.SourceContour),
                          ("resid", "Residuals", sherpa.plot.ResidContour),
                          ("ratio", "Ratio of Data to Model", sherpa.plot.RatioContour),
                          ("fit", "", sherpa.plot.FitContour),
                          ("fit_resid", None, None)])
def test_contour_xxx(plotfunc, title, pcls, session, requires_pylab):
    """Check we can call contour_xxx()/get_xxx_contour().

    There's no real way to test this without a backend.
    """

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    mdl = s.create_model_component('gauss2d', 'mdl')
    mdl.xpos = 10
    mdl.ypos = 0
    mdl.fwhm = 3
    mdl.ampl = 100

    s.set_source(mdl)

    getattr(s, "contour_" + plotfunc)()

    fig = plt.gcf()

    if plotfunc == 'fit_resid':
        assert len(fig.axes) == 2

        for idx, (ax, title) in enumerate(zip(fig.axes,
                                              ['', 'Residuals'])):

            assert ax.get_subplotspec().get_geometry() == (2, 1, idx, idx)
            assert ax.get_title() == title
            assert ax.xaxis.get_label().get_text() == 'x0'
            assert ax.yaxis.get_label().get_text() == 'x1'

    else:
        assert len(fig.axes) == 1

        ax = fig.axes[0]
        assert ax.get_subplotspec().get_geometry() == (1, 1, 0, 0)
        assert ax.get_title() == title
        assert ax.xaxis.get_label().get_text() == 'x0'
        assert ax.yaxis.get_label().get_text() == 'x1'

        plot = getattr(s, "get_{}_contour".format(plotfunc))()
        assert isinstance(plot, pcls)

        if plotfunc == 'fit':
            assert plot.datacontour.title == ''
            assert plot.modelcontour.title == 'Model'
        else:
            assert plot.title == title

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ctype", ["data", "model"])
def test_get_xxx_contour_prefs_pylab(ctype, session, requires_pylab):

    s = session()
    p = getattr(s, f"get_{ctype}_contour_prefs")()
    assert isinstance(p, dict)
    assert p == {'xlog': False, 'ylog': False,
                 'alpha': None, 'linewidths': None, 'colors': None,
                 'label': None, 'levels': None, 'linestyles': 'solid'}


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ctype", ["data", "kernel", "model" ,"psf", "ratio", "resid", "source"])
def test_get_contour_prefs_returns_something(session, ctype):
    """Check this returns something.

    We do not require a plotting backend so the returned dictionary
    can be empty.

    """

    s = session()
    p = s.get_contour_prefs(ctype)
    assert isinstance(p, dict)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_contour_prefs_does_not_like_fit(session):
    """Check this errors out.

    We could support "fit" but the contour_prefs settings are currently
    unused so it is likely to be confusing to users, so we error out.

    """

    # This is a different error to test_get_contour_prefs_fails
    s = session()
    with pytest.raises(ArgumentErr,
                       match="^Use 'data' or 'model' instead of 'fit'$"):
        s.get_contour_prefs("fit")


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ctype", ["", "not-a-plot", None, 1])
def test_get_contour_prefs_fails(session, ctype):
    """Check this call fails."""

    s = session()
    with pytest.raises(PlotErr,
                       match=r"^Plot type '.*' not found in \['data', "):
        s.get_contour_prefs(ctype)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_scatter_plot_empty(session):
    """Very basic check we can call get_scatter_plot

    This can be run even without a plotting backend available.
    """

    s = session()
    p = s.get_scatter_plot()
    assert isinstance(p, sherpa.plot.ScatterPlot)
    for f in ['x', 'y', 'xerr', 'yerr', 'xlabel', 'ylabel', 'title']:
        assert getattr(p, f) is None


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_scatter_plot(session):
    """Very basic check we can call plot_scatter/get_scatter_plot

    This can be run even without a plotting backend available.
    """

    s = session()

    x = [1, 5, 10]
    y = [-5, 2, -3]
    s.plot_scatter(x, y)

    p = s.get_scatter_plot()
    assert isinstance(p, sherpa.plot.ScatterPlot)

    assert p.x == pytest.approx(x)
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'x'
    assert p.ylabel == 'y'
    assert p.title == 'Scatter: (x,y)'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_scatter_plot_labels_noname(session):
    """Very basic check we can call plot_scatter/get_scatter_plot

    This can be run even without a plotting backend available.
    """

    s = session()

    x = [1, 5, 10]
    y = [-5, 2, -3]
    s.plot_scatter(x, y, xlabel='a x', ylabel='43')

    p = s.get_scatter_plot()
    assert isinstance(p, sherpa.plot.ScatterPlot)

    assert p.x == pytest.approx(x)
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'a x'
    assert p.ylabel == '43'
    assert p.title == 'Scatter: (x,y)'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_scatter_plot_labels(session):
    """Very basic check we can call plot_scatter/get_scatter_plot

    This can be run even without a plotting backend available.
    """

    s = session()

    x = [1, 5, 10]
    y = [-5, 2, -3]
    s.plot_scatter(x, y, xlabel='a x', ylabel='43', name='Fred Fred')

    p = s.get_scatter_plot()
    assert isinstance(p, sherpa.plot.ScatterPlot)

    assert p.x == pytest.approx(x)
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'a x'
    assert p.ylabel == '43'
    assert p.title == 'Scatter: Fred Fred'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_trace_plot_empty(session):
    """Very basic check we can call get_trace_plot

    This can be run even without a plotting backend available.
    """

    s = session()
    p = s.get_trace_plot()
    assert isinstance(p, sherpa.plot.TracePlot)
    for f in ['x', 'y', 'xerr', 'yerr', 'xlabel', 'ylabel', 'title']:
        assert getattr(p, f) is None


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_trace_plot(session):
    """Very basic check we can call get_trace_plot/plot_trace

    This can be run even without a plotting backend available.
    """

    s = session()

    y = [-5, 2, -3]
    s.plot_trace(y)

    p = s.get_trace_plot()
    assert isinstance(p, sherpa.plot.TracePlot)

    assert p.x == pytest.approx([0, 1, 2])
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'iteration'
    assert p.ylabel == 'x'
    assert p.title == 'Trace: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_trace_plot_labels_noname(session):
    """Very basic check we can call get_trace_plot/plot_trace

    Note that the xlabel setting doesn't seem to do anything.

    This can be run even without a plotting backend available.
    """

    s = session()

    y = [-5, 2, -3]
    s.plot_trace(y, xlabel='stats')

    p = s.get_trace_plot()
    assert isinstance(p, sherpa.plot.TracePlot)

    assert p.x == pytest.approx([0, 1, 2])
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'iteration'
    assert p.ylabel == 'x'
    assert p.title == 'Trace: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_trace_plot_labels(session):
    """Very basic check we can call get_trace_plot/plot_trace

    Note that the xlabel setting doesn't seem to do anything
    and the name is used interestingly.

    This can be run even without a plotting backend available.
    """

    s = session()

    y = [-5, 2, -3]
    s.plot_trace(y, xlabel='stats', name='Awesome sauce')

    p = s.get_trace_plot()
    assert isinstance(p, sherpa.plot.TracePlot)

    assert p.x == pytest.approx([0, 1, 2])
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'iteration'
    assert p.ylabel == 'Awesome sauce'
    assert p.title == 'Trace: Awesome sauce'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_cdf_plot_empty(session):

    s = session()
    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)
    for f in ["points", "x", "y", "median", "lower", "upper", "xlabel", "ylabel", "title"]:
        assert getattr(p, f) is None


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_replot_no_data(session, requires_pylab):
    """what does replot=True do for plot_cdf?

    The code doesn't check for this evantuality,
    so errors out (depends on the backend).
    """

    s = session()

    x = np.asarray([2, 8, 4, 6])

    # error can depend on matplotlib version
    with pytest.raises(ValueError):
        s.plot_cdf(x, replot=True)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_replot(session, requires_pylab):
    """what does replot=True do for plot_cdf?

    """

    from matplotlib import pyplot as plt

    s = session()

    z = np.asarray([1, 2, 3, 5])
    s.plot_cdf(z, xlabel='fo')

    s.plot_cdf(np.asarray([-10, -5, 0]), xlabel='ba', replot=True)

    fig = plt.gcf()

    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.xaxis.get_label().get_text() == 'fo'

    assert len(ax.lines) == 4
    line = ax.lines[0]

    assert line.get_xdata() == pytest.approx(z)
    assert line.get_ydata() == pytest.approx([0.25, 0.5, 0.75, 1])

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x)

    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)

    xsort = np.sort(x)
    assert p.points == pytest.approx(x)
    assert p.x == pytest.approx(xsort)
    assert p.y == pytest.approx(np.arange(1, 8) / 7)
    assert p.median == pytest.approx(15.0)
    assert p.lower == pytest.approx(11.903866)
    assert p.upper == pytest.approx(17.144201)
    assert p.xlabel == 'x'
    assert p.ylabel == 'p(<=x)'
    assert p.title == 'CDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_labels_noname(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x, xlabel='a b')

    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)

    assert p.xlabel == 'a b'
    assert p.ylabel == 'p(<=a b)'
    assert p.title == 'CDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_labels(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x, xlabel='a b', name='b a')

    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)

    assert p.xlabel == 'a b'
    assert p.ylabel == 'p(<=a b)'
    assert p.title == 'CDF: b a'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_cdf_plot_empty(session):

    s = session()

    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 10

    assert toks[0] == 'points = None'
    assert toks[1] == 'x      = None'
    assert toks[2] == 'y      = None'
    assert toks[3] == 'median = None'
    assert toks[4] == 'lower  = None'
    assert toks[5] == 'upper  = None'
    assert toks[6] == 'xlabel = None'
    assert toks[7] == 'ylabel = None'
    assert toks[8] == 'title  = None'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("use_numpy", [False, True])
def test_show_cdf_plot(session, use_numpy, old_numpy_printing):
    """This was to show issue #912 that has now been fixed.

    The display of the numeric values can depend on the
    NumPy version, so force the legacy output.
    """

    s = session()

    x = [20, 15, 25, 10]
    if use_numpy:
        x = np.asarray(x)

    s.plot_cdf(x)

    p = s.get_cdf_plot()
    assert isinstance(p, sherpa.plot.CDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 10

    assert toks[0] == 'points = [20,15,25,10]'
    assert toks[1] == 'x      = [10,15,20,25]'
    assert toks[2] == 'y      = [ 0.25, 0.5 , 0.75, 1.  ]'
    assert toks[3] == 'median = 17.5'
    assert toks[4].startswith('lower  = 12.37')
    assert toks[5].startswith('upper  = 22.62')
    assert toks[6] == 'xlabel = x'
    assert toks[7] == 'ylabel = p(<=x)'
    assert toks[8] == 'title  = CDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_pdf_plot_empty(session):

    s = session()
    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)
    for f in ["points", "xlo", "xhi", "y", "xlabel", "ylabel", "title"]:
        assert getattr(p, f) is None


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x)

    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)

    xgrid = np.arange(2, 8.5, 0.5)

    assert p.points == pytest.approx(x)
    assert p.xlo == pytest.approx(xgrid[:-1])
    assert p.xhi == pytest.approx(xgrid[1:])

    y = np.zeros(12)
    y[0] = 0.5
    y[4] = 0.5
    y[8] = 0.5
    y[11] = 0.5

    assert p.y == pytest.approx(y)

    assert p.xlabel == 'x'
    assert p.ylabel == 'probability density'
    assert p.title == 'PDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_replot_no_data(session, requires_pylab):
    """what does replot=True do for plot_pdf?

    The code doesn't check for this evantuality,
    so errors out (depends on the backend).
    """

    s = session()

    x = np.asarray([2, 8, 4, 6])

    # check on the error so we know when the code has changed
    with pytest.raises(TypeError) as exc:
        s.plot_pdf(x, replot=True)

    assert "'NoneType' has no len" in str(exc.value)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_replot(session, requires_pylab):
    """what does replot=True do for plot_pdf?

    """

    from matplotlib import pyplot as plt

    s = session()

    z = np.asarray([1, 2, 3, 5])
    s.plot_pdf(z, xlabel='fo', bins=5)

    s.plot_pdf(np.asarray([-10, -5, 0]), xlabel='ba', replot=True)

    fig = plt.gcf()

    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.xaxis.get_label().get_text() == 'fo'

    assert len(ax.lines) == 2
    line = ax.lines[0]

    x = np.asarray([1, 1.8, 1.8, 2.6, 2.6, 3.4, 3.4, 4.2, 4.2, 5])
    y = np.repeat([0.3125, 0.3125, 0.3125, 0, 0.3125], 2)

    assert line.get_xdata() == pytest.approx(x)
    assert line.get_ydata() == pytest.approx(y)

    pts = ax.lines[1]

    x = np.asarray([1.4, 2.2, 3, 3.8, 4.6])
    y = np.asarray([0.3125, 0.3125, 0.3125, 0, 0.3125])

    assert pts.get_xdata() == pytest.approx(x)
    assert pts.get_ydata() == pytest.approx(y)

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_labels_noname(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x, xlabel='x^2 x', bins=4)

    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)

    assert p.xlo.size == 4
    assert p.xhi.size == 4
    assert p.y.size == 4

    assert p.xlabel == 'x^2 x'
    assert p.ylabel == 'probability density'
    assert p.title == 'PDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_labels(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x, xlabel='x^2 x', name='no name', bins=4)

    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)

    assert p.xlo.size == 4
    assert p.xhi.size == 4
    assert p.y.size == 4

    assert p.xlabel == 'x^2 x'
    assert p.ylabel == 'probability density'
    assert p.title == 'PDF: no name'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_pdf_plot_empty(session):

    s = session()

    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 8

    assert toks[0] == 'points = None'
    assert toks[1] == 'xlo    = None'
    assert toks[2] == 'xhi    = None'
    assert toks[3] == 'y      = None'
    assert toks[4] == 'xlabel = None'
    assert toks[5] == 'ylabel = None'
    assert toks[6] == 'title  = None'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_pdf_plot(session, old_numpy_printing):
    """This is important as it also checks normed=False

    The display of the numeric values can depend on the
    NumPy version, so force the legacy output.
    """

    s = session()

    x = np.asarray([20, 15, 25, 10])
    s.plot_pdf(x, bins=3, normed=False)

    p = s.get_pdf_plot()
    assert isinstance(p, sherpa.plot.PDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 8

    assert toks[0] == 'points = [20,15,25,10]'
    assert toks[1] == 'xlo    = [ 10., 15., 20.]'
    assert toks[2] == 'xhi    = [ 15., 20., 25.]'
    assert toks[3] == 'y      = [1,1,2]'
    assert toks[4] == 'xlabel = x'
    assert toks[5] == 'ylabel = probability density'
    assert toks[6] == 'title  = PDF: x'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_data_plot_recalc(session):
    """Basic testing of get_data_plot(recalc=False)"""

    s = session()

    s.load_arrays(1, [1, 2], [1, 0])
    s.get_data_plot()

    s.load_arrays(1, [20, 30, 40], [25, 40, 60], [10, 12, 14], Data1DInt)

    p = s.get_data_plot(recalc=False)
    assert isinstance(p, sherpa.plot.DataHistogramPlot)
    assert p.xlo is None
    assert p.y is None

    p = s.get_data_plot(recalc=True)
    assert isinstance(p, sherpa.plot.DataHistogramPlot)
    assert p.xlo == pytest.approx([20, 30, 40])
    assert p.xhi == pytest.approx([25, 40, 60])
    assert p.y == pytest.approx([10, 12, 14])


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype,extraargs",
                         [('model', []), ('model_component', ['mdl']),
                          ('source', []), ('source_component', ['mdl'])])
def test_xxx_plot_nodata(ptype, extraargs, session):
    """Basic testing of get_xxx_plot when there's no data"""

    s = session()
    s._add_model_types(basic)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1
    s.set_source(mdl)

    func = getattr(s, 'get_{}_plot'.format(ptype))
    retval = func(*extraargs, recalc=False)
    assert retval.y is None


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype,extraargs",
                         [('model', []), ('model_component', ['mdl']),
                          ('source', []), ('source_component', ['mdl'])])
def test_xxx_plot_nodata_recalc(ptype, extraargs, session):
    """Basic testing of get_xxx_plot when there's no data and recalc=True"""

    s = session()
    s._add_model_types(basic)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1
    s.set_source(mdl)

    func = getattr(s, 'get_{}_plot'.format(ptype))
    with pytest.raises(IdentifierErr):
        func(*extraargs)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype,extraargs",
                         [('model', []), ('model_component', ['mdl']),
                          ('source', []), ('source_component', ['mdl'])])
def test_model_plot_recalc(ptype, extraargs, session):
    """Basic testing of get_model_plot(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    s.load_arrays(1, [1, 2], [1, 0])
    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1
    s.set_source(mdl)

    func = getattr(s, 'get_{}_plot'.format(ptype))

    # Seed the data for the recalc=False call
    func(*extraargs)

    s.load_arrays(1, [20, 30, 40], [25, 40, 60], [10, 12, 14], Data1DInt)
    s.set_source(mdl)

    # What data should be returned here? At the moment it uses the
    # current dataset to identify the value, but perhaps it should
    # just use the ModelPlot.
    #
    p = func(*extraargs, recalc=False)
    assert isinstance(p, sherpa.plot.ModelHistogramPlot)
    assert p.xlo is None
    assert p.xhi is None
    assert p.y is None

    p = func(*extraargs, recalc=True)
    assert isinstance(p, sherpa.plot.ModelHistogramPlot)
    assert p.xlo == pytest.approx([20, 30, 40])
    assert p.xhi == pytest.approx([25, 40, 60])
    assert p.y == pytest.approx([162.5, 450, 1200])


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype,pclass,y1,y2",
                         [('resid', sherpa.plot.ResidPlot, [-10, -12], [-20, -28, -36]),
                          ('ratio', sherpa.plot.RatioPlot, [1/11, 0], [1/3, 12/40, 14/50]),
                          ('delchi', sherpa.plot.DelchiPlot, [-5, -6], [-10, -14, -18]),
                          ('chisqr', sherpa.plot.ChisqrPlot, [25, 36], [100, 14*14, 18*18]),
                          ])
def test_xxx_plot_recalc(ptype, pclass, y1, y2, session):
    """Basic testing of get_xxx_plot(recalc=False)

    Unlike data/model this does not try changing the type
    of the dataset.
    """

    s = session()
    s._add_model_types(basic)

    s.load_arrays(1, [1, 2], [1, 0])
    s.set_staterror(1, [2, 2])
    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1
    s.set_source(mdl)

    # Set up the data to check in the recalc=False call below
    func = getattr(s, 'get_{}_plot'.format(ptype))
    func()

    s.load_arrays(1, [20, 30, 40], [10, 12, 14])
    s.set_staterror(1, [2, 2, 2])
    s.set_source(mdl)

    p = func(recalc=False)
    assert isinstance(p, pclass)

    assert p.x == pytest.approx([1, 2])
    assert p.y == pytest.approx(y1)

    p = func(recalc=True)
    assert isinstance(p, pclass)

    assert p.x == pytest.approx([20, 30, 40])
    assert p.y == pytest.approx(y2)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_fit_plot_recalc(session):
    """Basic testing of get_fit_plot(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    s.load_arrays(1, [1, 2], [1, 0])
    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1
    s.set_source(mdl)

    s.get_fit_plot()

    s.load_arrays(1, [20, 30, 40], [10, 12, 14])
    s.set_source(mdl)

    p = s.get_fit_plot(recalc=False)
    assert isinstance(p, sherpa.plot.FitPlot)

    assert p.dataplot.x == pytest.approx([1, 2])
    assert p.dataplot.y == pytest.approx([1, 0])

    assert p.modelplot.x == pytest.approx([1, 2])
    assert p.modelplot.y == pytest.approx([11, 12])

    p = s.get_fit_plot(recalc=True)
    assert isinstance(p, sherpa.plot.FitPlot)

    assert p.dataplot.x == pytest.approx([20, 30, 40])
    assert p.dataplot.y == pytest.approx([10, 12, 14])

    assert p.modelplot.x == pytest.approx([20, 30, 40])
    assert p.modelplot.y == pytest.approx([30, 40, 50])


def check_pvalue(caplog, plot):
    """Is the output as expected?"""

    assert len(caplog.records) == 1
    logname, loglvl, logmsg = caplog.record_tuples[0]
    assert logname == 'sherpa.ui.utils'
    assert loglvl == logging.INFO

    toks = logmsg.split('\n')
    assert len(toks) == 5
    assert toks[0] == 'Likelihood Ratio Test'
    assert toks[1].startswith('null statistic   =  -52.56')
    assert toks[2].startswith('alt statistic    =  -54.93')
    assert toks[3].startswith('likelihood ratio =  2.36')

    # The p-value is very sensitive so not really a useful check,
    # and can sometimes get a < 1/nruns answer (and nruns=100).
    #
    assert toks[4].startswith('p-value          =  0.') or \
        toks[4].startswith('p-value          <  0.01')

    assert isinstance(plot, sherpa.plot.LRHistogram)
    assert plot.xlabel == 'Likelihood Ratio'
    assert plot.ylabel == 'Frequency'
    assert plot.title == 'Likelihood Ratio Distribution'

    assert plot.lr == pytest.approx(2.3637744995453147)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_pvalue_plot(session, caplog):
    """Basic testing of get_pvalue_plot

    This is a made-up test, and I have no idea whether it's
    sensible.
    """

    s = session()
    s._add_model_types(basic)

    s.set_stat('cash')
    s.set_method('simplex')

    x = np.asarray([5, 7, 9, 11, 13, 20, 22])
    y = np.asarray([5, 4, 7, 9, 2, 6, 5])

    s.load_arrays(1, x, y)
    bgnd = s.create_model_component('const1d', 'bgnd')
    line = s.create_model_component('gauss1d', 'line')

    bgnd.c0 = 5
    line.pos.set(10.5, frozen=True)
    line.fwhm.set(0.2, frozen=True)
    line.ampl = 2.5

    s.set_source(bgnd)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        p = s.get_pvalue_plot(bgnd, bgnd+line, num=100, recalc=True)

    check_pvalue(caplog, p)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pvalue_requires_null_model(session):
    """We need a null_model argument"""

    s = session()
    with pytest.raises(TypeError) as te:
        # We add conv_model (to an invalid value) to avoid having
        # to set up a dataset
        s.plot_pvalue(None, None, conv_model=True)

    assert str(te.value) == 'null model cannot be None'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pvalue_requires_alt_model(session):
    """We need a alt_model argument"""

    s = session()
    with pytest.raises(TypeError) as te:
        # We add conv_model (to an invalid value) to avoid having
        # to set up a dataset. Similarly all we need is null_model
        # to be set, not that it's set properly
        s.plot_pvalue(False, None, conv_model=True)

    assert str(te.value) == 'alternative model cannot be None'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pvalue(session, caplog, plot_backends):
    """Basic testing of plot_pvalue

    This is a made-up test, and I have no idea whether it's
    sensible.
    """

    s = session()
    s._add_model_types(basic)

    s.set_stat('cash')
    s.set_method('simplex')

    x = np.asarray([5, 7, 9, 11, 13, 20, 22])
    y = np.asarray([5, 4, 7, 9, 2, 6, 5])

    s.load_arrays(1, x, y)
    bgnd = s.create_model_component('const1d', 'bgnd')
    line = s.create_model_component('gauss1d', 'line')

    bgnd.c0 = 5
    line.pos.set(10.5, frozen=True)
    line.fwhm.set(0.2, frozen=True)
    line.ampl = 2.5

    s.set_source(bgnd)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        s.plot_pvalue(bgnd, bgnd+line, num=100)

    p = s.get_pvalue_plot(recalc=False)
    check_pvalue(caplog, p)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_split_plot(session):
    """Check we can call get_split_plot.

    Limited checking of the result.
    """

    s = session()
    splot = s.get_split_plot()
    assert isinstance(splot, sherpa.plot.SplitPlot)
    assert splot.rows == 2
    assert splot.cols == 1


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_component_plot_invalid(session):
    """invalid model argument for get_model_component_plot"""

    # Fortunately we don't need to set up a dataset (this
    # could change)
    #
    s = session()
    with pytest.raises(ArgumentTypeErr) as exc:
        s.get_model_component_plot(id=1, model=3)

    emsg = "'model' must be a model object or model expression string"
    assert str(exc.value) == emsg


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_component_plot_string(session):
    """Check we can call get_model_component_plot with string model"""

    s = session()
    s._add_model_types(basic)

    x = np.asarray([3, 7, 12])
    y = np.asarray([4, 5, 8])
    s.load_arrays(1, x, y)

    s.create_model_component('const1d', 'gmdl')
    gmdl.c0 = 14

    plot = s.get_model_component_plot(id=1, model='gmdl')
    assert plot.x == pytest.approx(x)
    assert plot.y == pytest.approx([14, 14, 14])
    assert plot.title == 'Model component: const1d.gmdl'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_model_component_plot_model(session):
    """Check we can call get_model_component_plot with a model"""

    s = session()
    s._add_model_types(basic)

    x = np.asarray([3, 7, 12])
    y = np.asarray([4, 5, 8])
    s.load_arrays(1, x, y)

    s.create_model_component('const1d', 'gmdl')
    gmdl.c0 = 14

    plot = s.get_model_component_plot(id=1, model=gmdl)
    assert plot.x == pytest.approx(x)
    assert plot.y == pytest.approx([14, 14, 14])
    assert plot.title == 'Model component: const1d.gmdl'


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_scatter_empty_replot(session, requires_pylab):
    """plot_scatter with replot=False and no data

    Just check the current behavior
    """

    from matplotlib import pyplot as plt

    s = session()

    x = np.arange(3)
    y = x + 5
    s.plot_scatter(x, y, replot=True)

    fig = plt.gcf()

    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.xaxis.get_label().get_text() == ''
    assert ax.yaxis.get_label().get_text() == ''

    assert len(ax.lines) == 1
    line = ax.lines[0]

    assert line.get_xdata() == [None]
    assert line.get_ydata() == [None]

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_scatter(session, requires_pylab):
    """Simple test of plot_scatter"""

    from matplotlib import pyplot as plt

    s = session()

    x = np.arange(3)
    y = x + 5
    s.plot_scatter(x, y, marker='*', linestyle=':',
                   color='r', markerfacecolor='k')

    fig = plt.gcf()

    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.xaxis.get_label().get_text() == 'x'
    assert ax.yaxis.get_label().get_text() == 'y'

    assert len(ax.lines) == 1
    line = ax.lines[0]

    assert line.get_xdata() == pytest.approx(x)
    assert line.get_ydata() == pytest.approx(y)

    assert line.get_color() == 'r'
    assert line.get_markerfacecolor() == 'k'
    assert line.get_marker() == '*'
    assert line.get_linestyle() == ':'

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_trace_empty_replot(session, requires_pylab):
    """plot_trace with replot=False and no data

    Just check the current behavior
    """

    s = session()

    y = np.arange(100, 104)

    # error can depend on matplotlib version
    with pytest.raises(ValueError):
        s.plot_trace(y, replot=True)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_trace(session, requires_pylab):
    """Simple test of plot_trace"""

    from matplotlib import pyplot as plt

    s = session()

    y = np.asarray([100, 101, 99, 100])
    s.plot_trace(y, xlabel='stat', name='bob',
                 marker='*', linestyle=':',
                 color='r', markerfacecolor='k')

    fig = plt.gcf()

    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.xaxis.get_label().get_text() == 'iteration'
    assert ax.yaxis.get_label().get_text() == 'bob'
    assert ax.get_title() == 'Trace: bob'

    assert len(ax.lines) == 1
    line = ax.lines[0]

    assert line.get_xdata() == pytest.approx(np.arange(4))
    assert line.get_ydata() == pytest.approx(y)

    assert line.get_color() == 'r'
    assert line.get_markerfacecolor() == 'k'
    assert line.get_marker() == '*'
    assert line.get_linestyle() == ':'

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_data_contour_recalc(session):
    """Basic testing of get_data_contour(recalc=False)"""

    s = session()

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.get_data_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    p = s.get_data_contour(recalc=False)
    assert isinstance(p, sherpa.plot.DataContour)
    assert p.x0 == pytest.approx(x0)
    assert p.x1 == pytest.approx(x1)
    assert p.y == pytest.approx(y)

    p = s.get_data_contour(recalc=True)
    assert isinstance(p, sherpa.plot.DataContour)
    assert p.x0 == pytest.approx(nx0)
    assert p.x1 == pytest.approx(nx1)
    assert p.y == pytest.approx(ny)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_model_contour_recalc(session):
    """Basic testing of get_model_contour(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.create_model_component('gauss2d', 'gmdl')
    s.set_source('gmdl')

    s.get_model_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    s.create_model_component('const2d', 'cmdl')
    s.set_source('cmdl')

    p = s.get_model_contour(recalc=False)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(x0)
    assert p.x1 == pytest.approx(x1)
    # just check the model isn't flat
    assert p.y.min() < p.y.max()

    p = s.get_model_contour(recalc=True)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(nx0)
    assert p.x1 == pytest.approx(nx1)
    # just check the model is flat
    assert p.y.min() == p.y.max()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_source_contour_recalc(session):
    """Basic testing of get_source_contour(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.create_model_component('gauss2d', 'gmdl')
    s.set_source('gmdl')

    s.get_source_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    s.create_model_component('const2d', 'cmdl')
    s.set_source('cmdl')

    p = s.get_source_contour(recalc=False)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(x0)
    assert p.x1 == pytest.approx(x1)
    # just check the model isn't flat
    assert p.y.min() < p.y.max()

    p = s.get_source_contour(recalc=True)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(nx0)
    assert p.x1 == pytest.approx(nx1)
    # just check the model is flat
    assert p.y.min() == p.y.max()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_ratio_contour_recalc(session):
    """Basic testing of get_ratio_contour(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.create_model_component('gauss2d', 'gmdl')
    s.set_source('gmdl')

    s.get_ratio_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    s.create_model_component('const2d', 'cmdl')
    s.set_source('cmdl')

    p = s.get_ratio_contour(recalc=False)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(x0)
    assert p.x1 == pytest.approx(x1)
    # just check the model isn't flat
    assert p.y.min() < p.y.max()

    p = s.get_ratio_contour(recalc=True)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(nx0)
    assert p.x1 == pytest.approx(nx1)

    ygood = np.isfinite(p.y)
    assert ygood.sum() == ygood.size - 1
    assert not ygood[40]


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_resid_contour_recalc(session):
    """Basic testing of get_resid_contour(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.create_model_component('gauss2d', 'gmdl')
    s.set_source('gmdl')

    s.get_resid_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    s.create_model_component('const2d', 'cmdl')
    s.set_source('cmdl')

    p = s.get_resid_contour(recalc=False)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(x0)
    assert p.x1 == pytest.approx(x1)
    # just check the model isn't flat
    assert p.y.min() < p.y.max()

    p = s.get_resid_contour(recalc=True)
    assert isinstance(p, sherpa.plot.ModelContour)
    assert p.x0 == pytest.approx(nx0)
    assert p.x1 == pytest.approx(nx1)

    ygood = np.isfinite(p.y)
    assert ygood.sum() == ygood.size - 1
    assert not ygood[40]


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_fit_contour_recalc(session):
    """Basic testing of get_fit_contour(recalc=False)"""

    s = session()
    s._add_model_types(basic)

    x1, x0 = np.mgrid[-4:5, 6:15]

    x0 = x0.flatten()
    x1 = x1.flatten()
    y = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, x0, x1, y, Data2D)

    s.create_model_component('gauss2d', 'gmdl')
    s.set_source('gmdl')

    s.get_fit_contour()

    nx1, nx0 = np.mgrid[2:5, 12:15]

    nx0 = x0.flatten()
    nx1 = x1.flatten()
    ny = 100 / np.sqrt((x0 - 10)**2 + x1**2)

    s.load_arrays(1, nx0, nx1, ny, Data2D)

    s.create_model_component('const2d', 'cmdl')
    s.set_source('cmdl')

    p = s.get_fit_contour(recalc=False)
    assert isinstance(p, sherpa.plot.FitContour)
    pd = p.datacontour
    pm = p.modelcontour
    assert isinstance(pd, sherpa.plot.DataContour)
    assert isinstance(pm, sherpa.plot.ModelContour)
    assert pd.x0 == pytest.approx(x0)
    assert pd.x1 == pytest.approx(x1)
    assert pm.x0 == pytest.approx(x0)
    assert pm.x1 == pytest.approx(x1)
    assert pd.y == pytest.approx(y)
    # just check the model isn't flat
    assert pm.y.min() < pm.y.max()

    p = s.get_fit_contour(recalc=True)
    assert isinstance(p, sherpa.plot.FitContour)
    pd = p.datacontour
    pm = p.modelcontour
    assert isinstance(pd, sherpa.plot.DataContour)
    assert isinstance(pm, sherpa.plot.ModelContour)
    assert pd.x0 == pytest.approx(nx0)
    assert pd.x1 == pytest.approx(nx1)
    assert pm.x0 == pytest.approx(nx0)
    assert pm.x1 == pytest.approx(nx1)
    assert pd.y == pytest.approx(ny)
    # just check the model is flat
    assert pm.y.min() == pm.y.max()


@pytest.mark.parametrize("ptype",
                         ["resid", "ratio", "delchi"])
def test_plot_fit_xxx_pylab(ptype, clean_ui, requires_pylab):
    """Just ensure we can create a plot_fit_xxx call."""

    from matplotlib import pyplot as plt

    setup_example(1)
    pfunc = getattr(ui, 'plot_fit_{}'.format(ptype))
    pfunc(xlog=True, ylog=True)

    fig = plt.gcf()
    axes = fig.axes

    # This test occasionally fails because len(axes) == 3
    # but it's not obvious why - so let's print some
    # info in the hope it's informative
    print(plt.get_current_fig_manager())
    print(fig)
    print(axes)
    for ax in axes:
        print(ax.get_xlabel())
        print(ax.get_ylabel())
        print(ax.get_title())
        print(ax.lines)
        print(ax.get_xlim())
        print(ax.get_ylim())
        print('---')

    assert len(axes) == 2
    assert axes[0].xaxis.get_label().get_text() == ''

    assert axes[0].xaxis.get_scale() == 'log'
    assert axes[0].yaxis.get_scale() == 'log'

    assert axes[1].xaxis.get_scale() == 'log'
    assert axes[1].yaxis.get_scale() == 'linear'

    # Check we have the correct data (at least in the
    # number of data objects). The residual plot has
    # the data but also the axis line.
    #
    assert len(axes[0].lines) == 2
    assert len(axes[1].lines) == 2


@pytest.mark.parametrize("ptype",
                         ["resid", "ratio", "delchi"])
def test_plot_fit_xxx_overplot_pylab(ptype, caplog, clean_ui, requires_pylab):
    """Just ensure we can create a plot_fit_xxx(overplot=True) call."""

    from matplotlib import pyplot as plt

    setup_example(1)
    setup_example(2)

    pfunc = getattr(ui, 'plot_fit_{}'.format(ptype))
    pfunc(xlog=True, ylog=True)
    pfunc(2, overplot=True)

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 2
    assert axes[0].xaxis.get_label().get_text() == ''

    assert axes[0].xaxis.get_scale() == 'log'
    assert axes[0].yaxis.get_scale() == 'log'

    assert axes[1].xaxis.get_scale() == 'log'
    assert axes[1].yaxis.get_scale() == 'linear'

    # Check we have the correct data (at least in the
    # number of data objects). The residual plot has
    # the data but also the axis line.
    #
    assert len(axes[0].lines) == 4
    assert len(axes[1].lines) == 4

    # data is repeated so can check
    for idx in [0, 1]:
        l0 = axes[idx].lines
        assert l0[0].get_xydata() == pytest.approx(l0[2].get_xydata())
        assert l0[1].get_xydata() == pytest.approx(l0[3].get_xydata())


@pytest.mark.parametrize("idval", [None, "bob"])
def test_plot_fit_resid_handles_data_log(idval, clean_ui, requires_pylab):
    """Check that log handling is correct: data=log

    See also test_plot_fit_resid_handles_resid_log.

    I thought we had tests of this, but apparently not.
    """

    from matplotlib import pyplot as plt

    setup_example(idval)
    ui.set_xlog('data')
    ui.plot_fit_resid(idval)

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 2
    assert axes[0].xaxis.get_label().get_text() == ''

    assert axes[0].xaxis.get_scale() == 'log'
    assert axes[0].yaxis.get_scale() == 'linear'

    assert axes[1].xaxis.get_scale() == 'log'
    assert axes[1].yaxis.get_scale() == 'linear'


@pytest.mark.parametrize("idval", [None, "bob"])
def test_plot_fit_resid_handles_resid_log(idval, clean_ui, requires_pylab):
    """Check that log handling is correct: resid=log

    We need to decide whether we want the residual setting to override
    the linear display of the fit plot here. At present the code is
    that if resid has xlog set then both will be drawn logged (since
    the X axis is shared via a sharex=True argument to plt.subplots)
    but we may decide this should change.

    """

    from matplotlib import pyplot as plt

    setup_example(idval)
    ui.set_xlog('resid')
    ui.plot_fit_resid(idval)

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 2
    assert axes[0].xaxis.get_label().get_text() == ''

    assert axes[0].xaxis.get_scale() == 'log'
    assert axes[0].yaxis.get_scale() == 'linear'

    assert axes[1].xaxis.get_scale() == 'log'
    assert axes[1].yaxis.get_scale() == 'linear'


@pytest.mark.xfail
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("label", ["source", "model"])
@pytest.mark.parametrize("idval", [None, 1, "bob"])
def test_get_foo_component_plot_recalc(session, label, idval):
    """get_foo_component_plot(recalc=False) should not need a model.

    At the moment it does, hence the xfail mark. See issue #1512
    """

    s = session()
    s._add_model_types(basic)

    plot = getattr(s, f"plot_{label}_component")
    get = getattr(s, f"get_{label}_component_plot")

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 10
    mdl.c1 = 1

    if idval is None:
        s.load_arrays(1, [-1, 1, 2], [5, 12, 3])
        s.set_source(mdl)
    else:
        s.load_arrays(idval, [-1, 1, 2], [5, 12, 3])
        s.set_source(idval, mdl)

    if idval is None:
        plot(mdl)
        plotobj = get(recalc=False)
    else:
        plot(idval, mdl)
        plotobj = get(idval, recalc=False)

    # The main check is that the get call does not fail, but add some
    # basic check (fortunately as there is no response here the source
    # and model values are the same).
    #
    assert plotobj.y == pytest.approx([9, 11, 12])


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype,mclass,label",
                         [("source", sherpa.plot.ComponentSourcePlot,
                           "Source model"),
                          ("model", sherpa.plot.ComponentModelPlot,
                           "Model")
                          ])
def test_get_xxx_components_simple(session, ptype, mclass, label):
    """Simple check of get_xxx_components

    This does not require a backend
    """

    s = session()
    s._add_model_types(basic)

    s.load_arrays(3, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")

    c1.c0 = 0.5
    c2.c0 = 2
    c3.ampl = 4
    c3.xhi = 100
    c3.xlow = 7
    s.set_source(3, c1 * (c2 + c3))

    pfunc = getattr(s, f"get_{ptype}_components_plot")
    out = pfunc(3)

    assert isinstance(out, sherpa.plot.MultiPlot)
    assert len(out.plots) == 2
    assert isinstance(out.plots[0], mclass)
    assert isinstance(out.plots[1], mclass)

    assert out.plots[0].x == pytest.approx([1, 10, 50])
    assert out.plots[1].x == pytest.approx([1, 10, 50])

    # This should be c1 * c2 and c1 * c3
    assert out.plots[0].y == pytest.approx([1, 1, 1])
    assert out.plots[1].y == pytest.approx([0, 2, 2])

    assert out.plots[0].title == f"{label} component: scale1d.c1 * const1d.c2"
    assert out.plots[1].title == f"{label} component: scale1d.c1 * box1d.c3"


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_simple_mpl(session, ptype, requires_pylab):
    """Simple check of plot_xxx_components using matpotlib"""

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    s.set_default_id(3)
    s.load_arrays(3, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")

    c1.c0 = 0.5
    c2.c0 = 2
    c3.ampl = 4
    c3.xhi = 100
    c3.xlow = 7
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc()

    fig = plt.gcf()
    assert len(fig.axes) == 1

    axes = fig.axes[0]
    assert axes.get_xlabel() == 'x'
    assert axes.get_ylabel() == 'y'
    assert axes.get_title() == 'Component plot'

    assert len(axes.lines) == 2

    assert axes.lines[0].get_xdata() == pytest.approx([1, 10, 50])
    assert axes.lines[1].get_xdata() == pytest.approx([1, 10, 50])

    # This should be c1 * c2 and c1 * c3
    assert axes.lines[0].get_ydata() == pytest.approx([1, 1, 1])
    assert axes.lines[1].get_ydata() == pytest.approx([0, 2, 2])


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_simple_bokeh(session, ptype):
    """Simple check of plot_xxx_components using bokeh"""

    backend_mod = pytest.importorskip("sherpa.plot.bokeh_backend")
    backend = backend_mod.BokehBackend()

    s = session()
    s._add_model_types(basic)

    s.set_plot_backend(backend)

    s.set_default_id(3)
    s.load_arrays(3, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")

    c1.c0 = 0.5
    c2.c0 = 2
    c3.ampl = 4
    c3.xhi = 100
    c3.xlow = 7
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc()

    # What to check? At the moment *very* basic (less than
    # test_plot_xxx_components_simple_mpl).
    #
    fig = backend.current_fig
    assert len(fig.children) == 1

    fig0 = fig.children[0][0]
    assert len(fig0.axis) == 2
    assert fig0.xaxis.axis_label == "x"
    assert fig0.yaxis.axis_label == "y"
    assert fig0.title.text == "Component plot"

    # Unlike the matplotlib case we do not check the plotted data.


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_scalar_mpl(session, ptype, requires_pylab):
    """Can I set a kwarg to a scalar using matpotlib"""

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    s.set_default_id("x")
    s.load_arrays("x", [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")

    c1.c0 = 0.5
    c2.c0 = 2
    c3.ampl = 4
    c3.xhi = 100
    c3.xlow = 7
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc(color="black", alpha=0.5, linestyle="dashed")

    axes = plt.gca()
    assert len(axes.lines) == 2
    for l in axes.lines:
        assert l.get_color() == "black"
        assert l.get_alpha() == 0.5
        assert l.get_linestyle() == "--"


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_scalar_bokeh(session, ptype):
    """Can I set a kwarg to a scalar using bokeh"""

    backend_mod = pytest.importorskip("sherpa.plot.bokeh_backend")
    backend = backend_mod.BokehBackend()

    s = session()
    s._add_model_types(basic)

    s.set_plot_backend(backend)

    s.set_default_id("x")
    s.load_arrays("x", [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")

    c1.c0 = 0.5
    c2.c0 = 2
    c3.ampl = 4
    c3.xhi = 100
    c3.xlow = 7
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc(color="black", alpha=0.5, linestyle="dashed")

    # For now there is no check the kwargs were used


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_kwargs_mpl(session, ptype, requires_pylab):
    """Can we have per-plot kwargs? matplotlib"""

    from matplotlib import pyplot as plt

    s = session()
    s._add_model_types(basic)

    s.load_arrays(1, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc(color=["orange", "black"],
          alpha=0.5,
          linestyle=["dashed", "dotted"]
          )

    axes = plt.gca()
    assert len(axes.lines) == 2
    assert axes.lines[0].get_color() == "orange"
    assert axes.lines[1].get_color() == "black"

    assert axes.lines[0].get_alpha() == 0.5
    assert axes.lines[1].get_alpha() == 0.5

    assert axes.lines[0].get_linestyle() == "--"
    assert axes.lines[1].get_linestyle() == ":"


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
def test_plot_xxx_components_kwargs_bokeh(session, ptype):
    """Can we have per-plot kwargs? bokeh"""

    backend_mod = pytest.importorskip("sherpa.plot.bokeh_backend")
    backend = backend_mod.BokehBackend()

    s = session()
    s._add_model_types(basic)

    s.set_plot_backend(backend)

    s.load_arrays(1, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    pfunc(color=["orange", "black"],
          alpha=0.5,
          linestyle=["dashed", "dotted"]
          )

    # For now there is no check the kwargs were used


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["source", "model"])
@pytest.mark.parametrize("kwargs,badarg",
                         [({"color": ["orange", "black", "green"]}, "color"),
                          ({"color": ["black", "black"],
                            "alpha": [0.8]}, "alpha"),
                          ({"color": ["black", "black", "black"],
                            "alpha": [0.8, 0.7, 1.0]}, "color")
                          ])
def test_plot_xxx_components_kwargs_mismatch(kwargs, badarg, session, ptype, plot_backends):
    """Incorrect number of kwargs - what happens?

    It is okay to use plot_backends here because this is not
    the global Session but a temporary one.
    """

    s = session()
    s._add_model_types(basic)

    s.load_arrays(1, [1, 10, 50], [2, 3, 4])

    c1 = s.create_model_component("scale1d", "c1")
    c2 = s.create_model_component("const1d", "c2")
    c3 = s.create_model_component("box1d", "c3")
    s.set_source(c1 * (c2 + c3))

    pfunc = getattr(s, f"plot_{ptype}_components")
    msg = f"keyword '{badarg}': expected 2 elements but found "
    with pytest.raises(ValueError, match=f"^{msg}"):
        pfunc(**kwargs)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("name,val",
                         [("rows", 0), ("cols", 0),
                          ("rows", 1.4), ("cols", 1.9)])
def test_plot_invalid_size(session, name, val):
    """This errors out before accessing the data."""

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])

    kwargs = {}
    kwargs[name] = val
    with pytest.raises(ArgumentErr,
                       match=f"^{name} must be a positive integer, not {val}$"):
        s.plot("data", "data", **kwargs)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("nplots", [1, 2, 3, 4, 5])
def test_plot_overplot_too_many_plots(session, nplots):
    """Check we error out."""

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])

    # It's important to include non-square numbers in the test
    args1 = ["data"] * nplots
    args2 = ["data"] + args1

    s.plot(*args1)
    with pytest.raises(ArgumentErr,
                       match="^Too many elements for overplot$"):
        s.plot(*args2, overplot=True)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("nrows,ncols,nplots",
                         [(1, 1, 1),
                          (2, 1, 2),
                          (2, 2, 3),
                          (2, 2, 4),
                          (2, 3, 5),
                          (2, 3, 6),
                          (3, 3, 7),
                          (3, 3, 8),
                          (3, 3, 9),
                          (3, 4, 10),
                          (3, 4, 11),
                          (3, 4, 12)
                          ])
def test_plot_check_default_size_pylab(session, nrows, ncols, nplots, requires_pylab):
    """Check we get the expected size.

    We also check that xlog and alpha can be given as arrays.
    """

    from matplotlib import pyplot as plt

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])
    args = ["data"] * nplots
    s.plot(*args)

    fig = plt.gcf()
    assert len(fig.axes) == nplots
    for idx, ax in enumerate(fig.axes):
        assert ax.get_subplotspec().get_geometry() == (nrows, ncols, idx, idx)
        assert len(ax.lines) == 1

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("nrows,ncols,kwargs",
                         [(2, 3, {}),  # this is the default setting
                          (2, 3, {"rows": 1, "cols": 3}),  # also when requested nplots is too small
                          (2, 3, {"cols": 3}),
                          (2, 3, {"rows": 2}),
                          (2, 3, {"rows": 2, "cols": 3}),
                          (3, 2, {"cols": 2}),
                          (3, 2, {"rows": 3}),
                          (3, 2, {"rows": 3, "cols": 2}),
                          (1, 5, {"rows": 1}),
                          (1, 5, {"cols": 5}),
                          (5, 1, {"cols": 1}),
                          (5, 1, {"rows": 5})
                          ])
def test_plot_check_size_nrows_ncols_pylab(session, nrows, ncols, kwargs, requires_pylab):
    """Check we get the expected size.

    We also check that xlog and alpha can be given as arrays.
    """

    from matplotlib import pyplot as plt

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])
    args = ["data"] * 5
    kwargs["xlog"] = [False, True, False, True, False]
    kwargs["alpha"] = [0.2, 0.4, 0.6, 0.8, 1.0]
    kwargs["color"] = "black"
    s.plot(*args, **kwargs)

    fig = plt.gcf()
    assert len(fig.axes) == 5
    for idx, ax in enumerate(fig.axes):
        assert ax.get_subplotspec().get_geometry() == (nrows, ncols, idx, idx)
        assert len(ax.lines) == 1
        assert ax.lines[0].get_alpha() == pytest.approx((idx + 1) * 0.2)
        assert ax.lines[0].get_color() == "black"
        scale = "log" if idx % 2 else "linear"
        assert ax.xaxis.get_scale() == scale

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("rows,cols", [(None, 2), (1, 2),
                                       (2, None), (2, 1),
                                       (2, 2)])
def test_plot_singleton_ignores_sizes(session, rows, cols, requires_pylab):
    """Single plots used to be special cased to ignore size arguments.

    This is no-longer the case, so test the new behavior.
    """

    nrows = 1 if rows is None else rows
    ncols = 1 if cols is None else cols

    from matplotlib import pyplot as plt

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])
    s.plot("data", rows=rows, cols=cols)

    fig = plt.gcf()
    assert len(fig.axes) == 1
    assert fig.axes[0].get_subplotspec().get_geometry() == (nrows, ncols, 0, 0)

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_uses_split_plot_prefs(session, requires_pylab):
    """Check we can use get_split_plot to change the size."""

    from matplotlib import pyplot as plt

    s = session()
    s.load_arrays(1, [1, 2], [1, 2])

    sp = s.get_split_plot()
    assert sp.rows == 2
    assert sp.cols == 1
    sp.rows = 3
    sp.cols = 1
    s.plot("data", "data")

    fig = plt.gcf()
    assert len(fig.axes) == 2
    assert fig.axes[0].get_subplotspec().get_geometry() == (3, 1, 0, 0)
    assert fig.axes[1].get_subplotspec().get_geometry() == (3, 1, 1, 1)

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_overplot(session, requires_pylab):
    """What happens when overplot is set. See issue #2128.

    This is a regression test, to check the current behaviour.
    """

    from matplotlib import pyplot as plt

    x1 = [1, 2, 3]
    x2 = [21, 22, 23]
    y1 = [4, 5, 6]
    y2 = [-4, -5, -6]

    s = session()
    s.load_arrays(1, x1, y1)
    s.load_arrays(2, x2, y2)

    s.plot("data", "data", 2)

    # Validate the original plot
    fig = plt.gcf()
    assert len(fig.axes) == 2
    for ax in fig.axes:
        assert ax.get_title() == ""
        assert ax.xaxis.get_label().get_text() == "x"
        assert ax.yaxis.get_label().get_text() == "y"

        assert len(ax.lines) == 1

    l1 = fig.axes[0].lines[0]
    l2 = fig.axes[1].lines[0]
    assert l1.get_xdata() == pytest.approx(x1)
    assert l2.get_xdata() == pytest.approx(x2)
    assert l1.get_ydata() == pytest.approx(y1)
    assert l2.get_ydata() == pytest.approx(y2)

    # Now overplot the "other" data.
    #
    s.plot("data", 2, "data", overplot=True)

    fig = plt.gcf()
    assert len(fig.axes) == 2
    for ax in fig.axes:
        assert ax.get_title() == ""
        assert ax.xaxis.get_label().get_text() == "x"
        assert ax.yaxis.get_label().get_text() == "y"

        assert len(ax.lines) == 2

    l1 = fig.axes[0].lines[0]
    l2 = fig.axes[1].lines[0]
    assert l1.get_xdata() == pytest.approx(x1)
    assert l2.get_xdata() == pytest.approx(x2)
    assert l1.get_ydata() == pytest.approx(y1)
    assert l2.get_ydata() == pytest.approx(y2)

    l1 = fig.axes[0].lines[1]
    l2 = fig.axes[1].lines[1]
    assert l1.get_xdata() == pytest.approx(x2)
    assert l2.get_xdata() == pytest.approx(x1)
    assert l1.get_ydata() == pytest.approx(y2)
    assert l2.get_ydata() == pytest.approx(y1)

    plt.close()


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_overplot_smaller(session, requires_pylab):
    """Check we can overplot less plots than the original.

    This is a regression test, to check the current behaviour.
    """

    from matplotlib import pyplot as plt

    x1 = [1, 2, 3]
    x2 = [21, 22, 23]
    y1 = [4, 5, 6]
    y2 = [-4, -5, -6]

    s = session()
    s.load_arrays(1, x1, y1)
    s.load_arrays(2, x2, y2)

    s.plot("data", "data", 2)
    s.plot("data", 2, overplot=True, alpha=0.5)

    fig = plt.gcf()
    assert len(fig.axes) == 2
    assert len(fig.axes[0].lines) == 2
    assert len(fig.axes[1].lines) == 1

    l1 = fig.axes[0].lines[0]
    l2 = fig.axes[1].lines[0]
    assert l1.get_xdata() == pytest.approx(x1)
    assert l2.get_xdata() == pytest.approx(x2)
    assert l1.get_ydata() == pytest.approx(y1)
    assert l2.get_ydata() == pytest.approx(y2)

    l1 = fig.axes[0].lines[1]
    assert l1.get_xdata() == pytest.approx(x2)
    assert l1.get_ydata() == pytest.approx(y2)

    assert fig.axes[0].lines[0].get_alpha() is None
    assert fig.axes[1].lines[0].get_alpha() is None
    assert fig.axes[0].lines[1].get_alpha() == 0.5

    plt.close()
