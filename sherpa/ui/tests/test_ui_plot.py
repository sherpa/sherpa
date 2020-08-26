#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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

from sherpa.data import Data1DInt, Data2D
from sherpa.models import basic
from sherpa.plot import CDFPlot, DataPlot, FitPlot, ModelPlot, \
    PDFPlot, PSFPlot, PSFKernelPlot, ScatterPlot, TracePlot,\
    DataContour, ModelContour, SourceContour, ResidContour, \
    RatioContour, FitContour, PSFContour, LRHistogram

from sherpa.stats import Chi2Gehrels
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr
from sherpa.utils.testing import requires_plotting, requires_pylab


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

    assert isinstance(f, FitPlot)
    assert isinstance(f.dataplot, DataPlot)
    assert isinstance(f.modelplot, ModelPlot)

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
    assert isinstance(plotobj, PSFPlot)
    assert plotobj.title == 'gauss1d.psfmdl'
    assert plotobj.xlabel == 'x'
    assert plotobj.ylabel == 'y'
    assert plotobj.xerr is None
    assert plotobj.yerr is None
    assert plotobj.x == pytest.approx(x)

    yexp = psfmdl(x)
    assert plotobj.y == pytest.approx(yexp)


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_kernel(session, caplog):
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
    assert isinstance(plotobj, PSFKernelPlot)
    assert plotobj.title == 'PSF Kernel'
    assert plotobj.xlabel == 'x'
    assert plotobj.ylabel == 'y'
    assert plotobj.xerr is None
    assert plotobj.yerr is None
    assert plotobj.x == pytest.approx(x)

    yexp = psfmdl(x)
    yexp /= yexp.sum()
    assert plotobj.y == pytest.approx(yexp)


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("plotobj,plotfunc",
                         [("_residplot", ui.plot_resid),
                          ("_delchiplot", ui.plot_delchi)
                          ])
def test_plot_resid_ignores_ylog(plotobj, plotfunc):
    """Do the plot_resid-family of routines ignore the ylog setting?

    Note that plot_chisqr is not included in support for ignoring
    ylog (since the data should be positive in this case).
    """

    # access it this way to ensure have access to the actual
    # object used by the session object in this test (as the
    # clean call done by the clean_ui fixture will reset the
    # plot objects).
    #
    prefs = getattr(ui._session, plotobj).plot_prefs

    setup_example(None)

    ui.set_ylog()
    assert prefs['ylog']

    plotfunc(ylog=True)

    # Note that the ylog setting has been removed (to reflect
    # what was displayed).
    #
    assert not prefs['ylog']


@requires_plotting
@pytest.mark.usefixtures("clean_ui")
@pytest.mark.parametrize("plotobj,plotfunc",
                         [("_residplot", ui.plot_fit_resid),
                          ("_delchiplot", ui.plot_fit_delchi)
                          ])
def test_plot_fit_resid_ignores_ylog(plotobj, plotfunc):
    """Do the plot_resid-family of routines ignore the ylog setting?"""

    # access it this way to ensure have access to the actual
    # object used by the session object in this test (as the
    # clean call done by the clean_ui fixture will reset the
    # plot objects).
    #
    rprefs = getattr(ui._session, plotobj).plot_prefs
    dprefs = ui._session._dataplot.plot_prefs

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
    because of an invaid identifier.
    """

    s = session()
    func = getattr(s, plottype)
    with pytest.raises(ArgumentErr) as exc:
        func("fooflan flim flam")

    emsg = "'fooflan flim flam' is not a valid plot type"
    assert str(exc.value) == emsg


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_single(session):
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

    assert ax.get_geometry() == (1, 1, 1)
    assert ax.get_title() == ''
    assert ax.xaxis.get_label().get_text() == 'x'
    assert ax.yaxis.get_label().get_text() == 'y'

    plt.close()

    s.plot("model", 1)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_geometry() == (1, 1, 1)
    assert ax.get_title() == 'Model'
    assert ax.xaxis.get_label().get_text() == 'x'
    assert ax.yaxis.get_label().get_text() == 'y'

    plt.close()


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_multiple(session):
    """Can we call plot() with multiple plot types?

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
           "ratio", "tst")

    # Note: I wanted to try source_component but it is not
    # clear it works when the id is not the default.

    fig = plt.gcf()
    assert len(fig.axes) == 5

    for i, (ax, title, ylabel) in enumerate(zip(fig.axes,
                                                ['', 'Model', '', 'Source',
                                                 'Ratio of Data to Model'],
                                                ['y', 'y', 'y', 'y',
                                                 'Data / Model']),
                                    1):

        assert ax.get_geometry() == (2, 3, i)
        assert ax.get_title() == title
        assert ax.xaxis.get_label().get_text() == 'x'
        assert ax.yaxis.get_label().get_text() == ylabel

    plt.close()


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_contour_single(session):
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

    assert ax.get_geometry() == (1, 1, 1)
    assert ax.get_title() == ''
    assert ax.xaxis.get_label().get_text() == 'x0'
    assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()

    s.contour("model", 1)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = fig.axes[0]

    assert ax.get_geometry() == (1, 1, 1)
    assert ax.get_title() == 'Model'
    assert ax.xaxis.get_label().get_text() == 'x0'
    assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_contour_multiple(session):
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

    for i, (ax, title) in enumerate(zip(fig.axes,
                                        ['', 'Model', 'Source', '',
                                         'Ratio of Data to Model']),
                                    1):

        assert ax.get_geometry() == (2, 3, i)
        assert ax.get_title() == title
        assert ax.xaxis.get_label().get_text() == 'x0'
        assert ax.yaxis.get_label().get_text() == 'x1'

    plt.close()


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("plotfunc,title,pcls",
                         [("data", "", DataContour),
                          ("model", "Model", ModelContour),
                          ("source", "Source", SourceContour),
                          ("resid", "Residuals", ResidContour),
                          ("ratio", "Ratio of Data to Model", RatioContour),
                          ("fit", "", FitContour),
                          ("fit_resid", None, None)])
def test_contour_xxx(plotfunc, title, pcls, session):
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

        for i, (ax, title) in enumerate(zip(fig.axes,
                                            ['', 'Residuals']),
                                        1):

            assert ax.get_geometry() == (2, 1, i)
            assert ax.get_title() == title
            assert ax.xaxis.get_label().get_text() == 'x0'
            assert ax.yaxis.get_label().get_text() == 'x1'

    else:
        assert len(fig.axes) == 1

        ax = fig.axes[0]
        assert ax.get_geometry() == (1, 1, 1)
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


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
@pytest.mark.parametrize("ptype", ["data", "model"])
def test_get_xxx_contour_prefs_pylab(ptype, session):

    s = session()
    p = getattr(s, "get_{}_contour_prefs".format(ptype))()
    assert isinstance(p, dict)
    assert p == {'xlog': False, 'ylog': False,
                 'alpha': None, 'linewidths': None, 'colors': None}


@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_scatter_plot_empty(session):
    """Very basic check we can call get_scatter_plot

    This can be run even without a plotting backend available.
    """

    s = session()
    p = s.get_scatter_plot()
    assert isinstance(p, ScatterPlot)
    for f in ['x', 'y', 'xerr', 'yerr', 'xlabel', 'ylabel', 'title']:
        assert getattr(p, f) is None


@requires_plotting
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
    assert isinstance(p, ScatterPlot)

    assert p.x == pytest.approx(x)
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'x'
    assert p.ylabel == 'y'
    assert p.title == 'Scatter: (x,y)'


@requires_plotting
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
    assert isinstance(p, ScatterPlot)

    assert p.x == pytest.approx(x)
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'a x'
    assert p.ylabel == '43'
    assert p.title == 'Scatter: (x,y)'


@requires_plotting
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
    assert isinstance(p, ScatterPlot)

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
    assert isinstance(p, TracePlot)
    for f in ['x', 'y', 'xerr', 'yerr', 'xlabel', 'ylabel', 'title']:
        assert getattr(p, f) is None


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_get_trace_plot(session):
    """Very basic check we can call get_trace_plot/plot_trace

    This can be run even without a plotting backend available.
    """

    s = session()

    y = [-5, 2, -3]
    s.plot_trace(y)

    p = s.get_trace_plot()
    assert isinstance(p, TracePlot)

    assert p.x == pytest.approx([0, 1, 2])
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'iteration'
    assert p.ylabel == 'x'
    assert p.title == 'Trace: x'


@requires_plotting
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
    assert isinstance(p, TracePlot)

    assert p.x == pytest.approx([0, 1, 2])
    assert p.y == pytest.approx(y)
    assert p.xerr is None
    assert p.yerr is None
    assert p.xlabel == 'iteration'
    assert p.ylabel == 'x'
    assert p.title == 'Trace: x'


@requires_plotting
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
    assert isinstance(p, TracePlot)

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
    assert isinstance (p, CDFPlot)
    for f in ["points", "x", "y", "median", "lower", "upper", "xlabel", "ylabel", "title"]:
        assert getattr(p, f) is None


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_replot_no_data(session):
    """what does replot=True do for plot_cdf?

    The code doesn't check for this evantuality,
    so errors out (depends on the backend).
    """

    s = session()

    x = np.asarray([2, 8, 4, 6])

    # error can depend on matplotlib version
    with pytest.raises(ValueError):
        s.plot_cdf(x, replot=True)


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_replot(session):
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


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x)

    p = s.get_cdf_plot()
    assert isinstance (p, CDFPlot)

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


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_labels_noname(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x, xlabel='a b')

    p = s.get_cdf_plot()
    assert isinstance (p, CDFPlot)

    assert p.xlabel == 'a b'
    assert p.ylabel == 'p(<=a b)'
    assert p.title == 'CDF: x'


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_cdf_labels(session):

    s = session()

    x = np.asarray([20, 10, 14, 15, 12, 16, 17])
    s.plot_cdf(x, xlabel='a b', name='b a')

    p = s.get_cdf_plot()
    assert isinstance (p, CDFPlot)

    assert p.xlabel == 'a b'
    assert p.ylabel == 'p(<=a b)'
    assert p.title == 'CDF: b a'


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_cdf_plot_empty(session):

    s = session()

    p = s.get_cdf_plot()
    assert isinstance (p, CDFPlot)

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


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_cdf_plot(session, old_numpy_printing):
    """This is issue #912

    The display of the numeric values can depend on the
    NumPy version, so force the legacy output.
    """

    s = session()

    x = np.asarray([20, 15, 25, 10])
    s.plot_cdf(x)

    p = s.get_cdf_plot()
    assert isinstance (p, CDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 10

    assert toks[0] == 'points = [10,15,20,25]'
    assert toks[1] == 'x      = [ 0.25, 0.5 , 0.75, 1.  ]'
    assert toks[2] == 'y      = [20,15,25,10]'
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
    assert isinstance (p, PDFPlot)
    for f in ["points", "xlo", "xhi", "y", "xlabel", "ylabel", "title"]:
        assert getattr(p, f) is None


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x)

    p = s.get_pdf_plot()
    assert isinstance (p, PDFPlot)

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


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_replot_no_data(session):
    """what does replot=True do for plot_pdf?

    The code doesn't check for this evantuality,
    so errors out (depends on the backend).
    """

    s = session()

    x = np.asarray([2, 8, 4, 6])

    # check on the error so we know when the code has changed
    with pytest.raises(TypeError) as exc:
        s.plot_pdf(x, replot=True)

    assert str(exc.value) == "unsupported operand type(s) for +: 'NoneType' and 'NoneType'"


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_replot(session):
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

    assert len(ax.lines) == 1
    line = ax.lines[0]

    x = np.asarray([1.4, 2.2, 3, 3.8, 4.6])
    y = np.asarray([0.3125, 0.3125, 0.3125, 0, 0.3125])

    assert line.get_xdata() == pytest.approx(x)
    assert line.get_ydata() == pytest.approx(y)

    plt.close()


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_labels_noname(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x, xlabel='x^2 x', bins=4)

    p = s.get_pdf_plot()
    assert isinstance (p, PDFPlot)

    assert p.xlo.size == 4
    assert p.xhi.size == 4
    assert p.y.size == 4

    assert p.xlabel == 'x^2 x'
    assert p.ylabel == 'probability density'
    assert p.title == 'PDF: x'


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pdf_labels(session):

    s = session()

    x = np.asarray([2, 8, 4, 6])
    s.plot_pdf(x, xlabel='x^2 x', name='no name', bins=4)

    p = s.get_pdf_plot()
    assert isinstance (p, PDFPlot)

    assert p.xlo.size == 4
    assert p.xhi.size == 4
    assert p.y.size == 4

    assert p.xlabel == 'x^2 x'
    assert p.ylabel == 'probability density'
    assert p.title == 'PDF: no name'


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_show_pdf_plot_empty(session):

    s = session()

    p = s.get_pdf_plot()
    assert isinstance (p, PDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 8

    assert toks[0] == 'points = None'
    assert toks[1] == 'xlo    = None'
    assert toks[2] == 'xhi    = None'
    assert toks[3] == 'y      = None'
    assert toks[4] == 'xlabel = None'
    assert toks[5] == 'ylabel = None'
    assert toks[6] == 'title  = None'


@requires_plotting
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
    assert isinstance (p, PDFPlot)

    toks = str(p).split('\n')
    assert len(toks) == 8

    assert toks[0] == 'points = [20,15,25,10]'
    assert toks[1] == 'xlo    = [ 10., 15., 20.]'
    assert toks[2] == 'xhi    = [ 15., 20., 25.]'
    assert toks[3] == 'y      = [1,1,2]'
    assert toks[4] == 'xlabel = x'
    assert toks[5] == 'ylabel = probability density'
    assert toks[6] == 'title  = PDF: x'


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

    assert isinstance(plot, LRHistogram)
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


@requires_plotting
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_plot_pvalue(session, caplog):
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


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_scatter_empty_replot(session):
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


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_scatter(session):
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


@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_trace_empty_replot(session):
    """plot_trace with replot=False and no data

    Just check the current behavior
    """

    s = session()

    y = np.arange(100, 104)

    # error can depend on matplotlib version
    with pytest.raises(ValueError):
        s.plot_trace(y, replot=True)



@requires_pylab
@pytest.mark.parametrize("session", [BaseSession, AstroSession])
def test_pylab_plot_trace(session):
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
