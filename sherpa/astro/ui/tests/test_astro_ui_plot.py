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
Basic tests of the plot functionality in sherpa.astro.ui.

It is based on sherpa/ui/tests/test_ui_plot.py, but focusses on the
new routines, and "astro-specific" (currently PHA only) capabilities
that are in the astro layer.

There is very-limited checking that the plots are displaying the
correct data; it is more a check that the routines can be called.

"""

import logging
import numpy as np

from sherpa.astro import ui

from sherpa.astro.plot import ARFPlot, BkgDataPlot, ModelHistogram, \
    SourcePlot

from sherpa.utils.err import IdentifierErr
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_plotting, requires_pylab, requires_xspec

import pytest


_data_chan = np.linspace(1, 10, 10, dtype=np.int8)
_data_counts = np.asarray([0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
                          dtype=np.int8)

_data_bkg = np.asarray([1, 0, 0, 1, 0, 0, 2, 0, 0, 1],
                       dtype=np.int8)

_arf = np.asarray([0.8, 0.8, 0.9, 1.0, 1.1, 1.1, 0.7, 0.6, 0.6, 0.6])

# Using a "perfect" RMF, in that there's a one-to-one mapping
# from channel to energy. I use a non-uniform grid to make
# it obvious when the bin width is being used: when using a
# constant bin width like 0.1 keV the factor of 10 is too easy
# to confuse for other terms.
#
_energies = np.linspace(0.5, 1.5, 11)
_energies = np.asarray([0.5, 0.65, 0.75, 0.8, 0.9, 1. , 1.1, 1.12, 1.3, 1.4, 1.5])
_energies_lo = _energies[:-1]
_energies_hi = _energies[1:]
_energies_mid = (_energies_lo + _energies_hi) / 2
_energies_width = _energies_hi - _energies_lo

# How much longer is the background exposure compared to the source
# exposure; chose a non-integer value to make it more obvious when
# it is being applied (e.g. compared to the scaling ratios applied
# to the backscal values).
#
_bexpscale = 2.5

# Make sure the arrays can't be changed
for _array in [_data_chan, _data_counts, _data_bkg, _arf, _energies]:
    _array.flags.writeable = False

del _array

# Normalisation of the models.
#
MODEL_NORM = 1.02e2
BGND_NORM = 0.4


def example_pha_data():
    """Create an example data set."""

    etime = 1201.0
    d = ui.DataPHA('example', _data_chan.copy(),
                   _data_counts.copy(),
                   exposure=etime,
                   backscal=0.2)

    a = ui.create_arf(_energies_lo.copy(),
                      _energies_hi.copy(),
                      specresp=_arf.copy(),
                      exposure=etime)

    r = ui.create_rmf(_energies_lo.copy(),
                      _energies_hi.copy(),
                      e_min=_energies_lo.copy(),
                      e_max=_energies_hi.copy(),
                      startchan=1,
                      fname=None)

    d.set_arf(a)
    d.set_rmf(r)
    return d


def example_pha_with_bkg_data(direct=True):
    """Create an example data set with background

    There is no response for the background.

    Parameters
    ----------
    direct : bool, optional
        If True then the background is added to the source
        and only the source is returned. If False then the
        return is (source, bgnd) and the background is not
        associated with the source.

    """

    d = example_pha_data()

    b = ui.DataPHA('example-bkg', _data_chan.copy(),
                   _data_bkg.copy(),
                   exposure=1201.0 * _bexpscale,
                   backscal=0.4)

    if direct:
        d.set_background(b)
        return d

    return d, b


def example_model():
    """Create an example model."""

    ui.create_model_component('const1d', 'cpt')
    cpt = ui.get_model_component('cpt')
    cpt.c0 = MODEL_NORM
    return cpt


def example_bkg_model():
    """Create an example background model."""

    ui.create_model_component('powlaw1d', 'bcpt')
    bcpt = ui.get_model_component('bcpt')
    bcpt.gamma = 0.0  # use a flat model to make it easy to evaluate
    bcpt.ampl = BGND_NORM
    return bcpt


def setup_example(idval):
    """Set up a simple dataset for use in the tests.

    A *very basic* ARF is used, along with an ideal RMF. The
    way this is created means the analysis is in channel space
    by default.

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier.

    See Also
    --------
    setup_example_bkg
    """

    d = example_pha_data()
    m = example_model()
    if idval is None:
        ui.set_data(d)
        ui.set_source(m)

    else:
        ui.set_data(idval, d)
        ui.set_source(idval, m)


def setup_example_bkg(idval):
    """Set up a simple dataset + background for use in the tests.

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier.

    See Also
    --------
    setup_example, setup_example_bkg_model
    """

    d = example_pha_with_bkg_data()
    m = example_model()
    if idval is None:
        ui.set_data(d)
        ui.set_source(m)

    else:
        ui.set_data(idval, d)
        ui.set_source(idval, m)


def setup_example_bkg_model(idval, direct=True):
    """Set up a simple dataset + background for use in the tests.

    This includes a model for the background, unlike
    setup-example_bkg.

    Parameters
    ----------
    idval : None, int, str
        The dataset identifier.

    See Also
    --------
    setup_example_bkg
    """

    d = example_pha_with_bkg_data(direct=direct)
    m = example_model()
    bm = example_bkg_model()
    if idval is None:
        if direct:
            ui.set_data(d)
        else:
            ui.set_data(d[0])
            ui.set_bkg(d[1])

        ui.set_source(m)
        ui.set_bkg_model(bm)

    else:
        if direct:
            ui.set_data(idval, d)
        else:
            ui.set_data(idval, d[0])
            ui.set_bkg(idval, d[1])
        ui.set_source(idval, m)
        ui.set_bkg_model(idval, bm)


"""
Functions that could be tested:

plot_model
plot_source
plot_model_component
plot_source_component
plot_order

plot_bkg
plot_bkg_model
plot_bkg_resid
plot_bkg_ratio
plot_bkg_delchi
plot_bkg_chisqr
plot_bkg_fit
plot_bkg_source
plot_bkg_fit_ratio
plot_bkg_fit_resid
plot_bkg_fit_delchi

plot_arf

get_arf_plot                 X
get_bkg_chisqr_plot
get_bkg_delchi_plot
get_bkg_fit_plot             X
get_bkg_model_plot           X
get_bkg_plot                 X
get_bkg_ratio_plot
get_bkg_resid_plot           X
get_bkg_source_plot
get_model_component_plot
get_model_plot               X
get_order_plot
get_source_component_plot
get_source_plot              X

"""


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_arf_plot(idval, clean_astro_ui):
    """Basic testing of get_arf_plot
    """

    setup_example(idval)
    if idval is None:
        ap = ui.get_arf_plot()
    else:
        ap = ui.get_arf_plot(idval)

    assert isinstance(ap, ARFPlot)

    assert ap.xlo == pytest.approx(_energies_lo)
    assert ap.xhi == pytest.approx(_energies_hi)

    assert ap.y == pytest.approx(_arf)

    assert ap.title == 'test-arf'
    assert ap.xlabel == 'Energy (keV)'

    # the y label depends on the backend (due to LaTeX)
    # assert ap.ylabel == 'cm$^2$'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_plot(idval, clean_astro_ui):
    """Basic testing of get_bkg_plot
    """

    setup_example_bkg(idval)
    if idval is None:
        bp = ui.get_bkg_plot()
    else:
        bp = ui.get_bkg_plot(idval)

    assert isinstance(bp, BkgDataPlot)

    assert bp.x == pytest.approx(_data_chan)

    # normalise by exposure time and bin width, but bin width here
    # is 1 (because it is being measured in channels).
    #
    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'example-bkg'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_plot_energy(idval, clean_astro_ui):
    """Basic testing of get_bkg_plot: energy
    """

    setup_example_bkg(idval)
    if idval is None:
        ui.set_analysis('energy')
        bp = ui.get_bkg_plot()
    else:
        ui.set_analysis(idval, 'energy')
        bp = ui.get_bkg_plot(idval)

    assert bp.x == pytest.approx(_data_chan)

    # normalise by exposure time and bin width,
    # and since this is incorrectly in channels there's no normalisation
    #
    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'example-bkg'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("gfunc", [ui.get_bkg_plot,
                                   ui.get_bkg_model_plot,
                                   ui.get_bkg_fit_plot])
def test_get_bkg_plot_no_bkg(idval, gfunc, clean_astro_ui):
    """Basic testing of get_bkg_XXX_plot when there's no background
    """

    setup_example(idval)
    with pytest.raises(IdentifierErr):
        if idval is None:
            gfunc()
        else:
            gfunc(idval)


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_model_plot(idval, clean_astro_ui):
    """Basic testing of get_model_plot
    """

    setup_example(idval)
    if idval is None:
        mp = ui.get_model_plot()
    else:
        mp = ui.get_model_plot(idval)

    assert isinstance(mp, ModelHistogram)

    assert mp.xlo == pytest.approx(_data_chan)
    assert mp.xhi == pytest.approx(_data_chan + 1)

    # The model is a constant, but integrated across the energy bin,
    # so the energy width is important here to get the normalization
    # right. It should also be divided by the channel width, but in
    # this case each bin has a channel width of 1.
    #
    yexp = _arf * MODEL_NORM * _energies_width
    assert mp.y == pytest.approx(yexp)

    assert mp.title == 'Model'
    assert mp.xlabel == 'Channel'
    assert mp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_model_plot_energy(idval, clean_astro_ui):
    """Basic testing of get_model_plot: energy
    """

    setup_example(idval)
    if idval is None:
        ui.set_analysis('energy')
        mp = ui.get_model_plot()
    else:
        ui.set_analysis(idval, 'energy')
        mp = ui.get_model_plot(idval)

    assert mp.xlo == pytest.approx(_energies_lo)
    assert mp.xhi == pytest.approx(_energies_hi)

    # This should be normalized by the bin width, but it is cancelled
    # out by the fact that the model normalization has to be multiplied
    # by the bin width (both in energy).
    #
    yexp = _arf * MODEL_NORM
    assert mp.y == pytest.approx(yexp)

    assert mp.title == 'Model'
    assert mp.xlabel == 'Energy (keV)'
    assert mp.ylabel == 'Counts/sec/keV'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_source_plot_warning(idval, caplog, clean_astro_ui):
    """Does get_source_plot create a warning about channel space?

    This is a logged warning, not a UserWarning.
    """

    setup_example(idval)
    if idval is None:
        ui.get_source_plot()
    else:
        ui.get_source_plot(idval)

    emsg = 'Channel space is unappropriate for the PHA unfolded ' + \
           'source model,\nusing energy.'

    assert len(caplog.record_tuples) == 1
    rec = caplog.record_tuples[0]
    assert len(rec) == 3
    loc, lvl, msg = rec

    assert loc == 'sherpa.astro.plot'
    assert lvl == logging.WARNING
    assert msg == emsg


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_source_plot_energy(idval, clean_astro_ui):
    """Basic testing of get_source_plot: energy
    """

    setup_example(idval)
    if idval is None:
        ui.set_analysis('energy')
        sp = ui.get_source_plot()
    else:
        ui.set_analysis(idval, 'energy')
        sp = ui.get_source_plot(idval)

    assert isinstance(sp, SourcePlot)

    assert sp.xlo == pytest.approx(_energies_lo)
    assert sp.xhi == pytest.approx(_energies_hi)

    yexp = MODEL_NORM * np.ones(10)
    assert sp.y == pytest.approx(yexp)

    assert sp.title == 'Source Model of example'
    assert sp.xlabel == 'Energy (keV)'

    # y label depends on the backend
    # assert sp.ylabel == 'f(E)  Photons/sec/cm$^2$/keV'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("direct", [True, False])
def test_get_bkg_model_plot(idval, direct, clean_astro_ui):
    """Basic testing of get_bkg_model_plot

    We test ui.set_bkg as well as datapha.set_background to check
    issue #879.

    The same ARF is used as the source (by construction), which is
    likely to be a common use case.
    """

    setup_example_bkg_model(idval, direct=direct)
    if idval is None:
        bp = ui.get_bkg_model_plot()
    else:
        bp = ui.get_bkg_model_plot(idval)

    assert bp.xlo == pytest.approx(_data_chan)
    assert bp.xhi == pytest.approx(_data_chan + 1)

    yexp = _arf * BGND_NORM * _energies_width
    if direct:
        yexp /= _bexpscale

    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Model'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("direct", [True, False])
def test_get_bkg_model_plot_energy(idval, direct, clean_astro_ui):
    """Basic testing of get_bkg_model_plot: energy

    We test ui.set_bkg as well as datapha.set_background,
    since I have seen subtle differences due to the extra
    logic that set_bkg can do (issues #879 and #880)
    """

    setup_example_bkg_model(idval, direct=direct)
    if idval is None:
        ui.set_analysis('energy')
        bp = ui.get_bkg_model_plot()
    else:
        ui.set_analysis(idval, 'energy')
        bp = ui.get_bkg_model_plot(idval)

    if direct:
        assert bp.xlo == pytest.approx(_data_chan)
        assert bp.xhi == pytest.approx(_data_chan + 1)
    else:
        assert bp.xlo == pytest.approx(_energies_lo)
        assert bp.xhi == pytest.approx(_energies_hi)

    yexp = _arf * BGND_NORM
    if direct:
        yexp *= (_energies_width / _bexpscale)
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Model'
    if direct:
        assert bp.xlabel == 'Channel'
        assert bp.ylabel == 'Counts/sec/channel'
    else:
        assert bp.xlabel == 'Energy (keV)'
        assert bp.ylabel == 'Counts/sec/keV'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_resid_plot(idval, clean_astro_ui):
    """Basic testing of get_bkg_resid_plot
    """

    setup_example_bkg_model(idval)
    if idval is None:
        bp = ui.get_bkg_resid_plot()
    else:
        bp = ui.get_bkg_resid_plot(idval)

    assert bp.x == pytest.approx(_data_chan)

    # correct the counts by the bin width and exposure time
    #
    yexp = _data_bkg / (1201.0 * _bexpscale) - _arf * BGND_NORM * _energies_width / _bexpscale
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Residuals of example-bkg - Bkg Model'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_resid_plot_energy(idval, clean_astro_ui):
    """Basic testing of get_bkg_resid_plot: energy
    """

    setup_example_bkg_model(idval)
    if idval is None:
        ui.set_analysis('energy')
        bp = ui.get_bkg_resid_plot()
    else:
        ui.set_analysis(idval, 'energy')
        bp = ui.get_bkg_resid_plot(idval)

    assert bp.x == pytest.approx(_data_chan)

    # correct the counts by the bin width (which is 1) and exposure time
    #
    yexp = _data_bkg / (1201.0 * _bexpscale) - _arf * BGND_NORM * _energies_width / _bexpscale
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Residuals of example-bkg - Bkg Model'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_fit_plot(idval, clean_astro_ui):
    """Basic testing of get_bkg_fit_plot
    """

    setup_example_bkg_model(idval)
    if idval is None:
        fp = ui.get_bkg_fit_plot()
    else:
        fp = ui.get_bkg_fit_plot(idval)

    dp = fp.dataplot
    mp = fp.modelplot

    assert dp.title == 'example-bkg'
    assert mp.title == 'Background Model Contribution'

    for plot in [dp, mp]:
        assert plot.xlabel == 'Channel'
        assert plot.ylabel == 'Counts/sec/channel'
        assert plot.x == pytest.approx(_data_chan)

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dp.y == pytest.approx(dp.y)

    yexp = _arf * BGND_NORM * _energies_width / _bexpscale
    assert mp.y == pytest.approx(yexp)


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_bkg_fit_plot_energy(idval, clean_astro_ui):
    """Basic testing of get_bkg_fit_plot: energy
    """

    setup_example_bkg_model(idval)
    if idval is None:
        ui.set_analysis('energy')
        fp = ui.get_bkg_fit_plot()
    else:
        ui.set_analysis(idval, 'energy')
        fp = ui.get_bkg_fit_plot(idval)

    dp = fp.dataplot
    mp = fp.modelplot

    assert dp.title == 'example-bkg'
    assert mp.title == 'Background Model Contribution'

    for plot in [dp, mp]:
        assert plot.xlabel == 'Channel'
        assert plot.ylabel == 'Counts/sec/channel'
        assert plot.x == pytest.approx(_data_chan)

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dp.y == pytest.approx(dp.y)

    yexp = _arf * BGND_NORM * _energies_width / _bexpscale
    assert mp.y == pytest.approx(yexp)


def check_bkg_fit(plotfunc):
    """Is the background fit displayed?

    This only checks the plot object, not the plot "hardcopy" output
    (e.g. the pixel display/PNG output).
    """

    dplot = ui._session._bkgdataplot
    mplot = ui._session._bkgmodelplot

    # check the "source" plots are not set
    for plot in [ui._session._dataplot, ui._session._modelplot]:
        assert plot.x is None
        assert plot.y is None

    xlabel = 'Channel' if plotfunc == ui.plot_bkg_fit else ''

    # check plot basics
    for plot in [dplot, mplot]:
        assert plot.xlabel == xlabel
        assert plot.ylabel == 'Counts/sec/channel'

        assert plot.x == pytest.approx(_data_chan)

    assert dplot.title == 'example-bkg'
    assert mplot.title == 'Background Model Contribution'

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dplot.y == pytest.approx(yexp)

    yexp = _arf * BGND_NORM * _energies_width / _bexpscale
    assert mplot.y == pytest.approx(yexp)


def check_bkg_resid(plotfunc):
    """Is the background residual displayed?

    This only checks the plot object, not the plot "hardcopy" output
    (e.g. the pixel display/PNG output).

    There is limited checks of the actual values, since the
    assumption is that this test is just to check the plots were
    constructed from its components, and that other tests above
    will have checked that the components work.

    """

    # check the "other" background plots are not set
    plot = None
    for pf, pd in [(ui.plot_bkg_fit_delchi, ui._session._bkgdelchiplot),
                   (ui.plot_bkg_fit_ratio, ui._session._bkgratioplot),
                   (ui.plot_bkg_fit_resid, ui._session._bkgresidplot)]:
        if pf == plotfunc:
            assert plot is None  # a precaution
            plot = pd
            continue
        else:
            assert pd.x is None
            assert pd.y is None

    # very limited checks
    #
    assert plot.xlabel == 'Channel'
    assert plot.ylabel != ''  # depends on the plot type
    assert plot.title == ''

    assert plot.x == pytest.approx(_data_chan)
    assert plot.y is not None

    # the way the data and model are constructed, all residual values
    # should be negative, and the ratio values positive (or zero).
    #
    if plotfunc == ui.plot_bkg_fit_ratio:
        assert np.all(plot.y >= 0)
    else:
        assert np.all(plot.y < 0)


@requires_plotting
@pytest.mark.usefixtures("clean_astro_ui")
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,checkfuncs",
                         [(ui.plot_bkg_fit, [check_bkg_fit]),
                          (ui.plot_bkg_fit_delchi, [check_bkg_fit, check_bkg_resid]),
                          (ui.plot_bkg_fit_ratio, [check_bkg_fit, check_bkg_resid]),
                          (ui.plot_bkg_fit_resid, [check_bkg_fit, check_bkg_resid])])
def test_bkg_plot_xxx(idval, plotfunc, checkfuncs):
    """Test background plotting - channel space"""

    setup_example_bkg_model(idval)
    if idval is None:
        plotfunc()
    else:
        plotfunc(idval)

    # The X label of the plots may depend on whether there are 1
    # or two plots. The following isn't ideal but let's see how
    # it goes.
    #
    for checkfunc in checkfuncs:
        checkfunc(plotfunc)


# The following tests were added in a separate PR to those above, and
# rather than try to work out whether they do the same thing, they
# have been left separate.
#

@pytest.fixture
def basic_pha1(make_data_path):
    """Create a basic PHA-1 data set/setup"""

    ui.set_default_id('tst')
    ui.load_pha(make_data_path('3c273.pi'))
    ui.subtract()
    ui.notice(0.5, 7)
    ui.set_source(ui.powlaw1d.pl)
    pl = ui.get_model_component('pl')
    pl.gamma = 1.93
    pl.ampl = 1.74e-4


@pytest.fixture
def basic_img(make_data_path):
    """Create a basic image data set/setup"""

    ui.set_default_id(2)
    ui.load_image(make_data_path('img.fits'))
    ui.set_source(ui.gauss2d.gmdl)
    ui.guess()


_basic_plotfuncs = [ui.plot_data,
                    ui.plot_bkg,
                    ui.plot_model,
                    ui.plot_source,
                    ui.plot_resid,
                    ui.plot_delchi,
                    ui.plot_ratio,
                    ui.plot_fit,
                    ui.plot_fit_delchi,
                    ui.plot_fit_resid,
                    ui.plot_arf,
                    ui.plot_chisqr]


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_function(clean_astro_ui, basic_pha1):
    # can we call plot; do not try to be exhaustive
    ui.plot("data", "bkg", "fit", "arf")


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", _basic_plotfuncs)
def test_pha1_plot(clean_astro_ui, basic_pha1, plotfunc):
    plotfunc()


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", [ui.plot_model_component,
                                      ui.plot_source_component])
def test_pha1_plot_component(clean_astro_ui, basic_pha1, plotfunc):
    plotfunc("pl")


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", [ui.int_unc, ui.int_proj])
def test_pha1_int_plot(clean_astro_ui, basic_pha1, plotfunc):
    plotfunc('pl.gamma')


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", [ui.reg_unc, ui.reg_proj])
def test_pha1_reg_plot(clean_astro_ui, basic_pha1, plotfunc):
    plotfunc('pl.gamma', 'pl.ampl')


_img_plotfuncs = [ui.contour_data,
                  ui.contour_fit,
                  ui.contour_fit_resid,
                  ui.contour_model,
                  ui.contour_ratio,
                  ui.contour_resid,
                  ui.contour_source ]


@requires_plotting
@requires_fits
@requires_data
def test_img_contour_function(clean_astro_ui, basic_img):
    # can we call contour; do not try to be exhaustive
    ui.contour("data", "model", "source", "fit")


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", _img_plotfuncs)
def test_img_contour(clean_astro_ui, basic_img, plotfunc):
    plotfunc()


# Add in some pylab-specific tests to change default values
#

@requires_pylab
@requires_fits
@requires_data
def test_pha1_plot_data_options(clean_astro_ui, basic_pha1):
    """Test that the options have changed things, where easy to do so"""

    from matplotlib import pyplot as plt
    import matplotlib

    prefs = ui.get_data_plot_prefs()

    # check the preference are as expected for the boolean cases
    assert not prefs['xerrorbars']
    assert prefs['yerrorbars']
    assert not prefs['xlog']
    assert not prefs['ylog']

    prefs['xerrorbars'] = True
    prefs['yerrorbars'] = False
    prefs['xlog'] = True
    prefs['ylog'] = True

    prefs['color'] = 'orange'

    prefs['linecolor'] = 'brown'
    prefs['linestyle'] = '-.'

    prefs['marker'] = 's'
    prefs['markerfacecolor'] = 'cyan'
    prefs['markersize'] = 10

    ui.plot_data()

    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'log'

    assert ax.get_xlabel() == 'Energy (keV)'
    assert ax.get_ylabel() == 'Counts/sec/keV'

    # It is not clear whether an 'exact' check on the value, as
    # provided by pytest.approx, makes sense, or whether a "softer"
    # check - e.g.  just check whether it is less- or greater- than a
    # value - should be used. It depends on how often matplotlib
    # tweaks the axis settings and how sensitive it is to
    # platform/backend differences. Let's see how pytest.approx works
    #
    xmin, xmax = ax.get_xlim()
    assert xmin == pytest.approx(0.40110954270367555)
    assert xmax == pytest.approx(11.495805054836712)

    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(7.644069935298475e-05)
    assert ymax == pytest.approx(0.017031102671151491)

    assert len(ax.lines) == 1
    line = ax.lines[0]

    # Apparently color wins out over linecolor
    assert line.get_color() == 'orange'
    assert line.get_linestyle() == '-.'
    assert line.get_marker() == 's'
    assert line.get_markerfacecolor() == 'cyan'
    assert line.get_markersize() == pytest.approx(10.0)

    # assume error bars handled by a collection; test a subset
    # of values
    #
    assert len(ax.collections) == 1
    coll = ax.collections[0]

    assert len(coll.get_segments()) == 42

    # The return value depends on matplotlib version (>= 3.3
    # returns something). What has changed? Maybe this should
    # not be tested?
    #
    expected = [(None, None)]
    if matplotlib.__version__ >= '3.3.0':
        expected = [(0.0, None)]

    assert coll.get_linestyles() == expected

    # looks like the color has been converted to individual channels
    # - e.g. floating-point values for R, G, B, and alpha.
    #
    colors = coll.get_color()
    assert len(colors) == 1
    assert len(colors[0]) == 4
    r, g, b, a = colors[0]
    assert r == pytest.approx(1)
    assert g == pytest.approx(0.64705882)
    assert b == pytest.approx(0)
    assert a == pytest.approx(1)


@requires_pylab
@requires_fits
@requires_data
def test_pha1_plot_model_options(clean_astro_ui, basic_pha1):
    """Test that the options have changed things, where easy to do so

    In matplotlib 3.1 the plot_model call causes a MatplotlibDeprecationWarning
    to be created:

    Passing the drawstyle with the linestyle as a single string is deprecated since Matplotlib 3.1 and support will be removed in 3.3; please pass the drawstyle separately using the drawstyle keyword argument to Line2D or set_drawstyle() method (or ds/set_ds()).

    This warning is hidden by the test suite (sherpa/conftest.py) so that
    it doesn't cause the tests to fail. Note that a number of other tests
    in this module also cause this warning to be displayed.

    """

    from matplotlib import pyplot as plt

    # Note that for PHA data sets, the mode is drawn as a histogram,
    # so get_model_plot_prefs doesn't actually work. We need to change
    # the histogram prefs instead. See issue
    # https://github.com/sherpa/sherpa/issues/672
    #
    # prefs = ui.get_model_plot_prefs()
    prefs = ui.get_model_plot().histo_prefs

    # check the preference are as expected for the boolean cases
    assert not prefs['xlog']
    assert not prefs['ylog']

    # Only change the X axis here
    prefs['xlog'] = True

    prefs['color'] = 'green'

    prefs['linecolor'] = 'red'
    prefs['linestyle'] = 'dashed'

    prefs['marker'] = '*'
    prefs['markerfacecolor'] = 'yellow'
    prefs['markersize'] = 8

    ui.plot_model()

    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'linear'

    assert ax.get_xlabel() == 'Energy (keV)'
    assert ax.get_ylabel() == 'Counts/sec/keV'

    # It is not clear whether an 'exact' check on the value, as
    # provided by pytest.approx, makes sense, or whether a "softer"
    # check - e.g.  just check whether it is less- or greater- than a
    # value - should be used. It depends on how often matplotlib
    # tweaks the axis settings and how sensitive it is to
    # platform/backend differences. Let's see how pytest.approx works
    #
    xmin, xmax = ax.get_xlim()
    assert xmin == pytest.approx(0.40770789163447285)
    assert xmax == pytest.approx(11.477975806572461)

    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(-0.00045772936258082011, rel=0.01)
    assert ymax == pytest.approx(0.009940286575890335, rel=0.01)

    assert len(ax.lines) == 1
    line = ax.lines[0]

    # Apparently color wins out over linecolor
    assert line.get_color() == 'green'
    assert line.get_linestyle() == '--'  # note: input was dashed
    assert line.get_marker() == '*'
    assert line.get_markerfacecolor() == 'yellow'
    assert line.get_markersize() == pytest.approx(8.0)

    assert len(ax.collections) == 0


@requires_pylab
@requires_fits
@requires_data
def test_pha1_plot_fit_options(clean_astro_ui, basic_pha1):
    """Test that the options have changed things, where easy to do so"""

    from matplotlib import pyplot as plt
    import matplotlib

    dprefs = ui.get_data_plot_prefs()
    dprefs['xerrorbars'] = True
    dprefs['yerrorbars'] = False
    dprefs['xlog'] = True
    dprefs['ylog'] = True
    dprefs['color'] = 'orange'
    dprefs['linestyle'] = '-.'
    dprefs['marker'] = 's'
    dprefs['markerfacecolor'] = 'cyan'
    dprefs['markersize'] = 10

    mprefs = ui.get_model_plot().histo_prefs
    mprefs['color'] = 'green'
    mprefs['linestyle'] = 'dashed'
    mprefs['marker'] = '*'
    mprefs['markerfacecolor'] = 'yellow'
    mprefs['markersize'] = 8

    ui.plot_fit(alpha=0.7)

    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'log'

    assert ax.get_xlabel() == 'Energy (keV)'
    assert ax.get_ylabel() == 'Counts/sec/keV'

    xmin, xmax = ax.get_xlim()
    assert xmin == pytest.approx(0.40110954270367555)
    assert xmax == pytest.approx(11.495805054836712)

    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(7.644069935298475e-05)
    assert ymax == pytest.approx(0.017031102671151491)

    assert len(ax.lines) == 2

    # DATA
    #
    line = ax.lines[0]

    # Apparently color wins out over linecolor
    assert line.get_color() == 'orange'
    assert line.get_linestyle() == '-.'
    assert line.get_marker() == 's'
    assert line.get_markerfacecolor() == 'cyan'
    assert line.get_markersize() == pytest.approx(10.0)
    assert line.get_alpha() == pytest.approx(0.7)

    # MODEL
    #
    line = ax.lines[1]
    assert line.get_color() != 'green'  # option over-ridden
    assert line.get_linestyle() == '-'  # option over-ridden
    assert line.get_marker() == 'None'
    assert line.get_markerfacecolor() == line.get_color()  # option over-ridden
    assert line.get_markersize() == pytest.approx(6.0)
    assert line.get_alpha() == pytest.approx(0.7)

    # assume error bars handled by a collection; test a subset
    # of values
    #
    assert len(ax.collections) == 1
    coll = ax.collections[0]

    assert len(coll.get_segments()) == 42

    # The return value depends on matplotlib version (>= 3.3
    # returns something). What has changed? Maybe this should
    # not be tested?
    #
    expected = [(None, None)]
    if matplotlib.__version__ >= '3.3.0':
        expected = [(0.0, None)]

    assert coll.get_linestyles() == expected

    # looks like the color has been converted to individual channels
    # - e.g. floating-point values for R, G, B, and alpha.
    #
    colors = coll.get_color()
    assert len(colors) == 1
    assert len(colors[0]) == 4
    r, g, b, a = colors[0]
    assert r == pytest.approx(1)
    assert g == pytest.approx(0.64705882)
    assert b == pytest.approx(0)
    assert a == pytest.approx(0.7)


@requires_pylab
@requires_fits
@requires_data
@requires_xspec
def test_pha1_reg_proj(clean_astro_ui, basic_pha1):
    """This is potentially a time-consuming test to run, so simplify
    as much as possible.
    """

    from matplotlib import pyplot as plt

    pl = ui.get_model_component("pl")
    ui.set_source(ui.xsphabs.gal * pl)
    gal = ui.get_model_component("gal")

    ui.fit()

    ui.reg_proj("pl.gamma", "gal.nh", min=(1.6, 0), max=(2.5, 0.2),
                nloop=(3, 3))

    ax = plt.gca()
    assert ax.get_xscale() == 'linear'
    assert ax.get_yscale() == 'linear'

    assert ax.get_xlabel() == 'pl.gamma'
    assert ax.get_ylabel() == 'gal.nH'
    assert ax.get_title() == 'Region-Projection'

    xmin, xmax = ax.get_xlim()
    assert xmin == pytest.approx(1.6)
    assert xmax == pytest.approx(2.5)

    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(0.0)
    assert ymax == pytest.approx(0.2)

    assert len(ax.lines) == 1
    line = ax.lines[0]
    assert line.get_xdata().size == 1

    x0 = line.get_xdata()[0]
    y0 = line.get_ydata()[0]

    assert x0 == pytest.approx(pl.gamma.val)
    assert y0 == pytest.approx(gal.nh.val)

    # pylab get_confid_point_defaults() returns
    # {'symbol': '+', 'color': None}
    #
    assert line.get_marker() == '+'

    # the number depends on the matplotlib version: 2 for 2.2.3 and
    # 3 for 3.1.1; it's not clear what the "extra" one is in matplotlib 3
    # (it isn't obviously visible). DJB guesses that this would be
    # clearer if we ran with more bins along each axis, but this would
    # take more time.
    #
    ncontours = len(ax.collections)
    assert ncontours in [2, 3]


DATA_PREFS = {'alpha': None,
              'barsabove': False,
              'capsize': None,
              'color': None,
              'drawstyle': 'default',
              'ecolor': None,
              'linecolor': None,
              'linestyle': 'None',
              'marker': '.',
              'markerfacecolor': None,
              'markersize': None,
              'ratioline': False,
              'xaxis': False,
              'xerrorbars': False,
              'xlog': False,
              'yerrorbars': True,
              'ylog': False}

MODEL_PREFS = {'alpha': None,
               'barsabove': False,
               'capsize': None,
               'color': None,
               'drawstyle': 'default',
               'ecolor': None,
               'linecolor': None,
               'linestyle': '-',
               'marker': 'None',
               'markerfacecolor': None,
               'markersize': None,
               'ratioline': False,
               'xaxis': False,
               'xerrorbars': False,
               'xlog': False,
               'yerrorbars': False,
               'ylog': False}


CONTOUR_PREFS = {'alpha': None,
                 'colors': None,
                 'linewidths': None,
                 'xlog': False,
                 'ylog': False}


@requires_pylab
@pytest.mark.parametrize("funcname,expected",
                         [("get_data_plot_prefs", DATA_PREFS),
                          ("get_model_plot_prefs", MODEL_PREFS),
                          ("get_data_contour_prefs", CONTOUR_PREFS),
                          ("get_model_contour_prefs", CONTOUR_PREFS)])
def test_plot_defaults(funcname, expected):
    """What are the plot defaults?"""

    prefs = getattr(ui, funcname)()
    assert isinstance(prefs, dict)
    assert prefs == expected


@requires_fits
@requires_data
def test_pha1_get_model_plot_filtered(clean_astro_ui, basic_pha1):
    """Does get_model_plot register filters correctly?

    If there is an ignored range within a model, is it ignored?
    This also holds for plot_model and the model_component
    variants, but this is untested
    """

    # Ignore the 2-3 keV range
    ui.ignore(2, 3)

    mplot = ui.get_model_plot()

    assert mplot.xlo.size == mplot.xhi.size
    assert mplot.xlo.size == mplot.y.size

    assert mplot.y.size == 566

    assert mplot.xlo[0] == pytest.approx(0.46720001101493835)
    assert mplot.xhi[0] == pytest.approx(0.48179998993873596)

    assert mplot.xlo[-1] == pytest.approx(9.854999542236328)
    assert mplot.xhi[-1] == pytest.approx(9.869600296020508)

    assert mplot.xlo[100] == pytest.approx(1.9271999597549438)
    assert mplot.xhi[100] == pytest.approx(1.9417999982833862)

    assert mplot.xlo[101] == pytest.approx(3.0806000232696533)
    assert mplot.xhi[101] == pytest.approx(3.0952000617980957)
