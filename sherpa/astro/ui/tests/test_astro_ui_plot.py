#
#  Copyright (C) 2019, 2020, 2021, 2022
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
Basic tests of the plot functionality in sherpa.astro.ui.

It is based on sherpa/ui/tests/test_ui_plot.py, but focusses on the
new routines, and "astro-specific" (currently PHA only) capabilities
that are in the astro layer.

There is very-limited checking that the plots are displaying the
correct data; it is more a check that the routines can be called.

"""

import copy
import logging
import numpy as np

import pytest

from sherpa.astro import ui
import sherpa.plot

from sherpa.astro.data import DataPHA
from sherpa.astro.instrument import create_arf
from sherpa.astro.plot import ARFPlot, BkgDataPlot, FluxHistogram, ModelHistogram, \
    OrderPlot, SourcePlot, BkgSourcePlot, ComponentModelPlot, ComponentSourcePlot, \
    ModelPHAHistogram, BkgModelHistogram, BkgModelPHAHistogram, \
    DataPHAPlot
from sherpa.data import Data1D, Data1DInt
from sherpa.models import basic
from sherpa.models.template import create_template_model
from sherpa.plot import FitPlot, PlotErr

from sherpa.utils.err import DataErr, IdentifierErr, ModelErr
from sherpa.utils.testing import requires_data, requires_fits, \
    requires_plotting, requires_pylab, requires_xspec

import sherpa.ui.utils
import sherpa.astro.ui.utils

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
_energies = np.asarray([0.5, 0.65, 0.75, 0.8, 0.9, 1., 1.1, 1.12, 1.3, 1.4, 1.5])
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


@pytest.mark.parametrize("idval", [None, 'x'])
def test_get_arf_plot_no_arf(idval, clean_astro_ui):
    """Errors out if we are unthoughtful enough to have no ARF
    """

    setup_example(idval)
    ui.get_data(idval).set_arf(None)

    with pytest.raises(DataErr) as exc:
        if idval is None:
            idval = 1
            ap = ui.get_arf_plot()
        else:
            ap = ui.get_arf_plot(idval)

    emsg = "data set '{}' does not have an associated ARF".format(idval)
    assert str(exc.value) == emsg


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
def test_get_arf_plot_recalc(idval, clean_astro_ui):
    """Check recalc=False handling
    """

    setup_example(idval)
    if idval is None:
        ap = ui.get_arf_plot(recalc=False)
    else:
        ap = ui.get_arf_plot(idval, recalc=False)

    assert isinstance(ap, ARFPlot)
    assert ap.xlo is None
    assert ap.y is None
    assert ap.title is None


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_order_plot(idval, clean_astro_ui):
    """Basic testing of get_order_plot: orders=None
    """

    setup_example(idval)
    if idval is None:
        op = ui.get_order_plot()
    else:
        op = ui.get_order_plot(idval)

    assert isinstance(op, OrderPlot)

    # Why is this using analysis=channel?
    # Would this be changed by #884?
    #
    xaxis = _data_chan.reshape((1, _data_chan.size))
    assert op.xlo == pytest.approx(xaxis)
    assert op.xhi == pytest.approx(xaxis + 1)

    yexp = _arf * 1.02e2 * (_energies[1:] - _energies[:-1])
    yexp.resize((1, yexp.size))
    assert op.y == pytest.approx(yexp)

    assert op.title == 'Model Orders [1]'
    assert op.xlabel == 'Channel'
    assert op.ylabel == 'Counts/sec/channel'


@pytest.mark.parametrize("idval", [None, 1, "one", 23])
def test_get_order_plot_recalc(idval, clean_astro_ui):
    """Check recalc=False handling
    """

    setup_example(idval)
    if idval is None:
        op = ui.get_order_plot(recalc=False)
    else:
        op = ui.get_order_plot(idval, recalc=False)

    assert isinstance(op, OrderPlot)
    assert op.xlo is None
    assert op.y is None
    assert op.title == 'Model'


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

    assert bp.xlo == pytest.approx(_data_chan)

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

    assert bp.xlo == pytest.approx(_energies_lo)
    assert bp.x == pytest.approx(_energies_mid)

    # normalise by exposure time and bin width
    #
    yexp = _data_bkg / (1201.0 * _bexpscale) / _energies_width
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'example-bkg'
    assert bp.xlabel == 'Energy (keV)'
    assert bp.ylabel == 'Counts/sec/keV'


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
    issue #879 and #880 has been resolved.

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
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Background Model Contribution'
    assert bp.xlabel == 'Channel'
    assert bp.ylabel == 'Counts/sec/channel'


def test_get_bkg_plot_recalc(clean_astro_ui):

    setup_example_bkg_model(None)

    bp = ui.get_bkg_plot(recalc=False)
    assert bp.xlo is None
    assert bp.y is None
    assert bp.title is None


@pytest.mark.parametrize("func", [ui.get_bkg_resid_plot,
                                  ui.get_bkg_ratio_plot,
                                  ui.get_bkg_delchi_plot,
                                  ui.get_bkg_chisqr_plot])
def test_get_bkg_xxx_plot_recalc(func, clean_astro_ui):

    setup_example_bkg_model(None)

    bp = func(recalc=False)
    assert bp.x is None
    assert bp.y is None
    assert bp.title == 'Model'


def test_get_bkg_model_plot_recalc(clean_astro_ui):

    setup_example_bkg_model(None)

    bp = ui.get_bkg_model_plot(recalc=False)
    assert bp.xlo is None
    assert bp.y is None
    assert bp.title == 'Background Model Contribution'


def test_get_bkg_source_plot_recalc(clean_astro_ui):

    setup_example_bkg_model(None)

    bp = ui.get_bkg_source_plot(recalc=False)
    assert bp.xlo is None
    assert bp.y is None
    assert bp.title == 'Source'


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

    assert bp.xlo == pytest.approx(_energies_lo)
    assert bp.xhi == pytest.approx(_energies_hi)

    yexp = _arf * BGND_NORM
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Background Model Contribution'
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
    yexp = _data_bkg / (1201.0 * _bexpscale) - _arf * BGND_NORM * _energies_width
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

    assert bp.x == pytest.approx(_energies_mid)

    # correct the counts by the bin width and exposure time
    #
    yexp = _data_bkg / (1201.0 * _bexpscale * _energies_width) - _arf * BGND_NORM
    assert bp.y == pytest.approx(yexp)

    assert bp.title == 'Residuals of example-bkg - Bkg Model'
    assert bp.xlabel == 'Energy (keV)'
    assert bp.ylabel == 'Counts/sec/keV'


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
    assert isinstance(dp, BkgDataPlot)
    assert isinstance(mp, BkgModelPHAHistogram)

    assert dp.title == 'example-bkg'
    assert mp.title == 'Background Model Contribution'

    for plot in [dp, mp]:
        assert plot.xlabel == 'Channel'
        assert plot.ylabel == 'Counts/sec/channel'

    assert dp.xlo == pytest.approx(_data_chan)
    assert mp.xlo == pytest.approx(_data_chan)

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dp.y == pytest.approx(dp.y)

    yexp = _arf * BGND_NORM * _energies_width
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
    assert isinstance(dp, BkgDataPlot)
    assert isinstance(mp, BkgModelPHAHistogram)

    assert dp.title == 'example-bkg'
    assert mp.title == 'Background Model Contribution'

    for plot in [dp, mp]:
        assert plot.xlabel == 'Energy (keV)'
        assert plot.ylabel == 'Counts/sec/keV'

    assert dp.xlo == pytest.approx(_energies_lo)
    assert mp.xlo == pytest.approx(_energies_lo)

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dp.y == pytest.approx(dp.y)

    yexp = _arf * BGND_NORM
    assert mp.y == pytest.approx(yexp)


def check_bkg_fit(plotfunc, isfit=True):
    """Is the background fit displayed?

    This only checks the plot object, not the plot "hardcopy" output
    (e.g. the pixel display/PNG output).
    """

    dplot = ui._session._bkgdataplot
    mplot = ui._session._bkgmodelplot
    assert isinstance(dplot, BkgDataPlot)
    assert isinstance(mplot, BkgModelPHAHistogram)

    # check the "source" plots are not set
    for plot in [ui._session._dataplot, ui._session._modelplot]:
        assert plot.x is None
        assert plot.y is None

    xlabel = 'Channel' if plotfunc == ui.plot_bkg_fit else ''

    # check plot basics
    for plot in [dplot, mplot]:
        assert plot.xlabel == xlabel
        assert plot.ylabel == 'Counts/sec/channel'

    assert dplot.xlo == pytest.approx(_data_chan)
    assert mplot.xlo == pytest.approx(_data_chan)

    assert dplot.title == 'example-bkg'
    assert mplot.title == 'Background Model Contribution'

    yexp = _data_bkg / (1201.0 * _bexpscale)
    assert dplot.y == pytest.approx(yexp)

    yexp = _arf * BGND_NORM * _energies_width
    assert mplot.y == pytest.approx(yexp)


def check_bkg_resid(plotfunc, isfit=True):
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
    for pfs, pd in [([ui.plot_bkg_delchi, ui.plot_bkg_fit_delchi], ui._session._bkgdelchiplot),
                    ([ui.plot_bkg_ratio, ui.plot_bkg_fit_ratio], ui._session._bkgratioplot),
                    ([ui.plot_bkg_resid, ui.plot_bkg_fit_resid], ui._session._bkgresidplot)]:
        if plotfunc in pfs:
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

    if plotfunc in [ui.plot_bkg_delchi, ui.plot_bkg_fit_delchi]:
        assert plot.ylabel == 'Sigma'
    elif plotfunc in [ui.plot_bkg_ratio, ui.plot_bkg_fit_ratio]:
        assert plot.ylabel == 'Data / Model'
    elif plotfunc in [ui.plot_bkg_resid, ui.plot_bkg_fit_resid]:
        assert plot.ylabel == 'Counts/sec/channel'
    else:
        assert False  # check I've caught everything

    # Plot title depends on whether this is part of a fit or not.
    #
    if isfit:
        assert plot.title == ''
    elif plotfunc == ui.plot_bkg_delchi:
        assert plot.title == 'Sigma Residuals for example-bkg'
    elif plotfunc == ui.plot_bkg_ratio:
        assert plot.title == 'Ratio of example-bkg : Bkg Model'
    elif plotfunc == ui.plot_bkg_resid:
        assert plot.title == 'Residuals of example-bkg - Bkg Model'
    else:
        assert False  # check I've caught everything

    assert plot.x == pytest.approx(_data_chan)
    assert plot.y is not None

    # the way the data and model are constructed, all residual values
    # should be negative, and the ratio values positive (or zero).
    #
    if plotfunc in [ui.plot_bkg_ratio, ui.plot_bkg_fit_ratio]:
        assert np.all(plot.y >= 0)
    else:
        assert np.all(plot.y < 0)


def check_bkg_chisqr(plotfunc, isfit=True):
    """Is the background residual displayed?"""

    # check the "other" background plots are not set
    for pd in [ui._session._bkgdelchiplot,
               ui._session._bkgratioplot,
               ui._session._bkgresidplot]:
        assert pd.x is None
        assert pd.y is None

    plot = ui._session._bkgchisqrplot
    assert plot.x is not None
    assert plot.y is not None

    # Very limited checks. The y-axis label depends on the
    # LaTeX emulation of the backend so will need changing
    # once we have multiple backends.
    #
    assert plot.xlabel == 'Channel'
    assert plot.ylabel == '$\\chi^2$'
    assert plot.title == '$\\chi^2$ for example-bkg'

    assert np.all(plot.y >= 0)


def check_bkg_model(plotfunc, isfit=True):
    """Is the background model displayed?"""

    # check the "other" background plots are not set
    for pd in [ui._session._bkgdelchiplot,
               ui._session._bkgratioplot,
               ui._session._bkgresidplot,
               ui._session._bkgchisqrplot]:
        assert pd.x is None
        assert pd.y is None

    assert ui._session._bkgsourceplot.xlo is None
    assert ui._session._bkgsourceplot.xhi is None
    assert ui._session._bkgsourceplot.y is None

    plot = ui._session._bkgmodelhisto
    assert isinstance(plot, BkgModelHistogram)
    assert plot.xlo is not None
    assert plot.xhi is not None
    assert plot.y is not None

    assert plot.xlabel == 'Channel'
    assert plot.ylabel == 'Counts/sec/channel'
    assert plot.title == 'Background Model Contribution'

    assert np.all(plot.y >= 0)


def check_bkg_source(plotfunc, isfit=True):
    """Is the background source model displayed?"""

    # check the "other" background plots are not set
    for pd in [ui._session._bkgdelchiplot,
               ui._session._bkgratioplot,
               ui._session._bkgresidplot,
               ui._session._bkgchisqrplot]:
        assert pd.x is None
        assert pd.y is None

    # check the background model is not set
    assert ui._session._bkgmodelhisto.xlo is None
    assert ui._session._bkgmodelhisto.xhi is None
    assert ui._session._bkgmodelhisto.y is None

    plot = ui._session._bkgsourceplot
    assert isinstance(plot, BkgSourcePlot)
    assert plot.xlo is not None
    assert plot.xhi is not None
    assert plot.y is not None

    # Very limited checks. The y-axis label depends on the
    # LaTeX emulation of the backend so will need changing
    # once we have multiple backends.
    #
    assert plot.xlabel == 'Energy (keV)'
    assert plot.ylabel == 'f(E)  Photons/sec/cm$^2$/keV '
    assert plot.title == 'Source Model of example-bkg'

    assert np.all(plot.y >= 0)


@requires_plotting
@pytest.mark.parametrize("idval", [None, 1, "one", 23])
@pytest.mark.parametrize("plotfunc,checkfuncs",
                         [(ui.plot_bkg_delchi, [check_bkg_resid]),
                          (ui.plot_bkg_ratio, [check_bkg_resid]),
                          (ui.plot_bkg_resid, [check_bkg_resid]),
                          (ui.plot_bkg_chisqr, [check_bkg_chisqr]),
                          (ui.plot_bkg_model, [check_bkg_model]),
                          (ui.plot_bkg_source, [check_bkg_source]),
                          (ui.plot_bkg_fit, [check_bkg_fit]),
                          (ui.plot_bkg_fit_delchi, [check_bkg_fit, check_bkg_resid]),
                          (ui.plot_bkg_fit_ratio, [check_bkg_fit, check_bkg_resid]),
                          (ui.plot_bkg_fit_resid, [check_bkg_fit, check_bkg_resid])])
def test_bkg_plot_xxx(idval, plotfunc, checkfuncs, clean_astro_ui):
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
    isfit = checkfuncs[0] == check_bkg_fit
    for checkfunc in checkfuncs:
        checkfunc(plotfunc, isfit=isfit)


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
def basic_pha1_bg(make_data_path):
    """Create a basic PHA-1 data set/setup for the background.

    This is currently only used to get a dataset with no response.
    """

    ui.set_default_id('tst')
    ui.load_pha(make_data_path('3c273_bg.pi'))
    ui.notice(33, 676)
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
                    ui.plot_fit_ratio,
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
@pytest.mark.parametrize("plotfunc",
                         [ui.plot_data,
                          ui.plot_model,
                          ui.plot_source,
                          ui.plot_resid,
                          ui.plot_delchi,
                          ui.plot_ratio,
                          ui.plot_fit,
                          ui.plot_fit_delchi,
                          ui.plot_fit_resid,
                          ui.plot_fit_ratio,
                          ui.plot_chisqr])
def test_xxx_plot_clearwindow(hide_logging, clean_astro_ui, plotfunc):
    """If set clearwindow=False does it run sensibly for Data1D data?

    This does not have a great check for what "sensibly" means.
    """

    from matplotlib import pyplot as plt

    ui.load_arrays(1, [1, 5, 10], [4, 2, 6])
    ui.set_source(ui.const1d.mdl)

    plotfunc()
    axes = plt.gcf().axes
    nlines0 = [len(ax.lines) for ax in axes]

    plotfunc(clearwindow=False)
    nlines1 = [len(ax.lines) for ax in axes]

    # Check that there are twice the number of lines. This
    # should be sufficient for all the plot types to check
    # that new data has been added.
    #
    assert nlines1 == [2 * n for n in nlines0]


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", _basic_plotfuncs)
def test_pha1_plot_clearwindow(hide_logging, clean_astro_ui, basic_pha1, plotfunc):
    """If set clearwindow=False does it run sensibly for PHA data?

    This does not have a great check for what "sensibly" means.
    """

    from matplotlib import pyplot as plt

    plotfunc()
    axes = plt.gcf().axes
    nlines0 = [len(ax.lines) for ax in axes]

    plotfunc(clearwindow=False)
    nlines1 = [len(ax.lines) for ax in axes]

    # Check that there are twice the number of lines. This
    # should be sufficient for all the plot types to check
    # that new data has been added.
    #
    assert nlines1 == [2 * n for n in nlines0]


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", [ui.plot_bkg_model,
                                      ui.plot_bkg_source,
                                      ui.plot_bkg_fit,
                                      ui.plot_bkg_fit_resid,
                                      ui.plot_bkg_fit_ratio,
                                      ui.plot_bkg_fit_delchi])
def test_pha1_bkg_plot(plotfunc, clean_astro_ui, basic_pha1, hide_logging):
    """Test issue #943 and general check of plot_bkg_xxx"""

    ui.unsubtract()
    ui.set_bkg_source(ui.const1d.bmdl + ui.gauss1d.bgmdl)
    plotfunc()


@requires_plotting
@requires_fits
@requires_data
@pytest.mark.parametrize("plotfunc", [ui.plot_model_component,
                                      ui.plot_source_component])
def test_pha1_plot_component(clean_astro_ui, basic_pha1, plotfunc):
    plotfunc("pl")


def test_data1_get_model_component_plot(clean_astro_ui):
    """There is no response concept for non-PHA data.

    See test_pha1_get_model_component_plot_add_response
    for the PHA explanation.
    """

    x = np.asarray([1, 2, 3])
    ui.load_arrays(1, x, x)
    ui.set_source(ui.polynom1d.mdl)
    mdl.c0 = -3
    mdl.c1 = 2

    # Apparently we have sherpa.astro.plot.ComponentModelPlot
    # and sherpa.plot.ComponentModelPlot.
    #
    mplot = ui.get_model_component_plot('mdl')
    assert isinstance(mplot, sherpa.plot.ComponentModelPlot)
    assert mplot.title == 'Model component: polynom1d.mdl'

    # Unlike the PHA case it's easy to evaluate the model
    # for this test.
    assert mplot.x == pytest.approx([1, 2, 3])
    assert mplot.y == pytest.approx([-1, 1, 3])


def check_pha1_model_component_plot(mplot, mdl):
    """The checks for get_model_component_plot with PHA data

    This is for when the response is included.
    """

    assert isinstance(mplot, ComponentModelPlot)

    # The plot title includes the exposure time whose output could change
    # so do not use an equality check but something a bit-more forgiving
    # (could use a regexp but not worth it).
    #
    assert mplot.title.startswith('Model component: apply_rmf(apply_arf((38564.60')
    assert mplot.title.endswith(' * powlaw1d.pl)))')

    # This may change if we change the filtering/grouping code
    # Note that the model is evaluated on the un-grouped data.
    #
    assert len(mplot.xlo) == 644
    assert mplot.xlo[0] == pytest.approx(0.467200)
    assert mplot.xlo[-1] == pytest.approx(9.85500)
    assert mplot.xhi[0] == pytest.approx(0.481800)
    assert mplot.xhi[-1] == pytest.approx(9.869600)

    # The model is evaluated over the full channel range so we need to
    # restrict to the noticed channel range (of 33 to 676 inclusive).
    # I am 'cheating' with the argument to mdl() here since it is currently
    # ignored as long as it's specified.
    #
    yexp = mdl([1])
    assert len(yexp) == 1024
    yexp = yexp[32:676]

    dy = mplot.xhi - mplot.xlo
    texp = 38564.608926889
    assert mplot.y == pytest.approx(yexp / dy / texp)

    # manual check of one the y values
    #
    assert mplot.y[0] == pytest.approx(0.00501761)


# NOTE: the following model_component tests are just different
# enough in how you send in the data that it's not easy to parametrize
# them, so we have a bunch of similar routines rather than a single
# parametrized test.

@requires_fits
@requires_data
def test_pha1_get_model_component_plot_add_response(clean_astro_ui, basic_pha1):
    """Do we automatically add in the response? See issue #1020

    Note that this test does not need requires_plotting since
    we don't use anything that requires the plot backend.
    """

    # It's important to also test the "specify a string not
    # a model object" feature of the get_xxx_plot call here.
    # The explicit call in xxxx_with_response below does
    # not use a string.
    #
    mplot = ui.get_model_component_plot('pl')
    rsp = ui.get_response()
    mdl = rsp(pl)
    check_pha1_model_component_plot(mplot, mdl)


@requires_fits
@requires_data
def test_pha1_get_model_component_plot_with_response(clean_astro_ui, basic_pha1):
    """What happens if we explicitly include the response?

    Note that this test does not need requires_plotting since
    we don't use anything that requires the plot backend.
    """

    rsp = ui.get_response()
    mdl = rsp(pl)
    mplot = ui.get_model_component_plot(mdl)
    check_pha1_model_component_plot(mplot, mdl)


@requires_fits
@requires_data
def test_pha1_get_model_component_plot_no_response(clean_astro_ui, basic_pha1_bg):
    """What happens if we have no response?"""

    # safety check the file remains response "free".
    pha = ui.get_data('tst')
    assert pha.units == 'channel'
    assert pha.get_response() == (None, None)
    with pytest.raises(DataErr) as exc:
        ui.get_response()

    assert str(exc.value).startswith('No instrument response found for dataset')
    assert str(exc.value).endswith('3c273_bg.pi')

    # use the background dataset as it's response "free"
    mplot = ui.get_model_component_plot('pl')
    print(mplot)

    assert isinstance(mplot, ComponentModelPlot)
    assert mplot.title == 'Model component: powlaw1d.pl'

    xlo = np.arange(33, 677)
    assert mplot.xlo == pytest.approx(xlo)
    assert mplot.xhi == pytest.approx(xlo + 1)

    # The y values are currently pl(channels) / bin width / time
    # for the selected channel range. Fortunately the bin width is 1.
    #
    chans = np.arange(33, 677)
    yexp = pl(chans)
    texp = 38564.608926889
    assert mplot.y == pytest.approx(yexp / texp)

    # manual check of one the y values
    #
    assert np.log10(mplot.y[0]) == pytest.approx(-11.276372)


def check_pha1_plot_model_component_plot():
    """Check the model component plot.

    This is much-more limited than the test_pha1_get_model_component_plot_xxx
    variants, as all we do is check the y axis range makes sense
    (as it wouldn't have with #1020 unfixed).

    Any test calling this requires @requires_plotting
    """
    from matplotlib import pyplot as plt

    # Just check the Y range to be close to the expected range - this
    # could vary a bit across matplotlib versions
    #
    # The range is expected to be ~1e-5 to 1e-2
    ylim = plt.ylim()
    assert ylim[0] > 1e-6
    assert ylim[1] > 1e-2
    assert ylim[1] < 2e-2


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_model_component_add_response(clean_astro_ui, basic_pha1):
    """Do we automatically add in the response? See issue #1020.
    """
    ui.plot_model_component('pl', ylog=True)
    check_pha1_plot_model_component_plot()


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_with_model_component_add_response(clean_astro_ui, basic_pha1):
    """Do we automatically add in the response? See issue #1020.
    """
    ui.plot('model_component', 'pl', ylog=True)
    check_pha1_plot_model_component_plot()


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_model_component_with_response(clean_astro_ui, basic_pha1):
    """Plot is okay if we include the response
    """
    rsp = ui.get_response()
    ui.plot_model_component(rsp(pl), ylog=True)
    check_pha1_plot_model_component_plot()


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_with_model_component_with_response(clean_astro_ui, basic_pha1):
    """Plot is okay if we include the response
    """
    rsp = ui.get_response()
    ui.plot('model_component', rsp(pl), ylog=True)
    check_pha1_plot_model_component_plot()


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_model_component_no_response(clean_astro_ui, basic_pha1_bg):
    """PHA file but no response
    """
    from matplotlib import pyplot as plt

    ui.plot_model_component('pl', ylog=True)

    # The range is expected to be ~1e-14 to ~7e-12 which is very different
    # from the 'with response' values
    ylim = plt.ylim()
    assert ylim[1] < 1e-11


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_with_model_component_no_response(clean_astro_ui, basic_pha1_bg):
    """PHA file but no response
    """
    from matplotlib import pyplot as plt

    ui.plot('model_component', pl, ylog=True)

    # The range is expected to be ~1e-14 to ~7e-12 which is very different
    # from the 'with response' values
    ylim = plt.ylim()
    assert ylim[1] < 1e-11


@requires_plotting
@requires_fits
@requires_data
def test_pha1_plot_order(clean_astro_ui, basic_pha1):
    ui.plot_order()


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
                  ui.contour_source]


@requires_plotting
@requires_fits
@requires_data
def test_img_contour_function(clean_astro_ui, basic_img):
    # can we call contour; do not try to be exhaustive
    ui.contour("data", "model", "source", "fit")


@requires_pylab
@requires_fits
@requires_data
def test_img_contour_function_kwarg(clean_astro_ui, basic_img):
    """Check we can change the alpha setting."""

    from matplotlib import pyplot as plt

    ui.contour("data", "model", "source", "fit", alpha=0.2)

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 4

    for i, ax in enumerate(axes, 1):

        w = i - 1
        assert ax.get_subplotspec().get_geometry() == (2, 2, w, w)

        assert ax.get_xscale() == 'linear'
        assert ax.get_yscale() == 'linear'

        assert len(ax.lines) == 0
        assert len(ax.collections) > 0

        col = ax.collections[0]
        assert col.get_alpha() == 0.2

    plt.close()


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
def test_pha1_plot_data_options(caplog, clean_astro_ui, basic_pha1):
    """Test that the options have changed things, where easy to do so"""

    from matplotlib import pyplot as plt
    import matplotlib

    prefs = ui.get_data_plot_prefs('tst')

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

    # linecolor is unused, but check we can set it as existing code
    # may reference it.
    prefs['linecolor'] = 'brown'
    prefs['linestyle'] = '-.'

    prefs['marker'] = 's'
    prefs['markerfacecolor'] = 'cyan'
    prefs['markersize'] = 10

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.plot_data()

    # check for linecolor warning
    assert len(caplog.record_tuples) == 1
    rec = caplog.record_tuples[0]
    assert len(rec) == 3
    loc, lvl, msg = rec

    assert loc == 'sherpa.plot.pylab_backend'
    assert lvl == logging.WARNING
    assert msg == 'The linecolor attribute (brown) is unused.'

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

    assert len(ax.lines) == 2
    line = ax.lines[0]

    assert line.get_color() == 'orange'
    assert line.get_linestyle() == '-.'
    assert line.get_marker() == 'None'
    assert line.get_markerfacecolor() == 'orange'
    assert line.get_markersize() == pytest.approx(6.0)

    line = ax.lines[1]

    assert line.get_color() == 'orange'
    assert line.get_linestyle() == 'None'
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
def test_pha1_plot_model_options(caplog, clean_astro_ui, basic_pha1):
    """Test that the options have changed things, where easy to do so

    In matplotlib 3.1 the plot_model call causes a MatplotlibDeprecationWarning
    to be created:

    Passing the drawstyle with the linestyle as a single string is deprecated since Matplotlib 3.1 and support will be removed in 3.3; please pass the drawstyle separately using the drawstyle keyword argument to Line2D or set_drawstyle() method (or ds/set_ds()).

    This warning is hidden by the test suite (sherpa/conftest.py) so that
    it doesn't cause the tests to fail. Note that a number of other tests
    in this module also cause this warning to be displayed.

    """

    from matplotlib import pyplot as plt

    prefs = ui.get_model_plot_prefs('tst')

    # check the preference are as expected for the boolean cases
    assert not prefs['xlog']
    assert not prefs['ylog']

    # Only change the X axis here
    prefs['xlog'] = True

    prefs['color'] = 'green'

    # linecolor is unused, but check we can set it as existing code
    # may reference it.
    prefs['linecolor'] = 'red'
    prefs['linestyle'] = 'dashed'

    prefs['marker'] = '*'
    prefs['markerfacecolor'] = 'yellow'
    prefs['markersize'] = 8

    with caplog.at_level(logging.INFO, logger='sherpa'):
        ui.plot_model()

    # check for linecolor warning
    assert len(caplog.record_tuples) == 1
    rec = caplog.record_tuples[0]
    assert len(rec) == 3
    loc, lvl, msg = rec

    assert loc == 'sherpa.plot.pylab_backend'
    assert lvl == logging.WARNING
    assert msg == 'The linecolor attribute (red) is unused.'

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
    assert xmin == pytest.approx(0.4011095556520917)
    assert xmax == pytest.approx(11.495805328091725)

    ymin, ymax = ax.get_ylim()
    assert ymin == pytest.approx(-0.00045772936258082011, rel=0.01)
    assert ymax == pytest.approx(0.009940286575890335, rel=0.01)

    assert len(ax.lines) == 2

    # Apparently color wins out over linecolor
    line = ax.lines[0]
    assert line.get_color() == 'green'
    assert line.get_linestyle() == '--'  # note: input was dashed
    assert line.get_marker() == 'None'
    assert line.get_markerfacecolor() == 'green'
    assert line.get_markersize() == pytest.approx(6.0)

    pts = ax.lines[1]
    assert pts.get_color() == 'green'
    assert pts.get_linestyle() == 'None'
    assert pts.get_marker() == '*'
    assert pts.get_markerfacecolor() == 'yellow'
    assert pts.get_markersize() == pytest.approx(8.0)

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
    # TODO: why is the y range different now> have the errors changed
    assert ymin == pytest.approx(7.644069935298475e-05)
    assert ymax == pytest.approx(0.017031102671151491)

    assert len(ax.lines) == 4

    # DATA
    #
    line = ax.lines[0]

    assert line.get_color() == 'orange'
    assert line.get_linestyle() == '-.'
    assert line.get_marker() == 'None'
    assert line.get_markerfacecolor() == 'orange'
    assert line.get_markersize() == pytest.approx(6.0)
    assert line.get_alpha() == pytest.approx(0.7)

    # MODEL - line
    #
    line = ax.lines[1]
    assert line.get_color() != 'green'  # option over-ridden
    assert line.get_linestyle() == 'None'  # option over-ridden
    assert line.get_marker() == 's'
    assert line.get_markerfacecolor() == 'cyan'  # line.get_color()  # option over-ridden
    assert line.get_markersize() == pytest.approx(10.0)
    assert line.get_alpha() == pytest.approx(0.7)

    # MODEL - points
    #
    pts = ax.lines[2]
    assert pts.get_color() != 'green'  # option over-ridden
    assert pts.get_linestyle() == '-'
    assert pts.get_marker() == 'None'
    assert pts.get_markerfacecolor() == '#1f77b4'  # line.get_color()  # option over-ridden
    assert pts.get_markersize() == pytest.approx(6.0)
    assert pts.get_alpha() == pytest.approx(0.7)

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


@requires_fits
@requires_data
@pytest.mark.parametrize("units,xlabel,ylabel,xlo,xhi",
                         [("channel", 'Channel', 'Counts/sec/channel',
                           33, 677),
                          ("wavelength", 'Wavelength (Angstrom)', 'Counts/sec/Angstrom',
                           26.537710718511885,  1.2562229845315145)])
def test_bug920(units, xlabel, ylabel, xlo, xhi, clean_astro_ui, basic_pha1):
    """plot_model units appear to depend on setting.

    We check a few things, but it's the x axis value that is #920.
    """

    assert ui.get_analysis() == 'energy'
    mplot1 = copy.deepcopy(ui.get_model_plot())

    ui.set_analysis(units)
    assert ui.get_analysis() == units

    # You need to create the model plot to trigger the bug.
    mplot2 = copy.deepcopy(ui.get_model_plot())

    ui.set_analysis('energy')
    assert ui.get_analysis() == 'energy'

    mplot3 = copy.deepcopy(ui.get_model_plot())

    assert mplot1.xlabel == 'Energy (keV)'
    assert mplot2.xlabel == xlabel
    assert mplot3.xlabel == 'Energy (keV)'

    assert mplot1.ylabel == 'Counts/sec/keV'
    assert mplot2.ylabel == ylabel
    assert mplot3.ylabel == 'Counts/sec/keV'

    # mplot3 should be the same as mplot1
    assert mplot1.xlo[0] == pytest.approx(0.46720001101493835)
    assert mplot1.xhi[-1] == pytest.approx(9.869600296020508)

    assert mplot2.xlo[0] == pytest.approx(xlo)
    assert mplot2.xhi[-1] == pytest.approx(xhi)

    assert mplot1.xlo.size == mplot1.y.size
    assert mplot2.xlo.size == mplot2.y.size
    assert mplot3.xlo.size == mplot3.y.size

    # The number of bins are the same, it's just the numbers are different
    assert mplot1.xlo.size == 644
    assert mplot2.xlo.size == 644
    assert mplot3.xlo.size == 644

    # This should be equality, but allow small differences
    assert mplot3.xlo == pytest.approx(mplot1.xlo)
    assert mplot3.y == pytest.approx(mplot1.y)


def validate_flux_histogram(fhist, energy):
    """Limited checks for test_pha1_plot_foo_flux/test_pha1_get_foo_flux_hist

    histtype is 'Photon' or 'Energy'
    """

    assert fhist is not None
    assert isinstance(fhist, FluxHistogram)

    # Very minimal checks.
    #
    assert fhist.flux.shape == (200,)
    assert fhist.modelvals.shape == (200, 3)
    assert fhist.xlo.shape == (21,)
    assert fhist.xhi.shape == (21,)
    assert fhist.y.shape == (21,)

    # Do the fluxes and the histogram agree?
    #
    fmin = fhist.flux.min()
    fmax = fhist.flux.max()
    assert fmin == pytest.approx(fhist.xlo[0])
    assert fmax == pytest.approx(fhist.xhi[-1])

    # Check labels. The label depends on the plot backend through the
    # use of LaTeX.
    #
    if energy:
        assert fhist.xlabel.startswith('Energy flux (ergs cm')
    else:
        assert fhist.xlabel.startswith('Photon flux (Photons cm')

    assert fhist.xlabel.find(' sec') > 20
    assert fhist.xlabel.endswith(')')

    assert fhist.ylabel == 'Frequency'
    if energy:
        assert fhist.title == 'Energy flux distribution'
    else:
        assert fhist.title == 'Photon flux distribution'


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy,plotfunc,getfunc",
                         [(True, ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (False, ui.plot_photon_flux, ui.get_photon_flux_hist)])
@pytest.mark.parametrize("correlated", [False, True])
def test_pha1_plot_foo_flux(energy, plotfunc, getfunc, correlated, clean_astro_ui, basic_pha1):
    """Can we call plot_energy/photon_flux and then the get_ func (recalc=False)

    We extend the basic_pha1 test by including an XSPEC
    absorption model.
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    # At this point the return should be None
    res = getfunc(recalc=False)
    assert isinstance(res, FluxHistogram)
    assert res.flux is None
    assert res.xlo is None
    assert res.y is None

    # Since the results are not being inspected here, the "quality"
    # of the results isn't important, so we can use a relatively-low
    # number of iterations.
    #
    plotfunc(lo=0.5, hi=2, num=200, bins=20, correlated=correlated)

    # check we can access these results (relying on the fact that the num
    # and bins arguments have been changed from their default values).
    #
    res = getfunc(recalc=False)
    validate_flux_histogram(res, energy)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy,plotfunc,getfunc",
                         [(True, ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (False, ui.plot_photon_flux, ui.get_photon_flux_hist)])
def test_pha1_plot_foo_flux_recalc(energy, plotfunc, getfunc, clean_astro_ui, basic_pha1):
    """Just check we can call recalc on the routine

    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    # Since the results are not being inspected here, the "quality"
    # of the results isn't important, so we can use a relatively-low
    # number of iterations.
    #
    plotfunc(lo=0.5, hi=2, num=200, bins=20, correlated=True)
    plotfunc(lo=0.5, hi=2, num=200, bins=20, correlated=True, recalc=False)

    # check we can access these results (relying on the fact that the num
    # and bins arguments have been changed from their default values).
    #
    res = getfunc(recalc=False)
    validate_flux_histogram(res, energy)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy,getfunc", [(True, ui.get_energy_flux_hist),
                                            (False, ui.get_photon_flux_hist)])
@pytest.mark.parametrize("correlated", [False, True])
def test_pha1_get_foo_flux_hist(energy, getfunc, correlated, clean_astro_ui, basic_pha1):
    """Can we call get_energy/photon_flux_hist?

    See test_pha1_plot_foo_flux.
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    # Since the results are not being inspected here, the "quality"
    # of the results isn't important, so we can use a relatively-low
    # number of iterations.
    #
    res = getfunc(lo=0.5, hi=2, num=200, bins=20, correlated=correlated)
    validate_flux_histogram(res, energy)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("getfunc", [ui.get_energy_flux_hist,
                                     ui.get_photon_flux_hist])
def test_pha1_get_foo_flux_hist_no_data(getfunc, clean_astro_ui, basic_pha1):
    """What happens when there's no data?

    This is a regression test to check if this behavior gets
    changed: for example maybe an error gets raised.
    """

    empty = getfunc(recalc=False)
    assert empty is not None
    assert isinstance(empty, FluxHistogram)

    for field in ["modelvals", "flux", "xlo", "xhi", "y"]:
        assert getattr(empty, field) is None


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("getfunc", [ui.get_energy_flux_hist,
                                     ui.get_photon_flux_hist])
@pytest.mark.parametrize("scale", [None, 1.0, 2.0])
def test_pha1_get_foo_flux_hist_scales(getfunc, scale,
                                       clean_astro_ui, basic_pha1,
                                       hide_logging, reset_seed):
    """Can we call get_energy/photon_flux_hist with the scales argument.

    Run with the covar-calculated errors and then manually-specified
    errors and check that the parameter distributions change in the
    expected manner. The test is run for the covariance errors
    multiplied by scale. If scale is None then scales is not
    given.

    This is issue #801.
    """

    np.random.seed(42873)

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()
    ui.covar()
    covmat = ui.get_covar_results().extra_output
    errs = np.sqrt(covmat.diagonal())

    if scale is None:
        scales = None
    else:
        errs = scale * errs
        scales = errs

    # Run enough times that we can get a reasonable distribution
    res = getfunc(lo=0.5, hi=2, num=1000, bins=20, correlated=False,
                  scales=scales)

    pvals = res.modelvals

    if scale is None:
        scale = 1.0

    # rely on the fixed seed, but we still get a range of values
    # for the first parameter. Why is that?
    #
    std0 = np.std(pvals[:, 0]) / scale
    assert std0 > 0.032
    assert std0 < 0.040

    std1 = np.std(pvals[:, 1]) / scale
    assert std1 == pytest.approx(0.1328676086353561, rel=1e-3)

    std2 = np.log10(np.std(pvals[:, 2])) - np.log10(scale)
    assert std2 == pytest.approx(-4.520334184088533, rel=1e-3)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("plotfunc,getfunc",
                         [(ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (ui.plot_photon_flux, ui.get_photon_flux_hist)])
@pytest.mark.parametrize("scale", [1.0, 2.0])
def test_pha1_plot_foo_flux_scales(plotfunc, getfunc, scale,
                                   clean_astro_ui, basic_pha1,
                                   hide_logging, reset_seed):
    """Can we call plot_energy/photon_flux with the scales argument.

    Based on test_pha1_get_foo_flux_hist_scales. By using some
    non-standard arguments we can use the get routine (recalc=False)
    to check that the plot did use the scales.
    """

    np.random.seed(38259)

    # unlike test_pha1_get-foo_flux_hist_scales, only use a power-law
    # component
    #
    # Ensure near the minimum
    ui.fit()
    ui.covar()
    covmat = ui.get_covar_results().extra_output
    errs = scale * np.sqrt(covmat.diagonal())

    # Run enough times that we can get a reasonable distribution
    plotfunc(lo=0.5, hi=2, num=1000, bins=20, correlated=False,
             scales=errs)

    res = getfunc(recalc=False)
    pvals = res.modelvals

    assert res.flux.shape == (1000,)
    assert res.modelvals.shape == (1000, 2)
    assert res.xlo.shape == (21,)
    assert res.xhi.shape == (21,)
    assert res.y.shape == (21,)

    std0 = np.std(pvals[:, 0]) / scale
    std1 = np.log10(np.std(pvals[:, 1])) - np.log10(scale)

    assert std0 == pytest.approx(0.06927323647207109)
    assert std1 == pytest.approx(-4.951357667423949)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("plotfunc,getfunc,ratio",
                         [(ui.plot_energy_flux, ui.get_energy_flux_hist, 1.141794001904922),
                          (ui.plot_photon_flux, ui.get_photon_flux_hist, 1.1892431836751618)])
def test_pha1_plot_foo_flux_model(plotfunc, getfunc, ratio,
                                  clean_astro_ui, basic_pha1,
                                  hide_logging, reset_seed):
    """Can we call plot_energy/photon_flux with the model argument.

    Based on test_pha1_get_foo_flux_hist_scales. By using some
    non-standard arguments we can use the get routine (recalc=False)
    to check that the plot did use the scales.

    """

    np.random.seed(286728)

    # This time we want the absorbing component to make a difference
    # between the two plots.
    #
    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()
    ui.covar()
    covmat = ui.get_covar_results().extra_output
    errs = np.sqrt(covmat.diagonal())

    # Due to the way the get* routines work in the sherpa.astro.ui module,
    # the following will return the same object, so x1 and x2 will be
    # identical (and hence only contain values from the second plotfunc).
    #
    #   plotfunc()
    #   x1 = getfunc(recalc=False)
    #   plotfunc()
    #   x2 = getunc(recalc=False)
    #
    # This means that the tests below check the values from getfunc
    # before calling the new plotfunc.
    #
    # It is probably true that x1 would change contents as soon as
    # the second plotfunc() call is made above (unsure).
    #

    # Absorbed flux
    plotfunc(lo=0.5, hi=2, num=1000, bins=19, correlated=False)
    res = getfunc(recalc=False)

    avals = res.modelvals
    assert avals.shape == (1000, 3)

    std1 = np.std(avals[:, 1])
    std2 = np.log10(np.std(avals[:, 2]))
    assert std1 == pytest.approx(0.1330728846451271, rel=1e-3)
    assert std2 == pytest.approx(-4.54079387550295, rel=1e-3)

    assert res.clipped.shape == (1000,)

    assert res.y.shape == (20,)

    aflux = np.median(res.flux)

    # Unabsorbed flux
    plotfunc(lo=0.5, hi=2, model=orig_mdl, num=1000, bins=21,
             correlated=False)
    res = getfunc(recalc=False)

    uvals = res.modelvals
    assert uvals.shape == (1000, 3)

    std1 = np.std(uvals[:, 1])
    std2 = np.log10(np.std(uvals[:, 2]))
    assert std1 == pytest.approx(0.13648119989822335, rel=1e-3)
    assert std2 == pytest.approx(-4.550978251403581, rel=1e-3)

    assert res.y.shape == (22,)

    uflux = np.median(res.flux)

    # Is the unabsorbed to absorbed median flux close to the value
    # calculated from the best-fit solutions (specified as a number
    # rather than calculated here to act as a regression test).
    #
    got = uflux / aflux
    assert got == pytest.approx(ratio, rel=1e-3)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy,plotfunc,getfunc",
                         [(True, ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (False, ui.plot_photon_flux, ui.get_photon_flux_hist)])
def test_pha1_plot_foo_flux_soft(energy, plotfunc, getfunc, clean_astro_ui, basic_pha1):
    """Check we can send clip=soft
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    # Since the results are not being inspected here, the "quality"
    # of the results isn't important, so we can use a relatively-low
    # number of iterations.
    #
    plotfunc(lo=0.5, hi=2, num=200, bins=20, correlated=True, clip='soft')

    # check we can access these results (relying on the fact that the num
    # and bins arguments have been changed from their default values).
    #
    res = getfunc(recalc=False)
    validate_flux_histogram(res, energy)

    # check we have clip information and assume at least one bin is
    # clipped (expect ~ 40 of 200 from testing this)
    clip = res.clipped
    c0 = clip == 0
    c1 = clip == 1
    assert (c0 | c1).all()
    assert c0.any()
    assert c1.any()


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy, plotfunc,getfunc",
                         [(True, ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (False, ui.plot_photon_flux, ui.get_photon_flux_hist)])
def test_pha1_plot_foo_flux_none(energy, plotfunc, getfunc, clean_astro_ui, basic_pha1):
    """Check we can send clip=none

    Copy of test_pha1_plot_foo_flux_none
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    # Since the results are not being inspected here, the "quality"
    # of the results isn't important, so we can use a relatively-low
    # number of iterations.
    #
    plotfunc(lo=0.5, hi=2, num=200, bins=20, correlated=True, clip='none')

    # check we can access these results (relying on the fact that the num
    # and bins arguments have been changed from their default values).
    #
    res = getfunc(recalc=False)
    validate_flux_histogram(res, energy)

    # check we have clip information but that it's all zeros
    clip = res.clipped
    c0 = clip == 0
    c1 = clip == 1
    assert (c0 | c1).all()
    assert c0.all()
    assert not(c1.any())


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy,getfunc", [(True, ui.get_energy_flux_hist),
                                            (False, ui.get_photon_flux_hist)])
def test_pha1_get_foo_flux_soft(energy, getfunc, clean_astro_ui, basic_pha1):
    """Can we send clip=soft?
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    res = getfunc(lo=0.5, hi=2, num=200, bins=20, correlated=False,
                  clip='soft')
    validate_flux_histogram(res, energy)

    # check we have clip information and assume at least one bin is
    # clipped (expect ~ 40 of 200 from testing this)
    clip = res.clipped
    c0 = clip == 0
    c1 = clip == 1
    assert (c0 | c1).all()
    assert c0.any()
    assert c1.any()


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("energy, getfunc", [(True, ui.get_energy_flux_hist),
                                             (False, ui.get_photon_flux_hist)])
def test_pha1_get_foo_flux_none(energy, getfunc, clean_astro_ui, basic_pha1):
    """Can we send clip=none?

    Copy of test_pha1_get_foo_flux_soft
    """

    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()

    res = getfunc(lo=0.5, hi=2, num=200, bins=20, correlated=False,
                  clip='none')
    validate_flux_histogram(res, energy)

    # check we have clip information and all are 0
    clip = res.clipped
    c0 = clip == 0
    c1 = clip == 1
    assert (c0 | c1).all()
    assert c0.all()
    assert not(c1.any())


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("plotfunc,getfunc",
                         [(ui.plot_energy_flux, ui.get_energy_flux_hist),
                          (ui.plot_photon_flux, ui.get_photon_flux_hist)])
def test_pha1_plot_foo_flux_multi(plotfunc, getfunc,
                                  make_data_path, clean_astro_ui,
                                  hide_logging, reset_seed):
    """Can we call plot_energy/photon_flux with multiple datasets.
    """

    np.random.seed(7267239)

    ui.load_pha(1, make_data_path('obs1.pi'))
    ui.load_pha(3, make_data_path('obs1.pi'))

    ui.notice(0.5, 7)

    mdl = ui.xswabs.gal * ui.powlaw1d.pl
    gal = ui.get_model_component('gal')
    pl = ui.get_model_component('pl')
    gal.nh = 0.04
    gal.nh.freeze()
    pl.gamma = 1.93
    pl.ampl = 1.74e-4

    ui.set_source(1, mdl)
    ui.set_source(3, mdl)

    ui.set_stat('cstat')

    # Ensure near the minimum
    ui.fit()

    n = 200

    # Use all datasets since id=None
    plotfunc(lo=0.5, hi=7, num=n, bins=19, correlated=False)
    res = getfunc(recalc=False)
    assert res.y.shape == (20,)
    avals = res.modelvals.copy()

    # Use all datasets as explicit
    plotfunc(lo=0.5, hi=7, num=n, bins=19, correlated=False,
             id=1, otherids=(3,))
    res = getfunc(recalc=False)
    assert res.y.shape == (20,)
    bvals = res.modelvals.copy()

    # Use only dataset 1 (the shorter dataset)
    plotfunc(lo=0.5, hi=7, num=n, bins=19, correlated=False,
             id=1)
    res = getfunc(recalc=False)
    assert res.y.shape == (20,)
    cvals = res.modelvals.copy()

    assert avals.shape == (n, 2)
    assert bvals.shape == (n, 2)
    assert cvals.shape == (n, 2)

    # Let's just check the standard deviation of the gamma parameter,
    # which should be similar for avals and bvals, and larger for cvals.
    #
    s1 = np.std(avals[:, 0])
    s2 = np.std(bvals[:, 0])
    s3 = np.std(cvals[:, 0])

    assert s1 == pytest.approx(0.1774218722298197)
    assert s2 == pytest.approx(0.1688597911809269)
    assert s3 == pytest.approx(0.22703663059917598)


@requires_plotting
@requires_fits
@requires_data
@requires_xspec
@pytest.mark.parametrize("getfunc,ratio",
                         [(ui.get_energy_flux_hist, 1.149863517748443),
                          (ui.get_photon_flux_hist, 1.1880213994341908)])
def test_pha1_get_foo_flux_hist_model(getfunc, ratio,
                                      clean_astro_ui, basic_pha1,
                                      hide_logging, reset_seed):
    """Can we call get_energy/photon_flux_hist with the model argument.

    Very similar to test_pha1_plot_foo_flux_model.
    """

    np.random.seed(2731)

    # This time we want the absorbing component to make a difference
    # between the two plots.
    #
    orig_mdl = ui.get_source('tst')
    gal = ui.create_model_component('xswabs', 'gal')
    gal.nh = 0.04
    ui.set_source('tst', gal * orig_mdl)

    # Ensure near the minimum
    ui.fit()
    ui.covar()
    covmat = ui.get_covar_results().extra_output
    errs = np.sqrt(covmat.diagonal())

    # See commentary in test_pha1_plot_foo_flux_model about the
    # potentially-surprising behavior of the return value of the
    # get routines

    # Absorbed flux
    res = getfunc(lo=0.5, hi=2, num=1000, bins=19, correlated=False)

    avals = res.modelvals
    assert avals.shape == (1000, 3)

    std1 = np.std(avals[:, 1])
    std2 = np.log10(np.std(avals[:, 2]))
    assert std1 == pytest.approx(0.13478302893162564, rel=1e-3)
    assert std2 == pytest.approx(-4.518960679037794, rel=1e-3)

    assert res.clipped.shape == (1000,)

    assert res.y.shape == (20,)

    aflux = np.median(res.flux)

    # Unabsorbed flux
    res = getfunc(lo=0.5, hi=2, model=orig_mdl, num=1000, bins=21,
                  correlated=False)

    uvals = res.modelvals
    assert uvals.shape == (1000, 3)

    std1 = np.std(uvals[:, 1])
    std2 = np.log10(np.std(uvals[:, 2]))
    assert std1 == pytest.approx(0.13896034019912706, rel=1e-3)
    assert std2 == pytest.approx(-4.534244395838702, rel=1e-3)

    assert res.y.shape == (22,)

    uflux = np.median(res.flux)

    # Is the unabsorbed to absorbed median flux close to the value
    # calculated from the best-fit solutions (specified as a number
    # rather than calculated here to act as a regression test).
    #
    got = uflux / aflux
    assert got == pytest.approx(ratio, rel=1e-3)


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title,plotcls",
                         [("model", [], "Model",
                           sherpa.plot.ModelPlot),
                          ("model_component", ['mdl'],
                           "Model component: polynom1d.mdl",
                           sherpa.plot.ComponentModelPlot),
                          ("source", [], "Source",
                           sherpa.plot.SourcePlot),
                          ("source_component", ['mdl'],
                           "Source model component: polynom1d.mdl",
                           sherpa.plot.ComponentSourcePlot)])
def test_data1d_get_model_plot(cls, plottype, extraargs, title, plotcls):
    """Check we can get a motel plot of a Data1D model/source.

    Note that we test both Session classes here.
    """

    x = np.arange(10, 20, 2)
    y = np.asarray([100, 102, 104, 106, 108])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, x, y, Data1D)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 90
    mdl.c1 = 1

    s.set_source(mdl)

    plot = getattr(s, "get_{}_plot".format(plottype))(*extraargs)

    assert isinstance(plot, plotcls)
    assert plot.title == title
    assert plot.x == pytest.approx(x)
    assert plot.y == pytest.approx(y)


def create_template():
    """Create a simple template model"""

    # Evalation grid
    x = np.arange(5, 25, 1)

    # The parameter for the grid (FWHM in this case)
    grid = np.arange(1, 10, 1)

    mdl = basic.Gauss1D()
    mdl.pos = 15
    mdl.ampl = 100

    templates = []
    for fwhm in grid:
        mdl.fwhm = fwhm
        template = basic.TableModel()
        template.load(x, mdl(x))
        templates.append(template)

    grid.resize((grid.size, 1))
    return create_template_model('tmdl', ["fwhm"], grid, templates)


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title,plotcls",
                         [("model", [], "Model",
                           sherpa.plot.ModelPlot),
                          ("model_component", ['tmdl'],
                           "Model component: template.tmdl",
                           sherpa.plot.ComponentModelPlot),
                          ("source", [], "Source",
                           sherpa.plot.SourcePlot),
                          ("source_component", ['tmdl'],
                           "Source model component: template.tmdl",
                           sherpa.plot.ComponentSourcePlot)])
def test_data1d_get_model_plot_template(cls, plottype, extraargs, title, plotcls):
    """Template models are handled slightly differently.
    """

    x = np.arange(10, 20, 2)
    y = np.asarray([100, 102, 104, 106, 108])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, x, y, Data1D)

    mdl = create_template()
    s._tbl_models.append(mdl)
    s._add_model_component(mdl)

    s.set_source(mdl)

    # best to pick a value on the grid used by create_template
    mdl.fwhm = 3

    plot = getattr(s, "get_{}_plot".format(plottype))(*extraargs)

    assert isinstance(plot, plotcls)
    assert plot.title == title
    assert plot.x == pytest.approx(x)

    # This is the template value, so it should be close to this
    # (although there's interpolation both on the grid of parameter
    # values and the fata grid).
    #
    mdl = basic.Gauss1D()
    mdl.pos = 15
    mdl.fwhm = 3
    mdl.ampl = 100
    yexp = mdl(x)

    assert plot.y == pytest.approx(yexp)

    # Although this class is intended for use, it currently isn't
    # so check it isn't being used (only for the source_component
    # case but leave the test for everything).
    #
    assert not isinstance(plot, sherpa.plot.ComponentTemplateSourcePlot)


@requires_pylab
@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title",
                         [("model", [], "Model"),
                          ("model_component", ['mdl'],
                           "Model component: polynom1d.mdl"),
                          ("source", [], "Source"),
                          ("source_component", ['mdl'],
                           "Source model component: polynom1d.mdl")])
def test_data1d_plot_model(cls, plottype, extraargs, title):
    """Check we can plot a Data1D model/source.

    For the Data1D case source and model return the
    same plots apart for the title.

    Note that we test both Session classes here.
    """

    from matplotlib import pyplot as plt

    x = np.arange(10, 20, 2)
    y = np.asarray([100, 110, 105, 95, 120])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, x, y, Data1D)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 90
    mdl.c1 = 1

    s.set_source(mdl)

    getattr(s, "plot_{}".format(plottype))(*extraargs)

    ax = plt.gca()
    assert ax.get_xscale() == 'linear'
    assert ax.get_yscale() == 'linear'

    assert ax.get_xlabel() == 'x'
    assert ax.get_ylabel() == 'y'
    assert ax.get_title() == title

    assert len(ax.lines) == 1

    # The line
    line = ax.lines[0]
    assert line.get_xdata().size == 5

    xplot = line.get_xdata()
    yplot = line.get_ydata()

    assert xplot == pytest.approx(x)
    assert yplot == pytest.approx(np.asarray([100, 102, 104, 106, 108]))


@requires_pylab
@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title",
                         [("model", [], "Model"),
                          ("model_component", ['tmdl'],
                           "Model component: template.tmdl"),
                          ("source", [], "Source"),
                          ("source_component", ['tmdl'],
                           "Source model component: template.tmdl")])
def test_data1d_plot_model_template(cls, plottype, extraargs, title):
    """Template models are handled slightly differently.
    """

    from matplotlib import pyplot as plt

    x = np.arange(10, 20, 2)
    y = np.asarray([100, 102, 104, 106, 108])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, x, y, Data1D)

    mdl = create_template()
    s._tbl_models.append(mdl)
    s._add_model_component(mdl)

    s.set_source(mdl)

    # best to pick a value on the grid used by create_template
    mdl.fwhm = 3

    plot = getattr(s, "plot_{}".format(plottype))
    plot(*extraargs)

    ax = plt.gca()
    assert ax.get_xscale() == 'linear'
    assert ax.get_yscale() == 'linear'

    assert ax.get_xlabel() == 'x'
    assert ax.get_ylabel() == 'y'
    assert ax.get_title() == title

    assert len(ax.lines) == 1
    line = ax.lines[0]
    assert line.get_xdata().size == 5

    xplot = line.get_xdata()
    yplot = line.get_ydata()

    assert xplot == pytest.approx(x)

    # This is the template value, so it should be close to this
    # (although there's interpolation both on the grid of parameter
    # values and the fata grid).
    #
    mdl = basic.Gauss1D()
    mdl.pos = 15
    mdl.fwhm = 3
    mdl.ampl = 100
    yexp = mdl(x)

    assert yplot == pytest.approx(yexp)


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title,plotcls",
                         [("model", [], "Model",
                           sherpa.plot.ModelHistogramPlot),
                          ("model_component", ['mdl'],
                           "Model component: polynom1d.mdl",
                           sherpa.plot.ComponentModelHistogramPlot),
                          ("source", [], "Source",
                           sherpa.plot.SourceHistogramPlot),
                          ("source_component", ['mdl'],
                           "Source model component: polynom1d.mdl",
                           sherpa.plot.ComponentSourceHistogramPlot)])
def test_data1dint_get_model_plot(cls, plottype, extraargs, title, plotcls):
    """Check we can plot a Data1DInt model.

    This is plotted as a histogram, not as x/y values like Data1D

    For the Data1DInt case source and model return the
    same plots apart for the title.

    Note that we test both Session classes here.
    """

    xlo = np.asarray([10, 12, 16, 20, 22])
    xhi = np.asarray([12, 16, 19, 22, 26])
    y = np.asarray([50, 54, 58, 60, 64])
    yexp = np.asarray([22, 56, 52.5, 42, 96])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, xlo, xhi, y, Data1DInt)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 0
    mdl.c1 = 1

    s.set_source(mdl)

    plot = getattr(s, "get_{}_plot".format(plottype))(*extraargs)

    assert isinstance(plot, plotcls)
    assert plot.title == title
    # NOTE: we add an x mid-point
    assert plot.x == pytest.approx((xlo + xhi) / 2)
    assert plot.xlo == pytest.approx(xlo)
    assert plot.xhi == pytest.approx(xhi)
    assert plot.y == pytest.approx(yexp)


@requires_pylab
@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("plottype,extraargs,title",
                         [("model", [], "Model"),
                          ("model_component", ['mdl'],
                           "Model component: polynom1d.mdl"),
                          ("source", [], "Source"),
                          ("source_component", ['mdl'],
                           "Source model component: polynom1d.mdl")])
def test_data1dint_plot_model(cls, plottype, extraargs, title):
    """Check we can plot a Data1DInt model.

    This is plotted as a histogram, not as x/y values like Data1D

    For the Data1DInt case source and model return the
    same plots apart for the title.

    Note that we test both Session classes here.
    """

    from matplotlib import pyplot as plt

    xlo = np.asarray([10, 12, 16, 20, 22])
    xhi = np.asarray([12, 16, 19, 22, 26])
    y = np.asarray([50, 54, 58, 60, 64])

    s = cls()
    s._add_model_types(basic)

    s.load_arrays(1, xlo, xhi, y, Data1DInt)

    mdl = s.create_model_component('polynom1d', 'mdl')
    mdl.c0 = 0
    mdl.c1 = 1

    s.set_source(mdl)

    getattr(s, "plot_{}".format(plottype))(*extraargs)

    ax = plt.gca()
    assert ax.get_xscale() == 'linear'
    assert ax.get_yscale() == 'linear'

    assert ax.get_xlabel() == 'x'
    assert ax.get_ylabel() == 'y'
    assert ax.get_title() == title

    assert len(ax.lines) == 2

    # The line
    line = ax.lines[0]
    assert line.get_xdata().size == 11

    xplot = line.get_xdata()
    yplot = line.get_ydata()

    # do it a different way to the plot backend
    xexp = np.dstack((xlo, xhi)).flatten()
    yexp = np.asarray([22, 56, 52.5, 42, 96]).repeat(2)

    # Check for nan
    #
    good = np.ones(11, dtype=bool)
    good[6] = False
    assert np.isfinite(xplot) == pytest.approx(good)
    assert np.isfinite(yplot) == pytest.approx(good)

    assert xplot[good] == pytest.approx(xexp)
    assert yplot[good] == pytest.approx(yexp)

    # The points
    pts = ax.lines[1]
    assert pts.get_xdata().size == 5

    xplot = pts.get_xdata()
    yplot = pts.get_ydata()

    # do it a different way to the plot backend
    xexp = (xlo + xhi) / 2
    yexp = np.asarray([22, 56, 52.5, 42, 96])

    assert xplot == pytest.approx(xexp)
    assert yplot == pytest.approx(yexp)


# Tests from test_ui_plot but now checking out DataPHA
#
@requires_fits
@requires_data
def test_pha1_data_plot_recalc(clean_astro_ui, basic_pha1):
    """Basic testing of get_data_plot(recalc=False)"""

    ui.get_data_plot()

    ui.ignore(None, 1)
    ui.ignore(5, None)

    p = ui.get_data_plot(recalc=False)
    assert isinstance(p, DataPHAPlot)
    assert p.xlo.size == 42
    assert (p.xlo[0] + p.xhi[0]) / 2 == pytest.approx(0.5183)
    assert (p.xlo[-1] + p.xhi[-1]) / 2 == pytest.approx(8.2198)

    p = ui.get_data_plot(recalc=True)
    assert isinstance(p, DataPHAPlot)
    assert p.xlo.size == 26
    assert (p.xlo[0] + p.xhi[0]) / 2 == pytest.approx(1.0658)
    assert (p.xlo[-1] + p.xhi[-1]) / 2 == pytest.approx(4.4822)


@pytest.mark.parametrize("getfunc", [ui.get_data_plot, ui.get_model_plot])
@pytest.mark.parametrize("idval", [None, 'x'])
def test_get_xxx_plot_nodata(getfunc, idval, clean_astro_ui):

    with pytest.raises(IdentifierErr) as exc:
        getfunc(idval)

    iid = 1 if idval is None else idval
    assert str(exc.value) == 'data set {} has not been set'.format(iid)


@requires_fits
@requires_data
@pytest.mark.parametrize("ptype,extraargs,pclass,nbins1,nbins2",
                         [('model', [], ModelHistogram, 644, 252),
                          ('model_component', ['pl'], ComponentModelPlot, 644, 252),
                          ('source', [], SourcePlot, 1090, 1090),
                          ('source_component', ['pl'], ComponentSourcePlot, 1090, 1090)])
def test_pha1_model_plot_recalc(ptype, extraargs, pclass, nbins1, nbins2,
                                clean_astro_ui, basic_pha1):
    """Basic testing of get_model_plot(recalc=False)"""

    func = getattr(ui, 'get_{}_plot'.format(ptype))

    # Seed the data for the recalc=False call
    func(*extraargs)

    ui.ignore(None, 1)
    ui.ignore(5, None)

    p = func(*extraargs, recalc=False)
    assert isinstance(p, pclass)
    assert p.xlo.size == nbins1

    p = func(*extraargs, recalc=True)
    assert isinstance(p, pclass)
    assert p.xlo.size == nbins2


@requires_fits
@requires_data
def test_pha1_fit_plot_recalc(clean_astro_ui, basic_pha1):
    """Basic testing of get_fit_plot(recalc=False)"""

    ui.get_fit_plot('tst')

    ui.ignore(None, 1)
    ui.ignore(5, None)

    p = ui.get_fit_plot('tst', recalc=False)
    assert isinstance(p, FitPlot)
    assert isinstance(p.dataplot, DataPHAPlot)
    assert isinstance(p.modelplot, ModelPHAHistogram)
    assert p.dataplot.xlo.size == 42
    assert p.modelplot.xlo.size == 42

    p = ui.get_fit_plot('tst', recalc=True)
    assert isinstance(p, FitPlot)
    assert isinstance(p.dataplot, DataPHAPlot)
    assert isinstance(p.modelplot, ModelPHAHistogram)
    assert p.dataplot.xlo.size == 26
    assert p.modelplot.xlo.size == 26


@requires_fits
@requires_data
def test_pha1_bkg_fit_plot_no_model(clean_astro_ui, basic_pha1):
    """Error out if get_bkg_fit_plot has no source model"""

    ui.get_fit_plot('tst')
    with pytest.raises(ModelErr) as exc:
        ui.get_bkg_fit_plot('tst')

    emsg = 'background model 1 for data set tst has not been set'
    assert str(exc.value) == emsg


@requires_fits
@requires_data
def test_pha1_bkg_fit_plot_recalc(clean_astro_ui, make_data_path):
    """Basic testing of get_bkg_fit_plot(recalc=False)"""

    infile = make_data_path('3c273.pi')
    ui.load_pha(infile)
    ui.get_bkg().name = 'my-name.pi'

    ui.ignore(None, 1)
    ui.ignore(5, None)

    ui.set_bkg_model(ui.polynom1d.bpl)

    p = ui.get_bkg_fit_plot(recalc=False)
    assert isinstance(p, FitPlot)
    assert p.dataplot is None
    assert p.modelplot is None

    p = ui.get_bkg_fit_plot(recalc=True)
    assert isinstance(p, FitPlot)
    assert isinstance(p.dataplot, DataPHAPlot)
    assert isinstance(p.modelplot, ModelPHAHistogram)
    assert p.dataplot.xlo.size == 26
    assert p.modelplot.xlo.size == 26
    assert p.dataplot.title == 'my-name.pi'
    assert p.modelplot.title == 'Background Model Contribution'


@requires_pylab
@requires_fits
@requires_data
def test_pha1_plot_multiple_args(clean_astro_ui, basic_pha1):
    """Check plot('bkg', 1, 1, 'bkg', 1, 2)

    The dataid is actually 'tst', and the way we check
    it has worked is by getting it to check for an unknown
    background component: if the second argument is
    recognized then we get an error.

    Of course, just checking ('bkg', 'tst', 1) is in some
    ways enough.

    """

    with pytest.raises(IdentifierErr) as exc:
        ui.plot('bkg', 'tst', 1, 'bkg', 'tst', 'up')

    emsg = 'background data set up in PHA data set tst has not been set'
    assert str(exc.value) == emsg


def example_data1d():
    """Create a Data1D object"""
    x = np.asarray([2, 5, 20])
    y = np.asarray([2, 20, 200])
    return Data1D('name1d', x, y)


def example_data1dint():
    """Create a Data1DInt object"""
    x1 = np.asarray([2, 5, 20])
    x2 = np.asarray([5, 10, 40])
    y = np.asarray([2, 20, 200])
    return Data1DInt('name1d int', x1, x2, y)


def example_datapha():
    """Create a DataPHA object.

    The background is set to be equal to the source just
    because it is easier.
    """
    chans = np.arange(1, 6, dtype=np.int16)
    d = DataPHA('example.pha', chans, chans)

    ebins = np.arange(1, 7)
    arf = create_arf(ebins[:-1], ebins[1:])

    d.set_arf(arf)

    # re-create so we don't have a circular dataset
    b = DataPHA('background.pha', chans, chans)
    d.set_background(b)

    return d


def test_datapha_plot_after_clean():
    """Check we can get the DataPHA plot back after a clean."""

    d = example_datapha()

    s = sherpa.astro.ui.utils.Session()
    s.set_data(d)

    d1 = s.get_data_plot()

    # Change the yerrorbars setting. It depends whether we
    # have a plot backend or not.
    #
    prefs = s.get_data_plot_prefs()
    have_backend = sherpa.plot.backend.name != 'dummy'
    if have_backend:
        assert prefs['yerrorbars']
    else:
        assert 'yerrorbars' not in prefs

    prefs['yerrorbars'] = False

    assert isinstance(d1, DataPHAPlot)
    assert not d1.histo_prefs['yerrorbars']

    s.clean()
    s.set_data(d)

    # Check the yerrorbars setting is back to True
    #
    d2 = s.get_data_plot()
    assert isinstance(d2, DataPHAPlot)

    prefs = s.get_data_plot_prefs()

    if have_backend:
        assert prefs['yerrorbars']
        assert d2.histo_prefs['yerrorbars']
    else:
        assert 'yerrorbars' not in prefs


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("datafunc", [example_data1d,
                                      example_data1dint])
@pytest.mark.parametrize("plotfunc",
                         ['data',
                          'model',
                          'source',
                          'fit',
                          'resid',
                          'ratio',
                          'delchi'])
def test_set_plot_opt_x(cls, datafunc, plotfunc):
    """Does set_xlog/xlinear work?

    We run for both session types here, rather than duplicate
    the code across tests. We run both the log and linear
    versions in the same test since it makes it easy to check
    if the call worked.
    """

    s = cls()
    s._add_model_types(basic)

    s.set_xlog()

    # load the data *after* calling the log method, as it should
    # not matter what the data type is
    #
    s.set_data(datafunc())
    is_int = hasattr(s.get_data(), 'xlo')

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(mdl)

    plot = getattr(s, 'plot_{}'.format(plotfunc))
    pdata = getattr(s, 'get_{}_plot'.format(plotfunc))

    # Create the plot
    plot()
    p1 = pdata()
    if plotfunc == 'fit':
        if is_int:
            assert p1.dataplot.histo_prefs['xlog']
            assert p1.modelplot.histo_prefs['xlog']
        else:
            assert p1.dataplot.plot_prefs['xlog']
            assert p1.modelplot.plot_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi']:
        assert p1.plot_prefs['xlog']
    elif is_int:
        assert p1.histo_prefs['xlog']
    else:
        assert p1.plot_prefs['xlog']

    s.set_xlinear()
    plot()

    # Technically not needed as p1 is the same as p2
    p2 = pdata()
    if plotfunc == 'fit':
        if is_int:
            assert not p2.dataplot.histo_prefs['xlog']
            assert not p2.modelplot.histo_prefs['xlog']
        else:
            assert not p2.dataplot.plot_prefs['xlog']
            assert not p2.modelplot.plot_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi']:
        assert not p2.plot_prefs['xlog']
    elif is_int:
        assert not p2.histo_prefs['xlog']
    else:
        assert not p2.plot_prefs['xlog']


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("datafunc", [example_data1d,
                                      example_data1dint])
@pytest.mark.parametrize("plotfunc,answer",
                         [('data', True),
                          ('model', True),
                          ('source', True),
                          ('fit', True),
                          ('resid', False),
                          ('ratio', False),
                          ('delchi', False)])
def test_set_plot_opt_y(cls, datafunc, plotfunc, answer):
    """Does set_ylog/ylinear work?

    Unlike test_set_plot_opt_x, the Y axis setting of the plot
    does not necessarily follow the set_ylog setting (e.g.
    the residual plots).
    """

    s = cls()
    s._add_model_types(basic)

    s.set_ylog()

    # load the data *after* calling the log method, as it should
    # not matter what the data type is
    #
    s.set_data(datafunc())
    is_int = hasattr(s.get_data(), 'xlo')

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(mdl)

    plot = getattr(s, 'plot_{}'.format(plotfunc))
    pdata = getattr(s, 'get_{}_plot'.format(plotfunc))

    plot()
    p1 = pdata()
    if plotfunc == 'fit':
        if is_int:
            assert p1.dataplot.histo_prefs['ylog'] == answer
            assert p1.modelplot.histo_prefs['ylog'] == answer
        else:
            assert p1.dataplot.plot_prefs['ylog'] == answer
            assert p1.modelplot.plot_prefs['ylog'] == answer
    elif plotfunc in ['resid', 'ratio', 'delchi']:
        # We use plot_prefs even if is_int is True for
        # the residual-style plots. That is, the ordering
        # of the checks here is important.
        #
        assert p1.plot_prefs['ylog'] == answer
    elif is_int:
        assert p1.histo_prefs['ylog'] == answer
    else:
        assert p1.plot_prefs['ylog'] == answer

    s.set_ylinear()
    plot()

    # Technically not needed as p1 is the same as p2
    p2 = pdata()
    if plotfunc == 'fit':
        if is_int:
            assert not p2.dataplot.histo_prefs['ylog']
            assert not p2.modelplot.histo_prefs['ylog']
        else:
            assert not p2.dataplot.plot_prefs['ylog']
            assert not p2.modelplot.plot_prefs['ylog']
    elif plotfunc in ['resid', 'ratio', 'delchi']:
        assert not p2.plot_prefs['ylog']
    elif is_int:
        assert not p2.histo_prefs['ylog']
    else:
        assert not p2.plot_prefs['ylog']


@pytest.mark.parametrize("cls",
                         [sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("datafunc", [example_datapha])
@pytest.mark.parametrize("plotfunc",
                         ['data',
                          'model',
                          'source',
                          'fit',
                          'resid',
                          'ratio',
                          'delchi',
                          'bkg',
                          'bkg_model',
                          'bkg_source',
                          'bkg_fit',
                          'bkg_resid',
                          'bkg_ratio',
                          'bkg_delchi'])
def test_set_plot_opt_x_astro(cls, datafunc, plotfunc):
    """Does set_xlog/xlinear work? Astro data objects only.

    Given that cls and datafunc are single-element lists they could be
    hard-coded, but leave as is to mirror test_set_plot_opt.
    """

    s = cls()
    s._add_model_types(basic)

    s.set_xlog()

    # load the data *after* calling the log method, as it should
    # not matter what the data type is
    #
    s.set_data(datafunc())

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(mdl)
    s.set_bkg_source(mdl)

    plot = getattr(s, 'plot_{}'.format(plotfunc))
    pdata = getattr(s, 'get_{}_plot'.format(plotfunc))

    # Create the plot
    plot()
    p1 = pdata()
    if plotfunc in ['fit', 'bkg_fit']:
        assert p1.dataplot.histo_prefs['xlog']
        # Check the current behavior of the model plot in case it changes.
        # Note that for DataPHA datasets the modelplot object is
        # created on-the-fly using sherpa.astro.plot.ModelPHAHistogram,
        # rather than a session._plotobj value that is also set in
        # session._plot_types, so it doesn't get changed by set_xlog etc.
        #
        # Ideally this would match the dataplot setting but it isn't
        # actually required for the plot to work.
        #
        if sherpa.plot.backend.name != 'dummy':
            assert not p1.modelplot.histo_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi', 'bkg_resid', 'bkg_ratio', 'bkg_delchi']:
        assert p1.plot_prefs['xlog']
    else:
        assert p1.histo_prefs['xlog']

    s.set_xlinear()
    plot()

    # Technically not needed as p1 is the same as p2
    p2 = pdata()
    if plotfunc in ['fit', 'bkg_fit']:
        assert not p2.dataplot.histo_prefs['xlog']
        if sherpa.plot.backend.name != 'dummy':
            assert not p2.modelplot.histo_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi', 'bkg_resid', 'bkg_ratio', 'bkg_delchi']:
        assert not p2.plot_prefs['xlog']
    else:
        assert not p2.histo_prefs['xlog']


@pytest.mark.parametrize("cls",
                         [sherpa.astro.ui.utils.Session])
@pytest.mark.parametrize("datafunc", [example_datapha])
@pytest.mark.parametrize("plotfunc,answer",
                         [('data', True),
                          ('model', True),
                          ('source', True),
                          ('fit', True),
                          ('resid', False),
                          ('ratio', False),
                          ('delchi', False),
                          ('bkg', True),
                          ('bkg_model', True),
                          ('bkg_source', True),
                          ('bkg_fit', True),
                          ('bkg_resid', False),
                          ('bkg_ratio', False),
                          ('bkg_delchi', False)])
def test_set_plot_opt_y_astro(cls, datafunc, plotfunc, answer):
    """Does set_ylog/ylinear work?  Astro data objects only.
    """

    s = cls()
    s._add_model_types(basic)

    s.set_ylog()

    # load the data *after* calling the log method, as it should
    # not matter what the data type is
    #
    s.set_data(datafunc())

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(mdl)
    s.set_bkg_source(mdl)

    plot = getattr(s, 'plot_{}'.format(plotfunc))
    pdata = getattr(s, 'get_{}_plot'.format(plotfunc))

    plot()
    p1 = pdata()
    if plotfunc in ['fit', 'bkg_fit']:
        assert p1.dataplot.histo_prefs['ylog'] == answer
        # Check the current behavior of the model plot in case it changes
        if sherpa.plot.backend.name != 'dummy':
            assert not p1.modelplot.histo_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi', 'bkg_resid', 'bkg_ratio', 'bkg_delchi']:
        assert p1.plot_prefs['ylog'] == answer
    else:
        assert p1.histo_prefs['ylog'] == answer

    s.set_ylinear()
    plot()

    # Technically not needed as p1 is the same as p2
    p2 = pdata()
    if plotfunc in ['fit', 'bkg_fit']:
        assert not p2.dataplot.histo_prefs['ylog']
        if sherpa.plot.backend.name != 'dummy':
            assert not p2.modelplot.histo_prefs['xlog']
    elif plotfunc in ['resid', 'ratio', 'delchi', 'bkg_resid', 'bkg_ratio', 'bkg_delchi']:
        assert not p2.plot_prefs['ylog']
    else:
        assert not p2.histo_prefs['ylog']


@requires_pylab
def test_set_plot_opt_with_plot_x():
    """Does set_xlog/xlinear work with plot()  Astro data objects only.

    We could repeat the other tests - i.e. query the
    plot objects - but that is a bit messy, so require
    matplotlib.
    """

    from matplotlib import pyplot as plt

    s = sherpa.astro.ui.utils.Session()
    s._add_model_types(basic)

    s.set_xlog()

    s.set_data(1, example_data1d())
    s.set_data(2, example_data1dint())
    s.set_data(3, example_datapha())

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(1, mdl)
    s.set_source(2, mdl)
    s.set_source(3, mdl)

    s.plot('fit', 1, 'fit', 2, 'fit', 3)

    fig = plt.gcf()

    assert len(fig.axes) == 3
    for idx, ax in enumerate(fig.axes):
        assert ax.get_xscale() == 'log', idx
        assert ax.get_yscale() == 'linear', idx

    s.set_xlinear()

    s.plot('fit', 1, 'fit', 2, 'fit', 3)

    fig = plt.gcf()

    assert len(fig.axes) == 3
    for idx, ax in enumerate(fig.axes):
        assert ax.get_xscale() == 'linear', idx
        assert ax.get_yscale() == 'linear', idx


@requires_pylab
def test_set_plot_opt_with_plot_y():
    """Does set_ylog/ylinear work with plot()  Astro data objects only.
    """

    from matplotlib import pyplot as plt

    s = sherpa.astro.ui.utils.Session()
    s._add_model_types(basic)

    s.set_ylog()

    s.set_data(1, example_data1d())
    s.set_data(2, example_data1dint())
    s.set_data(3, example_datapha())

    # Set up a model, in case we need it.
    mdl = s.create_model_component('polynom1d', 'm1')
    mdl.c0 = 1
    mdl.c1 = 1
    s.set_source(1, mdl)
    s.set_source(2, mdl)
    s.set_source(3, mdl)

    s.plot('fit', 1, 'fit', 2, 'fit', 3)

    fig = plt.gcf()

    assert len(fig.axes) == 3
    for idx, ax in enumerate(fig.axes):
        assert ax.get_xscale() == 'linear', idx
        assert ax.get_yscale() == 'log', idx

    s.set_ylinear()

    s.plot('fit', 1, 'fit', 2, 'fit', 3)

    fig = plt.gcf()

    assert len(fig.axes) == 3
    for idx, ax in enumerate(fig.axes):
        assert ax.get_xscale() == 'linear', idx
        assert ax.get_yscale() == 'linear', idx


@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
def test_set_opt_invalid(cls):
    """Check we error out if called with an invalid option"""

    s = cls()
    with pytest.raises(PlotErr) as exc:
        s.set_xlog('notdata')

    msg = "Plot type 'notdata' not found in ["
    assert str(exc.value).startswith(msg)


@requires_pylab
@pytest.mark.parametrize("cls",
                         [sherpa.ui.utils.Session, sherpa.astro.ui.utils.Session])
def test_set_plot_opt_explicit(cls):
    """Check we can call set_xlog('data').

    We don't check all options (unlike the set_xlog/ylog
    tests above) since we assume they work. This is specific
    to setting just the data options.
    """

    from matplotlib import pyplot as plt

    s = cls()
    s._add_model_types(basic)

    s.set_xlog('data')

    d1 = example_data1d()
    d2 = example_data1dint()

    s.set_data(1, d1)
    s.set_data(2, d2)

    mdl = s.create_model_component('polynom1d', 'm1')
    s.set_source(1, mdl)
    s.set_source(2, mdl)

    s.plot('data', 'model', 'data', 2, 'model', 2)

    fig = plt.gcf()

    assert len(fig.axes) == 4

    for idx, ax in enumerate(fig.axes[0:4:2]):
        assert ax.get_xscale() == 'log', idx
        assert ax.get_yscale() == 'linear', idx

    for idx, ax in enumerate(fig.axes[1:4:2]):
        assert ax.get_xscale() == 'linear', idx
        assert ax.get_yscale() == 'linear', idx


@requires_pylab
def test_set_plot_opt_explicit_astro():
    """Check we can call set_xlog('data') with astro data.

    We don't check all options (unlike the set_xlog/ylog
    tests above) since we assume they work. This is specific
    to setting just the data options.
    """

    from matplotlib import pyplot as plt

    s = sherpa.astro.ui.utils.Session()
    s._add_model_types(basic)

    s.set_xlog('data')

    d1 = example_datapha()

    s.set_data(1, d1)

    mdl = s.create_model_component('polynom1d', 'm1')
    s.set_source(1, mdl)
    s.set_bkg_source(1, mdl)

    s.plot('data', 'model', 'bkg', 'bkgmodel')

    fig = plt.gcf()
    assert len(fig.axes) == 4

    # Only data has X axis drawn in log scale. We may decide to
    # make the 'bkg' plot act like a "data" plot in the future.
    #
    assert fig.axes[0].get_xscale() == 'log'
    assert fig.axes[0].get_yscale() == 'linear'

    assert fig.axes[1].get_xscale() == 'linear'
    assert fig.axes[1].get_yscale() == 'linear'

    assert fig.axes[2].get_xscale() == 'linear'
    assert fig.axes[2].get_yscale() == 'linear'

    assert fig.axes[3].get_xscale() == 'linear'
    assert fig.axes[3].get_yscale() == 'linear'

    plt.close(fig)

    # Now set the model y axis to log
    #
    s.set_ylog('model')

    s.plot('data', 'model', 'bkg', 'bkgmodel')

    fig = plt.gcf()
    assert len(fig.axes) == 4

    assert fig.axes[0].get_xscale() == 'log'
    assert fig.axes[0].get_yscale() == 'linear'

    assert fig.axes[1].get_xscale() == 'linear'
    assert fig.axes[1].get_yscale() == 'log'

    assert fig.axes[2].get_xscale() == 'linear'
    assert fig.axes[2].get_yscale() == 'linear'

    assert fig.axes[3].get_xscale() == 'linear'
    assert fig.axes[3].get_yscale() == 'linear'

    plt.close(fig)


def check_plot2_xscale(xscale):
    """Are there two plots, y-axis linear, x axis set?

    Any test using this needs @requires_pylab
    """

    from matplotlib import pyplot as plt

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 2
    assert axes[0].xaxis.get_label().get_text() == ''

    assert axes[0].xaxis.get_scale() == xscale
    assert axes[0].yaxis.get_scale() == 'linear'

    assert axes[1].xaxis.get_scale() == xscale
    assert axes[1].yaxis.get_scale() == 'linear'


@requires_pylab
@pytest.mark.parametrize("idval", [None, "bob"])
@pytest.mark.parametrize("plottype,xscale", [('data', 'log'),
                                             ('resid', 'log'),
                                             ('bkg', 'linear'),
                                             ('bkgresid', 'linear')])
def test_plot_fit_resid_set_xlog(idval, plottype, xscale, clean_astro_ui):
    """Check that set_xlog handling for plot_fit_resid.

    What is the X-axis scaling when you call set_xlog(plottype)?
    We use a range of plot types as the behavior is not always
    obvious, so let's ensure we test them.
    """

    setup_example(idval)
    ui.set_xlog(plottype)
    ui.plot_fit_resid(idval)
    check_plot2_xscale(xscale)


@requires_pylab
@pytest.mark.parametrize("idval", [None, "bob"])
@pytest.mark.parametrize("plottype,xscale", [('data', 'linear'),
                                             ('resid', 'linear'),
                                             ('bkg', 'log'),
                                             ('bkgresid', 'log')])
def test_plot_bkg_fit_resid_set_xlog(idval, plottype, xscale, clean_astro_ui):
    """Check that set_xlog handling for plot_bkg_fit_resid.

    This logic could be added to test_plot_fit_resid_set_xlog but it
    would complicate the setup, so we duplicate the code.

    """

    setup_example_bkg_model(idval)
    ui.set_xlog(plottype)
    ui.plot_bkg_fit_resid(idval)
    check_plot2_xscale(xscale)


@requires_pylab
@pytest.mark.parametrize("plot,yscale", [('data', 'linear'),
                                         ('ratio', 'linear'),
                                         ('fit', 'linear'),
                                         ('bkg', 'log'),
                                         ('bkg_resid', 'linear'),
                                         ('bkg_fit', 'log')])
def test_set_ylog_bkg(plot, yscale, clean_astro_ui):
    """Check y axis after_ylog('bkg').

    The idea is to check how separate the "background" from
    "non-background" plots are.  I use ylog since other tests have
    used xlog.

    """

    from matplotlib import pyplot as plt

    setup_example_bkg_model(1)
    ui.set_ylog('bkg')

    pfunc = getattr(ui, f'plot_{plot}')
    pfunc()

    fig = plt.gcf()
    axes = fig.axes
    assert len(axes) == 1
    assert axes[0].xaxis.get_scale() == 'linear'
    assert axes[0].yaxis.get_scale() == yscale


@requires_plotting
def test_pha_model_plot_filter_range_manual_1024(clean_astro_ui):
    """Check if issue #1024 is fixed.

    The problem was that the PHA model plot was resetting the
    range so that the notice call ended up not changing anything.
    Unfortunately this does not show up in the "manual" case,
    because the data is not grouped, so you have to rely on
    test_pha_model_plot_filter_range_1024.
    """

    setup_example(None)
    ui.set_analysis('energy')

    ui.plot_model()
    ui.notice(0.77, 1.125)

    f = ui.get_filter()
    assert f == '0.750000000000:1.300000000000'


@requires_fits
@requires_data
@requires_plotting
def test_pha_model_plot_filter_range_1024(make_data_path, clean_astro_ui):
    """Check if issue #1024 is fixed.

    Unlike test_pha_model_plot_filter_range_manual_1024
    this test does show issue #1024.
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_source(ui.powlaw1d.pl)

    ui.plot_model()
    ui.notice(0.5, 5)

    f = ui.get_filter()
    assert f == '0.467200011015:5.022399902344'


@requires_fits
@requires_data
@requires_plotting
@pytest.mark.parametrize("mask", [True, np.ones(46, dtype=bool)])
def test_pha_model_plot_filter_range_1024_true(mask, make_data_path, clean_astro_ui):
    """Special-case handling of mask: all selected.
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_source(ui.powlaw1d.pl)

    d = ui.get_data()
    d.mask = mask
    ui.plot_model()
    f = ui.get_filter()
    assert f == '0.001460000058:14.950400352478'


@requires_fits
@requires_data
@requires_plotting
@pytest.mark.parametrize("mask,expected",
                         [(False, 'No noticed bins'),
                          (np.zeros(46, dtype=bool), '')])
def test_pha_model_plot_filter_range_1024_false(mask, expected, make_data_path, clean_astro_ui):
    """Special-case handling of mask: all masked out.
    """

    ui.load_pha(make_data_path('3c273.pi'))
    ui.set_source(ui.powlaw1d.pl)

    d = ui.get_data()
    d.mask = mask

    # See #1220 for why we don't have a unique value for the filter
    assert ui.get_filter() == expected

    # We can not guarantee what error wil be raised here because of
    # issue #1220 so just pick a generic exception which will catch
    # both SherpaErr cases and general Pythonic ones.
    #
    # Note that thanks to #1024 the array of False values doesn't
    # cause an error, but instead drops all the filters on the
    # dataset.
    #
    with pytest.raises(Exception):
        ui.plot_model()

    assert ui.get_filter() == expected


@requires_fits
@requires_data
@pytest.mark.parametrize("coord", ["logical", "image", "physical", pytest.param("world", marks=pytest.mark.xfail), pytest.param("wcs", marks=pytest.mark.xfail)])
def test_1380_plot(coord, make_data_path, clean_astro_ui):
    """The contour data should ideally remain the same.

    See also sherpa/astro/tests/test_astro_data2.py::test_1380_data

    This is the actual bug report (it only fails when a valid plot
    backend is present but we try even if there's no backend just
    to check).

    """

    infile = make_data_path("image2.fits")
    ui.load_image(infile)

    img = ui.get_data()
    assert isinstance(img, ui.DataIMG)
    assert ui.get_coord() == "logical"

    # All we do here is check we can call contour_data, not that the
    # plot has actually done anything.
    #
    ui.contour_data()

    # We can not call contour_data when the world coordinate-system
    # is set, but fortunately this is not needed to trigger #1380;
    # we just need the set_coord call.
    #
    ui.set_coord(coord)

    ui.set_coord("logical")
    ui.contour_data()
