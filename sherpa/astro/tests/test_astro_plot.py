#
#  Copyright (C) 2007, 2015, 2018, 2019  Smithsonian Astrophysical Observatory
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
from sherpa.utils.testing import requires_data, requires_fits
from sherpa.astro.data import DataPHA
from sherpa.astro.plot import BkgDataPlot, DataPlot, SourcePlot
from sherpa.models.basic import Const1D, Gauss1D
from sherpa.astro.instrument import create_arf, create_delta_rmf
from sherpa import stats

import pytest
import logging


def test_sourceplot():

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "energy"

    # use a model that is "okay" to use with keV bins
    #
    m1 = Const1D('bgnd')
    m2 = Gauss1D('abs1')
    src = 100 * m1 * (1 - m2) * 10000

    m1.c0 = 0.01
    m2.pos = 5.0
    m2.fwhm = 4.0
    m2.ampl = 0.1

    sp = SourcePlot()
    sp.prepare(data, src)

    # add in several asserts to check that something has been
    # added to the object
    #
    assert sp.xlabel == 'Energy (keV)'

    # the following depends on the backend
    # assert sp.ylabel == 'f(E)  Photons/sec/cm$^2$/keV'

    assert sp.title == 'Source Model of '

    assert sp.xlo == pytest.approx(bins[:-1])
    assert sp.xhi == pytest.approx(bins[1:])

    # The check of the values is just to check that things are going
    # as expected, so the model values have been adjusted so that
    # an "integer" check can be used with enough precision to make
    # sure that the model is being evaluated correctly, but without
    # a very-high-precision check
    #
    yexp = np.asarray([9998, 9997, 9997, 9997, 9996, 9996, 9995, 9994,
                       9994, 9993, 9992, 9991, 9990, 9988, 9987, 9985,
                       9983, 9982, 9980, 9977, 9975, 9973, 9970, 9967,
                       9964, 9961, 9958, 9955, 9951, 9948, 9944, 9941,
                       9937, 9934, 9930, 9927, 9923, 9920, 9917, 9914,
                       9911, 9909, 9907, 9905, 9903, 9902, 9901, 9900,
                       9900, 9900, 9900, 9901, 9902, 9903, 9905, 9907,
                       9909, 9911, 9914, 9917, 9920, 9923, 9927, 9930,
                       9934, 9937, 9941, 9944, 9948, 9951, 9955, 9958,
                       9961, 9964, 9967, 9970, 9973, 9975, 9977, 9980,
                       9982, 9983, 9985, 9987, 9988, 9990, 9991, 9992,
                       9993, 9994, 9994, 9995, 9996, 9996, 9997, 9997,
                       9997, 9998, 9998])

    assert (sp.y.astype(np.int) == yexp).all()
    # sp.plot()


# Low-level test of the DataPlot prepare method for PHA style analysis
# with a range of statistics. Note that the results are not checked,
# just that the call to the prepare method can be called without
# error. This test can also be run when there is no plotting backend.
#
# Extra tests could be added to check the __str__ method of DataPlot(),
# since this does query the state of the data (e.g. filtering,
# background subtraction) when creating the arrays.
#
# The pytest.param calls seem to get recorded as 2 xfails; I think
# this is for the error and because of warning messages, but it is not
# clear.
#
@requires_data
@requires_fits
@pytest.mark.parametrize("stat",
                         [None,
                          stats.Chi2(),
                          stats.Chi2ConstVar(),
                          stats.Chi2DataVar(),
                          stats.Chi2Gehrels(),
                          stats.Chi2ModVar(),
                          stats.Chi2XspecVar(),
                          stats.LeastSq(),
                          stats.Cash(),
                          stats.CStat(),
                          stats.WStat(),
                         ])
def test_astro_data_plot_with_stat_simple(make_data_path, stat):

    from sherpa.astro import io

    infile = make_data_path('3c273.pi')
    pha = io.read_pha(infile)

    # tweak the data set so that we aren't using the default
    # options (it shouldn't matter for this test but just
    # in case).
    #
    # Note that background subtraction would normally be an issue
    # for some of the stats (e.g. WStat), but this shouldn't
    # trigger a problem here.
    #
    pha.set_analysis('energy')
    pha.subtract()
    pha.ignore(None, 0.5)
    pha.ignore(7.0, None)

    dplot = DataPlot()
    dplot.prepare(pha, stat=stat)


@pytest.fixture
def setup_pha():
    """Create a PHA dataset with the given label."""

    nchans = 10
    arf0 = 0.6

    e0 = 0.1
    de = 0.125

    def _return_data(label, exposure=None, areascal=None, backscal=None):

        egrid = e0 + np.arange(nchans + 1) * de
        elo = egrid[:-1]
        ehi = egrid[1:]

        x = np.arange(1, 11)
        y = np.ones(nchans)

        arf = create_arf(elo, ehi, arf0 * np.ones(nchans))
        rmf = create_delta_rmf(elo, ehi, e_min=elo, e_max=ehi)

        data = DataPHA(label, x, y, exposure=exposure)
        data.units = "energy"
        data.rate = True

        data.set_arf(arf)
        data.set_rmf(rmf)

        data.areascal = areascal
        data.backscal = backscal
        return data

    return _return_data


# as this is ungrouped data, the "rate=False" option is not per bin
@pytest.mark.parametrize("unit,xstart,xwidth,rate,xlabel,ylabel,ynorm",
                         [('energy', 0.1, 0.125, True, 'Energy (keV)',
                           'Counts/sec/keV', 100.0 * 0.125),
                          ('energy', 0.1, 0.125, False, 'Energy (keV)',
                           'Counts', 1.0),
                          ('channel', 0.5, 1.0, True, 'Channel',
                           'Counts/sec/channel', 100.0),
                          ('channel', 0.5, 1.0, False, 'Channel',
                           'Counts', 1.0)])
def test_pha_bkg_dataplot(setup_pha, unit, xstart, xwidth, rate, xlabel,
                          ylabel, ynorm):
    """Test out BkgDataPlot"""

    bkg = setup_pha('example-bkg', exposure=100.0)
    bkg.units = unit
    bkg.rate = rate

    bplot = BkgDataPlot()
    bplot.prepare(bkg)

    assert bplot.xlabel == xlabel
    assert bplot.ylabel == ylabel
    assert bplot.title == 'example-bkg'

    nchans = 10
    xmid = xstart + xwidth / 2 + np.arange(nchans) * xwidth

    assert bplot.x == pytest.approx(xmid)
    assert bplot.y == pytest.approx(np.ones(nchans) / ynorm)
    assert bplot.xerr == pytest.approx(np.ones(nchans) * xwidth)
    assert bplot.yerr is None


@pytest.mark.parametrize("unit", ["channel", "energy", "wave"])
def test_pha_bkg_dataplot_areascal(setup_pha, unit):
    """The areascal is applied for DataPHA obects but not backscal.

    This is really a check for the DataPHA class, but the behavior
    is important for the rescale behavior, so add an explicit check
    here. At present it is assumed that plots of DataPHA data apply
    the area scaling to the y axis, but do not apply the backscale
    value (since that requires both a source and background dataset
    to be meaningful).

    See Also
    --------
    test_pha_bkg_dataplot_exposure
    """

    pha0 = setup_pha('example-bkg', areascal=1.0, exposure=100, backscal=1e-6)
    pha0.units = unit
    pha0.rate = False

    bplot0 = BkgDataPlot()
    bplot0.prepare(pha0)

    pha = setup_pha('example-bkg', areascal=0.4, exposure=200, backscal=1e-5)
    pha.units = unit
    pha.rate = False

    bplot = BkgDataPlot()
    bplot.prepare(pha)

    # since rate is False the different exposure values are not relevant
    # here: see test_pha_bkg_dataplot_exposure
    assert bplot.y == pytest.approx(bplot0.y / 0.4)


@pytest.mark.parametrize("unit", ["channel", "energy", "wave"])
def test_pha_bkg_dataplot_exposure(setup_pha, unit):
    """Is exposure handled correctly?

    This is really a check for the DataPHA class, but the behavior
    is important for the rescale behavior, so add an explicit check
    here.

    See Also
    --------
    test_pha_bkg_dataplot_areascal
    """

    pha0 = setup_pha('example-bkg', areascal=1.0, exposure=100, backscal=1e-6)
    pha0.units = unit
    pha0.rate = True

    bplot0 = BkgDataPlot()
    bplot0.prepare(pha0)

    pha = setup_pha('example-bkg', areascal=1.0, exposure=200, backscal=1e-5)
    pha.units = unit
    pha.rate = True

    bplot = BkgDataPlot()
    bplot.prepare(pha)

    assert bplot.y == pytest.approx(bplot0.y / 2.0)


def test_pha_bkg_dataplot_rescale_vals(setup_pha):
    """Test out BkgDataPlot using the rescale option: compare to Sherpa

    This tests that the background scaling matches that done by the
    DataPHA class. The test_pha_bkg_rescale test checks that other
    settings are as expected.
    """

    src = setup_pha('example-src', exposure=100.0, areascal=1.0, backscal=1e-6)
    src.units = 'channel'
    src.rate = False

    bkg = setup_pha('example-bkg', exposure=500.0, areascal=0.8, backscal=2e-6)
    bkg.units = 'channel'
    bkg.rate = False

    bplot = BkgDataPlot()
    bplot.prepare(bkg, source=src)

    # by using ungrouped data, with channel analysis and displaying
    # counts, not a rate, the y value from the plot should match
    # the background returned by sum_background_data.
    #
    src.set_background(bkg)
    expected = src.sum_background_data()

    assert bplot.y == pytest.approx(expected)


@pytest.mark.parametrize("unit,xstart,xwidth,rate,xlabel,ylabel,norm",
                         [('energy', 0.1, 0.125, True, 'Energy (keV)',
                           'Counts/sec/keV', 500.0 * 0.125),
                          ('energy', 0.1, 0.125, False, 'Energy (keV)',
                           'Counts', 1.0),
                          ('channel', 0.5, 1.0, True, 'Channel',
                           'Counts/sec/channel', 500.0),
                          ('channel', 0.5, 1.0, False, 'Channel',
                           'Counts', 1.0)])
def test_pha_bkg_dataplot_rescale(setup_pha, unit, xstart, xwidth, rate,
                                  xlabel, ylabel, norm):
    """Test out BkgDataPlot using the rescale option.

    Here we use the same data for source and background, but adjust the
    scale factors to ensure that some scaling has happened.
    """

    src = setup_pha('example-src', exposure=100.0, areascal=1.0, backscal=1e-6)
    src.units = unit
    src.rate = rate

    bkg = setup_pha('example-bkg', exposure=500.0, areascal=0.8, backscal=2e-6)
    bkg.units = unit
    bkg.rate = rate

    bplot = BkgDataPlot()
    bplot.prepare(bkg, source=src)

    assert bplot.xlabel == xlabel
    assert bplot.ylabel == (ylabel + ' (scaled)')
    assert bplot.title == 'example-bkg'

    nchans = 10
    xmid = xstart + xwidth / 2 + np.arange(nchans) * xwidth

    assert bplot.x == pytest.approx(xmid)

    # the scaling factors could be included in norm, but leave that value
    # to be the same as the "non-rescaled" test, so add in the areascale,
    # backscale, and exposure-time ratios here where appropriate.
    #
    src_factors = 1.0 * 1e-6
    bkg_factors = 0.8 * 2e-6
    if not rate:
        src_factors *= 100.0
        bkg_factors *= 500.0

    expected = np.ones(nchans) * src_factors / (bkg_factors * norm)

    assert bplot.y == pytest.approx(expected)
    assert bplot.xerr == pytest.approx(np.ones(nchans) * xwidth)
    assert bplot.yerr is None


def test_pha_bkg_dataplot_nowarn(caplog, setup_pha):
    """Check we don't get warnings if the missing value is "unused".

    areascal is handled elsewhere, and exposure is only used when
    the rate flag is not set
    """

    bplot = BkgDataPlot()
    src = setup_pha('a', exposure=100.0, areascal=1.0, backscal=1e-6)
    bkg = setup_pha('b', exposure=None, areascal=None, backscal=2e-6)
    bplot.prepare(bkg, source=src)

    assert len(caplog.records) == 0


@pytest.mark.parametrize('change_bkg', [False, True])
@pytest.mark.parametrize('label,kw',
                         [('backscale', 'backscal'),
                          ('exposure', 'exposure')])
def test_pha_bkg_dataplot_warn(caplog, setup_pha, change_bkg, label, kw):
    "Do we get warnings when a required attribute is not set?"

    assert len(caplog.records) == 0

    kwargs = {'exposure': 1.0, 'backscal': 1.0}
    if change_bkg:
        srctype = 'background'
        src = setup_pha('example-src', **kwargs)
        del kwargs[kw]
        bkg = setup_pha('example-bkg', **kwargs)
    else:
        srctype = 'source'
        bkg = setup_pha('example-bkg', **kwargs)
        del kwargs[kw]
        src = setup_pha('example-src', **kwargs)

    src.rate = False
    bkg.rate = False

    bplot = BkgDataPlot()
    bplot.prepare(bkg, source=src)

    assert len(caplog.records) == 1
    log_name, log_level, message = caplog.record_tuples[0]
    assert log_name == 'sherpa.astro.plot'
    assert log_level == logging.WARNING

    wmsg = 'The {} value is missing from the '.format(label) + \
           '{} dataset, so the correction is ignored.'.format(srctype)
    assert message == wmsg
