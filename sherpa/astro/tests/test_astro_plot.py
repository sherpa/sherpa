#
#  Copyright (C) 2007, 2015, 2018, 2019, 2020, 2021
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

import logging

import numpy as np

import pytest

from sherpa.utils.testing import requires_data, requires_fits

from sherpa.astro.data import DataARF, DataPHA
from sherpa.astro.instrument import create_delta_rmf
from sherpa.astro.plot import SourcePlot, \
    DataPHAPlot, ModelPHAHistogram, OrderPlot, \
    EnergyFluxHistogram, PhotonFluxHistogram,  _check_hist_bins
from sherpa.astro import plot as aplot
from sherpa.astro import hc
from sherpa.data import Data1D
from sherpa.models.basic import Const1D, Gauss1D, Polynom1D, PowLaw1D
from sherpa import stats
from sherpa.utils.err import IOErr, PlotErr


def check_sourceplot_energy(sp, factor=0):
    """Check for test_sourceplot/test_sourceplot_channel

    Note that the rate setting does not change these tests.
    """

    assert sp.xlabel == 'Energy (keV)'

    # the following depends on the backend
    # assert sp.ylabel == 'f(E)  Photons/sec/cm$^2$/keV'
    if factor == 0:
        assert sp.ylabel.startswith('f(E)  Photons/sec/cm')
        assert sp.ylabel.endswith('/keV ')
    elif factor == 1:
        assert sp.ylabel.startswith('E f(E)  Photons/sec/cm')
        assert sp.ylabel.find('/keV ') == -1
    elif factor == 2:
        # This says E^2 f(E) ... but the exact format depends on
        # the back end
        assert sp.ylabel.startswith('E')
        assert sp.ylabel.find('^2') != -1
        assert sp.ylabel.find(' f(E)  Photons/sec/cm') != -1
        assert sp.ylabel.find('/keV ') == -1
    else:
        raise RuntimeError("unsupported factor")

    assert sp.title == 'Source Model of '

    bins = np.arange(0.1, 10.1, 0.1)
    assert sp.xlo == pytest.approx(bins[:-1])
    assert sp.xhi == pytest.approx(bins[1:])

    # The check of the values is just to check that things are going
    # as expected, rather than a test from first principles.
    #
    yexp = np.asarray([9998.300959, 9997.9935414, 9997.63869638,
                       9997.23070803, 9996.76346026, 9996.2304594, 9995.62486769,
                       9994.9395488, 9994.16712626, 9993.30005567, 9992.33071112,
                       9991.25148615, 9990.05490914, 9988.73377265, 9987.28127577,
                       9985.69117819, 9983.95796409, 9982.07701351, 9980.04477839,
                       9977.85895997, 9975.51868378, 9973.02466824, 9970.37938236,
                       9967.58718801, 9964.65446216, 9961.58969431, 9958.40355487,
                       9955.10893021, 9951.72092088, 9948.25679985, 9944.73592864,
                       9941.17962987, 9937.6110159, 9934.05477423, 9930.53691143,
                       9927.08445867, 9923.7251427, 9920.48702778, 9917.39813433,
                       9914.48604166, 9911.77748224, 9909.29793588, 9907.07123238,
                       9905.11917118, 9903.4611666, 9902.11392658, 9901.09117249,
                       9900.40340655, 9900.05773225, 9900.05773225, 9900.40340655,
                       9901.09117249, 9902.11392658, 9903.4611666, 9905.11917118,
                       9907.07123238, 9909.29793588, 9911.77748224, 9914.48604166,
                       9917.39813433, 9920.48702778, 9923.7251427, 9927.08445867,
                       9930.53691143, 9934.05477423, 9937.6110159, 9941.17962987,
                       9944.73592864, 9948.25679985, 9951.72092088, 9955.10893021,
                       9958.40355487, 9961.58969431, 9964.65446216, 9967.58718801,
                       9970.37938236, 9973.02466824, 9975.51868378, 9977.85895997,
                       9980.04477839, 9982.07701351, 9983.95796409, 9985.69117819,
                       9987.28127577, 9988.73377265, 9990.05490914, 9991.25148615,
                       9992.33071112, 9993.30005567, 9994.16712626, 9994.9395488,
                       9995.62486769, 9996.2304594, 9996.76346026, 9997.23070803,
                       9997.63869638, 9997.9935414, 9998.300959, 9998.5662520])

    if factor == 1:
        yexp *= 0.1
    elif factor == 2:
        yexp *= 0.01

    assert sp.y == pytest.approx(yexp)


def check_sourceplot_wavelength(sp, factor=0):
    """Check for test_sourceplot_wavelength.

    See check_sourceplot_energy.
    """

    assert sp.xlabel == 'Wavelength (Angstrom)'

    if factor == 0:
        assert sp.ylabel.startswith('f(lambda)  Photons/sec/cm')
        assert sp.ylabel.endswith('/Angstrom ')
    elif factor == 1:
        assert sp.ylabel.startswith('lambda f(lambda)  Photons/sec/cm')
        assert sp.ylabel.find('/Angstrom ') == -1
    elif factor == 2:
        # This says lambda^2 f(lambda) ... but the exact format depends on
        # the back end
        assert sp.ylabel.startswith('lambda')
        assert sp.ylabel.find('^2') != -1
        assert sp.ylabel.find(' f(lambda)  Photons/sec/cm') != -1
        assert sp.ylabel.find('/Angstrom ') == -1
    else:
        raise RuntimeError("unsupported factor")

    assert sp.title == 'Source Model of '

    # Use the code from DataPHA to convert from keV to Angstroms.
    #
    ebins = np.arange(0.1, 10.1, 0.1)
    lbins = hc / ebins
    assert sp.xlo == pytest.approx(lbins[:-1])
    assert sp.xhi == pytest.approx(lbins[1:])

    # Although the same model as check_sourceplot_energy the
    # interpretation of the model params means that there is no simple
    # conversion of the values used in that case, so we have a
    # different set of expected values.
    #
    yexp = np.asarray([10000., 10000., 10000., 10000., 10000., 10000.,
                       9999.99999864, 9999.99949354, 9999.97224156,
                       9999.55447151, 9996.86451659, 9987.51382359,
                       9966.695068 , 9933.27750675, 9891.33791156,
                       9847.93710785, 9809.78226046, 9781.16179805,
                       9763.58920097, 9756.46908383, 9758.02782873,
                       9766.08494609, 9778.54263095, 9793.6274431,
                       9809.96473056, 9826.55998249, 9842.73898966,
                       9858.07743665, 9872.33534744, 9885.40252063,
                       9897.25609947, 9907.92910978, 9917.48798344,
                       9926.01701814, 9933.60797703, 9940.3533812,
                       9946.34238774, 9951.65843565, 9956.37807112,
                       9960.57053782, 9964.29784539, 9967.61512156,
                       9970.57111803, 9973.20878565, 9975.56586528,
                       9977.67546184, 9979.56658298, 9981.26463295,
                       9982.79185821, 9984.16774487, 9985.40937017,
                       9986.5317114, 9987.54791625, 9988.46953851,
                       9989.30674326, 9990.06848503, 9990.76266241,
                       9991.39625226, 9991.97542597, 9992.50565038,
                       9992.99177532, 9993.43810959, 9993.84848703,
                       9994.22632392, 9994.57466896, 9994.89624683,
                       9995.19349616, 9995.46860277, 9995.72352864,
                       9995.96003731, 9996.17971622, 9996.3839962,
                       9996.57416871, 9996.75140095, 9996.9167492,
                       9997.07117061, 9997.21553359, 9997.35062702,
                       9997.47716842, 9997.59581118, 9997.70715101,
                       9997.81173162, 9997.9100499 , 9998.00256036,
                       9998.0896793 , 9998.17178835, 9998.24923778,
                       9998.32234937, 9998.39141906, 9998.45671928,
                       9998.51850107, 9998.57699596, 9998.63241771,
                       9998.68496386, 9998.73481709, 9998.78214653,
                       9998.82710888, 9998.86984944, 9998.91050307])

    if factor == 1:
        yexp *= (lbins[:-1] - lbins[1:])
    elif factor == 2:
        yexp *= (lbins[:-1] - lbins[1:])**2

    assert sp.y == pytest.approx(yexp)


def test_sourceplot(caplog):

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "energy"
    assert data.rate
    assert data.plot_fac == 0

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
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_energy(sp)


def test_sourceplot_counts(caplog):
    """test_sourceplot but when rate=False is chosen"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "energy"
    data.rate = False

    # Note that the model evaluation in done in Angstroms
    #
    m1 = Const1D('bgnd')
    m2 = Gauss1D('abs1')
    src = 100 * m1 * (1 - m2) * 10000

    m1.c0 = 0.01
    m2.pos = 5.0
    m2.fwhm = 4.0
    m2.ampl = 0.1

    sp = SourcePlot()
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_energy(sp)


@pytest.mark.parametrize("factor", [1, 2])
def test_sourceplot_facn(factor, caplog):
    """Change plot factor for test_sourceplot"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "energy"
    data.plot_fac = factor
    assert data.rate

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
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_energy(sp, factor=factor)


def test_sourceplot_channels(caplog):
    """Although we ask for channels we get energy units"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "channel"

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
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.astro.plot'
    assert lvl == logging.WARN
    assert msg == 'Channel space is unappropriate for the PHA unfolded source model,\nusing energy.'

    check_sourceplot_energy(sp)


def test_sourceplot_wavelength(caplog):
    """Check we get wavelength units"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "wave"

    # Note that the model evaluation in done in Angstroms
    #
    m1 = Const1D('bgnd')
    m2 = Gauss1D('abs1')
    src = 100 * m1 * (1 - m2) * 10000

    m1.c0 = 0.01
    m2.pos = 5.0
    m2.fwhm = 4.0
    m2.ampl = 0.1

    sp = SourcePlot()
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_wavelength(sp)


@pytest.mark.parametrize("factor", [1, 2])
def test_sourceplot_wavelength_facn(factor, caplog):
    """Change plot factor for test_sourceplot_wavelength"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "wavelength"
    data.plot_fac = factor
    assert data.rate

    m1 = Const1D('bgnd')
    m2 = Gauss1D('abs1')
    src = 100 * m1 * (1 - m2) * 10000

    m1.c0 = 0.01
    m2.pos = 5.0
    m2.fwhm = 4.0
    m2.ampl = 0.1

    sp = SourcePlot()
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_wavelength(sp, factor=factor)


def test_sourceplot_wavelength_counts(caplog):
    """test_sourceplot_wavelength but when rate=False is chosen"""

    bins = np.arange(0.1, 10.1, 0.1)
    data = DataPHA('', np.arange(10), np.ones(10),
                   bin_lo=bins[:-1].copy(),
                   bin_hi=bins[1:].copy())
    data.units = "wave"
    data.rate = False

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
    with caplog.at_level(logging.INFO, logger='sherpa'):
        sp.prepare(data, src)

    assert len(caplog.records) == 0
    check_sourceplot_wavelength(sp)


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
                          stats.WStat()])
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

    dplot = DataPHAPlot()
    dplot.prepare(pha, stat=stat)


@pytest.mark.parametrize("ptype",
                         [aplot.ModelPHAHistogram, aplot.SourcePlot])
def test_plot_fail_with_non_pha(ptype):
    """plots don't like Data1D objects"""

    x = np.arange(3)
    y = np.ones(3)
    d = Data1D('tst', x, y)

    m = Const1D()

    p = ptype()

    with pytest.raises(IOErr) as exc:
        p.prepare(d, m)

    assert str(exc.value) == "data set 'tst' does not contain a PHA spectrum"


@requires_data
@requires_fits
def test_dataphahistogram_prepare_wavelength(make_data_path):
    """Check we can use wavelength setting"""

    from sherpa.astro.io import read_pha

    # could fake a dataset but it's easier to use one
    infile = make_data_path('3c273.pi')
    pha = read_pha(infile)
    pha.name = 'my-name.pi'

    # Also check out the type='counts' option
    pha.set_analysis('wave', type='counts')
    pha.notice(3, 5)

    plot = aplot.DataPHAPlot()
    plot.prepare(pha)

    assert plot.xlabel == 'Wavelength (Angstrom)'
    assert plot.ylabel == 'Counts/Angstrom'
    assert plot.title == 'my-name.pi'

    # data is inverted
    assert plot.xlo[0] > plot.xlo[-1]

    # can we access the "pseudo" x attribute?
    assert plot.x[0] > plot.x[-1]

    assert np.all(plot.y > 0)

    # regression test
    yexp = np.asarray([39.44132393, 68.92072985, 64.39425057,
                       48.69954162, 41.17454851, 88.87982014,
                       74.70504415, 79.22094498, 94.57773635])
    assert plot.y == pytest.approx(yexp)


@requires_data
@requires_fits
def test_modelphahistogram_prepare_wavelength(make_data_path):
    """Check we can use wavelength setting"""

    from sherpa.astro.io import read_pha

    # could fake a dataset but it's easier to use one
    infile = make_data_path('3c273.pi')
    pha = read_pha(infile)
    pha.name = 'my-name.pi'

    pha.set_analysis('wave')
    pha.notice(3, 5)

    mdl = Const1D()

    # not bothered too much about the model (e.g. setting a response)
    #
    plot = aplot.ModelPHAHistogram()
    plot.prepare(pha, mdl)

    assert plot.xlabel == 'Wavelength (Angstrom)'
    assert plot.ylabel == 'Counts/sec/Angstrom'
    assert plot.title == 'Model'

    # data is inverted
    assert plot.xlo[0] > plot.xlo[-1]
    assert plot.xlo[0] > plot.xhi[0]
    assert np.all(plot.y > 0)
    assert plot.y.size == 9


@requires_data
@requires_fits
def test_sourceplot_prepare_wavelength(make_data_path):
    """Check we can use wavelength setting"""

    from sherpa.astro.io import read_pha

    # could fake a dataset but it's easier to use one
    infile = make_data_path('3c273.pi')
    pha = read_pha(infile)
    pha.name = 'my-name.pi'

    pha.set_analysis('wave')
    pha.notice(3, 5)

    mdl = Const1D()

    # not bothered too much about the model (e.g. setting a response)
    #
    plot = aplot.SourcePlot()
    plot.prepare(pha, mdl)

    assert plot.xlabel == 'Wavelength (Angstrom)'
    assert plot.ylabel.startswith('f(lambda)  Photons/sec/cm')
    assert plot.ylabel.find('^2') != -1
    assert plot.ylabel.endswith('/Angstrom ')
    assert plot.title == 'Source Model of my-name.pi'

    # data is inverted
    assert plot.xlo[0] > plot.xlo[-1]
    assert plot.xlo[0] > plot.xhi[0]
    assert np.all(plot.y > 0)
    assert plot.y.size == 1090


def test_pha_data_with_gaps_977():
    """If the lo/hi edges don't quite match what happens?

    Ideally it recognizes the data is the same (to float32
    "precision"). At the moment this is only done for PHA
    data and model plots, not for generic histogram plots.
    See issue #977
    """

    chans = np.arange(1, 6)
    vals = np.arange(1, 6)

    blo = np.asarray([100, 99, 98, 97, 96])
    bhi = np.asarray([101, 100.0000000001, 99, 98, 97])

    d = DataPHA('x', chans, vals, bin_lo=blo, bin_hi=bhi)
    d.set_analysis('wave')

    p = DataPHAPlot()
    p.prepare(d)

    assert p.y == pytest.approx([1, 2, 3, 4, 5])

    xlo = p.xlo
    xhi = p.xhi

    assert xlo.size == 5
    assert xlo[0] == pytest.approx(101)
    assert xhi[0] == pytest.approx(100)

    assert xlo[-1] == pytest.approx(97)
    assert xhi[-1] == pytest.approx(96)

    # This is an equality check, not with pytest.approx,
    # since this is enforced by the plot code. This fails
    # before #977 is fixed.
    #
    assert (xlo[1:] == xhi[:-1]).all()


def test_pha_model_with_gaps_977():
    """If the lo/hi edges don't quite match what happens?

    See test_pha_data_with_gaps_977.
    """

    chans = np.arange(1, 6)
    vals = np.arange(1, 6)

    blo = np.asarray([100, 99, 98, 97, 96])
    bhi = np.asarray([101, 100.0000000001, 99, 98, 97])

    d = DataPHA('x', chans, vals, bin_lo=blo, bin_hi=bhi)
    d.set_analysis('wave')

    mdl = Polynom1D()
    mdl.c0 = 0.1
    mdl.c1 = 1.1

    p = ModelPHAHistogram()
    p.prepare(d, mdl)

    assert p.y == pytest.approx([1.2, 2.3, 3.4, 4.5, 5.6])

    xlo = p.xlo
    xhi = p.xhi

    assert xlo.size == 5
    assert xlo[0] == pytest.approx(101)
    assert xhi[0] == pytest.approx(100)

    assert xlo[-1] == pytest.approx(97)
    assert xhi[-1] == pytest.approx(96)

    # This is an equality check, not with pytest.approx,
    # since this is enforced by the plot code. This fails
    # before #977 is fixed.
    #
    assert (xlo[1:] == xhi[:-1]).all()


@pytest.mark.parametrize("energy,cls",
                         [(True, EnergyFluxHistogram),
                          (False, PhotonFluxHistogram)])
def test_str_flux_histogram_empty(energy, cls):
    """Check str of an empty flux histogram"""

    obj = cls()
    out = str(obj).split('\n')

    assert out[0] == 'modelvals = None'
    assert out[1] == 'clipped = None'
    assert out[2] == 'flux = None'
    assert out[3] == 'xlo    = None'
    assert out[4] == 'xhi    = None'
    assert out[5] == 'y      = None'

    # the exact text depends on the plot backend
    if energy:
        assert out[6].startswith('xlabel = Energy flux ')
        assert out[8] == 'title  = Energy flux distribution'
    else:
        assert out[6].startswith('xlabel = Photon flux ')
        assert out[8] == 'title  = Photon flux distribution'

    assert out[7] == 'ylabel = Frequency'
    assert out[9].startswith('histo_prefs = ')

    assert len(out) == 10


@pytest.mark.parametrize("energy,cls",
                         [(True, EnergyFluxHistogram),
                          (False, PhotonFluxHistogram)])
def test_str_flux_histogram_full(energy, cls, old_numpy_printing):
    """Check str of a flux histogram"""

    # First column is flux, next two are pars, and the
    # last column is a clipped column.
    #
    args = np.asarray([[1.0, 0.1, 1.1, 1],
                       [1.5, 0.2, 1.1, 1],
                       [2.0, 1.0, 2.0, 0],
                       [0.5, 0.4, 0.9, 1]])

    obj = cls()
    obj.prepare(args, 3)

    out = str(obj).split('\n')
    print(out)

    # lines 1-4 are the modelvals array;assume they are
    # displayed correctly
    assert out[0] == 'modelvals = [[ 0.1, 1.1],'
    assert out[4] == 'clipped = [ 1., 1., 0., 1.]'
    assert out[5] == 'flux = [ 1. , 1.5, 2. , 0.5]'
    assert out[6] == 'xlo    = [ 0.5  , 0.875, 1.25 , 1.625]'
    assert out[7] == 'xhi    = [ 0.875, 1.25 , 1.625, 2.   ]'
    assert out[8] == 'y      = [ 1., 1., 1., 1.]'

    # the exact text depends on the plot backend
    if energy:
        assert out[9].startswith('xlabel = Energy flux ')
        assert out[11] == 'title  = Energy flux distribution'
    else:
        assert out[9].startswith('xlabel = Photon flux ')
        assert out[11] == 'title  = Photon flux distribution'

    assert out[10] == 'ylabel = Frequency'
    assert out[12].startswith('histo_prefs = ')

    assert len(out) == 13


# Based on sherpa/astro/ui/tests/test_astro_plot.py
#
_data_chan = np.linspace(1, 10, 10, dtype=np.int8)
_data_counts = np.asarray([0, 1, 2, 3, 4, 0, 1, 2, 3, 4],
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


def example_pha_data():
    """Create an example data set."""

    etime = 1201.0
    d = DataPHA('example', _data_chan.copy(),
                _data_counts.copy(),
                exposure=etime,
                backscal=0.2)

    a = DataARF('example-arf',
                _energies_lo.copy(),
                _energies_hi.copy(),
                _arf.copy(),
                exposure=etime)

    r = create_delta_rmf(_energies_lo.copy(),
                         _energies_hi.copy(),
                         e_min=_energies_lo.copy(),
                         e_max=_energies_hi.copy(),
                         offset=1, name='example-rmf')

    d.set_arf(a)
    d.set_rmf(r)
    return d


@pytest.mark.parametrize("orders", [None, [1]])
def test_orderplot_checks_colors_explicit_orders(orders):
    """Check we raise a PlotErr.

    We check both the default orders setting (None) and an explicit
    setting as, when writing this test, I found an error in the
    handling of the None case.

    """

    oplot = OrderPlot()

    pha = example_pha_data()
    model = PowLaw1D('example-pl')

    with pytest.raises(PlotErr) as pe:
        oplot.prepare(pha, model, orders=orders, colors=['green', 'orange'])

    assert str(pe.value) == "orders list length '1' does not match colors list length '2'"


def test_orderplot_check_title():
    """Is the title set?"""

    oplot = OrderPlot()

    pha = example_pha_data()
    model = PowLaw1D('example-pl')

    assert oplot.title is 'Model'
    oplot.prepare(pha, model)
    assert oplot.title == 'Model Orders [1]'


def test_orderplot_check_range():
    """Are the x/y values as expected?

    Note that this is not a proper test as the dataset only
    contains a single response. However, the code is not
    clear as what is meant to happen with multiple responses
    so do not spend too much time on this test here.
    """

    oplot = OrderPlot()

    pha = example_pha_data()
    model = Const1D('example-mdl')

    oplot.prepare(pha, model, orders=[2])
    assert len(oplot.xlo) == 1
    assert len(oplot.xhi) == 1
    assert len(oplot.y) == 1

    assert oplot.xlo[0] == pytest.approx(np.arange(1, 11))
    assert oplot.xhi[0] == pytest.approx(np.arange(2, 12))

    # constant model / exposure time
    yexp = np.ones(10) / 1201
    assert oplot.y[0] == pytest.approx(yexp)


def test_check_hist_bins():
    '''Should work for reversed order (which can happen for plotting
    wavelength) and xlo and xhi reversed (which can happen if pha/arf/rmf are
    not sorted in increasing energy)
    '''
    xlo = np.arange(5, dtype=float)
    xhi = xlo + 1
    xhi[2] = 2.9999999
    for x1, x2 in [(xlo, xhi), (xhi, xlo),
                   (xlo[::-1], xhi[::-1]), (xlo[::-1], xhi[::-1])]:
        out1, out2 = _check_hist_bins(x1.copy(), x2.copy())
        assert np.all(out1[1:] == out2[:-1])
