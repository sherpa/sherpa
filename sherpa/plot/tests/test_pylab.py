#
#  Copyright (C) 2016, 2018-2021, 2024-2025
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

from sherpa.data import Data1D, Data1DInt
from sherpa.models.basic import Const1D
from sherpa.stats import Chi2DataVar
from sherpa.plot import DataPlot, Histogram, ModelHistogramPlot, \
    DelchiPlot, RatioPlot, ResidPlot, Contour


@pytest.mark.parametrize("plottype", [DelchiPlot, RatioPlot, ResidPlot])
def test_ignore_ylog_prefs(plottype, requires_pylab):
    """Do the "residual" style plots ignore the ylog preference setting?"""

    from matplotlib import pyplot as plt

    data = Data1D('tst', np.asarray([1, 2, 3]), np.asarray([10, 12, 10.5]))
    mdl = Const1D('tst-model')
    mdl.c0 = 11.1

    plot = plottype()
    plot.plot_prefs['xlog'] = True
    plot.plot_prefs['ylog'] = True
    plot.prepare(data, mdl, stat=Chi2DataVar())
    plot.plot()

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'linear'


@pytest.mark.parametrize("plottype", [DelchiPlot, RatioPlot, ResidPlot])
def test_ignore_ylog_kwarg(plottype, requires_pylab):
    """Do the "residual" style plots ignore the ylog keyword argument?"""

    from matplotlib import pyplot as plt

    data = Data1D('tst', np.asarray([1, 2, 3]), np.asarray([10, 12, 10.5]))
    mdl = Const1D('tst-model')
    mdl.c0 = 11.1

    plot = plottype()
    plot.prepare(data, mdl, stat=Chi2DataVar())
    plot.plot(xlog=True, ylog=True)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = plt.gca()
    assert ax.get_xscale() == 'log'
    assert ax.get_yscale() == 'linear'


def test_warning_dataplot_linecolor(caplog, requires_pylab):
    """We get a warning when using linecolor: DataPlot"""

    data = Data1D('tst', np.asarray([1, 2, 3]), np.asarray([10, 12, 10.5]))
    plot = DataPlot()
    plot.prepare(data, stat=None)
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.plot(linecolor='mousey')

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.plot.pylab_backend'
    assert lvl == logging.WARNING
    assert msg == 'The linecolor attribute (mousey) is unused.'


def test_warning_plot_hist_linecolor(caplog, requires_pylab):
    """We get a warning when using linecolor: ModelHistogramPlot"""

    data = Data1DInt('tst', np.asarray([1, 2, 3]),
                     np.array([2, 2.5, 4]),
                     np.asarray([10, 12, 10.5]))
    mdl = Const1D()
    plot = ModelHistogramPlot()
    plot.prepare(data, mdl, stat=None)
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.plot(linecolor='mousey')

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.plot.pylab_backend'
    assert lvl == logging.WARNING
    assert msg == 'The linecolor attribute (mousey) is unused.'


def test_pylab_specific_option_DataPlot(requires_pylab):
    """Test a matplotlib-specific option (markerfacecolor)"""

    from matplotlib import pyplot as plt

    data = Data1D('tst', np.asarray([1, 2, 3]), np.asarray([10, 12, 10.5]))
    plot = DataPlot()
    plot.prepare(data, stat=None)
    plot.plot(marker='o', markerfacecolor='red')

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = plt.gca()
    lines = ax.get_lines()
    assert len(lines) == 1
    line = lines[0]
    mfc = line.get_markerfacecolor()
    # Depending on matplotlib version this can be a string or an RGBA tuple
    if isinstance(mfc, str):
        assert mfc == 'red'
    else:
        assert np.allclose(mfc, (1.0, 0.0, 0.0, 1.0))

    assert line.get_markeredgewidth() == pytest.approx(1.0)


def test_pylab_specific_histogram_overplot(requires_pylab):
    """Test a matplotlib-specific option (markeredgewidth) when overplotting a
    histogram

    We are not testing all the properties. The fact that it doesn't error
    out means it was passed through.
    `test_backend_option_does_not_exist` ensures that.
    """
    from matplotlib import pyplot as plt

    edges = np.asarray([-10, -5, 5, 12, 17, 20, 30, 56, 60])
    y = np.asarray([28, 62, 17, 4, 2, 4, 125, 55])

    d = Data1DInt('example histogram', edges[:-1], edges[1:], y)
    dplot = DataPlot()
    dplot.prepare(d)
    dplot.plot(markeredgewidth=2, markevery=2, fillstyle='none')

    hplot = Histogram()
    hplot.overplot(d.xlo, d.xhi, d.y, sketch_params=5)

    fig = plt.gcf()
    assert len(fig.axes) == 1

    ax = plt.gca()
    lines = ax.get_lines()
    assert len(lines) == 3
    line = lines[0]
    assert line.get_markeredgewidth() == 2


def test_backend_option_does_not_exist(requires_pylab):
    """Test that we get an error when a backend-specific option doesn't exist

    We use ModelHistogramPlot for a change.
    """
    edges = np.asarray([-10, -5, 5, 12, 17, 20, 30, 56, 60])
    y = np.asarray([28, 62, 17, 4, 2, 4, 125, 55])

    d = Data1DInt('example histogram', edges[:-1], edges[1:], y)
    mdl = Const1D()
    plot = ModelHistogramPlot()
    plot.prepare(d, mdl, stat=None)

    with pytest.raises(AttributeError, match="got an unexpected keyword argument 'mousey'"):
        plot.plot(mousey='green')


def test_backend_specific_option_contour(requires_pylab, caplog):
    """Test matplotlib-specific options (vmin, vmax) for contour plots"""

    from matplotlib import pyplot as plt

    x1, x0 = np.mgrid[20:30, 5:20]
    y = np.sqrt((x0 - 10)**2 + (x1 - 31)**2)
    contour = Contour()
    contour.contour(x0, x1, y, vmax=2, vmin=2)
    # Not testing the actual colors of the contour levels,
    # but no error was raised, so it worked.

    assert len(caplog.record_tuples) == 0