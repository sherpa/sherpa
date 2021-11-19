#
#  Copyright (C) 2016, 2018, 2019, 2020, 2021  Smithsonian Astrophysical Observatory
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

from sherpa.utils.testing import requires_pylab
from sherpa.data import Data1D, Data1DInt
from sherpa.models.basic import Const1D
from sherpa.stats import Chi2DataVar
from sherpa.plot import DataPlot, ModelHistogramPlot, \
    DelchiPlot, RatioPlot, ResidPlot


@requires_pylab
@pytest.mark.parametrize("plottype", [DelchiPlot, RatioPlot, ResidPlot])
def test_ignore_ylog_prefs(plottype):
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


@requires_pylab
@pytest.mark.parametrize("plottype", [DelchiPlot, RatioPlot, ResidPlot])
def test_ignore_ylog_kwarg(plottype):
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


@requires_pylab
def test_warning_dataplot_linecolor(caplog):
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
    assert msg == 'The linecolor attribute, set to mousey, is unused.'


@requires_pylab
def test_warning_plot_hist_linecolor(caplog):
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
