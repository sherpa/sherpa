#
#  Copyright (C) 2016, 2018, 2019  Smithsonian Astrophysical Observatory
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

from matplotlib import pyplot as plt

from sherpa.utils.testing import requires_pylab
from sherpa.data import Data1D
from sherpa.models.basic import Const1D
from sherpa.stats import Chi2DataVar
from sherpa.plot import DelchiPlot, RatioPlot, ResidPlot

import pytest


@requires_pylab
def test_axes_default():
    """Have we default values for all the axis settings?"""

    from sherpa.plot.pylab_backend import _errorbar_defaults

    # assert all needed defaults have been found, but do not
    # check the actual values. Split into multiple lines so that
    # if there's a failure we can see which is False.
    #
    fields = ('ecolor', 'capsize', 'barsabove')
    defs = [e in _errorbar_defaults for e in fields]
    assert all(defs)


@requires_pylab
@pytest.mark.parametrize("plottype", [DelchiPlot, RatioPlot, ResidPlot])
def test_ignore_ylog_prefs(plottype):
    """Do the "residual" style plots ignore the ylog preference setting?"""

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
