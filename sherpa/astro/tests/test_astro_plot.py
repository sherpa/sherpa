#
#  Copyright (C) 2007, 2015, 2018  Smithsonian Astrophysical Observatory
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
import numpy
from sherpa.utils.testing import SherpaTestCase, requires_data, requires_fits
from sherpa.astro.data import DataPHA
from sherpa.astro.plot import DataPlot, SourcePlot
from sherpa.models.basic import PowLaw1D
from sherpa.astro.optical import AbsorptionGaussian
from sherpa import stats

import pytest


import logging
logger = logging.getLogger('sherpa')


class test_plot(SherpaTestCase):

    def setUp(self):
        self.old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        self.data = DataPHA('', numpy.arange(10), numpy.ones(10),
                            bin_lo = numpy.arange(0.1, 10, 0.1),
                            bin_hi = numpy.arange(0.2, 10.1, 0.1) )
        self.data.units="energy"
        self.src = PowLaw1D('p1')*AbsorptionGaussian('abs1')

    def tearDown(self):
        logger.setLevel(self.old_level)

    def test_sourceplot(self):
        sp = SourcePlot()
        sp.prepare(self.data, self.src)
        #sp.plot()


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
