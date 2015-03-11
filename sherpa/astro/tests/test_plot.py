# 
#  Copyright (C) 2007  Smithsonian Astrophysical Observatory
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

from sherpa.utils import SherpaTestCase, needs_data
from sherpa.astro.data import DataPHA
from sherpa.all import *
from sherpa.astro.all import *
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
