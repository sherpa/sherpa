# 
#  Copyright (C) 2007, 2015  Smithsonian Astrophysical Observatory
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
from sherpa.utils import SherpaTestCase
from sherpa.astro.data import DataPHA
from sherpa.astro.plot import SourcePlot
from sherpa.models.basic import PowLaw1D, Polynom1D
from sherpa.astro.optical import AbsorptionGaussian
from sherpa.astro import ui
from sherpa.astro.instrument import DataRMF, DataARF

try:
    from unittest import mock
except ImportError:
    import mock

from six.moves import reload_module

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


def test_bug_106(mock_chips):
    # We need to make sure after the monkey patching
    # classes and functions are up to date.
    reload_module(ui)

    # We only need few points.
    channels = 32

    def apply_rmf(*args):
        """
        A simple rmf that simply divides the signal by 2
        """
        return args[0]/2

    def apply_arf(*args):
        """
        A simple arf that simply divides the signal by 5
        """
        return args[0]/5


    # Fake energy bins in 1024 channels
    e_lo = numpy.linspace(1, 20, channels, endpoint=False)
    e_hi = numpy.linspace(e_lo[1], 20, channels)

    # No errors
    ui.set_stat('chi2')

    def get_indep():
        return e_lo, e_hi

    # RMF stub definitions
    rmf_dict = {
        'detchans': channels,
        'energ_lo': e_lo,
        'energ_hi': e_hi,
        'get_indep': get_indep,
        'e_min': e_lo,
        'e_max': e_hi,
        'apply_rmf': apply_rmf,
    }

    # ARF stub definitions
    arf_dict = {
        'energ_lo': e_lo,
        'energ_hi': e_hi,
        'get_indep': get_indep,
        'apply_arf': apply_arf,
    }

    rmf = mock.MagicMock(spec=DataRMF, **rmf_dict)
    arf = mock.MagicMock(spec=DataARF, **arf_dict)

    source = Polynom1D()
    ui.set_source(source)
    source.c1.thaw()
    source.c0 = 0
    source.c1 = 1

    # Note noise is None, so the source is perfectly deterministic
    ui.fake_pha(1, arf=arf, rmf=rmf, exposure=50000, grouped=False, noise=None)

    # This works. We expect the source emission to be unaltered.
    ui.plot_source()
    expected = (e_lo+e_hi)/2
    chips = mock_chips[1]
    actual = chips.add_histogram.call_args[0][2]
    numpy.testing.assert_array_equal(expected, actual)

    # This does not work. Same test as above, but with `plot('source', 1)`
    ui.plot('source', 1)
    actual = chips.add_curve.call_args[0][1]
    numpy.testing.assert_array_equal(expected, actual)

    ##  Comment out above test to see the following tests run

    # This works. Now we expect the model to be attenuated by a factor 10
    ui.plot_model()
    expected = (e_lo+e_hi)/20
    actual = chips.add_histogram.call_args[0][2]
    numpy.testing.assert_array_almost_equal(expected, actual)

    # This does not work. Same test as above, but with `plot('model', 1)`
    ui.plot('model', 1)
    actual = chips.add_curve.call_args[0][1]
    numpy.testing.assert_array_almost_equal(expected, actual)