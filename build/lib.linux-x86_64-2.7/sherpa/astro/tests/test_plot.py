#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
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
        self.src = PowLaw1D('p1')*XSphabs('abs1')

    def tearDown(self):
        logger.setLevel(self.old_level)

    def test_sourceplot(self):
        sp = SourcePlot()
        sp.prepare(self.data, self.src)
        #sp.plot()
