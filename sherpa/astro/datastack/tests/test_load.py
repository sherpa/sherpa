#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import os
import tempfile
import numpy as np
from sherpa.astro.ui import *
from sherpa.astro.datastack import *
logger = logging.getLogger('sherpa')
logger.setLevel(logging.ERROR)


class test_load(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)
        self.create_files()

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        os.remove(self.lisname)
        os.remove(self.name1)
        os.remove(self.name2)
        set_stack_verbose(False)



    def test_case_1(self):
        load_ascii("@{}".format(self.lisname))
        assert len(ui._session._data) == 2
        assert len(DATASTACK.datasets) == 2

        load_data("@{}".format(self.lisname))
        assert len(ui._session._data) == 4
        assert len(DATASTACK.datasets) == 4

        load_data("@"+"/".join((self._this_dir, 'data', 'pha.lis')))
        assert len(ui._session._data) == 6
        assert len(DATASTACK.datasets) == 6


    def create_files(self):
        fd1, self.name1 = tempfile.mkstemp()
        fd2, self.name2 = tempfile.mkstemp()

        fdlis, self.lisname = tempfile.mkstemp()

        ascii1 = open(self.name1, 'w')
        ascii2 = open(self.name2, 'w')
        lis = open(self.lisname, 'w')

        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)

        for x, y in zip(x1, y1):
            ascii1.write("{} {}\n".format(x, y))
        ascii1.close()
        os.close(fd1)

        for x, y in zip(x2, y2):
            ascii2.write("{} {}\n".format(x, y))
        ascii2.close()
        os.close(fd2)

        lis.write("{}\n".format(self.name1))
        lis.write("{}\n".format(self.name2))
        lis.close()
        os.close(fdlis)


