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


class test_query(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")
        set_stack_verbose(False)

    def test_case_1(self):
        load_pha('@'+'/'.join((self._this_dir, 'data', 'pha.lis')))

        f = query_by_header_keyword('INSTRUME', 'ACIS')

        assert f == [1,2]

        f = query_by_obsid('7867')

        assert f == [2]

        ds = DataStack()

        ds.load_pha('@'+'/'.join((self._this_dir, 'data', 'pha.lis')))

        f = ds.query_by_obsid('4938')

        assert f == [3]

        f = query_by_obsid(ds, '7867')

        assert f == [4]

        



