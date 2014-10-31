#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import os
from sherpa.astro.datastack import *
logger = logging.getLogger('sherpa')


class test_design(SherpaTestCase):
    def setUp(self):
        clear_stack()
        ui.clean()
        set_template_id("ID")
        logger.setLevel(logging.ERROR)
        set_stack_verbosity(logging.ERROR)
        self._this_dir = os.path.dirname(sys.modules[self.__module__].__file__)

    def tearDown(self):
        clear_stack()
        ui.clean()
        set_template_id("__ID")



    def test_case_1(self):
        datadir = '/'.join((self._this_dir, 'data'))
        ls = '@'+'/'.join((datadir, '3c273.lis'))
        load_pha(ls, use_errors=True)

        assert 2 == len(DATASTACK.datasets)
        assert 2 == len(ui._session._data)

        load_pha("myid", '/'.join((datadir, "3c273.pi")))

        assert 2 == len(DATASTACK.datasets)
        assert 3 == len(ui._session._data)

        load_pha('/'.join((datadir, "3c273.pi")))

        assert 3 == len(DATASTACK.datasets)
        assert 4 == len(ui._session._data)

        load_pha([], '/'.join((datadir, "3c273.pi")))

        assert 4 == len(DATASTACK.datasets)
        assert 5 == len(ui._session._data)

        ds = DataStack()

        load_pha(ds, ls)

        assert 4 == len(DATASTACK.datasets)
        assert 7 == len(ui._session._data)
        assert 2 == len(ds.datasets)

        load_pha(ds, '/'.join((datadir, "3c273.pi")))

        assert 4 == len(DATASTACK.datasets)
        assert 8 == len(ui._session._data)
        assert 3 == len(ds.datasets)

        dids = DATASTACK.get_stack_ids()
        assert dids == [1,2,3,4]

        sids = ui._session._data.keys()
        assert sids == [1,2,3,4,5,6,7, "myid"]

        set_source([1,2], "powlaw1d.pID")
        set_source([3,4], "brokenpowerlaw.bpID")

        dsids = ds.get_stack_ids()
        assert dsids == [5,6,7]

        p1 = ui._session._model_components['p1']
        p2 = ui._session._model_components['p2']
        bp3 = ui._session._model_components['bp3']
        bp4 = ui._session._model_components['bp4']

        assert p1 is not None
        assert p2 is not None
        assert bp3 is not None
        assert bp4 is not None

        set_source(1, "polynom1d.poly1")
        set_source([2,3,4], "atten.attID")

        poly1 = ui._session._model_components['poly1']
        a2 = ui._session._model_components['att2']
        a3 = ui._session._model_components['att3']
        a4 = ui._session._model_components['att4']

        assert poly1 is not None
        assert a2 is not None
        assert a3 is not None
        assert a4 is not None

        clean()

        assert 0 == len(DATASTACK.datasets)
        assert 0 == len(ui._session._data)
        assert 3 == len(ds.datasets)