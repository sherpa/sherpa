#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import numpy as np
from sherpa.astro import datastack
from sherpa.astro import ui
import logging

logger = logging.getLogger('sherpa')
logger.setLevel(logging.ERROR)


class test_oo(SherpaTestCase):

    def setUp(self):
        datastack.clear_stack()
        datastack.set_template_id("__ID")
        ui.clean()
        self.ds = datastack.DataStack()

    def tearDown(self):
        self.ds.clear_stack()
        datastack.clear_stack()
        datastack.set_template_id("__ID")
        ui.clean()

    def test_case_1(self):
        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds = self.ds

        ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        ds.set_source('const1d.const * polynom1d.poly__ID')

        poly1 = ui._session._get_model_component('poly1')
        poly2 = ui._session._get_model_component('poly2')
        poly3 = ui._session._get_model_component('poly3')
        const = ui._session._get_model_component('const')

        ds.freeze('poly')
        ds.thaw('poly.c0')
        ds.thaw('poly.c1')
        ds.thaw('poly.c2')
        ds.thaw('const')

        assert poly1.c0.frozen is False
        assert poly1.c1.frozen is False
        assert poly1.c2.frozen is False
        assert poly1.c3.frozen is True

        assert poly2.c0.frozen is False
        assert poly2.c1.frozen is False
        assert poly2.c2.frozen is False
        assert poly2.c3.frozen is True

        assert poly3.c0.frozen is False
        assert poly3.c1.frozen is False
        assert poly3.c2.frozen is False
        assert poly3.c3.frozen is True

        ds.set_par('poly.c1', 0.45)

        assert poly1.c1.val == 0.45
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        ds[1].set_par('poly.c1', 0.1)

        assert poly1.c1.val == 0.1
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        ds.set_par('const.c0', 2)

        assert const.c0.val == 2

        ds.set_par('const.integrate', False)
        ds.freeze('const.c0')

        vals = ds.get_par('poly.c1.val')
        assert ([0.1, 0.45, 0.45] == vals).all()

        pars = ds.get_par('const.c0')

        ds.fit()

        assert round(poly1.c1.val) == 1
        assert round(poly1.c2.val) == 3
        assert round(poly2.c1.val) == 3
        assert round(poly2.c2.val) == 1
        assert round(poly3.c1.val) == 3
        assert round(poly3.c2.val) == 1

        ds.clear_stack()

        x1 = np.arange(50)+100
        y1 = 7*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds.load_arrays([[x1, y1], [x2, y2], [x3, y3]])

        datastack.set_template_id("foo")

        ds.set_source('const1d.constfoo * polynom1d.polyfoo')

        const1 = ui._session._get_model_component('const1')
        const2 = ui._session._get_model_component('const2')
        const3 = ui._session._get_model_component('const3')

        ds[2,3].link('const.c0')

        ds[2].set_par('const.c0', 3)
        ds[1].set_par('const.c0', 7)

        ds[1].freeze('const.c0')

        assert const2.c0.frozen is False

        ds.fit()

        assert const2.c0.val == const3.c0.val
        assert const3.c0._link is const2.c0

        ds.unlink("const.c0")

        assert const3.c0._link is not const2.c0

if __name__ == '__main__':

    from sherpa.utils import SherpaTest
    import sherpa.astro
    SherpaTest(sherpa.astro.datastack).test()

