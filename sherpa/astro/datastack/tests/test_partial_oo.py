#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_

from sherpa.utils import SherpaTestCase
import numpy as np
from sherpa.astro.datastack import *
logger = logging.getLogger('sherpa')
logger.setLevel(logging.ERROR)


class test_partial_oo(SherpaTestCase):

    def setUp(self):
        self.ds = DataStack()
        set_template_id("__ID")
        clear_stack()
        ui.clean()

    def tearDown(self):
        self.ds.clear_stack()
        clear_stack()
        set_template_id("__ID")
        ui.clean()

    def test_case_1(self):
        x1 = np.arange(50)+100
        y1 = 2*(3*x1**2 + x1)
        x2 = np.arange(50)
        y2 = 2*(x2**2 + 3*x2)
        x3 = np.arange(50)+200
        y3 = 2*(x3**2+3*x3)

        ds = self.ds

        load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

        set_source(ds, 'const1d.const * polynom1d.poly__ID')

        poly1 = ui._session._get_model_component('poly1')
        poly2 = ui._session._get_model_component('poly2')
        poly3 = ui._session._get_model_component('poly3')
        const = ui._session._get_model_component('const')

        freeze(ds, 'poly')
        thaw(ds, 'poly.c0')
        thaw(ds, 'poly.c1')
        thaw(ds, 'poly.c2')
        thaw(ds, 'const')

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

        set_par(ds, 'poly.c1', 0.45)

        assert poly1.c1.val == 0.45
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par(ds[1], 'poly.c1', 0.1)

        assert poly1.c1.val == 0.1
        assert poly2.c1.val == 0.45
        assert poly3.c1.val == 0.45

        set_par(ds, 'const.c0', 2)

        assert const.c0.val == 2

        set_par(ds, 'const.integrate', False)
        freeze(ds, 'const.c0')

        vals = get_par(ds, 'poly.c1.val')
        assert ([0.1, 0.45, 0.45] == vals).all()

        pars = get_par(ds, 'const.c0')

        fit(ds)

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

        load_arrays(ds, [[x1, y1], [x2, y2], [x3, y3]])

        set_template_id("foo")

        set_source(ds, 'const1d.constfoo * polynom1d.polyfoo')

        const1 = ui._session._get_model_component('const1')
        const2 = ui._session._get_model_component('const2')
        const3 = ui._session._get_model_component('const3')

        link(ds[2,3], 'const.c0')

        set_par(ds[2], 'const.c0', 3)
        set_par(ds[1], 'const.c0', 7)

        freeze(ds[1], 'const.c0')

        assert const2.c0.frozen is False

        fit(ds)

        assert const2.c0.val == const3.c0.val
        assert const3.c0._link is const2.c0

        unlink(ds, "const.c0")

        assert const3.c0._link is not const2.c0

if __name__ == '__main__':

    from sherpa.utils import SherpaTest
    import sherpa.astro
    SherpaTest(sherpa.astro.datastack).test()

