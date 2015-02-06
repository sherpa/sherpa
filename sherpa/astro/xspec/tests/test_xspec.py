#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
import numpy
import sherpa.astro.xspec as xs
from sherpa.astro import ui
from sherpa.utils import SherpaTestCase, needs_data

import logging
error = logging.getLogger(__name__).error

def is_proper_subclass(obj, cls):
    if type(cls) is not tuple:
        cls = (cls,)
    if obj in cls:
        return False
    return issubclass(obj, cls)


class test_xspec(SherpaTestCase):

    def test_create_model_instances(self):
        count = 0

        for cls in dir(xs):
            if not cls.startswith('XS'):
                continue

            cls = getattr(xs, cls)

            if is_proper_subclass(cls, (xs.XSAdditiveModel,
                                        xs.XSMultiplicativeModel)):
                m = cls()
                count += 1

        self.assertEqual(count, 164)

    def test_evaluate_model(self):
        m = xs.XSbbody()
        out = m([1,2,3,4])
        if m.calc.__name__.startswith('C_'):
            otype = numpy.float64
        else:
            otype = numpy.float32
        self.assert_(out.dtype.type is otype)
        self.assertEqual(int(numpy.flatnonzero(out == 0.0)), 3)


    def test_xspec_models(self):
        models = [model for model in dir(xs) if model[:2] == 'XS']
        models.remove('XSModel')
        models.remove('XSMultiplicativeModel')
        models.remove('XSAdditiveModel')
        models.remove('XSTableModel')

        xx = numpy.arange(0.1, 11.01, 0.01, dtype=float)
        xlo = numpy.array(xx[:-1])
        xhi = numpy.array(xx[1:])
        for model in models:
            cls = getattr(xs, model)
            foo = cls('foo')
            vals = foo(xlo,xhi)
            try:
                self.assert_(not numpy.isnan(vals).any() and
                             not numpy.isinf(vals).any() )
            except AssertionError:
                error('XS%s model evaluation failed' % model)
                raise

    @needs_data
    def test_set_analysis_wave_fabrizio(self):
        rmf = self.datadir + '/ciao4.3/fabrizio/Data/3c273.rmf'
        arf = self.datadir + '/ciao4.3/fabrizio/Data/3c273.arf'

        ui.set_model("fabrizio", "xspowerlaw.p1")
        ui.fake_pha("fabrizio", arf, rmf, 10000)

        model = ui.get_model("fabrizio")
        bare_model, _ = ui._session._get_model_status("fabrizio")
        y = bare_model.calc([1,1], model.xlo, model.xhi)
        y_m = numpy.mean(y)

        ui.set_analysis("fabrizio","wave")

        model2 = ui.get_model("fabrizio")
        bare_model2, _ = ui._session._get_model_status("fabrizio")
        y2 = bare_model2.calc([1,1], model2.xlo, model2.xhi)
        y2_m = numpy.mean(y2)

        self.assertAlmostEqual(y_m, y2_m)

    def test_xsxset_get(self):
	# TEST CASE #1 Case insentitive keys
	xs.set_xsxset('fooBar', 'somevalue')
	self.assertEqual('somevalue', xs.get_xsxset('Foobar'))


if __name__ == '__main__':

    from sherpa.utils import SherpaTest
    SherpaTest(xs).test()
