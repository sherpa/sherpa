#_PYTHON_INSERT_SAO_COPYRIGHT_HERE_(2007)_
#_PYTHON_INSERT_GPL_LICENSE_HERE_
from numpy import arange
import sherpa.astro.models as models
from sherpa.utils import SherpaFloat, SherpaTestCase
from sherpa.models.model import ArithmeticModel


class test_models(SherpaTestCase):

    def test_create_and_evaluate(self):
        x = arange(1.0, 5.0)
        count = 0

        for cls in dir(models):
            clsobj = getattr(models, cls)

            if ((not isinstance(clsobj, type)) or
                (not issubclass(clsobj, ArithmeticModel)) or
                (clsobj is ArithmeticModel)):
                continue

            # These have a very different interface than the others
            if cls in ('JDPileup', 'MultiResponseSumModel'):
                continue

            m = clsobj()
            self.assertEqual(type(m).__name__.lower(), m.name)
            count += 1

            if m.name == 'linebroad':
                m.vsini = 1e6

            try:
                if m.name.count('2d') or (m.name == 'hubblereynolds'):
                    pt_out  = m(x, x)
                    int_out = m(x, x, x, x)
                else:
                    pt_out  = m(x)
                    int_out = m(x, x)
            except ValueError:
                self.fail("evaluation of model '%s' failed" % cls)

            for out in (pt_out, int_out):
                self.assert_(out.dtype.type is SherpaFloat)
                self.assertEqual(out.shape, x.shape)

        self.assertEqual(count, 18)
