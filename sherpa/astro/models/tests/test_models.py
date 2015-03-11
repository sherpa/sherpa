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
