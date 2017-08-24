#
#  Copyright (C) 2007, 2016  Smithsonian Astrophysical Observatory
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

import logging
import operator
import numpy
from sherpa.utils import SherpaFloat, SherpaTestCase
from sherpa.utils.err import ModelErr
from sherpa.models.model import *
from sherpa.models.model import ArithmeticModel
from sherpa.models.parameter import Parameter, tinyval
from sherpa.models.basic import Sin, Const1D


def my_sin(pars, x):
    return (pars[2].val *
            numpy.sin(2.0*numpy.pi * (x - pars[1].val) / pars[0].val))


class test_model(SherpaTestCase):

    def setUp(self):
        self.m = Sin('m')

    def test_name(self):
        self.assertEqual(self.m.name, 'm')

    def test_iter(self):
        for part in self.m:
            self.assertIs(part, self.m)

    def test_num_pars(self):
        self.assertEqual(len(self.m.pars), 3)

    def test_par_names(self):
        self.assertEqual([p.name for p in self.m.pars],
                         ['period', 'offset', 'ampl'])

    def test_getpar(self):
        for par in (self.m.period, self.m.PerioD, self.m.PERIod):
            self.assertIs(par, self.m.pars[0])
        self.assertRaises(AttributeError, getattr, self.m, 'perio')

    def test_setpar(self):
        self.assertNotEqual(self.m.offset.val, 17.0)
        self.m.offset = 17
        self.assertEqual(self.m.offset.val, 17.0)
        self.m.ofFseT = 18
        self.assertEqual(self.m.offset.val, 18.0)
        self.assertRaises(AttributeError, setattr, self.m, 'ofset', 19)

    def test_calc_and_call(self):
        x = numpy.arange(10.0)
        refvals = my_sin(self.m.pars, x)
        pars = [p.val for p in self.m.pars]
        for vals in (self.m.calc(pars, x), self.m(x)):
            self.assertEqualWithinTol(vals, refvals, 1e-12)

    def test_get_thawed_pars(self):
        tp = [p.val for p in self.m.pars if not p.frozen]
        self.assertEqual(self.m.thawedpars, tp)

    def test_set_thawed_pars(self):
        pars = [7, 8, 9]
        self.assertNotEqual(pars, self.m.thawedpars)
        self.m.thawedpars = pars
        self.assertEqual(pars, self.m.thawedpars)
        self.assertRaises(ModelErr, setattr, self.m, 'thawedpars', pars[:2])
        self.assertRaises(ValueError, setattr, self.m, 'thawedpars',
                          [1, 2, 'ham'])

        pars[0] = self.m.pars[0].hard_min / 10
        pars[1] = self.m.pars[1].hard_max * 10

        logger = logging.getLogger('sherpa')
        old_level = logger.getEffectiveLevel()
        logger.setLevel(logging.ERROR)

        try:
            self.m.thawedpars = pars
        finally:
            logger.setLevel(old_level)

        self.assertEqual(self.m.pars[0].val, self.m.pars[0].min)
        self.assertEqual(self.m.pars[1].val, self.m.pars[1].max)

    def test_get_mins_maxes(self):
        self.assertEqual(self.m.thawedparmins,
                         [p.min for p in self.m.pars if not p.frozen])
        self.assertEqual(self.m.thawedparmaxes,
                         [p.max for p in self.m.pars if not p.frozen])
        self.assertEqual(self.m.thawedparhardmins,
                         [p.hard_min for p in self.m.pars if not p.frozen])
        self.assertEqual(self.m.thawedparhardmaxes,
                         [p.hard_max for p in self.m.pars if not p.frozen])


class test_composite_model(SherpaTestCase):

    def setUp(self):
        self.m  = Const1D('m')
        self.m2 = Const1D('m2')
        self.m.c0  = 2
        self.m2.c0 = 4
        self.s = Sin('s')
        self.x = 1.0
        self.xx = numpy.arange(-10.0, 10.0)

    def test_iter(self):
        m = 3 * self.m + self.m2
        parts = list(m)
        self.assertIs(type(parts[0]), BinaryOpModel)
        self.assertIs(type(parts[1]), ArithmeticConstantModel)
        self.assertIs(parts[2], self.m)
        self.assertIs(parts[3], self.m2)

    def test_unop(self):
        for op in (abs, operator.neg):
            m = op(self.m)
            self.assertTrue(isinstance(m, UnaryOpModel))
            self.assertEqual(m(self.x), op(self.m(self.x)))

    def test_binop(self):
        ops = [operator.add, operator.sub, operator.mul,
               operator.floordiv, operator.truediv, operator.mod,
               operator.pow]

        if hasattr(operator, 'div'): # Python 2
            ops.append(operator.div)

        for op in ops:
            for m in (op(self.m, self.m2.c0.val), op(self.m.c0.val, self.m2),
                      op(self.m, self.m2)):
                self.assertTrue(isinstance(m, BinaryOpModel))
                self.assertEqual(m(self.x), op(self.m.c0.val, self.m2.c0.val))

    def test_complex_expression(self):
        cmplx = (3 * self.m + self.m2) / (self.m ** 3.2)
        m = self.m(self.x)
        m2 = self.m2(self.x)
        self.assertEqual(cmplx(self.x), (3*m + m2) / (m ** 3.2))

    def test_filter(self):
        m = self.s[::2]
        self.assertIs(type(m), FilterModel)
        self.assertTrue(numpy.all(m(self.xx) == self.s(self.xx)[::2]))

    def test_nested(self):
        for func in (numpy.sin, self.s):
            m = self.m.apply(func)
            self.assertIs(type(m), NestedModel)
            self.assertEqual(m(self.x), func(self.m(self.x)))


# Test support for renamed parameters by sub-classing
# the Sin model (which lets the tests be re-used).
#
class RenamedPars(Sin):

    def __init__(self, name='renamedpars'):
        self._renamedpars = [('norm', 'ampl')]
        Sin.__init__(self, name)


class test_model_renamed(test_model):

    def setUp(self):
        self.m = RenamedPars('m')

    def test_getpar_rename(self):
        for par in (self.m.norm, self.m.NorM, self.m.NOrm):
            self.assertIs(par, self.m.pars[2])

    def test_setpar_rename(self):
        self.m.ampl = 1
        self.assertNotEqual(self.m.ampl.val, 12.0)
        self.m.norm = 12
        self.assertEqual(self.m.ampl.val, 12.0)
        self.m.NoRM = 18
        self.assertEqual(self.m.ampl.val, 18.0)
        self.m.ampl = 1
        self.assertEqual(self.m.ampl.val, 1.0)


# During testing an XSPEC model was found to refer to one of its
# parameters with a different case to how it was created, which
# caused problems. Replicate this problem case here.
#
class ParameterCase(ArithmeticModel):
    """Re-implemenent Sin model so can copy tests"""

    def __init__(self, name='parametercase'):

        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0)
        self._renamedpars = [('NORM', 'ampl')]
        pars = (self.perioD, self.oFFSEt, self.NORM)

        self._basemodel = Sin()

        ArithmeticModel.__init__(self, name, pars)

    def calc(self, *args, **kwargs):
        for par in self.pars:
            setattr(self._basemodel, par.name, par.val)

        self._basemodel.integrate = self.integrate
        return self._basemodel.calc(*args, **kwargs)


class test_model_parametercase_instance(test_model_renamed):

    def setUp(self):
        self.m = ParameterCase()

    def test_name(self):
        self.assertEqual(self.m.name, 'parametercase')
        self.assertEqual(self.m.NORM.name, 'ampl')
