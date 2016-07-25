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

import operator
from numpy import arange
from sherpa.utils import SherpaFloat, SherpaTestCase
from sherpa.models.parameter import *
from sherpa.utils.err import ParameterErr


class test_parameter(SherpaTestCase):

    def setUp(self):
        self.p = Parameter('model', 'name', 0, -10, 10, -100, 100, 'units')
        self.afp = Parameter('model', 'name', 0, alwaysfrozen=True)

    def test_name(self):
        self.assertEqual(self.p.modelname, 'model')
        self.assertEqual(self.p.name, 'name')
        self.assertEqual(self.p.fullname, 'model.name')

    def test_alwaysfrozen(self):
        self.assertTrue(self.afp.frozen)
        self.afp.frozen = True
        self.assertTrue(self.afp.frozen)
        self.afp.freeze()
        self.assertTrue(self.afp.frozen)
        self.assertRaises(ParameterErr, self.afp.thaw)
        self.assertRaises(ParameterErr, setattr, self.afp, 'frozen', 0)

    def test_readonly_attributes(self):
        self.assertEqual(self.p.alwaysfrozen, False)
        self.assertRaises(AttributeError, setattr, self.p, 'alwaysfrozen', 1)
        self.assertEqual(self.p.hard_min, -100.0)
        self.assertRaises(AttributeError, setattr, self.p, 'hard_min', -1000)
        self.assertEqual(self.p.hard_max, 100.0)
        self.assertRaises(AttributeError, setattr, self.p, 'hard_max', 1000)

    def test_val(self):
        self.p.val = -7
        self.assertEqual(self.p.val, -7)
        self.assertTrue(type(self.p.val) is SherpaFloat)
        self.assertRaises(ValueError, setattr, self.p, 'val', 'ham')
        self.assertRaises(ParameterErr, setattr, self.p, 'val', -101)
        self.assertRaises(ParameterErr, setattr, self.p, 'val', 101)

    def test_min_max(self):
        for attr, sign in (('min', -1), ('max', 1)):
            setattr(self.p, attr, sign*99)
            val = getattr(self.p, attr)
            self.assertEqual(val, sign*99)
            self.assertTrue(type(val) is SherpaFloat)
            self.assertRaises(ValueError, setattr, self.p, attr, 'ham')
            self.assertRaises(ParameterErr, setattr, self.p, attr, -101)
            self.assertRaises(ParameterErr, setattr, self.p, attr, 101)

    def test_frozen(self):
        self.p.frozen = 1.0
        self.assertTrue(self.p.frozen is True)
        self.p.frozen = []
        self.assertTrue(self.p.frozen is False)
        self.assertRaises(TypeError, setattr, self.p.frozen, arange(10))
        self.p.link = self.afp
        self.assertTrue(self.p.frozen is True)
        self.p.link = None
        self.p.freeze()
        self.assertTrue(self.p.frozen is True)
        self.p.thaw()
        self.assertTrue(self.p.frozen is False)

    def test_link(self):
        self.p.link = None
        self.assertTrue(self.p.link is None)
        self.assertNotEqual(self.p.val, 17.3)
        self.afp.val = 17.3
        self.p.link = self.afp
        self.assertEqual(self.p.val, 17.3)
        self.p.unlink()
        self.assertTrue(self.p.link is None)
        self.assertRaises(ParameterErr, setattr, self.afp, 'link', self.p)
        self.assertRaises(ParameterErr, setattr, self.p, 'link', 3)
        self.assertRaises(ParameterErr, setattr, self.p, 'link',
                          3*self.p + 2)

    def test_iter(self):
        for part in self.p:
            self.assertTrue(part is self.p)


class test_composite_parameter(SherpaTestCase):

    def setUp(self):
        self.p = Parameter('model', 'p', 2)
        self.p2 = Parameter('model', 'p2', 4)

    def test_unop(self):
        for op in (abs, operator.neg):
            p = op(self.p)
            self.assertTrue(isinstance(p, UnaryOpParameter))
            self.assertEqual(p.val, op(self.p.val))

    def test_binop(self):
        ops = [operator.add, operator.sub, operator.mul,
               operator.floordiv, operator.truediv, operator.mod,
               operator.pow]

        if hasattr(operator, 'div'): # Python 2
            ops.append(operator.div)

        for op in ops:
            for p in (op(self.p, self.p2.val), op(self.p.val, self.p2),
                      op(self.p, self.p2)):
                self.assertTrue(isinstance(p, BinaryOpParameter))
                self.assertEqual(p.val, op(self.p.val, self.p2.val))

    def test_iter(self):
        p = 3 * self.p + self.p2
        parts = list(p)
        self.assertTrue(type(parts[0]) is BinaryOpParameter)
        self.assertTrue(type(parts[1]) is ConstantParameter)
        self.assertTrue(parts[2] is self.p)
        self.assertTrue(parts[3] is self.p2)

    def test_complex_expression(self):
        cmplx = (3 * self.p + self.p2) / (self.p ** 3.2)
        p = self.p.val
        p2 = self.p2.val
        self.assertEqual(cmplx.val, (3*p + p2) / (p ** 3.2))
