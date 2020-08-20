#
#  Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

"""Do the unary and binary operators work for models?"""

from functools import reduce
import operator

import numpy as np

import pytest

from sherpa.models import basic
from sherpa.models.model import BinaryOpModel, UnaryOpModel


def test_basic_unop_neg_raw():

    cpt = basic.Polynom2D()
    mdl = UnaryOpModel(cpt, np.negative, '<->')

    assert mdl.name == '<->(polynom2d)'
    assert mdl.op == np.negative


def test_basic_unop_neg():

    cpt = basic.Polynom2D()
    mdl = -cpt

    assert mdl.name == '-(polynom2d)'
    assert mdl.op == np.negative


def test_basic_unop_abs_raw():

    cpt = basic.Polynom2D()
    mdl = UnaryOpModel(cpt, np.absolute, 'foo')

    assert mdl.name == 'foo(polynom2d)'
    assert mdl.op == np.absolute


def test_basic_unop_abs():

    cpt = basic.Polynom2D()
    mdl = abs(cpt)

    assert mdl.name == 'abs(polynom2d)'
    assert mdl.op == np.absolute


@pytest.mark.parametrize("op", [np.add, np.multiply, np.subtract, np.divide,
                                np.floor_divide, np.true_divide,
                                np.remainder, np.power])
def test_basic_binop_raw(op):

    l = basic.Polynom2D()
    r = basic.Gauss2D()
    mdl = BinaryOpModel(l, r, op, 'xOx')

    assert isinstance(mdl, BinaryOpModel)
    assert mdl.name == '(polynom2d xOx gauss2d)'
    assert mdl.op == op
    assert len(mdl.parts) == 2
    assert mdl.parts[0] == l
    assert mdl.parts[1] == r


@pytest.mark.parametrize("op,opstr",
                         [(np.add, '+'), (np.multiply, '*'),
                          (np.subtract, '-'), (np.divide, '/'),
                          (np.floor_divide, '//'), (np.true_divide, '/'),
                          (np.remainder, '%'), (np.power, '**')])
def test_basic_binop(op, opstr):

    l = basic.Polynom2D()
    r = basic.Gauss2D()
    mdl = op(l, r)

    assert isinstance(mdl, BinaryOpModel)
    assert mdl.name == '(polynom2d {} gauss2d)'.format(opstr)
    assert mdl.op == op
    assert len(mdl.parts) == 2
    assert mdl.parts[0] == l
    assert mdl.parts[1] == r


def test_eval_op():
    """Pick a "complicated" combination just to check."""

    x = np.asarray([2, 4, 5, 6, 7])

    m1 = basic.Const1D()
    m1.c0 = 10

    m2 = basic.Polynom1D()
    m2.c0 = 5
    m2.c1 = 1

    m3 = basic.Box1D()
    m3.xlow = 5
    m3.xhi = 6

    mdl = m1 + 2 * (m2 + (-m3))

    expected_m1 = 10 * np.ones(5)
    expected_m2 = 5 + np.asarray(x)
    expected_m3 = np.asarray([0, 0, 1, 1, 0])

    expected = expected_m1 + 2 * (expected_m2 - expected_m3)

    got = mdl(x)
    assert got == pytest.approx(expected)


def test_combine_models1d():
    """Check we can combine 1D models"""

    mdls = [basic.Box1D(), basic.Const1D(), basic.Gauss1D(), basic.NormGauss1D()]
    mdl = reduce(operator.add, mdls)
    assert isinstance(mdl, BinaryOpModel)

    # now multiply by a constant
    #
    mdl *= 2
    assert isinstance(mdl, BinaryOpModel)

    # Check we can call it as a 1D model; there is minimal checks
    # of the response.
    #
    bins = np.arange(2, 10, 2)
    y1 = mdl(bins)
    y2 = mdl(bins[:-1], bins[1:])

    assert y1.shape == (4, )
    assert y2.shape == (3, )


def test_combine_models2d():
    """Check we can combine 2D models"""

    mdls = [basic.Box2D(), basic.Const2D(), basic.Gauss2D(), basic.NormGauss2D()]
    mdl = reduce(operator.add, mdls)
    assert isinstance(mdl, BinaryOpModel)

    # now multiply by a constant
    #
    mdl *= 2
    assert isinstance(mdl, BinaryOpModel)

    # Check we can call it as a 2D model; there is minimal checks
    # of the response.
    #
    bins1 = np.arange(2, 10, 2)
    bins2 = bins1 - 4

    # only test non-integrated model
    y1 = mdl(bins1, bins2)

    assert y1.shape == (4, )


def test_combine_models_1d_2d():
    """What happens when we combine 1D and 2D models?"""

    m1 = basic.Box1D()
    m2 = basic.Box2D()
    mdl = m1 + m2

    assert isinstance(mdl, BinaryOpModel)

    # You can evaluate this, and has some sense, but it
    # unlikely to be what anyone wants.
    #
    y = mdl([1, 2, 3, 4, 5], [-5, -4, -3, -2, -1])
    assert y.shape == (5, )
