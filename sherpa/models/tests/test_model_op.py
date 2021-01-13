#
#  Copyright (C) 2020, 2021  Smithsonian Astrophysical Observatory
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

from sherpa.astro.ui.utils import Session
from sherpa.models import basic
from sherpa.models.model import ArithmeticConstantModel, BinaryOpModel, UnaryOpModel
from sherpa.utils.err import ModelErr
from sherpa.utils.testing import requires_data, requires_fits, requires_xspec


def test_basic_unop_neg_raw():

    cpt = basic.Polynom2D()
    mdl = UnaryOpModel(cpt, np.negative, '<->')

    assert mdl.name == '<->(polynom2d)'
    assert mdl.op == np.negative
    assert mdl.opstr == '<->'
    assert mdl.ndim == 2


def test_basic_unop_neg():

    cpt = basic.Polynom2D()
    mdl = -cpt

    assert mdl.name == '-(polynom2d)'
    assert mdl.op == np.negative
    assert mdl.opstr == '-'
    assert mdl.ndim == 2


def test_basic_unop_abs_raw():

    cpt = basic.Polynom2D()
    mdl = UnaryOpModel(cpt, np.absolute, 'foo')

    assert mdl.name == 'foo(polynom2d)'
    assert mdl.op == np.absolute
    assert mdl.opstr == 'foo'
    assert mdl.ndim == 2


def test_basic_unop_abs():

    cpt = basic.Polynom2D()
    mdl = abs(cpt)

    assert mdl.name == 'abs(polynom2d)'
    assert mdl.op == np.absolute
    assert mdl.opstr == 'abs'
    assert mdl.ndim == 2


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
    assert mdl.opstr == 'xOx'
    assert len(mdl.parts) == 2
    assert mdl.parts[0] == l
    assert mdl.parts[1] == r
    assert mdl.ndim == 2


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
    assert mdl.name == f'(polynom2d {opstr} gauss2d)'
    assert mdl.op == op
    assert mdl.opstr == opstr
    assert len(mdl.parts) == 2
    assert mdl.parts[0] == l
    assert mdl.parts[1] == r
    assert mdl.ndim == 2


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
    assert mdl.ndim == 1

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
    assert mdl.ndim == 1

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
    assert mdl.ndim == 2

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

    with pytest.raises(ModelErr) as exc:
        m1 + m2

    assert str(exc.value) == 'Models do not match: 1D (box1d) and 2D (box2d)'


@pytest.mark.parametrize("val", [2, [1, 3, 4]])
def test_arithmeticconstantmodel(val):
    """Does ArithmeticConstantModel handle ndim?

    With the current design we can not distinguish the model
    dimensionality from the value to be wrapped, so ndim is
    set to None for all arguments.
    """

    mdl = ArithmeticConstantModel(val)
    assert mdl.ndim is None


def test_arithmeticconstantmodel_dimerr():

    x = np.ones(6).reshape(2, 3)
    with pytest.raises(ModelErr) as exc:
        ArithmeticConstantModel(x)

    assert str(exc.value) == 'The constant must be a scalar or 1D, not 2D'


@requires_xspec
def test_xspec_expression():
    """Check we can enforce XSPEC expressions"""

    from sherpa.astro import xspec

    # pick an additive and multiplicative model
    a1 = xspec.XSpowerlaw()
    m1 = xspec.XSwabs()

    m = 2 * (a1 + m1)

    assert a1.ndim == 1
    assert m1.ndim == 1
    assert m.ndim == 1


@requires_data
@requires_fits
@requires_xspec
def test_instrument_model(make_data_path):
    """Check the full response model"""

    from sherpa.astro import xspec

    s = Session()
    s.load_pha(make_data_path('3c273.pi'))

    a1 = xspec.XSpowerlaw()
    m1 = xspec.XSwabs()

    s.set_source(m1 * a1)

    src = s.get_source()
    mdl = s.get_model()
    assert src.ndim == 1
    assert mdl.ndim == 1


@requires_data
@requires_fits
def test_load_table_model(make_data_path):
    """What does load_table_model do?

    Even though this is not a FITS file, the code appears to
    need the I/O module to work.
    """

    s = Session()
    s.load_table_model('tbl', make_data_path('double.dat'))
    tbl = s.get_model_component('tbl')
    assert tbl.ndim is None
