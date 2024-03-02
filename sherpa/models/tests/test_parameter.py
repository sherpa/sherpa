#
#  Copyright (C) 2007, 2016, 2018 - 2024
#  Smithsonian Astrophysical Observatory
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

import numpy as np

import pytest

from sherpa.models.basic import Gauss1D, Const1D, PowLaw1D
from sherpa.models.parameter import Parameter, UnaryOpParameter, \
    BinaryOpParameter, ConstantParameter, hugeval
from sherpa.utils.err import ParameterErr
from sherpa.utils.numeric_types import SherpaFloat
from sherpa import ui


def setUp_p():
    return Parameter('model', 'name', 0, -10, 10, -100, 100, 'units')


def setUp_afp():
    return Parameter('model', 'name', 0, alwaysfrozen=True)


def test_name():
    p = setUp_p()
    assert p.modelname == 'model'
    assert p.name == 'name'
    assert p.fullname == 'model.name'


def test_alwaysfrozen():
    afp = setUp_afp()
    assert afp.frozen
    afp.frozen = True
    assert afp.frozen
    afp.freeze()
    afp.frozen = True

    with pytest.raises(ParameterErr):
        afp.thaw()

    with pytest.raises(ParameterErr):
        afp.frozen = 0


def test_readonly_attributes():
    p = setUp_p()
    assert not p.alwaysfrozen
    with pytest.raises(AttributeError):
        p.alwaysfrozen = 1

    assert p.hard_min == -100.0
    with pytest.raises(AttributeError):
        p.hard_min = -1000

    assert p.hard_max == 100.0
    with pytest.raises(AttributeError):
        p.hard_max = 1000


def test_val():
    p = setUp_p()
    p.val = -7

    assert p.val == -7
    assert type(p.val) is SherpaFloat

    with pytest.raises(ValueError):
        p.val = 'ham'

    with pytest.raises(ParameterErr):
        p.val = -101

    with pytest.raises(ParameterErr):
        p.val = 101


def test_min_max():
    p = setUp_p()
    for attr, sign in (('min', -1), ('max', 1)):
        setattr(p, attr, sign * 99)
        val = getattr(p, attr)
        assert val == sign * 99
        assert type(val) is SherpaFloat

        with pytest.raises(ValueError):
            setattr(p, attr, 'ham')

        with pytest.raises(ParameterErr):
            setattr(p, attr, -101)

        with pytest.raises(ParameterErr):
            setattr(p, attr, 101)


def test_frozen():
    p = setUp_p()
    p.frozen = 1.0
    assert p.frozen is True

    p.frozen = []
    assert p.frozen is False

    with pytest.raises(ValueError):
        p.frozen = np.arange(10)

    afp = setUp_afp()
    p.link = afp
    assert p.frozen is True

    p.link = None
    p.freeze()
    assert p.frozen is True

    p.thaw()
    assert p.frozen is False


def test_link():
    p = setUp_p()
    assert p.link is None
    assert p.val == 0.0

    afp = setUp_afp()
    afp.val = 7.3
    p.link = afp
    assert p.val == 7.3

    p.unlink()
    assert p.link is None

    with pytest.raises(ParameterErr):
        afp.link = p

    with pytest.raises(ParameterErr):
        p.link = 3

    with pytest.raises(ParameterErr):
        p.link = 3 * p + 2


def test_iter():
    p = setUp_p()
    for part in p:
        assert part is p


def setUp_composite():
    p = Parameter('model', 'p', 2)
    p2 = Parameter('model', 'p2', 4)
    return p, p2


@pytest.mark.parametrize("op", [abs, operator.neg,
                                # Test a selection of numpy ufuncs
                                np.abs, np.sin, np.expm1
                                ])
def test_unop(op):
    porig = setUp_p()
    p = op(porig)
    assert isinstance(p, UnaryOpParameter)
    assert p.val == op(p.val)


def test_unop_string():
    '''Check that the string output is sensible.

    This implicitly also checks the logic to define a format string for
    printing the model representation.
    '''
    porig = setUp_p()
    assert repr(-porig) == "<UnaryOpParameter '-(model.name)'>"
    assert repr(np.cos(porig)) == "<UnaryOpParameter 'numpy.cos(model.name)'>"


def custom_func(a, b):
   return a + b

custom_ufunc = np.frompyfunc(custom_func, nin=2, nout=1)


@pytest.mark.parametrize("op", [operator.add, operator.sub, operator.mul,
                                operator.floordiv, operator.truediv, operator.mod,
                                operator.pow,
                                # Test a random selection of numpy ufuncs
                                np.add, np.divide, np.true_divide,
                                np.heaviside, np.greater, np.arctan2,
                                # and out own, custom made ufunc
                                custom_ufunc,
                                ])
def test_binop(op):
    p, p2 = setUp_composite()
    for comb in (op(p, p2.val), op(p.val, p2),
                 op(p, p2)):
        assert isinstance(comb, BinaryOpParameter)
        assert comb.val == op(p.val, p2.val)


def test_binop_string():
    '''Check that the string output is sensible.

    This implicitly also checks the logic to define the format string for
    printing the model representation.
    '''
    p, p2 = setUp_composite()
    comp = p + p2
    assert repr(comp) == "<BinaryOpParameter 'model.p + model.p2'>"
    comp = 1 + p * 2.7**p2
    assert repr(comp) == "<BinaryOpParameter '1 + model.p * 2.7 ** model.p2'>"
    comp = np.logaddexp2(p, p2)
    assert repr(comp) == "<BinaryOpParameter 'numpy.logaddexp2(model.p, (model.p2))'>"


def test_binop_string_with_custom_ufunc():
    '''Repeat the previous test, but with a custom ufunc.'''
    def func(a, b):
        return a + b

    uf = np.frompyfunc(func, nin=2, nout=1)
    p, p2 = setUp_composite()
    comp = uf(p, p2)
    assert repr(comp) == "<BinaryOpParameter 'func(model.p, (model.p2))'>"


class TestBrackets:
    """Provide a set of model instances for the tests."""

    a = Parameter('m', 'a', 2)
    b = Parameter('m', 'b', 3)
    c = Parameter('m', 'c', 4)
    d = Parameter('m', 'd', 5)

    # It would be nice to instead use a principled set of states,
    # but let's just try a somewhat-random set of expressions.
    #
    # This uses many of the same expressions from test_model_op.py
    # because why re-invent the wheel. Although there are differences.
    #
    @pytest.mark.parametrize("expr,expected",
                             [(a, "a"),  # TODO: should this be m.a?
                              (abs(a), "abs(m.a)"),
                              (abs(a) + b, "abs(m.a) + m.b"),
                              (b + abs(a), "m.b + abs(m.a)"),
                              (abs(a + b), "abs(m.a + m.b)"),
                              (abs(a + b * c), "abs(m.a + m.b * m.c)"),
                              (abs(a - b * c), "abs(m.a - m.b * m.c)"),
                              (abs((a + b) * c), "abs((m.a + m.b) * m.c)"),
                              (abs((a - b) * c), "abs((m.a - m.b) * m.c)"),
                              (abs((a - b) / c), "abs((m.a - m.b) / m.c)"),
                              (abs((a * b) - c), "abs(m.a * m.b - m.c)"),
                              (abs((a / b) - c), "abs(m.a / m.b - m.c)"),
                              (a * abs(b * (c + d)), "m.a * abs(m.b * (m.c + m.d))"),
                              (abs(b * (c + d)) * (a + d), "abs(m.b * (m.c + m.d)) * (m.a + m.d)"),
                              (-a, "-(m.a)"),
                              # (+a, "+(m.a)"),  currently not supported
                              (-a + b, "-(m.a) + m.b"),
                              (-a + 2, "-(m.a) + 2"),
                              # (+a + b, "+(m.a) + m.b"),  currenly not supported
                              (-(a + b), "-(m.a + m.b)"),
                              (-(a * b), "-(m.a * m.b)"),
                              (-(a - b), "-(m.a - m.b)"),
                              (-(a * b - c), "-(m.a * m.b - m.c)"),
                              (-(a - b * c), "-(m.a - m.b * m.c)"),
                              (a - a - b, "m.a - m.a - m.b"),
                              (a - (a - b), "m.a - (m.a - m.b)"),
                              (a - (b - (c - d)), 'm.a - (m.b - (m.c - m.d))'),
                              (a - (b + (c - d)), 'm.a - (m.b + m.c - m.d)'),
                              (b - (c + d), 'm.b - (m.c + m.d)'),
                              ((a - b) - (c + d), 'm.a - m.b - (m.c + m.d)'),
                              (a - (b - (c + d)), 'm.a - (m.b - (m.c + m.d))'),
                              (a - (b + (c + d)), 'm.a - (m.b + m.c + m.d)'),
                              (2 * (a + b) - c * 3, "2 * (m.a + m.b) - m.c * 3"),
                              (abs(2 * (a + b) - c * 3), "abs(2 * (m.a + m.b) - m.c * 3)"),
                              (a + a, "m.a + m.a"),
                              (a * b, "m.a * m.b"),
                              (a - a, "m.a - m.a"),
                              (a / b, "m.a / m.b"),
                              (a + b + c, "m.a + m.b + m.c"),
                              (a * b * c, "m.a * m.b * m.c"),
                              ((a * b) + c, "m.a * m.b + m.c"),
                              ((a + b) * c, "(m.a + m.b) * m.c"),
                              (a + (b * c), "m.a + m.b * m.c"),
                              (a * (b + c), "m.a * (m.b + m.c)"),
                              ((a + b) * (c + d), "(m.a + m.b) * (m.c + m.d)"),
                              ((a * b) * (c + d), "m.a * m.b * (m.c + m.d)"),
                              ((a + b) * (c * d), "(m.a + m.b) * m.c * m.d"),
                              ((a + (b * c) + d), "m.a + m.b * m.c + m.d"),
                              (100 * a * (b + c), "100 * m.a * (m.b + m.c)"),
                              (100 * (a * (b + c)), "100 * m.a * (m.b + m.c)"),
                              (a + b + 2 * c + d + a, "m.a + m.b + 2 * m.c + m.d + m.a"),
                              (a + b + c * 2 + d + a, "m.a + m.b + m.c * 2 + m.d + m.a"),
                              (a + b * (c - 2) + d + a, "m.a + m.b * (m.c - 2) + m.d + m.a"),
                              (a + b * (2 - c) + d + a, "m.a + m.b * (2 - m.c) + m.d + m.a"),
                              ((a + b + c) + (c + b + d + a), "m.a + m.b + m.c + m.c + m.b + m.d + m.a"),
                              ((a + b + c) + (c + b - d + a), "m.a + m.b + m.c + m.c + m.b - m.d + m.a"),
                              ((a + b + c) + (c + b - abs(d) + a), "m.a + m.b + m.c + m.c + m.b - abs(m.d) + m.a"),
                              ((a + b + c) * (c + b + d + a), "(m.a + m.b + m.c) * (m.c + m.b + m.d + m.a)"),
                              ((a + b + c) * (c + b - d + a), "(m.a + m.b + m.c) * (m.c + m.b - m.d + m.a)"),
                              ((a + b + c) * (c + b - abs(d) + a), "(m.a + m.b + m.c) * (m.c + m.b - abs(m.d) + m.a)"),
                              ((a * b * c) * (c + b + d + a), "m.a * m.b * m.c * (m.c + m.b + m.d + m.a)"),
                              ((a + b + c) * (c * b * d * a), "(m.a + m.b + m.c) * m.c * m.b * m.d * m.a"),
                              ((a + b + c) * (c * b + d * a), "(m.a + m.b + m.c) * (m.c * m.b + m.d * m.a)"),
                              (2 * a * 2, "2 * m.a * 2"),
                              (a * 2 * 2, "m.a * 2 * 2"),
                              (2 * a + 2 * (b + c - 4) * 3, "2 * m.a + 2 * (m.b + m.c - 4) * 3"),
                              # How about expressions with exponentation
                              (a**2, "m.a ** 2"),
                              (a**-2, "m.a ** -2"),
                              (a + b**2, "m.a + m.b ** 2"),
                              (a + (b + c)**2, "m.a + (m.b + m.c) ** 2"),
                              (a - b**(c - 2) - a, "m.a - m.b ** (m.c - 2) - m.a"),
                              ((a ** 2) ** b, "(m.a ** 2) ** m.b"),
                              (-a ** 2, "-(m.a ** 2)"),
                              ((-a)**2, "(-(m.a)) ** 2"),
                              # remainder and integer division, for fun
                              (a // 2, "m.a // 2"),
                              ((a + b) // 2, "(m.a + m.b) // 2"),
                              ((a * b) // 2, "(m.a * m.b) // 2"),
                              (a % 2, "m.a % 2"),
                              ((a + b) % 2, "(m.a + m.b) % 2"),
                              ((a * b) % 2, "(m.a * m.b) % 2")
                             ])
    def test_brackets(self, expr, expected):
        """Do we get the expected number of brackets?"""

        assert isinstance(expr, Parameter)
        assert expr.name == expected


def test_unary_op_precedence():
    """In checking precedences (adding TestBrackets case) DJB found
    this case so check the result."""

    p = Parameter("m", "x", 3)

    expr = p ** 2
    assert expr.name == "m.x ** 2"
    assert expr.val == pytest.approx(9)

    expr = (-p) ** 2
    assert expr.name == "(-(m.x)) ** 2"
    assert expr.val == pytest.approx(9)

    expr = -p ** 2
    assert expr.name == "-(m.x ** 2)"
    assert expr.val == pytest.approx(-9)

    expr = p ** -2
    assert expr.name == "m.x ** -2"
    assert expr.val == pytest.approx(1 / 9)


def test_iter_composite():
    p, p2 = setUp_composite()
    pnew = 3 * p + p2
    parts = list(pnew)
    assert len(parts) == 4
    assert type(parts[0]) is BinaryOpParameter
    assert type(parts[1]) is ConstantParameter
    assert parts[2] is p
    assert parts[3] is p2


def test_complex_expression():
    p, p2 = setUp_composite()
    cmplx = (3 * p + p2) / (p ** 3.2)
    p = p.val
    p2 = p2.val
    assert cmplx.val == (3 * p + p2) / (p ** 3.2)


def setUp_link():
    src1 = Gauss1D()
    src2 = Gauss1D()
    src1.pos.set(val=4, min=-10, max=10)
    src2.pos.set(val=5, min=-20, max=20)
    return src1, src2


def tst_pos(gauss, pos, minval, maxval, frozen=False,
            link=None):
    assert gauss.val == pos
    assert gauss.min == minval
    assert gauss.max == maxval
    assert gauss.frozen == frozen
    if gauss.link is None or link is None:
        assert gauss.link == link
    else:
        # These two asserts are just to point out that
        # the values we check against here are different
        # to the input values.
        #
        assert minval == -10
        assert maxval == 10
        tst_pos(gauss.link, pos, -20, 20)
    assert gauss.default_val == pos

    # Note: these are not minval/maxval but the defaults.
    assert gauss.default_min == -hugeval
    assert gauss.default_max == hugeval


def tst_unlink(src1, src2):
    # Are the min/max checks in tst_pos really adding anything
    # since we send in the parameter values.
    #
    src1.pos.unlink()
    tst_pos(src1.pos, src1.pos.default_val, -10, 10)

    src2.pos.unlink()
    tst_pos(src2.pos, 5, -20, 20)


def test_link_setup():
    src1, src2 = setUp_link()
    tst_pos(src1.pos, 4, -10, 10)
    tst_pos(src2.pos, 5, -20, 20)


def test_link_unlink_val():
    src1, src2 = setUp_link()
    tst_unlink(src1, src2)


def test_link_unlink_val_low_level():
    src1, src2 = setUp_link()

    src1.pos.val = src2.pos.val
    src1.pos.link = src2.pos

    tst_pos(src1.pos, 5, -10, 10, frozen=True, link=src1.pos)
    tst_pos(src2.pos, 5, -20, 20)

    tst_unlink(src1, src2)


def test_link_unlink_val_ui():
    src1, src2 = setUp_link()

    ui.link(src1.pos, src2.pos)

    tst_pos(src1.pos, 5, -10, 10, frozen=True, link=src1.pos)
    tst_pos(src2.pos, 5, -20, 20)

    tst_unlink(src1, src2)


def test_link_parameter_setting():
    """See https://github.com/sherpa/sherpa/issues/742

    See also test_link_parameter_evaluation.
    """

    mdl = PowLaw1D()
    lmdl = Const1D()

    # check we have the -10/10 range for the gamma
    # (so that if this changes we know to change this test)
    #
    assert mdl.gamma.min == pytest.approx(-10)
    assert mdl.gamma.max == pytest.approx(10)

    # Just check we can link the parameter
    lmdl.c0 = 2
    mdl.gamma = lmdl.c0
    assert mdl.gamma.val == pytest.approx(2)

    lmdl.c0 = 3
    assert mdl.gamma.val == pytest.approx(3)

    # What happens when we set lmdl.c0 to a value outside
    # the valid range for mdl?
    #
    # The error is raised when we try to access the value
    # of the parameter with the link, and *not* when we set
    # the parameter that is linked to.
    #
    lmdl.c0 = 23
    emsg = 'parameter powlaw1d.gamma has a maximum of 10'
    with pytest.raises(ParameterErr, match=emsg):
        mdl.gamma.val

    lmdl.c0 = -23
    emsg = 'parameter powlaw1d.gamma has a minimum of -10'
    with pytest.raises(ParameterErr, match=emsg):
        mdl.gamma.val

    lmdl.c0 = -2
    assert mdl.gamma.val == pytest.approx(-2)


def test_link_parameter_evaluation():
    """See also test_link_parameter_setting

    A version of this test, using an XSPEC table model, is found
    in sherpa/astro/xspec/tests/test_xspec_unit::test_xstbl_link_parameter_evaluation
    """

    # What happens when we try to evaluate a model whose
    # parameter is out-of-range thanks to a link?
    #
    mdl = PowLaw1D()
    lmdl = Const1D()

    grid = np.arange(1, 5)

    mdl.gamma = lmdl.c0
    lmdl.c0 = 2

    y2 = mdl(grid)
    assert (y2 > 0).all()

    lmdl.c0 = 12
    emsg = 'parameter powlaw1d.gamma has a maximum of 10'
    with pytest.raises(ParameterErr, match=emsg):
        mdl(grid)


@pytest.mark.parametrize("attr", ["val", "link"])
def test_link_manual(attr):
    """Check out parameter linking. Use .val or .link"""

    p = Parameter('model', 'eta', 2)
    q = Parameter('other', 'bob', 3)

    assert p.link is None
    assert q.link is None

    assert p.val == 2
    assert q.val == 3

    setattr(p, attr, 2 * q)
    assert p.val == 6

    assert q.link is None
    assert isinstance(p.link, BinaryOpParameter)
    assert isinstance(p.link.parts[0], ConstantParameter)
    assert p.link.parts[0].val == 2
    assert isinstance(p.link.parts[1], Parameter)
    assert p.link.parts[1] is q


def test_explicit_numpy_combination():
    """This was a question I wondered when developing test_brackets,
    so add a check.
    """

    p1 = Parameter("m1", "a", 4)
    p2 = Parameter("m2", "b", 2)
    p3 = Parameter("m4", "xx", 10)

    # These should be the same but check they are.
    #
    implicit = p1 * (p2 + p3)
    explicit = np.multiply(p1,
                           np.add(p2, p3))

    assert isinstance(implicit, BinaryOpParameter)
    assert isinstance(explicit, BinaryOpParameter)

    # The model version of this
    # test_model_op::test_explicit_numpy_combination()
    # has the names being the same, but they are not here.
    # Is this related to #1653
    #
    # assert explicit.name == implicit.name
    assert implicit.name == "m1.a * (m2.b + m4.xx)"
    assert explicit.name == "numpy.multiply(m1.a, (numpy.add(m2.b, m4.xx)))"

    # Check they evaluate to the same values.
    #
    yexp = 48
    assert implicit.val == pytest.approx(yexp)
    assert explicit.val == pytest.approx(yexp)
