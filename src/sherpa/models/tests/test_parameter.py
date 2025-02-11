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
    BinaryOpParameter, ConstantParameter, hugeval, expand_par
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
    assert repr(comp) == "<BinaryOpParameter 'numpy.logaddexp2(model.p, model.p2)'>"


def test_binop_string_with_custom_ufunc():
    '''Repeat the previous test, but with a custom ufunc.'''
    def func(a, b):
        return a + b

    uf = np.frompyfunc(func, nin=2, nout=1)
    p, p2 = setUp_composite()
    comp = uf(p, p2)
    assert repr(comp) == "<BinaryOpParameter 'func(model.p, model.p2)'>"


def test_multiop_string_with_custom_ufunc():
    '''Get an operator which we do not support.'''
    def func(a, b, c):
        return a + b + c

    uf = np.frompyfunc(func, nin=3, nout=1)
    p, p2 = setUp_composite()
    p3 = Parameter("model", "c", 15)

    # This errors out because we do not suport three parameters.  It
    # does not seem worth checking the error message as this comes
    # from upstream code that could change.
    #
    with pytest.raises(TypeError):
        uf(p, p2, p3)


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
                              # (+a + b, "+(m.a) + m.b"),  currently not supported
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
                              ((a * b) % 2, "(m.a * m.b) % 2"),
                              # abs is special, how about other functions
                              (np.sin(a), "numpy.sin(m.a)"),
                              (np.sin(a + b), "numpy.sin(m.a + m.b)"),
                              (np.sin(((a) * (b))), "numpy.sin(m.a * m.b)"),
                              (np.cos((np.sin(a) * ((2) + np.tan(b)))),
                               "numpy.cos(numpy.sin(m.a) * (2 + numpy.tan(m.b)))"),
                              (np.cos(a)**2 + np.sin(a)**2,
                               "numpy.cos(m.a) ** 2 + numpy.sin(m.a) ** 2"),
                              (np.logaddexp(a, b),
                               "numpy.logaddexp(m.a, m.b)"),
                              (2 + np.logaddexp(a, b * 2),
                               "2 + numpy.logaddexp(m.a, m.b * 2)"),
                              (np.logaddexp(2 / a, b * 2) - (2 * c),
                               "numpy.logaddexp(2 / m.a, m.b * 2) - 2 * m.c"),
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
    assert explicit.name == "numpy.multiply(m1.a, numpy.add(m2.b, m4.xx))"

    # Check they evaluate to the same values.
    #
    yexp = 48
    assert implicit.val == pytest.approx(yexp)
    assert explicit.val == pytest.approx(yexp)


def test_basic_reset():
    """Try to understand the default_val setting."""

    p1 = Parameter("m1", "a", 4)
    assert p1.val == 4
    assert p1.default_val == 4
    assert p1.link is None

    p1.reset()
    assert p1.val == 4
    assert p1.default_val == 4
    assert p1.link is None


def test_set_default_to_val(check_str):
    """Try to understand the default_val setting."""

    p1 = Parameter("m1", "a", 4)
    p1.default_val = 2

    # Note that p1.val has not changed, unlike test_set_default_to_parameter
    assert p1.val == 4
    assert p1.default_val == 2
    assert p1.link is None

    check_str(str(p1),
              [ "val         = 4.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = 2.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])

    p1.reset()
    assert p1.val == 2
    assert p1.default_val == 2
    assert p1.link is None

    check_str(str(p1),
              [ "val         = 2.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = 2.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])


def test_set_default_to_parameter(check_str):
    """We support this, so make sure it is tested.

    There is the question of whether it makes sense to support this,
    but that's a different issue.
    """

    p1 = Parameter("m1", "a", 4)
    p2 = Parameter("m2", "b", 12)

    p1.default_val = 2 + p2

    # Note that p1.val has changed, unlike test_set_default_to_val
    assert p1.val == 14
    assert p1.default_val == 14
    assert p1.link is not None

    check_str(str(p1),
              [ "val         = 14.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = True",
                "link        = 2 + m2.b",
                "default_val = 14.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])

    p1.reset()

    assert p1.val == 14
    assert p1.default_val == 14
    assert p1.link is not None

    check_str(str(p1),
              [ "val         = 14.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = True",
                "link        = 2 + m2.b",
                "default_val = 14.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])


def test_set_default_to_parameter_to_parameter(check_str):
    """We support this, so make sure it is tested.

   There is the question of
    whether it makes sense to support this, but that's a different
    issue.

    """

    p1 = Parameter("m1", "a", 4)
    p2 = Parameter("m2", "b", 12)
    p3 = Parameter("m2", "b", 200)

    p2.default_val = p3 / 100
    p1.default_val = p2 - 5

    assert p1.link is not None
    assert p1.val == -3
    assert p1.default_val == -3


def test_set_default_to_value_after_parameter(check_str):
    """Check what happens here."""

    p1 = Parameter("m1", "a", 4)
    p2 = Parameter("m2", "b", 12)

    p1.default_val = 2 + p2

    p1.default_val = 3

    # Note: this is now back to the original value
    assert p1.val == 4
    assert p1.default_val == 3
    assert p1.link is None

    check_str(str(p1),
              [ "val         = 4.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = 3.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])

    p1.reset()

    assert p1.val == 3
    assert p1.default_val == 3

    check_str(str(p1),
              [ "val         = 3.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = 3.0",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])


def test_set_value_to_none(check_str):
    """This really should not be allowed. See #1992

    Treat as a regression test for now, but we could change the
    behaviour.
    """

    p1 = Parameter("m1", "a", 4)
    p1.val = None

    assert np.isnan(p1.val)
    assert np.isnan(p1.default_val)

    check_str(str(p1),
              [ "val         = nan",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = nan",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])

    p1.reset()

    assert np.isnan(p1.val)
    assert np.isnan(p1.default_val)


def test_set_default_to_none(check_str):
    """This really should not be allowed. See #1992

    Treat as a regression test for now, but we could change the
    behaviour.
    """

    p1 = Parameter("m1", "a", 4)
    p1.default_val = None

    assert p1.val == 4
    assert np.isnan(p1.default_val)

    check_str(str(p1),
              [ "val         = 4.0",
                "min         = -3.4028234663852886e+38",
                "max         = 3.4028234663852886e+38",
                "units       = ",
                "frozen      = False",
                "link        = None",
                "default_val = nan",
                "default_min = -3.4028234663852886e+38",
                "default_max = 3.4028234663852886e+38"
               ])

    p1.reset()

    assert np.isnan(p1.val)
    assert np.isnan(p1.default_val)


@pytest.mark.parametrize("dval,expected",
                         [(0, "minimum of 1"),
                          (20, "maximum of 5")])
def test_set_default_to_invalid(dval, expected):
    """Pick a value outside the min/max range."""

    p1 = Parameter("m1", "a", 4, min=1, max=5)

    with pytest.raises(ParameterErr,
                       match=f"parameter m1.a has a {expected}"):
        p1.default_val = dval


def test_link_check_very_simple():
    """Check link cycles.

    This is related to the old (pre GitHub) issue 12287 which was
    related to cycles in links but for which no test was added
    (at least for the check_simple / check_involved cases).

    """

    p1 = Parameter("m1", "a", 400)
    p2 = Parameter("m2", "b", 12)

    p2.val = p1

    # Should this report some sort of linked cycle or just clean
    # out the old setting?
    #
    p1.val = p2

    assert p2.link is None
    assert p1.link is not None

    assert p1.val == 12
    assert p2.val == 12


def test_link_check_simple():
    """Check link cycles. See issue #1993."""

    p1 = Parameter("m1", "a", 400)
    p2 = Parameter("m2", "b", 12)

    p2.val = p1 / 100

    # This should report some sort of linked cycle
    p1.val = p2 - 5

    with pytest.raises(RecursionError):
        p1.val


def test_link_check_involved():
    """Check link cycles. See issue #1993."""

    p1 = Parameter("m1", "a", 4)
    p2 = Parameter("m2", "b", 12)
    p3 = Parameter("m2", "b", 200)

    p3.val = 1 + p1
    p2.val = p3 / 100

    # This should report some sort of linked cycle
    p1.val = p2 - 5

    with pytest.raises(RecursionError):
        p1.val


def test_set_for_default_params():
    """Check some corner cases."""

    p = Parameter("mdl", "a", 5, min=0, max=10)
    assert p.val == 5
    assert p.default_val == 5
    assert p.min == 0
    assert p.max == 10
    assert p.default_min == 0
    assert p.default_max == 10

    p.set(default_min=-5)
    assert p.val == 5
    assert p.default_val == 5
    assert p.min == 0
    assert p.max == 10
    assert p.default_min == -5
    assert p.default_max == 10

    p.set(default_max=20)
    assert p.val == 5
    assert p.default_val == 5
    assert p.min == 0
    assert p.max == 10
    assert p.default_min == -5
    assert p.default_max == 20

    p.set(default_val=2)
    assert p.val == 5
    assert p.default_val == 2
    assert p.min == 0
    assert p.max == 10
    assert p.default_min == -5
    assert p.default_max == 20


def test_link_expression_add():
    """Check that the model expression behaves as expected.

    This encodes a number of checks which should be split out.
    """

    p1 = Parameter("p", "a", 10)
    p2 = Parameter("p", "b", -5)
    p2.frozen = True

    expr = 2 * p1 + p2 + 40

    assert expr.name == "2 * p.a + p.b + 40"
    assert expr.val == 55
    assert not expr.frozen

    # Deconstruct the tree
    assert isinstance(expr, BinaryOpParameter)
    assert expr.op == np.add
    assert isinstance(expr.lhs, BinaryOpParameter)
    assert isinstance(expr.rhs, ConstantParameter)

    # The constant term (40) should be frozen
    assert expr.rhs.frozen
    assert expr.rhs.val == pytest.approx(40)

    # Do not use instance on Parameter as the objects we expect here
    # are all sub-classes of it, so it adds no real power to the test.
    #
    assert expr.lhs.op == np.add
    assert isinstance(expr.lhs.lhs, BinaryOpParameter)
    assert type(expr.lhs.rhs) == Parameter

    assert expr.lhs.rhs.name == "b"
    assert expr.lhs.rhs.val == pytest.approx(-5)
    assert expr.lhs.rhs.frozen

    assert expr.lhs.lhs.op == np.multiply
    assert isinstance(expr.lhs.lhs.lhs, ConstantParameter)
    assert type(expr.lhs.lhs.rhs) == Parameter

    # The constant term (2) should be frozen
    assert expr.lhs.lhs.lhs.frozen
    assert expr.lhs.lhs.lhs.val == pytest.approx(2)

    assert expr.lhs.lhs.rhs.name == "a"
    assert expr.lhs.lhs.rhs.val == pytest.approx(10)
    assert not expr.lhs.lhs.rhs.frozen


def test_link_expression_sub():
    """Check that the model expression behaves as expected.

    This encodes a number of checks which should be split out.
    """

    p1 = Parameter("p", "a", 10)
    p2 = Parameter("p", "b", -5)
    p2.frozen = True

    expr = 2 * p1 - p2 + 40

    assert expr.name == "2 * p.a - p.b + 40"
    assert expr.val == 65
    assert not expr.frozen

    # Deconstruct the tree
    assert isinstance(expr, BinaryOpParameter)
    assert expr.op == np.add
    assert isinstance(expr.lhs, BinaryOpParameter)
    assert isinstance(expr.rhs, ConstantParameter)

    # The constant term (40) should be frozen
    assert expr.rhs.frozen
    assert expr.rhs.val == pytest.approx(40)

    # Do not use instance on Parameter as the objects we expect here
    # are all sub-classes of it, so it adds no real power to the test.
    #
    assert expr.lhs.op == np.subtract
    assert isinstance(expr.lhs.lhs, BinaryOpParameter)
    assert type(expr.lhs.rhs) == Parameter

    assert expr.lhs.rhs.name == "b"
    assert expr.lhs.rhs.val == pytest.approx(-5)
    assert expr.lhs.rhs.frozen

    assert expr.lhs.lhs.op == np.multiply
    assert isinstance(expr.lhs.lhs.lhs, ConstantParameter)
    assert type(expr.lhs.lhs.rhs) == Parameter

    # The constant term (2) should be frozen
    assert expr.lhs.lhs.lhs.frozen
    assert expr.lhs.lhs.lhs.val == pytest.approx(2)

    assert expr.lhs.lhs.rhs.name == "a"
    assert expr.lhs.lhs.rhs.val == pytest.approx(10)
    assert not expr.lhs.lhs.rhs.frozen


a = Parameter("m", "a", 5)
b = Parameter("n", "b", 10)
fb = Parameter("fn", "b", 10)
fb.freeze()


@pytest.mark.parametrize("expr,expected",
                         [(a, [a]),
                          (2 + a, [a]),
                          (a / 2, [a]),
                          (a * b, [a, b]),
                          (a + b + 2 * a, [a, b]),
                          (fb, [fb]),
                          (2 + fb, [fb]),
                          (fb + a, [fb, a]),
                          (fb - 2 * (a + fb + b), [fb, a, b]),
                          (2 ** a, [a]),
                          (a ** 2, [a]),
                          (a ** b, [a, b]),
                          (fb ** (1 + b), [fb, b]),
                          (a ** fb, [a, fb]),
                          (fb ** fb, [fb])
                          ])
def test_expand_par(expr, expected):
    """check the expansion."""

    res = expand_par(expr)
    for gterm, eterm in zip(res, expected):
        assert gterm == eterm

    assert len(res) == len(expected)
