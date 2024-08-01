#
#  Copyright (C) 2020 - 2024
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

"""Do the unary and binary operators work for models?

Other tests related to model expressions are also made here.

"""

from functools import reduce
import operator
import re

import numpy as np

import pytest

from sherpa.astro.instrument import ARF1D, RMF1D, RSPModelNoPHA, create_arf, create_delta_rmf
from sherpa.astro.ui.utils import Session
from sherpa.instrument import PSFModel
from sherpa.models import basic
from sherpa.models.model import ArithmeticConstantModel, \
    ArithmeticFunctionModel, BinaryOpModel, UnaryOpModel, Model, \
    model_deconstruct
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


def test_basic_unop_pos():

    cpt = basic.Polynom2D()
    mdl = +cpt

    assert mdl.name == '+(polynom2d)'
    assert mdl.op == np.positive
    assert mdl.opstr == '+'
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
    assert mdl.name == 'polynom2d xOx gauss2d'
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
    assert mdl.name == f'polynom2d {opstr} gauss2d'
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


def test_eval_add_sub_op():
    """Another version of test_eval_op focussed on + and - unary ops"""

    x = np.asarray([2, 4, 5, 6, 7])

    m1 = basic.Const1D("c")
    m1.c0 = 10

    m2 = basic.Polynom1D("p")
    m2.c0 = 5
    m2.c1 = 1

    m3 = basic.Box1D("b")
    m3.xlow = 5
    m3.xhi = 6

    mdl = m1 - (+m2 - m3)
    assert mdl.ndim == 1

    assert mdl.name == "c - (+(p) - b)"

    bins = np.arange(2, 10, 2)
    yexp = m1(bins) + m3(bins) - m2(bins)
    y = mdl(bins)
    assert y == pytest.approx(yexp)


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

    with pytest.raises(ModelErr,
                       match=re.escape("Models do not match: 1D (box1d) and 2D (box2d)")):
        m1 + m2


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
    with pytest.raises(ModelErr,
                       match="The constant must be a scalar or 1D, not 2D"):
        ArithmeticConstantModel(x)


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


@pytest.mark.parametrize("flag", [True, False])
def test_unop_integrate(flag):
    """Is the integrate flag carried over"""

    mdl = basic.Const1D()
    mdl.integrate = flag

    umdl = -mdl
    with pytest.raises(AttributeError):
        assert umdl.integrate == flag


def test_aconstant_integrate():
    """Do we have an integrate setting?

    test_unop_integrate_unset and others use an
    ArithmeticConstantModel (implicitly) because at the time of
    writing it doesn't have an integrate setting. Make this a
    regression test (it is not a priori a problem if it does change
    behavior, we just want to know we have a test).

    """

    mdl = ArithmeticConstantModel(3.2)
    with pytest.raises(AttributeError):
        mdl.integrate


def test_unop_integrate_unset():
    """Is the integrate flag not set for an unary op model"""

    # UnaryOpModel converts the constant to an
    # ArithmeticConstantModel.
    #
    mdl = UnaryOpModel(4, np.negative, "-")

    with pytest.raises(AttributeError):
        mdl.integrate


@pytest.mark.parametrize("flag", [True, False])
def test_binop_integrate_same(flag):
    """Is the integrate flag carried over when the same"""

    mdl1 = basic.Const1D()
    mdl1.integrate = flag

    mdl2 = basic.Const1D()
    mdl2.integrate = flag

    mdl = mdl1 + mdl2
    with pytest.raises(AttributeError):
        assert mdl.integrate == flag


@pytest.mark.parametrize("flag", [True, False])
def test_binop_integrate_different(flag):
    """Is the integrate flag not carried over when different"""

    mdl1 = basic.Const1D()
    mdl1.integrate = flag

    mdl2 = basic.Const1D()
    mdl2.integrate = not flag

    mdl = mdl1 + mdl2
    with pytest.raises(AttributeError):
        mdl.integrate


def test_binop_integrate_unset():
    """Is the integrate flag not set for a binary op model"""

    # The ArithmeticConstantModel, which the binary-op model casts
    # the constant terms, has no integrate setting.
    #
    mdl = BinaryOpModel(4, 4, np.add, '+')
    with pytest.raises(AttributeError):
        mdl.integrate


@pytest.mark.parametrize("m1,m2", [(1, basic.Const1D()),
                                   (basic.Const1D(), 1),
                                   (1, basic.Scale1D()),
                                   (basic.Scale1D(), 1)])
def test_binop_arithmeticmodel_integrate(m1, m2):
    """Do we have an integrate setting?

    This is a regression test.
    """

    mdl = m1 + m2
    assert not hasattr(mdl, "integrate")


AMODEL = ArithmeticFunctionModel(np.sin)


@pytest.mark.parametrize("m1,m2", [(AMODEL, basic.Const1D()),
                                   (basic.Const1D(), AMODEL),
                                   (AMODEL, basic.Scale1D()),
                                   (basic.Scale1D(), AMODEL)])
def test_binop_arithmeticfunction_integrate(m1, m2):
    """What's the integrate setting when ArithmeticFunctionModel is present?

    This is a regression test.
    """

    mdl = m1 + m2
    assert not hasattr(mdl, "integrate")


@pytest.mark.parametrize("m1,m2", [(AMODEL, 2),
                                   (2, AMODEL)])
def test_binop_arithmeticfunction_works(m1, m2):
    """Check we can evaluate a function model in a binop"""

    grid = np.asarray([-0.2, -0.1, 0.1, 0.3])

    # There's no support for creating these models easily.
    #
    mdl = BinaryOpModel(m1, m2, np.multiply, '*')
    y = mdl(grid)

    expected = 2 * np.sin(grid)
    assert y == pytest.approx(expected)


class FakeResponse1D:
    """sherpa.astro.instrument.Response1D requires a PHA. This doesn't.

    This has limited functionality.
    """

    def __init__(self, arf, rmf):
        self.arf = arf
        self.rmf = rmf

    def __call__(self, model):
        return RSPModelNoPHA(self.arf, self.rmf,
                             self.arf.exposure * model)


class TestBrackets:
    """Provide a set of model instances for the tests."""

    a = basic.Const1D('a')
    b = basic.Const1D('b')
    c = basic.Const1D('c')
    d = basic.Const1D('d')

    # We don't need to 'load' the model data to use it here
    tm = basic.TableModel('tm')

    # Convolution-style model (PSF)
    cm = PSFModel('cm', basic.Const1D('cmflat'))

    # Convolution-style model (PHA)
    #
    egrid = np.arange(0.1, 0.5, 0.1)
    chans = np.arange(1, egrid.size)
    fake_arf = create_arf(egrid[:-1], egrid[1:], exposure=100.0)
    fake_rmf = create_delta_rmf(egrid[:-1], egrid[1:])

    arf = ARF1D(fake_arf)
    rmf = RMF1D(fake_rmf)
    rsp = FakeResponse1D(fake_arf, fake_rmf)

    # It would be nice to instead use a principled set of states,
    # but let's just try a somewhat-random set of expressions.
    #
    @pytest.mark.parametrize("model,expected",
                             [(a, "a"),
                              (abs(a), "abs(a)"),
                              (abs(a) + b, "abs(a) + b"),
                              (b + abs(a), "b + abs(a)"),
                              (abs(a + b), "abs(a + b)"),
                              (abs(a + b * c), "abs(a + b * c)"),
                              (abs(a - b * c), "abs(a - b * c)"),
                              (abs((a + b) * c), "abs((a + b) * c)"),
                              (abs((a - b) * c), "abs((a - b) * c)"),
                              (abs((a - b) / c), "abs((a - b) / c)"),
                              (abs((a * b) - c), "abs(a * b - c)"),
                              (abs((a / b) - c), "abs(a / b - c)"),
                              (a * abs(b * (c + d)), "a * abs(b * (c + d))"),
                              (abs(b * (c + d)) * (a + d), "abs(b * (c + d)) * (a + d)"),
                              (-a, "-(a)"),
                              (+a, "+(a)"),
                              ((-a) + ((b)), "-(a) + b"),
                              # the following is ugly but is valid Python
                              ((-(a)) + ((+b)), "-(a) + +(b)"),
                              (-a + b, "-(a) + b"),
                              (-a + 2, "-(a) + 2.0"),
                              (+a + b, "+(a) + b"),
                              (-(a + b), "-(a + b)"),
                              (-(a * b), "-(a * b)"),
                              (-(a - b), "-(a - b)"),
                              (-(a * b - c), "-(a * b - c)"),
                              (-(a - b * c), "-(a - b * c)"),
                              (a - a - b, "a - a - b"),
                              (a - (a - b), "a - (a - b)"),
                              (a - (b - (c - d)), 'a - (b - (c - d))'),
                              (a - (b + (c - d)), 'a - (b + c - d)'),
                              (b - (c + d), 'b - (c + d)'),
                              ((a - b) - (c + d), 'a - b - (c + d)'),
                              (a - (b - (c + d)), 'a - (b - (c + d))'),
                              (a - (b + (c + d)), 'a - (b + c + d)'),
                              (2 * (a + b) - c * 3, "2.0 * (a + b) - c * 3.0"),
                              (abs(2 * (a + b) - c * 3), "abs(2.0 * (a + b) - c * 3.0)"),
                              (a + a, "a + a"),
                              (a * b, "a * b"),
                              (a - a, "a - a"),
                              (a / b, "a / b"),
                              (a + b + c, "a + b + c"),
                              (a * b * c, "a * b * c"),
                              ((a * b) + c, "a * b + c"),
                              ((a + b) * c, "(a + b) * c"),
                              (a + (b * c), "a + b * c"),
                              (a * (b + c), "a * (b + c)"),
                              ((a + b) * (c + d), "(a + b) * (c + d)"),
                              ((a * b) * (c + d), "a * b * (c + d)"),
                              ((a + b) * (c * d), "(a + b) * c * d"),
                              ((a + (b * c) + d), "a + b * c + d"),
                              (100 * a * (b + c), "100.0 * a * (b + c)"),
                              (100 * (a * (b + c)), "100.0 * a * (b + c)"),
                              (a + b + 2 * c + d + a, "a + b + 2.0 * c + d + a"),
                              (a + b + c * 2 + d + a, "a + b + c * 2.0 + d + a"),
                              (a + b * (c - 2) + d + a, "a + b * (c - 2.0) + d + a"),
                              (a + b * (2 - c) + d + a, "a + b * (2.0 - c) + d + a"),
                              ((a + b + c) + (c + b + d + a), "a + b + c + c + b + d + a"),
                              ((a + b + c) + (c + b - d + a), "a + b + c + c + b - d + a"),
                              ((a + b + c) + (c + b - abs(d) + a), "a + b + c + c + b - abs(d) + a"),
                              ((a + b + c) * (c + b + d + a), "(a + b + c) * (c + b + d + a)"),
                              ((a + b + c) * (c + b - d + a), "(a + b + c) * (c + b - d + a)"),
                              ((a + b + c) * (c + b - abs(d) + a), "(a + b + c) * (c + b - abs(d) + a)"),
                              ((a * b * c) * (c + b + d + a), "a * b * c * (c + b + d + a)"),
                              ((a + b + c) * (c * b * d * a), "(a + b + c) * c * b * d * a"),
                              ((a + b + c) * (c * b + d * a), "(a + b + c) * (c * b + d * a)"),
                              (2 * a * 2, "2.0 * a * 2.0"),
                              (a * 2 * 2, "a * 2.0 * 2.0"),
                              (2 * a + 2 * (b + c - 4) * 3, "2.0 * a + 2.0 * (b + c - 4.0) * 3.0"),
                              (tm * (a + b) + tm * (a * b),
                               'tm * (a + b) + tm * a * b'),
                              (tm * (a + b) + tm * (a * (b + 3)),
                               'tm * (a + b) + tm * a * (b + 3.0)'),
                              (cm(a) + b, 'cm(a) + b'),
                              (a * cm(b + c), 'a * cm(b + c)'),
                              (a + cm(b + 2 * d + c),
                               'a + cm(b + 2.0 * d + c)'),
                              (arf(b * (c * d)) + d,
                               "apply_arf(100.0 * b * c * d) + d"),
                              (a + 2 * arf(b * (c + d)),
                               "a + 2.0 * apply_arf(100.0 * b * (c + d))"),
                              (a + 2 * rmf(b * (c + d)),
                               "a + 2.0 * apply_rmf(b * (c + d))"),
                              # Manually combining RMF1D and ARF1D is interesting as we would normally
                              # use RSPModelNoPHA
                              (a * rmf(arf(a * (b + c * d))) + d * arf(a + b),
                               'a * apply_rmf(apply_arf(100.0 * a * (b + c * d))) + d * apply_arf(100.0 * (a + b))'),
                              # Repeat but with rsp instead
                              (a * rsp(a * (b + c * d)) + d * arf(a + b),
                               'a * apply_rmf(apply_arf(100.0 * a * (b + c * d))) + d * apply_arf(100.0 * (a + b))'),
                              (arf(b * (c + d)),
                               "apply_arf(100.0 * b * (c + d))"),
                              # How about expressions with exponentation
                              (a**2, "a ** 2.0"),
                              (a**-2, "a ** -2.0"),
                              (a + b**2, "a + b ** 2.0"),
                              (a + (b + c)**2, "a + (b + c) ** 2.0"),
                              (a - b**(c - 2) - a, "a - b ** (c - 2.0) - a"),
                              ((a ** 2) ** b, "(a ** 2.0) ** b"),
                              (-a ** 2, "-(a ** 2.0)"),
                              ((-a)**2, "(-(a)) ** 2.0"),
                              # remainder and integer division, for fun
                              (a // 2, "a // 2.0"),
                              ((a + b) // 2, "(a + b) // 2.0"),
                              ((a * b) // 2, "(a * b) // 2.0"),
                              (a % 2, "a % 2.0"),
                              ((a + b) % 2, "(a + b) % 2.0"),
                              ((a * b) % 2, "(a * b) % 2.0")
                             ])
    def test_brackets(self, model, expected):
        """Do we get the expected number of brackets?"""

        assert model.name == expected

        # Can we check that the string expression, when evaluated
        # as a model, returns the same result?
        #
        # Now, thanks to how models are converted to strings, the
        # ARF/RMF cases do not work, so skip if any expected string
        # contains "apply".
        #
        if expected.find("apply") != -1:
            return

        got = eval(expected, None,
                   {"a": self.a,
                    "b": self.b,
                    "c": self.c,
                    "d": self.d,
                    "tm": self.tm,
                    "cm": self.cm})
        assert isinstance(got, Model)

        # Just because we can
        assert got.name == expected

        # Evaluate the two models and check they get the same. Since
        # most of the models are very simple (e.g. return 1 for each
        # bin) this does not test that much, but should be sufficient
        # here.
        #
        # However, those models that include convolution-style
        # expressions need to be folded through a data object which we
        # do not want to set up here, so we skip those tests. There is
        # a similar issue with the table model.
        #
        if expected.find("cm") != -1 or expected.find("tm") != -1:
            return

        x = [2, 5, 7]
        ymdl = model(x)
        ygot = got(x)

        assert ygot == pytest.approx(ymdl)


def test_explicit_numpy_combination():
    """This was a question I wondered when developing test_brackets,
    so add a check.
    """

    mdl1 = basic.Scale1D("mdl1")
    mdl2 = basic.Box1D("mdl2")
    mdl3 = basic.Gauss1D("mdl3")

    mdl1.c0 = 4
    mdl2.xlow = 5
    mdl2.xhi = 15
    mdl2.ampl = 2
    mdl3.pos = 10
    mdl3.fwhm = 5
    mdl3.ampl = 10

    # These should be the same but check they are.
    #
    implicit = mdl1 * (mdl2 + mdl3)
    explicit = np.multiply(mdl1,
                           np.add(mdl2, mdl3))

    assert isinstance(implicit, BinaryOpModel)
    assert isinstance(explicit, BinaryOpModel)

    # Check the names are the same.
    #
    assert explicit.name == implicit.name

    # Check they evaluate to the same values.
    #
    x = np.arange(4, 14, 2)
    y2 = mdl2(x)
    y3 = mdl3(x)
    yexp = 4 * (y2 + y3)

    assert implicit(x) == pytest.approx(yexp)
    assert explicit(x) == pytest.approx(yexp)

    # Have an actual test, just in case,
    assert yexp == pytest.approx([0.73812041, 14.78302164,
                                  33.66851795, 48, 33.66851795])


# Models for testing model_deconstruct. It is tempting to add these as
# attributes in a class, as is done above for TestBrackets, but we
# need to call methods and set attributes for these which would
# complicate things, so leave as module-level symbols.
#
BOX1 = basic.Box1D("b1")
BOX2 = basic.Box1D("b2")
GAUSS1 = basic.Gauss1D("g1")
GAUSS2 = basic.Gauss1D("g2")
SCALE1 = basic.Scale1D("s1")

BOX1.xlow = -10
BOX1.xhi = 10
BOX1.ampl.set(5, max=5)

BOX2.xlow = -10
BOX2.xhi = 10
BOX2.ampl.set(4, max=4)

GAUSS1.pos = -2
GAUSS1.fwhm = 5
GAUSS1.ampl = 8

GAUSS2.pos = 1
GAUSS2.fwhm = 5
GAUSS2.ampl = 6

# This isn't really needed, but do it to point out we expect these
# models not to change.
#
BOX1.freeze()
BOX2.freeze()
GAUSS1.freeze()
GAUSS2.freeze()
SCALE1.freeze()


@pytest.mark.parametrize("model,expecteds",
                         [(BOX1, ["b1"]),
                          # unary operator
                          (-BOX1, ["-(b1)"]),
                          (-(-BOX1), ["-(-(b1))"]),
                          # binary operator of singletons
                          (GAUSS1 + GAUSS2, ["g1", "g2"]),
                          (GAUSS1 - GAUSS2, ["g1", "-(g2)"]),
                          (GAUSS1 * GAUSS2, ["g1 * g2"]),
                          (GAUSS1 / GAUSS2, ["g1 / g2"]),
                          (GAUSS1 // GAUSS2, ["g1 // g2"]),
                          # sneaky test with a constant
                          (BOX1 + 2, ["b1", "2.0"]),
                          (2 + BOX1, ["2.0", "b1"]),
                          (BOX1 * 2, ["b1 * 2.0"]),
                          (2 * BOX1, ["2.0 * b1"]),
                          # more-complex binary operator, but still 1 term on one side
                          #
                          (BOX1 + (GAUSS1 + GAUSS2), ["b1", "g1", "g2"]),
                          ((BOX1 + GAUSS1) + GAUSS2, ["b1", "g1", "g2"]),
                          (BOX1 + (GAUSS1 * GAUSS2), ["b1", "g1 * g2"]),
                          ((BOX1 + GAUSS1) * GAUSS2, ["b1 * g2", "g1 * g2"]),
                          (BOX1 * (GAUSS1 + GAUSS2), ["b1 * g1", "b1 * g2"]),
                          ((BOX1 * GAUSS1) + GAUSS2, ["b1 * g1", "g2"]),
                          (BOX1 * (GAUSS1 * GAUSS2), ["b1 * g1 * g2"]),
                          ((BOX1 * GAUSS1) * GAUSS2, ["b1 * g1 * g2"]),
                          # add in negation and divison
                          (BOX1 - (GAUSS1 + GAUSS2), ["b1", "-(g1)", "-(g2)"]),
                          ((BOX1 - GAUSS1) + GAUSS2, ["b1", "-(g1)", "g2"]),
                          (BOX1 - (GAUSS1 * GAUSS2), ["b1", "-(g1 * g2)"]),
                          ((BOX1 * GAUSS1) - GAUSS2, ["b1 * g1", "-(g2)"]),
                          (BOX1 + (BOX2 / GAUSS1), ["b1", "b2 / g1"]),
                          ((BOX1 + BOX2) / GAUSS1, ["b1 / g1", "b2 / g1"]),
                          (GAUSS1 / (BOX1 * BOX2), ["g1 / (b1 * b2)"]),
                          ((GAUSS1 / BOX1) * BOX2, ["g1 / b1 * b2"]),
                          (GAUSS1 / (BOX1 + BOX2), ["g1 / (b1 + b2)"]),
                          # What happens with a unary term applied to a complex
                          # expression? There is no expansion, which is not
                          # ideal but safest. However, if written as a BinOp
                          # it does get expanded, which is a bit surprising!
                          (-(BOX1 + BOX2 * GAUSS1 + GAUSS2),
                           ["-(b1 + b2 * g1 + g2)"]),
                          (BOX2 - (BOX1 + BOX2 * GAUSS1 + GAUSS2),
                           ["b2", "-(b1)", "-(b2 * g1)", "-(g2)"]),
                          # note that we do not try to simplify the expressions
                          ((1 / GAUSS1) * (BOX1 - BOX2),
                           ["1.0 / g1 * b1", "1.0 / g1 * -(b2)"]),
                          # try combining binary operators
                          # - addition
                          ((BOX1 + BOX2) + (GAUSS1 + GAUSS2),
                           ["b1", "b2", "g1", "g2"]),
                          (BOX1 + (BOX2 + GAUSS1) + GAUSS2,
                           ["b1", "b2", "g1", "g2"]),
                          (BOX1 + ((BOX2 + GAUSS1) + GAUSS2),
                           ["b1", "b2", "g1", "g2"]),
                          (BOX1 + (BOX2 + GAUSS1) + GAUSS2,
                           ["b1", "b2", "g1", "g2"]),
                          # - addition and multiplication
                          ((BOX1 + BOX2) * (GAUSS1 + GAUSS2),
                           ["b1 * g1", "b1 * g2", "b2 * g1", "b2 * g2"]),
                          ((BOX1 - BOX2) * (GAUSS1 - GAUSS2),
                           ["b1 * g1", "b1 * -(g2)", "-(b2) * g1", "-(b2) * -(g2)"]),
                          ((BOX1 * BOX2) - (GAUSS1 * GAUSS2),
                           ["b1 * b2", "-(g1 * g2)"]),
                          # - division
                          (GAUSS1 * GAUSS2 / (BOX1 * BOX2), ["g1 * g2 / (b1 * b2)"]),
                          (GAUSS1 * GAUSS2 // (BOX1 * BOX2), ["(g1 * g2) // (b1 * b2)"]),
                          ((GAUSS1 + GAUSS2) / (BOX1 + BOX2), ["g1 / (b1 + b2)", "g2 / (b1 + b2)"]),
                          ((GAUSS1 + GAUSS2) // (BOX1 + BOX2), ["(g1 + g2) // (b1 + b2)"]),
                          ((GAUSS1 - GAUSS2) / (BOX1 + BOX2), ["g1 / (b1 + b2)", "-(g2) / (b1 + b2)"]),
                          ((GAUSS1 - GAUSS2) // (BOX1 + BOX2), ["(g1 - g2) // (b1 + b2)"]),
                          # - other combinations
                          (GAUSS1 * GAUSS2 ** (BOX1 * BOX2), ["g1 * g2 ** (b1 * b2)"]),
                          (GAUSS1 * GAUSS2 % (BOX1 * BOX2), ["(g1 * g2) % (b1 * b2)"]),
                          # - more realistic options
                          (BOX1 * (BOX2 * GAUSS1 + GAUSS2) * BOX1,
                           ["b1 * b2 * g1 * b1", "b1 * g2 * b1"]),
                          (BOX1 * (BOX2 * GAUSS1 - GAUSS2) + BOX1,
                           ["b1 * b2 * g1", "b1 * -(g2)", "b1"]),
                          (BOX1 * BOX2 * (GAUSS1 + BOX2 * GAUSS2),
                           ["b1 * b2 * g1", "b1 * b2 * b2 * g2"]),
                          ((GAUSS1 + BOX2 * GAUSS2) * BOX1 * BOX2,
                           ["g1 * b1 * b2", "b2 * g2 * b1 * b2"]),
                          ((GAUSS1 + GAUSS2 - BOX2) * (BOX1 + BOX2),
                           ["g1 * b1", "g1 * b2", "g2 * b1", "g2 * b2", "-(b2) * b1", "-(b2) * b2"]),
                          ((GAUSS1 + (GAUSS2 - BOX2)) * (BOX1 + BOX2),
                           ["g1 * b1", "g1 * b2", "g2 * b1", "g2 * b2", "-(b2) * b1", "-(b2) * b2"]),
                          ((GAUSS1 + GAUSS2 * GAUSS1 - BOX2) * (BOX1 + BOX2),
                           ["g1 * b1", "g1 * b2", "g2 * g1 * b1", "g2 * g1 * b2", "-(b2) * b1", "-(b2) * b2"]),
                          (((1 + GAUSS2) * GAUSS1 - BOX2) * (BOX1 + BOX2),
                           ["1.0 * g1 * b1", "1.0 * g1 * b2", "g2 * g1 * b1", "g2 * g1 * b2", "-(b2) * b1", "-(b2) * b2"]),
                          ((GAUSS1 + (1 + (GAUSS2 / BOX1))) * (BOX1 + BOX2),
                           ["g1 * b1", "g1 * b2", "1.0 * b1", "1.0 * b2", "g2 / b1 * b1", "g2 / b1 * b2"]),
                          ((GAUSS1 + (1 + (GAUSS2 / (BOX1 + BOX2)))) * (BOX1 + BOX2),
                           ["g1 * b1", "g1 * b2", "1.0 * b1", "1.0 * b2", "g2 / (b1 + b2) * b1", "g2 / (b1 + b2) * b2"]),
                          (((GAUSS1 + GAUSS2 - 4) / (BOX1 + BOX2) * SCALE1),
                           ["g1 / (b1 + b2) * s1",
                            "g2 / (b1 + b2) * s1",
                            "-(4.0) / (b1 + b2) * s1"]),
                           ((BOX1 + BOX2) ** 2 / 2 + GAUSS1,
                            ["(b1 + b2) ** 2.0 / 2.0",
                             "g1"]),
                           (GAUSS1 + (BOX1 - BOX2) ** 2 / 2,
                            ["g1",
                             "(b1 - b2) ** 2.0 / 2.0"])
                          ])
def test_model_deconstruct(model, expecteds):
    """Check model deconstruction using the model name and evaluation.

    This uses a set of known models in the assumption we have covered
    all the relevant code paths. Particular care is needed to ensure
    conditions like m1 / (m2 + m3) are included. At some point this
    becomes a regression test, as it's not neessarily important that
    the best possible deconstruction is created, just that we get
    consistent results (e.g. given that the code does not simplify
    expressions, in particular the handling of negation).

    The evaluation test checks that

        model(x) = sum_i term_i(x)

    for all the terms that model_deconstruct creates. The idea is that
    the models are set to non-zero values over the range they use (ie
    within -10 to 10 for the box components), and that there are some
    components which vary with x just to check everything passes
    through correctly.

    """

    terms = model_deconstruct(model)
    assert len(terms) == len(expecteds)
    for term, expected in zip(terms, expecteds):
        assert term.name == expected

    x = np.arange(-9, 9, 1)
    got_model = model(x)
    got_terms = np.zeros_like(got_model)
    for term in terms:
        got_terms += term(x)

    assert got_terms == pytest.approx(got_model)


@pytest.mark.parametrize("ntotal", [101, 1001])
def test_model_deconstruct_possible_recursion_error_lhs(ntotal):
    """Try and test the recursion-handling on the LHS.

    It is not obvious if the recursion-handling is always going to
    trigger when 1000 frames are hit, so we try a value which is
    known - for Python 3.11 - to trigger a recursion error, but it
    may not with other Python cases.  We also have a run where we
    know we don't hit the recursion error so we can check the
    fall-over path.

    We could error out if the 1001 case does not trigger a recursion
    error, but it seems poor form for the tests to fail if Python
    changes things in the future when Sherpa may not be being
    developed any more.

    """

    def mk(n):
        return basic.Scale1D(f"c{n}")

    mdl = mk(0)
    for i in range(1, ntotal):
        mdl += mk(i)

    # mdl is "(((...(c0 + c1) + .. ) + c999) + c1000)" where the
    # brackets are left in to point out the structure, when
    # ntotal=1001.
    #
    cpts = model_deconstruct(mdl)

    # Did we create ntotal individual components?
    #
    ncpts = len(cpts)
    if ncpts == ntotal:
        for cpt in cpts:
            assert isinstance(cpt, basic.Scale1D)

        return

    # Just check we aren't creating too-many components.
    #
    assert ncpts < ntotal

    # Check we have the expected break down.
    #
    assert isinstance(cpts[0], BinaryOpModel)
    for cpt in cpts[1:]:
        assert isinstance(cpt, basic.Scale1D)


@pytest.mark.parametrize("ntotal", [101, 1001])
def test_model_deconstruct_possible_recursion_error_rhs(ntotal):
    """Try and test the recursion-handling on the LHS.

    It is not obvious if the recursion-handling is always going to
    trigger when 1000 frames are hit, so we try a value which is
    known - for Python 3.11 - to trigger a recursion error, but it
    may not with other Python cases. We also have a run where we
    know we don't hit the recursion error so we can check the
    fall-over path.

    We could error out if the 1001 case does not trigger a recursion
    error, but it seems poor form for the tests to fail if Python
    changes things in the future when Sherpa may not be being
    developed any more.

    """

    def mk(n):
        return basic.Scale1D(f"c{n}")

    mdl = mk(0)
    for i in range(1, ntotal):
        mdl = BinaryOpModel(mk(i), mdl, np.add, "+")

    # mdl is "(c1000 + (c999 + (... + (c1 + c0)...)))" where the
    # brackets are left in to point out the structure, when
    # ntotal=1001.
    #
    cpts = model_deconstruct(mdl)

    # Did we create ntotal individual components?
    #
    ncpts = len(cpts)
    if ncpts == ntotal:
        for cpt in cpts:
            assert isinstance(cpt, basic.Scale1D)

        return

    # Just check we aren't creating too-many components.
    #
    assert ncpts < ntotal

    # Check we have the expected break down.
    #
    for cpt in cpts[:-1]:
        assert isinstance(cpt, basic.Scale1D)

    assert isinstance(cpts[-1], BinaryOpModel)
