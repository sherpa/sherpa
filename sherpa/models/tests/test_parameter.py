#
#  Copyright (C) 2007, 2016, 2018, 2019, 2020, 2021
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

from numpy import arange

import pytest

from sherpa.utils import SherpaFloat
from sherpa.models.parameter import Parameter, UnaryOpParameter, \
    BinaryOpParameter, ConstantParameter, hugeval
from sherpa.utils.err import ParameterErr
from sherpa.models.basic import Gauss1D, Const1D, PowLaw1D
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
        p.frozen = arange(10)

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


@pytest.mark.parametrize("op", (abs, operator.neg))
def test_unop(op):
    porig = setUp_p()
    p = op(porig)
    assert isinstance(p, UnaryOpParameter)
    assert p.val == op(p.val)


@pytest.mark.parametrize("op", [operator.add, operator.sub, operator.mul,
                                operator.floordiv, operator.truediv, operator.mod,
                                operator.pow])
def test_binop(op):
    p, p2 = setUp_composite()
    for comb in (op(p, p2.val), op(p.val, p2),
                 op(p, p2)):
        assert isinstance(comb, BinaryOpParameter)
        assert comb.val == op(p.val, p2.val)


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

    grid = arange(1, 5)

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
