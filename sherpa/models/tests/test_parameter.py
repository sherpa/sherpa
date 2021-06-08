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


# The p/afp names are from the original version of the tests.
#
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

    # with pytest.raises(TypeError):
    #     p.frozen = arange(10)
    #
    # At present the following raises a ValueError about the
    # truth value of an array with more than one element because
    # of the attempt to call bool(arange(10)). In the original
    # version of the test there was an error in the test that
    # caused an unrelated TypeError to be raised, making it look
    # like it had passed.
    #
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


class ParBase:

    def __init__(self, pos1, pos2):
        self.src1 = Gauss1D()
        self.src1.pos = pos1
        self.src1_pos = pos1
        self.src2 = Gauss1D()
        self.src2.pos = pos2
        self.src2_pos = pos2
        self.tst_pos(self.src1.pos, self.src1_pos)
        self.tst_pos(self.src2.pos, self.src2_pos)

    def tst_unlink(self):
        self.src1.pos.unlink()
        self.tst_pos(self.src1.pos, self.src1.pos.default_val,
                     min=self.src1.pos.min, max=self.src1.pos.max)
        self.src2.pos.unlink()
        self.tst_pos(self.src2.pos, self.src2_pos, min=self.src2.pos.min,
                     max=self.src2.pos.max)

    def tst_pos(self, gauss, pos, min=-hugeval, max=hugeval, frozen=False,
                link=None, link_min=None, link_max=None):
        assert gauss.val == pos
        assert gauss.min == min
        assert gauss.max == max
        assert gauss.frozen == frozen
        if gauss.link is None or link is None:
            assert gauss.link == link
        else:
            self.tst_pos(gauss.link, pos)
        assert gauss.default_val == pos
        assert gauss.default_min == min
        assert gauss.default_max == max


class ParVal(ParBase):

    def __init__(self, pos1, pos2):
        ParBase.__init__(self, pos1, pos2)

    def tst(self):
        self.tst_pos(self.src1.pos, self.src2_pos, frozen=True,
                     link=self.src1.pos)
        self.tst_pos(self.src2.pos, self.src2_pos)

    def tst_low_level_val_link(self):
        self.src1.pos.val = self.src2.pos.val
        self.src1.pos.link = self.src2.pos
        self.tst()
        self.tst_unlink()

    def tst_ui_val_link(self):
        ui.link(self.src1.pos, self.src2.pos)
        self.tst()
        self.tst_unlink()


def test_link_unlink_val():
    tst = ParVal(4, 5)
    tst.tst_unlink()


def test_link_unlink_val_low_level():
    tst = ParVal(4, 5)
    tst.tst_low_level_val_link()


def test_link_unlink_val_ui():
    tst = ParVal(4, 5)
    tst.tst_ui_val_link()


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
