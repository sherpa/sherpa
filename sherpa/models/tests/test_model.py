#
#  Copyright (C) 2007, 2016, 2017, 2020  Smithsonian Astrophysical Observatory
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

from collections import namedtuple
import logging
import operator
import warnings

import numpy

import pytest

from sherpa.utils.err import ModelErr
from sherpa.models.model import ArithmeticModel, ArithmeticConstantModel, \
    BinaryOpModel, FilterModel, NestedModel, UnaryOpModel
from sherpa.models.parameter import Parameter, tinyval
from sherpa.models.basic import Sin, Const1D


def validate_warning(warning_capturer, parameter_name="norm", model_name="ParameterCase", actual_name="Ampl", num=1):
    assert num == len(warning_capturer)
    for warning in warning_capturer:
        assert issubclass(warning.category, DeprecationWarning)
        expected_warning_message = 'Parameter name {} is deprecated for model {}, use {} instead'.format(
            parameter_name, model_name, actual_name
        )
        assert expected_warning_message == str(warning.message)


def my_sin(pars, x):
    return (pars[2].val *
            numpy.sin(2.0 * numpy.pi * (x - pars[1].val) / pars[0].val))


# Test support for renamed parameters by sub-classing
# the Sin model (which lets the tests be re-used).
#
class RenamedPars(Sin):
    # The only reason I am extending Sin is to inherit the method
    # implementations

    def __init__(self, name='renamedpars'):
        self.period = Parameter(name, 'period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'ampl', 1, 1e-05, hard_min=0, aliases=['norm'])
        ArithmeticModel.__init__(self, name,
                                 (self.period, self.offset, self.ampl))


# During testing an XSPEC model was found to refer to one of its
# parameters with a different case to how it was created, which
# caused problems. Replicate this problem case here.
#
class ParameterCase(ArithmeticModel):
    """Re-implemenent Sin model so can copy tests"""

    def __init__(self, name='parametercase'):
        self.period = Parameter(name, 'Period', 1, 1e-10, 10, tinyval)
        self.offset = Parameter(name, 'Offset', 0, 0, hard_min=0)
        self.ampl = Parameter(name, 'Ampl', 1, 1e-05, hard_min=0, aliases=["NORM"])

        with warnings.catch_warnings(record=True) as warn:
            warnings.simplefilter("always", DeprecationWarning)
            pars = (self.perioD, self.oFFSEt, self.NORM)
            validate_warning(warn)

        self._basemodel = Sin()
        ArithmeticModel.__init__(self, name, pars)

    def calc(self, *args, **kwargs):
        for par in self.pars:
            setattr(self._basemodel, par.name, par.val)

        self._basemodel.integrate = self.integrate
        return self._basemodel.calc(*args, **kwargs)


def setup_model():
    return Sin('m'), ['period', 'offset', 'ampl']


def setup_renamed():
    return RenamedPars('m'), ['period', 'offset', 'ampl']


def setup_parametercase():
    return ParameterCase(), ['Period', 'Offset', 'Ampl']


def setup_composite():
    out = namedtuple('composite', ['m', 'm2', 's', 'x', 'xx'])
    out.m  = Const1D('m')
    out.m2 = Const1D('m2')
    out.m.c0  = 2
    out.m2.c0 = 4
    out.s = Sin('s')
    out.x = 1.0
    out.xx = numpy.arange(-10.0, 10.0)
    return out


@pytest.mark.parametrize('setup,name',
                         [(setup_model, 'm'),
                          (setup_renamed, 'm'),
                          (setup_parametercase, 'parametercase')])
def test_name(setup, name):
    mdl, _ = setup()
    assert mdl.name == name


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_iter(setup):
    mdl, _ = setup()
    for part in mdl:
        assert part is mdl


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_num_pars(setup):
    mdl, _ = setup()
    assert len(mdl.pars) == 3


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_par_names(setup):
    mdl, names = setup()
    assert len(names) == len(mdl.pars)
    for name, par in zip(names, mdl.pars):
        assert par.name == name


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_getpar(setup):
    mdl, _ = setup()
    for par in (mdl.period, mdl.PerioD, mdl.PERIod):
        assert par is mdl.pars[0]

    with pytest.raises(AttributeError):
        mdl.perio

    with pytest.raises(AttributeError):
        getattr(mdl, 'perio')


@pytest.mark.parametrize('setup,warns',
                         [(setup_renamed, ["norm", "RenamedPars", "ampl"]),
                          (setup_parametercase, [])])
def test_getpar_rename(setup, warns):
    mdl, _ = setup()

    with warnings.catch_warnings(record=True) as warn:
        warnings.simplefilter("always", DeprecationWarning)

        for par in (mdl.norm, mdl.NorM, mdl.NOrm):
            assert par is mdl.pars[2]

        validate_warning(warn, *warns, num=3)


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_setpar(setup):
    mdl, _ = setup()

    assert mdl.offset.val != 17.0

    mdl.offset = 17
    assert mdl.offset.val == 17.0

    mdl.ofFseT = 18
    assert mdl.offset.val == 18.0

    with pytest.raises(AttributeError):
        mdl.ofset = 19

    with pytest.raises(AttributeError):
        setattr(mdl, 'ofset', 19)


@pytest.mark.parametrize('setup,warns',
                         [(setup_renamed, ["norm", "RenamedPars", "ampl"]),
                          (setup_parametercase, [])])
def test_setpar_rename(setup, warns):
    mdl, _ = setup()

    mdl.ampl = 1
    assert mdl.ampl.val != 12.0

    with warnings.catch_warnings(record=True) as warn:
        warnings.simplefilter("always", DeprecationWarning)

        mdl.norm = 12
        validate_warning(warn, *warns)

    assert mdl.ampl.val == 12.0

    with warnings.catch_warnings(record=True) as warn:
        warnings.simplefilter("always", DeprecationWarning)

        mdl.NoRM = 18
        validate_warning(warn, *warns)

    assert mdl.ampl.val == 18.0
    mdl.ampl = 1
    assert mdl.ampl.val == 1.0


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_calc_and_call(setup):
    mdl, _ = setup()

    x = numpy.arange(10.0)
    refvals = my_sin(mdl.pars, x)

    pars = [p.val for p in mdl.pars]

    y1 = mdl.calc(pars, x)
    assert y1 == pytest.approx(refvals, rel=1e-12)

    y2 = mdl(x)
    assert y2 == pytest.approx(refvals, rel=1e-12)


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_get_thawed_pars(setup):
    mdl, _names = setup()

    tp = [p.val for p in mdl.pars if not p.frozen]
    assert len(mdl.thawedpars) == len(tp)
    for a, b in zip(mdl.thawedpars, tp):
        assert a == b


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_set_thawed_pars(setup):
    mdl, _ = setup()
    pars = [7, 8, 9]

    assert pars != mdl.thawedpars

    mdl.thawedpars = pars

    assert pars == mdl.thawedpars

    with pytest.raises(ModelErr):
        mdl.thawedpars = [7, 8]

    with pytest.raises(ValueError):
        mdl.thawedpars = [1, 2, 'ham']

    pars[0] = mdl.pars[0].hard_min / 10
    pars[1] = mdl.pars[1].hard_max * 10
    pars[2] = 8

    logger = logging.getLogger('sherpa')
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)

    try:
        mdl.thawedpars = pars
    finally:
        logger.setLevel(old_level)

    assert mdl.pars[0].val == mdl.pars[0].min
    assert mdl.pars[1].val == mdl.pars[1].max
    assert mdl.pars[2].val == 8


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_get_mins_maxes(setup):
    mdl, _ = setup()

    assert mdl.thawedparmins == [p.min for p in mdl.pars if not p.frozen]
    assert mdl.thawedparmaxes == [p.max for p in mdl.pars if not p.frozen]

    assert mdl.thawedparhardmins == [p.hard_min for p in mdl.pars if not p.frozen]
    assert mdl.thawedparhardmaxes == [p.hard_max for p in mdl.pars if not p.frozen]


def test_composite_iter():
    out = setup_composite()
    m = 3 * out.m + out.m2
    parts = list(m)

    assert type(parts[0]) is BinaryOpModel
    assert type(parts[1]) is ArithmeticConstantModel
    assert parts[1].val == 3.0
    assert parts[2] is out.m
    assert parts[3] is out.m2
    assert len(parts) == 4


def test_composite_unop():
    out = setup_composite()
    for op in (abs, operator.neg):
        m = op(out.m)

        assert isinstance(m, UnaryOpModel)

        y1 = m(out.x)
        y2 = op(out.m(out.x))
        assert y1 == y2


def test_composite_binop():
    out = setup_composite()
    ops = [operator.add, operator.sub, operator.mul,
           operator.floordiv, operator.truediv, operator.mod,
           operator.pow]

    for op in ops:
        for m in (op(out.m, out.m2.c0.val),
                  op(out.m.c0.val, out.m2),
                  op(out.m, out.m2)):

            assert isinstance(m, BinaryOpModel)

            y1 = m(out.x)
            y2 = op(out.m.c0.val, out.m2.c0.val)
            assert y1 == y2


def test_composite_complex_expression():
    out = setup_composite()
    cmplx = (3 * out.m + out.m2) / (out.m ** 3.2)
    m = out.m(out.x)
    m2 = out.m2(out.x)

    y1 = cmplx(out.x)
    y2 = (3 * m + m2) / (m ** 3.2)
    assert y1 == y2


def test_filter():
    out = setup_composite()
    m = out.s[::2]

    assert type(m) is FilterModel

    assert numpy.all(m(out.xx) == out.s(out.xx)[::2])


def test_nested():
    out = setup_composite()

    for func in (numpy.sin, out.s):
        m = out.m.apply(func)

        assert type(m) is NestedModel

        y1 = m(out.x)
        y2 = func(out.m(out.x))
        assert y1 == y2
