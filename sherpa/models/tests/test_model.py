#
#  Copyright (C) 2007, 2016, 2017, 2020, 2021
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

from collections import namedtuple
import logging
import operator
import warnings

# Repeat the logic from sherpa/models/model.py
#
try:
    from hashlib import md5 as hashfunc
except ImportError:
    from hashlib import sha256 as hashfunc

import numpy

import pytest

from sherpa.utils.err import ModelErr
from sherpa.models.model import ArithmeticModel, ArithmeticConstantModel, \
    ArithmeticFunctionModel, BinaryOpModel, CompositeModel, FilterModel, \
    Model, NestedModel, UnaryOpModel, RegridWrappedModel, modelCacher1d, \
    expand, separate
from sherpa.models.parameter import Parameter, hugeval, tinyval
from sherpa.models.basic import Sin, Const1D, Box1D, Gauss1D, Polynom1D, \
    Scale1D, Integrate1D


def validate_warning(warning_capturer, parameter_name="norm",
                     model_name="ParameterCase", actual_name="Ampl", num=1):
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


class ShowableModel(ArithmeticModel):
    """Test out hidden parameter"""

    def __init__(self, name='thetestmodel'):
        self.notseen = Parameter(name, 'notseen', 10, min=5, max=15, hidden=True)
        self.p1 = Parameter(name, 'P1', 1, min=0, max=10)
        self.accent = Parameter(name, 'accent', 2, min=0, max=10, frozen=True)
        self.norm = Parameter(name, 'norm', 100, min=0, max=1e4)

        ArithmeticModel.__init__(self, name, (self.notseen, self.p1, self.accent, self.norm))


def setup_model():
    return Sin('m'), ['period', 'offset', 'ampl']


def setup_renamed():
    return RenamedPars('m'), ['period', 'offset', 'ampl']


def setup_parametercase():
    return ParameterCase(), ['Period', 'Offset', 'Ampl']


def setup_composite():
    out = namedtuple('composite', ['m', 'm2', 's', 'x', 'xx'])
    out.m = Const1D('m')
    out.m2 = Const1D('m2')
    out.m.c0 = 2
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


@pytest.mark.parametrize('setup',
                         [setup_model, setup_renamed, setup_parametercase])
def test_get_mins_maxes_errors(setup):
    """Check that we get errors when setting values incorrectly"""
    mdl, _ = setup()

    with pytest.raises(ModelErr):
        mdl.thawedparmins = [-5, -7]

    with pytest.raises(ModelErr):
        mdl.thawedparmaxes = [2, 3, 4, 5]


def test_get_mins_maxes_warns(caplog):
    """What happens if min/max set outside valid range"""

    mdl = Box1D()
    mdl.xlow.set(min=0, max=100, val=50)
    mdl.xhi.set(min=200, max=300, val=250)

    # Some tests can change the logging level, so make sure
    # we have the level set correctly.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        mdl.thawedparmins = [10, -1e40, -10]
        mdl.thawedparmaxes = [1e40, 290, 10]

    assert mdl.thawedparmins == [10, -hugeval, -10]
    assert mdl.thawedparmaxes == [hugeval, 290, 10]

    # The handling of setting min to greater > hard_max
    # (and vice versa) isn't as easily checked, as the
    # parameter value ends up causing the code to
    # error out, so we set the parameter value to the
    # limit beforehand.
    #
    mdl.xhi.set(max=hugeval, val=hugeval)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        mdl.thawedparmins = [10, 1e40, -10]
        mdl.thawedparmaxes = [1000, 1e41, 10]

    mdl.xlow.set(min=-hugeval, val=-hugeval)

    with caplog.at_level(logging.INFO, logger='sherpa'):
        mdl.thawedparmins = [-1e41, hugeval, -10]
        mdl.thawedparmaxes = [-1e40, hugeval, 10]

    assert len(caplog.records) == 6
    msgs = []
    for (log_name, log_level, message) in caplog.record_tuples:
        assert log_level == logging.WARNING
        assert log_name == 'sherpa.models.model'
        msgs.append(message)

    assert msgs[0] == 'value of parameter box1d.xhi minimum is below hard minimum; setting to hard minimum'
    assert msgs[1] == 'value of parameter box1d.xlow maximum is above hard maximum; setting to hard maximum'

    assert msgs[2] == 'value of parameter box1d.xhi minimum is above hard maximum; setting to hard maximum'
    assert msgs[3] == 'value of parameter box1d.xhi maximum is above hard maximum; setting to hard maximum'

    assert msgs[4] == 'value of parameter box1d.xlow minimum is below hard minimum; setting to hard minimum'
    assert msgs[5] == 'value of parameter box1d.xlow maximum is below hard minimum; setting to hard minimum'


def test_reset():
    """The reset method restores the last-changed values, so is a no-op here"""

    mdl, _ = setup_model()

    defvals = [p.val for p in mdl.pars]
    assert not mdl.ampl.frozen

    mdl.offset = 10
    mdl.ampl.freeze()

    newvals = [p.val for p in mdl.pars]
    assert newvals != defvals

    # Since this test does not adjust the parameter values using a low-level
    # routine it has no effect.
    #
    mdl.reset()

    vals = [p.val for p in mdl.pars]
    assert vals == newvals
    assert mdl.ampl.frozen


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


def test_show_model():
    """Can we show the model?"""

    m = ShowableModel()
    n = Const1D()
    n.c0 = 2

    # add a test of parameter linking
    m.norm = 10 + n.c0

    toks = str(m).split('\n')
    assert len(toks) == 6
    assert toks[0] == 'thetestmodel'
    assert toks[1].strip().split() == ['Param', 'Type', 'Value', 'Min', 'Max', 'Units']
    assert toks[2].strip().split() == ['-----', '----', '-----', '---', '---', '-----']
    assert toks[3].strip().split() == ['thetestmodel.P1', 'thawed', '1', '0', '10']
    assert toks[4].strip().split() == ['thetestmodel.accent', 'frozen', '2', '0', '10']
    assert toks[5].strip().split() == ['thetestmodel.norm', 'linked', '12',
                                       'expr:', '(10', '+', 'const1d.c0)']

    # Should hidden parameters be capable of being thawed?
    assert m.notseen.val == 10
    assert m.notseen.min == 5
    assert m.notseen.max == 15
    assert not m.notseen.frozen


def test_binop_checks_sizes():
    """We error out if the sides generate different sizes.

    It requires a bit of effort to trigger this behavior.
    """

    m1 = Const1D()
    m1.c0 = 2
    m2 = numpy.asarray([2, 3, 4])
    m = m1 + m2
    assert isinstance(m, BinaryOpModel)

    # check we can evaluate with 3 elements:
    ans = m([100, 200, 300])
    assert ans == pytest.approx([4, 5, 6])

    # If we use a different number we error out.
    #
    with pytest.raises(ValueError) as exc:
        m([10, 20])

    emsg = "shape mismatch between 'Const1D: 2' and 'ArithmeticConstantModel: 3'"
    assert str(exc.value) == emsg


def test_functionmodel_check():
    """ArithmeticFunctionModel errors out with invalid input"""

    with pytest.raises(ModelErr) as exc:
        ArithmeticFunctionModel(Const1D())

    assert str(exc.value) == 'ArithmeticFunctionModel instance cannot be created from another model'

    with pytest.raises(ModelErr) as exc:
        ArithmeticFunctionModel(23)

    assert str(exc.value) == 'attempted to create ArithmeticFunctionModel from non-callable object of type int'


@pytest.mark.parametrize('model,mtype',
                         [(ArithmeticConstantModel(23, name='the-23'),
                           ArithmeticConstantModel),
                          (ArithmeticFunctionModel(numpy.sin),
                           ArithmeticFunctionModel)])
def test_unop_arithmeticxxx(model, mtype):
    """Can we apply a function to an Arithmetic*Model object?

    Unlike the BinOpModel test we can't rely on Python numeric
    types, since there's no slots for __neg__ or __abs__ for
    the ArithmeticXXXModel classes, so we can not just say
    mneg = -<model>.

    """

    mneg = UnaryOpModel(model, numpy.negative, '-')
    assert isinstance(mneg, UnaryOpModel)
    assert mneg.op == numpy.negative
    assert isinstance(mneg.arg, mtype)

    x = numpy.linspace(0.1, 0.5, 5)
    y1 = -1 * model(x)
    y2 = mneg(x)
    assert y2 == pytest.approx(y1)


@pytest.mark.parametrize('model,mtype',
                         [(23, ArithmeticConstantModel),
                          (ArithmeticConstantModel(23, name='the-23'),
                           ArithmeticConstantModel),
                          (ArithmeticFunctionModel(numpy.sin),
                           ArithmeticFunctionModel)])
def test_binop_arithmeticxxx(model, mtype):
    """Can we create and combine Arithmetic*Model objects?"""

    m1 = Const1D()
    m1.c0 = 2

    mleft = model + m1
    assert isinstance(mleft, BinaryOpModel)
    assert mleft.op == numpy.add
    assert len(mleft.parts) == 2
    assert isinstance(mleft.parts[0], mtype)
    assert isinstance(mleft.parts[1], Const1D)

    mright = m1 + model
    assert isinstance(mright, BinaryOpModel)
    assert mright.op == numpy.add
    assert len(mright.parts) == 2
    assert isinstance(mright.parts[0], Const1D)
    assert isinstance(mright.parts[1], mtype)

    x = numpy.linspace(0.1, 0.5, 5)
    if mtype == ArithmeticConstantModel:
        y1 = 25 * numpy.ones(5)
    else:
        y1 = model(x) + 2

    y2 = mleft(x)
    y3 = mright(x)
    assert y2 == pytest.approx(y1)
    assert y3 == pytest.approx(y1)


@pytest.mark.parametrize('value,name,expected',
                         [(23, None, '23.0'),
                          (numpy.asarray([3, 5, -1]), None, 'float64[3]'),
                          (23, '24', '24'),  # this is not a good name
                          (numpy.asarray([3, 5, -1]), 'arrayval', 'arrayval')])
def test_constant_show(value, name, expected):
    """Does the ArithmeticConstantModel convert names as expected?"""

    m = ArithmeticConstantModel(value, name=name)
    assert m.name == expected


def test_integrate1d_show_no_model():
    """Check Integrate1D show"""

    imdl = Integrate1D()
    out = str(imdl).split('\n')
    assert out[0] == 'integrate1d'
    assert out[3].strip().startswith('integrate1d.epsabs frozen ')
    assert out[4].strip().startswith('integrate1d.epsrel frozen ')
    assert out[5].strip().startswith('integrate1d.maxeval frozen ')
    assert len(out) == 6


def test_integrate1d_show_model():
    """Check Integrate1D show when applied to a model

    Note that we do not show the integrate1d settings in this version.
    """

    imdl = Integrate1D()
    bmdl = Scale1D()
    mdl = imdl(bmdl)

    out = str(mdl).split('\n')
    assert out[0] == 'integrate1d(scale1d)'
    assert out[3].strip().startswith('scale1d.c0   thawed ')
    assert len(out) == 4


def test_integrate1d_fail_not_integrated():
    """Check Integrate1D requires an integrated grid"""

    imdl = Integrate1D()
    bmdl = Scale1D()
    mdl = imdl(bmdl)
    with pytest.raises(ModelErr) as exc:
        mdl([1.1, 1.2, 1.3, 1.4])

    assert str(exc.value).startswith('A non-overlapping integrated grid is required ')


def test_integrate1d_basic(caplog):
    """Check Integrate1D works

    There is no documentation on how it's supposed to work, so
    this is more a "this currently works, let's hope it continues
    to do so" approach than a test from first principles.
    """

    imdl = Integrate1D()
    bmdl = Scale1D()
    mdl = imdl(bmdl)
    bmdl.c0 = 4

    xlo = numpy.asarray([1.1, 1.2, 1.4, 1.8, 2.4])
    xhi = numpy.asarray([1.2, 1.3, 1.8, 2.0, 3.0])

    # I don't know what the rule is for creating this warning, but
    # check we see it.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        y = mdl(xlo, xhi)

    expected = 4 * numpy.asarray([0.1, 0.1, 0.4, 0.2, 0.6])
    assert y == pytest.approx(expected)

    assert len(caplog.records) == 1
    name, lvl, msg = caplog.record_tuples[0]
    assert name == 'sherpa.models.basic'
    assert lvl == logging.WARNING
    assert msg.startswith('Gauss-Kronrod integration failed with tolerance ')


def test_integrate1d_basic_epsabs(caplog):
    """Check Integrate1D works

    This time adjust epsabs so that we don't get the Gauss-Kronrod
    warning.
    """

    imdl = Integrate1D()
    bmdl = Scale1D()
    mdl = imdl(bmdl)
    bmdl.c0 = 4
    imdl.epsabs = numpy.finfo(numpy.float32).eps

    xlo = numpy.asarray([1.1, 1.2, 1.4, 1.8, 2.4])
    xhi = numpy.asarray([1.2, 1.3, 1.8, 2.0, 3.0])

    # I don't know what the rule is for creating this warning, but
    # check we see it.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        y = mdl(xlo, xhi)

    expected = 4 * numpy.asarray([0.1, 0.1, 0.4, 0.2, 0.6])
    assert y == pytest.approx(expected)

    assert len(caplog.records) == 0


def check_cache(mdl, expected, x, xhi=None):
    """Check the cache contents.

    We assume only one value is being cached at a time. The
    code matches that in sherpa.models.model.modelCacher1d,
    so all it does is check we are using this method.
    """

    cache = mdl._cache
    assert len(cache) == 1

    pars = [p.val for p in mdl.pars]
    data = [numpy.asarray(pars).tobytes(),
            b'1' if mdl.integrate else b'0',
            x.tobytes()]
    if xhi is not None:
        data.append(xhi.tobytes())

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_evaluate_no_cache1d():
    """Check we can turn off cacheing: 1d"""

    xgrid = numpy.arange(2, 10, 1.5)

    mdl = Polynom1D()
    mdl.integrate = False
    mdl._use_caching = False
    assert len(mdl._cache) == 0

    # Check the default values
    expected = numpy.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.c0 = 5
    mdl.c1 = 2

    expected = 5 + 2 * xgrid
    assert mdl(xgrid) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache1d():
    """Check we run with cacheing on: 1d"""

    xgrid = numpy.arange(2, 10, 1.5)

    mdl = Polynom1D()
    mdl.integrate = False
    mdl._use_caching = True
    assert len(mdl._cache) == 0

    # Check the default values
    expected = numpy.ones(6)
    assert mdl(xgrid) == pytest.approx(expected)
    check_cache(mdl, expected, xgrid)

    mdl.c0 = 5
    mdl.c1 = 2

    expected = 5 + 2 * xgrid
    assert mdl(xgrid) == pytest.approx(expected)
    check_cache(mdl, expected, xgrid)


def test_evaluate_no_cache1dint():
    """Check we can turn off cacheing: 1dint"""

    xgrid = numpy.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl._use_caching = False
    assert len(mdl._cache) == 0

    # Check the default values
    expected = numpy.ones(5) * 1.5
    assert mdl(xlo, xhi) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.c0 = 5
    mdl.c1 = 2

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2
    assert mdl(xlo, xhi) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache1dint():
    """Check we run with cacheing on: 1dint"""

    xgrid = numpy.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl._use_caching = True
    assert len(mdl._cache) == 0

    # Check the default values
    expected = numpy.ones(5) * 1.5
    assert mdl(xlo, xhi) == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)

    mdl.c0 = 5
    mdl.c1 = 2

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2
    assert mdl(xlo, xhi) == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)


def test_evaluate_cache_swap():
    """Check issue #959 when swapping integrate flag caused problems.

    Note that the problem causing #959 is actually tested in
    test_evaluate_cache1dint but it's nice to have this
    separate check in case things change.
    """

    xgrid = numpy.arange(2, 10, 1.5)
    xlo, xhi = xgrid[:-1], xgrid[1:]

    mdl = Polynom1D()
    mdl._use_caching = True

    mdl.c0 = 5
    mdl.c1 = 2

    mdl.integrate = False
    expected = 5 + 2 * xlo

    y1 = mdl(xlo, xhi)
    assert y1 == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)

    mdl.integrate = True

    def xval(x):
        return 5 + 2 * x

    dx = xhi - xlo
    ylo = xval(xlo)
    yhi = xval(xhi)

    expected = dx * (yhi + ylo) / 2

    y2 = mdl(xlo, xhi)
    assert y2 == pytest.approx(expected)
    check_cache(mdl, expected, xlo, xhi)


def test_evaluate_cache_arithmeticconstant():
    """Check we run with cacheing: ArihmeticConstant"""

    mdl = ArithmeticConstantModel(2.3)
    assert not hasattr(mdl, '_use_caching')


def test_evaluate_cache_unaryop():
    """UnaryOp has no cache"""

    mdl = Polynom1D()
    assert hasattr(mdl, '_use_caching')

    fmdl = -mdl
    assert isinstance(fmdl, UnaryOpModel)
    assert not hasattr(fmdl, '_use_caching')


def test_evaluate_cache_binaryop():
    """BinaryOp has no cache"""

    mdl = Polynom1D()
    assert hasattr(mdl, '_use_caching')

    fmdl = mdl + 2
    assert isinstance(fmdl, BinaryOpModel)
    assert not hasattr(fmdl, '_use_caching')


def test_evaluate_cache_regrid1d():
    """How about a regridded model?"""

    mdl = Polynom1D()

    x = numpy.arange(2, 20, 0.5)
    rmdl = mdl.regrid(x)

    assert isinstance(rmdl, RegridWrappedModel)
    assert not hasattr(rmdl, '_use_caching')


class DoNotUseModel(Model):
    """This is only used to get a integrate-free model.

    ArithmeticModel sets a integrate field, so this is derived
    from Model. This requires adding some fields from ArithmeticModel
    to support use with modelCacher1d.
    """

    # We need this for modelCacher1d
    _use_caching = True
    _cache = {}
    _cache_ctr = {'hits': 0, 'misses': 0, 'check': 0}
    _queue = ['']

    @modelCacher1d
    def calc(self, p, *args, **kwargs):
        """p is ignored."""

        return numpy.ones(args[0].size)


def test_cache_integrate_fall_through_no_integrate():
    """Try and test the fall-through of the integrate setting.

    This is a bit contrived.
    """

    mdl = DoNotUseModel('notme')
    x = numpy.asarray([2, 3, 7, 100])
    y = mdl(x)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [numpy.asarray(pars).tobytes(),
            b'0', # not integrated
            x.tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_integrate_fall_through_integrate_true():
    """See also test_cache_integrate_fall_through_no_integrate."""

    mdl = DoNotUseModel('notme')
    x = numpy.asarray([2, 3, 7, 100])
    y = mdl(x, integrate=True)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [numpy.asarray(pars).tobytes(),
            b'1', # integrated
            x.tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_integrate_fall_through_integrate_false():
    """See also test_cache_integrate_fall_through_no_integrate."""

    mdl = DoNotUseModel('notme')
    x = numpy.asarray([2, 3, 7, 100])
    y = mdl(x, integrate=False)

    expected = [1, 1, 1, 1]
    assert y == pytest.approx(expected)

    # A simplified version of check_cache (as we have no
    # integrate setting).
    #
    cache = mdl._cache
    assert len(cache) == 1

    pars = []
    data = [numpy.asarray(pars).tobytes(),
            b'0', # not integrated
            x.tobytes()]

    token = b''.join(data)
    digest = hashfunc(token).digest()
    assert digest in cache
    assert cache[digest] == pytest.approx(expected)


def test_cache_status_single(caplog):
    """Check cache_status for a single model."""

    p = Polynom1D()
    with caplog.at_level(logging.INFO, logger='sherpa'):
        p.cache_status()

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.models.model'
    assert lvl == logging.INFO
    toks = msg.split()
    assert toks[0] == 'polynom1d'
    assert toks[1] == 'size:'
    assert toks[2] == '1'
    assert toks[3] == 'hits:'
    assert toks[4] == '0'
    assert toks[5] == 'misses:'
    assert toks[6] == '0'
    assert toks[7] == 'check:'
    assert toks[8] == '0'
    assert len(toks) == 9


def test_cache_status_multiple(caplog):
    """Check cache_status for a multi-component model.

    Unlike test_cache_syayus_single we also have evaluated the model
    so we can check that the cache status has changed.
    """

    # The model expression includes an ArithmeticConstant model (the
    # term 2) which does not have a cache and so is ignored by
    # cache_status.
    #
    p = Polynom1D()
    b = Box1D()
    c = Const1D()
    mdl = c * (2 * p + b)

    # One model is not cached
    b._use_caching = False

    mdl([0.1, 0.2, 0.3])
    mdl([0.1, 0.2, 0.3])
    mdl([0.1, 0.2, 0.3, 0.4])

    with caplog.at_level(logging.INFO, logger='sherpa'):
        mdl.cache_status()

    assert len(caplog.records) == 3

    tokens = []
    for lname, lvl, msg in caplog.record_tuples:
        assert lname == 'sherpa.models.model'
        assert lvl == logging.INFO
        toks = msg.split()
        assert len(toks) == 9
        assert toks[1] == 'size:'
        assert toks[2] == '1'
        assert toks[3] == 'hits:'
        assert toks[5] == 'misses:'
        assert toks[7] == 'check:'
        assert toks[8] == '3'

        tokens.append(toks)

    toks = tokens[0]
    assert toks[0] == 'const1d'
    assert toks[4] == '1'
    assert toks[6] == '2'

    toks = tokens[1]
    assert toks[0] == 'polynom1d'
    assert toks[4] == '1'
    assert toks[6] == '2'

    toks = tokens[2]
    assert toks[0] == 'box1d'
    assert toks[4] == '0'
    assert toks[6] == '0'


def test_cache_clear_single(caplog):
    """Check cache_clear for a single model."""

    p = Polynom1D()

    # There's no official API for accessing the cache data,
    # so do it directly.
    #
    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    p([1, 2, 3])
    p([1, 2, 3])
    p([1, 2, 3, 4])

    assert len(p._cache) == 1
    assert p._cache_ctr['check'] == 3
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 2

    p.cache_clear()

    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0


def test_cache_clear_multiple(caplog):
    """Check cache_clear for a combined model."""

    p = Polynom1D()
    b = Box1D()
    c = Const1D()
    mdl = c * (p + 2 * b)

    # Ensure one component doesn't use the cache
    c._use_caching = False

    # There's no official API for accessing the cache data,
    # so do it directly.
    #
    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    assert len(b._cache) == 0
    assert b._cache_ctr['check'] == 0
    assert b._cache_ctr['hits'] == 0
    assert b._cache_ctr['misses'] == 0

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 0
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0

    mdl([1, 2, 3])
    mdl([1, 2, 3])
    mdl([1, 2, 3, 4])

    assert len(p._cache) == 1
    assert p._cache_ctr['check'] == 3
    assert p._cache_ctr['hits'] == 1
    assert p._cache_ctr['misses'] == 2

    assert len(b._cache) == 1
    assert b._cache_ctr['check'] == 3
    assert b._cache_ctr['hits'] == 1
    assert b._cache_ctr['misses'] == 2

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 3
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0

    mdl.cache_clear()

    assert len(p._cache) == 0
    assert p._cache_ctr['check'] == 0
    assert p._cache_ctr['hits'] == 0
    assert p._cache_ctr['misses'] == 0

    assert len(b._cache) == 0
    assert b._cache_ctr['check'] == 0
    assert b._cache_ctr['hits'] == 0
    assert b._cache_ctr['misses'] == 0

    assert len(c._cache) == 0
    assert c._cache_ctr['check'] == 0
    assert c._cache_ctr['hits'] == 0
    assert c._cache_ctr['misses'] == 0


@pytest.mark.parametrize("method", ["create", "expand"])
def test_compositemodel_not_implented(method):
    c = CompositeModel('bob', ())
    with pytest.raises(NotImplementedError):
        getattr(c, method)()


def test_unaryopmodel_create():
    b = Box1D()
    b.xlow = 10
    b.xhi = 20
    b.ampl.set(12, max=13)
    c = -b

    x = numpy.arange(5, 25, 2)

    n = c.create(b)

    assert n(x) == pytest.approx(c(x))


def test_unaryopmodel_create_invalid():
    c = -Box1D()
    with pytest.raises(TypeError):
        c.create()

    with pytest.raises(TypeError):
        c.create(c, c)


def test_binaryopmodel_create():
    b1 = Box1D()
    b1.xlow = 10
    b1.xhi = 20
    b1.ampl.set(12, max=13)

    b2 = Box1D()
    b2.xlow = 12
    b2.xhi = 22
    b2.ampl.set(-5, min=-5)

    c = b1 + b2

    x = numpy.arange(5, 25, 2)

    n = c.create(b1, b2)

    assert n(x) == pytest.approx(c(x))


def test_binaryopmodel_create_invalid():
    c = Box1D() + Box1D()
    with pytest.raises(TypeError):
        c.create()

    with pytest.raises(TypeError):
        c.create(c)

    with pytest.raises(TypeError):
        c.create(c, c, c)


BOX = Box1D()
LINE = Gauss1D()
LINE2 = Gauss1D('other')
BGND = Polynom1D()


@pytest.mark.parametrize("expression,expected",
                         [(BOX, 'box1d'),
                          (abs(BOX), 'abs(box1d)'),
                          (BOX + LINE, '(box1d + gauss1d)'),
                          (BOX - LINE, '(box1d - gauss1d)'),
                          (BOX * LINE, '(box1d * gauss1d)'),
                          (BOX / LINE, '(box1d / gauss1d)'),
                          (BOX * (LINE + BGND), '((box1d * gauss1d) + (box1d * polynom1d))'),
                          (BOX * (LINE - BGND), '((box1d * gauss1d) - (box1d * polynom1d))'),
                          ((LINE + BGND) * BOX, '((gauss1d * box1d) + (polynom1d * box1d))'),
                          ((LINE - BGND) * BOX, '((gauss1d * box1d) - (polynom1d * box1d))'),
                          (BOX - (LINE + BGND), '(box1d - (gauss1d + polynom1d))'),
                          (BOX - (LINE * BGND), '(box1d - (gauss1d * polynom1d))'),
                          (BOX * (LINE - (LINE2 + BGND)), '((box1d * gauss1d) - ((box1d * other) + (box1d * polynom1d)))'),
                          ((LINE - (LINE2 + BGND)) * BOX, '((gauss1d * box1d) - ((other * box1d) + (polynom1d * box1d)))'),
                          (BOX * (LINE - (LINE2 // BGND)), '((box1d * gauss1d) - (box1d * (other // polynom1d)))'),
                          ((LINE - (LINE2 // BGND)) * BOX, '((gauss1d * box1d) - ((other // polynom1d) * box1d))'),
                          (BGND + BOX + (LINE + LINE2) * BOX - 2 * BGND, '(((polynom1d + box1d) + ((gauss1d * box1d) + (other * box1d))) - (2.0 * polynom1d))')])
def test_model_expand(expression, expected):

    # Note we compare the name field as easiest thing to do
    answer = expand(expression)
    assert answer.name == expected


@pytest.mark.parametrize("expression,expected",
                         [(BOX, ['box1d']),
                          (abs(BOX), ['abs(box1d)']),
                          (BOX + LINE, ['box1d', 'gauss1d']),
                          (BOX - LINE, ['box1d', '-(gauss1d)']),
                          (BOX * LINE, ['(box1d * gauss1d)']),
                          (BOX / LINE, ['(box1d / gauss1d)']),
                          (BOX * (LINE + BGND), ['(box1d * gauss1d)', '(box1d * polynom1d)']),
                          (BOX * (LINE - BGND), ['(box1d * gauss1d)', '-((box1d * polynom1d))']),
                          ((LINE + BGND) * BOX, ['(gauss1d * box1d)', '(polynom1d * box1d)']),
                          ((LINE - BGND) * BOX, ['(gauss1d * box1d)', '-((polynom1d * box1d))']),
                          (BOX + (LINE + BGND), ['box1d', 'gauss1d', 'polynom1d']),
                          (BOX - (LINE + BGND), ['box1d', '-(gauss1d)', '-(polynom1d)']),
                          (BOX - (LINE - BGND), ['box1d', '-(gauss1d)', '-(-(polynom1d))']),
                          (BOX - (LINE * BGND), ['box1d', '-((gauss1d * polynom1d))']),
                          (BOX * (LINE + (LINE2 + BGND)), ['(box1d * gauss1d)', '(box1d * other)', '(box1d * polynom1d)']),
                          (BOX * (LINE - (LINE2 + BGND)), ['(box1d * gauss1d)', '-((box1d * other))', '-((box1d * polynom1d))']),
                          ((LINE - (LINE2 + BGND)) * BOX, ['(gauss1d * box1d)', '-((other * box1d))', '-((polynom1d * box1d))']),
                          (BOX * (LINE - (LINE2 // BGND)), ['(box1d * gauss1d)', '-((box1d * (other // polynom1d)))']),
                          ((LINE - (LINE2 // BGND)) * BOX, ['(gauss1d * box1d)', '-(((other // polynom1d) * box1d))']),
                          (BGND + BOX + (LINE + LINE2) * BOX - 2 * BGND, ['polynom1d', 'box1d', '(gauss1d * box1d)', '(other * box1d)', '-((2.0 * polynom1d))']),
                          ((BOX + LINE) * (LINE2 + BGND), ['(box1d * other)', '(gauss1d * other)', '(box1d * polynom1d)', '(gauss1d * polynom1d)']),
                          ((BOX + LINE) * (LINE2 - BGND), ['(box1d * other)', '(gauss1d * other)', '-((box1d * polynom1d))', '-((gauss1d * polynom1d))']),
                          ((BOX + LINE) * (LINE2 * BGND), ['(box1d * (other * polynom1d))', '(gauss1d * (other * polynom1d))']),
                          ((BOX - LINE) * (LINE2 * BGND), ['(box1d * (other * polynom1d))', '-((gauss1d * (other * polynom1d)))']),
                          ((BOX * LINE) * (LINE2 + BGND), ['((box1d * gauss1d) * other)', '((box1d * gauss1d) * polynom1d)']),
                          ((BOX * LINE) * (LINE2 - BGND), ['((box1d * gauss1d) * other)', '-(((box1d * gauss1d) * polynom1d))']),
                          ])
def test_model_separate(expression, expected):

    # Note we compare the name field as easiest thing to do
    answer = separate(expand(expression))
    print(expression.name)
    print(answer)
    for ans, eans in zip(answer, expected):
        assert ans.name == eans

    assert len(answer) == len(expected)
