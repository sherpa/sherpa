#
#  Copyright (C) 2007, 2016, 2018, 2020 - 2024
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

import numpy as np

import pytest

from sherpa.models import basic
from sherpa.models.parameter import hugeval
from sherpa.models.model import ArithmeticModel, RegriddableModel1D, \
    RegriddableModel2D
from sherpa.utils.err import ModelErr
from sherpa.utils.numeric_types import SherpaFloat


def userfunc(pars, x, *args, **kwargs):
    return x


EXCLUDED_MODELS = (ArithmeticModel, RegriddableModel1D, RegriddableModel2D,
                   basic.Const)

TESTABLE = []
for name in dir(basic):
    cls = getattr(basic, name)
    if not isinstance(cls, type) or \
       not issubclass(cls, ArithmeticModel) or \
       cls in EXCLUDED_MODELS:
        continue

    TESTABLE.append((name, cls))


def test_expected_number():
    assert len(TESTABLE) == 34


@pytest.mark.parametrize("name,cls", TESTABLE)
def test_create_and_evaluate(name, cls):

    # These have very different interfaces than the others
    if name in ['Integrator1D', 'Integrate1D']:
        return

    x = np.arange(1.0, 5.0)

    m = cls()
    assert type(m).__name__.lower() == m.name

    if isinstance(m, basic.TableModel):
        m.load(x, x)
    if isinstance(m, basic.UserModel):
        m.calc = userfunc

    if m.ndim == 2:
        pt_out = m(x, x)
        int_out = m(x, x, x, x)
    else:
        if m.name in ('log', 'log10'):
            xx = -x
        else:
            xx = x
        pt_out = m(xx)
        int_out = m(xx, xx)

    for out in (pt_out, int_out):
        assert out.dtype.type is SherpaFloat
        assert out.shape == x.shape


def test_polynom2d_guess():
    """Check guess works okay: polynom2d"""

    mdl = basic.Polynom2D()

    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(1.0)
        else:
            assert p.val == pytest.approx(0.0)

        assert p.min == pytest.approx(-hugeval)
        assert p.max == pytest.approx(hugeval)

    x = np.asarray([1, 2, 30])
    y = np.asarray([-2, 4, -2])
    z = np.asarray([10, 12, -2])

    mdl.guess(z, x, y)

    # This is a regression test
    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(5.0)
        else:
            assert p.val == pytest.approx(0.0)

    for par in [mdl.c, mdl.cx1y2, mdl.cx2y1, mdl.cx2y2]:
        assert par.min == pytest.approx(-2)
        assert par.max == pytest.approx(12)

    def check(par, maxval):
        assert par.min == pytest.approx(-maxval)
        assert par.max == pytest.approx(maxval)

    check(mdl.cy1, 233.33333333333334)
    check(mdl.cy2, 38.88888888888889)
    check(mdl.cx1, 48.275862068965516)
    check(mdl.cx2, 1.6646848989298455)
    check(mdl.cx1y1, 8.045977011494253)


def test_polynom2d_guess_ymin0():
    """Check guess works okay: polynom2d

    What happens when ymin=0
    """

    mdl = basic.Polynom2D()

    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(1.0)
        else:
            assert p.val == pytest.approx(0.0)

        assert p.min == pytest.approx(-hugeval)
        assert p.max == pytest.approx(hugeval)

    x = np.asarray([1, 2, 30])
    y = np.asarray([-2, 4, -2])
    z = np.asarray([10, 12, 0])

    mdl.guess(z, x, y)

    # This is a regression test
    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(6.0)
        else:
            assert p.val == pytest.approx(0.0)

    for par in [mdl.c, mdl.cx1y2, mdl.cx2y1, mdl.cx2y2]:
        assert par.min == pytest.approx(0)
        assert par.max == pytest.approx(12)

    def check(par, maxval):
        assert par.min == pytest.approx(-maxval)
        assert par.max == pytest.approx(maxval)

    check(mdl.cy1, 200)
    check(mdl.cy2, 33.33333333333333)
    check(mdl.cx1, 41.37931034482759)
    check(mdl.cx2, 1.426872770511296)
    check(mdl.cx1y1, 6.896551724137931)


def test_polynom2d_guess_ymax0():
    """Check guess works okay: polynom2d

    What happens when ymax=0
    """

    mdl = basic.Polynom2D()

    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(1.0)
        else:
            assert p.val == pytest.approx(0.0)

        assert p.min == pytest.approx(-hugeval)
        assert p.max == pytest.approx(hugeval)

    x = np.asarray([1, 2, 30])
    y = np.asarray([-2, 4, -2])
    z = np.asarray([-10, -12, 0])

    mdl.guess(z, x, y)

    # This is a regression test
    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(-6.0)
        else:
            assert p.val == pytest.approx(0.0)

    for par in [mdl.c, mdl.cx1y2, mdl.cx2y1, mdl.cx2y2]:
        assert par.min == pytest.approx(-12)
        assert par.max == pytest.approx(0)

    def check(par, maxval):
        assert par.min == pytest.approx(-maxval)
        assert par.max == pytest.approx(maxval)

    check(mdl.cy1, 200)
    check(mdl.cy2, 33.33333333333333)
    check(mdl.cx1, 41.37931034482759)
    check(mdl.cx2, 1.426872770511296)
    check(mdl.cx1y1, 6.896551724137931)


def test_polynom2d_guess_y0():
    """Check guess works okay: polynom2d

    What happens when y=0
    """

    mdl = basic.Polynom2D()

    for p in mdl.pars:
        if p.name == 'c':
            assert p.val == pytest.approx(1.0)
        else:
            assert p.val == pytest.approx(0.0)

        assert p.min == pytest.approx(-hugeval)
        assert p.max == pytest.approx(hugeval)

    x = np.asarray([1, 2, 30])
    y = np.asarray([-2, 4, -2])
    z = np.asarray([0, 0, 0])

    mdl.guess(z, x, y)

    # This is a regression test
    for p in mdl.pars:
        assert p.val == pytest.approx(0.0)

    def check(par, maxval):
        assert par.min == pytest.approx(-maxval)
        assert par.max == pytest.approx(maxval)

    for par in mdl.pars:
        check(par, 0)


def test_tablemodel_invalid_x():
    """We error out if X is invalid."""

    tbl = basic.TableModel()
    with pytest.raises(ModelErr,
                       match="Unable to treat array as numeric"):
        tbl.load([1, 2, "x"], [3, 4, 5])


def test_tablemodel_invalid_y():
    """We error out if y is invalid."""

    tbl = basic.TableModel()
    with pytest.raises(ModelErr,
                       match="Unable to treat array as numeric"):
        tbl.load([3, 4, 5], [1, 2, "x"])
