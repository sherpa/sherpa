#
#  Copyright (C) 2007, 2016, 2018, 2020, 2021, 2022, 2023
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

import sherpa.astro.models as models
from sherpa.models.model import ArithmeticModel, RegriddableModel2D, \
    RegriddableModel1D, boolean_to_byte
from sherpa.models import basic
from sherpa.utils.numeric_types import SherpaFloat


# See sherpa/models/tests/test_basic.py
#
EXCLUDED_MODELS = (ArithmeticModel, RegriddableModel1D,
                   RegriddableModel2D)


TESTABLE = []
for name in dir(models):
    cls = getattr(models, name)
    if not isinstance(cls, type) or \
       not issubclass(cls, ArithmeticModel) or \
       cls in EXCLUDED_MODELS:
        continue

    TESTABLE.append((name, cls))


def test_expected_number():
    assert len(TESTABLE) == 21


@pytest.mark.parametrize("name,cls", TESTABLE)
def test_create_and_evaluate(name, cls):

    # This has a very different interface than the others
    if name == 'JDPileup':
        return

    x = np.arange(1.0, 5.0)

    m = cls()
    assert type(m).__name__.lower() == m.name

    if m.name == 'linebroad':
        # for some reason with the default params and this x grid
        # the model will fail
        m.vsini = 1e6

    if m.ndim == 2:
        pt_out = m(x, x)
        int_out = m(x, x, x, x)
    else:
        pt_out = m(x)
        int_out = m(x, x)

    for out in (pt_out, int_out):
        assert out.dtype.type is SherpaFloat
        assert out.shape == x.shape


@pytest.mark.parametrize("cls", [models.Atten,
                                 models.BBody,
                                 models.BBodyFreq,
                                 models.BPL1D,
                                 models.Beta1D,
                                 models.Edge,
                                 models.LineBroad,
                                 models.Lorentz1D,
                                 models.NormBeta1D,
                                 models.PseudoVoigt1D,
                                 models.Schechter,
                                 models.Voigt1D])
def test_send_keyword_1d(cls):
    """What happens if we use an un-supported keyword?"""

    mdl = cls()
    mdl._use_caching = False

    if cls == models.LineBroad:
        mdl.vsini = 1e6

    # Not guaranteed to produce interesting results
    x = [10, 12, 13, 15]
    y1 = mdl(x)
    y2 = mdl(x, not_an_argument=True)
    assert y2 == pytest.approx(y1)


@pytest.mark.parametrize("cls", [models.Beta2D,
                                 models.DeVaucouleurs2D,
                                 models.Disk2D,
                                 models.HubbleReynolds,
                                 models.Lorentz2D,
                                 models.Sersic2D,
                                 models.Shell2D])
def test_send_keyword_2d(cls):
    """What happens if we use an un-supported keyword?"""

    mdl = cls()
    mdl._use_caching = False

    # Not guaranteed to produce interesting results
    x0 =  [10, 12, 13, 15]
    x1 =  [5, 7, 12, 13]
    y1 = mdl(x0, x1)
    y2 = mdl(x0, x1, not_an_argument=True)
    assert y2 == pytest.approx(y1)


@pytest.mark.parametrize("test_input, expected", [
    (True, b'1'),
    (False, b'0'),
    (None, b'0'),
    ("foo", b'0')
])
def test_boolean_to_byte(test_input, expected):
    assert boolean_to_byte(test_input) == expected


def test_voigt():
    """Regression test"""

    x = np.linspace(-0.8, 0.8, 10)
    voigt = models.Voigt1D()
    voigt.fwhm_g = 0.3
    voigt.fwhm_l = 0.3
    voigt_result = voigt(x)

    # it's a symmetric model with this grid
    expected = np.asarray([0.07779156, 0.13213921, 0.27123415, 0.67357158, 1.35326402])
    expected = np.hstack((expected, expected[::-1]))

    assert voigt_result == pytest.approx(expected)


def test_pseudovoigt():
    """Regression test"""

    x = np.linspace(-0.8, 0.8, 10)
    pvoigt = models.PseudoVoigt1D()
    pvoigt.fwhm = 0.3
    pvoigt_result = pvoigt(x)

    expected = np.asarray([0.03603509, 0.05828602, 0.11206344, 0.43013659, 2.0127256])
    expected = np.hstack((expected, expected[::-1]))

    assert pvoigt_result == pytest.approx(expected)


@pytest.mark.parametrize('cls', [models.Voigt1D, models.PseudoVoigt1D])
def test_get_center(cls):

    model = cls()
    model.pos = 12.1
    assert model.get_center() == (12.1, )


@pytest.mark.parametrize('cls', [models.Voigt1D, models.PseudoVoigt1D])
def test_set_center(cls):

    model = cls()
    model.set_center(12.1)
    assert model.pos.val == 12.1
