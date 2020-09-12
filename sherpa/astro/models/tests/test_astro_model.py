#
#  Copyright (C) 2007, 2016, 2018, 2020  Smithsonian Astrophysical Observatory
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
from sherpa.utils import SherpaFloat
from sherpa.utils.testing import SherpaTestCase
from sherpa.models.model import ArithmeticModel, RegriddableModel2D, RegriddableModel1D, boolean_to_byte
from sherpa.models.basic import Const


EXCLUDED_MODEL_CLASSES = (ArithmeticModel, RegriddableModel1D, RegriddableModel2D, Const)


def test_create_and_evaluate():

    x = np.arange(1.0, 5.0)
    count = 0

    for cls in dir(models):
        clsobj = getattr(models, cls)

        if not isinstance(clsobj, type) or \
           not issubclass(clsobj, ArithmeticModel) or \
           clsobj in EXCLUDED_MODEL_CLASSES:
            continue

        # These have a very different interface than the others
        if cls in ('JDPileup', 'MultiResponseSumModel'):
            continue

        m = clsobj()
        assert type(m).__name__.lower() == m.name
        count += 1

        if m.name == 'linebroad':
            m.vsini = 1e6

        try:
            if m.name.count('2d') or m.name == 'hubblereynolds':
                pt_out  = m(x, x)
                int_out = m(x, x, x, x)
            else:
                pt_out  = m(x)
                int_out = m(x, x)
        except ValueError:
            assert False, "evaluation of model '{}' failed".format(cls)

        for out in (pt_out, int_out):
            assert out.dtype.type is SherpaFloat
            assert out.shape == (4, )

    assert count == 19


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

    expected = np.asarray([0.07779156, 0.13213921, 0.27123415, 0.67357158, 1.35326402])
    expected = np.hstack((expected, expected[::-1]))

    assert voigt_result == pytest.approx(expected)
