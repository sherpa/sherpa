#
#  Copyright (C) 2007, 2016, 2018, 2020
#      Smithsonian Astrophysical Observatory
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
from pytest import approx

import sherpa.models.basic as basic
from sherpa.utils import SherpaFloat, _utils
from sherpa.models.model import ArithmeticModel, RegriddableModel1D, RegriddableModel2D


def userfunc(pars, x, *args, **kwargs):
    return x


class TestBasic:
    excluded_models = (ArithmeticModel, RegriddableModel1D, RegriddableModel2D, basic.Const)

    def test_create_and_evaluate(self):
        x = np.arange(1.0, 5.0)
        count = 0

        for cls in dir(basic):
            clsobj = getattr(basic, cls)

            if not isinstance(clsobj, type) \
                or not issubclass(clsobj, ArithmeticModel) \
                or clsobj in self.excluded_models:
                continue

            # These have very different interfaces than the others
            if cls == 'Integrator1D' or cls == 'Integrate1D':
                continue

            m = clsobj()
            if isinstance(m, basic.TableModel):
                m.load(x,x)
            if isinstance(m, basic.UserModel):
                m.calc = userfunc
            assert type(m).__name__.lower() == m.name
            count += 1

            try:
                if m.name.count('2d'):
                    pt_out  = m(x, x)
                    int_out = m(x, x, x, x)
                else:
                    if m.name in ('log', 'log10'):
                        xx = -x
                    else:
                        xx = x
                    pt_out  = m(xx)
                    int_out = m(xx, xx)
            except ValueError:
                self.fail("evaluation of model '%s' failed" % cls)

            for out in (pt_out, int_out):
                assert out.dtype.type is SherpaFloat
                assert out.shape == x.shape

        assert count == 33


def test_basic():
    my_test_basic = TestBasic()
    my_test_basic.test_create_and_evaluate()


class TestVoigt1D:

    def test_voigt(self):
        def myVoigt(x, alpha, gamma):
            sigma = alpha / np.sqrt(2 * np.log(2))
            arg = (x + 1j *gamma) / sigma / np.sqrt(2)
            wofz = _utils.calc_wofz
            return np.real(wofz(arg)) / sigma / np.sqrt( 2 * np.pi )
    
        x = np.linspace(-0.8, 0.8, 10)
        voigt = basic.Voigt1D()
        voigt_result = voigt(x)

        alpha, gamma = 0.1, 0.1
        myvoigt_result = myVoigt(x, alpha, gamma)

        scipy_result = [0.05065948, 0.08477187, 0.17116719, 0.50527797,
                        1.79706175, 1.79706175, 0.50527797, 0.17116719,
                        0.08477187, 0.05065948]

        assert voigt_result == approx(myvoigt_result)
        assert voigt_result == approx(scipy_result)

def test_voigt():
    my_test_voigt = TestVoigt1D()
    my_test_voigt.test_voigt()
