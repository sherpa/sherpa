#
#  Copyright (C) 2021  Smithsonian Astrophysical Observatory
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

"""
Tests for sherpa.astro.utils.xspec
"""

from io import StringIO
import os
import pathlib

import pytest

from sherpa.astro.utils import xspec
from sherpa.utils.testing import requires_xspec


@requires_xspec
@pytest.mark.parametrize("path", [False, True])
def test_parse_model_dat(path):
    """Can we parse the model.dat file (via string or pathlike)

    There's limited checking of the data since we can't
    guarantee it won't change significantly over time.
    """

    from sherpa.astro.xspec import get_xspath_manager

    path = get_xspath_manager()
    mfile = os.path.join(path, 'model.dat')

    if path:
        mfile = pathlib.PurePath(mfile)

    parsed = xspec.parse_xspec_model_description(mfile)
    assert len(parsed) > 0

    # Check we only have Add, Mul, Con, and Acn types
    #
    mtypes = set()
    for model in parsed:
        mtypes.add(model.modeltype)

    assert mtypes == set(["Add", "Mul", "Con", "Acn"])


def test_parse_additive_model():
    """Can we parse an Additive model"""

    model = StringIO("""apec           3  0.         1.e20           C_apec    add  0
kT      keV     1.    0.008   0.008   64.0      64.0      .01
Abundanc " "    1.    0.      0.      5.        5.        -0.001
Redshift " "    0.   -0.999  -0.999   10.       10.       -0.01

""")

    parsed = xspec.parse_xspec_model_description(model)
    assert len(parsed) == 1
    model = parsed[0]
    assert model.name == 'apec'
    assert model.clname == 'XSapec'
    assert model.modeltype == 'Add'
    assert model.flags == pytest.approx([0])
    assert model.elo == pytest.approx(0)
    assert model.ehi == pytest.approx(1e20)
    assert model.initString is None
    assert model.funcname == 'apec'
    assert model.language == 'C++ style'

    assert len(model.pars) == 4

    for par in model.pars:
        assert par.paramtype == 'Basic'

    assert model.pars[0].name == 'kT'
    assert not model.pars[0].frozen
    assert model.pars[0].units == 'keV'
    assert model.pars[0].default == pytest.approx(1)
    assert model.pars[0].softmin == pytest.approx(0.008)
    assert model.pars[0].hardmin == pytest.approx(0.008)
    assert model.pars[0].softmax == pytest.approx(64.0)
    assert model.pars[0].hardmax == pytest.approx(64.0)

    assert model.pars[1].name == 'Abundanc'
    assert model.pars[1].frozen
    assert model.pars[1].units is None
    assert model.pars[1].default == pytest.approx(1)
    assert model.pars[1].softmin == pytest.approx(0.0)
    assert model.pars[1].hardmin == pytest.approx(0.0)
    assert model.pars[1].softmax == pytest.approx(5.0)
    assert model.pars[1].hardmax == pytest.approx(5.0)

    assert model.pars[2].name == 'Redshift'
    assert model.pars[2].frozen
    assert model.pars[2].units is None
    assert model.pars[2].default == pytest.approx(0)
    assert model.pars[2].softmin == pytest.approx(-0.999)
    assert model.pars[2].hardmin == pytest.approx(-0.999)
    assert model.pars[2].softmax == pytest.approx(10.0)
    assert model.pars[2].hardmax == pytest.approx(10.0)

    assert model.pars[3].name == 'norm'
    assert not model.pars[3].frozen
    assert model.pars[3].units is None
    assert model.pars[3].default == pytest.approx(1)
    assert model.pars[3].softmin == pytest.approx(0)
    assert model.pars[3].hardmin == pytest.approx(0)
    assert model.pars[3].softmax == pytest.approx(1e24)
    assert model.pars[3].hardmax == pytest.approx(1e24)


def test_parse_multiplicative_model():
    """Can we parse a Multiplicative model"""

    model = StringIO("""zdust          4   0.         1.e20          mszdst    mul  0 1
$method   " "   1       1       1       3       3       -0.01
E_BmV    " "   0.1     0.0     0.0     100.    100.     0.01
*redshift " "     0.0
break  Hz   2.42E17 1.E10   1.E15     1.E19    1.E25    1E10
""")

    parsed = xspec.parse_xspec_model_description(model)
    assert len(parsed) == 1
    model = parsed[0]
    assert model.name == 'zdust'
    assert model.clname == 'XSzdust'
    assert model.modeltype == 'Mul'
    assert model.flags == pytest.approx([0, 1])
    assert model.elo == pytest.approx(0)
    assert model.ehi == pytest.approx(1e20)
    assert model.initString is None
    assert model.funcname == 'mszdst'
    assert model.language == 'Fortran - single precision'

    assert len(model.pars) == 4

    assert model.pars[0].paramtype == 'Switch'
    assert model.pars[0].name == 'method'
    assert model.pars[0].units is None
    assert model.pars[0].default == pytest.approx(1)
    assert model.pars[0].softmin == pytest.approx(1)
    assert model.pars[0].hardmin == pytest.approx(1)
    assert model.pars[0].softmax == pytest.approx(3)
    assert model.pars[0].hardmax == pytest.approx(3)

    assert model.pars[1].paramtype == 'Basic'
    assert model.pars[1].name == 'E_BmV'
    assert not model.pars[1].frozen
    assert model.pars[1].units is None
    assert model.pars[1].default == pytest.approx(0.1)
    assert model.pars[1].softmin == pytest.approx(0.0)
    assert model.pars[1].hardmin == pytest.approx(0.0)
    assert model.pars[1].softmax == pytest.approx(100.0)
    assert model.pars[1].hardmax == pytest.approx(100.0)

    assert model.pars[2].paramtype == 'Scale'
    assert model.pars[2].name == 'redshift'
    assert model.pars[2].units is None
    assert model.pars[2].default == pytest.approx(0)
    assert model.pars[2].softmin is None
    assert model.pars[2].hardmin is None
    assert model.pars[2].softmax is None
    assert model.pars[2].hardmax is None

    assert model.pars[3].paramtype == 'Basic'
    assert model.pars[3].name == 'break_'
    assert not model.pars[3].frozen
    assert model.pars[3].units == 'Hz'
    assert model.pars[3].default == pytest.approx(2.42e17)
    assert model.pars[3].softmin == pytest.approx(1e15)
    assert model.pars[3].hardmin == pytest.approx(1e10)
    assert model.pars[3].softmax == pytest.approx(1e19)
    assert model.pars[3].hardmax == pytest.approx(1e25)



def test_create_model_additive():
    """Fake up an additive model"""

    model = StringIO("""apec           1  0.         1.e20           C_apec    add  0
kT      keV     1.    0.008   0.008   64.0      64.0      .01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    # Very basic checks of the output
    #
    python = converted.python.split('\n')
    compiled = converted.compiled.split('\n')

    assert 'class XSapec(XSAdditiveModel):' in python
    assert '    _calc = _models.C_apec' in python
    assert "    def __init__(self, name='apec'):" in python
    assert "        self.kT = Parameter(name, 'kT', 1.0, min=0.008, max=64.0, hard_min=0.008, hard_max=64.0, units='keV')" in python
    assert "        self.norm = Parameter(name, 'norm', 1.0, min=0.0, max=1e+24, hard_min=0.0, hard_max=1e+24)" in python
    assert '        XSAdditiveModel.__init__(self, name, (self.kT,self.norm))' in python

    assert converted.compiled.find('\nextern "C" {\n\n}\n') > -1
    assert '  XSPECMODELFCT_C_NORM(C_apec, 2),' in compiled


def test_create_model_multiplicative():
    """Fake up a multplicative model"""

    model = StringIO("""abcd           1  0.         1.e20           foos    mul  0
nH      cm^-3   1.0   1.e-6  1.e-5  1.e19  1.e20   -0.01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    # Very basic checks of the output
    #
    python = converted.python.split('\n')
    compiled = converted.compiled.split('\n')

    assert 'class XSabcd(XSMultiplicativeModel):' in python
    assert '    _calc = _models.foos' in python
    assert "    def __init__(self, name='abcd'):" in python
    assert "        self.nH = Parameter(name, 'nH', 1.0, min=1e-05, max=1e+19, hard_min=1e-06, hard_max=1e+20, frozen=True, units='cm^-3')" in python
    assert '        XSMultiplicativeModel.__init__(self, name, (self.nH,))' in python

    assert '  void foos_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);' in compiled
    assert '  XSPECMODELFCT(foos, 1),' in compiled


def test_create_model_convolution():
    """Fake up a convolution model"""

    model = StringIO("""rgsxsrc        1  0.         1.e20           rgsxsrc   con  0
order    " "  -1.   -3.    -3.      -1.       -1.       -1

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    # Very basic checks of the output
    #
    python = converted.python.split('\n')
    compiled = converted.compiled.split('\n')

    assert 'class XSrgsxsrc(XSConvolutionKernel):' in python
    assert '    _calc = _models.rgsxsrc' in python
    assert "    def __init__(self, name='rgsxsrc'):" in python
    assert "        self.order = Parameter(name, 'order', -1.0, min=-3.0, max=-1.0, hard_min=-3.0, hard_max=-1.0, frozen=True)" in python
    assert '        XSConvolutionKernel.__init__(self, name, (self.order,))' in python

    assert '  void rgsxsrc_(float* ear, int* ne, float* param, int* ifl, float* photar, float* photer);' in compiled
    assert '  XSPECMODELFCT_CON_F77(rgsxsrc, 1),' in compiled
