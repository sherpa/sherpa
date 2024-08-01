#
#  Copyright (C) 2021, 2023, 2024
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

"""
Tests for sherpa.astro.utils.xspec
"""

from io import StringIO
import logging
import os
import pathlib

import pytest

from sherpa.astro.utils import xspec
from sherpa.utils.testing import requires_xspec


@pytest.mark.parametrize("namefunc", [None, True, lambda x: 3])
def test_fail_invalid_namefunc(namefunc):
    """We do basic checking on namefunc"""

    # The contents shouldn't be checked
    model = '/dev/null'

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model, namefunc)

    assert str(ve.value) == 'namefunc must be a callable which takes and returns a string'


def test_fail_invalid_model_line():
    """There is some check of a model line"""

    model = StringIO("""apec           3  0.         1.e20           C_apec    add
Abundanc " "    1.    0.      0.      5.        5.        -0.001

""")

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model)

    assert str(ve.value).startswith('Expected: modelname npars ')


def test_fail_missing_parameters():
    """We have less parameters than required"""

    model = StringIO("""apec           3  0.         1.e20           C_apec    add  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001

""")

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model)

    assert str(ve.value) == 'model=apec missing 2 parameters'


def test_fail_toomany_parameters():
    """We have more parameters than required

    This actually errors out as it tries to parse the extra
    parameter as a model (that is, this is only checked for
    indirectly).

    """

    model = StringIO("""apec           1  0.         1.e20           C_apec    add  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001
Redshift " "    0.   -0.999  -0.999   10.       10.       -0.01

""")

    # As we don't have an explicit check the actual error depends on the
    # next line so it's not worth checking the message.
    #
    with pytest.raises(ValueError):
        xspec.parse_xspec_model_description(model)


def test_fail_unknown_model():
    """There is some check of model types

    It mmight be better to just skip these.
    """

    model = StringIO("""apec           1  0.         1.e20           C_apec    sub  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001

""")

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model)

    assert str(ve.value).startswith('Unexpected model type sub in:\n')


def test_fail_invalid_basic_parameter():
    """There is a check of model types.

    It might be better to just skip these.
    """

    model = StringIO("""apec           1  0.         1.e20           C_apec    sub  0
Abundanc " "    1.

""")

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model)

    assert str(ve.value).startswith('Expected 6 values after units; model=apec\n')


def test_fail_periodic_parameter():
    """We may add support for these in the future.

    """

    model = StringIO("""apec           1  0.         1.e20           C_apec    add  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001 P
""")

    with pytest.raises(ValueError) as ve:
        xspec.parse_xspec_model_description(model)

    assert str(ve.value).startswith('Periodic parameters are unsupported; model=apec:\n')


def test_warn_parse_repeated_parname(caplog):
    """Check we warn if the parameter name is repeated"""

    model = StringIO("""apec           3  0.         1.e20           C_apec    add  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001
kT      keV     1.    0.008   0.008   64.0      64.0      .01
abundanc " "    0.   -0.999  -0.999   10.       10.       -0.01

""")

    assert len(caplog.records) == 0
    xspec.parse_xspec_model_description(model)

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.astro.utils.xspec'
    assert lvl == logging.WARNING
    assert msg == 'model=apec re-uses parameter name abundanc'


def test_skip_repeated_model(caplog):
    """Technically we could use the same model but with different
    parameters but there was a bug in a version of XSPEC where
    the same model was used incorrectly, so check we handle this.

    We only skip it when we come to creating the code.

    """

    model = StringIO("""apec           1  0.         1.e20           C_foo     mul  0
Abundanc " "    1.    0.      0.      5.        5.        -0.001

fred           0  0.         1.e20           C_foo     mul  0
""")

    parsed = xspec.parse_xspec_model_description(model)
    assert len(parsed) == 2
    assert len(caplog.records) == 0

    # As an extra benefit as we skip both models we can check the
    # 'no valid models' error.
    #
    with pytest.raises(ValueError) as ve:
        xspec.create_xspec_code(parsed)

    assert str(ve.value) == 'No supported models were found!'

    assert len(caplog.records) == 2
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == 'sherpa.astro.utils.xspec'
    assert lvl == logging.WARNING
    assert msg == 'Skipping model apec as it calls foo which is used by 2 different models'

    lname, lvl, msg = caplog.record_tuples[1]
    assert lname == 'sherpa.astro.utils.xspec'
    assert lvl == logging.WARNING
    assert msg == 'Skipping model fred as it calls foo which is used by 2 different models'


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


def check_compiled(got, suffix):
    """Check the compiled output.

    The output depends on the XSPEC version because it will cause

        #define XSPEC_major_minor_micro

    lines to be included. To support this the code splits on these
    markers and checks the text before and after this, and checks the
    define lines make sense.

    This test must be run with the requires_xspec fixture. This
    means that there must be one XSPEC define line present.

    Note that the text before the defines does not depend on
    the models being evaluated.

    """

    # This should be expressible with the re module, but it would
    # be a bit hard to read.
    #
    idx = got.find("\n#define XSPEC_")
    if idx == -1:
        assert False, f"No #define XSPEC fragment found."

    prefix = '''// Includes

#include <iostream>

#include <xsTypes.h>
#include <XSFunctions/Utilities/funcType.h>
'''

    assert got[:idx] == prefix

    # Check the XSPEC version lines
    #   - have form major_minor_micro
    #   - are in order
    #   - there's at least one of them.
    #
    toks = got[idx + 1:].split("\n")
    prev = (0, 0, 0)
    while True:
        tok = toks.pop(0)
        if not tok.startswith("#define XSPEC_"):
            break

        ntok = tok[14:].split("_")
        assert len(ntok) == 3

        # check we can convert to integers
        version = tuple([int(n) for n in ntok])
        assert version > prev
        prev = version

    # There should be a few fixed lines until we get to the Defines
    # header (which could have been included here but the tests
    # read better with the text there rather than here).
    #
    assert tok == ""

    tok = toks.pop(0)
    assert tok == '#include "sherpa/astro/xspec_extension.hh"'

    tok = toks.pop(0)
    assert tok == ""

    # Reconstruct the text to make it easy to compare to the user input.
    #
    got = "\n".join(toks)
    assert got == suffix


@requires_xspec
def test_create_model_additive(caplog):
    """Fake up an additive model"""

    model = StringIO("""apec           1  0.         1.e20           C_apec    add  0
kT      keV     1.    0.008   0.008   64.0      64.0      .01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    expected = '''
class XSapec(XSAdditiveModel):
    """XSPEC AdditiveModel: apec

    Parameters
    ----------
    kT
    norm

    """
    _calc = _models.C_apec

    def __init__(self, name='apec'):
        self.kT = XSParameter(name, 'kT', 1.0, min=0.008, max=64.0, hard_min=0.008, hard_max=64.0, units='keV')

        # norm parameter is automatically added by XSAdditiveModel
        pars = (self.kT,)
        XSAdditiveModel.__init__(self, name, pars)

'''
    assert converted.python == expected

    check_compiled(converted.compiled,
                   '''// Defines

void cppModelWrapper(const double* energy, int nFlux, const double* params,
  int spectrumNumber, double* flux, double* fluxError, const char* initStr,
  int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
  int, RealArray&, RealArray&, const string&));

extern "C" {
  XSCCall apec;
  void C_apec(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr) {
    const size_t nPar = 2;
    cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError, initStr, nPar, apec);
  }
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT_C_NORM(C_apec, 2),
  { NULL, NULL, 0, NULL }
};

// Module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
''')

    assert len(caplog.records) == 0


@requires_xspec
def test_create_model_multiplicative(caplog):
    """Fake up a multplicative model"""

    model = StringIO("""abcd           1  0.         1.e20           foos    mul  0
nH      cm^-3   1.0   1.e-6  1.e-5  1.e19  1.e20   -0.01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    expected = '''
class XSabcd(XSMultiplicativeModel):
    """XSPEC MultiplicativeModel: abcd

    Parameters
    ----------
    nH

    """
    _calc = _models.foos

    def __init__(self, name='abcd'):
        self.nH = XSParameter(name, 'nH', 1.0, min=1e-05, max=1e+19, hard_min=1e-06, hard_max=1e+20, frozen=True, units='cm^-3')

        pars = (self.nH,)
        XSMultiplicativeModel.__init__(self, name, pars)

'''
    assert converted.python == expected

    check_compiled(converted.compiled,
                   '''// Defines

extern "C" {
  xsf77Call foos_;
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT(foos, 1),
  { NULL, NULL, 0, NULL }
};

// Module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
''')

    assert len(caplog.records) == 0


@requires_xspec
def test_create_model_convolution(caplog):
    """Fake up a convolution model

    Also checks case handling
    """

    model = StringIO("""rgsxsrc        1  0.         1.e20           rgsxSRC   con  0
order    " "  -1.   -3.    -3.      -1.       -1.       -1

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    expected = '''
class XSrgsxsrc(XSConvolutionKernel):
    """XSPEC ConvolutionKernel: rgsxsrc

    Parameters
    ----------
    order

    """
    _calc = _models.rgsxsrc

    def __init__(self, name='rgsxsrc'):
        self.order = XSParameter(name, 'order', -1.0, min=-3.0, max=-1.0, hard_min=-3.0, hard_max=-1.0, frozen=True)

        pars = (self.order,)
        XSConvolutionKernel.__init__(self, name, pars)

'''
    assert converted.python == expected

    check_compiled(converted.compiled,
                   '''// Defines

extern "C" {
  xsf77Call rgsxsrc_;
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT_CON_F77(rgsxsrc, 1),
  { NULL, NULL, 0, NULL }
};

// Module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
''')

    assert len(caplog.records) == 0


@requires_xspec
def test_create_model_spectrum_specific(caplog):
    """Fake up an additive model that is marked as spectrum specific"""

    model = StringIO("""abcd           1  0.         1.e20           c_foos    add  0 1
xs   ""    10    1 2  20 30  0.01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    expected = '''import warnings

class XSabcd(XSAdditiveModel):
    """XSPEC AdditiveModel: abcd

    Parameters
    ----------
    xs
    norm

    """
    _calc = _models.foos

    def __init__(self, name='abcd'):
        self.xs = XSParameter(name, 'xs', 10.0, min=2.0, max=20.0, hard_min=1.0, hard_max=30.0)

        # norm parameter is automatically added by XSAdditiveModel
        pars = (self.xs,)
        XSAdditiveModel.__init__(self, name, pars)
        self._use_caching = False
        warnings.warn('support for models like xsabcd (recalculated per spectrum) is untested.')

'''

    assert converted.python == expected

    check_compiled(converted.compiled,
                   '''// Defines

extern "C" {
  xsccCall foos;
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT_C_NORM(foos, 2),
  { NULL, NULL, 0, NULL }
};

// Module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
''')

    # Do we get a warning to screen?
    #
    assert len(caplog.records) == 1
    assert caplog.records[0].msg == "abcd needs to be re-calculated per spectrum; this is untested."


@requires_xspec
def test_create_model_error_specific(caplog):
    """Fake up an additive model that is marked as needing the errors"""

    model = StringIO("""abcd           1  0.         1.e20           C_foos    add  1 0
xs   ""    10    1 2  20 30  0.01

""")
    parsed = xspec.parse_xspec_model_description(model)

    converted = xspec.create_xspec_code(parsed)

    expected = '''import warnings

class XSabcd(XSAdditiveModel):
    """XSPEC AdditiveModel: abcd

    Parameters
    ----------
    xs
    norm

    """
    _calc = _models.C_foos

    def __init__(self, name='abcd'):
        self.xs = XSParameter(name, 'xs', 10.0, min=2.0, max=20.0, hard_min=1.0, hard_max=30.0)

        # norm parameter is automatically added by XSAdditiveModel
        pars = (self.xs,)
        XSAdditiveModel.__init__(self, name, pars)
        warnings.warn('support for models like xsabcd (variances are calculated by the model) is untested.')

'''
    assert converted.python == expected

    check_compiled(converted.compiled,
                   '''// Defines

void cppModelWrapper(const double* energy, int nFlux, const double* params,
  int spectrumNumber, double* flux, double* fluxError, const char* initStr,
  int nPar, void (*cppFunc)(const RealArray&, const RealArray&,
  int, RealArray&, RealArray&, const string&));

extern "C" {
  XSCCall foos;
  void C_foos(const double* energy, int nFlux, const double* params, int spectrumNumber, double* flux, double* fluxError, const char* initStr) {
    const size_t nPar = 2;
    cppModelWrapper(energy, nFlux, params, spectrumNumber, flux, fluxError, initStr, nPar, foos);
  }
}

// Wrapper

static PyMethodDef Wrappers[] = {
  XSPECMODELFCT_C_NORM(C_foos, 2),
  { NULL, NULL, 0, NULL }
};

// Module

static struct PyModuleDef wrapper_module = {
  PyModuleDef_HEAD_INIT,
  "_models",
  NULL,
  -1,
  Wrappers,
};

PyMODINIT_FUNC PyInit__models(void) {
  import_array();
  return PyModule_Create(&wrapper_module);
}
''')

    # Do we get a warning to screen?
    #
    assert len(caplog.records) == 1
    assert caplog.records[0].msg == "abcd calculates model variances; this is untested/unsupported in Sherpa"


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


@requires_xspec
def test_convert_model_dat():
    """Can we convert the model.dat file?

    This is more just an existence test (i.e. does it work) with
    limited checking of the output. Checking the warnings that
    are created by parsing/converting the XSPEC model.dat is
    tricky to do reliably as it depends on the XSPEC version,
    so we do not check them (there are specific tests above).

    """

    from sherpa.astro.xspec import get_xspath_manager

    path = get_xspath_manager()
    mfile = os.path.join(path, 'model.dat')

    # We do not check the return value, just that it runs
    parsed = xspec.parse_xspec_model_description(mfile)
    xspec.create_xspec_code(parsed)
