#
#  Copyright (C) 2017, 2018, 2020, 2021, 2022
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
Should these tests be moved to test_session.py?

Note that this test is almost duplicated in
sherpa/astro/ui/tests/test_astro_ui_unit.py
"""

import logging

import numpy as np

import pytest

from sherpa import ui
from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel
from sherpa.utils.err import ArgumentTypeErr, IdentifierErr
from sherpa.utils.logging import SherpaVerbosity


# This is part of #397
#
def test_list_samplers():
    """Ensure list_samplers returns a list."""

    samplers = ui.list_samplers()

    assert isinstance(samplers, list)
    assert len(samplers) > 0


def test_list_samplers_contents():
    """Are the expected values included"""

    # Test that the expected values exist in this list,
    # but do not enforce these are the only values.
    #
    samplers = ui.list_samplers()
    for expected in ['mh', 'metropolismh']:
        assert expected in samplers


def test_all_has_no_repeated_elements():
    """Ensure __all__ does not contain repeated elements.

    It is not actually a problem if this happens, but it does
    indicate the possibility of confusion over what functionality
    the repeated symbol has (originally noticed with erf, which
    could be either sherpa.utils.erf or the ModelWrapper version
    of sherpa.models.basic.Erf). See
    https://github.com/sherpa/sherpa/issues/502
    """

    n1 = len(ui.__all__)
    n2 = len(set(ui.__all__))
    assert n1 == n2


@pytest.mark.parametrize("func", [ui.notice_id, ui.ignore_id])
@pytest.mark.parametrize("lo", [None, 1, "1:5"])
def test_check_ids_not_none(func, lo):
    """Check they error out when id is None

    There used to be the potential for different behavior depending
    if the lo argument was a string or not,hence the check.
    """

    with pytest.raises(ArgumentTypeErr) as exc:
        func(None, lo)

    assert str(exc.value) == "'ids' must be an identifier or list of identifiers"


@pytest.mark.parametrize("func", [ui.notice, ui.ignore])
@pytest.mark.parametrize("lo,hi", [(1, 5), (1, None), (None, 5), (None, None),
                                   ("1:5", None)])
def test_filter_no_data_is_an_error(func, lo, hi, clean_ui):
    """Does applying a filter lead to an error?

    This test was added because it was noted that an update for Python 3
    had lead to an error state not being reached.
    """

    with pytest.raises(IdentifierErr) as ie:
        func(lo, hi)

    assert str(ie.value) == 'No data sets found'


def test_save_filter_data1d(tmp_path, clean_ui):
    """Check save_filter [Data1D]"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x)

    ui.notice(2, 4)
    ui.notice(6, 8)

    outfile = tmp_path / "filter.dat"
    ui.save_filter(str(outfile))

    expected = [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]

    d = ui.unpack_data(str(outfile), colkeys=['X', 'FILTER'])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None


def test_set_iter_method_type_not_string():
    with pytest.raises(ArgumentTypeErr) as te:
        ui.set_iter_method(23)

    assert str(te.value) == "'23' must be a string"


def test_set_iter_method_type_not_enumeration():
    with pytest.raises(TypeError) as te:
        ui.set_iter_method('a random string')

    assert str(te.value) == "a random string is not an iterative fitting method"


class NonIterableObject:
    """Something that tuple(..) of will error out on"""

    pass


@pytest.mark.parametrize("func",
                         [ui.notice_id, ui.ignore_id])
def test_filter_errors_out_invalid_id(func):
    """Just check we create the expected error message.

    Somewhat contrived.
    """

    ids = NonIterableObject()
    with pytest.raises(ArgumentTypeErr) as te:
        func(ids)

    assert str(te.value) == "'ids' must be an identifier or list of identifiers"


def test_set_model_autoassign_func_type():

    with pytest.raises(ArgumentTypeErr) as te:
        ui.set_model_autoassign_func(23)

    assert str(te.value) == "'func' must be a function or other callable object"


class DummyModel(ArithmeticModel):
    pass


def test_guess_warns_no_guess_names_model(caplog, clean_ui):
    """Do we warn when the named model has no guess"""

    ui.load_arrays(1, [1, 2, 3], [-3, 4, 5])
    cpt = DummyModel('dummy')

    assert len(caplog.records) == 0
    ui.guess(cpt)

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "WARNING: No guess found for dummy"


def test_guess_warns_no_guess_no_argument(caplog, clean_ui):
    """Do we warn when the (implied) model has no guess"""

    ui.load_arrays(1, [1, 2, 3], [-3, 4, 5])
    cpt = DummyModel('dummy')
    ui.set_source(cpt + cpt)

    assert len(caplog.records) == 0
    ui.guess()

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "WARNING: No guess found for (dummy + dummy)"


class Parameter2(Parameter):
    """All we want is a sub-class of Parameter."""

    # What version of Python allows us to drop the pass statement?
    pass


class Const(ArithmeticModel):
    """A constant model"""

    def calc(self, pars, *args, **kwargs):
        return pars[0] * np.ones_like(args[0])


class Const1(Const):
    """A constant model using Parameter

    sherpa.models.basic.Const1D could have been used but here we can
    see that Const1/Const2 are the same, apart from the parameter
    class.

    """

    def __init__(self, name='const1'):
        self.con = Parameter(name, 'con', 1)
        Const.__init__(self, name, (self.con, ))


class Const2(Const):
    """A constant model using Parameter2"""

    def __init__(self, name='const2'):
        self.con = Parameter2(name, 'con', 1)
        Const.__init__(self, name, (self.con, ))


@pytest.mark.parametrize("mdlcls", [Const1, pytest.param(Const2, marks=pytest.mark.xfail)])
@pytest.mark.parametrize("method,getter",
                         [(ui.covar, ui.get_covar_results),
                          (ui.conf, ui.get_conf_results),
                          (ui.proj,ui.get_proj_results)])
def test_est_errors_works_single_parameter(mdlcls, method, getter, clean_ui):
    """This is issue #1397.

    Rather than require XSPEC, we create a subclass of the Parameter
    class to check it works. We are not too concerned with the actual
    results hence the relatively low tolerance on the numeric checks.

    """

    mdl = mdlcls()

    ui.load_arrays(1, [1, 2, 3, 4], [4, 2, 1, 3.5])
    ui.set_source(mdl)
    with SherpaVerbosity("ERROR"):
        ui.fit()

        # this is where #1397 fails with Const2
        method(mdl.con)

    atol = 1e-4
    assert ui.calc_stat() == pytest.approx(0.7651548418626658, abs=atol)

    results = getter()
    assert results.parnames == (f"{mdl.name}.con", )
    assert results.sigma == pytest.approx(1.0)

    assert results.parvals == pytest.approx((2.324060647544594, ), abs=atol)

    # The covar errors are -/+ 1.3704388763054511
    #     conf             -1.3704388763054511 / +1.3704388763054514
    #     proj             -1.3704388762971822 / +1.3704388763135826
    #
    err = 1.3704388763054511
    assert results.parmins == pytest.approx((-err, ), abs=atol)
    assert results.parmaxes == pytest.approx((err, ), abs=atol)
