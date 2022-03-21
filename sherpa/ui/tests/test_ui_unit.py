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

from io import StringIO
import logging

import numpy as np

import pytest

from sherpa import ui
from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel
from sherpa.utils.err import ArgumentTypeErr, IdentifierErr, ParameterErr
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


@pytest.mark.parametrize("mdlcls", [Const1, Const2])
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


@pytest.mark.parametrize("method", [ui.conf, ui.covar, ui.proj])
def test_err_estimate_errors_on_frozen(method, clean_ui):
    """Check we error out with frozen par with conf/proj/covar.

    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.set_source(ui.polynom1d.mdl)
    with pytest.raises(ParameterErr) as exc:
        method(mdl.c0, mdl.c1)

    assert str(exc.value) == "parameter 'mdl.c1' is frozen"


@pytest.mark.parametrize("method", [ui.conf, ui.covar, ui.proj])
def test_err_estimate_errors_model_all_frozen(method, clean_ui):
    """Check we error out with frozen model with conf/proj/covar.

    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.set_source(ui.polynom1d.mdl)
    for par in mdl.pars:
        par.freeze()

    with pytest.raises(ParameterErr) as exc:
        method(mdl)

    assert str(exc.value) == "Model 'polynom1d.mdl' has no thawed parameters"


@pytest.mark.parametrize("method", [ui.conf, ui.covar, ui.proj])
@pytest.mark.parametrize("id", [1, "xx"])
@pytest.mark.parametrize("otherids", [[2, 3], ["foo", "bar"]])
def test_err_estimate_errors_on_list_argument(method, id, otherids, clean_ui):
    """Check we error out with a list argument with conf/proj/covar.

    We had documented that you could say conf(1, [2, 3]) but this
    is not true, so check it does error out. Fortunately we can
    do this without setting up any dataset or model.

    """

    with pytest.raises(ArgumentTypeErr) as exc:
        method(id, otherids)

    assert str(exc.value) == "identifiers must be integers or strings"


def setup_err_estimate_multi_ids(strings=False):
    """Create the environment used in test_err_estimate_xxx tests.

    The model being fit is polynom1d with c0=50 c1=-2
    and was evaluated and passed through sherpa.utils.poisson_noise
    to create the datasets.

    Since we can have string or integer ids we allow either,
    but do not try to mix them.

    """

    if strings:
        id1 = "1"
        id2 = "2"
        id3 = "3"
    else:
        id1 = 1
        id2 = 2
        id3 = 3

    ui.load_arrays(id1, [1, 3, 7, 12], [50, 40,27, 20])
    ui.load_arrays(id2, [-3, 4, 5], [55, 34, 37])
    ui.load_arrays(id3, [10, 12, 20], [24, 26, 7])

    # NOTE: dataset "not-used" is not used in the fit and is not
    # drawn from the distributino used to create the other datasets.
    #
    ui.load_arrays("not-used", [2000, 2010, 2020], [10, 12, 14])

    mdl = ui.create_model_component("polynom1d", "mdl")
    mdl.c1.thaw()
    ui.set_source(id1, mdl)
    ui.set_source(id2, mdl)
    ui.set_source(id3, mdl)

    # apply the model to dataset not-used just so we can check we
    # don't end up using it
    mdl_not_used = ui.create_model_component("scale1d", "mdl_not_used")
    ui.set_source("not-used", mdl + mdl_not_used)

    # use cstat so we have an approximate goodness-of-fit just to
    # check we are getting sensible results.
    #
    ui.set_stat("cstat")
    ui.set_method("simplex")


# This is a regression test, there's no "correct" value to test against.
#
ERR_EST_C0_MIN = -2.8006781174692676
ERR_EST_C0_MAX = 2.9038373212762494
ERR_EST_C1_MIN = -0.19798615552319254
ERR_EST_C1_MAX = 0.21022843992924245


@pytest.mark.parametrize("strings", [False, True])
@pytest.mark.parametrize("idval,otherids",
                         [(1, (2, 3)),
                          (2, [3, 1]),
                          (3, [2, 1])])
def test_err_estimate_multi_ids(strings, idval, otherids, clean_ui):
    """Ensure we can use multiple ids with conf/proj/covar.

    Since this uses the same logic we only test the conf routine;
    ideally we'd use all but that's harder to test.

    The fit and error analysis should be the same however the ordering
    is done.
    """

    # This is a bit ugly
    if strings:
        idval = str(idval)
        if type(otherids) == tuple:
            otherids = (str(otherids[0]), str(otherids[1]))
        else:
            otherids = [str(otherids[0]), str(otherids[1])]

    datasets = tuple([idval] + list(otherids))

    setup_err_estimate_multi_ids(strings=strings)
    ui.fit(idval, *otherids)

    # The "reduced statistic" is ~0.42 for the fit.
    #
    res = ui.get_fit_results()
    assert res.datasets == datasets
    assert res.numpoints == 10  # sum of datasets 1, 2, 3
    assert res.statval == pytest.approx(3.379367979541458)

    # since there's a model assigned to dataset not-used the
    # overall statistic is not the same as res.statval.
    #
    assert ui.calc_stat() == pytest.approx(4255.615602052843)

    assert mdl.c0.val == pytest.approx(46.046607302070015)
    assert mdl.c1.val == pytest.approx(-1.9783953989993386)

    ui.conf(*datasets)
    res = ui.get_conf_results()

    assert res.datasets == datasets
    assert res.parnames == ("mdl.c0", "mdl.c1")

    assert res.parmins == pytest.approx([ERR_EST_C0_MIN, ERR_EST_C1_MIN])
    assert res.parmaxes == pytest.approx([ERR_EST_C0_MAX, ERR_EST_C1_MAX])


@pytest.mark.parametrize("strings", [False, True])
@pytest.mark.parametrize("idval,otherids",
                         [(1, (2, 3)),
                          (2, [3, 1]),
                          (3, [2, 1])])
def test_err_estimate_model(strings, idval, otherids, clean_ui):
    """Ensure we can use model with conf/proj/covar.

    This is test_err_estimate_multi_ids but

      - added an extra model to each source (that evaluates to 0)
      - we include the model expression in the call.

    The fit and error analysis should be the same however the ordering
    is done.
    """

    # This is a bit ugly
    if strings:
        idval = str(idval)
        if type(otherids) == tuple:
            otherids = (str(otherids[0]), str(otherids[1]))
        else:
            otherids = [str(otherids[0]), str(otherids[1])]

    datasets = tuple([idval] + list(otherids))

    setup_err_estimate_multi_ids(strings=strings)

    zero = ui.create_model_component("scale1d", "zero")
    zero.c0 = 0
    zero.c0.freeze()

    for id in datasets:
        # In this case we have
        #   orig == mdl
        # but let's be explicit in case the code changes
        #
        orig = ui.get_source(id)
        ui.set_source(id, orig + zero)

    ui.fit(idval, *otherids)

    res = ui.get_fit_results()
    assert res.datasets == datasets
    assert res.numpoints == 10
    assert res.statval == pytest.approx(3.379367979541458)
    assert ui.calc_stat() == pytest.approx(4255.615602052843)
    assert mdl.c0.val == pytest.approx(46.046607302070015)
    assert mdl.c1.val == pytest.approx(-1.9783953989993386)

    # I wanted to have zero.co thawed at this stage, but then we can not
    # use the ERR_EST_C0/1_xxx values as the fit has changed (and mdl.c0
    # and zero.c0 are degenerate to boot).
    #
    ui.conf(*datasets, mdl)
    res = ui.get_conf_results()

    assert res.datasets == datasets
    assert res.parnames == ("mdl.c0", "mdl.c1")

    assert res.parmins == pytest.approx([ERR_EST_C0_MIN, ERR_EST_C1_MIN])
    assert res.parmaxes == pytest.approx([ERR_EST_C0_MAX, ERR_EST_C1_MAX])


@pytest.mark.parametrize("strings", [False, True])
@pytest.mark.parametrize("idval,otherids",
                         [(1, (2, 3)),
                          (2, [3, 1]),
                          (3, [2, 1])])
def test_err_estimate_single_parameter(strings, idval, otherids, clean_ui):
    """Ensure we can fti a single parameter with conf/proj/covar.

    Since this uses the same logic we only test the conf routine;
    ideally we'd use all but that's harder to test.

    We use the same model as test_err_estimate_multi_ids but
    here we only want to evaluate the error for the mdl.c1 component.

    The fit and error analysis should be the same however the ordering
    is done.
    """

    # This is a bit ugly
    if strings:
        idval = str(idval)
        if type(otherids) == tuple:
            otherids = (str(otherids[0]), str(otherids[1]))
        else:
            otherids = [str(otherids[0]), str(otherids[1])]

    datasets = tuple([idval] + list(otherids))
    setup_err_estimate_multi_ids(strings=strings)
    ui.fit(idval, *otherids)

    # pick an odd ordering just to check we pick it up
    ui.conf(datasets[0], mdl.c1, datasets[1], datasets[2])
    res = ui.get_conf_results()

    assert res.datasets == datasets
    assert res.parnames == ("mdl.c1", )

    assert res.parmins == pytest.approx([ERR_EST_C1_MIN])
    assert res.parmaxes == pytest.approx([ERR_EST_C1_MAX])


def test_show_all_empty(clean_ui):
    """Checks several routines at once!"""

    out = StringIO()
    ui.show_all(outfile=out)
    assert out.getvalue() == "\n"


def test_show_all_basic(clean_ui):
    """Set up a very basic data/model/fit"""

    ui.load_arrays(1, [1, 2, 4], [3, 5, 5])
    ui.set_source(ui.scale1d.mdl)
    ui.fit()
    ui.conf()
    ui.proj()
    ui.covar()

    def get(value):
        out = StringIO()
        getattr(ui, f"show_{value}")(outfile=out)
        ans = out.getvalue()
        assert len(ans) > 1

        # trim the trailing "\n"
        return ans[:-1]

    # All we are really checking is that the show_all output is the
    # comppsite of the following. We are not checking that the
    # actual output makes sense for any command.
    #
    expected = get("data") + get("model") + get("fit") + get("conf") + \
        get("proj") + get("covar")

    got = get("all")

    assert expected == got


def test_show_conf_basic(clean_ui):
    """Set up a very basic data/model/fit"""

    ui.load_arrays(1, [1, 2, 4], [3, 5, 5])
    ui.set_source(ui.scale1d.mdl)
    ui.fit()
    ui.conf()

    out = StringIO()
    ui.show_conf(outfile=out)
    got = out.getvalue().split('\n')

    assert len(got) == 12
    assert got[0] == "Confidence:Dataset               = 1"
    assert got[1] == "Confidence Method     = confidence"
    assert got[2] == "Iterative Fit Method  = None"
    assert got[3] == "Fitting Method        = levmar"
    assert got[4] == "Statistic             = chi2gehrels"
    assert got[5] == "confidence 1-sigma (68.2689%) bounds:"
    assert got[6] == "   Param            Best-Fit  Lower Bound  Upper Bound"
    assert got[7] == "   -----            --------  -----------  -----------"
    assert got[8] == "   mdl.c0            4.19798     -1.85955      1.85955"
    assert got[9] == ""
    assert got[10] == ""
    assert got[11] == ""


@pytest.mark.parametrize("string", [True, False])
def test_freeze_parameter(string, clean_ui):
    """Can we freeze a parameter?"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    assert not mdl.c1.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    if string:
        ui.freeze("mdl.c1")
    else:
        ui.freeze(mdl.c1)

    assert mdl.c1.frozen
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 2


@pytest.mark.parametrize("string", [True, False])
def test_freeze_frozen_parameter(string, clean_ui):
    """Can we freeze a frozen_parameter?"""

    mdl = ui.create_model_component("polynom1d", "mdl")
    ui.set_source(mdl)
    assert mdl.c1.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 9

    if string:
        ui.freeze("mdl.c1")
    else:
        ui.freeze(mdl.c1)

    assert mdl.c1.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 9


@pytest.mark.parametrize("string", [True, False])
def test_freeze_alwaysfrozen_parameter(string, clean_ui):
    """Can we freeze an always-frozen parameter?"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    assert mdl.ref.alwaysfrozen
    assert mdl.ref.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    if string:
        ui.freeze("mdl.ref")
    else:
        ui.freeze(mdl.ref)

    assert mdl.ref.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1


@pytest.mark.parametrize("string", [True, False])
def test_thaw_parameter(string, clean_ui):
    """Can we thaw a parameter?"""

    mdl = ui.create_model_component("polynom1d", "mdl")
    ui.set_source(mdl)
    assert mdl.c1.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 9

    if string:
        ui.thaw("mdl.c1")
    else:
        ui.thaw(mdl.c1)

    assert not mdl.c1.frozen
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 8


@pytest.mark.parametrize("string", [True, False])
def test_thaw_thawed_parameter(string, clean_ui):
    """Can we thaw a thawed parameter? String argument"""

    mdl = ui.create_model_component("polynom1d", "mdl")
    ui.set_source(mdl)
    assert not mdl.c0.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 9

    if string:
        ui.thaw("mdl.c0")
    else:
        ui.thaw(mdl.c0)

    assert not mdl.c0.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 9


@pytest.mark.parametrize("string", [True, False])
def test_thaw_alwaysfrozen_parameter(string, clean_ui):
    """Can we thaw an always-frozen parameter?"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    assert mdl.ref.alwaysfrozen
    assert mdl.ref.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    with pytest.raises(ParameterErr) as pe:
        if string:
            ui.thaw("mdl.ref")
        else:
            ui.thaw(mdl.ref)

    assert mdl.ref.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    assert str(pe.value) == "parameter mdl.ref is always frozen and cannot be thawed"


@pytest.mark.parametrize("string", [True, False])
def test_freeze_model(string, clean_ui):
    """Can we freeze a model?

    We use a model with an alwaysfrozen parameter.
    """

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    assert mdl.ref.frozen
    assert not mdl.c1.frozen
    assert not mdl.c2.frozen
    assert not mdl.ampl.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    if string:
        ui.freeze("mdl")
    else:
        ui.freeze(mdl)

    assert mdl.ref.frozen
    assert mdl.c1.frozen
    assert mdl.c2.frozen
    assert mdl.ampl.frozen
    assert ui.get_num_par_thawed() == 0
    assert ui.get_num_par_frozen() == 4


@pytest.mark.parametrize("string", [True, False])
def test_thaw_model(string, clean_ui):
    """Can we thaw a model?

    We use a model with an alwaysfrozen parameter.
    """

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    mdl.c1.freeze()
    mdl.ampl.freeze()

    assert mdl.ref.frozen
    assert mdl.c1.frozen
    assert not mdl.c2.frozen
    assert mdl.ampl.frozen
    assert ui.get_num_par_thawed() == 1
    assert ui.get_num_par_frozen() == 3

    if string:
        ui.thaw("mdl")
    else:
        ui.thaw(mdl)

    assert mdl.ref.frozen
    assert not mdl.c1.frozen
    assert not mdl.c2.frozen
    assert not mdl.ampl.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1


@pytest.mark.parametrize("string", [True, False])
def test_freeze_multi_arguments(string, clean_ui):
    """Check we can combine model and parameters"""

    mdl1 = ui.create_model_component("logparabola", "mdl1")
    mdl2 = ui.create_model_component("polynom1d", "mdl2")
    ui.set_source(mdl1 + mdl2)
    assert not mdl1.c2.frozen
    assert not mdl2.c0.frozen
    assert ui.get_num_par_thawed() == 4
    assert ui.get_num_par_frozen() == 10

    if string:
        ui.freeze("mdl2", "mdl1.c2")
    else:
        ui.freeze(mdl2, mdl1.c2)

    assert mdl1.c2.frozen
    assert mdl2.c0.frozen
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 12


@pytest.mark.parametrize("string", [True, False])
def test_thaw_multi_arguments(string, clean_ui):
    """Check we can combine model and parameters"""

    mdl1 = ui.create_model_component("logparabola", "mdl1")
    mdl2 = ui.create_model_component("polynom1d", "mdl2")
    ui.set_source(mdl1 + mdl2)
    mdl1.c1.freeze()
    mdl1.ampl.freeze()
    assert mdl1.c1.frozen
    assert mdl1.ampl.frozen
    assert mdl2.c2.frozen
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 12

    if string:
        ui.thaw("mdl1", "mdl2.c2")
    else:
        ui.thaw(mdl1, mdl2.c2)

    assert not mdl1.c1.frozen
    assert not mdl1.ampl.frozen
    assert not mdl2.c2.frozen
    assert ui.get_num_par_thawed() == 5
    assert ui.get_num_par_frozen() == 9


def test_freeze_no_arguments(clean_ui):
    """This is a no-op"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1

    ui.freeze()
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1


def test_thaw_no_arguments(clean_ui):
    """This is a no-op"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)
    mdl.c1.freeze()
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 2

    ui.thaw()
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 2


@pytest.mark.parametrize("string", [True, False])
def test_freeze_invalid_arguments(string, clean_ui):
    """We error out with an invalid argument"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)

    with pytest.raises(ArgumentTypeErr) as ae:
        if string:
            ui.freeze("1")
        else:
            ui.freeze(1)

    assert str(ae.value) == "'par' must be a parameter or model object or expression string"


@pytest.mark.parametrize("string", [True, False])
def test_thaw_invalid_arguments(string, clean_ui):
    """We error out with an invalid argument"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)

    with pytest.raises(ArgumentTypeErr) as ae:
        if string:
            ui.thaw("1")
        else:
            ui.thaw(1)

    assert str(ae.value) == "'par' must be a parameter or model object or expression string"
