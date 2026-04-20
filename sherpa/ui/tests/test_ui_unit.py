#
#  Copyright (C) 2017, 2018, 2020-2025
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

from sherpa.fit import ErrorEstResults, FitResults
from sherpa.models.parameter import Parameter
from sherpa.models.model import ArithmeticModel
from sherpa import ui
from sherpa.utils.err import ArgumentTypeErr, DataErr, IdentifierErr, ParameterErr
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

    with pytest.raises(ArgumentTypeErr,
                       match="'ids' must be an identifier or list of identifiers"):
        func(None, lo)


@pytest.mark.parametrize("func", [ui.notice, ui.ignore])
@pytest.mark.parametrize("lo,hi", [(1, 5), (1, None), (None, 5), (None, None),
                                   ("1:5", None)])
def test_filter_no_data_is_an_error(func, lo, hi, clean_ui):
    """Does applying a filter lead to an error?

    This test was added because it was noted that an update for Python 3
    had lead to an error state not being reached.
    """

    with pytest.raises(IdentifierErr,
                       match="No data sets found"):
        func(lo, hi)


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
    with pytest.raises(ArgumentTypeErr,
                       match="^'meth' must be a string$"):
        ui.set_iter_method(23)


def test_set_iter_method_type_not_enumeration():
    with pytest.raises(TypeError,
                       match="^a random string is not an iterative fitting method$"):
        ui.set_iter_method('a random string')


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
    with pytest.raises(ArgumentTypeErr,
                       match="'ids' must be an identifier or list of identifiers"):
        func(ids)


def test_set_model_autoassign_func_type():

    with pytest.raises(ArgumentTypeErr,
                       match="'func' must be a function or other callable object"):
        ui.set_model_autoassign_func(23)


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
    assert lvl == logging.WARNING
    assert msg == "No guess found for dummy"


def test_guess_warns_no_guess_no_argument(caplog, clean_ui):
    """Do we warn when the (implied) model has no guess"""

    ui.load_arrays(1, [1, 2, 3], [-3, 4, 5])
    cpt = DummyModel('dummy')
    ui.set_source(cpt + cpt)

    assert len(caplog.records) == 0
    ui.guess()

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.models.model"
    assert lvl == logging.WARNING
    assert msg == "No guess found for dummy"


def test_guess_single_model(caplog, clean_ui):
    """What happens when call guess on a single model?"""

    # Pick a simple model in terms of it's guess behaviour.
    ui.load_arrays(1, [10, 20, 40], [5, 7, 3])
    ui.set_source(ui.scale1d.mdl)

    bound = 100000

    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min < -bound # we just care these change
    assert mdl.c0.max > bound

    ui.guess()

    assert mdl.c0.val == pytest.approx(5)
    assert mdl.c0.min == pytest.approx(0.03)
    assert mdl.c0.max == pytest.approx(700)

    assert len(caplog.records) == 0


def test_guess_multiple_model_single(caplog, clean_ui):
    """What happens when call guess on a single model from a composite?"""

    # Pick a simple model in terms of it's guess behaviour.
    ui.load_arrays(1, [10, 20, 40], [5, 7, 3])
    ui.set_source(ui.scale1d.mdl + ui.box1d.bmdl)

    bound = 100000

    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min < -bound # we just care these change
    assert mdl.c0.max > bound

    # check this doesn't get changhed
    assert bmdl.ampl.val == pytest.approx(1)
    assert bmdl.ampl.min < -bound
    assert bmdl.ampl.max > bound

    ui.guess(mdl)

    assert mdl.c0.val == pytest.approx(5)
    assert mdl.c0.min == pytest.approx(0.03)
    assert mdl.c0.max == pytest.approx(700)

    assert bmdl.ampl.val == pytest.approx(1)
    assert bmdl.ampl.min < -1000
    assert bmdl.ampl.max > 1000

    assert len(caplog.records) == 0


def test_guess_multiple_model_multiple(caplog, clean_ui):
    """What happens when call guess on multiple models from a composite?"""

    # Pick a simple model in terms of it's guess behaviour.
    ui.load_arrays(1, [10, 20, 40], [5, 7, 3])
    ui.set_source(ui.scale1d.mdl + ui.box1d.bmdl)

    bound = 100000

    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min < -bound # we just care these change
    assert mdl.c0.max > bound

    # check this doesn't get changhed
    assert bmdl.ampl.val == pytest.approx(1)
    assert bmdl.ampl.min < -bound
    assert bmdl.ampl.max > bound

    ui.guess()

    assert mdl.c0.val == pytest.approx(5)
    assert mdl.c0.min == pytest.approx(0.03)
    assert mdl.c0.max == pytest.approx(700)

    assert bmdl.ampl.val == pytest.approx(7)
    assert bmdl.ampl.min == pytest.approx(0.007)
    assert bmdl.ampl.max == pytest.approx(7000)

    assert len(caplog.records) == 0


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
    with pytest.raises(ParameterErr,
                       match="parameter 'mdl.c1' is frozen"):
        method(mdl.c0, mdl.c1)


@pytest.mark.parametrize("method", [ui.conf, ui.covar, ui.proj])
def test_err_estimate_errors_model_all_frozen(method, clean_ui):
    """Check we error out with frozen model with conf/proj/covar.

    """

    ui.load_arrays(1, [1, 2, 3], [1, 2, 3])
    ui.set_source(ui.polynom1d.mdl)
    for par in mdl.pars:
        par.freeze()

    with pytest.raises(ParameterErr,
                       match="Model 'polynom1d.mdl' has no thawed parameters"):
        method(mdl)


@pytest.mark.parametrize("method", [ui.conf, ui.covar, ui.proj])
@pytest.mark.parametrize("id", [1, "xx"])
@pytest.mark.parametrize("otherids", [[2, 3], ["foo", "bar"]])
def test_err_estimate_errors_on_list_argument(method, id, otherids, clean_ui):
    """Check we error out with a list argument with conf/proj/covar.

    We had documented that you could say conf(1, [2, 3]) but this
    is not true, so check it does error out. Fortunately we can
    do this without setting up any dataset or model.

    """

    with pytest.raises(ArgumentTypeErr,
                       match="identifiers must be integers or strings"):
        method(id, otherids)


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


def test_show_conf_basic(clean_ui, check_str):
    """Set up a very basic data/model/fit"""

    ui.load_arrays(1, [1, 2, 4], [3, 5, 5])
    ui.set_source(ui.scale1d.mdl)
    ui.fit()
    ui.conf()

    out = StringIO()
    ui.show_conf(outfile=out)
    got = out.getvalue()

    check_str(got,
              ["Confidence:Dataset               = 1",
               "Confidence Method     = confidence",
               "Fitting Method        = levmar",
               "Statistic             = chi2gehrels",
               "confidence 1-sigma (68.2689%) bounds:",
               "   Param            Best-Fit  Lower Bound  Upper Bound",
               "   -----            --------  -----------  -----------",
               "   mdl.c0            4.19798     -1.85955      1.85955",
               "",
               "",
               ""
               ])


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

    with pytest.raises(ParameterErr,
                       match="^parameter mdl.ref is always frozen and cannot be thawed$"):
        if string:
            ui.thaw("mdl.ref")
        else:
            ui.thaw(mdl.ref)

    assert mdl.ref.frozen
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1


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

    with pytest.raises(ArgumentTypeErr,
                       match="^'par' must be a parameter or model object or expression string$"):
        if string:
            ui.freeze("1")
        else:
            ui.freeze(1)


@pytest.mark.parametrize("string", [True, False])
def test_thaw_invalid_arguments(string, clean_ui):
    """We error out with an invalid argument"""

    mdl = ui.create_model_component("logparabola", "mdl")
    ui.set_source(mdl)

    with pytest.raises(ArgumentTypeErr,
                       match="^'par' must be a parameter or model object or expression string$"):
        if string:
            ui.thaw("1")
        else:
            ui.thaw(1)


def test_set_dep_none(clean_ui):
    """What happens if set_dep is called with None?"""

    x = np.arange(10, 30, 5)
    y = np.ones(4)
    ui.load_arrays(1, x, y)

    ui.set_dep(None)

    # We get [nan, nan, nan, nan]
    assert np.isnan(ui.get_dep()) == pytest.approx([1, 1, 1, 1])


def test_set_dep_scalar(clean_ui):
    """What happens if set_dep is called with a scalar?"""

    x = np.arange(10, 30, 5)
    ones = np.ones(4)
    ui.load_arrays(1, x, ones)

    ui.set_dep(3)
    assert ui.get_dep() == pytest.approx(3 * ones)


def test_set_dep_scalar_empty(clean_ui):
    """What happens if set_dep is called with a scalar but no data

    This is a corner case (since we may decide to remove Data arrays
    with no data).
    """

    ui.set_data(ui.Data1D("x", None, None))
    with pytest.raises(TypeError):
        # The message is some variant of:
        # "can't multiply sequence by non-int of type 'NoneType'"
        ui.set_dep(3)


def test_set_dep_array(clean_ui):
    """What happens if set_dep is called with an array?"""

    x = np.arange(10, 30, 5)
    ones = np.ones(4)
    ui.load_arrays(1, x, ones)

    ui.set_dep([2, 4, 5, 19])
    assert ui.get_dep() == pytest.approx([2, 4, 5, 19])


def test_set_dep_array_wrong(clean_ui):
    """What happens if set_dep is called with an array with the wrong length?"""

    x = np.arange(10, 30, 5)
    ones = np.ones(4)
    ui.load_arrays(1, x, ones)

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and y: 4 vs 3"):
        ui.set_dep([2, 4, 5])


def test_set_staterror_none(clean_ui):
    """What happens when we set the staterror to None?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    ui.load_arrays(1, np.arange(3), np.ones(3), staterror, syserror)
    ui.set_stat('cstat')

    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)

    # removing the statistical error means that the statistic is used;
    # for the likelihood stats we just get 1's
    #
    ui.set_staterror(None)
    assert ui.get_staterror() == pytest.approx(np.ones(3))
    assert ui.get_syserror() == pytest.approx(syserror)

    combo = np.sqrt(1 + 0.25) * np.ones(3)
    assert ui.get_error() == pytest.approx(combo)


def test_set_syserror_none(clean_ui):
    """What happens when we set the syserror to None?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    ui.load_arrays(1, np.arange(3), np.ones(3), staterror, syserror)
    ui.set_stat('cstat')

    ui.set_syserror(None)
    assert ui.get_staterror() == pytest.approx(staterror)
    with pytest.raises(DataErr,
                       match="data set '1' does not specify systematic errors"):
        ui.get_syserror()

    combo = staterror
    assert ui.get_error() == pytest.approx(combo)


def test_set_staterror_scalar_no_fractional(clean_ui):
    """What happens when we set the staterror to a scalar fractional=False?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    ui.load_arrays(1, np.arange(3), np.ones(3), staterror, syserror)
    ui.set_stat('cstat')

    ui.set_staterror(3)
    assert ui.get_staterror() == pytest.approx(3 * np.ones(3))
    assert ui.get_syserror() == pytest.approx(syserror)

    combo = np.sqrt(9 + 0.25) * np.ones(3)
    assert ui.get_error() == pytest.approx(combo)


def test_set_syserror_scalar_no_fractional(clean_ui):
    """What happens when we set the syserror to a scalar fractional=False?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    ui.load_arrays(1, np.arange(3), np.ones(3), staterror, syserror)
    ui.set_stat('cstat')

    ui.set_syserror(3)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(3 * np.ones(3))

    combo = np.sqrt(0.01 + 9) * np.ones(3)
    assert ui.get_error() == pytest.approx(combo)


def test_set_staterror_scalar_fractional(clean_ui):
    """What happens when we set the staterror to a scalar fractional=True?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, np.arange(3), y, staterror, syserror)
    ui.set_stat('cstat')

    ui.set_staterror(0.4, fractional=True)
    assert ui.get_staterror() == pytest.approx(0.4 * y)
    assert ui.get_syserror() == pytest.approx(syserror)

    combo = np.sqrt(0.16 * y * y + 0.25)
    assert ui.get_error() == pytest.approx(combo)


def test_set_syserror_scalar_fractional(clean_ui):
    """What happens when we set the syserror to a scalar fractional=True?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, np.arange(3), y, staterror, syserror)
    ui.set_stat('cstat')

    ui.set_syserror(0.4, fractional=True)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(0.4 * y)

    combo = np.sqrt(0.01 + 0.16 * y * y)
    assert ui.get_error() == pytest.approx(combo)


def test_set_staterror_array(clean_ui):
    """What happens when we set the staterror to an array?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, np.arange(3), y, None, syserror)
    ui.set_stat('cstat')

    ui.set_staterror(staterror)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)


def test_set_syserror_array(clean_ui):
    """What happens when we set the syserror to an array?"""

    staterror = 0.1 * np.ones(3)
    syserror = 0.5 * np.ones(3)
    combo = np.sqrt(0.01 + 0.25) * np.ones(3)

    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, np.arange(3), y, staterror, None)
    ui.set_stat('cstat')

    ui.set_syserror(syserror)
    assert ui.get_staterror() == pytest.approx(staterror)
    assert ui.get_syserror() == pytest.approx(syserror)
    assert ui.get_error() == pytest.approx(combo)


@pytest.mark.parametrize("field", ["staterror", "syserror"])
def test_set_error_array_wrong(field, clean_ui):
    """What happens when we set the stat/syserror to an array of the wrong length?"""

    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, np.arange(3), y)
    ui.set_stat('cstat')

    setfunc = getattr(ui, f"set_{field}")

    with pytest.raises(DataErr,
                       match=f"size mismatch between independent axis and {field}: 3 vs 4"):
        setfunc(np.asarray([1, 2, 3, 4]))


@pytest.mark.parametrize("fractional", [False, True])
def test_set_staterr_scalar_nodep(fractional, clean_ui):
    """Corner case: set_staterror with scalar and no data."""

    ui.set_data(ui.Data1D("x", None, None))
    with pytest.raises(TypeError):
        # The errors are some variant of:
        # "can't multiply sequence by non-int of type 'NoneType'"
        # "unsupported operand type(s) for *: 'float' and 'NoneType'"
        ui.set_staterror(2.3, fractional=fractional)


@pytest.mark.parametrize("ignore", [False, True])
def test_set_filter_unmasked(ignore, clean_ui):
    """What happens when we call set_filter to an unfiltered dataset?"""

    x = np.asarray([10, 20, 30])
    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, x, y)

    data = ui.get_data()
    assert data.mask is True

    if ignore:
        expected = [False, True, False]
    else:
        expected = [True, False, True]

    ui.set_filter(np.asarray([True, False, True]), ignore=ignore)
    assert data.mask == pytest.approx(np.asarray(expected))


def test_set_filter_unmasked_wrong(clean_ui):
    """What happens when we call set_filter to an unfiltered dataset with the wrong size?"""

    x = np.asarray([10, 20, 30])
    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, x, y)

    data = ui.get_data()

    with pytest.raises(DataErr,
                       match="size mismatch between independent axis and mask: 3 vs 2"):
        ui.set_filter(np.asarray([True, False]))


@pytest.mark.parametrize("ignore", [False, True])
def test_set_filter_masked(ignore, clean_ui):
    """What happens when we call set_filter to a filtered dataset?"""

    x = np.asarray([10, 20, 30, 40, 50])
    y = np.asarray([2, 3, 4, 5, 6])
    ui.load_arrays(1, x, y)

    ui.ignore(lo=15, hi=45)

    data = ui.get_data()

    orig = np.asarray([True, False, False, False, True])
    assert data.mask == pytest.approx(orig)

    # Unlike test_set_filter_masked the two expected values are not
    # logical inverses, since we need to consider the existing mask.
    #
    if ignore:
        expected = [False, False, False, False, True]
    else:
        expected = [True, False, True, False, True]

    ui.set_filter(np.asarray([True, False, True, False, False]), ignore=ignore)
    assert data.mask == pytest.approx(np.asarray(expected))


def test_set_filter_masked_wrong(clean_ui):
    """What happens when we call set_filter to a filtered dataset with the wrong size?"""

    x = np.asarray([10, 20, 30])
    y = np.asarray([2, 3, 4])
    ui.load_arrays(1, x, y)

    ui.ignore(lo=15, hi=45)

    with pytest.raises(DataErr,
                       match="size mismatch between data and filter: 3 vs 2"):
        ui.set_filter(np.asarray([True, False]))


def test_fit_output_single(clean_ui, caplog, check_str):
    """This is essentially the test from #1804 but for a single dataset"""

    ui.load_arrays(1, [1, 2, 3], [12, 10, 4])
    ui.set_source(ui.const1d.mdl)
    ui.set_stat("cash")

    assert len(caplog.records) == 0
    ui.fit()
    assert len(caplog.records) == 1

    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    check_str(msg,
              ["Dataset               = 1",
               "Method                = levmar",
               "Statistic             = cash",
               "Initial fit statistic = 6",
               "Final fit statistic   = -60.2932 at function evaluation 11",
               "Data points           = 3",
               "Degrees of freedom    = 2",
               "Change in statistic   = 66.2932",
               "   mdl.c0         8.66665      +/- 1.70862     "])


def test_calc_stat_info_output_single(clean_ui, caplog, check_str):
    """This is essentially the test from #1804 but for a single dataset
    and converted to check calc_stat

    """

    ui.load_arrays(1, [1, 2, 3], [12, 10, 4])
    ui.set_source(ui.const1d.mdl)
    ui.set_stat("cash")

    assert len(caplog.records) == 0
    ui.fit()
    assert len(caplog.records) == 1
    ui.calc_stat_info()
    assert len(caplog.records) == 2

    lname, lvl, msg = caplog.record_tuples[1]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    check_str(msg,
              ["Dataset               = 1",
               "Statistic             = cash",
               "Fit statistic value   = -60.2932",
               "Data points           = 3",
               "Degrees of freedom    = 2"])


def test_fit_output_multi(clean_ui, caplog, check_str):
    """This is essentially the test from #1804"""

    ui.load_arrays(1, [1, 2, 3], [12, 10, 4])
    ui.load_arrays(2, [1, 2, 3], [7, 8, 6])
    ui.set_source(ui.const1d.mdl)
    ui.set_source(2, mdl)
    ui.set_stat("cash")

    assert len(caplog.records) == 0
    ui.fit()
    assert len(caplog.records) == 1

    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    check_str(msg,
              ["Datasets              = 1, 2",
               "Method                = levmar",
               "Statistic             = cash",
               "Initial fit statistic = 12",
               "Final fit statistic   = -99.4885 at function evaluation 11",
               "Data points           = 6",
               "Degrees of freedom    = 5",
               "Change in statistic   = 111.488",
               "   mdl.c0         7.83332      +/- 1.14642     "])


def test_calc_stat_info_output_multi(clean_ui, caplog, check_str):
    """This is essentially the test from #1804 but for calc_stat_info"""

    ui.load_arrays(1, [1, 2, 3], [12, 10, 4])
    ui.load_arrays(2, [1, 2, 3], [7, 8, 6])
    ui.set_source(ui.const1d.mdl)
    ui.set_source(2, mdl)
    ui.set_stat("cash")

    assert len(caplog.records) == 0
    ui.fit()
    assert len(caplog.records) == 1
    ui.calc_stat_info()
    assert len(caplog.records) == 2

    lname, lvl, msg = caplog.record_tuples[1]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    check_str(msg,
              ["Dataset               = 1",
               "Statistic             = cash",
               "Fit statistic value   = -60.0362",
               "Data points           = 3",
               "Degrees of freedom    = 2",
               "",
               "Dataset               = 2",
               "Statistic             = cash",
               "Fit statistic value   = -39.4523",
               "Data points           = 3",
               "Degrees of freedom    = 2",
               "",
               "Datasets              = [1, 2]",
               "Statistic             = cash",
               "Fit statistic value   = -99.4885",
               "Data points           = 6",
               "Degrees of freedom    = 5"])


def test_num_pars_single(clean_ui):
    """Related to issue #777 and #1982.

    Let's just check what we get
    """

    gmdl = ui.create_model_component("gauss1d", "gmdl")
    gmdl.ampl.freeze()
    ui.set_source(gmdl)

    assert ui.get_num_par() == 3
    assert ui.get_num_par_thawed() == 2
    assert ui.get_num_par_frozen() == 1


def test_num_pars_combined(clean_ui):
    """Related to issue #777 and #1982.

    Let's just check what we get
    """

    gmdl = ui.create_model_component("gauss1d", "gmdl")
    gmdl.ampl.freeze()
    bmdl = ui.create_model_component("scale1d", "bmdl")
    ui.set_source(gmdl + bmdl)

    assert ui.get_num_par() == 4
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 1


def test_num_pars_duplicated(clean_ui):
    """Related to issue #777 and #1982.

    Let's just check what we get
    """

    # For fun we repeat the gmdl component in the model expression. As
    # written it does not make much sense but it is a simplified form
    # for more-complex cases where it might end up having the same
    # component appear multiple times.
    #
    gmdl = ui.create_model_component("gauss1d", "gmdl")
    gmdl.ampl.freeze()
    bmdl = ui.create_model_component("scale1d", "bmdl")
    ui.set_source(gmdl + bmdl + gmdl)

    assert ui.get_num_par() == 7
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 4


def test_num_pars_links(clean_ui):
    """Related to issue #777 and #1982.

    Let's just check what we get
    """

    # Have a linked parameter that is not part of the model
    # expression.
    #
    gmdl = ui.create_model_component("gauss1d", "gmdl")
    gmdl.ampl.freeze()
    bmdl = ui.create_model_component("scale1d", "bmdl")
    smdl = ui.create_model_component("scale1d", "sigma")
    gmdl.fwhm = 2.4 * smdl.c0  # approximate the sigma/FWHM scaling here
    ui.set_source(gmdl + bmdl)

    # The linked parameter is included here.
    assert ui.get_num_par() == 5
    assert ui.get_num_par_thawed() == 3
    assert ui.get_num_par_frozen() == 2


SEEDVAL = 1327

def set_rng_global():
    # Ensure no RNG is set
    ui.set_rng(None)
    np.random.seed(SEEDVAL)


def set_rng_local():
    # Set the global RNG state to a random state to make sure it isn't
    # being used (the assumption is that this will not match SEEDVAL).
    np.random.seed()
    ui.set_rng(np.random.RandomState(SEEDVAL))


# The uniform sampling creates a model value < 0 which is allowed,
# but unphysical. Really should have changed the parameter range.
#
N_SAMPLE = np.asarray([[2.9013215,  4.4862895 ],
                       [1.55735165, 3.86432178],
                       [0.42363933, 2.35716076]])
U_SAMPLE = np.asarray([[ 5.97118671e-01,  3.20748789e+00],
                       [ 9.39240153e-01,  3.48835118e+00],
                       [ 9.21034037e+02, -2.67949604e-01]])
T_SAMPLE = np.asarray([[9.75190011, 6.69601813],
                       [8.60285233, 6.37242034],
                       [0.46611989, 2.30776257]])

@pytest.mark.parametrize("xxx,expected,setrng",
                         [("normal", N_SAMPLE, set_rng_global),
                          ("normal", N_SAMPLE, set_rng_local),
                          ("uniform", U_SAMPLE, set_rng_global),
                          ("uniform", U_SAMPLE, set_rng_local),
                          ("t", T_SAMPLE, set_rng_global),
                          ("t", T_SAMPLE, set_rng_local),
                         ])
def test_xxx_sample_random(xxx, expected, setrng, clean_ui):
    """Check if xxx_sample is repeatable.

    At the moment we allow the global RNG (np.random.seed)

    """

    ui.load_arrays(1, [2, 3, 10], [3, 4, 1])
    ui.set_stat('cash')
    ui.set_source(ui.const1d.mdl)
    ui.fit()

    setrng()

    func = getattr(ui, f"{xxx}_sample")
    answer = func(num=3)

    # clear out the RNG state
    np.random.seed()
    ui.set_rng(None)

    assert answer[:, 0:2] == pytest.approx(expected)
    assert answer[:, 2] == pytest.approx([0] * 3)


@pytest.mark.parametrize("setrng", [set_rng_global, set_rng_local])
def test_normal_sample_correlate(setrng, clean_ui):
    """Does the correlate setting change anything with normal_sample?"""

    ui.load_arrays(1, [1, 2, 3, 4, 5], [2, 4, 7, 15, 19])
    ui.set_source(ui.polynom1d.mdl)
    mdl.c1.thaw()
    ui.fit()

    bestfit = np.asarray([mdl.c0.val, mdl.c1.val])

    setrng()
    e1t = ui.normal_sample(num=3, scale=1, correlate=True)

    setrng()
    e1f = ui.normal_sample(num=3, scale=1, correlate=False)

    setrng()
    e2t = ui.normal_sample(num=3, scale=2, correlate=True)

    setrng()
    e2f = ui.normal_sample(num=3, scale=2, correlate=False)

    # Are the sigma=2 values twice those of sigma=1, after subtracting
    # off the best-fit? The first column is dropped as it is the
    # statistic column and the last column as it is the clipped column.
    #
    scalevalue = 2.0
    expected = np.full((3, 2), scalevalue)

    def check_ratio(v1, v2):
        """Check sigma=1 and sigma=2 results scale"""

        d1 = v1 - bestfit
        d2 = v2 - bestfit
        r = d2 / d1
        assert r == pytest.approx(expected)

    check_ratio(e1t[:, 1:-1], e2t[:, 1:-1])
    check_ratio(e1f[:, 1:-1], e2f[:, 1:-1])

    # Assume that the correlated=true/false results are different,
    # which can be checked by comparing the statistic columns.
    #
    stat_t = e1t[:, 0]
    stat_f = e1f[:, 0]
    assert not stat_f == pytest.approx(stat_t)


def setup_linked_pars(swap: bool = False
                      ) -> tuple[ArithmeticModel, ArithmeticModel]:
    """Create simple setup for test_normal_sample_linked_par_xxx

    If swap is False then the linked par version and the original
    have the same order of thawedpars, otherwise they are swapped.
    """

    # Fit the same data with the same model.
    #
    # This data is best-described with a poisson process. It was
    # generated from a model with c0=-3, c1=4.
    #
    x = [1, 2, 3, 4, 5]
    y = [1, 4, 12, 19, 12]

    ui.load_arrays(1, x, y)
    ui.set_source(1, ui.polynom1d.mdl1)
    mdl1.c1.thaw()

    ui.load_arrays(2, x, y)
    ui.set_source(2, ui.polynom1d.mdl2)
    mdl2.c1.thaw()

    linkedpar = ui.create_model_component("scale1d", "linkedpar")
    if swap:
        mdl2.c0 = linkedpar.c0
    else:
        mdl2.c1 = linkedpar.c0

    return mdl1, mdl2

def basic_linked_par_validation(mdl1: ArithmeticModel,
                                mdl2: ArithmeticModel,
                                f1: FitResults,
                                f2: FitResults,
                                c1: ErrorEstResults,
                                c2: ErrorEstResults,
                                swap: bool = False) -> None:
    """Simple checks that the fit results match for linked pars.

    This is just pulled out of the test to avoid complication, and is
    just a test that the linked-parameter fit is working as expected.

    swap indicates whether the thawed pars are in the same order
    in the two cases (False) or swapped (True).
    """

    # We have the same fit results
    assert f1.datasets == (1, )
    assert f2.datasets == (2, )
    assert f2.statval == pytest.approx(f1.statval)
    assert f2.dof == f1.dof

    assert len(mdl2.thawedpars) == len(mdl1.thawedpars)
    assert len(mdl1.get_thawed_pars()) == len(mdl2.get_thawed_pars())

    assert mdl1.get_thawed_pars()[0].fullname == "mdl1.c0"
    assert mdl1.get_thawed_pars()[1].fullname == "mdl1.c1"

    if swap:
        assert mdl2.get_thawed_pars()[0].fullname == "mdl2.c1"
    else:
        assert mdl2.get_thawed_pars()[0].fullname == "mdl2.c0"

    assert mdl2.get_thawed_pars()[1].fullname == "linkedpar.c0"

    if swap:
        assert mdl2.thawedpars[0] == pytest.approx(mdl1.c1.val)
        assert mdl2.thawedpars[1] == pytest.approx(mdl1.c0.val)
        assert linkedpar.c0.val == pytest.approx(mdl1.c0.val)
    else:
        assert mdl2.thawedpars[0] == pytest.approx(mdl1.c0.val)
        assert mdl2.thawedpars[1] == pytest.approx(mdl1.c1.val)
        assert linkedpar.c0.val == pytest.approx(mdl1.c1.val)

    assert c1.parnames == ('mdl1.c0', 'mdl1.c1')
    if swap:
        assert c2.parnames == ('mdl2.c1', 'linkedpar.c0')
    else:
        assert c2.parnames == ('mdl2.c0', 'linkedpar.c0')

    if swap:
        assert c2.parmins[0] == pytest.approx(c1.parmins[1])
        assert c2.parmins[1] == pytest.approx(c1.parmins[0])
        assert c2.parmaxes[0] == pytest.approx(c1.parmaxes[1])
        assert c2.parmaxes[1] == pytest.approx(c1.parmaxes[0])
    else:
        assert c2.parmins == pytest.approx(c1.parmins)
        assert c2.parmaxes == pytest.approx(c1.parmaxes)


@pytest.mark.parametrize("correlate", [False, True])
@pytest.mark.parametrize("method", ["levmar", "simplex"])
@pytest.mark.parametrize("swap", [False, True])
def test_normal_sample_linked_par_chisq(swap, method, correlate, clean_ui):
    """Does it work with a linked parameter?

    Although the data is best described by a Poisson process,
    use chi-square statistic for the fit.
    """

    mdl1, mdl2 = setup_linked_pars(swap=swap)

    # Since the error estimation may depend on the optimizer, make
    # sure we check.
    #
    ui.set_stat("chi2datavar")
    ui.set_method(method)

    ui.fit(1)
    f1 = ui.get_fit_results()

    ui.fit(2)
    f2 = ui.get_fit_results()

    ui.covar(1)
    c1 = ui.get_covar_results()

    ui.covar(2)
    c2 = ui.get_covar_results()

    basic_linked_par_validation(mdl1, mdl2, f1, f2, c1, c2,
                                swap=swap)

    set_rng_local()
    s1 = ui.normal_sample(id=1, num=5, correlate=correlate)

    set_rng_local()
    s2 = ui.normal_sample(id=2, num=5, correlate=correlate)

    if not swap:
        # Since the thwedpars should have the same order for the
        # original and linked-par version, the results should match.
        #
        assert s2 == pytest.approx(s1)
        return

    # Since the parameter order is swapped we do not expect
    # the same results for the two calls. There are some things
    # we can check however.
    #
    # TODO: somehow the statistic values actually match, which DJB
    # does not understand! Although we want to be rng_i * sigma_j away
    # from the best-fit location with the sampled parameter values,
    # the sigma_j values (which are per-model parameter) have nothing
    # to do with the calculation of the statistic, which is from
    # chi2datavar in this case, and calculated per-data point.
    #
    assert s2[:, 0] == pytest.approx(s1[:, 0])

    if not correlate:
        # Since the thawed parameters are swapped for the two datasets,
        # the sampled values will be different. However, since we used
        # correlate=False, so the samples are drawn from the sigma value
        # from the covar run - i.e. the cX.parmaxes values - we can
        # reverse calculate what the random values would have been,
        # and so we can check that these are as expected.
        #
        d1 = (s1[:, 1:-1] - mdl1.thawedpars) / c1.parmaxes
        d2 = (s2[:, 1:-1] - mdl2.thawedpars) / c2.parmaxes

        assert d2 == pytest.approx(d1)

    # Safety check: we should be at the best-fit location here.
    #
    orig1 = ui.calc_stat(1)
    orig2 = ui.calc_stat(2)
    assert orig2 == pytest.approx(orig1)
    assert orig1 == pytest.approx(f1.statval)

    # Apply the first row of the sample values and check the
    # statistic.
    #
    # This is to check that we get the actual statistic values and
    # that, even though the parameter values are different, the stats
    # are the same. As mentioned, DJB does not understand why this is
    # so, so if this starts to fail in the future this may not be a
    # bad thing.
    #
    mdl2.thawedpars = s2[0, 1:-1]
    check2 = ui.calc_stat(2)

    mdl1.thawedpars = s1[0, 1:-1]
    check1 = ui.calc_stat(1)

    # Just check they are not the same.
    assert mdl2.thawedpars != pytest.approx(mdl1.thawedpars)

    assert check2 == pytest.approx(check1)
    assert check1 == pytest.approx(s1[0, 0])

    # This location should be worse than the best-fit location.
    #
    assert check1 > f1.statval


class StraightLineU(ArithmeticModel):
    """y = mx + c where c >= 0"""

    def __init__(self, name="straight") -> None:
        self.c = Parameter(name, 'c', 1, 0, 100, hard_min=0, hard_max=1000)
        self.m = Parameter(name, 'm', 0, -10, 10, hard_min=-10, hard_max=10)
        super().__init__(name, (self.c, self.m))

    def calc(self, p, *args, **kwargs):
        """returm m * args[0] + c"""
        return self.m.val * (args[0] - 35) + self.c.val


class StraightLineL(ArithmeticModel):
    """y = mx + c where c <= 0"""

    def __init__(self, name="straight") -> None:
        self.c = Parameter(name, 'c', -1, -100, 0, hard_min=-1000, hard_max=0)
        self.m = Parameter(name, 'm', 0, -10, 10, hard_min=-10, hard_max=10)
        super().__init__(name, (self.c, self.m))

    def calc(self, p, *args, **kwargs):
        """returm m * args[0] + c"""
        return self.m.val * (args[0] - 35) + self.c.val


def setup_upper_limit() -> None:
    """Set up a fit where a parameter has an upper limit with covar."""

    # Create a dataset where we expect the c0 value to be an upper
    # limit. This is done by fitting y = mx + c to data where c has a
    # 1-sigma error that includes 0 and then setting the hard minimum
    # to 0.
    #
    # The data was simulated from polynom1d with
    #  c0 = 1  c1 = 0.5  offset=35
    #
    ui.load_arrays(1, [5, 10, 20, 30], [15.49, 14.02, 7.25, 3.60])
    ui.set_staterror(1, 1.2)

    # Note that the model has not been "asssociated" with the UI code,
    # but it should be okay to use.
    #
    mdl = StraightLineU(name="mdl")
    ui.set_source(mdl)

    ui.set_stat('chi2datavar')
    ui.set_method('levmar')

    ui.fit()


def setup_lower_limit() -> None:
    """Invert/negate the setup_upper_limit so we have a lower-limit parameter"""

    # This uses the negative values for the data used in lower_limit
    # so that c is now a lower-limit.
    #
    ui.load_arrays(1, [5, 10, 20, 30], [-15.49, -14.02, -7.25, -3.60])
    ui.set_staterror(1, 1.2)

    mdl = StraightLineL(name="mdl")
    ui.set_source(mdl)

    ui.set_stat('chi2datavar')
    ui.set_method('levmar')

    ui.fit()


def check_c(val: float, is_upper: bool) -> None:
    """Check the c value is as expected."""

    expected = 0.6730508474511385
    if not is_upper:
        expected *= -1

    assert val == pytest.approx(expected)


def check_c_limits(lo: float | None,
                   hi: float | None,
                   is_upper: bool) -> None:
    """Check the c limits."""

    sigma = 1.3163895152773046
    if is_upper:
        assert lo is None
        assert hi == pytest.approx(sigma)
    else:
        assert lo == pytest.approx(-sigma)
        assert hi is None


def check_m(val: float, is_upper: bool) -> None:
    """Check the m value is as expected."""

    expected = -0.5022372881358701
    if not is_upper:
        expected *= -1

    assert val == pytest.approx(expected)


def check_m_limits(lo: float | None,
                   hi: float | None) -> None:
    """Check the m limits."""

    sigma = 0.06249067727078889
    assert lo == pytest.approx(-sigma)
    assert hi == pytest.approx(sigma)


def validate_limits(is_upper: bool) -> None:
    """Check the upper-limit fit and errors are as expected.

    It is possible that these will change as the optimiser code
    changes, or numerical changes in OTS code.

    This code has been pulled out of the test since this is not
    the main part of the test, so we can separate out these
    checks logically without making the test code too complex.

    """

    # Check in expected location
    fr = ui.get_fit_results()
    assert fr.succeeded
    assert len(fr.parvals) == 2
    check_c(fr.parvals[0], is_upper)
    check_m(fr.parvals[1], is_upper)

    cr = ui.get_covar_results()
    assert len(cr.parmins) == 2
    check_c_limits(cr.parmins[0], cr.parmaxes[0], is_upper)
    check_m_limits(cr.parmins[1], cr.parmaxes[1])


def reset_seed() -> None:
    """Ensure the Sherpa RNG state is at a known point."""

    # Use randomState for repeatability as it is possible that
    # default_rng() could change over time.
    #
    ui.set_rng(np.random.RandomState(234032))


def test_normal_sample_2335_upper_limit(clean_ui, caplog):
    """Are hard limits respected in the normal_sample output?

    This is a regression test, just to see what happens.

    See issue #2335.
    """

    setup_upper_limit()
    ui.covar()
    validate_limits(True)

    mdl = ui.get_source()
    bestfit = mdl.thawedpars
    stat_bestfit = ui.calc_stat()

    # correlate=False
    reset_seed()
    nlog1 = len(caplog.records)
    res_false = ui.normal_sample(num=5, correlate=False, clip="none")
    nlog2 = len(caplog.records)

    # Check we get a message about switching to covariance and the
    # hard limits. The logged message may or may not contain "below
    # minimum" nessages, and it's not clear exactly why.
    #
    ngot = nlog2 - nlog1
    assert ngot in [5, 7]

    r = caplog.records[nlog1]
    assert r.getMessage() == "hard minimum hit for parameter mdl.c"
    r = caplog.records[nlog1 + 1]
    assert r.getMessage() == "Covariance failed for 'mdl.c', trying Confidence..."
    r = caplog.records[nlog1 + 2]
    assert r.getMessage().startswith("mdl.c lower bound:\t")
    r = caplog.records[nlog1 + 3]
    assert r.getMessage().startswith("mdl.c upper bound:\t")
    r = caplog.records[nlog1 + 4]
    assert r.getMessage() == "hard minimum hit for parameter mdl.c"

    if ngot == 7:
        msg = "value of parameter mdl.c is below minimum; setting to minimum"
        r = caplog.records[nlog1 + 5]
        assert r.getMessage() == msg
        r = caplog.records[nlog1 + 6]
        assert r.getMessage() == msg

    # Check the statistic values.
    #
    stat_false = [1.96756011,
                  5.17975686,
                  38.494826,
                  3.72584152,
                  12.42435527]
    assert res_false[:, 0] == pytest.approx(stat_false)

    # We expect two values where c is less than the hard limit.
    #
    c_false = [0.3972328,
               -1.18409209,
               2.39903626,
               -0.59492707,
               0.42605331]
    assert res_false[:, 1] == pytest.approx(c_false)

    # Check that the second statistic value was calculated with c=0
    # and not c=-1.18.
    #
    stat = res_false[1, 0]
    assert stat > stat_bestfit  # this is a worse fit than the best fit

    # Setting c to a value < 0 gets reset to 0.
    mdl.thawedpars = res_false[1, 1:-1]
    assert mdl.c.val == pytest.approx(0.0)

    # Using c=0 gets us the expected statistic
    assert ui.calc_stat() == pytest.approx(stat)

    # Reset the fit
    mdl.thawedpars = bestfit

    # correlate=True
    reset_seed()
    nlog3 = len(caplog.records)
    res_true = ui.normal_sample(num=5, correlate=True, clip="none")
    nlog4 = len(caplog.records)

    # Should have extra message(s) added to the log.
    #
    ngot43 = nlog4 - nlog3
    assert nlog3 == (nlog2 + 1)
    assert ngot43 in [1, 2]

    r = caplog.records[nlog3]
    assert r.getMessage() == "hard minimum hit for parameter mdl.c"

    assert res_true[:, 0] == pytest.approx([3.26772618,
                                            7.86906547,
                                            1.58810658,
                                            5.40581264,
                                            3.12139264])

    c_true = [0.95056442,
              -1.05177614,
              0.92072754,
              2.42574202,
              0.37431964]
    assert res_true[:, 1] == pytest.approx(c_true)

    # Since the statistic calculation follows the same code as with
    # the correlate=False case assume that it is correct and do not
    # test it here.

    # Now check with the clip=hard setting. The checks
    # are limited:
    #  - do we get the same statistic values
    #  - are the values clipped?
    #
    reset_seed()
    res_false2 = ui.normal_sample(num=5, correlate=False, clip="hard")

    assert res_false2[:, 0] == pytest.approx(stat_false)

    c_false2 = np.clip(c_false, 0, None)
    assert res_false2[:, 1] == pytest.approx(c_false2)

    reset_seed()
    res_true2 = ui.normal_sample(num=5, correlate=True, clip="hard")

    c_true2 = np.clip(c_true, 0, None)
    assert res_true2[:, 1] == pytest.approx(c_true2)

    assert res_false2[:, 3] == pytest.approx([0, 1, 0, 1, 0])
    assert res_true2[:, 3] == pytest.approx([0, 1, 0, 0, 0])


def test_normal_sample_2335_lower_limit(clean_ui, caplog):
    """Are hard limits respected in the normal_sample output?

    This is a regression test, just to see what happens.

    See test_normal_sample_2335_upper_limit. The test results are
    similar to the upper-limit version but not enough that it's worth
    bundling the checks into a common routine.

    """

    setup_lower_limit()
    ui.covar()
    validate_limits(False)

    # correlate=False
    reset_seed()
    nlog1 = len(caplog.records)
    res_false = ui.normal_sample(num=5, correlate=False, clip="none")
    nlog2 = len(caplog.records)

    # Check we get a message about switching to covariance and the
    # hard limits.
    #
    ngot = nlog2 - nlog1
    assert ngot in [6, 7]

    r = caplog.records[nlog1]
    assert r.getMessage() == "hard maximum hit for parameter mdl.c"
    r = caplog.records[nlog1 + 1]
    assert r.getMessage() == "Covariance failed for 'mdl.c', trying Confidence..."
    r = caplog.records[nlog1 + 2]
    assert r.getMessage().startswith("mdl.c lower bound:\t")
    r = caplog.records[nlog1 + 3]
    assert r.getMessage().startswith("mdl.c upper bound:\t")
    r = caplog.records[nlog1 + 4]
    assert r.getMessage() == "hard maximum hit for parameter mdl.c"
    r = caplog.records[nlog1 + 5]
    assert r.getMessage() == "1 sigma bounds for parameter mdl.c could not be found, using soft limit minimum"

    if ngot == 7:
        r = caplog.records[nlog1 + 6]
        assert r.getMessage() == "value of parameter mdl.c is above maximum; setting to maximum"

    # Check the statistic values. These are significantly worse
    # than the upper-limit case.
    #
    stat_false = [1145.16881,
                  54074.8610,
                  20.8524903,
                  25915.2235,
                  1153.55855]
    assert res_false[:, 0] == pytest.approx(stat_false)

    # We expect one value where c is larger than the hard limit.  Note
    # that these c values are a lot "worse" than in the upper-limit
    # case because the soft-limit minimum was used.
    #
    c_false = [-21.62566763,
               -141.75157795,
               130.44204801,
               -96.99544622,
               -19.43630701]
    assert res_false[:, 1] == pytest.approx(c_false)

    # correlate=True
    reset_seed()
    nlog3 = len(caplog.records)
    res_true = ui.normal_sample(num=5, correlate=True, clip="none")
    nlog4 = len(caplog.records)

    # Should have been one extra message added to the log.
    #
    assert nlog3 == nlog2
    ngot43 = nlog4 - nlog3
    assert ngot43 in [1, 2]

    r = caplog.records[nlog3]
    assert r.getMessage() == "hard maximum hit for parameter mdl.c"

    # Note that these are significantly "better" than the
    # correlate=False results, as the scaling used is
    # more-appropriate when using the covariance matrix.
    #
    assert res_true[:, 0] == pytest.approx([3.26772618,
                                            3.88042719,
                                            1.58810658,
                                            1.49768585,
                                            3.12139264])

    c_true = [-0.39553727,
              -2.39787783,
              -0.42537416,
              1.07964032,
              -0.97178205]
    assert res_true[:, 1] == pytest.approx(c_true)

    # Now check with the clip=hard setting. The checks
    # are limited:
    #  - do we get the same statistic values
    #  - are the values clipped?
    #
    reset_seed()
    res_false2 = ui.normal_sample(num=5, correlate=False, clip="hard")

    assert res_false2[:, 0] == pytest.approx(stat_false)

    c_false2 = np.clip(c_false, None, 0)
    assert res_false2[:, 1] == pytest.approx(c_false2)

    reset_seed()
    res_true2 = ui.normal_sample(num=5, correlate=True, clip="hard")

    c_true2 = np.clip(c_true, None, 0)
    assert res_true2[:, 1] == pytest.approx(c_true2)


def test_normal_sample_sigma_warning_message_upper(clean_ui, caplog):
    """Check we correctly report the sigma limit.

    Note that test_normal_sample_sigma_warning_message_upper
    checks the upper-limit version, as it is subtly-different.

    """

    setup_upper_limit()
    mdl = ui.get_source()

    # correlate=False
    reset_seed()
    nlog1 = len(caplog.records)
    resf_sigma1 = ui.normal_sample(num=5, correlate=False,
                                   clip="hard")
    nlog2 = len(caplog.records)

    reset_seed()
    resf_sigma16 = ui.normal_sample(num=5, scale=1.6, correlate=False,
                                    clip="hard")
    nlog3 = len(caplog.records)

    reset_seed()
    resf_sigma2 = ui.normal_sample(num=5, scale=2, correlate=False,
                                   clip="hard")
    nlog4 = len(caplog.records)

    reset_seed()
    rest_sigma1 = ui.normal_sample(num=5, correlate=True, clip="hard")
    nlog5 = len(caplog.records)

    reset_seed()
    rest_sigma16 = ui.normal_sample(num=5, scale=1.6, correlate=True,
                                    clip="hard")
    nlog6 = len(caplog.records)

    reset_seed()
    rest_sigma2 = ui.normal_sample(num=5, scale=2, correlate=True,
                                   clip="hard")
    nlog7 = len(caplog.records)

    # This does not have the "x sigma bounds not found" message,
    # unlike the _lower version.
    #
    ngot21 = nlog2 - nlog1
    ngot32 = nlog3 - nlog2
    ngot43 = nlog4 - nlog3
    ngot54 = nlog5 - nlog4
    ngot65 = nlog6 - nlog5
    ngot76 = nlog7 - nlog6
    assert ngot21 in [5, 7]
    assert ngot32 in [5, 7]
    assert ngot43 in [5, 7]
    assert ngot54 in [1, 2]
    assert ngot65 in [1, 2]
    assert ngot76 in [1, 2]

    # Check how the statistic values vary. As sigma increases you
    # might expect the statistic to increase. At the moment just
    # treat as a regression test.
    #
    assert resf_sigma1[:, 0] == pytest.approx([1.96756011, 5.17975686, 38.494826, 3.72584152, 12.42435527])
    assert resf_sigma16[:, 0] == pytest.approx([2.21202788, 8.96141072, 76.96112604, 4.15693761, 20.58014045])
    assert resf_sigma2[:, 0] == pytest.approx([1.93025625, 8.39272471, 95.45773131, 4.09849966, 20.96955128])

    assert rest_sigma1[:, 0] == pytest.approx([3.26772618, 7.86906547, 1.58810658, 5.40581264, 3.12139264])
    assert rest_sigma16[:, 0] == pytest.approx([7.31972735, 14.19471891, 2.47293905, 19.31488558, 6.98017826])
    assert rest_sigma2[:, 0] == pytest.approx([11.77709675, 12.75671961, 3.57858166, 39.53542143, 11.29312376])

    # The c values are hard to check for correlate=False, so just
    # check the offsets from the expected value.
    #
    d1f = resf_sigma1[:, 1] - mdl.c.val
    d16f = resf_sigma16[:, 1] - mdl.c.val
    d2f = resf_sigma2[:, 1] - mdl.c.val

    assert d1f == pytest.approx([-0.27581805, -0.6730508474511385, 1.72598541, -0.6730508474511385, -0.24699754])
    assert d16f == pytest.approx([-0.44130888, -0.6730508474511385, 2.76157666, -0.6730508474511385, -0.39519606])
    assert d2f == pytest.approx([-0.5516361, -0.6730508474511385, 3.45197083, -0.6730508474511385, -0.49399507])

    def check_m_ratio(vsigma, vsigma1, expected):
        r = vsigma / vsigma1
        assert r == pytest.approx(np.full(5, expected))

    d1f = resf_sigma1[:, 2] - mdl.m.val
    d16f = resf_sigma16[:, 2] - mdl.m.val
    d2f = resf_sigma2[:, 2] - mdl.m.val

    # TODO: Why are these ratios not 1.6 and 2, and why is sigma=2
    # smaller than sigma=1.6?
    #
    check_m_ratio(d16f, d1f, 1.2785750682547592)
    check_m_ratio(d2f, d1f, 1.2422080062609844)

    # For correlate=True, there is no need to special case the
    # c value. Since this uses the covariance matrix we check
    # whether the separation from the best-fit location has
    # increased with sigma.
    #
    d1t = abs(rest_sigma1[:, 1:-1] - mdl.thawedpars)
    d16t = abs(rest_sigma16[:, 1:-1] - mdl.thawedpars)
    d2t = abs(rest_sigma2[:, 1:-1] - mdl.thawedpars)

    def check_ratio(vsigma, vsigma1, minval, maxval):
        r = vsigma / vsigma1
        assert r.min() == pytest.approx(minval)
        assert r.max() == pytest.approx(maxval)
        # check the ratio for the c column
        assert (r[[0, 2, 3, 4], 0] > 1.6).all()
        # This is the out-of-bounds value
        assert r[1, 0] == pytest.approx(1.0)

    check_ratio(d16t, d1t, 0.5592601525645855, 2.38351205346464)
    check_ratio(d2t, d1t, 0.6046847691090992, 3.9596727228691004)


def test_normal_sample_sigma_warning_message_lower(clean_ui, caplog):
    """Check we correctly report the sigma limit.

    This is a tricky behaviour to test, so as
    test_normal_sample_2335_lower_limit triggers the warning
    message, use it for this test.

    Note that test_normal_sample_sigma_warning_message_upper
    checks the upper-limit version, as it is subtly-different.

    """

    setup_lower_limit()
    mdl = ui.get_source()

    # correlate=False
    reset_seed()
    nlog1 = len(caplog.records)
    resf_sigma1 = ui.normal_sample(num=5, correlate=False,
                                   clip="hard")
    nlog2 = len(caplog.records)

    reset_seed()
    resf_sigma16 = ui.normal_sample(num=5, scale=1.6, correlate=False,
                                    clip="hard")
    nlog3 = len(caplog.records)

    reset_seed()
    resf_sigma2 = ui.normal_sample(num=5, scale=2, correlate=False,
                                   clip="hard")
    nlog4 = len(caplog.records)

    reset_seed()
    rest_sigma1 = ui.normal_sample(num=5, correlate=True,
                                   clip="hard")
    nlog5 = len(caplog.records)

    reset_seed()
    rest_sigma16 = ui.normal_sample(num=5, scale=1.6, correlate=True,
                                    clip="hard")
    nlog6 = len(caplog.records)

    reset_seed()
    rest_sigma2 = ui.normal_sample(num=5, scale=2, correlate=True,
                                   clip="hard")
    nlog7 = len(caplog.records)

    ngot21 = nlog2 - nlog1
    ngot32 = nlog3 - nlog2
    ngot43 = nlog4 - nlog3
    ngot54 = nlog5 - nlog4
    ngot65 = nlog6 - nlog5
    ngot76 = nlog7 - nlog6
    assert ngot21 in [6, 7]
    assert ngot32 in [6, 7]
    assert ngot43 in [6, 7]
    assert ngot54 in [1, 2]
    assert ngot65 in [1, 2]
    assert ngot76 in [1, 2]

    r = caplog.records[nlog1 + 5]
    assert r.getMessage() == "1 sigma bounds for parameter mdl.c could not be found, using soft limit minimum"

    r = caplog.records[nlog2 + 5]
    assert r.getMessage() == "1.6 sigma bounds for parameter mdl.c could not be found, using soft limit minimum"

    r = caplog.records[nlog3 + 5]
    assert r.getMessage() == "2 sigma bounds for parameter mdl.c could not be found, using soft limit minimum"

    # Check how the statistic values vary.
    #
    assert resf_sigma1[:, 0] == pytest.approx([1.14516881e+03, 5.40748610e+04, 2.08524903e+01, 2.59152235e+04, 1.15355855e+03])
    assert resf_sigma16[:, 0] == pytest.approx([1.12467203e+03, 5.37400095e+04, 3.00878135e+01, 2.59548008e+04, 1.20528623e+03])
    assert resf_sigma2[:, 0] == pytest.approx([1.12733428e+03, 5.37836483e+04, 2.87802714e+01, 2.59496320e+04, 1.19845547e+03])

    assert rest_sigma1[:, 0] == pytest.approx([3.26772618, 3.88042719, 1.58810658, 1.49768585, 3.12139264])
    assert rest_sigma16[:, 0] == pytest.approx([7.31972735, 5.95784881, 2.47293905, 1.66356551, 6.98017826])
    assert rest_sigma2[:, 0] == pytest.approx([11.77709675, 7.91490779, 3.57858165, 4.16485202, 11.29312376])

    # The c values should be the same, as they are drawn from the same
    # distribution as this is a special case, when correlate=False.
    # We should be able to compare the m values with sigma (again for
    # the correlate=False).
    #
    assert resf_sigma16[:, 1] == pytest.approx(resf_sigma1[:, 1])
    assert resf_sigma2[:, 1] == pytest.approx(resf_sigma1[:, 1])

    def check_m_ratio(vsigma, vsigma1, expected):
        r = vsigma / vsigma1
        assert r == pytest.approx(np.full(5, expected))

    d1f = resf_sigma1[:, 2] - mdl.m.val
    d16f = resf_sigma16[:, 2] - mdl.m.val
    d2f = resf_sigma2[:, 2] - mdl.m.val

    # TODO: Why are these ratios not 1.6 and 2, and why is sigma=2
    # smaller than sigma=1.6?
    #
    check_m_ratio(d16f, d1f, 1.2785750685947328)
    check_m_ratio(d2f, d1f, 1.24220800667757)

    # For correlate=True, there is no need to special case the
    # c value. Since this uses the covariance matrix we check
    # whether the separation from the best-fit location has
    # increased with sigma.
    #
    d1t = abs(rest_sigma1[:, 1:-1] - mdl.thawedpars)
    d16t = abs(rest_sigma16[:, 1:-1] - mdl.thawedpars)
    d2t = abs(rest_sigma2[:, 1:-1] - mdl.thawedpars)

    def check_ratio(vsigma, vsigma1, minval, maxval):
        r = vsigma / vsigma1
        assert r.min() == pytest.approx(minval)
        assert r.max() == pytest.approx(maxval)
        # check the ratio for the c column
        assert (r[[0, 1, 2, 4], 0] > 1.6).all()
        # This is the out-of-bounds value
        assert r[3, 0] == pytest.approx(1.0)

    check_ratio(d16t, d1t, 0.5592601537154354, 2.3835120525988995)
    check_ratio(d2t, d1t, 0.6046847674045083, 3.9596727215871725)


def test_uniform_sample_upper_limit(clean_ui, caplog):
    """Are hard limits respected in the uniform_sample output?

    This is a regression test, just to see what happens.

    """

    setup_upper_limit()

    # clip=none
    reset_seed()
    res = ui.uniform_sample(num=5, clip="none")

    stat = [23.09504763, 248.14278362, 3.84258011, 4.10371269, 2.52533989]
    c = [-1.52105109,  5.79701772, -2.45911342, 0.20786994, -0.93543076]
    m = [-0.39829977, -0.71934196, -0.57431753, -0.56908599, -0.55958875]

    assert res[:, 0] == pytest.approx(stat)
    assert res[:, 1] == pytest.approx(c)
    assert res[:, 2] == pytest.approx(m)
    assert res[:, 3] == pytest.approx([0] * 5)

    # clip=hard
    reset_seed()
    res = ui.uniform_sample(num=5, clip="hard")

    c = np.clip(c, 0, None)
    assert res[:, 0] == pytest.approx(stat)
    assert res[:, 1] == pytest.approx(c)
    assert res[:, 2] == pytest.approx(m)
    assert res[:, 3] == pytest.approx([1, 0, 1, 0, 1])


def test_t_sample_upper_limit(clean_ui, caplog):
    """Are hard limits respected in the t_sample output?

    This is a regression test, just to see what happens.

    """

    setup_upper_limit()

    # clip=none
    reset_seed()
    res = ui.t_sample(num=5, clip="none")

    stat = [8.6840743, 50.73089515, 2.06343103, 10.54354892, 1.6621687]
    c = [1.20415568, -3.25097801, 1.05196097, 3.29119311, 0.53070364]
    m = [-0.55668802, -0.73053488, -0.51084555, -0.4575194, -0.48985982]

    assert res[:, 0] == pytest.approx(stat)
    assert res[:, 1] == pytest.approx(c)
    assert res[:, 2] == pytest.approx(m)

    # clip=hard
    reset_seed()
    res = ui.t_sample(num=5, clip="hard")

    c = np.clip(c, 0, None)
    assert res[:, 0] == pytest.approx(stat)
    assert res[:, 1] == pytest.approx(c)
    assert res[:, 2] == pytest.approx(m)


def setup_multiple_datasets():
    """Two datasets with the same model to test ordering."""

    ui.load_arrays(1, [10, 12, 14, 16, 18],
                   [34, 38, 40, 41, 39])
    ui.load_arrays(2, [20, 22, 24, 26],
                   [54, 51, 37, 64])

    ui.set_stat("cash")

    # Fit y = mx + c to the two datasets, with m = 2, c = 10
    #
    pol = ui.create_model_component("polynom1d", "pol")
    pol.c1.thaw()

    ui.set_source(1, pol)
    ui.set_source(2, pol)


def setup_multiple_datasets_linked():
    """Two datasets with linked parameters to test ordering."""

    ui.load_arrays(1, [10, 12, 14, 16, 18],
                   [41, 42, 37, 56, 52])
    ui.load_arrays(2, [20, 22, 24, 26],
                   [28, 22, 28, 38])

    ui.set_stat("cash")

    # Fit y = mx + c to the two datasets, with m being the same but c
    # different for the two datasets: m = 2, c_1 = 20, c_2 = -20.
    # However, the exact values do not matter for this test, as long
    # as they are separate enough that generating values for them
    # results in obviously-different values.
    #
    # The c_1 parameter is set via a link so that the thawed parameters
    # are in a different order for dataset 1 and 2.
    #
    # >>> ui.get_source(1).get_thawed_pars()
    # [<Parameter 'c1' of model 'pol1'>, <Parameter 'c0' of model 'scale'>]
    # >>> ui.get_source(2).get_thawed_pars()
    # [<Parameter 'c0' of model 'pol2'>, <Parameter 'c1' of model 'pol1'>]
    #
    pol1 = ui.create_model_component("polynom1d", "pol1")
    pol2 = ui.create_model_component("polynom1d", "pol2")
    scale = ui.create_model_component("scale1d", "scale")
    pol1.c0 = scale.c0
    pol1.c1.thaw()
    pol2.c1 = pol1.c1

    ui.set_source(1, pol1)
    ui.set_source(2, pol2)


@pytest.mark.parametrize("order", [(1, 2), (2, 1)])
def test_fit_multiple_datasets(order, clean_ui):
    """Check fit works with multiple datasets."""

    setup_multiple_datasets()

    ui.fit(*order)
    fres = ui.get_fit_results()
    assert fres.datasets == order
    assert fres.parnames == ("pol.c0", "pol.c1")

    # Note if the fit finds a different solution (the exact values are
    # not the point of this test).
    #
    assert fres.parvals == pytest.approx([21.756126722384383,
                                          1.248116121054563])


@pytest.mark.parametrize("order,parnames",
                         [((1, 2), ('pol1.c1', 'pol2.c0', 'scale.c0')),
                          ((2, 1), ('pol2.c0', 'pol1.c1', 'scale.c0'))
                          ])
def test_fit_multiple_datasets_linked(order, parnames, clean_ui):
    """Check fit works with multiple datasets."""

    setup_multiple_datasets_linked()

    ui.fit(*order)
    fres = ui.get_fit_results()
    assert fres.datasets == order
    assert fres.parnames == parnames

    # Note if the fit finds a different solution (the exact values are
    # not the point of this test).
    #
    vals = {"scale.c0": 22.002871808482844,
            "pol1.c1": 1.6842614178514361,
            "pol2.c0": -9.716190831799699}

    assert fres.parvals[0] == pytest.approx(vals[parnames[0]])
    assert fres.parvals[1] == pytest.approx(vals[parnames[1]])
    assert fres.parvals[2] == pytest.approx(vals[parnames[2]])


@pytest.mark.parametrize("order", [(1, 2), (2, 1)])
def test_covar_multiple_datasets(order, clean_ui):
    """Check covar works with multiple datasets."""

    setup_multiple_datasets()
    pol = ui.get_model_component("pol")

    ui.fit(*order)
    ui.covar(*order)
    cres = ui.get_covar_results()
    assert cres.datasets == order
    assert cres.parnames == ("pol.c0", "pol.c1")
    assert cres.parmaxes == pytest.approx([7.613218314229104,
                                           0.42324434548533496])

    # What happens if only a single parameter is asked for?
    #
    ui.covar(*order, pol.c1)
    cres_s = ui.get_covar_results()
    assert cres_s.datasets == order
    assert cres_s.parnames == ("pol.c1", )
    assert cres_s.parmaxes[0] == pytest.approx(0.42324434548533496)

    # The covariance matrix does not know about the parameter restriction.
    #
    assert cres_s.extra_output == pytest.approx(cres.extra_output)


@pytest.mark.parametrize("order,parnames",
                         [((1, 2), ('pol1.c1', 'pol2.c0', 'scale.c0')),
                          ((2, 1), ('pol2.c0', 'pol1.c1', 'scale.c0'))
                          ])
def test_covar_multiple_datasets_linked(order, parnames, clean_ui):
    """Check covar works with multiple datasets."""

    setup_multiple_datasets_linked()
    scale = ui.get_model_component("scale")

    ui.fit(*order)
    ui.covar(*order)
    cres = ui.get_covar_results()
    assert cres.datasets == order
    assert cres.parnames == parnames

    # Just check the maxes are in order
    maxes = {"scale.c0": 10.996504758773307,
             "pol1.c1": 0.7720581246708447,
             "pol2.c0": 17.6988390328967
             }

    assert cres.parmaxes[0] == pytest.approx(maxes[parnames[0]])
    assert cres.parmaxes[1] == pytest.approx(maxes[parnames[1]])
    assert cres.parmaxes[2] == pytest.approx(maxes[parnames[2]])

    # What happens if only a single parameter is asked for?
    #
    ui.covar(*order, scale.c0)
    cres_s = ui.get_covar_results()
    assert cres_s.datasets == order
    assert cres_s.parnames == ("scale.c0", )
    assert cres_s.parmaxes[0] == pytest.approx(maxes["scale.c0"])

    # The covariance matrix does not know about the parameter restriction.
    #
    assert cres_s.extra_output == pytest.approx(cres.extra_output)


@pytest.mark.parametrize("order", [(1, 2), (2, 1)])
def test_get_draws_multiple_datasets(order, clean_ui):
    """Check get_draws works with multiple datasets."""

    setup_multiple_datasets()

    ui.fit(*order)
    # ui.covar(*order)  as of 4.18.0 this is not needed

    # Ensure repeatability (unlike the _linked version the covariance
    # matrix is the same, so this is repeatable).
    #
    ui.set_rng(np.random.RandomState(3287))
    stats, accept, params = ui.get_draws(order[0], [order[1]], niter=100)

    # Although this is repeatable, just check with the mean value
    # rather than a per-iteration check. This is partly to be similar
    # to the _linked tests.
    #
    means = params.mean(axis=1)
    assert len(means) == 2
    assert means[0] == pytest.approx(21.219479992931166)
    assert means[1] == pytest.approx(1.3153808885552025)


@pytest.mark.parametrize("order,parnames",
                         [((1, 2), ('pol1.c1', 'pol2.c0', 'scale.c0')),
                          ((2, 1), ('pol2.c0', 'pol1.c1', 'scale.c0'))
                          ])
def test_get_draws_multiple_datasets_linked(order, parnames, clean_ui):
    """Check get_draws works with multiple datasets."""

    setup_multiple_datasets_linked()

    ui.fit(*order)
    # ui.covar(*order)  as of 4.18.0 this is not needed

    # Ensure repeatability. However, unlike the fit and covar results,
    # we do not get the same results.
    #
    ui.set_rng(np.random.RandomState(3287))
    stats, accept, params = ui.get_draws(order[0], [order[1]], niter=100)

    # Check the mean values of the parameters, as we assume they will
    # be "close" to the best-fit values as long as niter is large
    # enough.
    #
    # Convert to a dict to make the tests easier to read.
    means = dict(zip(parnames, params.mean(axis=1)))
    assert means['scale.c0'] == pytest.approx(22, abs=1)
    assert means['pol1.c1'] == pytest.approx(1.7, abs=0.1)
    assert means['pol2.c0'] == pytest.approx(-10, abs=1)


def test_get_draws_mismatched_covar_linked(clean_ui, caplog):
    """What happens when the covar call doesn't match the fit?"""

    setup_multiple_datasets_linked()
    scale = ui.get_model_component("scale")

    # Not only are the datasets swapped in the covar call, only one
    # component is asked for, and not all of them. This does not
    # actually matter since get_draws recalculates the covariance
    # matrix.
    #
    ui.fit(1, 2)
    ui.covar(2, 1, scale.c0)

    ui.set_rng(np.random.RandomState(3287))
    n1 = len(caplog.records)
    stats, accept, params = ui.get_draws(1, [2], niter=100)
    n2 = len(caplog.records)

    # Check the expected extra messages.
    #
    assert n2 == (n1 + 3)
    assert caplog.records[n1].getMessage() == "Using Priors:"
    assert caplog.records[n1 + 1].getMessage().startswith("pol1.c1: ")
    assert caplog.records[n1 + 2].getMessage().startswith("pol2.c0: ")

    # Unlike test_get_draws_multiple_datasets_linked we know the
    # ordering and can be much more precise in the tests.
    #
    means = params.mean(axis=1)
    assert means[0] == pytest.approx(1.73371932)
    assert means[1] == pytest.approx(-10.46427734)
    assert means[2] == pytest.approx(22.00661571)

    # There's at least one accept.
    assert any(accept)
