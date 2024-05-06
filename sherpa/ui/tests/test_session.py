#
#  Copyright (C) 2016, 2017, 2019 - 2024
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
There is a copy of this test in sherpa/ui/astro/tests/ that needs
to be kept up to date.
"""

from io import StringIO
import logging
import pickle
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_array_equal

import pytest

from sherpa.models import Model, parameter
import sherpa.models.basic
from sherpa import estmethods as est
from sherpa import optmethods as opt
from sherpa import stats
from sherpa.ui.utils import Session, ModelWrapper
from sherpa.utils import poisson_noise
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, \
    IdentifierErr, IOErr, SessionErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_xspec

TEST = [1, 2, 3]
TEST2 = [4, 5, 6]


# bug #303
#
# Note: this test could be run even if plotting is not available,
# since in that case the initial get_data_plot_prefs() call will
# return an empty dictionary rather than a filled one. However,
# this seems to cause problems so keep the fixture.
#
def test_set_log():
    session = Session()
    assert not session.get_data_plot_prefs()['xlog']
    assert not session.get_data_plot_prefs()['ylog']
    session.set_xlog()
    assert session.get_data_plot_prefs()['xlog']
    session.set_ylog()
    assert session.get_data_plot_prefs()['ylog']
    session.set_xlinear()
    assert not session.get_data_plot_prefs()['xlog']
    session.set_ylinear()
    assert not session.get_data_plot_prefs()['ylog']


def test_set_log_clean():
    """Cleaning the session resets the plot preferences."""

    # relies on the tests in test_set_log
    session = Session()
    session.set_xlog()
    session.set_ylog()

    session.clean()
    assert not session.get_data_plot_prefs()['xlog']
    assert not session.get_data_plot_prefs()['ylog']


def test_set_log_does_not_change_other_sessions():
    """The plot preferences in different sessions are distinct.
    """

    session1 = Session()
    session2 = Session()
    session1.set_xlog()
    session2.set_ylog()

    assert session1.get_data_plot_prefs()['xlog']
    assert not session1.get_data_plot_prefs()['ylog']
    assert not session2.get_data_plot_prefs()['xlog']
    assert session2.get_data_plot_prefs()['ylog']


# bug #262
def test_list_ids():
    session = Session()
    session.load_arrays(1, TEST, TEST)
    session.load_arrays("1", TEST, TEST2)

    # order of 1 and "1" is not determined
    assert {1, "1"} == set(session.list_data_ids())
    assert_array_equal(TEST2, session.get_data('1').get_dep())
    assert_array_equal(TEST, session.get_data(1).get_dep())


# bug #297
def test_save_restore(tmp_path):
    outfile = tmp_path / "sherpa.save"
    session = Session()
    session.load_arrays(1, TEST, TEST2)
    session.save(str(outfile), clobber=True)
    session.clean()
    assert set() == set(session.list_data_ids())

    session.restore(str(outfile))
    assert {1, } == set(session.list_data_ids())
    assert_array_equal(TEST, session.get_data(1).get_indep()[0])
    assert_array_equal(TEST2, session.get_data(1).get_dep())


def test_save_clobber_check(tmp_path):
    """Check we handle clobber case"""

    tmp = tmp_path / "x.sherpa"
    tmp.write_text("not-empty")

    s = Session()
    with pytest.raises(IOErr,
                       match="' exists and clobber is not set"):
        s.save(str(tmp), clobber=False)


def test_restore_not_a_session(tmp_path):
    """Check we error out when it's not a Session object"""

    # Argh tmppath and pickle don't play well together
    tmp = tmp_path / "not-sherpa.state"
    with open(tmp, 'wb') as ofh:
        pickle.dump({"x": 23}, ofh)

    s = Session()
    with pytest.raises(ArgumentErr,
                       match="' does not contain a saved Sherpa session"):
        s.restore(str(tmp))


def test_models_models():
    """There are no models available by default"""

    s = Session()
    assert s.list_models() == []


def test_models_all():
    """We can get a list of all models"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    models = s.list_models()
    assert type(models) == list
    assert len(models) > 1
    assert all([type(m) == str for m in models])


# XSPEC is only needed for the xspec test but easiest to just
# mark it for all cases
#
@requires_xspec
@pytest.mark.parametrize("mtype", ["1d", "2d", "xspec"])
def test_models_filtered(mtype):
    """We can get a list of a subset of models

    Check the fix for issue #749
    """

    import sherpa.astro.xspec

    s = Session()
    s._add_model_types(sherpa.models.basic)
    s._add_model_types(sherpa.astro.xspec)

    models = s.list_models(mtype)
    assert type(models) == list
    assert len(models) > 1
    assert all([type(m) == str for m in models])


@requires_xspec
def test_models_all_xspec():
    """Check all models are returned"""

    import sherpa.astro.xspec

    s = Session()
    s._add_model_types(sherpa.models.basic)
    s._add_model_types(sherpa.astro.xspec)

    # Actually, it's not guaranteed that these models are disjoint
    # but they are at the moment.
    #
    mall = s.list_models('all')
    m1d = s.list_models('1d')
    m2d = s.list_models('2d')
    mxs = s.list_models('xspec')

    nall = len(mall)
    assert nall == len(m1d) + len(m2d) + len(mxs)

    mall = set(mall)
    assert nall == len(mall)

    mcomb = set(m1d).union(set(m2d)).union(set(mxs))
    assert nall == len(mcomb)


def test_get_functions():
    """Limited check of get_functions"""

    s = Session()
    fns = s.get_functions()
    assert type(fns) == list
    assert len(fns) > 1
    assert all([type(f) == str for f in fns])
    assert 'load_data' in fns


def test_list_functions():
    """Limited check of list_functions"""

    s = Session()
    store = StringIO()
    s.list_functions(outfile=store)
    txt = store.getvalue()
    fns = txt.split('\n')
    assert len(fns) > 1
    assert 'calc_stat' in fns


def test_paramprompt_function():
    """Does paramprompt toggle the state setting?"""

    s = Session()
    assert not s._paramprompt

    s.paramprompt(True)
    assert s._paramprompt

    s.paramprompt(False)
    assert not s._paramprompt


def test_paramprompt():
    """Does paramprompt work?"""

    s = Session()
    assert s.list_model_ids() == []
    assert s.list_model_components() == []

    # Add in some models
    s._add_model_types(sherpa.models.basic)

    models = s.list_models()
    assert models != []
    assert 'const1d' in models

    s.create_model_component('const1d', 'm1')
    assert s.list_model_ids() == []

    expected = ['m1']
    assert s.list_model_components() == expected

    # Now there are multiple components do not rely on any ordering
    # provided by list_model_components()
    #
    s.paramprompt(True)

    # paramprompt doesn't affect create_model_component (as shown
    # by the fact that we don't have to mock stdin)
    #
    s.create_model_component('const1d', 'm2')
    assert s.list_model_ids() == []
    expected = set(expected)
    expected.add('m2')
    assert set(s.list_model_components()) == expected

    # it does affect set_model
    #
    # pytest errors out if you try to read from stdin in a test, so
    # try to work around this here using
    # https://stackoverflow.com/questions/13238566/python-equivalent-of-input-using-sys-stdin
    # An alternative would be just to check that the error is
    # raised, and the text, but then we don't get to test the
    # behaviour of the paramprompt code.
    #
    with patch("sys.stdin", StringIO("2.1")):
        s.set_model('const1d.mx')

    assert s.list_model_ids() == [1]
    expected.add('mx')
    assert set(s.list_model_components()) == expected

    mx = s.get_model_component('mx')
    assert mx.c0.val == pytest.approx(2.1)

    with patch("sys.stdin", StringIO("2.1e-3 , 2.0e-3, 1.2e-2")):
        s.set_model('x', 'const1d.my')

    assert set(s.list_model_ids()) == set([1, 'x'])
    expected.add('my')
    assert set(s.list_model_components()) == expected

    my = s.get_model_component('my')
    assert my.c0.val == pytest.approx(2.1e-3)
    assert my.c0.min == pytest.approx(2.0e-3)
    assert my.c0.max == pytest.approx(1.2e-2)

    with patch("sys.stdin", StringIO("2.1e-3 ,, 1.2e-2")):
        s.set_model(2, 'const1d.mz')

    assert set(s.list_model_ids()) == set([1, 2, 'x'])
    expected.add('mz')
    assert set(s.list_model_components()) == expected

    mz = s.get_model_component('mz')
    assert mz.c0.val == pytest.approx(2.1e-3)
    assert mz.c0.min == pytest.approx(- parameter.hugeval)
    assert mz.c0.max == pytest.approx(1.2e-2)

    with patch("sys.stdin", StringIO("-12.1,-12.1")):
        s.set_model('const1d.f1')

    assert set(s.list_model_ids()) == set([1, 2, 'x'])
    expected.add('f1')
    assert set(s.list_model_components()) == expected

    f1 = s.get_model_component('f1')
    assert f1.c0.val == pytest.approx(-12.1)
    assert f1.c0.min == pytest.approx(-12.1)
    assert f1.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(" ,-12.1")):
        s.set_model('const1d.f2')

    # stop checking list_model_ids
    expected.add('f2')
    assert set(s.list_model_components()) == expected

    f2 = s.get_model_component('f2')
    assert f2.c0.val == pytest.approx(1.0)
    assert f2.c0.min == pytest.approx(-12.1)
    assert f2.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(" ,")):
        s.set_model('const1d.f3')

    f3 = s.get_model_component('f3')
    assert f3.c0.val == pytest.approx(1.0)
    assert f3.c0.min == pytest.approx(- parameter.hugeval)
    assert f3.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(" ,, ")):
        s.set_model('const1d.f4')

    f4 = s.get_model_component('f4')
    assert f4.c0.val == pytest.approx(1.0)
    assert f4.c0.min == pytest.approx(- parameter.hugeval)
    assert f4.c0.max == pytest.approx(parameter.hugeval)


def test_list_model_ids():
    """Does list_model_ids work?"""

    # This issue was found when developing the paramprompt test.
    s = Session()

    s._add_model_types(sherpa.models.basic)

    assert s.list_model_ids() == []

    s.set_model('ma', 'const1d.mdla')
    assert s.list_model_ids() == ['ma']

    s.set_model('mb', 'const1d.mdla')
    assert set(s.list_model_ids()) == set(['ma', 'mb'])

    s.set_model('const1d.mdla')
    assert set(s.list_model_ids()) == set([1, 'ma', 'mb'])


def test_list_methods():
    """we have some methods"""

    s = Session()
    methods = s.list_methods()
    assert type(methods) == list
    assert len(methods) > 1  # do not check exact number
    assert 'levmar' in methods


def test_list_iter_methods():
    """we have some iteration methods"""

    s = Session()
    methods = s.list_iter_methods()
    assert type(methods) == list
    assert len(methods) > 1  # do not check exact number
    assert 'none' in methods
    assert 'sigmarej' in methods


def test_get_method_default():
    """get_method returns default instance (LevMar)"""

    s = Session()
    method = s.get_method()
    assert isinstance(method, opt.LevMar)


@pytest.mark.parametrize("name,req",
                         [("GridSearch", opt.GridSearch),
                          ("levmar", opt.LevMar),
                          ("MONCAR", opt.MonCar),
                          ("neldermead", opt.NelderMead),
                          ("simplex", opt.NelderMead)])
def test_get_method_named(name, req):
    """get_method returns requested instance"""

    s = Session()
    method = s.get_method(name)
    assert isinstance(method, req)


def test_get_method_invalid():
    """Errors out with invalid argument"""

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'name' must be a string$"):
        s.get_method(opt.MonCar)


@pytest.mark.parametrize("name,req",
                         [("moncar", opt.MonCar),
                          (opt.GridSearch(), opt.GridSearch)])
def test_set_method(name, req):
    """We can set the method"""

    s = Session()
    s.set_method(name)
    ans = s.get_method()
    assert isinstance(ans, req)


def test_set_method_invalid():
    """Errors out with invalid argument"""

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'meth' must be a method name or object$"):
        s.set_method(sherpa.models.basic.Const1D)


def test_get_stat_default():
    """get_stat returns default instance (Chi2Gehrels)"""

    s = Session()
    stat = s.get_stat()
    assert isinstance(stat, stats.Chi2Gehrels)


@pytest.mark.parametrize("name,req",
                         [("chi2gehrels", stats.Chi2Gehrels),
                          ("Chi2DataVar", stats.Chi2DataVar),
                          ("leastsq", stats.LeastSq),
                          ("Wstat", stats.WStat)])
def test_get_stat_named(name, req):
    """get_stat returns requested instance"""

    s = Session()
    stat = s.get_stat(name)
    assert isinstance(stat, req)


def test_get_stat_invalid():
    """Errors out with invalid argument"""

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'name' must be a string$"):
        s.get_stat(stats.Cash)


@pytest.mark.parametrize("name,req",
                         [("chi2modvar", stats.Chi2ModVar),
                          (stats.CStat(), stats.CStat)])
def test_set_stat(name, req):
    """We can set the statistic"""

    s = Session()
    s.set_stat(name)
    ans = s.get_stat()
    assert isinstance(ans, req)


def test_set_stat_invalid():
    """Errors out with invalid argument"""

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'stat' must be a statistic name or object$"):
        s.set_stat(sherpa.models.basic.Const1D)


def test_get_default_id():
    """Does the default id react correctly?"""
    s = Session()
    assert s.get_default_id() == 1

    s.set_default_id('alpha')
    assert s.get_default_id() == 'alpha'


@pytest.mark.parametrize("name,req",
                         [('covar', est.Covariance),
                          ('conf', est.Confidence),
                          ('proj', est.Projection)])
def test_get_error_estimator(name, req):
    """Can we get the error estimator?"""

    s = Session()
    func = getattr(s, f'get_{name}')
    ans = func()
    assert isinstance(ans, req)


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_get_error_opt(name):
    """Can we get an option for an error estimator?

    Fortunately they all have a sigma option.
    """

    s = Session()
    func = getattr(s, f'get_{name}_opt')
    ans = func('sigma')
    assert ans == pytest.approx(1.0)


@pytest.mark.parametrize("name,fullname",
                         [('covar', 'covariance'),
                          ('conf', 'confidence'),
                          ('proj', 'projection')])
def test_get_error_opt_invalid(name, fullname):
    """We can not ask for an error option that does not exist"""

    s = Session()
    func = getattr(s, f'get_{name}_opt')

    emsg = f"^'the-real-sigma' is not a valid option for method {fullname}$"
    with pytest.raises(ArgumentErr, match=emsg):
        func('the-real-sigma')


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_set_error_opt(name):
    """Can we set an option for an error estimator?

    Fortunately they all have a sigma option.
    """

    s = Session()
    gfunc = getattr(s, f'get_{name}_opt')
    sfunc = getattr(s, f'set_{name}_opt')

    sfunc('sigma', 1.6)
    ans = gfunc('sigma')
    assert ans == pytest.approx(1.6)


@pytest.mark.parametrize("name,fullname",
                         [('covar', 'covariance'),
                          ('conf', 'confidence'),
                          ('proj', 'projection')])
def test_set_error_opt_invalid(name, fullname):
    """We can not set an error option that does not exist"""

    s = Session()
    func = getattr(s, f'set_{name}_opt')

    emsg = f"^'the-real-sigma' is not a valid option for method {fullname}$"
    with pytest.raises(ArgumentErr, match=emsg):
        func('the-real-sigma', 2.4)


@pytest.mark.parametrize("name,fullname",
                         [('covar', 'covariance'),
                          ('conf', 'confidence'),
                          ('proj', 'projection')])
def test_error_estimate_not_set(name, fullname):
    """Error out if the error estimate is not set"""

    s = Session()
    func = getattr(s, f'get_{name}_results')

    with pytest.raises(SessionErr,
                       match=f"^{fullname} has not been performed"):
        func()


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_error_estimate_not_run(name):
    """Can not rn the error estimate if data is not set up"""

    s = Session()
    s.set_default_id('bob')
    func = getattr(s, name)

    with pytest.raises(IdentifierErr,
                       match="^data set bob has not been set"):
        func()


def test_set_source_invalid():

    s = Session()
    with pytest.raises(ArgumentErr,
                       match="^invalid model expression: name 'made_up' is not defined$"):
        s.set_source('2 * made_up.foo')


def test_paramprompt_default():
    """The default setting is False."""

    s = Session()
    assert not s._paramprompt


@pytest.mark.parametrize("flag", [True, "yes", "yeS", 1, "1", "on", "ON"])
def test_paramprompt_set_true(flag):
    """The default setting is False."""

    s = Session()
    s.paramprompt(flag)
    assert s._paramprompt


@pytest.mark.parametrize("flag", [False, "no", "No", 0, "0", "off", "OFF"])
def test_paramprompt_set_false(flag):
    """The default setting is False."""

    s = Session()
    s._paramprompt = True
    s.paramprompt(flag)
    assert not s._paramprompt


def test_paramprompt_single_parameter_works(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO("223")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 0
    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(223)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max > 3e38

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_min_works(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO(",-100")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 0
    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min == pytest.approx(-100)
    assert mdl.c0.max > 3e38

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_max_works(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO(",,100")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 0
    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max == pytest.approx(100)

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_combo_works(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO("-2,-10,10")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 0
    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(-2)
    assert mdl.c0.min == pytest.approx(-10)
    assert mdl.c0.max == pytest.approx(10)

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_check_invalid(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO("typo\n-200")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "Please provide a float value; could not convert string to float: 'typo'"

    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(-200)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max > 3e38

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_check_invalid_min(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO(",typo\n,-200")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "Please provide a float value; could not convert string to float: 'typo'"

    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(1)
    assert mdl.c0.min == pytest.approx(-200)
    assert mdl.c0.max > 3e38

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_check_invalid_max(caplog):

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO("-2,,typo\n-200,,-2")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "Please provide a float value; could not convert string to float: 'typo'"

    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(-200)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max == pytest.approx(-2)

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_check_invalid_max_out_of_bound(caplog):
    """Note this creates two warnings"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO(",,typo\n,,-200")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 2
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "Please provide a float value; could not convert string to float: 'typo'"

    lname, lvl, msg = caplog.record_tuples[1]
    assert lname == "sherpa.models.parameter"
    assert lvl == logging.WARN
    assert msg == "parameter bob.c0 greater than new maximum; bob.c0 reset to -200"

    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(-200)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max == pytest.approx(-200)

    # remove the bob symbol from the global table
    s.clean()


def test_paramprompt_single_parameter_check_too_many_commas(caplog):
    """Check we tell users there was a problem"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO(",,,,\n12")):
            s.set_source("scale1d.bob")

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.INFO
    assert msg == "Error: Please provide a comma-separated list of floats; e.g. val,min,max"

    mdl = s.get_model_component('bob')
    assert mdl.c0.val == pytest.approx(12)
    assert mdl.c0.min < -3e38
    assert mdl.c0.max > 3e38

    # remove the bob symbol from the global table
    s.clean()


def test_add_user_pars_modelname_not_a_string():

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'model name' must be a string$"):
        s.add_user_pars(23, ['x'])


def test_add_user_pars_modelname_not_a_model1():

    s = Session()
    with pytest.raises(ArgumentTypeErr,
                       match="^'not a model' must be a user model$"):
        s.add_user_pars('not a model', ['x'])


def test_add_user_pars_modelname_not_a_model2():
    """Use an actual model, but not a user model"""

    s = Session()
    s._add_model_types(sherpa.models.basic)
    s.create_model_component('scale1d', 'foo')
    with pytest.raises(ArgumentTypeErr,
                       match="^'foo' must be a user model$"):
        s.add_user_pars('foo', ['x'])

    # remove the foo symbol from the global table
    s.clean()


def test_paramprompt_multi_parameter(caplog):
    """Check that paramprompt works with multiple parameters"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    cpt1 = s.create_model_component('const1d', 'bob')
    cpt2 = s.create_model_component('gauss1d', 'fred')

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with SherpaVerbosity('INFO'):
        with patch("sys.stdin", StringIO("\n5,1,5\n\n-5, -5, 0")):
            s.set_source(bob + fred)

    assert len(caplog.records) == 0

    assert cpt1.c0.val == pytest.approx(1)
    assert cpt1.c0.min < -3e38
    assert cpt1.c0.max > 3e38

    assert cpt2.fwhm.val == pytest.approx(5)
    assert cpt2.fwhm.min == pytest.approx(1)
    assert cpt2.fwhm.max == pytest.approx(5)

    assert cpt2.pos.val == pytest.approx(0)
    assert cpt2.pos.min < -3e38
    assert cpt2.pos.max > 3e38

    assert cpt2.ampl.val == pytest.approx(-5)
    assert cpt2.ampl.min == pytest.approx(-5)
    assert cpt2.ampl.max == pytest.approx(0)

    # remove the bob and fred symbols from the global table
    s.clean()


def test_paramprompt_eof(caplog):
    """What happens when we end early?"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    cpt1 = s.create_model_component('const1d', 'bob')
    cpt2 = s.create_model_component('gauss1d', 'fred')

    s.paramprompt(True)
    assert len(caplog.records) == 0

    with pytest.raises(EOFError):
        with SherpaVerbosity('INFO'):
            with patch("sys.stdin", StringIO("\n5,1,5\n2\n")):
                s.set_source(bob + fred)

    assert len(caplog.records) == 0

    assert cpt1.c0.val == pytest.approx(1)
    assert cpt1.c0.min < -3e38
    assert cpt1.c0.max > 3e38

    assert cpt2.fwhm.val == pytest.approx(5)
    assert cpt2.fwhm.min == pytest.approx(1)
    assert cpt2.fwhm.max == pytest.approx(5)

    assert cpt2.pos.val == pytest.approx(2)
    assert cpt2.pos.min < -3e38
    assert cpt2.pos.max > 3e38

    assert cpt2.ampl.val == pytest.approx(1)
    assert cpt2.ampl.min < -3e38
    assert cpt2.ampl.max > 3e38

    # remove the bob and fred symbols from the global table
    s.clean()


def test_delete_model_component_invalid_argument():
    """We could allow the component to be deleted this way"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    tst = s.create_model_component('gauss1d', 'gmdl')

    with pytest.raises(ArgumentTypeErr,
                       match="'name' must be a string"):
        s.delete_model_component(tst)


def test_delete_model_component_not_a_component():
    """Check correct error message for non-existent model"""

    s = Session()
    with pytest.raises(IdentifierErr,
                       match="model component 'tst' does not exist"):
        s.delete_model_component('tst')


def test_delete_model_component_warning(caplog):
    """Check we get a warning (which ends up being issue #16)"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.set_source('const1d.mdl + gauss1d.mdl2')
    assert s.list_model_components() == ['mdl', 'mdl2']

    assert len(caplog.records) == 0
    s.delete_model_component('mdl2')

    assert len(caplog.records) == 1
    lname, lvl, msg = caplog.record_tuples[0]
    assert lname == "sherpa.ui.utils"
    assert lvl == logging.WARNING
    assert msg == "the model component 'gauss1d.mdl2' is found in model 1 and cannot be deleted"

    assert s.list_model_components() == ['mdl', 'mdl2']

    # remove the mdl1 and mdl2 symbols from the global table
    s.clean()


def test_issue_16():
    """Check one of the examples from #16"""

    s = Session()
    s._add_model_types(sherpa.models.basic)

    s.dataspace1d(0.01, 11, 0.01, id=1)
    s.dataspace1d(2, 5, 0.1, id="tst")
    s.set_source(1, 'powlaw1d.pl1')
    s.set_source('tst', 'powlaw1d.pltst')

    assert s.list_data_ids() == [1, 'tst']
    assert s.list_model_ids() == [1, 'tst']
    assert s.list_model_components() == ['pl1', 'pltst']

    s.delete_model(id='tst')
    s.delete_model_component("pltst")
    s.delete_data(id='tst')

    assert s.list_data_ids() == [1]
    assert s.list_model_ids() == [1]
    assert s.list_model_components() == ['pl1']

    s.delete_model(id=1)
    s.delete_model_component("pl1")
    s.delete_data(id=1)

    assert s.list_data_ids() == []
    assert s.list_model_ids() == []
    assert s.list_model_components() == []


@pytest.mark.parametrize("value", ["data", "model", "source", "fit", "resid", "ratio",
                                   "delchi", "chisqr", "psf", "kernel", "compsource",
                                   "compmodel", "source_component", "model_component",])
def test_set_default_id_check_invalid(value):
    """Check we error out with an invalid id.

    There are other checks of this logic, but not as simple
    as this. The test in
    sherpa/astro/ui/tests/test_astro_session.py::test_id_checks_astro_session
    adds a test for astro-specific keywords.

    """

    s = Session()

    emsg = f"^identifier '{value}' is a reserved word$"
    with pytest.raises(IdentifierErr, match=emsg):
        s.set_default_id(value)


class ModelWithDoc(Model):
    """This has a doc string

    This line is not included in the model-wrapped doc string.
    """

    _hidden = False


class ModelWithNoDoc(Model):
    # This does not have a doc string

    pass


def test_modelwrapper_init():
    """Basic ModelWrapper check"""

    s = Session()
    wrap = ModelWrapper(s, ModelWithDoc)
    assert repr(wrap) == "<ModelWithDoc model type>"


def test_modelwrapper_checks_session():
    """Check we error out"""

    with pytest.raises(ValueError,
                       match="session=.* is not a Session instance"):
        ModelWrapper(ModelWithDoc, None)


def test_modelwrapper_checks_model():
    """Check we error out"""

    s = Session()
    with pytest.raises(ValueError,
                       match="modeltype=Gauss1D is not a Model class"):
        ModelWrapper(s, "Gauss1D")


def test_modelwrapper_str_with_doc():
    """Basic ModelWrapper check"""

    s = Session()
    wrap = ModelWrapper(s, ModelWithDoc)
    assert str(wrap) == "This has a doc string\n\n    This line is not included in the model-wrapped doc string.\n    "


def test_modelwrapper_str_no_doc():
    """Basic ModelWrapper check"""

    s = Session()
    wrap = ModelWrapper(s, ModelWithNoDoc)
    assert str(wrap) == "<ModelWithNoDoc model type>"


def test_modelwrapper_what_is_the_docstring_has_doc():
    """What do we want help(wrappedmodel) to return"""

    # Check that the class and wrapped docstrings are different.
    #
    assert ModelWrapper.__doc__ is not None

    s = Session()
    wrap = ModelWrapper(s, ModelWithDoc)
    assert wrap.__doc__ != ModelWrapper.__doc__

    # Check the wrapped docstring references the model and description.
    #
    assert wrap.__doc__.startswith("Create a modelwithdoc model instance.\n\n    This has a doc string\n\n    Instances can")


def test_modelwrapper_what_is_the_docstring_no_doc():
    """What do we want help(wrappedmodel) to return"""

    # Check that the class and wrapped docstrings are different.
    #
    assert ModelWrapper.__doc__ is not None

    s = Session()
    wrap = ModelWrapper(s, ModelWithNoDoc)
    assert wrap.__doc__ != ModelWrapper.__doc__

    # Check the wrapped docstring references the model but no description.
    assert wrap.__doc__.startswith("Create a modelwithnodoc model instance.\n\n    Instances can")


@pytest.mark.parametrize("attr", ["_hidden", "_foo_bar"])
def test_modelwrapper_getattr_no_hidden(attr):
    """How does attribute access work?

    We can't access attributes of the original model or unknown
    attributes if they start with _.

    """

    s = Session()
    wrap = ModelWrapper(s, ModelWithDoc)

    with pytest.raises(AttributeError,
                       match=f"^'ModelWrapper' object has no attribute '{attr}'$"):
        getattr(wrap, attr)


def test_modelwrapper_getattr_instance():
    """The access works for the model instance"""

    s = Session()
    wrap = ModelWrapper(s, ModelWithDoc)
    mdl = wrap("bob")
    assert not mdl._hidden


def test_modelwrapper_getattr_instance_unknown():
    """The access fails for the model instance when field is unknown"""

    s = Session()
    wrap = ModelWrapper(s, ModelWithNoDoc)
    mdl = wrap("bob")
    with pytest.raises(AttributeError,
                       match="^'ModelWithNoDoc' object has no attribute '_hidden'$"):
        got = mdl._hidden


def test_modelwrapper_direct():
    """Compare behavior with test_modelwrapper_call"""

    s = Session()

    wrap = ModelWrapper(s, ModelWithNoDoc)
    assert s.list_model_components() == []

    mdl = wrap("foo")
    assert isinstance(mdl, ModelWithNoDoc)
    assert mdl.name == "modelwithnodoc.foo"

    assert s.list_model_components() == ["foo"]


def test_modelwrapper_call():
    """The reason for the class"""

    s = Session()

    wrap = ModelWrapper(s, ModelWithNoDoc)
    assert s.list_model_components() == []

    mdl = wrap.foo
    assert isinstance(mdl, ModelWithNoDoc)
    assert mdl.name == "modelwithnodoc.foo"

    assert s.list_model_components() == ["foo"]


def test_show_data_explicit_id():
    """Check we can call show_data with an id"""

    s = Session()
    s.load_arrays(1, [1, 2], [10, 20])
    s.load_arrays("bob", [12, 14, 20], [3, 1, 2])

    out = StringIO()
    s.show_data("bob", outfile=out)

    # not a full check of the output
    assert out.getvalue().startswith("Data Set: bob\nname      = \n")


def test_show_kernel_explicit_id(tmp_path):
    """Check we can call show_kernel with an id"""

    s = Session()

    s.load_arrays("bob", [1, 2, 3, 4, 5], [10, 5, 3, 2, 1])

    tmp = tmp_path / "tmp.dat"
    tmp.write_text("! X Y\n1 10\n2 5\n 3 0\n4 1\n")

    s.load_psf("foo", str(tmp), comment="!")
    s.set_psf("bob", "foo")

    out = StringIO()
    s.show_kernel("bob", outfile=out)

    # Since we do not have many other tests of this, do a more-complete
    # check than some of the other routines here.
    #
    toks = out.getvalue().split("\n")
    assert toks[0] == "PSF Kernel: bob"
    assert toks[1] == "psfmodel.foo"
    assert " Param " in toks[2]
    assert " ----- " in toks[3]

    # Can the file path contain spaces?
    words = toks[4].split()
    assert len(words) >= 3
    assert words[0] == "foo.kernel"
    assert words[1] == "frozen"

    def check(lstr, name, val, minval, maxval):
        words = lstr.split()
        assert len(words) == 5
        assert words[0] == name
        assert words[1] == "frozen"
        assert words[2] == val
        assert words[3] == minval
        assert words[4] == maxval

    check(toks[5], "foo.size", "4", "4", "4")
    check(toks[6], "foo.center", "2", "2", "2")
    check(toks[7], "foo.radial", "0", "0", "1")
    check(toks[8], "foo.norm", "1", "0", "1")

    assert toks[9] == ""
    assert toks[10] == ""
    assert toks[11] == ""
    assert len(toks) == 12


def test_show_psf_explicit_id(tmp_path):
    """Check we can call show_psf with an id"""

    s = Session()

    s.load_arrays("bob", [1, 2, 3, 4, 5], [10, 5, 3, 2, 1])

    tmp = tmp_path / "tmp.dat"
    tmp.write_text("! X Y\n1 10\n2 5\n 3 0\n4 1\n")

    s.load_psf("foo", str(tmp), comment="!")
    s.set_psf("bob", "foo")

    out = StringIO()
    s.show_psf("bob", outfile=out)

    print(out.getvalue())

    # Since we do not have many other tests of this, do a more-complete
    # check than some of the other routines here.
    #
    toks = out.getvalue().split("\n")
    assert toks[0] == "PSF Model: bob"
    assert toks[1].startswith("name      = ")
    assert toks[2] == "x         = Float64[4]"
    assert toks[3] == "y         = Float64[4]"
    assert toks[4] == "staterror = None"
    assert toks[5] == "syserror  = None"

    assert toks[6] == ""
    assert toks[7] == ""
    assert toks[8] == ""
    assert len(toks) == 9


def test_show_kernel_multiple(tmp_path):
    """Check we can call show_kernel with multiple datasets, not all with a PSF"""

    s = Session()

    s.load_arrays("bob", [1, 2, 3, 4, 5], [10, 5, 3, 2, 1])
    s.load_arrays(1, [1, 2, 3], [4, 5, 6])

    tmp = tmp_path / "tmp.dat"
    tmp.write_text("! X Y\n1 10\n2 5\n 3 0\n4 1\n")

    s.load_psf("foo", str(tmp), comment="!")
    s.set_psf("bob", "foo")

    out = StringIO()
    s.show_kernel("bob", outfile=out)

    print(out.getvalue())

    # Since we do not have many other tests of this, do a more-complete
    # check than some of the other routines here.
    #
    toks = out.getvalue().split("\n")
    assert toks[0] == "PSF Kernel: bob"
    assert toks[1] == "psfmodel.foo"
    assert toks[2][2:].startswith(" Param ")
    assert toks[3][2:].startswith(" ----- ")
    assert toks[4].startswith("   foo.kernel   frozen ")
    assert toks[5].startswith("   foo.size     frozen ")
    assert toks[6].startswith("   foo.center   frozen ")
    assert toks[7].startswith("   foo.radial   frozen ")
    assert toks[8].startswith("   foo.norm     frozen ")
    assert toks[9] == ""
    assert toks[10] == ""
    assert toks[11] == ""
    assert len(toks) == 12


@pytest.mark.parametrize("idval", [None, 1])
def test_load_filter_simple(idval, tmp_path):
    """Although there is a version in astro we only check the non-astro
    version as the behavior is different (e.g. the astro version needs
    I/O support).

    """

    s = Session()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2])

    infile = tmp_path / "filter.dat"
    infile.write_text("0\n0\n1\n")

    if idval is None:
        s.load_filter(str(infile), ncols=1, ignore=True)
    else:
        s.load_filter(idval, str(infile), ncols=1, ignore=True)

    assert s.get_data().mask == pytest.approx([True, True, False])


@pytest.mark.parametrize("idval", [None, 1])
def test_save_filter_simple(idval, tmp_path):
    """Although there is a version in astro we only check the non-astro
    version as the behavior is different (e.g. the astro version needs
    I/O support).

    """

    s = Session()
    s.load_arrays(idval, [1, 2, 3], [12, 15, 2])
    s.get_data(idval).mask = [True, False, True]

    outfile = tmp_path / "filter.dat"
    if idval is None:
        s.save_filter(str(outfile), comment="P ", linebreak="+")
    else:
        s.save_filter(idval, str(outfile), comment="P ", linebreak="+")

    assert outfile.read_text() == "P X FILTER+1 1+2 0+3 1+"


def setup_linked_data(s):
    """A simple 'random' dataset."""

    # We don't care about randomness as much as repeatability here
    rng = np.random.RandomState(123)

    base = s.create_model_component("gauss1d", "base")
    base.pos = 100
    base.ampl = 50
    base.fwhm = 20

    x = np.arange(70, 130, 5)
    y_base = base(x)
    y = poisson_noise(y_base, rng=rng)

    s.load_arrays(1, x, y)


LPAR_STAT = -949.8639549430894
LPAR_FWHM = 19.910836513919275
LPAR_SIGMA = 19.910836513919275 / 2.3548200450309493
LPAR_AMPL = 47.4433673308311


def test_linked_parameter_explicit():
    """See sherpa/tests/test_fit_unit.py::test_linked_parameter_explicit

    We do not include test_fit_unit.py::test_linked_parameter_base as
    we expect this to work correctly thanks to the other tests.
    """

    s = Session()
    s._add_model_types(sherpa.models.basic)

    setup_linked_data(s)
    s.set_stat("cash")

    m2 = s.create_model_component("gauss1d", "fit2")
    s2 = s.create_model_component("scale1d", "sigma2")

    m2.pos = 100
    m2.pos.freeze()

    convert = 2 * np.sqrt(2 * np.log(2))
    m2.fwhm = convert * sigma2.c0

    # Adjust c0 to better match the FWHM=10 value (just so we can
    # check we end up in the same location as the
    # sherpa/tests/test_fit_unit.py tests). It's important to use the
    # .val field here so we are not accidentally setting up another
    # link.
    #
    sigma2.c0.val = 4.2466090014400955

    s.set_model(m2 + 0 * s2)

    assert s.get_num_par() == 4
    assert s.get_num_par_thawed() == 2

    s.fit()
    res = s.get_fit_results()

    assert s.calc_stat() == pytest.approx(LPAR_STAT)
    assert m2.fwhm.val == pytest.approx(LPAR_FWHM)
    assert s2.c0.val == pytest.approx(LPAR_SIGMA)
    assert m2.ampl.val == pytest.approx(LPAR_AMPL)

    assert res.statval == pytest.approx(LPAR_STAT)
    assert res.numpoints == 12
    assert res.dof == 10
    assert res.parnames == ('fit2.ampl', 'sigma2.c0')
    assert res.parvals[0] == pytest.approx(LPAR_AMPL)
    assert res.parvals[1] == pytest.approx(LPAR_SIGMA)


def test_linked_parameter_implicit():
    """See test_linked_parameter_explicit

    What happens if we do not link the parameter to the model expression.
    """

    s = Session()
    s._add_model_types(sherpa.models.basic)

    setup_linked_data(s)
    s.set_stat("cash")

    m3 = s.create_model_component("gauss1d", "fit3")
    s3 = s.create_model_component("scale1d", "sigma3")

    m3.pos = 100
    m3.pos.freeze()

    convert = 2 * np.sqrt(2 * np.log(2))
    m3.fwhm = convert * sigma3.c0

    # Adjust c0 to better match the FWHM=10 value (just so we can
    # check we end up in the same location as the
    # sherpa/tests/test_fit_unit.py tests). It's important to use the
    # .val field here so we are not accidentally setting up another
    # link.
    #
    sigma3.c0.val = 4.2466090014400955

    s.set_model(m3)

    assert s.get_num_par() == 4
    assert s.get_num_par_thawed() == 2

    s.fit()
    res = s.get_fit_results()

    assert s.calc_stat() == pytest.approx(LPAR_STAT)
    assert m3.fwhm.val == pytest.approx(LPAR_FWHM)
    assert s3.c0.val == pytest.approx(LPAR_SIGMA)
    assert m3.ampl.val == pytest.approx(LPAR_AMPL)

    assert res.statval == pytest.approx(LPAR_STAT)
    assert res.numpoints == 12
    assert res.dof == 10
    assert res.parnames == ('fit3.ampl', 'sigma3.c0')
    assert res.parvals[0] == pytest.approx(LPAR_AMPL)
    assert res.parvals[1] == pytest.approx(LPAR_SIGMA)
