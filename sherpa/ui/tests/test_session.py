#
#  Copyright (C) 2016, 2017, 2019, 2020, 2021
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
from unittest.mock import patch

from numpy.testing import assert_array_equal

import pytest

from sherpa.models import parameter
import sherpa.models.basic
from sherpa import estmethods as est
from sherpa import optmethods as opt
from sherpa import stats
from sherpa.ui.utils import Session
from sherpa.utils.err import ArgumentErr, ArgumentTypeErr, \
    IdentifierErr, SessionErr
from sherpa.utils.logging import SherpaVerbosity
from sherpa.utils.testing import requires_plotting, requires_xspec

TEST = [1, 2, 3]
TEST2 = [4, 5, 6]


# bug #303
#
# Note: this test could be run even if plotting is not available,
# since in that case the initial get_data_plot_prefs() call will
# return an empty dictionary rather than a filled one.
#
@requires_plotting
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
def test_save_restore(tmpdir):
    outfile = tmpdir.join("sherpa.save")
    session = Session()
    session.load_arrays(1, TEST, TEST2)
    session.save(str(outfile), clobber=True)
    session.clean()
    assert set() == set(session.list_data_ids())

    session.restore(str(outfile))
    assert {1, } == set(session.list_data_ids())
    assert_array_equal(TEST, session.get_data(1).get_indep()[0])
    assert_array_equal(TEST2, session.get_data(1).get_dep())


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
    with patch("sys.stdin", StringIO(u"2.1")):
        s.set_model('const1d.mx')

    assert s.list_model_ids() == [1]
    expected.add('mx')
    assert set(s.list_model_components()) == expected

    mx = s.get_model_component('mx')
    assert mx.c0.val == pytest.approx(2.1)

    with patch("sys.stdin", StringIO(u"2.1e-3 , 2.0e-3, 1.2e-2")):
        s.set_model('x', 'const1d.my')

    assert set(s.list_model_ids()) == set([1, 'x'])
    expected.add('my')
    assert set(s.list_model_components()) == expected

    my = s.get_model_component('my')
    assert my.c0.val == pytest.approx(2.1e-3)
    assert my.c0.min == pytest.approx(2.0e-3)
    assert my.c0.max == pytest.approx(1.2e-2)

    with patch("sys.stdin", StringIO(u"2.1e-3 ,, 1.2e-2")):
        s.set_model(2, 'const1d.mz')

    assert set(s.list_model_ids()) == set([1, 2, 'x'])
    expected.add('mz')
    assert set(s.list_model_components()) == expected

    mz = s.get_model_component('mz')
    assert mz.c0.val == pytest.approx(2.1e-3)
    assert mz.c0.min == pytest.approx(- parameter.hugeval)
    assert mz.c0.max == pytest.approx(1.2e-2)

    with patch("sys.stdin", StringIO(u"-12.1,-12.1")):
        s.set_model('const1d.f1')

    assert set(s.list_model_ids()) == set([1, 2, 'x'])
    expected.add('f1')
    assert set(s.list_model_components()) == expected

    f1 = s.get_model_component('f1')
    assert f1.c0.val == pytest.approx(-12.1)
    assert f1.c0.min == pytest.approx(-12.1)
    assert f1.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(u" ,-12.1")):
        s.set_model('const1d.f2')

    # stop checking list_model_ids
    expected.add('f2')
    assert set(s.list_model_components()) == expected

    f2 = s.get_model_component('f2')
    assert f2.c0.val == pytest.approx(1.0)
    assert f2.c0.min == pytest.approx(-12.1)
    assert f2.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(u" ,")):
        s.set_model('const1d.f3')

    f3 = s.get_model_component('f3')
    assert f3.c0.val == pytest.approx(1.0)
    assert f3.c0.min == pytest.approx(- parameter.hugeval)
    assert f3.c0.max == pytest.approx(parameter.hugeval)

    with patch("sys.stdin", StringIO(u" ,, ")):
        s.set_model('const1d.f4')

    f4 = s.get_model_component('f4')
    assert f4.c0.val == pytest.approx(1.0)
    assert f4.c0.min == pytest.approx(- parameter.hugeval)
    assert f4.c0.max == pytest.approx(parameter.hugeval)

    # TODO: test error cases
    #


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
    assert 'primini' in methods


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
    with pytest.raises(ArgumentTypeErr) as exc:
        s.get_method(opt.MonCar)

    assert str(exc.value) == "'name' must be a string"


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
    with pytest.raises(ArgumentTypeErr) as exc:
        s.set_method(sherpa.models.basic.Const1D)

    assert str(exc.value) == "'meth' must be a method name or object"


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
    with pytest.raises(ArgumentTypeErr) as exc:
        s.get_stat(stats.Cash)

    assert str(exc.value) == "'name' must be a string"


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
    with pytest.raises(ArgumentTypeErr) as exc:
        s.set_stat(sherpa.models.basic.Const1D)

    assert str(exc.value) == "'stat' must be a statistic name or object"


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
    func = getattr(s, 'get_{}'.format(name))
    ans = func()
    assert isinstance(ans, req)


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_get_error_opt(name):
    """Can we get an option for an error estimator?

    Fortunately they all have a sigma option.
    """

    s = Session()
    func = getattr(s, 'get_{}_opt'.format(name))
    ans = func('sigma')
    assert ans == pytest.approx(1.0)


@pytest.mark.parametrize("name,fullname",
                         [('covar', 'covariance'),
                          ('conf', 'confidence'),
                          ('proj', 'projection')])
def test_get_error_opt_invalid(name, fullname):
    """We can not ask for an error option that does not exist"""

    s = Session()
    func = getattr(s, 'get_{}_opt'.format(name))

    with pytest.raises(ArgumentErr) as exc:
        func('the-real-sigma')

    assert str(exc.value) == "'the-real-sigma' is not a valid option for method {}".format(fullname)


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_set_error_opt(name):
    """Can we set an option for an error estimator?

    Fortunately they all have a sigma option.
    """

    s = Session()
    gfunc = getattr(s, 'get_{}_opt'.format(name))
    sfunc = getattr(s, 'set_{}_opt'.format(name))

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
    func = getattr(s, 'set_{}_opt'.format(name))

    with pytest.raises(ArgumentErr) as exc:
        func('the-real-sigma', 2.4)

    assert str(exc.value) == "'the-real-sigma' is not a valid option for method {}".format(fullname)


@pytest.mark.parametrize("name,fullname",
                         [('covar', 'covariance'),
                          ('conf', 'confidence'),
                          ('proj', 'projection')])
def test_error_estimate_not_set(name, fullname):
    """Error out if the error estimate is not set"""

    s = Session()
    func = getattr(s, 'get_{}_results'.format(name))

    with pytest.raises(SessionErr) as exc:
        func()

    assert str(exc.value) == "{} has not been performed".format(fullname)


@pytest.mark.parametrize("name", ['covar', 'conf', 'proj'])
def test_error_estimate_not_run(name):
    """Can not rn the error estimate if data is not set up"""

    s = Session()
    s.set_default_id('bob')
    func = getattr(s, name)

    with pytest.raises(IdentifierErr) as exc:
        func()

    assert str(exc.value) == "data set bob has not been set"


def test_set_source_invalid():

    s = Session()
    with pytest.raises(ArgumentErr) as ae:
        s.set_source('2 * made_up.foo')

    assert str(ae.value) == "invalid model expression: name 'made_up' is not defined"


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


def test_add_user_pars_modelname_not_a_string():

    s = Session()
    with pytest.raises(ArgumentTypeErr) as exc:
        s.add_user_pars(23, ['x'])

    assert str(exc.value) == "'model name' must be a string"


def test_add_user_pars_modelname_not_a_model1():

    s = Session()
    with pytest.raises(ArgumentTypeErr) as exc:
        s.add_user_pars('not a model', ['x'])

    assert str(exc.value) == "'not a model' must be a user model"


def test_add_user_pars_modelname_not_a_model2():
    """Use an actual model, but not a user model"""

    s = Session()
    s._add_model_types(sherpa.models.basic)
    s.create_model_component('scale1d', 'foo')
    with pytest.raises(ArgumentTypeErr) as exc:
        s.add_user_pars('foo', ['x'])

    assert str(exc.value) == "'foo' must be a user model"
