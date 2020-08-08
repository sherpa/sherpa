#
#  Copyright (C) 2016, 2017, 2019, 2020  Smithsonian Astrophysical Observatory
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
import numpy

from sherpa.astro.data import DataPHA
from sherpa.utils.testing import requires_plotting, requires_xspec
from sherpa.ui.utils import Session
from sherpa.astro.ui.utils import Session as AstroSession
from numpy.testing import assert_array_equal
from sherpa.models import parameter, Const1D
import sherpa.models.basic

import pytest

from unittest.mock import patch

from io import StringIO


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
def test_models_all():
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


# Fix 476 - this should be in sherpa/ui/tests/test_session.py
def test_zero_division_calc_stat():
    ui = AstroSession()
    x = numpy.arange(100)
    y = numpy.zeros(100)
    ui.load_arrays(1, x, y, DataPHA)
    ui.group_counts(1, 100)
    ui.set_full_model(1, Const1D("const"))

    # in principle I wouldn't need to call calc_stat_info(), I could just
    # use _get_stat_info to reproduce the issue, However, _get_stat_info is not a public
    # method, so I want to double check that calc_stat_info does not throw an exception.
    # So, first we try to run calc_stat_info and make sure there are no exceptions.
    # Then, since calc_stat_info only logs something and doesn't return anything, we use
    # a white box approach to get the result from _get_stat_info.
    ui.calc_stat_info()
    assert ui._get_stat_info()[0].rstat is numpy.nan
