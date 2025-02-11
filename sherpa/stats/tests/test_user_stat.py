# -*- coding: utf-8 -*-
#
#  Copyright (C) 2017, 2019, 2025
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
import numpy as np
import pytest

from pytest import approx

from sherpa.astro import ui
from sherpa.data import Data1D
from sherpa.stats import UserStat
from sherpa.utils.err import IdentifierErr, StatErr


def test_user_stat_unit(clean_astro_ui):
    given_stat_error = [1.1, 2.2, 3.3]
    given_sys_error = [10.1, 10.2, 10.3]

    def calc_stat(data, _model, staterror, syserror=None, weight=None):
        # Make sure values are being injected correctly
        np.testing.assert_array_equal(given_stat_error, staterror)
        np.testing.assert_array_equal(given_sys_error, syserror)
        return 3.235, np.ones_like(data)

    xdata = [1, 2, 3]
    ydata = xdata

    ui.load_arrays(1, xdata, ydata, None, given_sys_error, Data1D)

    ui.set_model(1, 'polynom1d.p')

    ui.load_user_stat('customstat', calc_stat, lambda x: given_stat_error)
    ui.set_stat(eval('customstat'))

    try:
        ui.fit(1)
    except StatErr:
        pytest.fail("Call should not be throwing any exception (bug #341)")

    # Test the result is what we made the user stat return
    assert 3.235 == ui.get_fit_results().statval


@pytest.mark.parametrize("use_string", [False, True])
def test_user_model_stat_docs(use_string, clean_astro_ui):
    """This test reproduces the documentation shown at:
    http://cxc.harvard.edu/sherpa4.4/statistics/#userstat

    and:
    http://cxc.harvard.edu/sherpa/threads/user_model/

    I tried to be as faithful as possible to the original, although the examples in thedocs
    are not completely self-contained, so some changes were necessary. I changed the numpy
    reference, as it is imported as `np` here, and added a clean up of the environment
    before doing anything.

    For the model, the difference is that I am not importing the function from an
    external module, plus the dataset is different.

    Also, the stats docs do not perform a fit.

    The use_string parameter is to allow issue #2225 to be tested.
    When False, the original code is used, using eval to access the
    statistic name. When True we use a string (issue #2225).

    """
    def my_stat_func(data, model, staterror, syserror=None, weight=None):
        # A simple function to replicate χ2
        fvec = ((data - model) / staterror)**2
        stat = fvec.sum()
        return (stat, fvec)

    def my_staterr_func(data):
        # A simple staterror function
        return np.sqrt(data)

    def myline(pars, x):
        return pars[0]*x + pars[1]

    x = [1, 2, 3]
    y = [4, 5, 6.01]

    ui.clean()
    ui.load_arrays(1, x, y)
    ui.load_user_stat("mystat", my_stat_func, my_staterr_func)
    ui.load_user_model(myline, "myl")
    ui.add_user_pars("myl", ["m", "b"])

    if use_string:
        ui.set_stat("mystat")
    else:
        ui.set_stat(eval("mystat"))

    ui.set_model(eval("myl"))

    ui.fit()

    assert ui.get_par("myl.m").val == approx(1, abs=0.01)
    assert ui.get_par("myl.b").val == approx(3, abs=0.01)


def test_list_stats(clean_astro_ui):
    """Can we list the stats after adding a new one?"""

    # We do not hard code the statistic names here, in case we ever
    # add (or remove) any frmo the defaut set.
    #

    ostats = ui.list_stats()
    assert len(ostats) > 0
    assert len(ostats) == len(set(ostats))  # no duplicates
    assert ui.get_stat_name() in ostats

    def delme(*args, **kwargs):
        raise RuntimeError()

    statname = "not_A_Good_Stat"
    assert statname not in ostats          # just in case
    assert statname.lower() not in ostats  # just in case

    ui.load_user_stat(statname, delme)

    nstats = ui.list_stats()
    assert len(nstats) == len(ostats) + 1
    assert statname.lower() in nstats
    assert set(nstats).difference(set(ostats)) == set([statname.lower()])


@pytest.mark.parametrize("use_string", [False, True])
def test_get_stat_name_user_stat(use_string, clean_astro_ui):
    """What does get_stat_name return?"""

    def delme(*args, **kwargs):
        raise RuntimeError()

    statname = "What_A_Name"
    ui.load_user_stat(statname, delme)

    if use_string:
        ui.set_stat(statname)
    else:
        ui.set_stat(eval(statname))

    assert ui.get_stat_name() == statname.lower()


@pytest.mark.parametrize("use_string", [False, True])
def test_get_stat_user_stat(use_string, clean_astro_ui):
    """What does get_stat return?"""

    def delme(*args, **kwargs):
        raise RuntimeError()

    statname = "aName"
    ui.load_user_stat(statname, delme)

    if use_string:
        ui.set_stat(statname)
    else:
        ui.set_stat(eval(statname))

    s = ui.get_stat()

    assert isinstance(s, UserStat)
    assert s.name == statname
    assert s.statfunc == delme
    assert s.errfunc is None


@pytest.mark.parametrize("name", ["map", "Map"])
def test_load_user_stat_invalid_name(name, clean_astro_ui):
    """Check we error out"""

    with pytest.raises(IdentifierErr,
                       match=f"^'{name}' is reserved for the native Python function$"):
        ui.load_user_stat(name, None)


def test_341(clean_astro_ui):
    """
    The original reporter of bug #341 had a special implementation that should be captured
    by this test. The implementation has a proxy model that takes care of updating the actual
    model when it is evaluated. During a recent refactoring of the Stat and Fit code
    (PR #287) a regression was introduced by short-circuiting the evaluation of the model.

    """
    class ExampleModel():
        """ Class to define model
        """
        def __init__(self, x, y):
            self.x = np.array(x)
            self.y = np.array(y)
            self.parvals = [1, 2]
            self.parnames = ("m", "b")

        def calc_stat(self):
            return float(np.sum(np.abs(self.y - self.model())))

        def model(self):
            return self.parvals[0] * self.x + self.parvals[1]

    class CalcModel():
        """ Class to update model parameters
        """
        def __init__(self, model):
            self.model = model

        def __call__(self, pars, x):
            self.model.parvals = pars
            return np.ones_like(x)

    class CalcStat():
        """ Class to determine fit statistic
        """
        def __init__(self, model):
            self.model = model

        def __call__(self, _data, _model, *args, **kwargs):
            fit_stat = self.model.calc_stat()

            return fit_stat, np.ones(1)

    xdata = [1, 2, 3]
    ydata = [4, 5, 6]
    newmodel = ExampleModel(xdata, ydata)

    dummy_data = np.zeros(1)
    dummy_times = np.arange(1)
    ui.load_arrays(1, dummy_times, dummy_data)

    method = 'simplex'
    ui.set_method(method)

    ui.load_user_model(CalcModel(newmodel), 'simplemodel')
    ui.add_user_pars('simplemodel', newmodel.parnames)
    ui.set_model(1, 'simplemodel')

    calc_stat = CalcStat(newmodel)
    ui.load_user_stat('customstat', calc_stat, lambda x: np.ones_like(x))
    ui.set_stat(eval('customstat'))

    ui.fit(1)

    assert ui.get_par("simplemodel.m").val == approx(1, abs=0.00001)
    assert ui.get_par("simplemodel.b").val == approx(3, abs=0.00001)
