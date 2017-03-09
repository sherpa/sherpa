#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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

from sherpa.astro import ui
from sherpa.data import Data1D
from sherpa.utils.err import StatErr


def test_user_stat_unit():
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
    ui.set_stat(customstat)

    try:
        ui.fit(1)
    except StatErr:
        pytest.fail("Call should not be throwing and exception (bug #341)")

    # Test the result is what we made the user stat return
    assert 3.235 == ui.get_fit_results().statval


def test_user_model_stat_docs():
    """
    This test reproduces the documentation shown at:
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
    ui.set_stat(mystat)
    ui.load_user_model(myline, "myl")
    ui.add_user_pars("myl", ["m", "b"])
    ui.set_model(myl)

    ui.fit()

    assert abs(1 - ui.get_par("myl.m").val) < 0.01
    assert abs(3 - ui.get_par("myl.b").val) < 0.01


def test_341():
    """
    The original reporter of bug #341 had a special implementation that should be captured
    by this test. The implementation has a proxy model that takes care of updating the actual
    model when it is evaluated. During a recent refactoring of the Stat and Fit code
    (PR #287) a regression was introduced by short-circuiting the evaluation of the model.

    """
    class ExampleModel(object):
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

    class CalcModel(object):
        """ Class to update model parameters
        """
        def __init__(self, model):
            self.model = model

        def __call__(self, pars, x):
            self.model.parvals = pars
            return np.ones_like(x)

    class CalcStat(object):
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
    ui.set_stat(customstat)

    ui.fit(1)
    assert abs(1 - ui.get_par("simplemodel.m").val) < 0.00001
    assert abs(3 - ui.get_par("simplemodel.b").val) < 0.00001
