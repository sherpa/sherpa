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
from sherpa.utils.err import StatErr


def test_user_stat_unit():
    def calc_stat(data, _model):
        return 3.235, np.ones_like(data)

    xdata = [1, 2, 3]
    ydata = xdata

    ui.load_arrays(1, xdata, ydata)

    ui.set_model(1, 'polynom1d.p')

    ui.load_user_stat('customstat', calc_stat, lambda x: np.ones_like(x))
    ui.set_stat(customstat)

    try:
        ui.fit(1)
    except StatErr:
        pytest.fail("Call should not be throwing and exception (bug #341)")

    # Test the result is what we made the user stat return
    assert 3.235 == ui.get_fit_results().statval


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
