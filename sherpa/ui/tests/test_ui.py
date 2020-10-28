#
#  Copyright (C) 2012, 2015, 2016, 2018, 2019, 2020
#      Smithsonian Astrophysical Observatory
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

from collections import namedtuple
import logging

import numpy as np

import pytest

from sherpa.utils.testing import requires_data
from sherpa.data import Data1D
from sherpa.fit import Fit
from sherpa.stats import LeastSq
from sherpa.models import ArithmeticModel, Parameter
from sherpa.models.basic import PowLaw1D
from sherpa.models.parameter import hugeval
from sherpa.utils.err import IdentifierErr, StatErr, SessionErr
from sherpa import ui

# As the user model is called UserModel, refer to the Sherpa version
# by its full module name
import sherpa.models.basic


logger = logging.getLogger("sherpa")


@pytest.fixture(autouse=True)
def hide_logging():
    "hide INFO-level logging in all these tests"

    olevel = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(olevel)


class UserModel(ArithmeticModel):
    def __init__(self, name='usermodel'):
        self.param1 = Parameter(name, 'param1', 1, min=0, max=100)
        self.param2 = Parameter(name, 'param2', 1, min=-100, max=100)

        ArithmeticModel.__init__(self, name, (self.param1,
                                              self.param2))

    def calc(self, p, x, *args, **kwargs):
        return p[0] * x + p[1]


# A function-based version of this
#
def um_line(pars, xlo, *args, **kwargs):
    """A test user model: straight line with slope and intercept"""
    return pars[0] * np.asarray(xlo) + pars[1]


"""
Tests for PR #155

TODO: Some test cases cannot be implemented (easily) without mock
In particular one cannot test that the correct matrix is used when one is not
provided and that an error is thrown when, for any reason, the computed covariance matrix
is None.
"""


WRONG_STAT_MSG = "Fit statistic must be cash, cstat or wstat, not {}"
WSTAT_ERR_MSG = "No background data has been supplied. Use cstat"
NO_COVAR_MSG = "covariance has not been performed"
RIGHT_STATS = {'cash', 'cstat', 'wstat'}


@pytest.fixture
def setup_covar(make_data_path):

    print("A")
    ui.load_data(make_data_path('sim.poisson.1.dat'))
    ui.set_model(PowLaw1D("p1"))


@requires_data
@pytest.mark.parametrize('stat', sorted(set(ui.list_stats()) - RIGHT_STATS))
def test_covar_wrong_stat(stat, clean_ui, setup_covar):
    "Test an exception is thrown is the proper stat is not set"

    ui.covar()
    ui.set_stat(stat)
    with pytest.raises(ValueError) as exc:
        ui.get_draws()

    assert WRONG_STAT_MSG.format(stat) == str(exc.value)


@requires_data
def test_covar_wstat_no_background(clean_ui, setup_covar):
    "Test an exception is thrown when wstat is used without background"

    ui.covar()
    ui.set_stat("wstat")
    with pytest.raises(StatErr) as exc:
        ui.get_draws()

    assert WSTAT_ERR_MSG == str(exc.value)


@requires_data
@pytest.mark.parametrize("stat", sorted(RIGHT_STATS))
def test_no_covar(stat, clean_ui, setup_covar):
    "Test an exception is thrown if covar is not run"

    ui.set_stat(stat)
    with pytest.raises(SessionErr) as exc:
        ui.get_draws()

    assert NO_COVAR_MSG == str(exc.value)


# Test get_draws returns a valid response when the covariance matrix is provided
# Note the accuracy of the returned values is not assessed here
@requires_data
@pytest.mark.parametrize('stat', sorted(RIGHT_STATS - {'wstat'}))
def test_covar_as_argument(stat, clean_ui, setup_covar):

    ui.set_stat(stat)
    ui.fit()

    matrix = [[0.00064075, 0.01122127], [0.01122127, 0.20153251]]
    niter = 10
    stat, accept, params = ui.get_draws(niter=niter, covar_matrix=matrix)

    n = niter + 1
    assert stat.size == n
    assert accept.size == n
    assert params.shape == (2, n)
    assert np.any(accept)


# Test get_draws returns a valid response when the covariance matrix is not provided
# Note the accuracy of the returned values is not assessed here
@requires_data
@pytest.mark.parametrize('stat', sorted(RIGHT_STATS - {'wstat'}))
def test_covar_as_none(stat, clean_ui, setup_covar):

    ui.set_stat(stat)
    ui.fit()
    ui.covar()

    niter = 10
    stat, accept, params = ui.get_draws(niter=niter)

    n = niter + 1
    assert stat.size == n
    assert accept.size == n
    assert params.shape == (2, n)
    assert np.any(accept)


def identity(x):
    return x


@pytest.fixture
def setup_ui(make_data_path):

    ui.dataspace1d(1, 1000, dstype=ui.Data1D)

    out = namedtuple('setup_ui_data', ['ascii', 'single', 'double',
                                       'filter'])

    out.ascii = make_data_path('sim.poisson.1.dat')
    out.single = make_data_path('single.dat')
    out.double = make_data_path('double.dat')
    out.filter = make_data_path('filter_single_integer.dat')
    return out


@requires_data
def test_ui_ascii(clean_ui, setup_ui):
    ui.load_data(1, setup_ui.ascii)
    ui.load_data(1, setup_ui.ascii, 2)
    ui.load_data(1, setup_ui.ascii, 2, ("col2", "col1"))


# Test table model
@requires_data
def test_ui_table_model_ascii_table(clean_ui, setup_ui):
    ui.load_table_model('tbl', setup_ui.single)
    ui.load_table_model('tbl', setup_ui.double)


# Test user model
@requires_data
def test_ui_user_model_ascii_table(clean_ui, setup_ui):
    ui.load_user_model(identity, 'mdl', setup_ui.single)
    ui.load_user_model(identity, 'mdl', setup_ui.double)


@requires_data
def test_ui_filter_ascii(clean_ui, setup_ui):
    ui.load_filter(setup_ui.filter)
    ui.load_filter(setup_ui.filter, ignore=True)


@requires_data
def test_ui_add_model(clean_ui, setup_ui):
    ui.add_model(UserModel)
    ui.set_model('usermodel.user1')


@requires_data
def test_ui_set_full_model(clean_ui, setup_ui):
    ui.load_psf('psf1', 'gauss2d.g1')
    ui.set_full_model('psf1(gauss2d.g2)+const2d.c1')
    ui.get_model()


# Bug 12644
@requires_data
def test_ui_source_methods_with_full_model(clean_ui, setup_ui):

    ui.load_data('full', setup_ui.ascii)
    ui.set_full_model('full', 'powlaw1d.p1')

    # Test Case 1
    with pytest.raises(IdentifierErr) as exc:
        ui.get_source('full')

    emsg = "Convolved model\n'powlaw1d.p1'\n is set for dataset full. You should use get_model instead."
    assert str(exc.value) == emsg

    with pytest.raises(IdentifierErr) as exc:
        ui.plot_source('full')

    emsg = "Convolved model\n'powlaw1d.p1'\n is set for dataset full. You should use plot_model instead."
    assert str(exc.value) == emsg

    with pytest.raises(IdentifierErr) as exc:
        ui.get_source_plot('full')

    emsg = "Convolved model\n'powlaw1d.p1'\n is set for dataset full. You should use get_model_plot instead."
    assert str(exc.value) == emsg

    # Test Case 2
    ui.set_source('full', 'powlaw1d.p2')
    ui.get_source('full')

    # Test Case 3
    ui.load_data('not_full', setup_ui.ascii)
    with pytest.raises(IdentifierErr) as exc:
        ui.get_source('not_full')

    emsg = 'source not_full has not been set, consider using set_source() or set_model()'
    assert emsg == str(exc.value)


@pytest.mark.parametrize('model', ['gauss1d', 'delta1d', 'normgauss1d'])
def test_psf_model1d(model, clean_ui):
    ui.dataspace1d(1, 10)
    ui.load_psf('psf1d', model + '.mdl')
    ui.set_psf('psf1d')
    mdl = ui.get_model_component('mdl')
    assert mdl.get_center() == (4.0, )


@pytest.mark.parametrize('model', ['gauss2d', 'delta2d', 'normgauss2d'])
def test_psf_model2d(model):
    ui.dataspace2d([216, 261])
    ui.load_psf('psf2d', model + '.mdl')
    ui.set_psf('psf2d')
    mdl = ui.get_model_component('mdl')
    assert mdl.get_center() == (108.0, 130.0)


def test_does_user_model_get_cleaned(clean_ui):
    """Do user models get removed from the session by clean?"""

    mname = "test_model"
    with pytest.raises(IdentifierErr):
        ui.get_model_component(mname)

    ui.load_user_model(um_line, mname)
    mdl = ui.get_model_component(mname)
    assert mdl.name == "usermodel.{}".format(mname)
    assert isinstance(mdl, sherpa.models.basic.UserModel)

    ui.clean()

    with pytest.raises(IdentifierErr):
        ui.get_model_component(mname)


# Test simple use of load_user_model and add_user_pars
#
def test_user_model_create_pars_default(clean_ui):

    mname = "test_model"
    ui.load_user_model(um_line, mname)

    mdl = ui.get_model_component(mname)
    assert len(mdl.pars) == 1
    par = mdl.pars[0]
    assert par.name == "ampl"
    assert par.val == pytest.approx(1.0)
    assert not par.frozen
    assert par.units == ''
    assert par.min == pytest.approx(-1 * hugeval)
    assert par.max == pytest.approx(hugeval)


def test_user_model_create_pars_names(clean_ui):

    mname = "test_model"
    ui.load_user_model(um_line, mname)

    mdl = ui.get_model_component(mname)
    assert len(mdl.pars) == 1

    # add user pars doesn't change the existing instance, you have
    # to "get" the new version to see the change
    #
    ui.add_user_pars(mname, ['X1', 'x'])

    mdl = ui.get_model_component(mname)
    assert len(mdl.pars) == 2
    p0 = mdl.pars[0]
    p1 = mdl.pars[1]

    assert p0.name == 'X1'
    assert p0.val == pytest.approx(0.0)
    assert p0.units == ''
    assert not p0.frozen
    assert p0.min == pytest.approx(-1 * hugeval)
    assert p0.max == pytest.approx(hugeval)

    assert p1.name == 'x'
    assert p1.val == pytest.approx(0.0)
    assert p1.units == ''
    assert not p1.frozen
    assert p1.min == pytest.approx(-1 * hugeval)
    assert p1.max == pytest.approx(hugeval)


def test_user_model_create_pars_full(clean_ui):

    mname = "test_model"
    ui.load_user_model(um_line, mname)
    ui.add_user_pars(mname, ['pAr1', '_p'], [23.2, 3.1e2],
                     parunits=['', 'cm^2 s'],
                     parfrozen=[True, False],
                     parmins=[0, -100],
                     parmaxs=[100, 1e5])

    mdl = ui.get_model_component(mname)
    assert len(mdl.pars) == 2
    p0 = mdl.pars[0]
    p1 = mdl.pars[1]

    assert p0.name == 'pAr1'
    assert p0.val == pytest.approx(23.2)
    assert p0.units == ''
    assert p0.frozen
    assert p0.min == pytest.approx(0)
    assert p0.max == pytest.approx(100)

    assert p1.name == '_p'
    assert p1.val == pytest.approx(3.1e2)
    assert p1.units == 'cm^2 s'
    assert not p1.frozen
    assert p1.min == pytest.approx(-100)
    assert p1.max == pytest.approx(1e5)


# diagnose and test out GitHub issue #609; the inability to change
# parameter values for a user model using direct access
#
def test_user_model_change_par(clean_ui):

    mname = "test_model"
    ui.load_user_model(um_line, mname)
    ui.add_user_pars(mname, ['xXx', 'Y2'])

    mdl = ui.get_model_component(mname)
    assert len(mdl.pars) == 2
    p0 = mdl.pars[0]
    p1 = mdl.pars[1]

    assert p0.name == 'xXx'
    assert p1.name == 'Y2'
    assert p0.val == pytest.approx(0.0)
    assert p1.val == pytest.approx(0.0)

    # Use the user-supplied names:
    #
    mdl.xXx = 2.0
    assert p0.val == pytest.approx(2.0)

    mdl.Y2 = 3.0
    assert p1.val == pytest.approx(3.0)

    # Now all lower case
    #
    mdl.xxx = 4.0
    assert p0.val == pytest.approx(4.0)

    mdl.y2 = 12.0
    assert p1.val == pytest.approx(12.0)

    # Try with the set_par function
    #
    ui.set_par('test_model.xxx', 12.2)
    assert p0.val == pytest.approx(12.2)

    ui.set_par('test_model.y2', 14.0, frozen=True)
    assert p1.val == pytest.approx(14.0)
    assert p1.frozen

    ui.clean()


def test_user_model1d_eval_fail(clean_ui):
    """This is expected to fail as the number of pars does not match."""

    mname = "test_model"
    ui.load_user_model(um_line, mname)
    mdl = ui.get_model_component(mname)
    with pytest.raises(IndexError):
        mdl([2.3, 5.4, 8.7])


def test_user_model1d_eval(clean_ui):
    """Simple evaluation check for 1D case."""

    mname = "test_model"
    ui.load_user_model(um_line, mname)
    ui.add_user_pars(mname, ["slope", "intercept"])

    m = 2.1
    c = -4.8

    mdl = ui.get_model_component(mname)
    mdl.slope = m
    mdl.intercept = c

    x = np.asarray([2.3, 5.4, 8.7])
    y = mdl(x)

    yexp = x * m + c

    # This check require pytest >= 3.2.0
    #
    assert y == pytest.approx(yexp)


def test_user_model1d_fit(clean_ui):
    """Check can use in a fit."""

    mname = "test_model"
    ui.load_user_model(um_line, mname)
    ui.add_user_pars(mname, ["slope", "intercept"],
                     parvals=[1.0, 1.0])

    mdl = ui.get_model_component(mname)

    x = np.asarray([-2.4, 2.3, 5.4, 8.7, 12.3])

    # Set up the data to be scattered around y = -0.2 x + 2.8
    # Pick the deltas so that they sum to 0 (except for central
    # point)
    #
    slope = -0.2
    intercept = 2.8

    dy = np.asarray([0.1, -0.2, 0.14, -0.1, 0.2])
    ydata = x * slope + intercept + dy

    ui.load_arrays(1, x, ydata)

    ui.set_source(mname)
    ui.ignore(5.0, 6.0)  # drop the central bin

    ui.set_stat('leastsq')
    ui.set_method('simplex')
    ui.fit()

    fres = ui.get_fit_results()
    assert fres.succeeded
    assert fres.parnames == ('test_model.slope', 'test_model.intercept')
    assert fres.numpoints == 4
    assert fres.dof == 2

    # Tolerance has been adjusted to get the tests to pass on my
    # machine. It's really just to check that the values have chanegd
    # from their default values.
    #
    assert fres.parvals[0] == pytest.approx(slope, abs=0.01)
    assert fres.parvals[1] == pytest.approx(intercept, abs=0.05)

    # Thse should be the same values, so no need to use pytest.approx
    # (unless there's some internal translation between types done
    # somewhere?).
    #
    assert mdl.slope.val == fres.parvals[0]
    assert mdl.intercept.val == fres.parvals[1]


class MyCacheTestModel(ArithmeticModel):

    def calc(self, par, x):
        A = par[0]
        mylambda = par[1]
        b = par[2]
        fvec = A * np.exp(- mylambda * x) + b
        if self.counter == 0:
            assert self._use_caching
        else:
            assert not self._use_caching
        self.counter += 1
        return fvec

    def __init__(self, name='myexp'):
        self.A = Parameter(name, 'A', 1)
        self.mylambda = Parameter(name, 'mylambda', 2)
        self.b = Parameter(name, 'b', 3)
        self.counter = 0
        ArithmeticModel.__init__(self, name, (self.A, self.mylambda, self.b))


def test_cache():
    """To make sure that the runtime fit(cache=???) works"""

    x = np.array([1.0, 2.0, 3.0])
    model = MyCacheTestModel()
    par = np.array([1.1, 2.0, 3.0])
    y = model.calc(par, x)

    data = Data1D('tmp', x, y)
    fit = Fit(data, model, LeastSq())
    fit.fit(cache=False)


def test_list_psf_ids_empty(clean_ui):
    assert ui.list_psf_ids() == []


@pytest.mark.parametrize("id", [None, 2, "2"])
def test_list_psf_ids_single(id, clean_ui):
    ui.dataspace1d(1, 10, id=id)
    ui.load_psf('psf1d', 'gauss1d.mdl')
    ui.set_psf(id, 'psf1d')

    if id is None:
        id = 1

    assert ui.list_psf_ids() == [id]


def test_list_psf_ids_multi(clean_ui):

    ui.dataspace1d(1, 10)
    ui.dataspace1d(1, 10, id="2")
    ui.load_psf('psf1d', 'gauss1d.mdl')
    ui.set_psf('psf1d')
    ui.set_psf("2", 'psf1d')

    ans = ui.list_psf_ids()
    assert ans == [1, "2"]
