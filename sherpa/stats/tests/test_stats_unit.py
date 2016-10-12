#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

# Unit tests for changes to the stats module
#
# The chi-square tests are similar to those in
# sherpa/tests/test_fit_unit.py

import pytest
import numpy as np

from numpy.testing import assert_almost_equal, assert_equal

from sherpa.data import Data1D, DataSimulFit
from sherpa.models.model import SimulFitModel
from sherpa.models.basic import Polynom1D
from sherpa.utils.err import StatErr

from sherpa.stats import LeastSq, Chi2, Chi2Gehrels, Chi2DataVar, \
    Chi2ConstVar, Chi2ModVar, Chi2XspecVar, Cash, CStat, WStat, UserStat


def setup_single(stat, sys):
    """Return a single data set and model (as SimulFit objects).

    Parameters
    ----------
    stat, sys : bool
        Should statistical and systematic errors be explicitly set
        (True) or taken from the statistic (False)?

    Returns
    -------
    data, model
        DataSimulFit and SimulFitModel objects.

    """

    x = np.asarray([-10, -5, 3, 4])
    y = np.asarray([16.3, 2.4, 4.3, 5.6])
    if stat:
        # pick values close to sqrt(y)
        staterr = np.asarray([4.0, 1.6, 2.1, 2.4])
    else:
        staterr = None
    if sys:
        syserr = 0.05 * y
    else:
        syserr = None

    data = Data1D('tst1', x, y, staterror=staterr, syserror=syserr)

    mdl = Polynom1D('mdl1')
    mdl.c0 = 0
    mdl.c2 = 0.2
    mdl.offset = -1

    return DataSimulFit('tst', (data,)), SimulFitModel('mdl', (mdl,))


@pytest.mark.parametrize("stat", [Cash, CStat, WStat, UserStat])
def test_stats_calc_chisqr_missing(stat):
    """non chi-quare statistics do not have a calc_chisqr."""

    statobj = stat()
    assert not hasattr(statobj, 'calc_chisqr')


@pytest.mark.parametrize("usesys", [True, False])
def test_stats_calc_chisqr_chi2_nostat(usesys):
    """chi-quare statistic with no staterror fails"""

    data, models = setup_single(False, usesys)
    statobj = Chi2()
    with pytest.raises(StatErr):
        statobj.calc_chisqr(data, models)


# Numeric answers calculated using CIAO 4.8 (sherpa 4.8.0)
# on a 64-bit Linux machine.
#
exp_delta = np.asarray([0.1, 0.8, 1.1, 0.6])
exp_lsq = exp_delta * exp_delta

exp_chi2_tt = np.asarray([0.00060009, 0.24860162, 0.27153028, 0.06166073])
exp_chi2_tf = np.asarray([0.000625, 0.25, 0.27437642, 0.0625])

exp_dvar_ft = np.asarray([0.00058948, 0.26507621, 0.27840252, 0.06339814])
exp_dvar_ff = np.asarray([0.0006135, 0.26666667, 0.28139535, 0.06428571])

exp_cvar_ft = np.asarray([0.00127972, 0.08933058, 0.16814371, 0.04980355])
exp_cvar_ff = np.asarray([0.0013986, 0.08951049, 0.16923077, 0.05034965])

exp_mvar_xt = np.asarray([0.00059297, 0.19910403, 0.37274064, 0.07088847])
exp_mvar_xf = np.asarray([0.00061728, 0.2, 0.378125, 0.072])


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (Chi2, True, True, exp_chi2_tt),
    (Chi2, True, False, exp_chi2_tf),
    (LeastSq, True, True, exp_lsq),
    (LeastSq, True, False, exp_lsq),
    (LeastSq, False, True, exp_lsq),
    (LeastSq, False, False, exp_lsq),
    (Chi2Gehrels, True, True, exp_chi2_tt),
    (Chi2DataVar, True, True, exp_chi2_tt),
    (Chi2DataVar, True, False, exp_chi2_tf),
    (Chi2DataVar, False, True, exp_dvar_ft),
    (Chi2DataVar, False, False, exp_dvar_ff),
    (Chi2XspecVar, True, True, exp_chi2_tt),
    (Chi2XspecVar, True, False, exp_chi2_tf),
    (Chi2XspecVar, False, True, exp_dvar_ft),
    (Chi2XspecVar, False, False, exp_dvar_ff),
    (Chi2ConstVar, True, True, exp_chi2_tt),
    (Chi2ConstVar, True, False, exp_chi2_tf),
    (Chi2ConstVar, False, True, exp_cvar_ft),
    (Chi2ConstVar, False, False, exp_cvar_ff),
    (Chi2ModVar, True, True, exp_mvar_xt),
    (Chi2ModVar, True, False, exp_mvar_xf),
    (Chi2ModVar, False, True, exp_mvar_xt),
    (Chi2ModVar, False, False, exp_mvar_xf)
])
def test_stats_calc_chisqr(stat, usestat, usesys, expected):
    """chi-quare statistics calculates expected values"""

    data, models = setup_single(usestat, usesys)
    statobj = stat()
    answer = statobj.calc_chisqr(data, models)
    assert_almost_equal(answer, expected)


@pytest.mark.parametrize("stat,usestat,usesys,expected1", [
    (Chi2, True, True, exp_chi2_tt),
    (Chi2, True, False, exp_chi2_tf),
    (LeastSq, True, True, exp_lsq),
    (LeastSq, True, False, exp_lsq),
    (LeastSq, False, True, exp_lsq),
    (LeastSq, False, False, exp_lsq),
    (Chi2Gehrels, True, True, exp_chi2_tt),
    (Chi2DataVar, True, True, exp_chi2_tt),
    (Chi2DataVar, True, False, exp_chi2_tf),
    (Chi2DataVar, False, True, exp_dvar_ft),
    (Chi2DataVar, False, False, exp_dvar_ff),
    (Chi2XspecVar, True, True, exp_chi2_tt),
    (Chi2XspecVar, True, False, exp_chi2_tf),
    (Chi2XspecVar, False, True, exp_dvar_ft),
    (Chi2XspecVar, False, False, exp_dvar_ff),
    (Chi2ConstVar, True, True, exp_chi2_tt),
    (Chi2ConstVar, True, False, exp_chi2_tf),
    # See test_stats_calc_chisqr_multi_cvar for the next two
    # (Chi2ConstVar, False, True, exp_cvar_ft),
    # (Chi2ConstVar, False, False, exp_cvar_ff),
    (Chi2ModVar, True, True, exp_mvar_xt),
    (Chi2ModVar, True, False, exp_mvar_xf),
    (Chi2ModVar, False, True, exp_mvar_xt),
    (Chi2ModVar, False, False, exp_mvar_xf)
])
def test_stats_calc_chisqr_multi(stat, usestat, usesys, expected1):
    """chi-quare statistics calculates expected values

    Here multiple datasets are used. To save time, the first
    dataset is re-used, excluding the third point.
    """

    x, y = setup_single(usestat, usesys)

    data1 = x.datasets[0]
    model1 = y.parts[0]

    x, y = setup_single(usestat, usesys)

    data2 = x.datasets[0]
    data2.ignore(1, 3.5)

    # not an essential part of the stats code; more a check that
    # things are working correctly
    assert_equal(data2.mask, np.asarray([True, True, False, True]))

    mdata = DataSimulFit('simul', (data1, data2))
    mmodel = SimulFitModel('simul', (model1, model1))

    statobj = stat()
    answer = statobj.calc_chisqr(mdata, mmodel)

    expected = np.concatenate((expected1,
                               expected1[data2.mask]))
    assert_almost_equal(answer, expected)


@pytest.mark.parametrize("usestat,usesys,expected", [
    (False, True, np.concatenate((exp_cvar_ft,
                                  np.asarray([0.001141, 0.07887213,
                                              0.04401839])))),
    (False, False, np.concatenate((exp_cvar_ff,
                                   np.asarray([0.00123457, 0.07901235,
                                               0.04444444])))),
])
def test_stats_calc_chisqr_multi_cvar(usestat, usesys, expected):
    """chi-quare statistics calculates expected values: constvar

    Here multiple datasets are used. To save time, the first
    dataset is re-used, excluding the third point. This means that
    the constvar results are a little-bit different to the other
    statistics.
    """

    x, y = setup_single(usestat, usesys)

    data1 = x.datasets[0]
    model1 = y.parts[0]

    x, y = setup_single(usestat, usesys)

    data2 = x.datasets[0]
    data2.ignore(1, 3.5)

    # not an essential part of the stats code; more a check that
    # things are working correctly
    assert_equal(data2.mask, np.asarray([True, True, False, True]))

    mdata = DataSimulFit('simul', (data1, data2))
    mmodel = SimulFitModel('simul', (model1, model1))

    statobj = Chi2ConstVar()
    answer = statobj.calc_chisqr(mdata, mmodel)

    assert_almost_equal(answer, expected)
