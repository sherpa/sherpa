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

# Unit tests for changes to the stats module. Since this provides
# the primary interface to the stats C++ code, checks are made here
# of the results of that code (in part to make sure the data is being
# sent to these routines correctly).
#
# The calc_info/chisqr tests are similar to those in
# sherpa/tests/test_fit_unit.py; once there are enough here it may
# be possible to reduce the number of tests in test_fit_unit.py.
#
# TODO:
#   add calc_stat_from_data tests for wstat
#
#   add some tests with Data2D and/or DataIMG (to check nD where n > 1)
#
#   add some error checks for mixed data types - e.g. Data1D and Data1DInt
#   as well as those with and without errors (e.g. to check that they are
#   caught and raise an error if necessary)
#

import pytest
import numpy as np

from numpy.testing import assert_almost_equal, assert_equal

from sherpa.astro.data import DataPHA
from sherpa.data import Data1D, Data1DInt, DataSimulFit
from sherpa.models.model import SimulFitModel
from sherpa.models.basic import Const1D, Polynom1D
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
        DataSimulFit and SimulFitModel objects. The data is a Data1D
        object.

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


def setup_single_1dint(stat, sys):
    """Return a single data set and model (as SimulFit objects).

    Parameters
    ----------
    stat, sys : bool
        Should statistical and systematic errors be explicitly set
        (True) or taken from the statistic (False)?

    Returns
    -------
    data, model
        DataSimulFit and SimulFitModel objects. The data is a
        Data1DInt object.

    """

    xlo = np.asarray([-10, -5, 3, 4])
    xhi = np.asarray([-5, 2, 4, 7])

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

    data = Data1DInt('tst1', xlo, xhi, y, staterror=staterr, syserror=syserr)

    # As the model is integrated, the normalization values need to
    # be divided by the bin width, but the bins are not constant
    # width, so pick a value which gives reasonable values.
    #
    mdl = Polynom1D('mdl1')
    mdl.c0 = 0
    mdl.c2 = 0.08
    mdl.offset = -1

    mdl.offset.frozen = False
    mdl.c2.frozen = False

    return DataSimulFit('tst', (data,)), SimulFitModel('mdl', (mdl,))


def setup_multiple(usestat, usesys):
    """Return multiple data sets and model (as SimulFit objects).

    To save time, the first dataset is re-used, excluding the third point.

    Parameters
    ----------
    stat, sys : bool
        Should statistical and systematic errors be explicitly set
        (True) or taken from the statistic (False)?

    Returns
    -------
    data, model
        DataSimulFit and SimulFitModel objects. The data sets
        are Data1D objects.
    """

    x, y = setup_single(usestat, usesys)
    data1 = x.datasets[0]
    model1 = y.parts[0]

    x, y = setup_single(usestat, usesys)
    data2 = x.datasets[0]
    data2.ignore(1, 3.5)

    # not an essential part of the stats code; more a check that
    # things are working correctly (i.e. that this invariant hasn't
    # been broken by some change to the test).
    assert_equal(data2.mask, np.asarray([True, True, False, True]))

    mdata = DataSimulFit('simul', (data1, data2))
    mmodel = SimulFitModel('simul', (model1, model1))
    return mdata, mmodel


def setup_multiple_1dint(stat, sys):
    """Return multiple data sets and models (as SimulFit objects).

    To save time, the first dataset is re-used, excluding the third point.

    Parameters
    ----------
    stat, sys : bool
        Should statistical and systematic errors be explicitly set
        (True) or taken from the statistic (False)?

    Returns
    -------
    data, model
        DataSimulFit and SimulFitModel objects. The data sets are
        Data1DInt objects.

    """

    x, y = setup_single_1dint(stat, sys)
    data1 = x.datasets[0]
    model1 = y.parts[0]

    x, y = setup_single_1dint(stat, sys)
    # The bins cover (-10,-5), (-5,2), (3,4), (4,7)
    data2 = x.datasets[0]
    data2.ignore(2.5, 3.5)

    # not an essential part of the stats code; more a check that
    # things are working correctly
    assert_equal(data2.mask, np.asarray([True, True, False, True]))

    mdata = DataSimulFit('simul', (data1, data2))
    mmodel = SimulFitModel('simul', (model1, model1))
    return mdata, mmodel


def setup_single_pha(stat, sys, background=True):
    """Return a single data set and model (as SimulFit objects).

    This is aimed at wstat calculation, and so the DataPHA object has
    no attached response. The data set is grouped.

    Parameters
    ----------
    stat, sys : bool
        Should statistical and systematic errors be explicitly set
        (True) or taken from the statistic (False)?
    background : bool
        Should a background data set be included (True) or not (False)?
        The background is *not* subtracted when True.

    Returns
    -------
    data, model
        DataSimulFit and SimulFitModel objects. The data is a DataPHA
        object.

    """

    # If used the same bins as setup_single_1dint then could
    # re-use the results, but the bins are different, and it
    # is useful for the Data1DInt case to test non-consecutive
    # histogram bins.
    #
    channels = np.arange(1, 6, dtype=np.int16)
    counts = np.asarray([10, 13, 9, 17, 21], dtype=np.int16)

    if stat:
        staterror = np.asarray([3.0, 4.0, 3.0, 4.0, 5.0])
    else:
        staterror = None

    if sys:
        syserror = 0.2 * counts
    else:
        syserror = None

    grouping = np.asarray([1, -1, 1, -1, 1], dtype=np.int16)
    # quality = np.asarray([0, 0, 0, 0, 0], dtype=np.int16)
    quality = None

    exposure = 150.0
    backscal = 0.01

    # does not set areascal or header
    data = DataPHA(name='tstpha', channel=channels,
                   counts=counts, staterror=staterror,
                   syserror=syserror, grouping=grouping,
                   quality=quality, exposure=exposure,
                   backscal=backscal)

    if background:
        bgcounts = np.asarray([2, 1, 0, 2, 2], dtype=np.int16)

        if stat:
            bgstaterror = np.asarray([0.2, 0.4, 0.5, 0.3, 0.2])
        else:
            bgstaterror = None

        if sys:
            bgsyserror = 0.3 * bgcounts
        else:
            bgsyserror = None

        bggrouping = None
        bgquality = None

        bgexposure = 550.0
        bgbackscal = np.asarray([0.05, 0.06, 0.04, 0.04, 0.07])

        bgdata = DataPHA(name='bgpha', channel=channels,
                         counts=bgcounts, staterror=bgstaterror,
                         syserror=bgsyserror, grouping=bggrouping,
                         quality=bgquality, exposure=bgexposure,
                         backscal=bgbackscal)

        data.set_background(bgdata)

    # Trying a multi-component model, even though this actual
    # model is degenerate (cnst.c0 and poly.c0)
    cnst = Const1D('cnst')
    poly = Polynom1D('poly')

    cnst.c0 = 1.2
    poly.c0 = 7.9
    poly.c1 = 2.1
    poly.c1.frozen = False

    mdl = cnst + poly
    return DataSimulFit('tst', (data,)), SimulFitModel('mdl', (mdl,))


@pytest.mark.parametrize("stat", [Cash, CStat, WStat, UserStat])
def test_stats_calc_chisqr_missing(stat):
    """non chi-quare statistics do not have a calc_chisqr."""

    statobj = stat()
    assert not hasattr(statobj, 'calc_chisqr')


@pytest.mark.parametrize("usesys", [True, False])
def test_stats_calc_chisqr_chi2_nostat(usesys):
    """chi-quare statistic with no staterror fails"""

    data, model = setup_single(False, usesys)
    statobj = Chi2()
    with pytest.raises(StatErr):
        statobj.calc_chisqr(data, model)


@pytest.mark.parametrize("usesys", [True, False])
def test_stats_calc_stat_chi2_nostat(usesys):
    """chi-quare statistic with no staterror fails"""

    data, model = setup_single(False, usesys)
    statobj = Chi2()
    with pytest.raises(StatErr):
        statobj.calc_stat_from_data(data, model)


@pytest.mark.parametrize("usestat,usesys", [
    (True, True),
    (True, False),
    (False, True),
    (False, False),
])
def test_stats_calc_stat_wstat_nobg(usestat, usesys):
    """wstat statistic fails with no background"""

    statobj = WStat()

    data, model = setup_single(usestat, usesys)
    with pytest.raises(StatErr):
        statobj.calc_stat_from_data(data, model)

    data, model = setup_single_pha(False, False, background=False)
    with pytest.raises(StatErr):
        statobj.calc_stat_from_data(data, model)


def test_stats_calc_stat_wstat_diffbins():
    """wstat statistic fails when src/bg bin sizes do not match"""

    statobj = WStat()

    data, model = setup_single_pha(True, False, background=True)

    # Tweak data to have one-less bin than the background
    d = data.datasets[0]
    d.channel = d.channel[:-1]
    d.counts = d.channel[:-1]
    for attr in ['staterror', 'syserror', 'grouping', 'quality',
                 'backscal']:
        val = getattr(d, attr)
        if val is not None:
            try:
                setattr(d, attr, val[:-1])
            except TypeError:
                # assume a scalar, so leave be
                pass

    # There is no Sherpa error for this, which seems surprising
    with pytest.raises(TypeError):
        statobj.calc_stat_from_data(data, model)


# Numeric answers calculated using CIAO 4.8 (sherpa 4.8.0)
# on a 64-bit Linux machine.
#
exp_delta = np.asarray([0.1, 0.8, 1.1, 0.6])
exp_lsq = exp_delta * exp_delta

exp_chi2_tt = np.asarray([0.00060009, 0.24860162, 0.27153028, 0.06166073])
exp_chi2_tf = np.asarray([0.000625, 0.25, 0.27437642, 0.0625])

exp_gehrels_ft = np.asarray([0.00037075, 0.08296552, 0.11425155, 0.02887336])
exp_gehrels_ff = np.asarray([0.00038011, 0.08312068, 0.11475241, 0.02905606])

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
    (Chi2Gehrels, True, False, exp_chi2_tf),
    (Chi2Gehrels, False, True, exp_gehrels_ft),
    (Chi2Gehrels, False, False, exp_gehrels_ff),
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

    data, model = setup_single(usestat, usesys)
    statobj = stat()
    answer = statobj.calc_chisqr(data, model)
    assert_almost_equal(answer, expected)


exp1dint_delta = np.asarray([1.43333333, 0.02666667, 2.67333333, 4.72])
exp1dint_lsq = exp1dint_delta * exp1dint_delta

exp1dint_chi2_tt = np.asarray([0.123284728, 0.000276224018,
                               1.60375904, 3.81583996])
exp1dint_chi2_tf = np.asarray([0.128402778, 0.000277777778,
                               1.62056941, 3.86777778])

exp1dint_gehrels_ft = np.asarray([0.0761679608, 0.0000921839121,
                                  0.674812245, 1.78681175])
exp1dint_gehrels_ff = np.asarray([0.0780910272, 0.0000923563159,
                                  0.677770505, 1.79811827])

exp1dint_dvar_ft = np.asarray([0.121104527, 0.000294529122,
                               1.64434909, 3.92335869])
exp1dint_dvar_ff = np.asarray([0.126039536, 0.000296296296,
                               1.66202584, 3.97828571])

exp1dint_cvar_ft = np.asarray([0.262910838, 0.0000992561989,
                               0.993119463, 3.08206519])
exp1dint_cvar_ff = np.asarray([0.287334887, 0.0000994560995,
                               0.999540016, 3.11586014])

exp1dint_mvar_xt = np.asarray([0.111669408, 0.000291311631,
                               4.27207048, 2.14248346])
exp1dint_mvar_xf = np.asarray([0.115852130, 0.000293040293,
                               4.39346995, 2.15875969])


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (Chi2, True, True, exp1dint_chi2_tt),
    (Chi2, True, False, exp1dint_chi2_tf),
    (LeastSq, True, True, exp1dint_lsq),
    (LeastSq, True, False, exp1dint_lsq),
    (LeastSq, False, True, exp1dint_lsq),
    (LeastSq, False, False, exp1dint_lsq),
    (Chi2Gehrels, True, True, exp1dint_chi2_tt),
    (Chi2Gehrels, True, False, exp1dint_chi2_tf),
    (Chi2Gehrels, False, True, exp1dint_gehrels_ft),
    (Chi2Gehrels, False, False, exp1dint_gehrels_ff),
    (Chi2DataVar, True, True, exp1dint_chi2_tt),
    (Chi2DataVar, True, False, exp1dint_chi2_tf),
    (Chi2DataVar, False, True, exp1dint_dvar_ft),
    (Chi2DataVar, False, False, exp1dint_dvar_ff),
    (Chi2XspecVar, True, True, exp1dint_chi2_tt),
    (Chi2XspecVar, True, False, exp1dint_chi2_tf),
    (Chi2XspecVar, False, True, exp1dint_dvar_ft),
    (Chi2XspecVar, False, False, exp1dint_dvar_ff),
    (Chi2ConstVar, True, True, exp1dint_chi2_tt),
    (Chi2ConstVar, True, False, exp1dint_chi2_tf),
    (Chi2ConstVar, False, True, exp1dint_cvar_ft),
    (Chi2ConstVar, False, False, exp1dint_cvar_ff),
    (Chi2ModVar, True, True, exp1dint_mvar_xt),
    (Chi2ModVar, True, False, exp1dint_mvar_xf),
    (Chi2ModVar, False, True, exp1dint_mvar_xt),
    (Chi2ModVar, False, False, exp1dint_mvar_xf)
])
def test_stats_calc_chisqr_1dint(stat, usestat, usesys, expected):
    """chi-quare statistics calculates expected values"""

    data, model = setup_single_1dint(usestat, usesys)
    statobj = stat()
    answer = statobj.calc_chisqr(data, model)
    assert_almost_equal(answer, expected)


@pytest.mark.parametrize("stat,usestat,usesys,expected1", [
    (Chi2, True, True, exp_chi2_tt),
    (Chi2, True, False, exp_chi2_tf),
    (LeastSq, True, True, exp_lsq),
    (LeastSq, True, False, exp_lsq),
    (LeastSq, False, True, exp_lsq),
    (LeastSq, False, False, exp_lsq),
    (Chi2Gehrels, True, True, exp_chi2_tt),
    (Chi2Gehrels, True, False, exp_chi2_tf),
    (Chi2Gehrels, False, True, exp_gehrels_ft),
    (Chi2Gehrels, False, False, exp_gehrels_ff),
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

    Here multiple datasets are used.
    """

    data, model = setup_multiple(usestat, usesys)
    statobj = stat()
    answer = statobj.calc_chisqr(data, model)

    data2 = data.datasets[1]
    expected = np.concatenate((expected1,
                               expected1[data2.mask]))
    assert_almost_equal(answer, expected)


@pytest.mark.parametrize("stat,usestat,usesys,expected1", [
    (Chi2, True, True, exp1dint_chi2_tt),
    (Chi2, True, False, exp1dint_chi2_tf),
    (LeastSq, True, True, exp1dint_lsq),
    (LeastSq, True, False, exp1dint_lsq),
    (LeastSq, False, True, exp1dint_lsq),
    (LeastSq, False, False, exp1dint_lsq),
    (Chi2Gehrels, True, True, exp1dint_chi2_tt),
    (Chi2Gehrels, True, False, exp1dint_chi2_tf),
    (Chi2Gehrels, False, True, exp1dint_gehrels_ft),
    (Chi2Gehrels, False, False, exp1dint_gehrels_ff),
    (Chi2DataVar, True, True, exp1dint_chi2_tt),
    (Chi2DataVar, True, False, exp1dint_chi2_tf),
    (Chi2DataVar, False, True, exp1dint_dvar_ft),
    (Chi2DataVar, False, False, exp1dint_dvar_ff),
    (Chi2XspecVar, True, True, exp1dint_chi2_tt),
    (Chi2XspecVar, True, False, exp1dint_chi2_tf),
    (Chi2XspecVar, False, True, exp1dint_dvar_ft),
    (Chi2XspecVar, False, False, exp1dint_dvar_ff),
    (Chi2ConstVar, True, True, exp1dint_chi2_tt),
    (Chi2ConstVar, True, False, exp1dint_chi2_tf),
    # See test_stats_calc_chisqr_1dint_multi_cvar for the next two
    # (Chi2ConstVar, False, True, exp1dint_cvar_ft),
    # (Chi2ConstVar, False, False, exp1dint_cvar_ff),
    (Chi2ModVar, True, True, exp1dint_mvar_xt),
    (Chi2ModVar, True, False, exp1dint_mvar_xf),
    (Chi2ModVar, False, True, exp1dint_mvar_xt),
    (Chi2ModVar, False, False, exp1dint_mvar_xf)
])
def test_stats_calc_chisqr_1dint_multi(stat, usestat, usesys, expected1):
    """chi-quare statistics calculates expected values

    Here multiple datasets are used. To save time, the first
    dataset is re-used, excluding the third point.
    """

    data, model = setup_multiple_1dint(usestat, usesys)
    statobj = stat()
    answer = statobj.calc_chisqr(data, model)

    data2 = data.datasets[1]
    expected = np.concatenate((expected1,
                               expected1[data2.mask]))
    assert_almost_equal(answer, expected)


exp_cvar_2_ft = np.asarray([0.001141, 0.07887213, 0.04401839])
exp_cvar_2_ff = np.asarray([0.00123457, 0.07901235, 0.04444444])


@pytest.mark.parametrize("usestat,usesys,expected", [
    (False, True, np.concatenate((exp_cvar_ft, exp_cvar_2_ft))),
    (False, False, np.concatenate((exp_cvar_ff, exp_cvar_2_ff))),
])
def test_stats_calc_chisqr_multi_cvar(usestat, usesys, expected):
    """chi-quare statistics calculates expected values: constvar

    Here multiple datasets are used. To save time, the first
    dataset is re-used, excluding the third point. This means that
    the constvar results are a little-bit different to the other
    statistics.
    """

    data, model = setup_multiple(usestat, usesys)
    statobj = Chi2ConstVar()
    answer = statobj.calc_chisqr(data, model)
    assert_almost_equal(answer, expected)


@pytest.mark.parametrize("usestat,usesys,expected", [
    (False, True, np.concatenate((exp1dint_cvar_ft,
                                  np.asarray([0.234412563, 0.0000876356984,
                                              2.72405360])))),
    (False, False, np.concatenate((exp1dint_cvar_ff,
                                   np.asarray([0.253635117, 0.0000877914952,
                                               2.75041975])))),
])
def test_stats_calc_chisqr_1dint_multi_cvar(usestat, usesys, expected):
    """chi-quare statistics calculates expected values: constvar

    Here multiple datasets are used. To save time, the first
    dataset is re-used, excluding the third point. This means that
    the constvar results are a little-bit different to the other
    statistics.
    """

    data, model = setup_multiple_1dint(usestat, usesys)
    statobj = Chi2ConstVar()
    answer = statobj.calc_chisqr(data, model)
    assert_almost_equal(answer, expected)


# Numeric answers calculated using CIAO 4.8 (sherpa 4.8.0)
# on a 64-bit Linux machine.
#
stat_lsq = exp_lsq.sum()

stat_chi2_tt = exp_chi2_tt.sum()
stat_chi2_tf = exp_chi2_tf.sum()

stat_gehrels_ft = exp_gehrels_ft.sum()
stat_gehrels_ff = exp_gehrels_ff.sum()

stat_dvar_ft = exp_dvar_ft.sum()
stat_dvar_ff = exp_dvar_ff.sum()

stat_cvar_ft = exp_cvar_ft.sum()
stat_cvar_ff = exp_cvar_ff.sum()

stat_mvar_xt = exp_mvar_xt.sum()
stat_mvar_xf = exp_mvar_xf.sum()

stat_cash = -69.2032919676
stat_cstat = 0.630015576282


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (Chi2, True, True, stat_chi2_tt),
    (Chi2, True, False, stat_chi2_tf),
    (LeastSq, True, True, stat_lsq),
    (LeastSq, True, False, stat_lsq),
    (LeastSq, False, True, stat_lsq),
    (LeastSq, False, False, stat_lsq),
    (Chi2Gehrels, True, True, stat_chi2_tt),
    (Chi2Gehrels, True, False, stat_chi2_tf),
    (Chi2Gehrels, False, True, stat_gehrels_ft),
    (Chi2Gehrels, False, False, stat_gehrels_ff),
    (Chi2DataVar, True, True, stat_chi2_tt),
    (Chi2DataVar, True, False, stat_chi2_tf),
    (Chi2DataVar, False, True, stat_dvar_ft),
    (Chi2DataVar, False, False, stat_dvar_ff),
    (Chi2XspecVar, True, True, stat_chi2_tt),
    (Chi2XspecVar, True, False, stat_chi2_tf),
    (Chi2XspecVar, False, True, stat_dvar_ft),
    (Chi2XspecVar, False, False, stat_dvar_ff),
    (Chi2ConstVar, True, True, stat_chi2_tt),
    (Chi2ConstVar, True, False, stat_chi2_tf),
    (Chi2ConstVar, False, True, stat_cvar_ft),
    (Chi2ConstVar, False, False, stat_cvar_ff),
    (Chi2ModVar, True, True, stat_mvar_xt),
    (Chi2ModVar, True, False, stat_mvar_xf),
    (Chi2ModVar, False, True, stat_mvar_xt),
    (Chi2ModVar, False, False, stat_mvar_xf),

    (Cash, True, True, stat_cash),
    (Cash, True, False, stat_cash),
    (Cash, False, True, stat_cash),
    (Cash, False, False, stat_cash),
    (CStat, True, True, stat_cstat),
    (CStat, True, False, stat_cstat),
    (CStat, False, True, stat_cstat),
    (CStat, False, False, stat_cstat),

])
def test_stats_calc_stat(stat, usestat, usesys, expected):
    """statistic calculates expected values: single dataset"""

    data, model = setup_single(usestat, usesys)
    statobj = stat()
    # do not check fvec
    answer, _ = statobj.calc_stat_from_data(data, model)
    assert_almost_equal(answer, expected)


delta_lsq = exp_lsq[2]

delta_chi2_tt = exp_chi2_tt[2]
delta_chi2_tf = exp_chi2_tf[2]

delta_gehrels_ft = exp_gehrels_ft[2]
delta_gehrels_ff = exp_gehrels_ff[2]

delta_dvar_ft = exp_dvar_ft[2]
delta_dvar_ff = exp_dvar_ff[2]

# These do not follow the "expected" rule
delta_cvar_ft = exp_cvar_ft.sum() - exp_cvar_2_ft.sum()
delta_cvar_ff = exp_cvar_ff.sum() - exp_cvar_2_ff.sum()

delta_mvar_xt = exp_mvar_xt[2]
delta_mvar_xf = exp_mvar_xf[2]

delta_cash = -3.60309696433
delta_cstat = 0.340992230887


@pytest.mark.parametrize("stat,usestat,usesys,expected1,delta", [
    (Chi2, True, True, stat_chi2_tt, delta_chi2_tt),
    (Chi2, True, False, stat_chi2_tf, delta_chi2_tf),
    (LeastSq, True, True, stat_lsq, delta_lsq),
    (LeastSq, True, False, stat_lsq, delta_lsq),
    (LeastSq, False, True, stat_lsq, delta_lsq),
    (LeastSq, False, False, stat_lsq, delta_lsq),
    (Chi2Gehrels, True, True, stat_chi2_tt, delta_chi2_tt),
    (Chi2Gehrels, True, False, stat_chi2_tf, delta_chi2_tf),
    (Chi2Gehrels, False, True, stat_gehrels_ft, delta_gehrels_ft),
    (Chi2Gehrels, False, False, stat_gehrels_ff, delta_gehrels_ff),
    (Chi2DataVar, True, True, stat_chi2_tt, delta_chi2_tt),
    (Chi2DataVar, True, False, stat_chi2_tf, delta_chi2_tf),
    (Chi2DataVar, False, True, stat_dvar_ft, delta_dvar_ft),
    (Chi2DataVar, False, False, stat_dvar_ff, delta_dvar_ff),
    (Chi2XspecVar, True, True, stat_chi2_tt, delta_chi2_tt),
    (Chi2XspecVar, True, False, stat_chi2_tf, delta_chi2_tf),
    (Chi2XspecVar, False, True, stat_dvar_ft, delta_dvar_ft),
    (Chi2XspecVar, False, False, stat_dvar_ff, delta_dvar_ff),
    (Chi2ConstVar, True, True, stat_chi2_tt, delta_chi2_tt),
    (Chi2ConstVar, True, False, stat_chi2_tf, delta_chi2_tf),
    # Unlike test_stats_calcchisqr_multi, the following two
    # can be included in this test
    (Chi2ConstVar, False, True, stat_cvar_ft, delta_cvar_ft),
    (Chi2ConstVar, False, False, stat_cvar_ff, delta_cvar_ff),
    (Chi2ModVar, True, True, stat_mvar_xt, delta_mvar_xt),
    (Chi2ModVar, True, False, stat_mvar_xf, delta_mvar_xf),
    (Chi2ModVar, False, True, stat_mvar_xt, delta_mvar_xt),
    (Chi2ModVar, False, False, stat_mvar_xf, delta_mvar_xf),

    (Cash, True, True, stat_cash, delta_cash),
    (Cash, True, False, stat_cash, delta_cash),
    (Cash, False, True, stat_cash, delta_cash),
    (Cash, False, False, stat_cash, delta_cash),
    (CStat, True, True, stat_cstat, delta_cstat),
    (CStat, True, False, stat_cstat, delta_cstat),
    (CStat, False, True, stat_cstat, delta_cstat),
    (CStat, False, False, stat_cstat, delta_cstat),

])
def test_stats_calc_stat_multi(stat, usestat, usesys, expected1, delta):
    """statistic calculates expected values: multiple dataset"""

    data, model = setup_multiple(usestat, usesys)
    statobj = stat()
    # do not check fvec
    answer, _ = statobj.calc_stat_from_data(data, model)

    # correct for the one missing bin in the second dataset
    expected = 2 * expected1 - delta
    assert_almost_equal(answer, expected)


# These values were created on a 64-bit Linux machine using
# CIAO 4.8
#
stat_pha_lsq = 51.82

stat_pha_chi2_tt = 1.30511684776
stat_pha_chi2_tf = 2.0728

# with background subtraction
stat_pha_lsq_bg = 53.9958239163

stat_pha_chi2_bg_tt = 1.36147410407
stat_pha_chi2_bg_tf = 2.15970543269

stat_pha_gehrels_ff = 1.37612618038
stat_pha_gehrels_ft = 0.989239848562

stat_pha_gehrels_bg_ff = 1.43234664248
stat_pha_gehrels_bg_ft = 1.03105762001

# TODO: should the background subtraction actually raise an errror
#       for the likelihood statistics?
#
stat_pha_cash = -299.771783393
stat_pha_cash_bg = -297.366618448

stat_pha_cstat = 1.75191290087
stat_pha_cstat_bg = 1.8269289994


# This is not as extensive as some of the earlier checks as the
# assumption is that if it works for these variants it should be
# okay with the others.
#
@pytest.mark.parametrize("stat,usestat,usesys,havebg,usebg,expected", [
    (LeastSq, True, True, False, False, stat_pha_lsq),
    (LeastSq, True, True, True, False, stat_pha_lsq),
    (LeastSq, True, True, True, True, stat_pha_lsq_bg),

    (LeastSq, False, False, False, False, stat_pha_lsq),
    (LeastSq, False, False, True, False, stat_pha_lsq),
    (LeastSq, False, False, True, True, stat_pha_lsq_bg),

    (Chi2, True, True, False, False, stat_pha_chi2_tt),
    (Chi2, True, True, True, False, stat_pha_chi2_tt),
    (Chi2, True, True, True, True, stat_pha_chi2_bg_tt),

    (Chi2, True, False, False, False, stat_pha_chi2_tf),
    (Chi2, True, False, True, False, stat_pha_chi2_tf),
    (Chi2, True, False, True, True, stat_pha_chi2_bg_tf),

    (Chi2Gehrels, True, True, False, False, stat_pha_chi2_tt),
    (Chi2Gehrels, True, True, True, False, stat_pha_chi2_tt),
    (Chi2Gehrels, True, True, True, True, stat_pha_chi2_bg_tt),

    (Chi2Gehrels, True, False, False, False, stat_pha_chi2_tf),
    (Chi2Gehrels, True, False, True, False, stat_pha_chi2_tf),
    (Chi2Gehrels, True, False, True, True, stat_pha_chi2_bg_tf),

    (Chi2Gehrels, False, False, False, False, stat_pha_gehrels_ff),
    (Chi2Gehrels, False, False, True, False, stat_pha_gehrels_ff),
    (Chi2Gehrels, False, False, True, True, stat_pha_gehrels_bg_ff),

    (Chi2Gehrels, False, True, False, False, stat_pha_gehrels_ft),
    (Chi2Gehrels, False, True, True, False, stat_pha_gehrels_ft),
    (Chi2Gehrels, False, True, True, True, stat_pha_gehrels_bg_ft),

    (Cash, True, True, False, False, stat_pha_cash),
    (Cash, False, False, False, False, stat_pha_cash),
    (Cash, False, False, True, True, stat_pha_cash_bg),

    (CStat, True, True, False, False, stat_pha_cstat),
    (CStat, False, False, False, False, stat_pha_cstat),
    (CStat, False, False, True, True, stat_pha_cstat_bg),

])
def test_stats_calc_stat_pha(stat, usestat, usesys,
                             havebg, usebg, expected):
    """statistic calculates expected values: single PHA dataset"""

    data, model = setup_single_pha(usestat, usesys, background=havebg)
    if usebg:
        data.datasets[0].subtract()

    statobj = stat()
    # do not check fvec
    answer, _ = statobj.calc_stat_from_data(data, model)
    assert_almost_equal(answer, expected)
