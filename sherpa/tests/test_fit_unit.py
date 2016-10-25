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

# Unit tests for fit to support the statistic redesign (part of
# investigating #248 and #227). It is possible that many of these
# tests can be removed after these issues have been addressed (i.e.
# the functionality should be tested elsewhere).
#
# It includes very-basic tests of the fit/error calculations
# (at present fit only).
#
# Notes
#
# Since the multiple-dataset setup code returns Fit objects for the
# individual data as well as the combined set, they could be used
# for the "single dataset" test. This might clean up the test code
# a bit. However, for now leave as is. The data values and model
# parameter values should also be adjusted to provide a better test of
# the "q value" calculated for chi-square like statistics (at present
# the value is often << 1e-7, which is the limit for
# assert_almost_equal).
#
# Some of these tests may be superfluous (e.g. testing both the single
# and multiple dataset handling for calc_chisqr with wstat), but given
# the current design, where the functionality is embedded within several
# different classes, it is not clear how much we can rely on testing
# separate behavior and then assuming they work together well. This is
# the reason for some of the suggestions in the possible improvements
# section below: really we should be able to rely on the filter/group
# behavior of the data objects, so that we don't need to check as
# many combinations as we do), but the current logic (and API of
# some of the classes) means that they need to be tested here.
#
# As these tests are really a combination of the behavior of the
# Fit, Data, Model, and Stats classes, perhaps the file name shoule be
# changed.
#
# Possible improvements:
#
# Ideally would test grouping/filtering for non-wstat statistics, to check
# that, when passing through fit, the mapping gets handled correctly.
#
# Similarly, add tests of DataPHA objects (with and without background
# subtraction) for non-wstat statistics.
#
# The existing filter tests are very simple, in that the filtered data is
# a single, contiguous chunk. A test like:
#     notice(None, None)
#     ignore(a, b)
# could be used, where a >> xlo and b << xhi (i.e. a-b is a range
# of one-or-more bins within the dataset that does not extend to
# the edges).
#
# To dos:
#
# *) wstat and rstat/qval
#
# Shouldn't WStat report both a rstat and qval?
#
# *) weights (and chi-square stistices)
#
# The tests do not include any ability to change the weights (e.g. for
# the chi-square statistic). I (DougBurke) may have missed it, but it's
# not clear how the Data class handles weights, so is this functionality
# that the stats class has but isn't available to the Fit class?
#
# *) chi2datavar and chi2xspecvar
#
# There are several places where it is not clear to me (DougBurke) why
# the Chi2XspecVar result agrees or disagrees with the Chi2DataVar result.
# Note that the stats test should be probing those cases, as here they
# are just representative tests.
#

import pytest
import numpy as np

from numpy.testing import assert_almost_equal

from sherpa.fit import Fit, StatInfoResults
from sherpa.data import Data1D, DataSimulFit
from sherpa.astro.data import DataPHA
from sherpa.models.model import SimulFitModel
from sherpa.models.basic import Const1D, Gauss1D, Polynom1D, StepLo1D
from sherpa.utils.err import DataErr, EstErr, FitErr, StatErr

from sherpa.stats import LeastSq, Chi2, Chi2Gehrels, Chi2DataVar, \
    Chi2ConstVar, Chi2ModVar, Chi2XspecVar, Cash, CStat, WStat, UserStat

from sherpa.optmethods import LevMar, NelderMead
from sherpa.estmethods import Covariance, Confidence


def setup_stat_single(stat, usestat, usesys):
    """Set up a single dataset with the given statistic.

    A sherpa.data.Data1D instance is used.

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
        The statistic object to use
    usestat, usesys : bool
        Should statistical or systematic errors be used?

    Returns
    -------
    fit : sherpa.fit.Fit instance
    """

    # Ensure there's a bin with a dependent-axis value of 0
    x = np.asarray([10, 20, 30, 35, 42])
    y = np.asarray([12, 15, 10, 17, 0])
    if usestat:
        ystat = np.asarray([1.2, 1.3, 2.1, 0.8, 0.4])
    else:
        ystat = None
    if usesys:
        ysys = np.asarray([1.3, 1.2, 2.4, 0.2, 1.4])
    else:
        ysys = None

    data = Data1D("test", x, y, staterror=ystat, syserror=ysys)

    # Use a model which varies with position, to check that
    # the independent axis is being passed through correctly.
    # A simple model (single component) is sufficient, but
    # this is done to better-handle the last y value (of 0),
    # to make sure the data and model values are not too far
    # apart.
    #
    mdl = (Polynom1D("poly") * StepLo1D("step")) + Const1D("bg")
    pars = {p.fullname: p for p in mdl.pars}
    pars['poly.c0'].val = 8
    pars['poly.c1'].val = 0.15
    pars['step.xcut'].val = 40
    pars['bg.c0'].val = 2

    return Fit(data, mdl, stat=stat)


def setup_stat_multiple(stat, usestat, usesys):
    """Set up multiple datasets with the given statistic.

    The data is stored in sherpa.data.Data1D instances.

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
        The statistic object to use
    usestat, usesys : bool
        Should statistical or systematic errors be used?

    Returns
    -------
    fit, fits : sherpa.fit.Fit instance, list of sherpa.fit.Fit instances
        The first fit is for the combined dataset, and then fit objects
        for the individual data sets are returned in a list.

    """

    x1 = np.asarray([10, 20, 35])
    y1 = np.asarray([12, 15, 17])
    x2 = np.asarray([-400, -350, -325, -300])
    y2 = np.asarray([150, 30, 180, 98])
    x3 = np.asarray([100, 200, 300, 400, 500])
    y3 = np.asarray([2002, 2009, 2021, 2033, 2038])
    if usestat:
        ystat1 = np.asarray([1.2, 1.3, 0.8])
        ystat2 = np.asarray([30, 40, 25, 12])
        ystat3 = np.asarray([0.2, 0.3, 0.35, 0.4, 0.5])
    else:
        ystat1 = None
        ystat2 = None
        ystat3 = None
    if usesys:
        ysys1 = np.asarray([1.3, 1.2, 0.2])
        ysys2 = np.asarray([5, 6, 4, 3])
        ysys3 = np.asarray([1, 2, 2, 1, 2])
    else:
        ysys1 = None
        ysys2 = None
        ysys3 = None

    data1 = Data1D("test1", x1, y1, staterror=ystat1, syserror=ysys1)
    data2 = Data1D("test2", x2, y2, staterror=ystat2, syserror=ysys2)
    data3 = Data1D("test3", x3, y3, staterror=ystat3, syserror=ysys3)

    data = DataSimulFit("multidata", (data1, data2, data3))

    mdl1 = Const1D("mdl1")
    mdl1.c0 = 14
    mdl2 = Const1D("mdl2")
    mdl2.c0 = 100
    mdl3 = Polynom1D("mdl3")
    mdl3.c0 = 1990
    mdl3.c1 = 0.1

    mdl = SimulFitModel("multimodel", (mdl1, mdl2, mdl3))

    fit = Fit(data, mdl, stat=stat)
    fit1 = Fit(data1, mdl1, stat=stat)
    fit2 = Fit(data2, mdl2, stat=stat)
    fit3 = Fit(data3, mdl3, stat=stat)

    return fit, [fit1, fit2, fit3]


def setup_pha_single(scalar, usestat, usesys, flo, fhi, stat=None):
    """Set up a single PHA dataset.

    A sherpa.data.DataPHA instance is used, with an associated
    background data set.

    Parameters
    ----------
    scalar : bool
        Should a scalar BACKSCAL value be used for the background
        region (True) or an array (False)? The source always uses
        a scalar BACKSCAL.
    usestat, usesys : bool
        Should statistical or systematic errors be used?
    flo, fhi : None or number
        The channel filters (used in a notice call) that is applied
        to the source only.
    stat : sherpa.stats.Stat instance, optional
        The statistic object to use. If not set uses a WStat object.

    Returns
    -------
    fit : sherpa.fit.Fit instance
    """

    # ran numpy.random.poisson([10, 10, 10, 10, 10]) to get
    # the src_counts values below, and with values of 2 for
    # the background counts (just to get poisson-distributed
    # data; it doesn't actually matter how it is distributed
    # here as just evaluating at these values and not fitting).
    #
    # The channel starts at 1 just to follow the expected PHA
    # behavior, but it should not matter here.
    #
    channels = np.arange(1, 6, dtype=np.int)
    src_counts = np.asarray([14, 15, 11, 3, 8], dtype=np.int)
    bg_counts = np.asarray([2, 0, 2, 3, 4], dtype=np.int)

    # TODO: can add a grouping flag; if given use a larger set of bins and set
    # up grouping

    src = DataPHA('tst', channels, src_counts,
                  exposure=100.0, backscal=0.01)

    if scalar:
        backscal = 0.03
    else:
        backscal = np.asarray([0.03, 0.02, 0.04, 0.03, 0.03])

    bg = DataPHA('tstbg', channels, bg_counts,
                 exposure=200.0, backscal=backscal)

    if usestat:
        src.staterror = np.sqrt(src.counts)
        bg.staterror = np.sqrt(bg.counts)

    if usesys:
        src.syserror = 0.1 * src.counts
        bg.syserror = 0.1 * bg.counts

    src.set_background(bg)

    # Filtering is *only* applied to the source
    src.notice(lo=flo, hi=fhi)

    # Pick something that varies with the channel number, to
    # make sure that filtering is handled correctly.
    mdl = Polynom1D('mdl')
    mdl.c0 = 8
    mdl.c1 = 1.2

    # For now restrict to the default statistic values
    if stat is None:
        statobj = WStat()
    else:
        statobj = stat

    return Fit(src, mdl, stat=statobj)


def setup_pha_multiple(flo, fhi):
    """Set up multiple PHA datasets.

    A sherpa.data.DataPHA instance is used, with an associated
    background data set.

    Parameters
    ----------
    flo, fhi : None or number
        The channel filters (used in a notice call). These filters are
        applied to source only, both, and background only.

    Returns
    -------
    fit, fits : sherpa.fit.Fit instance, list of sherpa.fit.Fit instances
        The fit object for the combined data set and a list of those
        for the individual data sets.
    """

    # Using three data sets means that can have scalar and array
    # backscale values.
    #
    x1 = np.arange(1, 10)
    y1 = np.asarray([2, 3, 2, 5, 3, 0, 0, 2, 1])
    b1 = np.asarray([4, 2, 1, 2, 1, 0, 1, 3, 2])

    x2 = np.arange(1, 6)
    y2 = np.asarray([10, 12, 14, 12, 14])
    b2 = np.asarray([2, 4, 3, 2, 1])

    x3 = np.arange(1, 12)
    y3 = np.asarray([0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0])
    b3 = np.asarray([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])

    e1 = 100.0
    be1 = 125.0
    bscal1 = 0.012
    bbscal1 = np.asarray([0.024, 0.024, 0.023, 0.025, 0.03, 0.025,
                          0.04, 0.035, 0.031])

    e2 = 1000.0
    be2 = 1000.0
    bscal2 = np.asarray([0.12, 0.13, 0.12, 0.09, 0.08])
    bbscal2 = 0.24

    e3 = 800.0
    be3 = 700.0
    bscal3 = 0.0145 * np.ones(y3.size)
    bbscal3 = 0.0232 * np.ones(y3.size)
    bscal3[0] = 0.008
    bbscal3[-1] = 0.012

    src1 = DataPHA('p1', x1, y1, exposure=e1, backscal=bscal1)
    src2 = DataPHA('p2', x2, y2, exposure=e2, backscal=bscal2)
    src3 = DataPHA('p3', x3, y3, exposure=e3, backscal=bscal3)

    bg1 = DataPHA('b1', x1, b1, exposure=be1, backscal=bbscal1)
    bg2 = DataPHA('b2', x2, b2, exposure=be2, backscal=bbscal2)
    bg3 = DataPHA('b3', x3, b3, exposure=be3, backscal=bbscal3)

    src1.set_background(bg1)
    src2.set_background(bg2)
    src3.set_background(bg3)

    # Apply filtering to
    #  1: source only
    #  2: source and background (same filter)
    #  3: background only
    #
    # The background filters should be ignored
    #
    src1.notice(lo=flo, hi=fhi)
    src2.notice(lo=flo, hi=fhi)
    bg2.notice(lo=flo, hi=fhi)
    bg3.notice(lo=flo, hi=fhi)

    # Model components
    m1 = Gauss1D('m1')
    m1.fwhm = 0.8
    m1.pos = 4.05
    m1.ampl = 4.5

    m2 = Const1D('m2')
    m2.c0 = 12

    m3 = Polynom1D('m3')
    m3.c0 = 1.1
    m3.c1 = -0.05

    # Models for the datasets; at present no component is re-used,
    # since support for that functionality should be tested elsewhere,
    # and it should not be directly relevant to the statistics code
    # being tested.
    mexpr1 = m1
    mexpr2 = m2
    mexpr3 = m3

    data = DataSimulFit("multi", (src1, src2, src3))
    model = SimulFitModel("multi", (mexpr1, mexpr2, mexpr3))

    statobj = WStat()

    fit = Fit(data, model, stat=statobj)

    fit1 = Fit(src1, mexpr1, stat=statobj)
    fit2 = Fit(src2, mexpr2, stat=statobj)
    fit3 = Fit(src3, mexpr3, stat=statobj)

    return fit, [fit1, fit2, fit3]


@pytest.mark.parametrize("stat", [Cash, CStat, WStat])
def test_fit_raises_error_on_bgsubtraction(stat):
    """Check that fits using likelihood stats fail when bkg subtracted."""

    # It is not worth looping through all the combinations to the setup
    # routine here.
    statobj = stat()
    fit = setup_pha_single(False, False, False, None, None, stat=statobj)
    fit.data.subtract()

    with pytest.raises(FitErr) as excinfo:
        fit.fit()

    emsg = '{} statistics cannot be used with '.format(statobj.name) + \
           'background subtracted data'
    assert str(excinfo.value) == emsg

expected_leastsq = 31.5625

expected_chi2_tf = 36.9174680011
expected_chi2_tt = 9.73944684193

expected_gehrels_ff = 2.5415403715
expected_gehrels_ft = 1.85309576866

# The datavar_ft result was calculated using Sherpa 4.8.2 since this
# combination raised an error in earlier releases (see PR #153)
#
expected_chi2_datavar_ff = 6.49264705882
expected_chi2_datavar_ft = 3.76700948927

expected_chi2_constvar_ff = 2.9224537037
expected_chi2_constvar_ft = 2.1656375594

expected_chi2_xspecvar_ff = expected_chi2_datavar_ff
expected_chi2_xspecvar_ft = 3.07754451409

expected_mod_statonly = 3.9268028344
expected_mod_stat_sys = 2.50586379978

expected_cash = -169.183485665
expected_cstat = 6.07673552072


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (LeastSq, False, False, expected_leastsq),
    (LeastSq, True, False, expected_leastsq),
    (LeastSq, False, True, expected_leastsq),
    (LeastSq, True, True, expected_leastsq),

    # Note that Chi2 needs a staterror column but the other Chi2xxx don't
    (Chi2, True, False, expected_chi2_tf),
    (Chi2, True, True, expected_chi2_tt),

    (Chi2Gehrels, False, False, expected_gehrels_ff),
    (Chi2Gehrels, True, False, expected_chi2_tf),
    (Chi2Gehrels, False, True, expected_gehrels_ft),
    (Chi2Gehrels, True, True, expected_chi2_tt),

    (Chi2DataVar, False, False, expected_chi2_datavar_ff),
    (Chi2DataVar, True, False, expected_chi2_tf),
    (Chi2DataVar, False, True, expected_chi2_datavar_ft),
    (Chi2DataVar, True, True, expected_chi2_tt),

    (Chi2ConstVar, False, False, expected_chi2_constvar_ff),
    (Chi2ConstVar, True, False, expected_chi2_tf),
    (Chi2ConstVar, False, True, expected_chi2_constvar_ft),
    (Chi2ConstVar, True, True, expected_chi2_tt),

    # Chi2ModVar does not behave similarly to the other
    # Chi2 statistics; it is assumed that the following
    # is correct.
    (Chi2ModVar, False, False, expected_mod_statonly),
    (Chi2ModVar, True, False, expected_mod_statonly),
    (Chi2ModVar, False, True, expected_mod_stat_sys),
    (Chi2ModVar, True, True, expected_mod_stat_sys),

    (Chi2XspecVar, False, False, expected_chi2_xspecvar_ff),
    (Chi2XspecVar, True, False, expected_chi2_tf),
    (Chi2XspecVar, False, True, expected_chi2_xspecvar_ft),
    (Chi2XspecVar, True, True, expected_chi2_tt),

    (Cash, False, False, expected_cash),
    (Cash, True, False, expected_cash),
    (Cash, False, True, expected_cash),
    (Cash, True, True, expected_cash),

    (CStat, False, False, expected_cstat),
    (CStat, True, False, expected_cstat),
    (CStat, False, True, expected_cstat),
    (CStat, True, True, expected_cstat)
])
def test_fit_calc_stat_single(stat, usestat, usesys, expected):
    """Test the results from the calc_stat method of Fit for 1 data set.

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    # For now restrict the tests to using the default statistic values
    #
    fit = setup_stat_single(stat(), usestat, usesys)
    assert_almost_equal(fit.calc_stat(), expected)


# Unfortunately the parameter choces mean that the qval values,
# when calculated, are often very small.
#
q_chi2_tf = 0.0  # actually ~ 1e-9 but check os to 7 dp
q_chi2_tt = 0.0018035516

q_mod_statonly = 0.0475222107
q_mod_sysonly = 0.1134232786

q_chi2_datavar_ff = 0.0108321566
q_chi2_datavar_ft = 0.0522730121

q_chi2_xspecvar_ff = q_chi2_datavar_ff
q_chi2_xspecvar_ft = 0.07938028609

q_cstat = 0.0136973660


@pytest.mark.parametrize("stat,usestat,usesys,expected,qval", [
    (LeastSq, False, False, expected_leastsq, None),
    (LeastSq, True, False, expected_leastsq, None),
    (LeastSq, False, True, expected_leastsq, None),
    (LeastSq, True, True, expected_leastsq, None),

    (Chi2, True, False, expected_chi2_tf, q_chi2_tf),
    (Chi2, True, True, expected_chi2_tt, q_chi2_tt),

    (Chi2Gehrels, False, False, expected_gehrels_ff, 0.1108865604),
    (Chi2Gehrels, True, False, expected_chi2_tf, q_chi2_tf),
    (Chi2Gehrels, False, True, expected_gehrels_ft, 0.1734237670),
    (Chi2Gehrels, True, True, expected_chi2_tt, q_chi2_tt),

    (Chi2DataVar, False, False, expected_chi2_datavar_ff, q_chi2_datavar_ff),
    (Chi2DataVar, True, False, expected_chi2_tf, q_chi2_tf),
    (Chi2DataVar, False, True, expected_chi2_datavar_ft, q_chi2_datavar_ft),
    (Chi2DataVar, True, True, expected_chi2_tt, q_chi2_tt),

    (Chi2ConstVar, False, False, expected_chi2_constvar_ff, 0.0873549369),
    (Chi2ConstVar, True, False, expected_chi2_tf, q_chi2_tf),
    (Chi2ConstVar, False, True, expected_chi2_constvar_ft, 0.1411260797),
    (Chi2ConstVar, True, True, expected_chi2_tt, q_chi2_tt),

    (Chi2ModVar, False, False, expected_mod_statonly, q_mod_statonly),
    (Chi2ModVar, True, False, expected_mod_statonly, q_mod_statonly),
    (Chi2ModVar, False, True, expected_mod_stat_sys, q_mod_sysonly),
    (Chi2ModVar, True, True, expected_mod_stat_sys, q_mod_sysonly),

    (Chi2XspecVar, False, False, expected_chi2_xspecvar_ff,
     q_chi2_xspecvar_ff),
    (Chi2XspecVar, True, False, expected_chi2_tf, q_chi2_tf),
    (Chi2XspecVar, False, True, expected_chi2_xspecvar_ft,
     q_chi2_xspecvar_ft),
    (Chi2XspecVar, True, True, expected_chi2_tt, q_chi2_tt),

    (Cash, False, False, expected_cash, None),
    (Cash, True, False, expected_cash, None),
    (Cash, False, True, expected_cash, None),
    (Cash, True, True, expected_cash, None),

    (CStat, False, False, expected_cstat, q_cstat),
    (CStat, True, False, expected_cstat, q_cstat),
    (CStat, False, True, expected_cstat, q_cstat),
    (CStat, True, True, expected_cstat, q_cstat)

])
def test_fit_calc_stat_info_single(stat, usestat, usesys, expected, qval):
    """See test_fit_calc_stat_single for the origin of the
    expected values.
    """

    statobj = stat()
    fit = setup_stat_single(statobj, usestat, usesys)
    ans = fit.calc_stat_info()

    assert isinstance(ans, StatInfoResults)
    assert ans.name == ""
    assert ans.ids is None
    assert ans.bkg_ids is None
    assert ans.numpoints == 5
    assert ans.dof == 1
    assert_almost_equal(ans.statval, expected)

    # The choice of ans.statname is interesting for chi-squared
    # cases, since it depends on whether usestat is True.
    #
    if statobj.name.startswith('chi2') and usestat:
        assert ans.statname == 'chi2'
    else:
        assert ans.statname == statobj.name

    if qval is None:
        assert ans.qval is None
        assert ans.rstat is None

    else:
        assert_almost_equal(ans.qval, qval)
        # as there's only 1 dof, statval == rstat
        assert_almost_equal(ans.rstat, ans.statval)


chisqr_leastsq = [0.25, 4.0, 20.25, 3.0625, 4.0]
chisqr_chi2_tf = [0.173611111, 2.366863905, 4.591836735,
                  4.785156250, 25.0]
chisqr_chi2_tt = [0.079872204, 1.277955272, 1.991150442,
                  4.503676471, 1.886792453]

chisqr_gehrels_ff = [0.011966630, 0.162026931, 1.106107770,
                     0.112690724, 1.148748316]
chisqr_gehrels_ft = [0.011071045, 0.153096839, 0.841385758,
                     0.112525101, 0.735017026]

chisqr_datavar_ff = [0.020833333, 0.266666667, 2.025000000,
                     0.180147059, 4.0]
chisqr_datavar_ft = [0.018261505, 0.243309002, 1.284898477,
                     0.179724178, 2.040816327]

chisqr_const_ff = [0.023148148, 0.370370370, 1.875,
                   0.283564815, 0.370370370]
chisqr_const_ft = [0.020016013, 0.326797386, 1.222826087,
                   0.282518450, 0.313479624]

chisqr_mod_xf = [0.021739130, 0.307692308, 1.396551724,
                 0.200819672, 2.0]
chisqr_mod_xt = [0.018953753, 0.277008310, 0.999506417,
                 0.200294310, 1.010101010]

chisqr_xspecvar_ff = chisqr_datavar_ff
chisqr_xspecvar_ft = [0.018261505, 0.243309002, 1.284898477,
                      0.179724178, 1.351351351]


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (LeastSq, False, False, chisqr_leastsq),
    (LeastSq, True, False, chisqr_leastsq),
    (LeastSq, False, True, chisqr_leastsq),
    (LeastSq, True, True, chisqr_leastsq),

    (Chi2, True, False, chisqr_chi2_tf),
    (Chi2, True, True, chisqr_chi2_tt),

    (Chi2Gehrels, False, False, chisqr_gehrels_ff),
    (Chi2Gehrels, True, False, chisqr_chi2_tf),
    (Chi2Gehrels, False, True, chisqr_gehrels_ft),
    (Chi2Gehrels, True, True, chisqr_chi2_tt),

    (Chi2DataVar, False, False, chisqr_datavar_ff),
    (Chi2DataVar, True, False, chisqr_chi2_tf),
    (Chi2DataVar, False, True, chisqr_datavar_ft),
    (Chi2DataVar, True, True, chisqr_chi2_tt),

    (Chi2ConstVar, False, False, chisqr_const_ff),
    (Chi2ConstVar, True, False, chisqr_chi2_tf),
    (Chi2ConstVar, False, True, chisqr_const_ft),
    (Chi2ConstVar, True, True, chisqr_chi2_tt),

    (Chi2ModVar, False, False, chisqr_mod_xf),
    (Chi2ModVar, True, False, chisqr_mod_xf),
    (Chi2ModVar, False, True, chisqr_mod_xt),
    (Chi2ModVar, True, True, chisqr_mod_xt),

    (Chi2XspecVar, False, False, chisqr_xspecvar_ff),
    (Chi2XspecVar, True, False, chisqr_chi2_tf),
    (Chi2XspecVar, False, True, chisqr_xspecvar_ft),
    (Chi2XspecVar, True, True, chisqr_chi2_tt),

    (Cash, False, False, None),
    (Cash, True, False, None),
    (Cash, False, True, None),
    (Cash, True, True, None),

    (CStat, False, False, None),
    (CStat, True, False, None),
    (CStat, False, True, None),
    (CStat, True, True, None)
])
def test_fit_calc_chisqr_single(stat, usestat, usesys, expected):
    """Test the results from the calc_chisqr method of Fit for 1 data set.

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    # For now restrict the tests to using the default statistic values
    #
    fit = setup_stat_single(stat(), usestat, usesys)
    ans = fit.calc_chisqr()
    if expected is None:
        assert ans is None
    else:
        assert_almost_equal(ans, expected)


expected_multi_leastsq = 13837.0
expected_multi_chi2 = 225.06442572689826
expected_multi_chi2_tt = 43.574115140290814

expected_multi_chi2_gehrels_ff = 159.83814050225632
expected_multi_chi2_gehrels_ft = 103.64602941062347

expected_multi_chi2_datavar_ff = 216.53516387719822
expected_multi_chi2_datavar_ft = 122.10913451612139

expected_multi_chi2_constvar_ff = 121.5229005671909
expected_multi_chi2_constvar_ft = 100.48181690303166

expected_multi_chi2_xspecvar_ff = expected_multi_chi2_datavar_ff
expected_multi_chi2_xspecvar_ft = expected_multi_chi2_datavar_ft

expected_multi_chi2_mod_stat = 139.04938684379346
expected_multi_chi2_mod_sys = 112.20863259556128

expected_multi_cash = -137151.91981923723
expected_multi_cstat = 142.0254918602011


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (LeastSq, False, False, expected_multi_leastsq),
    (LeastSq, True, False, expected_multi_leastsq),
    (LeastSq, False, True, expected_multi_leastsq),
    (LeastSq, True, True, expected_multi_leastsq),

    (Chi2, True, False, expected_multi_chi2),
    (Chi2, True, True, expected_multi_chi2_tt),

    (Chi2Gehrels, False, False, expected_multi_chi2_gehrels_ff),
    (Chi2Gehrels, True, False, expected_multi_chi2),
    (Chi2Gehrels, False, True, expected_multi_chi2_gehrels_ft),
    (Chi2Gehrels, True, True, expected_multi_chi2_tt),

    (Chi2DataVar, False, False, expected_multi_chi2_datavar_ff),
    (Chi2DataVar, True, False, expected_multi_chi2),
    (Chi2DataVar, False, True, expected_multi_chi2_datavar_ft),
    (Chi2DataVar, True, True, expected_multi_chi2_tt),

    (Chi2ConstVar, False, False, expected_multi_chi2_constvar_ff),
    (Chi2ConstVar, True, False, expected_multi_chi2),
    (Chi2ConstVar, False, True, expected_multi_chi2_constvar_ft),
    (Chi2ConstVar, True, True, expected_multi_chi2_tt),

    (Chi2ModVar, False, False, expected_multi_chi2_mod_stat),
    (Chi2ModVar, True, False, expected_multi_chi2_mod_stat),
    (Chi2ModVar, False, True, expected_multi_chi2_mod_sys),
    (Chi2ModVar, True, True, expected_multi_chi2_mod_sys),

    (Chi2XspecVar, False, False, expected_multi_chi2_xspecvar_ff),
    (Chi2XspecVar, True, False, expected_multi_chi2),
    (Chi2XspecVar, False, True, expected_multi_chi2_xspecvar_ft),
    (Chi2XspecVar, True, True, expected_multi_chi2_tt),

    (Cash, False, False, expected_multi_cash),
    (Cash, True, False, expected_multi_cash),
    (Cash, False, True, expected_multi_cash),
    (Cash, True, True, expected_multi_cash),

    (CStat, False, False, expected_multi_cstat),
    (CStat, True, False, expected_multi_cstat),
    (CStat, False, True, expected_multi_cstat),
    (CStat, True, True, expected_multi_cstat)
])
def test_fit_calc_stat_multiple(stat, usestat, usesys, expected):
    """Test the results from the calc_stat method of Fit for 3 data sets.

    The idea is to check that the data sets are being sent around
    correctly, particularly when they have different sizes.

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    # For now restrict to the default statistic values
    #
    statobj = stat()
    fit, fits = setup_stat_multiple(statobj, usestat, usesys)

    assert_almost_equal(fit.calc_stat(), expected)

    # Test that the sum of the individual cases matches the expected
    # value.
    stats = [f.calc_stat() for f in fits]
    assert_almost_equal(sum(stats), expected)


# The choice of model parameter values does not create "good" fits,
# hence the q values are all << 1e-7, which means that the current
# checks can't distinguish them from zero.
#
q_multi_chi2_tf = 0.0  # ~ 1.8e-43
q_multi_chi2_tt = 1.69072e-06

q_multi_mod_statonly = 0.0
q_multi_mod_sysonly = 0.0

q_multi_xspecvar_ff = 0.0  # ~ 5e-22
q_multi_xspecvar_ft = 0.0

q_multi_cstat = 0.0  # ~ 4e-26


@pytest.mark.parametrize("stat,usestat,usesys,expected,qval", [
    (LeastSq, False, False, expected_multi_leastsq, None),
    (LeastSq, True, False, expected_multi_leastsq, None),
    (LeastSq, False, True, expected_multi_leastsq, None),
    (LeastSq, True, True, expected_multi_leastsq, None),

    (Chi2, True, False, expected_multi_chi2, q_multi_chi2_tf),
    (Chi2, True, True, expected_multi_chi2_tt, q_multi_chi2_tt),

    (Chi2Gehrels, False, False, expected_multi_chi2_gehrels_ff, 0.0),
    (Chi2Gehrels, True, False, expected_multi_chi2, q_multi_chi2_tf),
    (Chi2Gehrels, False, True, expected_multi_chi2_gehrels_ft, 0.0),
    (Chi2Gehrels, True, True, expected_multi_chi2_tt, q_multi_chi2_tt),

    (Chi2DataVar, False, False, expected_multi_chi2_datavar_ff, 0.0),
    (Chi2DataVar, True, False, expected_multi_chi2, q_multi_chi2_tf),
    (Chi2DataVar, False, True, expected_multi_chi2_datavar_ft, 0.0),
    (Chi2DataVar, True, True, expected_multi_chi2_tt, q_multi_chi2_tt),

    (Chi2ConstVar, False, False, expected_multi_chi2_constvar_ff, 0.0),
    (Chi2ConstVar, True, False, expected_multi_chi2, q_multi_chi2_tf),
    (Chi2ConstVar, False, True, expected_multi_chi2_constvar_ft, 0.0),
    (Chi2ConstVar, True, True, expected_multi_chi2_tt,
     q_multi_chi2_tt),

    (Chi2ModVar, False, False, expected_multi_chi2_mod_stat, 0.0),
    (Chi2ModVar, True, False, expected_multi_chi2_mod_stat, 0.0),
    (Chi2ModVar, False, True, expected_multi_chi2_mod_sys, 0.0),
    (Chi2ModVar, True, True, expected_multi_chi2_mod_sys, 0.0),

    (Chi2XspecVar, False, False, expected_multi_chi2_xspecvar_ff,
     q_multi_xspecvar_ff),
    (Chi2XspecVar, True, False, expected_multi_chi2, q_multi_chi2_tf),
    (Chi2XspecVar, False, True, expected_multi_chi2_xspecvar_ft,
     q_multi_xspecvar_ft),
    (Chi2XspecVar, True, True, expected_multi_chi2_tt,
     q_multi_chi2_tt),

    (Cash, False, False, expected_multi_cash, None),
    (Cash, True, False, expected_multi_cash, None),
    (Cash, False, True, expected_multi_cash, None),
    (Cash, True, True, expected_multi_cash, None),

    (CStat, False, False, expected_multi_cstat, q_multi_cstat),
    (CStat, True, False, expected_multi_cstat, q_multi_cstat),
    (CStat, False, True, expected_multi_cstat, q_multi_cstat),
    (CStat, True, True, expected_multi_cstat, q_multi_cstat)
])
def test_fit_calc_stat_info_multiple(stat, usestat, usesys, expected, qval):
    """Test the results from the calc_stat_info method of Fit for 3 data sets.

    The idea is to check that the data sets are being sent around
    correctly, particularly when they have different sizes.

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    # channels, per dataset
    nbins = [3, 4, 5]
    nchannels = sum(nbins)

    # number of free parameters, per model
    nfree = [1, 1, 1]

    ndof = nchannels - sum(nfree)

    # This test does not use the individual data set results
    statobj = stat()
    fit, _ = setup_stat_multiple(statobj, usestat, usesys)

    ans = fit.calc_stat_info()

    assert isinstance(ans, StatInfoResults)
    assert ans.name == ""
    assert ans.ids is None
    assert ans.bkg_ids is None
    assert ans.numpoints == nchannels
    assert ans.dof == ndof
    assert_almost_equal(ans.statval, expected)

    # The choice of ans.statname is interesting for chi-squared
    # cases, since it depends on whether usestat is True.
    #
    if statobj.name.startswith('chi2') and usestat:
        assert ans.statname == 'chi2'
    else:
        assert ans.statname == statobj.name

    if qval is None:
        assert ans.qval is None
        assert ans.rstat is None

    else:
        assert_almost_equal(ans.qval, qval)
        assert_almost_equal(ans.rstat, ans.statval / ndof)


chisqr_multi_leastsq = [4.0, 1.0, 9.0, 2500.0, 4900.0, 6400.0,
                        4.0, 4.0, 1.0, 1.0, 9.0, 4.0]

chisqr_multi_chi2 = [2.777777778, 0.591715976, 14.062500000,
                     2.777777778, 3.062500000, 10.240000000,
                     0.027777778, 100.000000000, 11.111111111,
                     8.163265306, 56.250000000, 16.000000000]
chisqr_multi_chi2_tt = [1.277955272, 0.319488818, 13.235294118,
                        2.702702703, 2.995110024, 9.984399376,
                        0.026143791, 3.846153846, 0.244498778,
                        0.242571255, 7.758620690, 0.941176471]

chisqr_multi_gehrels_ff = [0.191466084, 0.040506733, 0.331172741,
                           14.179887089, 114.377652560, 30.675043738,
                           0.033437936, 0.001910900, 0.000476097,
                           0.000473333, 0.004235405, 0.001877886]
chisqr_multi_gehrels_ft = [0.177136713, 0.038274210, 0.330686011,
                           12.418900956, 62.150769352, 28.490198523,
                           0.031098249, 0.001909988, 0.000475192,
                           0.000472438, 0.004233413, 0.001874366]

chisqr_multi_datavar_ff = [0.333333333, 0.066666667, 0.529411765,
                           16.666666667, 163.333333333, 35.555555556,
                           0.040816327, 0.001998002, 0.000497760,
                           0.000494805, 0.004426955, 0.001962709]
chisqr_multi_datavar_ft = [0.292184076, 0.060827251, 0.528169014,
                           14.285714286, 74.242424242, 32.653061224,
                           0.037383178, 0.001997004, 0.000496771,
                           0.000493827, 0.004424779, 0.001958864]

chisqr_multi_constvar_ff = [0.272727273, 0.068181818, 0.613636364,
                            21.834061135, 42.794759825, 55.895196507,
                            0.034934498, 0.001979610, 0.000494903,
                            0.000494903, 0.004454123, 0.001979610]
chisqr_multi_constvar_ft = [0.244548604, 0.062086093, 0.611967362,
                            17.921146953, 32.558139535, 49.042145594,
                            0.032388664, 0.001978631, 0.000493925,
                            0.000493925, 0.004451919, 0.001975699]

chisqr_multi_mod_stat = [0.285714286, 0.071428571, 0.642857143,
                         25.000000000, 49.000000000, 64.000000000,
                         0.040000000, 0.002000000, 0.000497512,
                         0.000495050, 0.004433498, 0.001960784]
chisqr_multi_mod_sys = [0.254939452, 0.064766839, 0.641025641,
                        20.000000000, 36.029411765, 55.172413793,
                        0.036697248, 0.001999000, 0.000496524,
                        0.000494071, 0.004431315, 0.001956947]

chisqr_multi_xspecvar_ff = [0.333333333, 0.066666667, 0.529411765,
                            16.666666667, 163.333333333, 35.555555556,
                            0.040816327, 0.001998002, 0.000497760,
                            0.000494805, 0.004426955, 0.001962709]
chisqr_multi_xspecvar_ft = [0.292184076, 0.060827251, 0.528169014,
                            14.285714286, 74.242424242, 32.653061224,
                            0.037383178, 0.001997004, 0.000496771,
                            0.000493827, 0.004424779, 0.001958864]

chisqr_multi_cash = None
chisqr_multi_cstat = None


@pytest.mark.parametrize("stat,usestat,usesys,expected", [
    (LeastSq, False, False, chisqr_multi_leastsq),
    (LeastSq, True, False, chisqr_multi_leastsq),
    (LeastSq, False, True, chisqr_multi_leastsq),
    (LeastSq, True, True, chisqr_multi_leastsq),

    (Chi2, True, False, chisqr_multi_chi2),
    (Chi2, True, True, chisqr_multi_chi2_tt),

    (Chi2Gehrels, False, False, chisqr_multi_gehrels_ff),
    (Chi2Gehrels, True, False, chisqr_multi_chi2),
    (Chi2Gehrels, False, True, chisqr_multi_gehrels_ft),
    (Chi2Gehrels, True, True, chisqr_multi_chi2_tt),

    (Chi2DataVar, False, False, chisqr_multi_datavar_ff),
    (Chi2DataVar, True, False, chisqr_multi_chi2),
    (Chi2DataVar, False, True, chisqr_multi_datavar_ft),
    (Chi2DataVar, True, True, chisqr_multi_chi2_tt),

    (Chi2ConstVar, False, False, chisqr_multi_constvar_ff),
    (Chi2ConstVar, True, False, chisqr_multi_chi2),
    (Chi2ConstVar, False, True, chisqr_multi_constvar_ft),
    (Chi2ConstVar, True, True, chisqr_multi_chi2_tt),

    (Chi2ModVar, False, False, chisqr_multi_mod_stat),
    (Chi2ModVar, True, False, chisqr_multi_mod_stat),
    (Chi2ModVar, False, True, chisqr_multi_mod_sys),
    (Chi2ModVar, True, True, chisqr_multi_mod_sys),

    (Chi2XspecVar, False, False, chisqr_multi_xspecvar_ff),
    (Chi2XspecVar, True, False, chisqr_multi_chi2),
    (Chi2XspecVar, False, True, chisqr_multi_xspecvar_ft),
    (Chi2XspecVar, True, True, chisqr_multi_chi2_tt),

    (Cash, False, False, chisqr_multi_cash),
    (Cash, True, False, chisqr_multi_cash),
    (Cash, False, True, chisqr_multi_cash),
    (Cash, True, True, chisqr_multi_cash),

    (CStat, False, False, chisqr_multi_cstat),
    (CStat, True, False, chisqr_multi_cstat),
    (CStat, False, True, chisqr_multi_cstat),
    (CStat, True, True, chisqr_multi_cstat)
])
def test_fit_calc_chisqr_multiple(stat, usestat, usesys, expected):
    """Test the results from the calc_chisqr method of Fit for 3 data sets.

    The idea is to check that the data sets are being sent around
    correctly, particularly when they have different sizes.

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    # For now restrict the tests to using the default statistic values
    #
    fit, fits = setup_stat_multiple(stat(), usestat, usesys)
    ans = fit.calc_chisqr()

    allans = [f.calc_chisqr() for f in fits]

    if expected is None:
        assert ans is None
        assert all([a is None for a in allans])
    else:
        assert_almost_equal(ans, expected)
        combined = np.concatenate(allans)
        assert_almost_equal(combined, expected)


wstat_single_scalar_stat = 18.907726418126835
wstat_single_scalar_lo_stat = 17.094447413086726
wstat_single_scalar_hi_stat = 15.29434361363597
wstat_single_scalar_mid_stat = 13.481064608595863

wstat_single_array_stat = 18.895406024373457

# The expected values have been calculated by
# manually filtering the backscal array, and so should be
# correct (i.e. the tests should fail once #248 is fixed).
#
wstat_single_array_lo_stat = 17.08212701933335
wstat_single_array_hi_stat = 15.282023219882593
wstat_single_array_mid_stat = 13.468744214842486


# There are some other combinations that we currently do not
# exercise for wstat:
#   - grouped vs ungrouped (WSTAT is primarily for ungrouped data,
#     but it is possible to use with grouped data)
#   - multiple backgrounds
#
@pytest.mark.parametrize("scalar,usestat,usesys,flo,fhi,expected", [
    (True, False, False, None, None, wstat_single_scalar_stat),
    (True, True, False, None, None, wstat_single_scalar_stat),
    (True, False, True, None, None, wstat_single_scalar_stat),
    (True, True, True, None, None, wstat_single_scalar_stat),
    # I am not going to iterate through every usestat/sys combination
    # with flo/fhi options
    #
    # Start with some edge cases
    (True, False, False, -10, 10, wstat_single_scalar_stat),
    (True, False, False, 1, 5, wstat_single_scalar_stat),
    #
    # Now some more-general choices
    (True, False, False, 2, None, wstat_single_scalar_lo_stat),
    (True, False, False, 1.2, None, wstat_single_scalar_lo_stat),
    (True, True, True, 2, None, wstat_single_scalar_lo_stat),
    (True, False, False, None, 4.1, wstat_single_scalar_hi_stat),
    (True, True, True, None, 4.1, wstat_single_scalar_hi_stat),
    (True, False, False, 1.2, 4.1, wstat_single_scalar_mid_stat),
    (True, False, False, 2, 4.1, wstat_single_scalar_mid_stat),
    (True, True, True, 1.2, 4.1, wstat_single_scalar_mid_stat),

    # switch to an array for the background backscal
    (False, False, False, None, None, wstat_single_array_stat),
    (False, True, False, None, None, wstat_single_array_stat),
    (False, False, True, None, None, wstat_single_array_stat),
    (False, True, True, None, None, wstat_single_array_stat),

    (False, False, False, -1, 6, wstat_single_array_stat),
    (False, False, False, 1, 5, wstat_single_array_stat),

    (False, True, False, 1.2, None, wstat_single_array_lo_stat),
    (False, False, True, 2, None, wstat_single_array_lo_stat),

    (False, False, False, None, 4.1, wstat_single_array_hi_stat),
    (False, False, False, 1.2, 4.1, wstat_single_array_mid_stat),
    (False, False, False, 2, 4.1, wstat_single_array_mid_stat),
    (False, False, False, 2, 4, wstat_single_array_mid_stat)
])
def test_fit_calc_stat_wstat_single(scalar, usestat, usesys,
                                    flo, fhi, expected):
    """Test the results from the calc_stat method of Fit for 1 data set: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit = setup_pha_single(scalar, usestat, usesys, flo, fhi)
    assert_almost_equal(fit.calc_stat(), expected)


@pytest.mark.parametrize("scalar,usestat,usesys,flo,fhi,nbin,expected", [
    (True, False, False, None, None, 5, wstat_single_scalar_stat),
    (True, True, False, None, None, 5, wstat_single_scalar_stat),
    (True, False, True, None, None, 5, wstat_single_scalar_stat),
    (True, True, True, None, None, 5, wstat_single_scalar_stat),

    (True, False, False, -10, 10, 5, wstat_single_scalar_stat),
    (True, False, False, 1, 5, 5, wstat_single_scalar_stat),

    (True, False, False, 2, None, 4, wstat_single_scalar_lo_stat),
    (True, False, False, 1.2, None, 4, wstat_single_scalar_lo_stat),
    (True, True, True, 2, None, 4, wstat_single_scalar_lo_stat),
    (True, False, False, None, 4.1, 4, wstat_single_scalar_hi_stat),
    (True, True, True, None, 4.1, 4, wstat_single_scalar_hi_stat),
    (True, False, False, 1.2, 4.1, 3, wstat_single_scalar_mid_stat),
    (True, False, False, 2, 4.1, 3, wstat_single_scalar_mid_stat),
    (True, True, True, 1.2, 4.1, 3, wstat_single_scalar_mid_stat),

    (False, False, False, None, None, 5, wstat_single_array_stat),
    (False, True, False, None, None, 5, wstat_single_array_stat),
    (False, False, True, None, None, 5, wstat_single_array_stat),
    (False, True, True, None, None, 5, wstat_single_array_stat),

    (False, False, False, -1, 6, 5, wstat_single_array_stat),
    (False, False, False, 1, 5, 5, wstat_single_array_stat),

    (False, True, False, 1.2, None, 4, wstat_single_array_lo_stat),
    (False, False, True, 2, None, 4, wstat_single_array_lo_stat),

    (False, False, False, None, 4.1, 4, wstat_single_array_hi_stat),
    (False, False, False, 1.2, 4.1, 3, wstat_single_array_mid_stat),
    (False, False, False, 2, 4.1, 3, wstat_single_array_mid_stat),
    (False, False, False, 2, 4, 3, wstat_single_array_mid_stat)
])
def test_fit_calc_stat_info_wstat_single(scalar, usestat, usesys,
                                         flo, fhi, nbin, expected):
    """Test the results from the calc_stat_info method of Fit for 1 data set: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit = setup_pha_single(scalar, usestat, usesys, flo, fhi)
    ans = fit.calc_stat_info()

    # mdl.c0 is the only thawed parameter
    dof = nbin - 1

    assert_almost_equal(fit.calc_stat(), expected)
    assert isinstance(ans, StatInfoResults)
    assert ans.name == ""
    assert ans.statname == "wstat"
    assert ans.ids is None
    assert ans.bkg_ids is None
    assert ans.numpoints == nbin
    assert ans.dof == dof
    assert_almost_equal(ans.statval, expected)

    # TODO: as this is derived from cstat, shouldn't it calculate
    #       a q value and reduced stat?
    assert ans.qval is None
    assert ans.rstat is None

    # assert_almost_equal(ans.qval, qval)
    # assert_almost_equal(ans.rstat, ans.statval / dof)


# Do not really need to go through all these options, but it's easy
# to keep them.
@pytest.mark.parametrize("scalar,usestat,usesys,flo,fhi", [
    (True, False, False, None, None),
    (True, True, False, None, None),
    (True, False, True, None, None),
    (True, True, True, None, None),
    (True, False, False, -10, 10),
    (True, False, False, 1, 5),
    (True, False, False, 2, None),
    (True, False, False, 1.2, None),
    (True, True, True, 2, None),
    (True, False, False, None, 4.1),
    (True, True, True, None, 4.1),
    (True, False, False, 1.2, 4.1),
    (True, False, False, 2, 4.1),
    (True, True, True, 1.2, 4.1),
    (False, False, False, None, None),
    (False, True, False, None, None),
    (False, False, True, None, None),
    (False, True, True, None, None),
    (False, False, False, -1, 6),
    (False, False, False, 1, 5),
    (False, True, False, 1.2, None),
    (False, False, True, 2, None),
    (False, False, False, None, 4.1),
    (False, False, False, 1.2, 4.1),
    (False, False, False, 2, 4.1),
    (False, False, False, 2, 4)
])
def test_fit_calc_chisqr_wstat_single(scalar, usestat, usesys,
                                      flo, fhi):
    """Test the results from the calc_chisqr method of Fit for 1 data set: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit = setup_pha_single(scalar, usestat, usesys, flo, fhi)
    ans = fit.calc_chisqr()
    assert ans is None


# TODO: is the filtering applied to both src and background or are
#       they different?
#       ditto grouping and quality
#
# QUS: how is background subtraction handled, as that is essentially the same
# as wstat calculation in terms of source + bgnd handling for those
# settings which can be different.

@pytest.mark.parametrize("flo,fhi,expected", [
    # The following four fail in CIAO 4.8 so the expected values
    # were created with an early version of a fix for #248
    # (using a slightly different approach to the final fix,
    # see #269, commit 332feb2ff70595d28fc4d6cd9e2f158a95ee46d9)
    #
    (None, None, 179.593636089),
    (3, None, 177.863206820),
    (None, 12, 77.7333219274),
    (5, 15, 97.7590037051),

])
def test_fit_calc_stat_wstat_grouped_single(flo, fhi, expected):
    """Test the results from the calc_stat method of Fit for 1 data set: wstat

    This uses grouped data; it could be included in
    test_fit_calc_stat_wstat_single but it's getting a bit messy
    so pull out into a separate test.

    """

    nbins = 20
    channels = np.arange(1, nbins + 1, dtype=np.int)
    src_counts = (10 + 5 * np.sin(channels / 2.0)).astype(np.int8)
    bg_counts = np.ones(nbins, dtype=np.int8)

    # grouping: 1 to start a bin, -1 to continue a bin
    # quality: 0 for good, 1 or more for bad
    #
    # The code behaves differently if a NumPy array or Python
    # list is used for the grouping column, so force a numpy
    # array for now.
    #
    grouping = np.asarray([1, -1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1,
                           1, 1, 1, 1, 1, -1, -1, -1])
    assert len(grouping) == nbins, 'internal error; can not count'

    src = DataPHA('tst', channels, src_counts, grouping=grouping,
                  exposure=100.0, backscal=0.01)

    # The backscal array should be variable so that it has to be
    # "averaged" in some way by the grouping.
    #
    backscal = 0.04 - 0.00015 * channels
    bg = DataPHA('tstbg', channels, bg_counts,
                 exposure=200.0, backscal=backscal)

    src.set_background(bg)

    # Filtering is *only* applied to the source
    src.notice(lo=flo, hi=fhi)

    # Pick something that varies with the channel number, to
    # make sure that filtering is handled correctly.
    mdl = Polynom1D('mdl')
    mdl.c0 = 8
    mdl.c1 = 1.2

    # For now restrict to the default statistic values
    fit = Fit(src, mdl, stat=WStat())
    assert_almost_equal(fit.calc_stat(), expected)


wstat_multi_all = 28.114709948

# The following three fail in CIAO 4.8 so the expected values
# were created with an early version of a fix for #248
# (using a slightly different approach to the final fix,
# see #269, commit 332feb2ff70595d28fc4d6cd9e2f158a95ee46d9)
#
wstat_multi_lo = 27.321455143
wstat_multi_hi = 27.970661573
wstat_multi_mid = 27.177406768


@pytest.mark.parametrize("flo,fhi,expected", [
    (None, None, wstat_multi_all),
    (None, 1000, wstat_multi_all),
    (0, None, wstat_multi_all),
    (0, 1000, wstat_multi_all),
    (1, None, wstat_multi_all),
    (1, 1000, wstat_multi_all),

    (2, None, wstat_multi_lo),
    (None, 8, wstat_multi_hi),
    (2, 8, wstat_multi_mid),
])
def test_fit_calc_stat_wstat_multiple(flo, fhi, expected):
    """Test the results from the calc_stat method of Fit for 3 data sets: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit, fits = setup_pha_multiple(flo, fhi)
    assert_almost_equal(fit.calc_stat(), expected)

    # Test that the sum of the individual cases matches the expected
    # value.
    stats = [f.calc_stat() for f in fits]
    assert_almost_equal(sum(stats), expected)


# Remember, the filter is only applied to the background for the third
# dataset, so it is ignored (as only fitler/grouping applied to the
# source is important for WStat).
#
wstat_multi_nbins = 9 + 5 + 11
wstat_multi_lo_nbins = wstat_multi_nbins - 2
wstat_multi_hi_nbins = wstat_multi_nbins - 1
wstat_multi_mid_nbins = wstat_multi_nbins - 3


@pytest.mark.parametrize("flo,fhi,nbin,expected", [
    (None, None, wstat_multi_nbins, wstat_multi_all),
    (None, 1000, wstat_multi_nbins, wstat_multi_all),
    (0, None, wstat_multi_nbins, wstat_multi_all),
    (0, 1000, wstat_multi_nbins, wstat_multi_all),
    (1, None, wstat_multi_nbins, wstat_multi_all),
    (1, 1000, wstat_multi_nbins, wstat_multi_all),

    (2, None, wstat_multi_lo_nbins, wstat_multi_lo),
    (None, 8, wstat_multi_hi_nbins, wstat_multi_hi),
    (2, 8, wstat_multi_mid_nbins, wstat_multi_mid),
])
def test_fit_calc_stat_info_wstat_multiple(flo, fhi, nbin, expected):
    """Test the results from the calc_stat_info method of Fit for 3 data sets: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit, _ = setup_pha_multiple(flo, fhi)
    ans = fit.calc_stat_info()

    dof = nbin - 5

    assert_almost_equal(fit.calc_stat(), expected)
    assert isinstance(ans, StatInfoResults)
    assert ans.name == ""
    assert ans.statname == "wstat"
    assert ans.ids is None
    assert ans.bkg_ids is None
    assert ans.numpoints == nbin
    assert ans.dof == dof
    assert_almost_equal(ans.statval, expected)

    # TODO: as this is derived from cstat, shouldn't it calculate
    #       a q value and reduced stat?
    assert ans.qval is None
    assert ans.rstat is None

    # assert_almost_equal(ans.qval, qval)
    # assert_almost_equal(ans.rstat, ans.statval / dof)


@pytest.mark.parametrize("flo,fhi,expected", [
    (None, None, wstat_multi_all),
    (None, 1000, wstat_multi_all),
    (0, None, wstat_multi_all),
    (0, 1000, wstat_multi_all),
    (1, None, wstat_multi_all),
    (1, 1000, wstat_multi_all),
    (2, None, wstat_multi_lo),
    (None, 8, wstat_multi_hi),
    (2, 8, wstat_multi_mid),
])
def test_fit_calc_chisqr_wstat_multiple(flo, fhi, expected):
    """Test the results from the calc_chisqr method of Fit for 3 data sets: wstat

    This was created using the CIAO 4.8 version of Sherpa
    5ba2a7acdba6ae23bdb9462a7a57f6a60cbed685 on Linux 64 to
    calculate the expected values.

    """

    fit, fits = setup_pha_multiple(flo, fhi)

    assert fit.calc_chisqr() is None
    assert all([f.calc_chisqr() is None for f in fits])


def test_fit_calc_stat_error_on_chi2():
    """Check that an error message is raised when staterror is None."""

    d = Data1D("x", np.asarray([1, 2]), np.asarray([3, 4]),
               staterror=None, syserror=None)
    m = Const1D()
    f = Fit(d, m, stat=Chi2())
    with pytest.raises(StatErr):
        f.calc_stat()

    d.syserror = [2, 3]
    with pytest.raises(StatErr):
        f.calc_stat()


def test_fit_calc_stat_error_on_wstat():
    """Check that raise errors when using wstat incorrectly"""

    # wstat will error out if there is no background; that is,
    # at present it only works with DataPHA data sets with an
    # associated background.
    #
    dbasic = Data1D("x", np.asarray([1, 2]), np.asarray([3, 4]),
                    staterror=None, syserror=None)
    dpha = DataPHA("x", np.asarray([1, 2]), np.asarray([3, 4]),
                   exposure=1.0, backscal=1.0)

    m = Const1D()

    for d in [dbasic, dpha]:
        f = Fit(d, m, stat=Chi2())
        with pytest.raises(StatErr):
            f.calc_stat()


# The Fit constructor should probably check, and error out if there
# is a size mis-match. These should really be a unit test of the
# stat class.
#
@pytest.mark.parametrize("stat", [
    LeastSq, Chi2, Chi2DataVar, Cash, CStat
])
def test_fit_calc_stat_error_on_mismatch(stat):
    """Check that an error message is raised when ndata != nmodel"""

    # Do not need the fit objects returned by setup_stat_multiple,
    # but easier to extract the information (data and models) from
    # them.
    #
    statobj = stat()
    _, fits = setup_stat_multiple(statobj, False, False)

    data = [f.data for f in fits]
    models = [f.model for f in fits]

    dobj = DataSimulFit('too-many-data', data)
    mobj = SimulFitModel('not-enough-models', models[0:2])
    f = Fit(dobj, mobj, stat=statobj)
    with pytest.raises(StatErr):
        f.calc_stat()

    dobj = DataSimulFit('not-enough-data', data[0:2])
    mobj = SimulFitModel('too-many-models', models)
    f = Fit(dobj, mobj, stat=statobj)
    with pytest.raises(StatErr):
        f.calc_stat()


def test_fit_calc_stat_error_on_mismatch_wstat():
    """Check that an error message is raised when ndata != nmodel + wstat"""

    statobj = WStat()
    _, fits = setup_pha_multiple(None, None)

    data = [f.data for f in fits]
    models = [f.model for f in fits]

    dobj = DataSimulFit('too-many-data', data)
    mobj = SimulFitModel('not-enough-models', models[0:2])
    f = Fit(dobj, mobj, stat=statobj)
    with pytest.raises(StatErr):
        f.calc_stat()

    dobj = DataSimulFit('not-enough-data', data[0:2])
    mobj = SimulFitModel('too-many-models', models)
    f = Fit(dobj, mobj, stat=statobj)
    with pytest.raises(StatErr):
        f.calc_stat()


def test_fit_calc_stat_error_no_cache():
    """Check that data/errors/models are not cached.

    At present the fit does not cache the results, which
    lets things be changed and the statistic re-evaluated.
    Add a check to ensure there's no cache-ing of
    some relevant values.
    """

    d = Data1D("x", np.asarray([1, 2]), np.asarray([3, 4]),
               staterror=[1.2, 1.3], syserror=None)
    m = Const1D()
    m.c0 = 3
    f = Fit(d, m, stat=Chi2())

    def ans(delta, dy):
        return (delta / dy) ** 2

    assert_almost_equal(f.calc_stat(), ans(1.0, 1.3))

    d.y = np.asarray([3.1, 4.1])
    d.staterror = np.asarray([1.1, 1.5])
    m.c0 = 4.1

    assert_almost_equal(f.calc_stat(), ans(1.0, 1.1))


@pytest.mark.parametrize("stat", [
    LeastSq, Chi2, Chi2DataVar, Cash, CStat, WStat, UserStat
])
def test_fit_str_single(stat):
    """Test the str method of the fit object: single dataset.

    Not an ideal test, since it essentially just repeats the code from
    the __str__ method. Note that the test also doesn't check that
    the data is "valid" (e.g. the wstat statistic can not be
    evaluated for this fit object since there's no background
    component).

    """

    statobj = stat()
    fit = setup_stat_single(statobj, False, False)
    out = str(fit)

    expected = [("data", "test"),
                ("model", "((poly * step) + bg)"),
                ("stat", stat.__name__),
                ("method", "LevMar"),
                ("estmethod", "Covariance")]
    expected = "\n".join(["{:9s} = {}".format(*e) for e in expected])

    assert out == expected


@pytest.mark.parametrize("stat", [
    LeastSq, Chi2, Chi2DataVar, Cash, CStat, WStat, UserStat
])
def test_fit_str_multiple(stat):
    """Test the str method of the fit object: multiple datasets.

    Not an ideal test, since it essentially just repeats the code from
    the __str__ method. Note that the test also doesn't check that
    the data is "valid" (e.g. the wstat statistic can not be
    evaluated for this fit object since there's no background
    component).

    """

    statobj = stat()
    fit, _ = setup_stat_multiple(statobj, True, True)
    out = str(fit)

    expected = [("data", "multidata"),
                ("model", "multimodel"),
                ("stat", stat.__name__),
                ("method", "LevMar"),
                ("estmethod", "Covariance")]
    expected = "\n".join(["{:9s} = {}".format(*e) for e in expected])

    assert out == expected


# Can include WStat in this check since the check for nobins is
# done before validating the data type matches the statistic.
#
@pytest.mark.parametrize("stat",
                         [LeastSq, Chi2, Chi2Gehrels, Chi2DataVar,
                          Chi2ConstVar, Chi2ModVar, Chi2XspecVar,
                          Cash, CStat, WStat])
def test_fit_single_fails_nobins(stat):
    """Check that the fit method fails: single dataset, no valid bins
    """

    # The statistical/systematic errors should not be relevant
    # to this check, so pick something that allows Chi2 to be
    # used.
    #
    statobj = stat()
    fit = setup_stat_single(statobj, True, False)
    fit.data.ignore(None, None)
    with pytest.raises(DataErr) as excinfo:
        fit.fit()

    emsg = 'mask excludes all data'
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("stat",
                         [LeastSq, Chi2, Chi2Gehrels, Chi2DataVar,
                          Chi2ConstVar, Chi2ModVar, Chi2XspecVar,
                          Cash, CStat, WStat])
def test_fit_single_fails_pha_nobins(stat):
    """Check that the fit method fails: single dataset, PHA, no valid bins

    This probably does not add much to test_fit_single_fails_nobins.
    """

    # The statistical/systematic errors should not be relevant
    # to this check, so pick something that allows Chi2 to be
    # used.
    #
    statobj = stat()
    fit = setup_pha_single(True, True, False, None, None, stat=statobj)
    fit.data.ignore(None, None)
    with pytest.raises(DataErr) as excinfo:
        fit.fit()

    emsg = 'mask excludes all data'
    assert str(excinfo.value) == emsg


# Note that prior to Sherpa 4.8.1 the error here would have been a
# ValueError with the message 'calculation of errors has failed
# using current statistic' coming from a call to self.get_staterror().
# This was changed by PR #153.
#
@pytest.mark.parametrize("usesys", [True, False])
def test_fit_single_fails_dvar(usesys):
    """Check that the fit method fails: single dataset, datavar, 0's

    Does the Chi2DataVar stat error out with an error message
    about the inability to calculate the statistic when there
    are 0 values in the dependent axis (statistical error)
    """

    statobj = Chi2DataVar()
    fit = setup_stat_single(statobj, False, usesys)
    with pytest.raises(FitErr) as excinfo:
        fit.fit()

    emsg = 'zeros found in uncertainties, consider using ' + \
           'calculated uncertainties'
    assert str(excinfo.value) == emsg


# Calculated using LevMar
#
fit_lsq = 26.1355932203

fit_chi2_tt = 3.76425157747
fit_chi2_tf = 8.32136015756

fit_gehrels_ff = 1.22882590966
fit_gehrels_ft = 1.04190448781

fit_xvar_ff = 2.1192140563
fit_xvar_ft = 1.61852598

fit_cvar_ff = 2.41996233522
fit_cvar_ft = 1.78836616916

fit_mvar_xt = 1.42493962338
fit_mvar_xf = 1.84372887153

fit_cash = -173.299239257
fit_cstat = 1.96098192899


@pytest.mark.parametrize("stat,usestat,usesys,finalstat", [
    (LeastSq, False, False, fit_lsq),
    (LeastSq, True, True, fit_lsq),
    (Chi2, True, True, fit_chi2_tt),
    (Chi2, True, False, fit_chi2_tf),
    (Chi2Gehrels, True, True, fit_chi2_tt),
    (Chi2Gehrels, True, False, fit_chi2_tf),
    (Chi2Gehrels, False, False, fit_gehrels_ff),
    (Chi2Gehrels, False, True, fit_gehrels_ft),
    (Chi2DataVar, True, True, fit_chi2_tt),
    (Chi2DataVar, True, False, fit_chi2_tf),
    (Chi2XspecVar, True, True, fit_chi2_tt),
    (Chi2XspecVar, True, False, fit_chi2_tf),
    (Chi2XspecVar, False, False, fit_xvar_ff),
    (Chi2XspecVar, False, True, fit_xvar_ft),
    (Chi2ConstVar, True, True, fit_chi2_tt),
    (Chi2ConstVar, True, False, fit_chi2_tf),
    (Chi2ConstVar, False, False, fit_cvar_ff),
    (Chi2ConstVar, False, True, fit_cvar_ft),
    (Chi2ModVar, True, True, fit_mvar_xt),
    (Chi2ModVar, True, False, fit_mvar_xf),
    (Chi2ModVar, False, False, fit_mvar_xf),
    (Chi2ModVar, False, True, fit_mvar_xt),

    (Cash, True, True, fit_cash),
    (Cash, True, False, fit_cash),
    (Cash, False, False, fit_cash),
    (Cash, False, True, fit_cash),

    (CStat, True, True, fit_cstat),
    (CStat, True, False, fit_cstat),
    (CStat, False, False, fit_cstat),
    (CStat, False, True, fit_cstat),
])
def test_fit_single(stat, usestat, usesys, finalstat):
    """Check that the fit method works: single dataset, successful fit

    This is a minimal test, in that it just checks the
    final statistic of the fit (not the model parameters).
    The data and models are not designed to give good fit
    results.
    """

    statobj = stat()
    fit = setup_stat_single(statobj, usestat, usesys)
    assert fit.method.name == 'levmar'
    fr = fit.fit()
    assert fr.succeeded is True
    assert_almost_equal(fr.statval, finalstat)


@pytest.mark.parametrize("stat,usestat,usesys,finalstat", [
    (LeastSq, False, False, fit_lsq),
    (Chi2, True, True, fit_chi2_tt),
    (Chi2Gehrels, True, True, fit_chi2_tt),
    (Chi2Gehrels, False, False, fit_gehrels_ff),
    (Chi2XspecVar, False, True, fit_xvar_ft),
    (Chi2ConstVar, False, True, fit_cvar_ft),
    (Chi2ModVar, True, True, fit_mvar_xt),
    (Chi2ModVar, True, False, fit_mvar_xf),
    (Chi2ModVar, False, False, fit_mvar_xf),
    (Chi2ModVar, False, True, fit_mvar_xt),

    (Cash, True, True, fit_cash),
    (CStat, True, True, fit_cstat),
])
def test_fit_single_nm(stat, usestat, usesys, finalstat):
    """Check that the fit method works: single dataset, successful fit, simplex

    Test several cases from test_fit_single but using the
    NelderMead optimiser instead of the default LevMar.
    """

    statobj = stat()
    fit = setup_stat_single(statobj, usestat, usesys)
    fit.method = NelderMead()
    fr = fit.fit()
    assert fr.succeeded is True
    assert_almost_equal(fr.statval, finalstat)


# Calculated using LevMar
#
fit_multi_lsq = 12992.8666667

fit_multi_chi2_tt = 27.6895333147
fit_multi_chi2_tf = 156.932316029

fit_multi_gehrels_ff = 136.866794565
fit_multi_gehrels_ft = 100.624578034

# The data sets used in the multiple case are different to
# the single case, and do not contain the problem values that
# caused the single case to error out.
#
fit_multi_dvar_ff = 171.96397815
fit_multi_dvar_ft = 117.858874425

fit_multi_xvar_ff = fit_multi_dvar_ff
fit_multi_xvar_ft = fit_multi_dvar_ft

fit_multi_cvar_ff = 114.086122486
fit_multi_cvar_ft = 92.2492728722

fit_multi_mvar_xt = 87.8410045128
fit_multi_mvar_xf = 107.836151226

fit_multi_cash = -137160.045124
fit_multi_cstat = 133.900187108


@pytest.mark.parametrize("stat,usestat,usesys,finalstat", [
    (LeastSq, False, False, fit_multi_lsq),
    (LeastSq, True, True, fit_multi_lsq),
    (Chi2, True, True, fit_multi_chi2_tt),
    (Chi2, True, False, fit_multi_chi2_tf),
    (Chi2Gehrels, True, True, fit_multi_chi2_tt),
    (Chi2Gehrels, True, False, fit_multi_chi2_tf),
    (Chi2Gehrels, False, False, fit_multi_gehrels_ff),
    (Chi2Gehrels, False, True, fit_multi_gehrels_ft),
    (Chi2DataVar, True, True, fit_multi_chi2_tt),
    (Chi2DataVar, True, False, fit_multi_chi2_tf),
    (Chi2DataVar, False, False, fit_multi_dvar_ff),
    (Chi2DataVar, False, True, fit_multi_dvar_ft),
    (Chi2XspecVar, True, True, fit_multi_chi2_tt),
    (Chi2XspecVar, True, False, fit_multi_chi2_tf),
    (Chi2XspecVar, False, False, fit_multi_xvar_ff),
    (Chi2XspecVar, False, True, fit_multi_xvar_ft),
    (Chi2ConstVar, True, True, fit_multi_chi2_tt),
    (Chi2ConstVar, True, False, fit_multi_chi2_tf),
    (Chi2ConstVar, False, False, fit_multi_cvar_ff),
    (Chi2ConstVar, False, True, fit_multi_cvar_ft),
    (Chi2ModVar, True, True, fit_multi_mvar_xt),
    (Chi2ModVar, True, False, fit_multi_mvar_xf),
    (Chi2ModVar, False, False, fit_multi_mvar_xf),
    (Chi2ModVar, False, True, fit_multi_mvar_xt),

    (Cash, True, True, fit_multi_cash),
    (Cash, True, False, fit_multi_cash),
    (Cash, False, False, fit_multi_cash),
    (Cash, False, True, fit_multi_cash),

    (CStat, True, True, fit_multi_cstat),
    (CStat, True, False, fit_multi_cstat),
    (CStat, False, False, fit_multi_cstat),
    (CStat, False, True, fit_multi_cstat),
])
def test_fit_multiple(stat, usestat, usesys, finalstat):
    """Check that the fit method works: multiple datasets, successful fit

    This is a minimal test, in that it just checks the
    final statistic of the fit (not the model parameters).
    The data and models are not designed to give good fit
    results.
    """

    statobj = stat()
    fit, _ = setup_stat_multiple(statobj, usestat, usesys)
    fr = fit.fit()
    assert fr.succeeded is True
    assert_almost_equal(fr.statval, finalstat)


# TODO: add explicit wstat fit tests

@pytest.mark.parametrize("method,estmethod,usestat,usesys", [
    (LevMar, Covariance, True, True),
    (NelderMead, Covariance, True, False),
    (LevMar, Covariance, False, False),
    (NelderMead, Covariance, False, True),
    (LevMar, Confidence, True, True),
    (NelderMead, Confidence, True, False),
    (LevMar, Confidence, False, False),
    (NelderMead, Confidence, False, True),
])
def test_est_errors_single_fails_lsq(method, estmethod, usestat, usesys):
    """Check that the est_errors method fails: single dataset, LeastSq

    The method setting should not affect this, but include as a check
    (do not check all combinations).
    """

    statobj = LeastSq()
    fit = setup_stat_single(statobj, usestat, usesys)
    fit.method = method()
    fit.estmethod = estmethod()
    fit.fit()

    with pytest.raises(EstErr) as excinfo:
        fit.est_errors()

    emsg = 'cannot estimate confidence limits with LeastSq'
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("method, estmethod", [
    (LevMar, Covariance),
    (LevMar, Confidence),
    (NelderMead, Covariance),
    (NelderMead, Confidence),
])
def test_est_errors_single_fails_nobins(method, estmethod):
    """Check that the est_errors method fails: single dataset, no bins

    Check that, if a user has changed the bins after a fit, that the
    code errors out when there are no bins.
    """

    # Assume do not need to loop over the statistic or stat/sys errors
    statobj = Chi2()
    fit = setup_stat_single(statobj, True, False)
    fit.method = method()
    fit.estmethod = estmethod()
    fit.fit()

    fit.data.ignore(None, None)
    with pytest.raises(DataErr) as excinfo:
        fit.est_errors()

    emsg = 'mask excludes all data'
    assert str(excinfo.value) == emsg


# Using Covariance does not trigger this error.
@pytest.mark.parametrize("stat, method, estmethod", [
    (Chi2, LevMar, Confidence),
    (Chi2, NelderMead, Confidence),
    # Since usestat=True in the call to setup_stat_single then
    # all the Chi2XXX statistics (other than ModVar) are going
    # to have the same error values
    (Chi2Gehrels, LevMar, Confidence),
    (Chi2Gehrels, NelderMead, Confidence),
    (Chi2ModVar, LevMar, Confidence),
    (Chi2ModVar, NelderMead, Confidence),
    (CStat, LevMar, Confidence),
    (CStat, NelderMead, Confidence),
])
def test_est_errors_single_fails_badfit(stat, method, estmethod):
    """Check that the est_errors method fails: single dataset, bad fit

    Check that if this is a bad fit the error estimate call fails.
    """

    statobj = stat()
    fit = setup_stat_single(statobj, True, False)
    fit.method = method()
    fit.estmethod = estmethod()

    with pytest.raises(EstErr) as excinfo:
        fit.est_errors()

    emsg = 'reduced statistic larger than 3'
    assert str(excinfo.value) == emsg


# Restrict the choices run here to save time; ideally the parameter
# values should be set rather than rely on a call to fit.
#
@pytest.mark.parametrize("stat,usestat,usesys", [
    (Chi2, True, True),
    (Chi2Gehrels, True, True),
    (Chi2Gehrels, True, False),
    (Chi2Gehrels, False, False),
    (Chi2Gehrels, False, True),
    (Chi2ModVar, True, True),
    (Chi2ModVar, True, False),
    (Cash, False, True),
    (CStat, False, True),
])
def test_est_errors_single(stat, usestat, usesys):
    """Check that the est_errors method works: single dataset, successful fit

    This is a minimal test, in that it just checks that the
    error routine runs without raising an error.
    """

    statobj = stat()
    fit = setup_stat_single(statobj, usestat, usesys)
    fit.estmethod = Covariance()
    fit.fit()

    result = fit.est_errors()

    # check that a new best-fit location was not found during error
    # analysis
    assert result.nfits == 0
