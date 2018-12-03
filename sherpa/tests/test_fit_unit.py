#
#  Copyright (C) 2016, 2018  Smithsonian Astrophysical Observatory
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
    Chi2ConstVar, Chi2ModVar, Chi2XspecVar, Likelihood, \
    Cash, CStat, WStat, UserStat

from sherpa.optmethods import LevMar, NelderMead, MonCar
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


def setup_pha_single(scalar, usestat, usesys, flo, fhi,
                     stat=None, itermethod=None):
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
    itermethod : {None, True, False}
        None means no iterated fit, True means sigmarej, False means primini.

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

    if itermethod is None:
        iopts = None
    elif itermethod:
        iopts = {'name': 'sigmarej', 'maxiters': 5,
                 'hrej': 3, 'lrej': 3, 'grow': 0}
    else:
        iopts = {'name': 'primini', 'maxiters': 10, 'tol': 1.0e-3}

    return Fit(src, mdl, stat=statobj, itermethod_opts=iopts)


def setup_pha_multiple(flo, fhi, stat=None):
    """Set up multiple PHA datasets.

    A sherpa.data.DataPHA instance is used, with an associated
    background data set.

    Parameters
    ----------
    flo, fhi : None or number
        The channel filters (used in a notice call). These filters are
        applied to source only, both, and background only.
    stat : sherpa.stats.Stat instance, optional
        The statistic object to use. If not set uses a WStat object.

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

    # For now restrict to the default statistic values
    if stat is None:
        statobj = WStat()
    else:
        statobj = stat

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
q_chi2_tf = 1.23237e-09  # needed for check of format() output
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

    # The logic for the statistical name is a bit complicated.
    statname = statobj.name
    if usestat and isinstance(statobj, Chi2) and type(statobj) != Chi2 and \
       type(statobj) != LeastSq:
        statname = 'chi2'

    assert ans.statname == statname

    if qval is None:
        assert ans.qval is None
        assert ans.rstat is None

    else:
        assert_almost_equal(ans.qval, qval)
        # as there's only 1 dof, statval == rstat
        assert_almost_equal(ans.rstat, ans.statval)

    # It is easiest just to do a direct string comparison, but
    # differences in precision/formatting can make this annoying.
    #
    def tostr_short(k, v):
        return "{:9s} = {}".format(k, v)

    def tostr_long(k, v):
        return "{:21s} = {}".format(k, v)

    # The __str__ and format methods use different formats when
    # displaying the numeric values. Unfortunately the format
    # method uses %g which - for the numbers here - can be less
    # than 7 decimal places, so assert_almost_equal is not
    # usable.
    #
    def checkval_short(s, k, v):
        """Check that s == 'k = <v>'"""

        assert s.startswith(tostr_short(k, ''))
        sval = s.split(' = ')[1]
        assert_almost_equal(float(sval), v)

    def checkval_long(s, k, v):
        """Check that s == 'k = <v>'"""

        assert s.startswith(tostr_long(k, ''))
        sval = s.split(' = ')[1]
        assert sval == "%g" % (v, )

    # Validate __str__
    #
    sans = str(ans).split("\n")
    assert sans[0] == tostr_short('name', '')
    assert sans[1] == tostr_short('ids', 'None')
    assert sans[2] == tostr_short('bkg_ids', 'None')

    assert sans[3] == tostr_short('statname', statname)

    checkval_short(sans[4], 'statval', expected)
    assert sans[5] == tostr_short('numpoints', '5')
    assert sans[6] == tostr_short('dof', '1')
    if qval is None:
        assert sans[7] == tostr_short('qval', 'None')
        assert sans[8] == tostr_short('rstat', 'None')
    else:
        checkval_short(sans[7], 'qval', qval)
        checkval_short(sans[8], 'rstat', expected)

    assert len(sans) == 9

    # Validate format
    #
    sans = ans.format().split("\n")

    assert sans[0] == tostr_long('Statistic', statname)
    checkval_long(sans[1], 'Fit statistic value', expected)
    assert sans[2] == tostr_long('Data points', '5')
    assert sans[3] == tostr_long('Degrees of freedom', '1')

    nlines = 4
    if qval is not None:
        checkval_long(sans[4], 'Probability [Q-value]', qval)
        checkval_long(sans[5], 'Reduced statistic', expected)
        nlines += 2

    assert len(sans) == nlines


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


# The qval values were calculated using changeset
# b47c94ec7a28b7572b0665ac735db404d0dad13b
qval_scalar = 0.000819438167689
qval_scalar_lo = 0.000675824858126
qval_scalar_hi = 0.00158163032114
qval_scalar_mid = 0.00118201779597

qval_array = 0.000824015790193
qval_array_lo = 0.000679780437499
qval_array_hi = 0.00159083325587
qval_array_mid = 0.00118932173205


@pytest.mark.parametrize("scalar,usestat,usesys,flo,fhi,nbin,expected,qval", [
    (True, False, False, None, None, 5, wstat_single_scalar_stat, qval_scalar),
    (True, True, False, None, None, 5, wstat_single_scalar_stat, qval_scalar),
    (True, False, True, None, None, 5, wstat_single_scalar_stat, qval_scalar),
    (True, True, True, None, None, 5, wstat_single_scalar_stat, qval_scalar),

    (True, False, False, -10, 10, 5, wstat_single_scalar_stat, qval_scalar),
    (True, False, False, 1, 5, 5, wstat_single_scalar_stat, qval_scalar),

    (True, False, False, 2, None, 4, wstat_single_scalar_lo_stat,
     qval_scalar_lo),
    (True, False, False, 1.2, None, 4, wstat_single_scalar_lo_stat,
     qval_scalar_lo),
    (True, True, True, 2, None, 4, wstat_single_scalar_lo_stat,
     qval_scalar_lo),

    (True, False, False, None, 4.1, 4, wstat_single_scalar_hi_stat,
     qval_scalar_hi),
    (True, True, True, None, 4.1, 4, wstat_single_scalar_hi_stat,
     qval_scalar_hi),

    (True, False, False, 1.2, 4.1, 3, wstat_single_scalar_mid_stat,
     qval_scalar_mid),
    (True, False, False, 2, 4.1, 3, wstat_single_scalar_mid_stat,
     qval_scalar_mid),
    (True, True, True, 1.2, 4.1, 3, wstat_single_scalar_mid_stat,
     qval_scalar_mid),

    (False, False, False, None, None, 5, wstat_single_array_stat, qval_array),
    (False, True, False, None, None, 5, wstat_single_array_stat, qval_array),
    (False, False, True, None, None, 5, wstat_single_array_stat, qval_array),
    (False, True, True, None, None, 5, wstat_single_array_stat, qval_array),

    (False, False, False, -1, 6, 5, wstat_single_array_stat, qval_array),
    (False, False, False, 1, 5, 5, wstat_single_array_stat, qval_array),

    (False, True, False, 1.2, None, 4, wstat_single_array_lo_stat,
     qval_array_lo),
    (False, False, True, 2, None, 4, wstat_single_array_lo_stat,
     qval_array_lo),

    (False, False, False, None, 4.1, 4, wstat_single_array_hi_stat,
     qval_array_hi),

    (False, False, False, 1.2, 4.1, 3, wstat_single_array_mid_stat,
     qval_array_mid),
    (False, False, False, 2, 4.1, 3, wstat_single_array_mid_stat,
     qval_array_mid),
    (False, False, False, 2, 4, 3, wstat_single_array_mid_stat,
     qval_array_mid),
])
def test_fit_calc_stat_info_wstat_single(scalar, usestat, usesys,
                                         flo, fhi, nbin, expected, qval):
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

    assert_almost_equal(ans.qval, qval)
    assert_almost_equal(ans.rstat, expected / dof)


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

# The qval values were calculated using changeset
# b47c94ec7a28b7572b0665ac735db404d0dad13b
qval_multi = 0.106711635305
qval_multi_lo = 0.0731604652765
qval_multi_hi = 0.0839948772314
qval_multi_mid = 0.0555167652247


@pytest.mark.parametrize("flo,fhi,nbin,expected, qval", [
    (None, None, wstat_multi_nbins, wstat_multi_all, qval_multi),
    (None, 1000, wstat_multi_nbins, wstat_multi_all, qval_multi),
    (0, None, wstat_multi_nbins, wstat_multi_all, qval_multi),
    (0, 1000, wstat_multi_nbins, wstat_multi_all, qval_multi),
    (1, None, wstat_multi_nbins, wstat_multi_all, qval_multi),
    (1, 1000, wstat_multi_nbins, wstat_multi_all, qval_multi),

    (2, None, wstat_multi_lo_nbins, wstat_multi_lo, qval_multi_lo),
    (None, 8, wstat_multi_hi_nbins, wstat_multi_hi, qval_multi_hi),
    (2, 8, wstat_multi_mid_nbins, wstat_multi_mid, qval_multi_mid),
])
def test_fit_calc_stat_info_wstat_multiple(flo, fhi, nbin, expected, qval):
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

    assert_almost_equal(ans.qval, qval)
    assert_almost_equal(ans.rstat, expected / dof)


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
    assert fr.succeeded
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
    assert fr.succeeded
    assert_almost_equal(fr.statval, finalstat)


# Since the background is being ignored in this fit (except for
# wstat), the status of the BACKSCAL value (whether it's an array
# of values or not) should make no difference to the fit for all
# but WStat.
#
fit_pha_lsq = 166.8
fit_pha_filt_lsq = 78.6

fit_pha_chi2 = 20.4880257038
fit_pha_filt_chi2 = 6.76891612899
fit_pha_gehrels = 11.0170341689
fit_pha_filt_gehrels = 4.15421342256
fit_pha_mvar = 14.425605409
fit_pha_filt_mvar = 5.73844232605

fit_pha_cash = -128.007204397
fit_pha_filt_cash = -138.731060074

fit_pha_cstat = 17.7443412653
fit_pha_filt_cstat = 6.42881185701

fit_pha_wstat_t = 18.5533195269
fit_pha_filt_wstat_t = 6.91786326729
fit_pha_wstat_f = 18.552838794

# This value is not from CIAO 4.8 since the code there does not
# calculate the correct value.
fit_pha_filt_wstat_f = 6.89372285597


@pytest.mark.parametrize("stat,scalar,usestat,usesys,filtflag,finalstat", [
    (LeastSq, True, False, False, False, fit_pha_lsq),
    (LeastSq, True, False, False, True, fit_pha_filt_lsq),
    (LeastSq, False, True, True, False, fit_pha_lsq),

    (Chi2, True, True, True, False, fit_pha_chi2),
    (Chi2, True, True, True, True, fit_pha_filt_chi2),
    (Chi2, False, True, True, False, fit_pha_chi2),

    (Chi2Gehrels, True, False, False, False, fit_pha_gehrels),
    (Chi2Gehrels, False, False, False, False, fit_pha_gehrels),
    (Chi2Gehrels, False, False, False, True, fit_pha_filt_gehrels),

    (Chi2ModVar, True, False, False, False, fit_pha_mvar),
    (Chi2ModVar, True, False, False, True, fit_pha_filt_mvar),
    (Chi2ModVar, False, False, False, False, fit_pha_mvar),

    (Cash, True, True, False, False, fit_pha_cash),
    (Cash, False, True, False, False, fit_pha_cash),
    (Cash, False, False, False, True, fit_pha_filt_cash),

    (CStat, True, True, False, False, fit_pha_cstat),
    (CStat, False, True, False, False, fit_pha_cstat),
    (CStat, False, False, False, True, fit_pha_filt_cstat),

    (WStat, True, True, False, False, fit_pha_wstat_t),
    (WStat, True, True, False, True, fit_pha_filt_wstat_t),
    (WStat, False, True, False, False, fit_pha_wstat_f),
    (WStat, False, True, False, True, fit_pha_filt_wstat_f),
])
def test_fit_single_pha(stat, scalar, usestat, usesys, filtflag, finalstat):
    """Check that the fit method works: single dataset, PHA, successful fit

    This is a minimal test, in that it just checks the
    final statistic of the fit (not the model parameters).
    The data and models are not designed to give good fit
    results.

    This combines two tests to avoid typing: all data and a simple
    filtered case.
    """

    statobj = stat()
    fit = setup_pha_single(scalar, usestat, usesys, None, None,
                           stat=statobj)
    assert fit.method.name == 'levmar'
    assert fit.data.subtracted is False

    numpoints = 5
    if filtflag:
        fit.data.ignore(3.8, 4.5)
        numpoints -= 1

    fr = fit.fit()
    assert fr.succeeded
    assert fr.numpoints == numpoints
    assert fr.dof == (numpoints - 1)
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
    assert fr.succeeded
    assert_almost_equal(fr.statval, finalstat)

    # As the datasets and models are independent, the fit
    # to just the third dataset should return the same parameter
    # values. Actually, this is incorrect, as the convergence is
    # not guaranteed to be the same in the single-fit versus
    # multiple-fit cases, so there could be small differences.
    # This is particularly the cases for the likelihood stats
    # since the LevMar optimiser in use here is not ideal in
    # these cases. This is why the number of decimal places is
    # much smaller for the likelihood cases.
    #
    # Recall the setup code to ensure separate models are used
    _, fits = setup_stat_multiple(statobj, usestat, usesys)
    fr2 = fits[2].fit()

    # assume only a single free parameter
    assert fr.parnames[-1] == fr2.parnames[-1]

    if isinstance(statobj, Likelihood):
        ndp = 0
    else:
        ndp = 7

    assert_almost_equal(fr2.parvals[-1], fr.parvals[-1],
                        decimal=ndp)


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


@pytest.mark.parametrize("stat", [Chi2, Chi2Gehrels, Cash, CStat])
def test_est_errors_multiple(stat):
    """Check that the est_errors method works: multiple datasets, successful fit

    This is a minimal test, in that it just checks that the
    error routine runs without raising an error.
    """

    statobj = stat()
    fit, _ = setup_stat_multiple(statobj, True, True)
    fit.estmethod = Covariance()
    fit.fit()

    result = fit.est_errors()

    # check that a new best-fit location was not found during error
    # analysis
    assert result.nfits == 0


@pytest.mark.parametrize("stat,scalar,usestat,usesys,filtflag", [
    (Chi2, False, True, True, False),
    (Chi2, False, True, True, True),
    (Chi2, True, True, True, True),
    (Chi2Gehrels, False, False, False, False),
    (Chi2Gehrels, False, False, False, True),
    (Chi2ModVar, True, True, True, True),
    (Chi2ModVar, True, True, True, False),
    (Cash, False, False, False, False),
    (Cash, False, False, False, True),
    (Cash, True, False, False, False),
    (Cash, True, False, False, True),
    (CStat, False, False, False, False),
    (CStat, False, False, False, True),
    (CStat, True, False, False, False),
    (CStat, True, False, False, True),
    (WStat, False, False, False, False),
    (WStat, False, False, False, True),
    (WStat, True, False, False, False),
    (WStat, True, False, False, True),
])
def test_est_errors_single_pha(stat, scalar, usestat, usesys, filtflag):
    """Check that the est_errors method works: single dataset, PHA, successful fit

    This is a minimal test, in that it just checks that the
    error routine runs without raising an error.
    """

    statobj = stat()
    fit = setup_pha_single(scalar, usestat, usesys, None, None,
                           stat=statobj)
    assert fit.method.name == 'levmar'
    assert fit.data.subtracted is False

    if filtflag:
        fit.data.ignore(3.8, 4.5)

    fit.estmethod = Covariance()
    fit.fit()
    result = fit.est_errors()

    # check that a new best-fit location was not found during error
    # analysis
    assert result.nfits == 0


# Due to the way the errors are set in setup_pha_multiple(), can not
# use Chi2 here.
#
@pytest.mark.parametrize("stat,flo,fhi", [
    (Chi2Gehrels, None, None),
    (Chi2Gehrels, 2, None),
    (Chi2Gehrels, None, 8),
    (Chi2Gehrels, 2, 8),
    (Chi2ModVar, None, None),
    (Chi2ModVar, 2, 8),
    (Cash, None, None),
    (Cash, 2, None),
    (Cash, None, 8),
    (Cash, 2, 8),
    (CStat, None, None),
    (CStat, 2, None),
    (CStat, None, 8),
    (CStat, 2, 8),
    (WStat, None, None),
    (WStat, 2, None),
    (WStat, None, 8),
    (WStat, 2, 8)
])
def test_est_errors_multiple_pha(stat, flo, fhi):
    """Check that the est_errors method works: multiple datasets, PHA, successful fit

    This is a minimal test, in that it just checks that the
    error routine runs without raising an error.
    """

    statobj = stat()
    fit, _ = setup_pha_multiple(flo, fhi, stat=statobj)
    fit.estmethod = Covariance()
    fit.fit()

    result = fit.est_errors()

    # check that a new best-fit location was not found during error
    # analysis
    assert result.nfits == 0


# Test iterated-fit methods.
#
def setup_single_iter(stat, sigmarej=True):
    """Create a data set and model for testing iterated-fit methods.

    A sherpa.data.Data1D instance is used.

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
        The statistic object to use
    sigmarej : bool, optional
        If True use sigmarej, otherwise Primini

    Returns
    -------
    fit : sherpa.fit.Fit instance
    """

    # data created from y = 10 + 2 * x with gaussian noise (scale=1.2)
    # added and then one point manually changed
    #
    x = np.linspace(0, 30, 7)
    y = np.asarray([8.722, 21.012, 30.780, 40.920, 22.84, 59.34, 68.20])
    dy = np.ones(x.size) * 1.2

    data = Data1D('idata', x, y, staterror=dy)

    mdl = Polynom1D("poly")
    mdl.c1.frozen = False

    # Default values taken from sherpa.ui.utils.clean()
    #
    if sigmarej:
        iopts = {'name': 'sigmarej', 'maxiters': 5,
                 'hrej': 3, 'lrej': 3, 'grow': 0}
    else:
        iopts = {'name': 'primini', 'maxiters': 10, 'tol': 1.0e-3}

    return Fit(data, mdl, stat=stat, itermethod_opts=iopts)


@pytest.mark.parametrize("stat", [LeastSq, Cash, CStat])
@pytest.mark.parametrize("sigmarej", [True, False])
def test_fit_iterfit_fails_nonchi2(stat, sigmarej):
    """Check that iterfit fails with non-chi2/leastsq stats"""

    statobj = stat()
    fit = setup_single_iter(statobj, sigmarej=sigmarej)
    with pytest.raises(FitErr) as excinfo:
        fit.fit()

    if sigmarej:
        emsg = "Sigma-rejection"
    else:
        emsg = "Primini's"

    emsg += " method requires a deviates array; use a chi-square  statistic"
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("stat", [LeastSq, Cash, CStat, WStat])
@pytest.mark.parametrize("sigmarej", [True, False])
def test_fit_iterfit_fails_nonchi2_wstat(stat, sigmarej):
    """Check that iterfit fails with wstat

    This needs a "PHA" data set to get to the error with WSTAT.
    Throw in the other stats just as a check.
    """

    statobj = stat()
    fit = setup_pha_single(True, False, False, None, None,
                           stat=statobj, itermethod=sigmarej)
    with pytest.raises(FitErr) as excinfo:
        fit.fit()

    if sigmarej:
        emsg = "Sigma-rejection"
    else:
        emsg = "Primini's"

    emsg += " method requires a deviates array; use a chi-square  statistic"
    assert str(excinfo.value) == emsg


@pytest.mark.parametrize("stat", [Chi2, Chi2Gehrels])
def test_fit_iterfit_single_sigmarej_chi2(stat):
    """Very limited test of iterated-fit code.

    Since setup_single_iter creates a staterror column then
    the Chi2-based statistics (module Chi2ModVar) should all
    give the same result.
    """

    statobj = stat()
    fit = setup_single_iter(statobj, sigmarej=True)

    # be explicit here since the result is not guaranteed to be a bool
    assert fit.data.mask is True

    fr = fit.fit()
    assert fr.succeeded

    # The results are much "worse" than the Chi2Gehrels case
    expected_mask = [True, True, False, False, False, False, False]
    assert np.all(fit.data.mask == expected_mask)

    assert_almost_equal(fr.statval, 0.0)

    mdl = fit.model
    assert_almost_equal(mdl.c0.val, 8.722)
    assert_almost_equal(mdl.c1.val, 2.458)

    assert fr.numpoints == 2
    assert fr.dof == 0


def test_fit_iterfit_single_sigmarej_chi2gehrels():
    """Very limited test of iterated-fit code."""

    # If remove the staterror column and use the data values
    # the fit is "better".
    #
    statobj = Chi2Gehrels()
    fit = setup_single_iter(statobj, sigmarej=True)

    fit.data.staterror = None

    # be explicit here since the result is not guaranteed to be a bool
    assert fit.data.mask is True

    fr = fit.fit()
    assert fr.succeeded

    # the discrepant point should be excluded
    expected_mask = [True, True, True, True, False, True, True]
    assert np.all(fit.data.mask == expected_mask)

    assert_almost_equal(fr.statval, 0.1914790757)

    mdl = fit.model
    assert_almost_equal(mdl.c0.val, 9.69168196604)
    assert_almost_equal(mdl.c1.val, 2.00600360315)

    assert fr.numpoints == 6
    assert fr.dof == 4


def test_fit_iterfit_single_sigmarej_ignore_chi2gehrels():
    """Very limited test of iterated-fit code.

    This ignores some data before the fit since this checks
    logic that is not tested above.
    """

    statobj = Chi2Gehrels()
    fit = setup_single_iter(statobj, sigmarej=True)

    fit.data.ignore(4, 6)
    fit.data.staterror = None

    # be explicit here since the result is not guaranteed to be a bool
    start_mask = [True, False, True, True, True, True, True]
    assert np.all(fit.data.mask == start_mask)

    fr = fit.fit()
    assert fr.succeeded

    # the discrepant point should be excluded
    expected_mask = [True, False, True, True, False, True, True]
    assert np.all(fit.data.mask == expected_mask)

    assert_almost_equal(fr.statval, 0.1245627587)

    mdl = fit.model
    assert_almost_equal(mdl.c0.val, 9.25537857670)
    assert_almost_equal(mdl.c1.val, 2.01845980545)

    assert fr.numpoints == 5
    assert fr.dof == 3


def test_wstat_rstat_qval_fields_not_none():
    """It turns out there are other tests that check if rstat/qval
    are populated, but leave these in.

    Note the tests do not check the numerical values (these are
    handled elsewhere), just that the qval and rstat fields are
    not None.
    """

    fit = setup_pha_single(True, False, False, None, None)

    s1 = fit.calc_stat_info()
    assert s1.rstat is not None
    assert s1.qval is not None

    fr = fit.fit()
    assert fr.rstat is not None
    assert fr.qval is not None

    s2 = fit.calc_stat_info()
    assert s2.rstat is not None
    assert s2.qval is not None
    assert s2.rstat < s1.rstat
    assert s2.qval > s1.qval


def validate(expected, val):
    """Compare numeric values which may also be NaN/None/pytest.approx values

    Try and make it easier for the test writer.

    Parameters
    ----------
    expected : float, None, NaN, or pytest.approx instance
        The expected value. If a float then the standard pytest.approx
        settings are used, but you can send in an instance with its own
        values.
    value : float, None, or NaN
        The value to test.

    Notes
    -----
    The test for whether we have a pytest.approx instance relies on
    internals of pytest.approx, so could be unreliable long term.

    """

    if expected is None:
        assert val is None
    elif hasattr(expected, 'expected'):
        # Assyme a pytest.approx value
        assert val == expected
    elif np.isfinite(expected):
        assert val == pytest.approx(expected)
    else:
        assert np.isnan(val)


# Add some basic tests to see how fits where the number of degrees-of-freedom
# is 0 or negative are handled.
#
# A number of earlier routines could be updated to use this, but
# leave that for a later update.
#
def assert_stat_info(statinfo, npoints, dof, statval, qval, rstat):
    """Ensure that the stat info matches the expected values.

    Parameters
    ----------
    statinfo : sherpa.stats.StatInfo instance
    npoints : int
        The expected number of data points.
    dof : float
        The expected number of degrees of freedom.
    statval : float
        The expected statistic value.
    qval : float or None or np.nan
        The expected qval value
    rstat : float or None or np.nan
        The expected rstat value.

    """

    assert statinfo.numpoints == npoints
    assert statinfo.dof == dof
    assert statinfo.statval == pytest.approx(statval)

    validate(qval, statinfo.qval)
    validate(rstat, statinfo.rstat)


def assert_fit_results(fitres, rstat):
    """Ensure that the fit results match the expected values.

    Limited checks.

    Parameters
    ----------
    fitres : sherpa.fit.FitResults instance
        This must have fitres.succeeded be True.
    rstat : None or np.nan or float or pytest.approx instance
        The expected rstat value.
    """

    validate(rstat, fitres.rstat)
    assert fitres.succeeded


# A "canary" to note when issue #563 (WStat assumes arrays are ndarrays
# but we do not ensure this, so how should this be fixed) is addressed
#
# Note that once #563 is fixed the call is expected to raise a
# sherpa.utils.err.StatErr error with the message
#   "No background data has been supplied. Use cstat"
#
def test_563_still_exists():
    d = Data1D('test', [1, 2, 3], [4, 5, 6])
    mdl = Polynom1D()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=WStat())
    with pytest.raises(AttributeError) as excinfo:
        fit.calc_stat_info()

    assert "'list' object has no attribute 'size'" in str(excinfo.value)


# Note that the dof=1 tests don't really add any coverage to other
# tests, since they are no different than dof>>1, but are included so that
# we have some "documentation" of the behavior for dof > 0, dof=0, dof < 0.
#

# These values were found by running the code (so are regression tests)
# rather than being calculated from first principles. It does not seem
# worth running all statistics, so I have chosen those with different
# behaviors:
#  cash  - likelihood-based, which doesn't have rstat/dof
#  cstat - likelihood-based but with a measure of goodness of fit
#  chi2  - our good old friendly chi-square statistic
#
dof1_cash = -18.2772161919084

dof1_cstat = 0.40863145212838614
dof1_cstat_qval = 0.522664944712

dof1_chi2 = 2.03
dof1_chi2_qval = 0.154220607327


@pytest.mark.parametrize("method", [LevMar, NelderMead, MonCar])
@pytest.mark.parametrize("stat,statargs", [(Cash, (dof1_cash, None, None)),
                                           (CStat, (dof1_cstat, dof1_cstat_qval, dof1_cstat)),
                                           (Chi2, (dof1_chi2, dof1_chi2_qval, dof1_chi2))])
def test_dof_1(method, stat, statargs):
    """DOF is 1"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    statinfo = fit.calc_stat_info()
    assert_stat_info(statinfo, 3, 1, statargs[0], statargs[1], statargs[2])


@pytest.mark.parametrize("method", [LevMar, NelderMead, MonCar])
@pytest.mark.parametrize("stat,statargs", [(Cash, (dof1_cash, None, None)),
                                           (CStat, (dof1_cstat, np.nan, np.nan)),
                                           (Chi2, (dof1_chi2, np.nan, np.nan))])
def test_dof_0(method, stat, statargs):
    """DOF is 0"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c2.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    statinfo = fit.calc_stat_info()
    assert_stat_info(statinfo, 3, 0, statargs[0], statargs[1], statargs[2])


@pytest.mark.parametrize("method", [LevMar, NelderMead, MonCar])
@pytest.mark.parametrize("stat,statargs",
                         [(Cash, (dof1_cash, None, None)),
                          (CStat, (dof1_cstat, np.nan, np.nan)),
                          (Chi2, (dof1_chi2, np.nan, np.nan))])
def test_dof_neg1(method, stat, statargs):
    """DOF is -1"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c2.thaw()
    mdl.c3.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    statinfo = fit.calc_stat_info()
    assert_stat_info(statinfo, 3, -1, statargs[0], statargs[1], statargs[2])


# Results from CIAO 4.10 on Linux Ubuntu. Expect to end up with
# c0 = 3.0, c1 = 1.0 to about 3dp and the rstat values, when not
# None, are ~ e-14 to e-30 - i.e. 0 to within 13 dp.
#
@pytest.mark.parametrize("method,stat,rstat",
                         [(LevMar, Cash, None),
                          (LevMar, CStat, pytest.approx(0.0, abs=1e-7)),
                          (LevMar, Chi2, 0.0),
                          (NelderMead, Cash, None),
                          (NelderMead, CStat, 0.0),
                          (NelderMead, Chi2, 0.0),
                          (MonCar, Cash, None),
                          (MonCar, CStat, 0.0),
                          (MonCar, Chi2, 0.0)])
def test_fit_dof_1(method, stat, rstat):
    """DOF is 1"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    fres = fit.fit()

    assert mdl.c0.val == pytest.approx(3.0, abs=0.001)
    assert mdl.c1.val == pytest.approx(1.0, abs=0.001)

    assert_fit_results(fres, rstat)


# Results from CIAO 4.10 on Linux Ubuntu. Expect to end up with
# c0 = 3.0, c1 = 1.0, both to about 3dp, and c2 = 0 to within
# about 1e-6, the rstat values are NaN
#
@pytest.mark.parametrize("method,stat,rstat",
                         [(LevMar, Cash, None),
                          (LevMar, CStat, np.nan),
                          (LevMar, Chi2, np.nan),
                          (NelderMead, Cash, None),
                          (NelderMead, CStat, np.nan),
                          (NelderMead, Chi2, np.nan),
                          (MonCar, Cash, None),
                          (MonCar, CStat, np.nan),
                          (MonCar, Chi2, np.nan)])
def test_fit_dof_0(method, stat, rstat):
    """DOF is 0"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c2.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    fres = fit.fit()

    assert mdl.c0.val == pytest.approx(3.0, abs=0.001)
    assert mdl.c1.val == pytest.approx(1.0, abs=0.001)
    assert mdl.c2.val == pytest.approx(0.0, abs=1e-6)

    assert_fit_results(fres, rstat)


# Results from CIAO 4.10 on Linux Ubuntu. The fit results are
# not consistent between different options (e.g. method and statistic),
# which is not surprising given that the fit is not-well constrained.
#
# The only check here is whether the fit succeeds or not, since this
# depends on the method.
#
@pytest.mark.parametrize("stat", [Cash, CStat, Chi2])
@pytest.mark.parametrize("method,success",
                         [(LevMar, False),
                          (NelderMead, True),
                          (MonCar, True)])
def test_fit_dof_neg1(stat, method, success):
    """DOF is -1"""

    d = Data1D('test', [1, 2, 3], [4, 5, 6], [1, 1, 1])
    mdl = Polynom1D()
    mdl.c1.thaw()
    mdl.c2.thaw()
    mdl.c3.thaw()
    mdl.c0 = 5.1

    fit = Fit(d, mdl, stat=stat(), method=method())
    fres = fit.fit()
    assert fres.succeeded == success
