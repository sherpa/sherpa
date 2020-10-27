#
#  Copyright (C) 2019, 2020  Smithsonian Astrophysical Observatory
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

# Test the appearance of warnings when plots are created. At present
# this is restricted to warnings about whether the displayed error bars
# are used in a fit or not.
#
# The tests require a minimum of pytest 3.3 (and may require 3.4, as it
# is not clear if the changes in the logging handling between 3.3 and 3.4
# are an issue here; I think they are not but have not tested it).
#

import logging

import numpy as np

import pytest

from sherpa.data import Data1D
from sherpa.models.basic import Const1D
from sherpa.plot import DataPlot, ModelPlot, ResidPlot, RatioPlot, \
    FitPlot, JointPlot
from sherpa.stats import LeastSq, Cash, Chi2DataVar, Chi2ModVar


def example_data():
    """Create a 1D example dataset"""

    x = np.array([0, 10, 20])
    y = np.array([10, 20, 40])
    return Data1D('exmp', x, y)


def example_model():
    """Create a 1D example model"""

    mdl = Const1D('exmp')
    mdl.c0 = 25
    return mdl


def check_for_warning(caplog, nwarn, statname=None):
    """Check that the warning was created ntimes.

    Parameters
    ----------
    caplog
        The caplog fixture sent to the test.
    nwarn : int
        The number of times the warning is created (>= 0).
    statname : None or str
        The statistic name (must be given if nwarn > 0).
    """

    assert len(caplog.records) == nwarn
    if nwarn == 0:
        return

    # How useful is it to check the exact message?
    expected = 'The displayed errorbars have been supplied with the ' + \
               'data or calculated using chi2xspecvar; the errors ' + \
               'are not used in fits with {}'.format(statname)

    for (log_name, log_level, message) in caplog.record_tuples:
        assert log_level == logging.WARNING
        assert log_name == 'sherpa.plot'
        assert message == expected


# NOTE: the following tests contain a lot of repeated code but I don't
#       want to refactor them until things have stabilized
#
_stats_with_flags = [(Chi2DataVar, False),
                     (Chi2ModVar, False),
                     (Cash, True),
                     (LeastSq, True)]
_stats_no_flags = [s[0] for s in _stats_with_flags]


@pytest.mark.parametrize("statClass,flag", _stats_with_flags)
def test_data_plot_see_errorbar_warnings(caplog, statClass, flag):
    """Do we see the warning when expected - data plot?

    This looks for the 'The displayed errorbars have been supplied with
    the data or calculated using chi2xspecvar; the errors are not used in
    fits with <>' message. These are messages displayed to the Sherpa
    logger at the warning level, rather than using the warnings module,
    so the Sherpa capture_all_warnings test fixture does not come into
    play.

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
    flag : bool
        True if the warning should be created, False otherwise

    """

    d = example_data()
    plot = DataPlot()

    # Internal check: this test requires that either yerrorbars is set
    # to True, or not included, in the plot preferences. So check this
    # assumption.
    #
    prefs = plot.plot_prefs
    prefname = 'yerrorbars'
    assert (prefname not in prefs) or prefs[prefname]

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.prepare(d, stat)

    if flag:
        nwarn = 1
    else:
        nwarn = 0

    check_for_warning(caplog, nwarn, stat.name)


@pytest.mark.parametrize("statClass", _stats_no_flags)
def test_data_plot_no_errors_no_errorbar_warnings(caplog, statClass):
    """Should not see warnings when no error bars are drawn (See #621).

    This is a copy of test_data_plot_see_errorbar_warnings except that
    the 'yerrorbars' preference setting is set to False, so we should
    not be generating any warnings.

    Parameters
    ----------
    statClass : sherpa.stats.Stat instance

    """

    d = example_data()
    plot = DataPlot()

    prefs = plot.plot_prefs
    prefname = 'yerrorbars'
    prefs[prefname] = False

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.prepare(d, stat)

    check_for_warning(caplog, 0, stat.name)


@pytest.mark.parametrize("plotClass", [ResidPlot, RatioPlot])
@pytest.mark.parametrize("statClass,flag", _stats_with_flags)
def test_residstyle_plot_see_errorbar_warnings(caplog, plotClass, statClass, flag):
    """Do we see the warning when expected - residual/ratio plots?

    This looks for the 'The displayed errorbars have been supplied with
    the data or calculated using chi2xspecvar; the errors are not used in
    fits with <>' message. These are messages displayed to the Sherpa
    logger at the warning level, rather than using the warnings module,
    so the Sherpa capture_all_warnings test fixture does not come into
    play.

    Parameters
    ----------
    plotClass : {sherpa.plot.ResidPlot, sherpa.plot.RatioPlot}
        The plot to test.
    statClass : sherpa.stats.Stat instance
    flag : bool
        True if the warning should be created, False otherwise

    """

    d = example_data()
    m = example_model()
    plot = plotClass()

    # Internal check: this test requires that either yerrorbars is set
    # to True, or not included, in the plot preferences. So check this
    # assumption.
    #
    prefs = plot.plot_prefs
    prefname = 'yerrorbars'
    assert (prefname not in prefs) or prefs[prefname]

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.prepare(d, m, stat)

    if flag:
        nwarn = 1
    else:
        nwarn = 0

    check_for_warning(caplog, nwarn, stat.name)


@pytest.mark.parametrize("plotClass", [ResidPlot, RatioPlot])
@pytest.mark.parametrize("statClass", _stats_no_flags)
def test_residstyle_plot_no_errors_no_errorbar_warnings(caplog, plotClass, statClass):
    """Should see no warnings (see #621)

    This is a copy of test_residstyle_plot_see_errorbar_warnings
    but with yerrorbars preference set to False.

    Parameters
    ----------
    plotClass : {sherpa.plot.ResidPlot, sherpa.plot.RatioPlot}
        The plot to test.
    statClass : sherpa.stats.Stat instance

    """

    d = example_data()
    m = example_model()
    plot = plotClass()

    prefs = plot.plot_prefs
    prefname = 'yerrorbars'
    prefs[prefname] = False

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):
        plot.prepare(d, m, stat)

    check_for_warning(caplog, 0, stat.name)


@pytest.mark.parametrize("statClass,flag", _stats_with_flags)
def test_fit_plot_see_errorbar_warnings(caplog, statClass, flag):
    """Do we see the warning when expected - fit plot?

    This looks for the 'The displayed errorbars have been supplied with
    the data or calculated using chi2xspecvar; the errors are not used in
    fits with <>' message. These are messages displayed to the Sherpa
    logger at the warning level, rather than using the warnings module,
    so the Sherpa capture_all_warnings test fixture does not come into
    play.

    Parameters
    ----------
    stat : sherpa.stats.Stat instance
    flag : bool
        True if the warning should be created, False otherwise

    """

    d = example_data()
    m = example_model()

    dplot = DataPlot()
    mplot = ModelPlot()
    fplot = FitPlot()

    # Internal check: this test requires that either yerrorbars is set
    # to True, or not included, in the plot preferences. So check this
    # assumption.
    #
    # I am skipping model plot here, since it is assumed that there
    # are no errors on the model.
    #
    prefname = 'yerrorbars'
    for plot in [dplot, fplot]:
        prefs = plot.plot_prefs
        assert (prefname not in prefs) or prefs[prefname]

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):

        dplot.prepare(d, stat)
        mplot.prepare(d, m, stat)
        fplot.prepare(dplot, mplot)

    if flag:
        nwarn = 1
    else:
        nwarn = 0

    check_for_warning(caplog, nwarn, stat.name)


@pytest.mark.parametrize("plotClass", [ResidPlot, RatioPlot])
@pytest.mark.parametrize("statClass,flag", _stats_with_flags)
def test_fit_residstyle_plot_see_errorbar_warnings(caplog, plotClass, statClass, flag):
    """Do we see the warning when expected - fit + resid/ratio plot?

    This looks for the 'The displayed errorbars have been supplied with
    the data or calculated using chi2xspecvar; the errors are not used in
    fits with <>' message. These are messages displayed to the Sherpa
    logger at the warning level, rather than using the warnings module,
    so the Sherpa capture_all_warnings test fixture does not come into
    play.

    Parameters
    ----------
    plotClass : {sherpa.plot.ResidPlot, sherpa.plot.RatioPlot}
        The plot to test.
    statClass : sherpa.stats.Stat instance
    flag : bool
        True if the warning should be created, False otherwise

    Notes
    -----
    Is this an accurate example of how 'plot_fit_resid' is created?
    """

    d = example_data()
    m = example_model()

    dplot = DataPlot()
    mplot = ModelPlot()
    fplot = FitPlot()
    rplot = plotClass()

    jplot = JointPlot()

    # Internal check: this test requires that either yerrorbars is set
    # to True, or not included, in the plot preferences. So check this
    # assumption.
    #
    # I am skipping model plot here, since it is assumed that there
    # are no errors on the model.
    #
    prefname = 'yerrorbars'
    for plot in [dplot, fplot, rplot]:
        prefs = plot.plot_prefs
        assert (prefname not in prefs) or prefs[prefname]

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):

        dplot.prepare(d, stat)
        mplot.prepare(d, m, stat)
        fplot.prepare(dplot, mplot)

        rplot.prepare(d, m, stat)

        jplot.plottop(fplot)
        jplot.plotbot(rplot)

    if flag:
        nwarn = 2
    else:
        nwarn = 0

    check_for_warning(caplog, nwarn, stat.name)


@pytest.mark.parametrize("plotClass", [ResidPlot, RatioPlot])
@pytest.mark.parametrize("statClass", _stats_no_flags)
def test_fit_residstyle_plot_no_errors_no_errorbar_warnings(caplog, plotClass, statClass):
    """Should not see warnings when no error bars are drawn (See #621).

    This is a copy of test_fit_residstyle_plot_see_errorbar_warnings
    except that the 'yerrorbars' preference setting for all plots is
    'False'.

    Parameters
    ----------
    plotClass : {sherpa.plot.ResidPlot, sherpa.plot.RatioPlot}
        The plot to test.
    statClass : sherpa.stats.Stat instance

    Notes
    -----
    Is this an accurate example of how 'plot_fit_resid' is created?
    """

    d = example_data()
    m = example_model()

    dplot = DataPlot()
    mplot = ModelPlot()
    fplot = FitPlot()
    rplot = plotClass()

    jplot = JointPlot()

    prefname = 'yerrorbars'
    for plot in [dplot, rplot]:
        prefs = plot.plot_prefs
        prefs[prefname] = False

    stat = statClass()

    # Ensure that the logging is set to WARNING since there
    # appears to be some test that changes it to ERROR.
    #
    with caplog.at_level(logging.INFO, logger='sherpa'):

        dplot.prepare(d, stat)
        mplot.prepare(d, m, stat)
        fplot.prepare(dplot, mplot)

        rplot.prepare(d, m, stat)

        jplot.plottop(fplot)
        jplot.plotbot(rplot)

    check_for_warning(caplog, 0, stat.name)
