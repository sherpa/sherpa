#
# Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

"""Very-basic tests of the HTML representation of objects.

"""

import pytest

from sherpa import fit
from sherpa.data import Data1D, DataSimulFit
from sherpa.models import SimulFitModel
from sherpa.models.basic import Const1D
from sherpa.optmethods import LevMar, NelderMead
from sherpa.stats import Chi2, LeastSq


def test_statinfo_basic():
    # Leave model argument as None as isn't used
    s = fit.StatInfoResults('chi2foo', 1.2, 12, None, 10)
    r = s._repr_html_()

    assert r is not None

    assert '<summary>Statistics summary (4)</summary>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2foo</div>' in r
    assert '<div class="dataname">Value</div><div class="dataval">1.2</div>' in r
    assert '<div class="dataname">Number of points</div><div class="dataval">12</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">10</div>' in r


def test_statinfo_chisq():
    """Fake rstat and qval"""
    s = fit.StatInfoResults('chi2foo', 1.2, 12, None, 10)
    s.qval= 0.9
    s.rstat = 1.1
    r = s._repr_html_()

    assert r is not None

    assert '<summary>Statistics summary (6)</summary>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2foo</div>' in r
    assert '<div class="dataname">Value</div><div class="dataval">1.2</div>' in r
    assert '<div class="dataname">Number of points</div><div class="dataval">12</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">10</div>' in r
    assert '<div class="dataname">Reduced statistic</div><div class="dataval">1.1</div>' in r
    assert '<div class="dataname">Probability (Q-value)</div><div class="dataval">0.9</div>' in r


def test_statinfo_single():
    s = fit.StatInfoResults('chi2foo', 1.2, 12, None, 10)
    s.ids = ["foo"]
    r = s._repr_html_()

    assert r is not None

    assert '<summary>Statistics summary (5)</summary>' in r
    assert '<div class="dataname">Dataset</div><div class="dataval">foo</div>' in r

    assert '<div class="dataname">Statistic</div><div class="dataval">chi2foo</div>' in r
    assert '<div class="dataname">Value</div><div class="dataval">1.2</div>' in r
    assert '<div class="dataname">Number of points</div><div class="dataval">12</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">10</div>' in r


def test_statinfo_single_bkg():
    s = fit.StatInfoResults('chi2foo', 1.2, 12, None, 10)
    s.ids = ["foo"]
    s.bkg_ids = [2]
    r = s._repr_html_()

    assert r is not None

    assert '<summary>Statistics summary (6)</summary>' in r
    assert '<div class="dataname">Dataset</div><div class="dataval">foo</div>' in r
    assert '<div class="dataname">Background</div><div class="dataval">2</div>' in r

    assert '<div class="dataname">Statistic</div><div class="dataval">chi2foo</div>' in r
    assert '<div class="dataname">Value</div><div class="dataval">1.2</div>' in r
    assert '<div class="dataname">Number of points</div><div class="dataval">12</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">10</div>' in r


def test_statinfo_multi():

    s = fit.StatInfoResults('chi2foo', 1.2, 12, None, 10)
    s.ids = [1, 2]
    s.bkg_ids = ["x", "y"]

    r = s._repr_html_()

    assert r is not None

    assert '<summary>Statistics summary (6)</summary>' in r
    assert '<div class="dataname">Datasets</div><div class="dataval">[1, 2]</div>' in r
    assert '<div class="dataname">Backgrounds</div><div class="dataval">[\'x\', \'y\']</div>' in r

    assert '<div class="dataname">Statistic</div><div class="dataval">chi2foo</div>' in r
    assert '<div class="dataname">Value</div><div class="dataval">1.2</div>' in r
    assert '<div class="dataname">Number of points</div><div class="dataval">12</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">10</div>' in r


@pytest.mark.parametrize('method', [LevMar, NelderMead])
def test_fitresults(method):
    d = Data1D('dx', [1, 2, 3], [4, 2, 2])
    m = Const1D()
    m.c0 = 3
    fr = fit.Fit(d, m, method=method(), stat=LeastSq()).fit()
    r = fr._repr_html_()

    assert r is not None

    assert '<summary>Fit parameters</summary>' in r
    assert '<summary>Summary (8)' in r
    assert '<td>const1d.c0</td>' in r

    assert '<div class="dataname">Method</div><div class="dataval">{}</div>'.format(fr.methodname) in r
    assert '<div class="dataname">Statistic</div><div class="dataval">leastsq</div>' in r

    assert '<div class="dataname">&#916; statistic</div><div class="dataval">0.333333</div>' in r
    assert '<div class="dataname">Number of data points</div><div class="dataval">3</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">2</div>' in r


def test_fitresults_named():
    d = Data1D('dx', [1, 2, 3], [4, 2, 2])
    m = Const1D()
    m.c0 = 3
    fr = fit.Fit(d, m, method=LevMar(), stat=LeastSq()).fit()
    fr.datasets = [1]
    r = fr._repr_html_()

    assert r is not None

    assert '<summary>Fit parameters</summary>' in r
    assert '<summary>Summary (9)' in r
    assert '<td>const1d.c0</td>' in r

    assert '<div class="dataname">Dataset</div><div class="dataval">1</div>' in r
    assert '<div class="dataname">Method</div><div class="dataval">{}</div>'.format(fr.methodname) in r
    assert '<div class="dataname">Statistic</div><div class="dataval">leastsq</div>' in r

    assert '<div class="dataname">&#916; statistic</div><div class="dataval">0.333333</div>' in r
    assert '<div class="dataname">Number of data points</div><div class="dataval">3</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">2</div>' in r


def test_fitresults_chisq():
    d = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 1.4, 1.4])
    m = Const1D()
    m.c0 = 3
    fr = fit.Fit(d, m, method=LevMar(), stat=Chi2()).fit()
    r = fr._repr_html_()

    assert r is not None

    assert '<summary>Fit parameters</summary>' in r
    assert '<summary>Summary (10)</summary>' in r
    assert '<td>const1d.c0</td>' in r

    assert '<div class="dataname">Method</div><div class="dataval">{}</div>'.format(fr.methodname) in r
    assert '<div class="dataname">Final statistic</div><div class="dataval">1.65289</div>' in r

    assert '<div class="dataname">Reduced statistic</div><div class="dataval">0.826446</div>' in r
    assert '<div class="dataname">Probability (Q-value)</div><div class="dataval">0.437602</div>' in r

    assert '<div class="dataname">Number of data points</div><div class="dataval">3</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">2</div>' in r


@pytest.mark.parametrize('method', [LevMar, NelderMead])
def test_fitresults_multi(method):
    """Fit multiple datasets"""

    d1 = Data1D('dx', [1, 2, 3], [4, 2, 2])
    d2 = Data1D('dx', [4, 5, 6, 10], [4, 4, 2, 4])
    d = DataSimulFit('combined', (d1, d2))

    m1 = Const1D()
    m1.c0 = 3
    m = SimulFitModel('silly', (m1, m1))

    fr = fit.Fit(d, m, method=method(), stat=LeastSq()).fit()
    fr.datasets = ['ddx', 'ddy']
    r = fr._repr_html_()

    assert r is not None

    assert '<summary>Summary (9)</summary>' in r
    assert '<td>const1d.c0</td>' in r

    assert '<div class="dataname">Datasets</div><div class="dataval">ddx,ddy</div>' in r
    assert '<div class="dataname">Method</div><div class="dataval">{}</div>'.format(fr.methodname) in r
    assert '<div class="dataname">Statistic</div><div class="dataval">leastsq</div>' in r

    assert '<div class="dataname">&#916; statistic</div><div class="dataval">0.142857</div>' in r
    assert '<div class="dataname">Number of data points</div><div class="dataval">7</div>' in r
    assert '<div class="dataname">Degrees of freedom</div><div class="dataval">6</div>' in r


def test_errresults():
    d = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    m = Const1D()
    m.c0 = 3
    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()
    r = er._repr_html_()

    assert r is not None

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (2)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>   -0.562226</td><td>    0.562226</td></tr>' in r




def test_errresults_limits_none():
    """Missing an error limit"""
    d = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    m = Const1D()
    m.c0 = 3
    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()

    # perhaps should just fake this instead?
    assert len(er.parmins) == 1
    er.parmins = (None, )

    r = er._repr_html_()

    assert r is not None

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (2)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>-----</td><td>    0.562226</td></tr>' in r


def test_errresults_limits_interval():
    """Missing an error limit"""
    d = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    m = Const1D()
    m.c0 = 3
    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()

    # perhaps should just fake this instead?
    assert len(er.parmaxes) == 1
    er.parmaxes = ([0.1, 0.2], )

    r = er._repr_html_()

    assert r is not None

    print(r)

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (2)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>   -0.562226</td><td>(1.000000e-01, 2.000000e-01)</td></tr>' in r


def test_errresults_named():
    d = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    m = Const1D()
    m.c0 = 3
    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()
    er.datasets = [2]
    r = er._repr_html_()

    assert r is not None

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (3)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Dataset</div><div class="dataval">2</div>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>   -0.562226</td><td>    0.562226</td></tr>' in r


def test_errresults_multi():
    d1 = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    d2 = Data1D('dx', [10, 11, 12, 13], [4, 4, 2, 4], [0.8, 1.1, 1.1, 0.9])
    d = DataSimulFit('combined', (d1, d2))

    m1 = Const1D()
    m1.c0 = 3
    m = SimulFitModel('silly', (m1, m1))

    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()
    r = er._repr_html_()

    assert r is not None

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (2)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>   -0.362415</td><td>    0.362415</td></tr>' in r


def test_errresults_multi_named():
    d1 = Data1D('dx', [1, 2, 3], [4, 2, 2], [1.2, 0.9, 0.9])
    d2 = Data1D('dx', [10, 11, 12, 13], [4, 4, 2, 4], [0.8, 1.1, 1.1, 0.9])
    d = DataSimulFit('combined', (d1, d2))

    m1 = Const1D()
    m1.c0 = 3
    m = SimulFitModel('silly', (m1, m1))

    f = fit.Fit(d, m, stat=Chi2())
    er = f.est_errors()
    er.datasets = ['dx', 'dy']
    r = er._repr_html_()

    assert r is not None

    assert '<summary>covariance 1&#963; (68.2689%) bounds</summary>' in r
    assert '<summary>Summary (3)' in r
    assert '<td>const1d.c0</td>' in r
    assert '<div class="dataname">Datasets</div><div class="dataval">dx,dy</div>' in r
    assert '<div class="dataname">Fitting Method</div><div class="dataval">levmar</div>' in r
    assert '<div class="dataname">Statistic</div><div class="dataval">chi2</div>' in r

    assert '<tr><td>const1d.c0</td><td>           3</td><td>   -0.362415</td><td>    0.362415</td></tr>' in r
