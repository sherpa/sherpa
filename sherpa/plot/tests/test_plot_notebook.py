#
# Copyright (C) 2020, 2021  Smithsonian Astrophysical Observatory
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

The testing depends on the backend, but for now this is all
kept with each test (i.e. the checks runin a test depend
on the backend).

"""

import numpy as np

from sherpa import plot
from sherpa.data import Data1D, Data2D
from sherpa.fit import Fit
from sherpa.models.basic import Const1D, Gauss2D, Polynom1D
from sherpa.stats import Chi2


def check_empty(r, summary, nsummary=0):
    """Is this an 'empty' response?"""

    if plot.backend.name == 'pylab':
        assert r is None
        return

    assert r is not None

    assert '<summary>{} ({})</summary>'.format(summary, nsummary) in r


def check_full(r, summary, label, title, nsummary=0):
    """Is this a 'full' response?"""

    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>{}</summary>'.format(summary) in r
        assert '<svg ' in r
        return

    assert '<summary>{} ({})</summary>'.format(summary, nsummary) in r
    assert '<div class="dataval">{}</div>'.format(label) in r
    assert '<div class="dataval">{}</div>'.format(title) in r


def test_histogram(override_plot_backend):
    p = plot.HistogramPlot()
    r = p._repr_html_()
    check_empty(r, 'HistogramPlot', nsummary=6)

    p.xlo = np.asarray([1, 2, 3])
    p.xhi = np.asarray([2, 2.5, 4])
    p.y = np.asarray([4, 2, 2])
    p.xlabel = 'X <sup>2</sup>'
    p.ylabel = 'Y'
    p.title = 'Title String'
    r = p._repr_html_()
    check_full(r, 'HistogramPlot', 'X <sup>2</sup>',
               'Title String', nsummary=6)


def test_pdfplot(override_plot_backend):
    p = plot.PDFPlot()
    r = p._repr_html_()
    check_empty(r, 'PDFPlot', nsummary=7)

    x = np.arange(10) / 2
    p.prepare(x)
    r = p._repr_html_()
    check_full(r, 'PDFPlot', 'probability density', 'PDF: x',
               nsummary=7)


def test_cdfplot(override_plot_backend):
    p = plot.CDFPlot()
    r = p._repr_html_()
    check_empty(r, 'CDFPlot', nsummary=9)

    x = np.arange(10) / 2
    p.prepare(x)
    r = p._repr_html_()
    # TODO: the <= should be converted to &lt;=
    check_full(r, 'CDFPlot', 'p(<=x)', 'CDF: x', nsummary=9)


def test_lrhist(override_plot_backend):
    p = plot.LRHistogram()
    r = p._repr_html_()
    check_empty(r, 'LRHistogram', nsummary=8)

    # no idea what sensible values should be
    x = np.arange(10) / 2
    p.prepare(x, 5, 2, 0.1, 0.2)
    r = p._repr_html_()
    check_full(r, 'LRHistogram', 'Frequency',
               'Likelihood Ratio Distribution', nsummary=8)


def test_data(override_plot_backend):
    p = plot.DataPlot()
    r = p._repr_html_()
    check_full(r, 'DataPlot', 'None', 'None', nsummary=7)  # NOT empty

    x = np.arange(5, 8, 0.5)
    y = np.ones(x.size)
    dr = y / 0.1
    d = Data1D('n n', x, y, staterror=dr)
    p.prepare(d)
    r = p._repr_html_()
    check_full(r, 'DataPlot', 'y', 'n n', nsummary=7)


def test_datacontour(override_plot_backend):
    p = plot.DataContour()
    r = p._repr_html_()
    check_empty(r, 'DataContour', nsummary=7)

    x0 = np.asarray([[10, 20, 30], [10, 20, 30]])
    x1 = np.asarray([[5, 5, 5], [15, 15, 15]])
    y = np.arange(x0.size)
    d = Data2D('n n', x0.flatten(), x1.flatten(), y.flatten())

    p.prepare(d)
    r = p._repr_html_()
    check_full(r, 'DataContour', 'x1', 'n n', nsummary=7)


def test_model(override_plot_backend):
    p = plot.ModelPlot()
    r = p._repr_html_()
    check_empty(r, 'ModelPlot', nsummary=7)

    x = np.arange(5, 8, 0.5)
    y = np.ones(x.size)
    d = Data1D('n n', x, y)

    m = Const1D()

    p.prepare(d, m)
    r = p._repr_html_()
    check_full(r, 'ModelPlot', 'y', 'Model', nsummary=7)


def test_modelcontour(override_plot_backend):
    p = plot.ModelContour()
    r = p._repr_html_()
    check_empty(r, 'ModelContour', nsummary=7)

    x0 = np.asarray([[10, 20, 30], [10, 20, 30]])
    x1 = np.asarray([[5, 5, 5], [15, 15, 15]])
    y = np.arange(x0.size)
    d = Data2D('n n', x0.flatten(), x1.flatten(), y.flatten())

    # need somethig we can contour
    m = Gauss2D()
    m.xpos = 21
    m.ypos = 9
    m.fwhm = 5

    # we need a stat argument (as of 4.12.1) but it isn't used so
    # send in None
    p.prepare(d, m, stat=None)
    r = p._repr_html_()
    check_full(r, 'ModelContour', 'x1', 'Model', nsummary=7)


def test_fit(override_plot_backend):
    p = plot.FitPlot()
    r = p._repr_html_()
    assert r is None  # note: always None

    x = np.arange(5, 8, 0.5)
    y = np.ones(x.size)
    d = Data1D('n n', x, y)

    m = Const1D()

    dplot = plot.DataPlot()
    dplot.prepare(d)

    mplot = plot.ModelPlot()
    mplot.prepare(d, m)

    p.prepare(dplot, mplot)
    r = p._repr_html_()

    # different to previous checks
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>FitPlot</summary>' in r
        assert '<svg ' in r
        return

    assert '<summary>DataPlot (' in r
    assert '<summary>ModelPlot (' in r
    assert '<div class="dataval">n n</div>' in r
    assert '<div class="dataval">Model</div>' in r


def test_fitcontour(override_plot_backend):
    p = plot.FitContour()
    r = p._repr_html_()
    assert r is None  # note: always None

    x0 = np.asarray([[10, 20, 30], [10, 20, 30]])
    x1 = np.asarray([[5, 5, 5], [15, 15, 15]])
    y = np.arange(x0.size)
    d = Data2D('n n', x0.flatten(), x1.flatten(), y.flatten())

    # need somethig we can contour
    m = Gauss2D()
    m.xpos = 21
    m.ypos = 9
    m.fwhm = 5

    dplot = plot.DataContour()
    dplot.prepare(d)

    mplot = plot.ModelContour()
    mplot.prepare(d, m, stat=None)

    p.prepare(dplot, mplot)
    r = p._repr_html_()

    # different to previous checks
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>FitContour</summary>' in r
        assert '<svg ' in r
        return

    assert '<summary>DataContour (' in r
    assert '<summary>ModelContour (' in r
    assert '<div class="dataval">n n</div>' in r
    assert '<div class="dataval">Model</div>' in r


def test_intproj(old_numpy_printing, override_plot_backend):
    p = plot.IntervalProjection()
    r = p._repr_html_()

    check_empty(r, 'IntervalProjection', nsummary=8)

    x = np.arange(5, 8, 0.5)
    y = np.asarray([2, 3, 4, 5, 4, 3])
    dy = y / 2
    d = Data1D('n n', x, y, staterror=dy)

    m = Const1D()

    fit = Fit(d, m, stat=Chi2())
    fr = fit.fit()
    assert fr.succeeded

    p.prepare(min=1, max=6, nloop=10)
    p.calc(fit, m.c0)

    r = p._repr_html_()
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>IntervalProjection</summary>' in r
        assert '<svg ' in r
        return

    assert '<summary>IntervalProjection (8)</summary>' in r

    assert '<div class="dataname">x</div><div class="dataval">[ 1.        1.555556  2.111111  2.666667  3.222222  3.777778  4.333333  4.888889\n  5.444444  6.      ]</div>' in r
    assert '<div class="dataname">nloop</div><div class="dataval">10</div>' in r


def test_regproj(old_numpy_printing, override_plot_backend):
    p = plot.RegionProjection()
    r = p._repr_html_()

    check_empty(r, 'RegionProjection', nsummary=13)

    x = np.arange(5, 8, 0.5)
    y = np.asarray([2, 3, 4, 5, 4, 3])
    dy = y / 2
    d = Data1D('n n', x, y, staterror=dy)

    m = Polynom1D()
    m.c1.thaw()

    fit = Fit(d, m, stat=Chi2())
    fr = fit.fit()
    assert fr.succeeded

    p.prepare(min=(-2, -1), max=(2, 2), nloop=(10, 20))
    p.calc(fit, m.c0, m.c1)

    r = p._repr_html_()
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>RegionProjection</summary>' in r
        assert '<svg ' in r
        return

    print(r)
    assert '<summary>RegionProjection (13)</summary>' in r

    assert '<div class="dataname">parval0</div><div class="dataval">-0.5315772076542427</div>' in r
    assert '<div class="dataname">parval1</div><div class="dataval">0.5854611101216837</div>' in r
    assert '<div class="dataname">sigma</div><div class="dataval">(1, 2, 3)</div>' in r

    assert '<div class="dataname">y</div><div class="dataval">[ 306.854444  282.795953  259.744431  237.699877  216.662291  196.631674\n' in r

    assert '<div class="dataname">levels</div><div class="dataval">[  3.606863   7.491188  13.140272]</div>' in r
    assert '<div class="dataname">min</div><div class="dataval">[-2, -1]</div>' in r
    assert '<div class="dataname">max</div><div class="dataval">[2, 2]</div>' in r

    assert '<div class="dataname">nloop</div><div class="dataval">(10, 20)</div>' in r
