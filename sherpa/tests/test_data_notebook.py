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

import numpy as np

from sherpa import data
from sherpa import plot


def check(r, summary, title, vals, nsummary):
    assert r is not None

    if plot.backend.name == 'pylab':
        assert '<summary>{} Plot</summary>'.format(summary) in r
        assert '<svg ' in r
        return

    assert '<summary>{} Summary ({})</summary>'.format(summary, nsummary) in r
    assert '<div class="dataval">{}</div>'.format(title) in r
    assert '<div class="dataval">{}</div>'.format(vals) in r
    assert '<svg ' not in r


def test_data1d(override_plot_backend):
    d = data.Data1D("x x",
                    np.asarray([1, 3, 6]),
                    np.asarray([3, 4, 7]))
    r = d._repr_html_()
    check(r, 'Data1D', 'x x', '[3 4 7]', nsummary=5)


def test_data1d_errs(old_numpy_printing, override_plot_backend):
    d = data.Data1D("x x",
                    np.asarray([1, 3, 6]),
                    np.asarray([3, 4, 7]),
                    staterror=np.asarray([0.25, 0.52, 0.2]),
                    syserror=np.asarray([2, 3, 1]))
    r = d._repr_html_()
    check(r, 'Data1D', 'x x', '[3 4 7]', nsummary=7)

    if plot.backend.name != 'dummy':
        return

    assert '<div class="dataname">Identifier</div><div class="dataval">x x</div>' in r
    assert '<div class="dataname">Number of bins</div><div class="dataval">3</div>' in r

    assert '<div class="dataname">Using</div><div class="dataval">1.0000-6.0000 x with 3 bins</div>' in r
    assert '<div class="dataname">X</div><div class="dataval">[1 3 6]</div>' in r
    assert '<div class="dataname">Y</div><div class="dataval">[3 4 7]</div>' in r
    assert '<div class="dataname">Statistical error</div><div class="dataval">[ 0.25  0.52  0.2 ]</div>' in r
    assert '<div class="dataname">Systematic error</div><div class="dataval">[2 3 1]</div>' in r


def test_data1dint(override_plot_backend):
    d = data.Data1DInt("x x",
                       np.asarray([1, 3, 6]),
                       np.asarray([2, 6, 7]),
                       np.asarray([3, 4, 8]))
    r = d._repr_html_()
    check(r, 'Data1DInt', 'x x', '[3 4 8]', nsummary=6)


# We still test with both backends to note if we ever do add
# pyplot support.
#
def test_data2d(override_plot_backend):
    d = data.Data2D("x x",
                    np.asarray([1, 3, 1, 2, 3]),
                    np.asarray([2, 2, 4, 4, 4]),
                    np.asarray([4, 5, 6, 7, 8]))
    r = d._repr_html_()

    assert r is not None

    assert '<summary>{} Summary (5)</summary>'.format('Data2D') in r
    assert '<div class="dataval">{}</div>'.format('x x') in r
    assert '<div class="dataval">{}</div>'.format('[2 2 4 4 4]') in r
    assert '<svg ' not in r
