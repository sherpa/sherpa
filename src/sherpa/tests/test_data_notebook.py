#
#  Copyright (C) 2020, 2021, 2022, 2023
#  Smithsonian Astrophysical Observatory
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
from sherpa.plot.testing import check_full


def test_data1d(old_numpy_printing, all_plot_backends):
    d = data.Data1D("x x",
                    np.asarray([1, 3, 6]),
                    np.asarray([3, 4, 7]))
    r = d._repr_html_()

    check_full(r, 'Data1D Plot', test_other=[
        '<div class="sherpa-text-fallback">&lt;Data1D data set instance &#x27;x x&#x27;&gt;</div>',
        '<summary>Data1D Summary (5)</summary>',
        '<div class="dataname">Identifier</div><div class="dataval">x x</div>',
        '<div class="dataname">Number of bins</div><div class="dataval">3</div>',
        '<div class="dataname">Using</div><div class="dataval">1.0000-6.0000 x with 3 bins</div>',
        '<div class="dataname">X</div><div class="dataval">[1 3 6]</div>',
        '<div class="dataname">Y</div><div class="dataval">[3 4 7]</div>',
    ])



def test_data1d_errs(old_numpy_printing, all_plot_backends):
    d = data.Data1D("x x",
                    np.asarray([1, 3, 6]),
                    np.asarray([3, 4, 7]),
                    staterror=np.asarray([0.25, 0.52, 0.2]),
                    syserror=np.asarray([2, 3, 1]))
    r = d._repr_html_()
    check_full(r, 'Data1D Plot', test_other=[
        '<div class="sherpa-text-fallback">&lt;Data1D data set instance &#x27;x x&#x27;&gt;</div>',
        '<summary>Data1D Summary (7)</summary>',
        '<div class="dataname">Identifier</div><div class="dataval">x x</div>',
        '<div class="dataname">Number of bins</div><div class="dataval">3</div>',

        '<div class="dataname">Using</div><div class="dataval">1.0000-6.0000 x with 3 bins</div>',
        '<div class="dataname">X</div><div class="dataval">[1 3 6]</div>',
        '<div class="dataname">Y</div><div class="dataval">[3 4 7]</div>',
        '<div class="dataname">Statistical error</div><div class="dataval">[ 0.25  0.52  0.2 ]</div>',
        '<div class="dataname">Systematic error</div><div class="dataval">[2 3 1]</div>',
    ])


def test_data1dint(old_numpy_printing, all_plot_backends):
    d = data.Data1DInt("x x",
                       np.asarray([1, 3, 6]),
                       np.asarray([2, 6, 7]),
                       np.asarray([3, 4, 8]))
    r = d._repr_html_()

    check_full(r, 'Data1DInt Plot', test_other=[
        '<div class="sherpa-text-fallback">&lt;Data1DInt data set instance &#x27;x x&#x27;&gt;</div>',
        '<summary>Data1DInt Summary (6)</summary>',
        '<div class="dataname">Identifier</div><div class="dataval">x x</div>',
        '<div class="dataname">Number of bins</div><div class="dataval">3</div>',
        '<div class="dataname">Using</div><div class="dataval">1.0000-7.0000 x with 3 bins</div>',
        '<div class="dataname">XLO</div><div class="dataval">[1 3 6]</div>',
        '<div class="dataname">XHI</div><div class="dataval">[2 6 7]</div>',
        '<div class="dataname">Y</div><div class="dataval">[3 4 8]</div>',
    ])


# We still test with both backends to note if we ever do add
# pyplot support.
#
def test_data2d(old_numpy_printing, all_plot_backends):
    d = data.Data2D("x x",
                    np.asarray([1, 3, 1, 2, 3]),
                    np.asarray([2, 2, 4, 4, 4]),
                    np.asarray([4, 5, 6, 7, 8]))
    r = d._repr_html_()

    assert r is not None

    assert '<div class="sherpa-text-fallback">&lt;Data2D data set instance &#x27;x x&#x27;&gt;</div>' in r

    assert '<summary>Data2D Summary (5)</summary>' in r
    assert '<div class="dataname">Identifier</div><div class="dataval">x x</div>' in r
    assert '<div class="dataname">Number of bins</div><div class="dataval">5</div>' in r
    assert '<div class="dataname">X0</div><div class="dataval">[1 3 1 2 3]</div>' in r
    assert '<div class="dataname">X1</div><div class="dataval">[2 2 4 4 4]</div>' in r
    assert '<div class="dataname">Y</div><div class="dataval">[4 5 6 7 8]</div>' in r

    assert '<svg ' not in r
