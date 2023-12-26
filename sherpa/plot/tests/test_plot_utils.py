#
#  Copyright (C) 2021, 2023
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
import numpy as np

import pytest

from sherpa.plot.utils import histogram_line


def test_histrogram_line():
    '''Histrogram line manually puts together the (x,y) pairs to draw a
    line that looks list the "steps" option in matplotlib, but works better
    for our input formats and deal with missing or non-continuous bins by
    inserting nans. In this test, we check that the function works even if
    the input hi/lo arrays are interchanged to reversed. That can happen when
    converting from energy to wavelength and it's easier to just allow the
    flexibility in the plotting code than to enforce a certain ordering
    convention.
    '''
    xlow = np.arange(5, dtype=float)
    xhigh = xlow + 1
    xhigh[2] = 2.9
    y = np.arange(5) / .4

    for xlo, xhi in [(xlow, xhigh), (xhigh, xlow)]:
        x, y2 = histogram_line(xlo, xhi, y)
        assert np.all(np.ma.masked_invalid(x) ==
                      np.ma.masked_invalid([0., 1., 1., 2., 2., 2.9,
                                            np.nan, 3., 4., 4., 5.]))
        assert np.all(np.ma.masked_invalid(y2) ==
                      np.ma.masked_invalid([0.,  0.,  2.5,  2.5,  5.,  5.,
                                            np.nan, 7.5, 7.5, 10., 10.]))

        x, y2 = histogram_line(xlo[::-1], xhi[::-1], y)
        assert np.all(np.ma.masked_invalid(x) ==
                      np.ma.masked_invalid([5., 4., 4., 3., np.nan, 2.9,
                                            2., 2., 1., 1., 0.]))
        assert np.all(np.ma.masked_invalid(y2) ==
                      np.ma.masked_invalid([0.,  0.,  2.5,  2.5, np.nan,  5.,
                                            5., 7.5, 7.5, 10., 10.]))


def test_histogram_line_integers_no_gap():
    """Test based on issue #1838

    This is just a regression test to document the current behavior.
    """

    x1 = np.asarray([10, 20, 30])
    x2 = np.asarray([20, 30, 40])
    y = np.asarray([12, 3, 15])

    # based on https://github.com/numpy/numpy/issues/17325
    assert np.issubdtype(x1.dtype, np.integer)
    assert np.issubdtype(x2.dtype, np.integer)
    assert np.issubdtype(y.dtype, np.integer)

    xh, yh = histogram_line(x1, x2, y)
    assert np.issubdtype(xh.dtype, np.integer)
    assert np.issubdtype(yh.dtype, np.integer)

    assert xh == pytest.approx([10, 20, 20, 30, 30, 40])
    assert yh == pytest.approx([12, 12, 3, 3, 15, 15])
    

def test_histogram_line_floats_no_gap():
    """Test based on issue #1838

    This is just a regression test to document the current behavior.
    """

    x1 = np.asarray([1, 2.5, 3])
    x2 = np.asarray([2.5, 3, 4])
    y = np.asarray([12, 3, 15])

    # based on https://github.com/numpy/numpy/issues/17325
    assert np.issubdtype(x1.dtype, np.floating)
    assert np.issubdtype(x2.dtype, np.floating)
    assert np.issubdtype(y.dtype, np.integer)

    xh, yh = histogram_line(x1, x2, y)
    assert np.issubdtype(xh.dtype, np.floating)
    assert np.issubdtype(yh.dtype, np.integer)

    assert xh == pytest.approx([1, 2.5, 2.5, 3, 3, 4])
    assert yh == pytest.approx([12, 12, 3, 3, 15, 15])
    

def test_histogram_line_integers_single_gap():
    """Test based on issue #1838

    This is just a regression test to document the current behavior.
    """

    x1 = np.asarray([10, 20, 30])
    x2 = np.asarray([20, 25, 40])
    y = np.asarray([12, 3, 15])

    # based on https://github.com/numpy/numpy/issues/17325
    assert np.issubdtype(x1.dtype, np.integer)
    assert np.issubdtype(x2.dtype, np.integer)
    assert np.issubdtype(y.dtype, np.integer)

    xh, yh = histogram_line(x1, x2, y)
    assert np.issubdtype(xh.dtype, np.floating)
    assert np.issubdtype(yh.dtype, np.floating)

    good = np.isfinite(xh)
    assert good == pytest.approx(np.isfinite(yh))
    assert good == pytest.approx([1, 1, 1, 1, 0, 1, 1])

    assert xh[good] == pytest.approx([10, 20, 20, 25, 30, 40])
    assert yh[good] == pytest.approx([12, 12, 3, 3, 15, 15])


def test_histogram_line_floats_single_gap():
    """Test based on issue #1838

    This is just a regression test to document the current behavior.
    """

    x1 = np.asarray([1, 2, 3])
    x2 = np.asarray([2, 2.5, 4])
    y = np.asarray([12, 3, 15])

    # based on https://github.com/numpy/numpy/issues/17325
    assert np.issubdtype(x1.dtype, np.integer)
    assert np.issubdtype(x2.dtype, np.floating)
    assert np.issubdtype(y.dtype, np.integer)

    xh, yh = histogram_line(x1, x2, y)
    assert np.issubdtype(xh.dtype, np.floating)
    assert np.issubdtype(yh.dtype, np.floating)

    good = np.isfinite(xh)
    assert good == pytest.approx(np.isfinite(yh))
    assert good == pytest.approx([1, 1, 1, 1, 0, 1, 1])

    assert xh[good] == pytest.approx([1, 2, 2, 2.5, 3, 4])
    assert yh[good] == pytest.approx([12, 12, 3, 3, 15, 15])
    
