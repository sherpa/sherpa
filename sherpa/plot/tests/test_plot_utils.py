#
#  Copyright (C) 2021
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

from sherpa.plot.utils import histogram_line


def test_histrogram_line():
    '''Histrogram line manually puts together the (x,y) pairs to draw a
    line that looks list the "steps" option in matplotlib, but works better
    for our input formats and deal with missing or non-continuous bins by
    inserting nans. In this test, we check that the fuction works even if
    the input hi/lo arrays are intercahnged to reversed. That can happen when
    converting from energy to wavelength and it's eaiser to just allow the
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
