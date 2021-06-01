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

__all__ = ('intersperse', 'histogram_line')


def intersperse(a, b):
    '''Interleave two arrays a and b

    See https://stackoverflow.com/questions/5347065/interweaving-two-numpy-arrays#5347492
    for interleaving arrays.
    '''
    out = np.empty((a.size + b.size, ), dtype=a.dtype)
    out[0::2] = a
    out[1::2] = b
    return out


def histogram_line(xlo, xhi, y):
    '''Manually create x, y arrays for draw histrogram line

    Draw the data as a histogram, manually creating the lines
    from the low to high edge of each bin. In the case on non-consequtive
    bins, the line will have some nan values, so that disjoint lines are
    drawn.

    In matplotlib, an alternative would be to create RectanglePatches,
    one for each bin, but I don't want each bin to go down to 0. I do
    not find the existing drawstyle options to be sufficient.

    Parameters
    ----------
    xlo, xhi : array
        Lower and upper bin boundaries. Typically, ``xlo`` will contain the
        lower boundary and ``xhi`` the upper boundary, but this function can
        deal with situations where that is reversed. Both arrays have to be
        monotonically increasing or decreasing.
    y : array
        Dependent values for each histrogram bin_hi

    Returns
    -------
    x, y2 : array
        x and y arrays for plotting the histogram line.
    '''
    if (len(xlo) != len(xhi)) or (len(y) != len(xlo)):
        raise ValueError('All input arrays have to have the same length.')
    # Deal with xhi <-> xlo switches. Those can occor when converting
    # from energy to wavelength
    # Deal with reversed order. Can happen when converting from energy
    # to wavelength, or if input PHA is not ordered in increasing energy.
    # But is both are happening at the same time, need to switch twice, which
    # is a no-op. So, we get to use the elusive Python XOR operator.
    if (xlo[0] > xhi[0]) ^ (xhi[0] > xhi[-1]):
        xlo, xhi = xhi, xlo
    idxs, = np.where(xhi[:-1] != xlo[1:])
    x = intersperse(xlo, xhi)

    y2 = intersperse(y, y)

    if idxs.size > 0:
        idxs = 2 * (idxs + 1)
        nans = [np.nan] * idxs.size

        # ensure the arrays are floats so we can add nan values
        #
        x = np.insert(x.astype(np.float64), idxs, nans)
        y2 = np.insert(y2.astype(np.float64), idxs, nans)

    return x, y2
