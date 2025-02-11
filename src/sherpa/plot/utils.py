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
'''Helper functions for plotting

'''
import numpy as np


__all__ = ('histogram_line', )


def histogram_line(xlo, xhi, y):
    '''Manually create x, y arrays to replicate a histogram

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

    # Deal with reversed order. Can happen when converting from energy
    # to wavelength, or if input PHA is not ordered in increasing energy.
    # But if both are happening at the same time, need to switch twice, which
    # is a no-op. So, we get to use the elusive Python XOR operator.
    if (xlo[0] > xhi[0]) ^ (xhi[0] > xhi[-1]):
        xlo, xhi = xhi, xlo

    # Combine the edges and replicate the y array. This ensures that x
    # has the "correct" type (e.g. if xlo is int but xhi is float the
    # result will be float), which was the cause of issue #1838.
    #
    x = np.vstack((xlo, xhi)).T.flatten()
    y2 = np.vstack((y, y)).T.flatten()

    # Where are the edges?
    idxs, = np.where(xhi[:-1] != xlo[1:])
    nedges = idxs.size
    if nedges == 0:
        return x, y2

    # nan values need to be added where the edges are in the
    # duplicated array
    nanidxs = 2 * (idxs + 1)
    nans = [np.nan] * nedges

    # ensure the arrays are floats so we can add nan values
    #
    x = np.insert(x.astype(np.float64), nanidxs, nans)
    y2 = np.insert(y2.astype(np.float64), nanidxs, nans)

    return x, y2
