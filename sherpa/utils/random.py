#
#  Copyright (C) 2023
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

"""
Basic routines for handling random numbers.
"""

import numpy as np

# Copy of sherpa.utils.SherpaFloat
#
SherpaFloat = np.float_


def poisson_noise(x):
    """Draw samples from a Poisson distribution.

    Parameters
    ----------
    x : scalar or array
       The expectation value for the distribution.

    Returns
    -------
    out : scalar or array
       A random realisation of the input array, drawn from
       the Poisson distribution, as a `SherpaFloat`.

    Notes
    -----
    The distribution is calculated by `np.poisson.poisson`.

    Examples
    --------
    >>> poisson_noise([10, 20, 5])
    array([ 13.,  21.,   6.])

    """

    x = np.asarray(x)

    # Using np.where() and indexing doesn't work with 0-d arrays, so
    # handle them separately
    if x.shape == ():
        x = SherpaFloat(x)
        if x <= 0.0:
            x = 0.0
        else:
            x = np.random.poisson(x)
        return SherpaFloat(x)

    x_out = np.zeros(x.shape, SherpaFloat)
    good = np.where(x > 0.0)
    x_out[good] = np.random.poisson(x[good])

    return x_out
