#
#  Copyright (C) 2007, 2015, 2016, 2018, 2019, 2020, 2023
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

"""Catch all for optimiser tests"""

import numpy as np

import pytest

from sherpa.optmethods.opt import SimplexRandom


def rosenbrock(x):
    """Rosenbrock(1, 1) = 0
    """
    val = np.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0, axis=0)
    return val


SEED = 2354


def test_is_simplexbase_repeatable():
    """The SimplexBase* classes uses RNG, so can we make it repeatable?

    This is based on the 'if __name__ == "__main__"' section of
    opt.py, (see issue #1737 in case the code gets cleaned up) and
    only one of the classes based on SimplexBase is used.

    The point of the test is to pick a npop size large enough so that
    the Simplex class seeds the simplex with randomly-distributed
    parameter values. For the SimpexRandom case this is npop>1.

    """

    x0 = np.array([-1.2, 1.0])
    xmin = [-1000, -1000]
    xmax = [1000, 1000]
    simp = SimplexRandom(func=rosenbrock, npop=4, xpar=x0, xmin=xmin,
                         xmax=xmax, step=x0 + 1.2, seed=SEED,
                         factor=10)

    # After the first line, which is the x0 array + function value,
    # the next three rows have 2 random numbers and then the function
    # value at that coordinate, and then sorted on the last column.
    # So this tests how repeatable the RNG us.
    #
    expected = np.asarray([[-1.2, 1.0, 24.2],
                           [4.21466190e-01, -1.72015071e+00, 3.60493288e+02],
                           [-3.66908476e+00, 6.96971374e+00, 4.23701601e+03],
                           [-6.71378027e-01, 1.07842741e+01, 1.06809687e+04]])
    assert simp.simplex == pytest.approx(expected)
