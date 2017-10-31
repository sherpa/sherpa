#
#  Copyright (C) 2017  Smithsonian Astrophysical Observatory
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
This adds additional tests to test_optical but uses pytest.
"""

import numpy as np
from numpy.testing import assert_allclose

import sherpa.astro.optical as models


def test_emissiongaussian_skew1():
    mdl = models.EmissionGaussian()

    # select an x-axis range which leads to non-zero
    # values in multiple bins
    x = np.arange(4995.0, 5006.0, 1.0)
    y = mdl(x)

    # This is symmetrical since the grid is symmetrical
    # about the pos parameter (default is 5000).
    #
    expected = np.zeros(11)
    expected[3:8] = [0.01045223, 0.2078923, 0.5632678, 0.2078923,
                     0.01045223]

    # The same absolute tolerance is used as in test_optical.py;
    # perhaps a smaller value could be used?
    assert_allclose(expected, y, atol=1.0e-4, rtol=0)


def test_emissiongaussian_skew3():
    # This currently applies several expected properties of the
    # result (i.e. is it peaked at the right place and skewed)
    # and does not do an explicit check of the expected values as
    # it is not clear what they are.
    #

    mdl = models.EmissionGaussian()
    mdl.skew = 3

    # select an x-axis range which leads to non-zero
    # values in multiple bins
    x = np.arange(4995.0, 5006.0, 1.0)
    y = mdl(x)

    # sanity check
    assert (y >= 0.0).all()

    # Check subset is working (i.e. that the outer bins
    # are zero).
    #
    idx0 = np.where((x < 4998) | (x > 5002))
    assert_allclose(np.zeros(6), y[idx0])

    # This should be asymmetrical about x=5000, with x>5000
    # being larger than x<5000, and the peak at 5000.
    #
    midval = y[5]
    assert (midval > y[x != 5000]).all()

    lsum = y[x < 5000].sum()
    hsum = y[x > 5000].sum()
    assert lsum < hsum


def test_emissiongaussian_skew025():
    # Since the skew is < 1 the lower-x values should now be
    # higher.
    #

    mdl = models.EmissionGaussian()
    mdl.skew = 0.25

    # select an x-axis range which leads to non-zero
    # values in multiple bins
    x = np.arange(4995.0, 5006.0, 1.0)
    y = mdl(x)

    # sanity check
    assert (y >= 0.0).all()

    # Check subset is working (i.e. that the outer bins
    # are zero).
    #
    idx0 = np.where((x < 4998) | (x > 5002))
    assert_allclose(np.zeros(6), y[idx0])

    # This should be asymmetrical about x=5000, with x>5000
    # being smaller than x<5000, and the peak at 5000.
    #
    midval = y[5]
    assert (midval > y[x != 5000]).all()

    lsum = y[x < 5000].sum()
    hsum = y[x > 5000].sum()
    assert lsum > hsum
