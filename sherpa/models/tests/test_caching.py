#
#  Copyright (C) 2025
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
"""Tests related to caching model evaluations

More tests with 1D data are in test_model.py and might be moved here at a later point.
"""
import numpy as np
import pytest

from sherpa.data import Data2D
from sherpa.models import Polynom2D, Gauss2D
from sherpa.astro.models import Beta2D, Lorentz2D
from sherpa.models.model import UnaryOpModel, BinaryOpModel, RegridWrappedModel


@pytest.mark.parametrize("modelclass", [Polynom2D, Gauss2D, Beta2D, Lorentz2D])
def test_default_cache_size(modelclass):
    """Check the default cache size for some models."""

    mdl = modelclass()
    assert mdl.cache == 0


def test_evaluate_no_cache2d():
    """Check we can turn off caching: 2d"""

    x0grid = np.array([2, 10, 1.5])
    x1grid = np.array([-2, 0, 3.])

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 0
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x0grid, x1grid) == pytest.approx(expected)
    assert len(mdl._cache) == 0

    mdl.cx1 = 5
    mdl.cy1 = 2

    expected = 1 + 2 * x0grid + 5 * x1grid
    assert mdl(x1grid, x0grid) == pytest.approx(expected)
    assert len(mdl._cache) == 0


def test_evaluate_cache2d():
    """Check we run with caching on: 1d"""

    x0grid = np.array([2, 10, 1.5])
    x1grid = np.array([-2, 0, 3.])

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 5
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x1grid, x0grid) == pytest.approx(expected)
    # Model has been run once, so there is one entry in the cache
    assert len(mdl._cache) == 1
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)

    mdl.cx1 = 5
    mdl.cy1 = 2

    expected2 = 1 + 2 * x0grid + 5 * x1grid
    assert mdl(x1grid, x0grid) == pytest.approx(expected2)
    assert len(mdl._cache) == 2

    # Now, the model has run a second time, so there are two entries in the cache
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)
    assert mdl._cache[list(mdl._cache.keys())[1]] == pytest.approx(expected2)


def test_cache_is_actually_used():
    """Most other tests check that the cache has values in it,
    but not that those values are actually returned.
    Here, we manipulate the cached value and then call the model
    to check that the cached value is used.
    """
    x0grid = [2, 10, 1.5]
    x1grid = [-2, 0, 3.]

    mdl = Polynom2D()
    mdl.integrate = False
    mdl.cache = 5

    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.ones(3)
    assert mdl(x0grid, x1grid) == pytest.approx(expected)

    # Manipulate the values in the cache
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * expected
    assert mdl(x0grid, x1grid) == pytest.approx(2 * expected)


def test_evaluate_cache_many_parameters():
    """Check that we cache correctly with > 2 parameters, e.g.
    for a Data2DInt (with x0lo, x1lo, x0hi, x1hi) or an ND model.
    Since Data2DInt is used more often, there are models defined for
    that in Sherpa already and that's an easier approach to write this test.
    """
    x0lo = np.arange(-20, -10, 1)
    x0hi = x0lo + 1
    x1lo = np.arange(-3.2, 8.4, 1.2)
    x1hi = x1lo + 1

    mdl = Polynom2D()
    mdl.integrate = True
    mdl.cache = 5

    assert len(mdl._cache) == 0
    out = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 1
    # Manipulate the values in the cache to check if it's used
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * out

    # second call should return value from cache
    out1 = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 1
    assert pytest.approx(out1) == 2 * out
    # Call with third array different, so should not return cached value
    x0hi2 = x1hi + 0.1
    out2 = mdl(x0lo, x1lo, x0hi2, x1hi)
    assert len(mdl._cache) == 2
    assert pytest.approx(out2) != out1
    assert pytest.approx(out2) != 2 * out
    # And just for good measure, try the first call again
    out3 = mdl(x0lo, x1lo, x0hi, x1hi)
    assert len(mdl._cache) == 2
    assert pytest.approx(out3) == 2 * out


def test_evaluate_cache_regrid2d():
    """How about a regridded model?"""

    mdl = Lorentz2D()

    x0 = np.arange(0, 10, 0.1)
    x1 = np.arange(2, 20, 0.1)
    rmdl = mdl.regrid(x0, x1)

    assert isinstance(rmdl, RegridWrappedModel)
    assert not hasattr(rmdl, 'cache')


def test_evaluate_cache2dint():
    """Check we run with caching on: 2dint"""

    x0lo = np.array([-20, -1, 10])
    x1lo = np.array([-25, -1, -10])
    x0hi = np.array([-10, 1, 15])
    x1hi = np.array([-10, 1, -8])

    mdl = Polynom2D()
    mdl.cache = 12
    assert len(mdl._cache) == 0

    # Check the default values
    expected = np.array([150, 4, 10])
    assert mdl(x0lo, x1lo, x0hi, x1hi) == pytest.approx(expected)
    assert len(mdl._cache) == 1
    assert mdl._cache[list(mdl._cache.keys())[0]] == pytest.approx(expected)

    # Manipulate the values in the cache to check if it's used
    mdl._cache[list(mdl._cache.keys())[0]] = 2 * expected
    assert mdl(x0lo, x1lo, x0hi, x1hi) == pytest.approx(2 * expected)

