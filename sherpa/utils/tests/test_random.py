#
#  Copyright (C) 2010, 2016, 2018 - 2024
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

from sherpa.utils.random import poisson_noise, \
    standard_normal, chisquare, choice
from sherpa.utils.numeric_types import SherpaFloat


def test_poisson_noise_checks_dtype():

    with pytest.raises(ValueError,
                       # pre NumPy 2.0 this ends with
                       #    'ham'
                       # and for >= 2.0 it is
                       #    np.str_('ham')
                       #
                       match="could not convert string to float: .*'ham'.*"):
        poisson_noise('ham')


def test_poisson_noise_checks_number_args():

    with pytest.raises(TypeError,
                       match=r"^poisson_noise\(\) takes from 1 to 2 positional arguments but 3 were given$"):
        poisson_noise(1, 2, 'ham')


def test_poisson_noise_checks_rng():

    with pytest.raises(AttributeError,
                       match="'int' object has no attribute 'poisson'"):
        poisson_noise(1, rng=123)


def test_poisson_noise_scalar_none():
    """Send in a scalar, get back a scalar"""

    out = poisson_noise(1000)
    assert type(out) == SherpaFloat
    assert out > 0.0


@pytest.mark.parametrize("x,expected", [(1000, 956), (4.3, 2)])
def test_poisson_noise_scalar_rng(x, expected):
    """Send in a scalar, get back a scalar"""

    rng = np.random.RandomState(43)
    out = poisson_noise(x, rng=rng)
    assert type(out) == SherpaFloat
    assert out == pytest.approx(expected)


@pytest.mark.parametrize("x", [-1000, 0])
@pytest.mark.parametrize("rng", [None, np.random.RandomState(23)])
def test_poisson_noise_scalar_not_positive(x, rng):

    out = poisson_noise(x, rng=rng)
    assert type(out) == SherpaFloat
    assert out == 0.0


def test_poisson_noise_array_none():
    """Send in an array, get back an array"""

    out = poisson_noise([1001, 1002, 0.0, 1003, -1004])
    assert type(out) == np.ndarray
    assert out.dtype.type == SherpaFloat

    # As the RNG is not set here, technically we could get zeros for
    # the positive expected values, but this is very unlikely.
    #
    ans = np.flatnonzero(out > 0.0)
    assert (ans == np.array([0, 1, 3])).all()


def test_poisson_noise_array_rng():
    """Send in an array, get back an array"""

    rng = np.random.RandomState(43)
    out = poisson_noise([12, 4.5, 0.0, 8.2, -1], rng=rng)
    assert type(out) == np.ndarray
    assert out.dtype.type == SherpaFloat

    assert out == pytest.approx([7, 2, 0, 8, 0])


def test_standard_normal_rng_is_none():
    """Coverage shows this was not tested, so check it

    It is not obvious if this is a real coverage failure,
    or just an issue with how the tests are run, but may
    as well add the check.
    """

    SEED = 1234
    SIZE = 5
    try:
        np.random.seed(SEED)
        rng = np.random.RandomState(SEED)
        a = standard_normal(rng=None, size=SIZE)
        b = standard_normal(rng=rng, size=SIZE)
    finally:
        np.random.seed()

    assert a == pytest.approx(b)


def test_chisquare_rng_is_none():
    """Coverage shows this was not tested, so check it

    It is not obvious if this is a real coverage failure,
    or just an issue with how the tests are run, but may
    as well add the check.
    """

    SEED = 1234
    SIZE = 5
    DF = 2
    try:
        np.random.seed(SEED)
        rng = np.random.RandomState(SEED)
        a = chisquare(rng=None, df=DF, size=SIZE)
        b = chisquare(rng=rng, df=DF, size=SIZE)
    finally:
        np.random.seed()

    assert a == pytest.approx(b)


def test_choice_rng_is_none():
    """Coverage shows this was not tested, so check it

    It is not obvious if this is a real coverage failure,
    or just an issue with how the tests are run, but may
    as well add the check.
    """

    SEED = 1234
    SIZE = 3
    VALS = [1, "xx", 45, " X "]
    try:
        np.random.seed(SEED)
        rng = np.random.RandomState(SEED)
        a = choice(rng=None, xs=VALS, n=SIZE)
        b = choice(rng=rng, xs=VALS, n=SIZE)
    finally:
        np.random.seed()

    assert a == pytest.approx(b)
