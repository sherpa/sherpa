#
#  Copyright (C) 2023, 2024
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

A number of routines have been added to support writing code that
allows both the legacy (pre 1.17) NumPy random API to be used
(with rng set to None) or the new generator approach.

"""

from typing import Literal, Optional, Sequence, SupportsFloat, Union, \
    overload

import numpy as np

from .numeric_types import SherpaFloat

# This should probably be an explicit type alias but for now leave it
# like this. Users are expected to use the Generator rather than
# RandomState form, but both are kept in (e.g. for testing).
#
# Once Python 3.10 or 3.12 is the minimum we can make this an actual
# alias that we can export. For now this "alias" can be used by other
# modules but the documentation is not ideal and the name may change.
#
RandomType = Union[np.random.Generator, np.random.RandomState]


__all__ = ("chisquare", "choice", "integers",
           "multivariate_normal", "normal", "poisson_noise",
           "random", "standard_normal", "uniform")


@overload
def poisson_noise(x: SupportsFloat,
                  rng: Optional[RandomType] = None
                  ) -> SherpaFloat:
    ...

@overload
def poisson_noise(x: Sequence[SupportsFloat],
                  rng: Optional[RandomType] = None
                  ) -> np.ndarray:
    ...

def poisson_noise(x, rng=None):
    """Draw samples from a Poisson distribution.


    Parameters
    ----------
    x : scalar or array
       The expectation value for the distribution.
    rng : numpy.random.Generator, numpy.random.RandomState, or None, optional
       Determines how the random numbers are created. If set to
       None then the `numpy.random.poisson` routine is used.

    Returns
    -------
    out : scalar or array
       A random realisation of the input array, drawn from
       the Poisson distribution, as a `SherpaFloat`.

    Notes
    -----
    All input values less than zero are replaced by zero.

    Examples
    --------

    When the rng parameter is left as None then the legacy Numpy
    random API is used:

    >>> np.random.seed(17389)
    >>> poisson_noise([10, 20, 5])
    array([ 7.,  20.,   6.])

    Note that the seed used by the legacy Numpy is not guaranteed to
    match the behavior of the numpy generators:

    >>> rng = np.random.default_rng(17389)
    >>> poisson_noise([10, 20, 5], rng=rng)
    array([12., 31., 7.])

    """

    if rng is None:
        simulate = np.random.poisson
    else:
        simulate = rng.poisson

    x = np.asarray(x)
    if x.shape == ():
        x = SherpaFloat(x)
        if x <= 0.0:
            x = 0.0
        else:
            x = simulate(x)

        return SherpaFloat(x)

    x_out = np.zeros(x.shape, SherpaFloat)
    good = np.where(x > 0.0)
    x_out[good] = simulate(x[good])
    return x_out


def random(rng: Optional[RandomType]) -> float:
    """Create a random value [0, 1.0)

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.random_sample` routine is used.

    Returns
    -------
    value : number

    """

    if rng is None:
        return np.random.random_sample()

    return rng.random()


# Try to come up with a sensible type for the size field.
SizeType = Union[int, Sequence[int], np.ndarray]


@overload
def uniform(rng: Optional[RandomType],
            low: float,
            high: float,
            size: Literal[None]
            ) -> float:
    ...

@overload
def uniform(rng: Optional[RandomType],
            low: float,
            high: float,
            size: SizeType
            ) -> np.ndarray:
    ...

def uniform(rng, low, high, size=None):
    """Create a random value within a uniform range.

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.uniform` routine is used.
    low, high
       The range [low, high].
    size
       The shape and size of the return.

    Returns
    -------
    value : number or ndarray

    """

    if rng is None:
        return np.random.uniform(low, high, size=size)

    return rng.uniform(low, high, size=size)


def integers(rng: Optional[RandomType],
             high: int
             ) -> int:
    """Create a random integer from [0, high).

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.randint` routine is used.
    high
       The upper limit (not inclusive).

    Returns
    -------
    value : number

    """

    if rng is None:
        return np.random.randint(high)

    # It appears that this is one place where we can not use the same
    # method for Generator and RandomState.
    #
    try:
        return rng.integers(high)  # type: ignore[union-attr]
    except AttributeError:
        # Fall back on RandomState
        return rng.randint(high)  # type: ignore[union-attr]


@overload
def normal(rng: Optional[RandomType],
           loc: float,
           scale: float,
           size: Literal[None]
           ) -> float:
    ...

@overload
def normal(rng: Optional[RandomType],
           loc: float,
           scale: float,
           size: SizeType
           ) -> np.ndarray:
    ...

def normal(rng, loc=0.0, scale=1.0, size=None):
    """Create a random value from a normal distribution.

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.normal` routine is used.
    loc : number, optional
       The location of the normal.
    scale : number, optional
       The scale of the normal.
    size : optional
       The shape and number of elements to create.

    Returns
    -------
    value : number or ndarray

    """

    if rng is None:
        return np.random.normal(loc=loc, scale=scale, size=size)

    return rng.normal(loc=loc, scale=scale, size=size)


@overload
def standard_normal(rng: Optional[RandomType],
                    size: Optional[None]
                    ) -> float:
    ...

@overload
def standard_normal(rng: Optional[RandomType],
                    size: SizeType
                    ) -> np.ndarray:
    ...

def standard_normal(rng, size=None):
    """Create a random value from a normal distribution (mean=0, stdev=1).

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.standard_normal` routine is used.
    size : optional
       The shape and number of elements to create.

    Returns
    -------
    value : number or ndarray

    """

    if rng is None:
        return np.random.standard_normal(size=size)

    return rng.standard_normal(size=size)


def multivariate_normal(rng: Optional[RandomType],
                        mean,
                        cov,
                        size=None
                        ) -> np.ndarray:
    """Create a random value from a multivariate normal distribution.

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.multivariate_normal` routine is used.
    mean, cov
       The mean and covariance matrix.
    size : optional
       The shape and number of elements to create.

    Returns
    -------
    value : number or ndarray

    """

    if rng is None:
        return np.random.multivariate_normal(mean, cov, size=size)

    return rng.multivariate_normal(mean, cov, size=size)


def chisquare(rng: Optional[RandomType],
              df,
              size=None
              ) -> np.ndarray:
    """Create a random value from a multivariate normal distribution.

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.chisquare` routine is used.
    df
       The degrees of freedom.
    size : optional
       The shape and number of elements to create.

    Returns
    -------
    value : number

    """

    if rng is None:
        return np.random.chisquare(df, size=size)

    return rng.chisquare(df, size=size)


def choice(rng: Optional[RandomType],
           xs: Sequence,
           n
           ) -> np.ndarray:
    """Create a subset of elements from xs with no duplication.

    Parameters
    ----------
    rng : numpy.random.Generator, numpy.random.RandomState, or None
       Determines how the random numbers are created. If set to
       None then the `numpy.random.choice` routine is used.
    xs
       Sequence of values.
    n
       The number of values to select from xs.

    Returns
    -------
    values : ndarray

    """

    if rng is None:
        return np.random.choice(xs, n, replace=False)

    return rng.choice(xs, n, replace=False)
