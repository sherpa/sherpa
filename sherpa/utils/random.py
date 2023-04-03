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

A number of routines have been added to support writing code that
allows both the legacy (pre 1.17) NumPy random API to be used
(with rng set to None) or the new generator approach.

"""

import numpy as np

# Copy of sherpa.utils.SherpaFloat
#
SherpaFloat = np.float_


def poisson_noise(x, rng=None):
    """Draw samples from a Poisson distribution.


    Parameters
    ----------
    x : scalar or array
       The expectation value for the distribution.
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.

    Returns
    -------
    out : scalar or array
       A random realisation of the input array, drawn from
       the Poisson distribution, as a `SherpaFloat`.

    Notes
    -----
    The distribution is calculated by the NumPy poisson generator
    and all values zero or less are replaced by zero.

    Examples
    --------

    >>> poisson_noise([10, 20, 5])
    array([ 13.,  21.,   6.])

    >>> rng = np.random.default_rng(17389)
    >>> poisson_noise([1, 2, 2, 3], rng=rng)
    array([2., 2., 5., 3.])

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
            x = x = simulate(x)

        return SherpaFloat(x)

    x_out = np.zeros(x.shape, SherpaFloat)
    good = np.where(x > 0.0)
    x_out[good] = simulate(x[good])
    return x_out


def random(rng):
    """Create a random value [0, 1.0)

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.

    Returns
    -------
    value : number

    """

    if rng is None:
        return np.random.random_sample()

    return rng.random()


def uniform(rng, low, high, size=None):
    """Create a random value within a uniform range.

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
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


def integers(rng, high):
    """Create a random integer from [0, high).

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
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
        return rng.integers(high)
    except AttributeError:
        return rng.randint(high)


def normal(rng, loc=0, scale=1, size=None):
    """Create a random value from a normal distribution.

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
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


def standard_normal(rng, size=None):
    """Create a random value from a normal distribution (mean=0, stdev=1).

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
    size : optional
        The shape and number of elements to create.

    Returns
    -------
    value : number or ndarray

    """

    if rng is None:
        return np.random.standard_normal(size=size)

    return rng.standard_normal(size=size)


def multivariate_normal(rng, mean, cov, size=None):
    """Create a random value from a multivariate normal distribution.

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
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


def chisquare(rng, df, size=None):
    """Create a random value from a multivariate normal distribution.

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
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


def choice(rng, xs, n):
    """Create a subset of elements from xs with no duplication.

    Parameters
    ----------
    rng : np.random.Generator, np.random.RandomState, or None, optional
        If set, the generator is used to create the random numbers. If
        not set then the legacy numpy RandomState instance is used.
    xs
        Sequence of values.
    n
        The number of values to select from xs

    Returns
    -------
    values : ndarray

    """

    if rng is None:
        return np.random.choice(xs, n, replace=False)

    return rng.choice(xs, n, replace=False)
