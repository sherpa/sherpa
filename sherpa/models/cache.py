#
#  Copyright (C) 2010, 2016 - 2024
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
from __future__ import annotations
import functools

import numpy as np
from typing import Callable
from sherpa.astro import hc

# What routine do we use for the hash in modelCacher1d?  As we do not
# need cryptographic security go for a "quick" algorithm, but md5 is
# not guaranteed to always be present.  There has been no attempt to
# check the run times of these routines for the expected data sizes
# they will be used with.
#
try:
    from hashlib import md5 as hashfunc
except ImportError:
    from hashlib import sha256 as hashfunc


__all__ = ['modelCacher1d', 'modelCacher1d_exp']


def boolean_to_byte(boolean_value: bool) -> bytes:
    """Convert a boolean to a byte value.

    Parameters
    ----------
    boolean_value : bool
        The value to convert. If not a boolean then it is
        treated as `False`.

    Returns
    -------
    val : bytes
        b'1' if `True` otherwise b'0'.
    """

    bmap = {True: b'1', False: b'0'}
    return bmap.get(boolean_value, b'0')


def modelCacher1d(func: Callable) -> Callable:
    """A decorator to cache 1D ArithmeticModel evaluations.

    Apply to the `calc` method of a 1D model to allow the model
    evaluation to be cached. The decision is based on the
    `_use_caching` attribute of the cache along with the `integrate`
    setting, the evaluation grid, parameter values, and the keywords
    sent to the model.

    Notes
    -----
    The keywords are included in the hash calculation even if they are
    not relevant for the model (as there's no easy way to find this
    out).

    Example
    -------

    Allow `MyModel` model evaluations to be cached::

        def MyModel(ArithmeticModel):
            ...
            @modelCacher1d
            def calc(self, p, *args, **kwargs):
                ...

    """

    @functools.wraps(func)
    def cache_model(cls, pars, xlo, *args, **kwargs):
        # Counts all accesses, even those that do not use the cache.
        cache_ctr = cls._cache_ctr
        cache_ctr['check'] += 1

        # Short-cut if the cache is not being used.
        #
        if not cls._use_caching:
            return func(cls, pars, xlo, *args, **kwargs)

        # Up until Sherpa 4.12.2 we used the kwargs to define the
        # integrate setting, with
        # boolean_to_byte(kwargs.get('integrate', False)) but
        # unfortunately this is used in code like
        #
        #    @modelCacher1d
        #    def calc(..):
        #        kwargs['integrate'] = self.integrate
        #        return somefunc(... **kwargs)
        #
        # and the decorator is applied to calc, which is not
        # called with a integrate kwarg, rather than the call to
        # somefunc, which was sent an integrate setting.
        #
        try:
            integrate = cls.integrate
        except AttributeError:
            # Rely on the integrate kwarg as there's no
            # model setting.
            #
            integrate = kwargs.get('integrate', False)

        data = [np.array(pars).tobytes(),
                boolean_to_byte(integrate),
                np.asarray(xlo).tobytes()]
        if args:
            data.append(np.asarray(args[0]).tobytes())

        # Add any keyword arguments to the list. This will
        # include the xhi named argument if given. Can the
        # value field fail here?
        #
        for k, v in kwargs.items():
            data.extend([k.encode(), np.asarray(v).tobytes()])

        # Is the value cached?
        #
        token = b''.join(data)
        digest = hashfunc(token).digest()
        cache = cls._cache
        if digest in cache:
            cache_ctr['hits'] += 1
            return cache[digest].copy()

        # Evaluate the model.
        #
        vals = func(cls, pars, xlo, *args, **kwargs)

        # remove first item in queue and remove from cache
        queue = cls._queue
        key = queue.pop(0)
        cache.pop(key, None)

        # append newest model values to queue
        queue.append(digest)
        cache[digest] = vals.copy()

        cache_ctr['misses'] += 1

        return vals

    return cache_model


# The inner function shares a lot of the code with modelCacher1d, but
# it's annoyingly hard to factor out the common code, so it's not
# done here.
def modelCacher1d_exp(canonical_nh : list[float]=[0.1],
                              nh_bound : list[float]=[np.inf]) -> Callable:
    def modelCacher1d_exp_inner(func: Callable) -> Callable:
        r"""A decorator to cache XSPEC absorption models.

        Apply to the `calc` method of a 1D model to allow the model
        evaluation to be cached. The decision is based on the
        `_use_caching` attribute of the cache along with the `integrate`
        setting, the evaluation grid, parameter values, and the keywords
        sent to the model.

        Unlike the more general `sherpa.models.model.modelCacher1d`, this
        version is designed for models with a single parameter and a model of
        the form :math:`y(x) = \exp(p * a(x))`.

        where `p` is the parameter. The cacher will evaluate the model once
        with `p=0.1` and then store that. This is useful in particular for
        absorption models, where `a(x)` is an expensive calculation with
        atomic cross-sections etc.

        Use this only with XSPEC models
        -------------------------------
        This cacher has some specific assumption built in how XSPEC
        decides if it receives a grid in wavelength or in energy space.
        Those assumptions will fail for normal Sherpa models.

        Notes
        -----
        The keywords are included in the hash calculation even if they are
        not relevant for the model (as there's no easy way to find this
        out).

        Example
        -------

        Allow `MyModel` model evaluations to be cached::

            def MyModel(ArithmeticModel):
                ...
                @modelCacher1d_exp()
                def calc(self, p, *args, **kwargs):
                    ...

        """
        @functools.wraps(func)
        def cache_model(cls, pars, xlo, *args, **kwargs):

            # Assert, not exception since this is a developer error.
            assert len(pars) == 1, "Only one parameter is allowed"

            # Counts all accesses, even those that do not use the cache.
            cache_ctr = cls._cache_ctr
            cache_ctr["check"] += 1

            # Short-cut if the cache is not being used.
            #
            if not cls._use_caching:
                return func(cls, pars, xlo, *args, **kwargs)

            try:
                integrate = cls.integrate
            except AttributeError:
                # Rely on the integrate kwarg as there's no
                # model setting.
                #
                integrate = kwargs.get("integrate", False)

            data = [
                boolean_to_byte(integrate),
                np.asarray(xlo).tobytes(),
            ]
            if args:
                data.append(np.asarray(args[0]).tobytes())

            # Add any keyword arguments to the list. This will
            # include the xhi named argument if given. Can the
            # value field fail here?
            #
            for k, v in kwargs.items():
                data.extend([k.encode(), np.asarray(v).tobytes()])

            # Is the value cached?
            #
            token = b"".join(data)
            digest = hashfunc(token).digest()
            cache = cls._cache

            # XSPEC assumes that the grid is in keV if xlo is in increasing order
            # and in Ang if xlo is in decreasing order.
            # For us, it's simpler to assume that the grid is in keV.
            # We are not checking that the grid is monotonically increasing,
            # because we're assuming that's done elsewhere.

            # It's not clear if xlo can be anything but an array in real live
            # but we have tests that pass in plain Python lists.
            if xlo[1] > xlo[0]:
                en_lo = np.asanyarray(xlo)
                if args:
                    en_hi = np.asanyarray(args[0])
                else:
                    en_hi = ()
            else:
                if args:
                    en_lo = hc / np.asanyarray(args[0])
                    en_hi = hc / np.asanyarray(xlo)
                else:
                    en_lo = hc / np.asanyarray(xlo)
                    en_hi = ()

            nhbin = np.digitize(en_lo, nh_bound)
            nh_array = np.array(canonical_nh)[nhbin]

            vals = np.zeros_like(xlo, dtype=float)

            if digest in cache:
                cache_ctr["hits"] += 1

            else:
                # Evaluate the model.
                for i, this_nh in enumerate(canonical_nh):
                    ind = nhbin == i
                    if ind.sum() == 0:
                        continue
                    # XSPEC models needs to be called with at least two energy bins.
                    # If the split is such that we have only one bin, we need to
                    # work around that.
                    if ind.sum() == 1:
                        # The easiest is if we have xlo, xhi as input
                        if args:
                            this_xlo = np.append(en_lo[ind], en_hi[ind])
                            this_args = (np.append(en_hi[ind], en_hi[ind] + 0.01), )
                        else:
                            this_xlo = np.append(xlo[ind], xlo[ind] + 0.01)
                            this_args = ( )
                        # Either way, we'll only use the value from the first bin.
                        use_xpec_out = [True, False]
                    if ind.sum() > 1:
                        this_xlo = en_lo[ind]
                        if args:
                            this_args = (en_hi[ind], )
                        else:
                            this_args = ( )
                        use_xpec_out = np.ones_like(this_xlo, dtype=bool)
                    xspecout = func(cls, (this_nh,), this_xlo, *this_args, **kwargs)
                    vals[ind] = xspecout[use_xpec_out]

                # remove first item in queue and remove from cache
                queue = cls._queue
                key = queue.pop(0)
                cache.pop(key, None)

                # append newest model values to queue
                queue.append(digest)
                cache[digest] = np.log(vals)

                cache_ctr["misses"] += 1

            return np.exp(pars[0] / nh_array * cache[digest])

        return cache_model
    return modelCacher1d_exp_inner