#
#  Copyright (C) 2025
#  MIT
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
"""Speed tests for the XSPEC additive models.

These tests check if additive XSPEC models are evaluated faster with
the norm decorator than without it. For complex models that is likely
the case, but for very simple models (like a constant) looking up a value
in the cache is comparable or even faster than skipping the cache-lookup
and evaluating the model every time.
Of course, NOT caching also reduced the memory footprint, so it should
only be done if it is beneficial on runtime.

The tests in this file are not run by default, because there are too many
false positive and false negatives from just other random processes in the
operating system that can affect the runtime of the tests.
(Of course, one could be more sophisticated and run the tests multiple times,
but the runtime is already > 2 min on my machine, so I don't want to do that.)

To run tests in this file use `pytest --runspeed`.

In tests on a Mac M1 system on Sherpa 4.17 we find the following models
are consistently faster without the decorator:

- bbody, bbodyrad, zbbody
- diskpn
- expdec
- nsa, nsagrav, nsatmos, nsmax
- optxagn, optxagnf
- pregpwrlw
- redge, refsch

For larger grids also: xposm
but in ns-type models, it's only nsa, nsagrav, not the other two.
"""
import copy
import inspect
import time
import types

import numpy as np
import pytest

from sherpa.utils.testing import requires_xspec
try:
    from sherpa.astro.xspec.tests.test_xspec import get_xspec_models, _hc
except ImportError:
    # The tests won't run anyway without XSPEC, so just define a dummy
    get_xspec_models = lambda: []


def make_grid_dense():
    """Return the 'standard' contiguous grid used in these tests.

    Returns elo, ehi, wlo, whi where the e-xxx values are in keV and
    the w-xxx values in Angstrom. The *lo/*hi values are this grid
    separated out into bin edges.

    The following condition holds: *hi[i] > *lo[i], and the
    values in the two arrays are either monotonically increasing
    or decreasing.

    """
    # XSgaussian has peak at 6.5 keV and sigma ~ 0.1 keV
    # so need to get to around 7 keV
    egrid = np.linspace(0.1, 10.01, num=1024, endpoint=True)
    elo = egrid[:-1]
    ehi = egrid[1:]

    wgrid = _hc / egrid
    whi = wgrid[:-1]
    wlo = wgrid[1:]

    return elo, ehi, wlo, whi


@pytest.mark.speed
@requires_xspec
@pytest.mark.parametrize("modelcls", get_xspec_models())
def test_evaluate_additive_xspec_model_normwrapper(modelcls):
    """Does the `sherpa.astro.xspec.eval_xspec_with_fixed_norm`
    wrapper change the results?

    It should not (within numerical tolerances).
    The hard part in this test is to obtain a "clean" version of the
    model without the wrapper.
    """

    from sherpa.astro import xspec

    # use an identifier in case there is an error
    mdl = modelcls()
    # Skip everything that is not an additive model
    # Easier to filter out here given the current design.
    if not isinstance(mdl, xspec.XSAdditiveModel):
        return

    elo, ehi, wlo, whi = make_grid_dense()

    # run once in case there is a startup time for loading
    mdl.norm = 1
    mdl.cache = 5
    evals = mdl(elo, ehi)

    # The norm=1 case has been run for initialization.
    # So, this should pull the value from the cache and be very fast.
    start_time_decorated = time.perf_counter()
    evals = mdl(elo, ehi)
    end_time_decorated = time.perf_counter()
    assert mdl._cache_ctr['hits'] == 1
    # Now, modify the model to remove `sherpa.astro.xspec.eval_xspec_with_fixed_norm`
    # decorator and empty the cache to make sure we are not using the cached value.
    mdl.calc = types.MethodType(inspect.unwrap(mdl.calc), mdl)
    mdl.cache = 0
    start_time_clean = time.perf_counter()
    evals_no_wrapper = mdl(elo, ehi)
    end_time_clean = time.perf_counter()
    assert mdl._cache_ctr['hits'] == 0

    assert evals == pytest.approx(evals_no_wrapper)

    run_time_decorated = end_time_decorated - start_time_decorated
    run_time_clean = end_time_clean - start_time_clean

    if modelcls().cache == 0:
        assert run_time_decorated > run_time_clean
    else:
        assert run_time_decorated < run_time_clean
