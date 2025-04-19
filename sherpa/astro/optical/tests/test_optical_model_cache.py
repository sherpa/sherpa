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
"""Cache tests for the models.

These tests check if models are evaluated faster without caching than with
caching. For very simple models (like a constant) looking up a value
in the cache is comparable or even faster than skipping the cache-lookup
and evaluating the model every time.
Of course, NOT caching also reduced the memory footprint, so it should
only be done if it is beneficial on runtime.

The tests in this file are not run by default, because there are too many
false positive and false negatives from just other random processes in the
operating system that can affect the runtime of the tests.

To run tests in this file use `pytest --runspeed`.

"""
import pytest
from sherpa.astro import optical
from sherpa.models import ArithmeticModel
from sherpa.conftest import get_runtime_with_cache


modellist = optical.__all__
modellist = list(set(modellist) - set(['AbsorptionVoigt', 'EmissionVoigt']))
modellist.sort()


@pytest.mark.speed
@pytest.mark.parametrize("modelcls", modellist)
def test_runtime_with_cache(modelcls):
    """Check if a model is evaluated faster with cache than without.

    This test will fail if the runtime with cache is slower than without;
    independent of the default setting of the cache. It is meant for developers
    to interactively try out things, and should not be run by default.
    """

    # use an identifier in case there is an error
    mdl = getattr(optical, modelcls)()
    # Skip everything that is not an Arithmetic model
    # Easier to filter out here given the current design.
    if not isinstance(mdl, ArithmeticModel):
        return

    # Cache for 2D models still under discussion
    if mdl.ndim > 1:
        return

    get_runtime_with_cache(mdl, timing=True)