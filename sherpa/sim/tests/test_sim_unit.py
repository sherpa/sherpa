#
#  Copyright (C) 2017, 2021  Smithsonian Astrophysical Observatory
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
See sherpa/astro/sim/tests_astro_sim_unit.py for the astro-specific
version of this.
"""

import pytest

from sherpa import sim
from sherpa.stats import LeastSq


# This is part of #397
#
def test_list_samplers():
    """Ensure list_samplers returns a list."""

    mcmc = sim.MCMC()
    samplers = mcmc.list_samplers()

    assert isinstance(samplers, list)
    assert len(samplers) > 0


def test_list_samplers_contents():
    """Are the expected values included"""

    # Test that the expected values exist in this list,
    # but do not enforce these are the only values.
    #
    samplers = sim.MCMC().list_samplers()
    for expected in ['mh', 'metropolismh']:
        assert expected in samplers


def test_results_repr():
    """Basic check of string output of LikelihoodRatioResults"""

    # We don't care that the results are meaningless
    ratios = [0.2, 0.4, 0.3]
    stats = [[2.1, 3], [2.2, 2.8], [4.2, 1.3]]
    samples = [1, 2, 3]  # what should this be?
    x = sim.LikelihoodRatioResults(ratios, stats, samples, 0.2, 0.4, 12.3, 14.4)

    out = repr(x)
    assert out == '<Likelihood ratio results instance>'


def test_results_str():
    """Basic check of string output of LikelihoodRatioResults"""

    # We don't care that the results are meaningless, which does make
    # the string output a bit annoying to validate.
    #
    ratios = [0.2, 0.4, 0.3]
    stats = [[2.1, 3], [2.2, 2.8], [4.2, 1.3]]
    samples = [1, 2, 3]  # what should this be?
    x = sim.LikelihoodRatioResults(ratios, stats, samples, 0.2, 0.4, 12.3, 14.4)

    out = str(x)
    toks = out.split('\n')
    assert len(toks) == 9
    assert toks[0].startswith('samples = [1')
    assert toks[1].startswith('stats   = [[2')
    assert toks[4].startswith('ratios  = [0')
    assert toks[5] == 'null    = 12.3'
    assert toks[6] == 'alt     = 14.4'
    assert toks[7] == 'lr      = 0.2'
    assert toks[8] == 'ppp     = 0.4'


def test_lrt_dows_not_like_gaussian():
    """A very basic check that LRT errors out for gaussian.

    We are taking advantage of the code to not have to set up very
    much infrastructure (e.g. calling run with non-sensical
    values like None).
    """

    lrt = sim.LikelihoodRatioTest()
    with pytest.raises(TypeError) as exc:
        lrt.run(None, None, None, stat=LeastSq())

    emsg = 'Sherpa fit statistic must be Cash or CStat for likelihood ratio test'
    assert str(exc.value) == emsg
