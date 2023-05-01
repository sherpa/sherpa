#
#  Copyright (C) 2017, 2021, 2023, 2024
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
See sherpa/astro/sim/tests_astro_sim_unit.py for the astro-specific
version of this.
"""

import numpy as np

import pytest

from sherpa import sim
from sherpa.sim.mh import dmvnorm, dmvt, rmvt
from sherpa.stats import Chi2DataVar, LeastSq


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


def make_sim():
    """Create a LikelihoodRatioResults object for testing.

    We don't care that the results are meaningless as this is
    mainly to check the display output.
    """

    ratios = [0.2, 0.4, 0.3]
    stats = [[2.1, 3], [2.2, 2.8], [4.2, 1.3]]
    samples = [1, 2, 3]  # what should this be?
    pnames = ["a.a", "a.b", "b.a"]
    pvals = [[1, 2, 3], [2, 3, 4], [1.5, 2.5, 3.5]]
    return sim.LikelihoodRatioResults(ratios, stats, samples,
                                      0.2, 0.4, 12.3, 14.4,
                                      parnames=pnames,
                                      parvals=pvals)


def test_results_repr():
    """Basic check of string output of LikelihoodRatioResults"""

    x = make_sim()
    assert repr(x) == '<Likelihood ratio results instance>'


def test_results_str():
    """Basic check of string output of LikelihoodRatioResults"""

    x = make_sim()
    toks = str(x).split('\n')
    assert len(toks) == 9
    assert toks[0].startswith('samples = [1,')
    assert toks[1].startswith('stats   = [[2.1,3. ')
    assert toks[4] == 'ratios  = [0.2,0.4,0.3]'
    assert toks[5] == 'null    = 12.3'
    assert toks[6] == 'alt     = 14.4'
    assert toks[7] == 'lr      = 0.2'
    assert toks[8] == 'ppp     = 0.4'


@pytest.mark.parametrize("cls", [LeastSq, Chi2DataVar])
def test_lrt_does_not_like_gaussian(cls):
    """A very basic check that LRT errors out for gaussian.

    We are taking advantage of the code to not have to set up very
    much infrastructure (e.g. calling run with non-sensical
    values like None).
    """

    lrt = sim.LikelihoodRatioTest()
    with pytest.raises(TypeError) as exc:
        lrt.run(None, None, None, stat=cls())

    emsg = 'Sherpa fit statistic must be Cash, CStat, or WStat for likelihood ratio test'
    assert str(exc.value) == emsg


def test_lrt_checks_argument_size():
    """Check we error out"""

    with pytest.raises(ValueError,
                       match=r"^len\(parnames\) = 3  parvals.shape = \(2, 2\)$"):
        # These values don't make much sense as all we care about are
        # the last two.
        sim.LikelihoodRatioResults([2, 1.5],
                                   [[14, 12], [13, 12.2]],
                                   [[4], [3]],
                                   1.87, 0.05, 11.9, 12.1,
                                   ["a.x", "b.y", "b.z"],
                                   [[12, 13], [14, 15]])


def test_rmvt_checks_dof():
    """Error check"""

    with pytest.raises(ValueError,
                       match="^The degrees of freedom must be > 0$"):
        rmvt(14.3, [[0.1, 0.05], [0.07, 0.08]], 0)


def test_dmvt_checks_symmetry():
    """Error check"""

    sigma = np.asarray([[0.1, 0.05], [0.07, 0.08]])
    with pytest.raises(ValueError,
                       match="^Error: sigma is not symmetric$"):
        dmvt(14.3, 2.3, sigma, 1)


def test_dmvnorm_checks_symmetry():
    """Error check"""

    sigma = np.asarray([[0.1, 0.05], [0.07, 0.08]])
    with pytest.raises(ValueError,
                       match="^Error: sigma is not symmetric$"):
        dmvnorm(14.3, 2.3, sigma)
