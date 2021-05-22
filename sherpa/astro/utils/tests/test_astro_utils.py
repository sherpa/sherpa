#
#  Copyright (C) 2008, 2016, 2018, 2020, 2021
#       Smithsonian Astrophysical Observatory
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

from sherpa.astro.utils import is_in, compile_energy_grid

LONGVALS = [100, 249, 400, 450, 500, 601, 1024]
SHORTVALS = [100, 249, 601, 1024]


# test_is_in_long and test_not_in_short were originally a
# single test case called test_response_filter_logic
#
@pytest.mark.parametrize("lo, hi",
                         [(50, 2400),  # outside case
                          (100, 200),  # lo case
                          (50, 1024),  # hi case
                          (250, 2000),  # 'hidden' lo case
                          (50, 250),  # 'hidden' hi case
                          (250, 600)])  # 'hidden' interval case w/ noticed channels inside
def test_is_in_long(lo, hi):
    assert is_in(LONGVALS, lo, hi)


def test_not_in_short():
    # 'hidden' interval case w/ *no* noticed channels inside
    assert not is_in(SHORTVALS, 250, 600)


def test_compile_energy_grid():
    '''This test just executes the doc string. Once we activate doctests,
    this test can be removed.'''
    elo, ehi, htable = compile_energy_grid([([1, 3, 5], [3, 5, 7]),
                                            ([0, 1, 2], [1, 2, 3])])
    assert np.all(elo == [0, 1, 2, 3, 5])
    assert np.all(ehi == [1, 2, 3, 5, 7])
    assert np.all(htable[0][0] == [1, 3, 4])


def test_compile_energy_grid_3inputs():
    '''Test more than two input grids'''
    elo, ehi, htable = compile_energy_grid([([1, 3, 5], [3, 5, 7]),
                                            ([0, 1, 2], [1, 2, 3]),
                                            ([0.5, 5.5, 7.5], [5.5, 7.5, 9])])
    assert np.allclose(elo, [0., 0.5, 1., 2., 3., 5., 5.5, 7., 7.5])
    assert np.allclose(ehi, [0.5, 1., 2., 3., 5., 5.5, 7., 7.5, 9.])
    assert np.all(htable[0][0] == [2, 4, 5])
    assert np.all(htable[2][1] == [5, 7, 8])
