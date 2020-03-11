#
#  Copyright (C) 2008, 2016, 2018, 2020
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

import pytest

from sherpa.astro.utils import is_in

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
                          (250, 600)  # 'hidden' interval case w/ noticed channels inside
                         ])
def test_is_in_long(lo, hi):
    assert is_in(LONGVALS, lo, hi)


def test_not_in_short():
    # 'hidden' interval case w/ *no* noticed channels inside
    assert not is_in(SHORTVALS, 250, 600)
