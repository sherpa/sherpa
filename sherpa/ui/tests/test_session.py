#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

from sherpa.ui.utils import Session
from numpy.testing import assert_array_equal


# bug #262
def test_list_ids():
    session = Session()
    session.load_arrays(1, [1, 2, 3], [1, 2, 3])
    session.load_arrays("1", [1, 2, 3], [4, 5, 6])
    assert [1, "1"] == session.list_data_ids()
    assert_array_equal([4, 5, 6], session.get_data('1').get_dep())
    assert_array_equal([1, 2, 3], session.get_data(1).get_dep())
