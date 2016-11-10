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

TEST = [1, 2, 3]
TEST2 = [4, 5, 6]


# bug #262
def test_list_ids():
    session = Session()
    session.load_arrays(1, TEST, TEST)
    session.load_arrays("1", TEST, TEST2)

    # order of 1 and "1" is not determined
    assert {1, "1"} == set(session.list_data_ids())
    assert_array_equal(TEST2, session.get_data('1').get_dep())
    assert_array_equal(TEST, session.get_data(1).get_dep())


# bug #297
def test_save_restore(tmpdir):
    outfile = tmpdir.join("sherpa.save")
    session = Session()
    session.load_arrays(1, TEST, TEST2)
    session.save(str(outfile), clobber=True)
    session.clean()
    assert set() == set(session.list_data_ids())

    session.restore(str(outfile))
    assert {1, } == set(session.list_data_ids())
    assert_array_equal(TEST, session.get_data(1).get_indep()[0])
    assert_array_equal(TEST2, session.get_data(1).get_dep())
