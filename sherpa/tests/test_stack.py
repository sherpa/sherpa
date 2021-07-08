#
#  Copyright (C) 2007, 2016, 2018, 2021  Smithsonian Astrophysical Observatory
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

import os

from sherpa.utils.testing import requires_stk


_this_dir = os.path.dirname(__file__)


def get_name(name):
    return '/'.join((_this_dir, name))


@requires_stk
def test_build_stack():
    import stk

    names = ['a', 'a1', 'a2', 'b', 'b1', 'b2']
    expected = [get_name(n) for n in names]

    out = stk.build('@+{}/{}'.format(_this_dir, 'a.lis'))
    for outval, expval in zip(out, expected):
        assert outval == expval

    assert len(out) == len(expected)
