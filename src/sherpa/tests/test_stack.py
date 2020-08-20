#
#  Copyright (C) 2007, 2016, 2018, 2020  Smithsonian Astrophysical Observatory
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


from sherpa.testing import SherpaTestCase, requires_stk

import os
_this_dir = os.path.dirname(__file__)


@requires_stk
class test_stack(SherpaTestCase):

    def test_build_stack(self):
        import stk

        def get_name(name):
            return '/'.join((_this_dir, name))

        out = stk.build('@+{}/{}'.format(_this_dir, 'a.lis'))
        self.assertEqual([get_name('a'), get_name('a1'), get_name('a2'),
                          get_name('b'), get_name('b1'), get_name('b2')], out)
