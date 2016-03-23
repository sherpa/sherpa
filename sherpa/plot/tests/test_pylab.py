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

from sherpa.utils.test import SherpaTestCase, requires_pylab


@requires_pylab
class pylab_test(SherpaTestCase):

    def test_axes_default(self):
        from sherpa.plot.pylab_backend import _errorbar_defaults
        # assert all needed defaults have been found
        assert all(e in _errorbar_defaults for e in ('ecolor', 'capsize', 'barsabove'))
