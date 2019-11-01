#
#  Copyright (C) 2016, 2018, 2019  Smithsonian Astrophysical Observatory
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

from sherpa.utils.testing import requires_pylab


@requires_pylab
def test_axes_default():
    """Have we default values for all the axis settings?"""

    from sherpa.plot.pylab_backend import _errorbar_defaults

    # assert all needed defaults have been found, but do not
    # check the actual values. Split into multiple lines so that
    # if there's a failure we can see which is False.
    #
    fields = ('ecolor', 'capsize', 'barsabove')
    defs = [e in _errorbar_defaults for e in fields]
    assert all(defs)
