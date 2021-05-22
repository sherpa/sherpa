#
#  Copyright (C) 2019  Smithsonian Astrophysical Observatory
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

# Although the Meta module does not need a FITS backend, you can not
# import it unless there's a FITS backend present, thanks to the
# design of sherpa.astro.io. So all the tests are behind a
# requires_fits decorator. The import could be done at the top of the
# file (within a check for ImportError) but leave as is.
#

from sherpa.utils.testing import requires_fits


@requires_fits
def test_str_singleton():
    """stringification: single value"""

    from sherpa.astro.io.meta import Meta
    store = Meta()

    # Could loop over this, but a bit easier to check the string
    # output this way.
    #
    store['key'] = ""
    assert str(store) == '\n key           = '

    store['key'] = "  "
    assert str(store) == '\n key           =   '

    store['key'] = " a string "
    assert str(store) == '\n key           =  a string '

    store['key'] = 23
    assert str(store) == '\n key           = 23'

    store['key'] = 23.0
    assert str(store) == '\n key           = 23.0'

    store['key'] = False
    assert str(store) == '\n key           = False'

    # Now some special cases
    for val in [None, "None", "NONE", "none"]:
        store['key'] = val
        assert str(store) == ''


@requires_fits
def test_str_multi():
    """Multiple keys are displayed as expected"""

    from sherpa.astro.io.meta import Meta
    store = Meta()
    store['Xkey'] = 'X X'
    store['xkey'] = ' y  y'
    store['a'] = 23
    store['INFILE'] = 'none'
    store['outfile'] = '/tmp/b.fits'

    lines = str(store).split('\n')
    assert len(lines) == 5
    assert lines[0] == ''
    assert lines[1] == ' Xkey          = X X'
    assert lines[2] == ' a             = 23'
    assert lines[3] == ' outfile       = /tmp/b.fits'
    assert lines[4] == ' xkey          =  y  y'
