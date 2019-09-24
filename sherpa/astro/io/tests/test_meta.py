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

import pytest

from sherpa.astro.io.meta import Meta


def test_empty():
    """An empty store is empty."""

    store = Meta()
    assert len(store.keys()) == 0


def test_key_does_not_exist():
    """We know when a key does not exist"""

    store = Meta()
    assert not store.has_key("key")


@pytest.mark.parametrize("value", ["", " a string ", 23, 23.0, True, None])
def test_add_keyword(value):
    """Can add a keyword"""

    store = Meta()
    store['key'] = value
    assert store.keys() == ['key']
    assert store.has_key('key')
    assert len(store.values()) == 1

    svalue = store['key']
    assert svalue == value


def test_error_if_missing_key():
    """Raise the correct error on missing key"""

    store = Meta()
    store['key'] = 1

    with pytest.raises(KeyError):
        store['keykey']


def test_get_with_known_key():
    """get works on a known key"""

    store = Meta()
    store['key1'] = 2
    store['key2'] = 23

    assert store.get('key1') == 2


def test_get_with_unknown_key():
    """the default value is used if no key exists"""

    store = Meta()
    store['key1'] = 2
    store['key2'] = 23

    assert store.get('key3', "unknown") == "unknown"


def test_overwrite_keyword():
    """We can change a keyword"""

    store = Meta()
    store['val1'] = True
    store['val1'] = False
    assert not store['val1']


def test_del_keyword():
    """We can delete a keyword"""

    store = Meta()
    store['x1'] = 1
    store['u1'] = 2
    del store['x1']

    assert store.keys() == ['u1']
    assert store.values() == [2]


def test_del_missing_keyword():
    """We can not delete a keyword that doesn't exist"""

    store = Meta()
    store['x1'] = 1
    store['u1'] = 2
    with pytest.raises(KeyError):
        del store['y1']

    assert store.keys() == ['u1', 'x1']
    assert store.values() == [2, 1]


def test_keyword_order():
    """keys are sorted, not based on input order"""

    store = Meta()
    store['key2'] = 1
    store['key1'] = 2
    store['a'] = True
    store['A'] = False

    keys = store.keys()
    assert keys == ['A', 'a', 'key1', 'key2']

    vals = store.values()
    assert vals == [False, True, 2, 1]


def test_key_pop():
    """Can pop a known value"""

    store = Meta()
    store['key1'] = True
    store['key2'] = False

    assert store.pop('key1')
    assert store.keys() == ['key2']


def test_key_pop():
    """Can pop an unknown"""

    store = Meta()
    store['key1'] = True
    store['key2'] = False

    assert store.pop('key') is None
    assert store.keys() == ['key1', 'key2']


def test_str_empty():
    """stringification: empty"""

    s = str(Meta())
    assert s == ''


def test_str_singleton():
    """stringification: single value"""

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


def test_str_multi():
    """Multiple keys are displayed as expected"""

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
