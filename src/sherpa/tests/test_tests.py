#
#  Copyright (C) 2023  MIT
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
'''Make sure we tests the test fixtures.'''
import pytest

def test_check_str_fails(check_str):
    '''The tests that use check_str make sure that it passes,
    but what if it were to pass everything?
    So, add a test here to make sure it actually fails when it should.'''
    with pytest.raises(AssertionError, match="^assert 'foo = "):
        check_str('foo = 1.23456', ['foo = 1.230000'])


def test_check_str_fails_with_float_cmp(check_str):
    # The match value may need updating if pytest.approx changes it's
    # assertion error.
    with pytest.raises(AssertionError, match="^assert 1.23456 "):
        check_str('foo = 1.23456', ['foo = 1.230000   # doctest: +FLOAT_CMP'])


def test_check_str_fails_with_float_cmp_number_only(check_str):
    # The match value may need updating if pytest.approx changes it's
    # assertion error.
    with pytest.raises(AssertionError, match="^assert 1.23456 "):
        check_str('1.23456', ['1.230000   # doctest: +FLOAT_CMP  '])
