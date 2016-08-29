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
import pytest
from sherpa.astro.utils._region import region_mask, Region


region_string = 'Circle(256 ,256 , 25)'
normalized_string = 'Circle(256,256,25)'
x = [1, 2, 3]
y = [1, 2, 3]
xm = [254, 255, 256]
ym = [254, 255, 256]


def test_region_string():
    r = Region(region_string)
    assert normalized_string == str(r)


def test_region_mask_string():
    r, m = region_mask(None, region_string, x, y, False, False)
    assert normalized_string == str(r)
    assert all(val == 0 for val in m)


def test_region_mask_string_match():
    r, m = region_mask(None, region_string, xm, ym, False, False)
    assert normalized_string == str(r)
    assert all(val != 0 for val in m)


def test_region_mask_file(tmpdir):
    f = tmpdir.join("reg")
    f.write(normalized_string)
    r, m = region_mask(None, str(f), x, y, True, False)
    assert normalized_string == str(r)
    assert all(val == 0 for val in m)


def test_region_non_string():
    with pytest.raises(TypeError):
        Region(3)


def test_region_null():
    with pytest.raises(TypeError):
        Region(None)


def test_region_mask_non_string():
    with pytest.raises(TypeError):
        region_mask(None, 3, x, y, False, False)


def test_region_mask_null():
    with pytest.raises(TypeError):
        region_mask(None, None, x, y, False, False)
