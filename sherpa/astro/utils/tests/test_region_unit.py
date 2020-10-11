#
#  Copyright (C) 2016, 2020  Smithsonian Astrophysical Observatory
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

import numpy as np

from sherpa.astro.utils._region import Region


# Use a non-rotational shape so we can swap the
# xm and ym values and see a difference in the mask.
region_string = 'ellipse(256 ,256 , 25, 10, 0)'
normalized_string = 'Ellipse(256,256,25,10,0)'
x = [1, 2, 3]
y = [1, 2, 3]
xm = [254, 275, 256]
ym = [254, 255, 256]


def test_region_empty():
    """The empty region"""
    r = Region()
    assert '' == str(r)


def test_region_string():
    """We can create with a string.

    Note the regon gets normalized.
    """
    r = Region(region_string)
    assert normalized_string == str(r)


def test_region_string_invalid():
    """Check we error out with invalid input"""
    wrong = 'circle 100 200 30'
    with pytest.raises(ValueError) as exc:
        Region(wrong)

    emsg = 'unable to parse region string: ' + wrong
    assert str(exc.value) == emsg


def test_region_file(tmpdir):
    """Can we create a region from a file?"""
    f = tmpdir.join('reg')
    f.write(normalized_string)

    r = Region(str(f), True)
    assert normalized_string == str(r)


def test_region_file_invalid(tmpdir):
    """What happens if the file is invalid.

    Note that the file parser apparently fails if
    there are spaces, so I use region_string here
    in case the parser ever gets improved, in which
    case this test will need a new invalid file.
    """
    f = tmpdir.join('reg')
    f.write(region_string)

    with pytest.raises(ValueError) as exc:
        Region(str(f), True)

    emsg = 'unable to read region from: {}'.format(f)
    assert str(exc.value) == emsg


@pytest.mark.parametrize("arg", [None, 3])
def test_region_invalid(arg):
    """What happens with some invalid arguments?"""
    with pytest.raises(TypeError):
        Region(arg)


def test_region_mask_nomatch():
    """All points are excluded"""
    r = Region(region_string)
    m = r.mask(x, y)
    assert all(val == 0 for val in m)


def test_region_mask_match():
    """All points are included"""
    r = Region(region_string)
    m = r.mask(xm, ym)
    assert all(val == 1 for val in m)


def test_region_mask_match_swapped():
    """All points are included

    This is just to check the ordering is handled
    correctly for test_region_mask_keywords
    """
    r = Region(region_string)
    m = r.mask(ym, xm)
    assert m == pytest.approx([1, 0, 1])


def test_region_mask_keywords1():
    """Check can use x0/x1 names"""
    r = Region(region_string)
    m = r.mask(x1=ym, x0=xm)
    assert m == pytest.approx([1, 1, 1])


def test_region_mask_keywords2():
    """Check can use x0/x1 names"""
    r = Region(region_string)
    m = r.mask(x1=xm, x0=ym)
    assert m == pytest.approx([1, 0, 1])


def test_region_mask_nomatch_file(tmpdir):
    """It shoudn't matter where the region comes from"""
    f = tmpdir.join("reg")
    f.write(normalized_string)

    r = Region(str(f), True)
    m = r.mask(x, y)
    assert all(val == 0 for val in m)


def test_region_mask_match_file(tmpdir):
    """It shoudn't matter where the region comes from"""
    f = tmpdir.join("reg")
    f.write(normalized_string)

    r = Region(str(f), True)
    m = r.mask(xm, ym)
    assert all(val == 1 for val in m)


def test_region_mask_null():
    """What happens if the region is empty?"""
    empty = Region()
    m = empty.mask(x, y)
    assert all(val == 0 for val in m)


def test_region_combine_nothing1():
    """Adding a region to nothing == region"""

    empty = Region()
    r = empty.combine(Region(region_string))
    assert normalized_string == str(r)


def test_region_ignore_nothing1():
    """Ignoring a region frm nothing == !region"""

    empty = Region()
    r = empty.combine(Region(region_string), True)
    assert '!' + normalized_string == str(r)


@pytest.mark.parametrize('ignore', [False, True])
def test_region_combine_nothing2(ignore):
    """Adding or ignoring nothing to a region == region"""

    orig = Region(region_string)
    r = orig.combine(Region(), ignore)
    assert normalized_string == str(r)


def test_region_combine():
    """Adding two regions together"""

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r = Region(r1).combine(Region(r2))

    r2 = r2.replace('rect', 'rectangle')
    expected = r1.capitalize() + '&' + r2.capitalize()
    assert expected == str(r)


def test_region_ignore_no_overlap():
    """Ignoring a region which doesn't overlap.

    It looks like the region lib doesn't try to simplify
    the combination of shapes.
    """

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r = Region(r1).combine(Region(r2), True)

    r2 = r2.replace('rect', 'rectangle')
    expected = r1.capitalize() + '&!' + r2.capitalize()
    assert expected == str(r)


def test_region_ignore_overlap():
    """Ignoring a region which does overlap."""

    r1 = 'circle(10,10,10)'
    r2 = 'rect(5,5,7,10)'
    r = Region(r1).combine(Region(r2), True)

    r2 = r2.replace('rect', 'rectangle')
    expected = r1.capitalize() + '&!' + r2.capitalize()
    assert expected == str(r)


def test_region_combine_combined():
    """Check multiple calls to region_combine"""

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r3 = 'ellipse(70,70,10,12,45)'

    rcomb = Region(r1).combine(Region(r2))
    rcomb = rcomb.combine(Region(r3))

    r2 = r2.replace('rect', 'rectangle')
    expected = '&'.join([x.capitalize() for x in [r1, r2, r3]])
    assert expected == str(rcomb)


def test_region_ignore_combined():
    """Check multiple calls to region_combine

    I am not sure I agree with the result but I just want
    to check the results.
    """

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r3 = 'ellipse(70,70,10,12,45)'

    rcomb = Region(r1).combine(Region(r2), True)
    rcomb = rcomb.combine(Region(r3))

    r2 = '!' + r2.replace('rect', 'rectangle').capitalize()
    expected = '&'.join([r1.capitalize(), r2, r3.capitalize()])
    assert expected == str(rcomb)


def test_region_combine_kwargs():
    """Adding two regions together with keyword arguments"""

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r = Region(r1).combine(exclude=0, region=Region(r2))

    r2 = r2.replace('rect', 'rectangle')
    expected = r1.capitalize() + '&' + r2.capitalize()
    assert expected == str(r)


def test_region_ignore_kwargs():
    """Removing a regions with keyword arguments"""

    r1 = 'circle(10,10,10)'
    r2 = 'rect(100,100,200,200)'
    r = Region(r1).combine(exclude=1, region=Region(r2))

    r2 = r2.replace('rect', 'rectangle')
    expected = r1.capitalize() + '&!' + r2.capitalize()
    assert expected == str(r)


@pytest.mark.parametrize("arg", [None, "circle(0,0,1)"])
def test_region_combine_invalid(arg):
    """Check we error out if not sent a region."""

    with pytest.raises(TypeError):
        Region().combine(arg)


def test_region_invert():
    """Can we invert a region?"""

    x = np.asarray([50, 95, 105])
    y = np.asarray([50, 105, 115])

    r = Region("box(100,100,10,20)")

    assert r.mask(x, y) == pytest.approx([0, 1, 0])

    r.invert()
    assert str(r) == "!Box(100,100,10,20)"

    assert r.mask(x, y) == pytest.approx([1, 0, 1])


def test_region_invert_empty():
    """Can we invert an empty region?"""

    x = np.asarray([50, 95, 105])
    y = np.asarray([50, 105, 115])

    r = Region()

    assert r.mask(x, y) == pytest.approx([0, 0, 0])

    r.invert()
    assert str(r) == ''

    assert r.mask(x, y) == pytest.approx([0, 0, 0])
