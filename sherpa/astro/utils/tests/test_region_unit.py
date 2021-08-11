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

import numpy as np

import pytest

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

    Note the region gets normalized.
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


def test_region_file(tmp_path):
    """Can we create a region from a file?"""
    f = tmp_path / 'reg'
    f.write_text(normalized_string)

    r = Region(str(f), True)
    assert normalized_string == str(r)


def test_region_file_invalid(tmp_path):
    """What happens if the file is invalid.

    Note that the file parser apparently fails if
    there are spaces, so I use region_string here
    in case the parser ever gets improved, in which
    case this test will need a new invalid file.
    """
    f = tmp_path / 'src.reg'
    f.write_text(region_string)

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


def test_region_mask_nomatch_file(tmp_path):
    """It shoudn't matter where the region comes from"""
    f = tmp_path / 'reg.reg'
    f.write_text(normalized_string)

    r = Region(str(f), True)
    m = r.mask(x, y)
    assert all(val == 0 for val in m)


def test_region_mask_match_file(tmp_path):
    """It shoudn't matter where the region comes from"""
    f = tmp_path / 'reg.src'
    f.write_text(normalized_string)

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


def check_region_invert_compare(reg):
    """Check that [field()]-rect(..) is filtering correctly."""

    y, x = np.mgrid[99:108, 100:106]
    x = x.flatten()
    y = y.flatten()

    m = reg.mask(x, y)

    # For boxes it's <= for the upper limit, not <
    expected = (x >= 100) & (x <= 104) & (y >= 100) & (y <= 106)
    expected = ~expected
    assert np.all(m == expected.flatten())

    # just check we are actually doing a filter...
    assert expected.size == 54
    assert expected.sum() == 19


def test_region_invert_compare_invert():
    """-r and field()-r should ideally filter the same way
    but they are potentially different.
    """

    shape = "rect(100, 100, 104, 106)"
    r = Region(shape)
    r.invert()
    check_region_invert_compare(r)


def test_region_invert_compare_steps():
    """-r and field()-r should ideally filter the same way
    but they are potentially different.
    """

    shape = "rect(100, 100, 104, 106)"
    r = Region("field()").combine(Region(shape), exclude=True)
    check_region_invert_compare(r)


def test_region_invert_compare_inone():
    """-r and field()-r should ideally filter the same way
    but they are potentially different.
    """

    shape = "rect(100, 100, 104, 106)"
    r = Region(f"field()-{shape}")
    check_region_invert_compare(r)


COMPLEX_REGION = ["circle(50,50,30)",
                  "ellipse(40,75,30,20,320)",
                  "-rotbox(30,30,10,5,45)"]


def check_complex_region_1245(reg):
    """Check the region filtering

    We do not check the string representation as different
    strings can represent the same filter
    """

    # make it rectangular just in case
    y, x = np.mgrid[0:101, 0:102]
    z = np.arange(x.size)
    mask = reg.mask(x.flatten(), y.flatten()) == 1

    # I have not checked this is correct, but it was based on data
    # created by image_data and analyzed by dax that seemed to make sense.
    #
    assert mask.sum() == 3722

    filt = z[mask]
    assert filt.min() == 2090
    assert filt.max() == 10133
    assert filt.sum() == 22558471


def check_complex_region_1245_no_exclude(reg):
    """check_complex_region_1245 but when the -rotbox is not applied
    """

    # make it rectangular just in case
    y, x = np.mgrid[0:101, 0:102]
    z = np.arange(x.size)
    mask = reg.mask(x.flatten(), y.flatten()) == 1

    # I have not checked this is correct, but it was based on data
    # created by image_data and analyzed by dax that seemed to make sense.
    #
    assert mask.sum() == 3757

    filt = z[mask]
    assert filt.min() == 2090
    assert filt.max() == 10133
    assert filt.sum() == 22671256


def test_complex_region_file_explicit(tmp_path):
    """Based on #1245.

    Here we are explicit in our region shape - that is
       a - c
       b - c

    although technically we could have used
       a - c
       b

    as c does not overlap a
    """

    regfile = tmp_path / 'test.reg'
    regstr = "\n".join([COMPLEX_REGION[0] + COMPLEX_REGION[2],
                        COMPLEX_REGION[1] + COMPLEX_REGION[2]])
    regfile.write_text(regstr)
    rfile = Region(str(regfile), True)

    assert str(rfile) == 'Circle(50,50,30)&!RotBox(30,30,10,5,45)|Ellipse(40,75,30,20,320)&!RotBox(30,30,10,5,45)'
    check_complex_region_1245(rfile)


def test_complex_region_file_implicit(tmp_path):
    """Based on #1245.

    Unlike test_complex_region_file_explicit the region
    file is just

       a
       b
       -c

    Note that this actualy behaves as if we've just ignored the "-c"
    filter (e.g. a dmcopy with this filter file ignores the -c
    expression) because c does not overlap b [or whatever logic
    is used by CIAO region expressions]

    """

    regfile = tmp_path / 'test.reg'
    regstr = "\n".join(COMPLEX_REGION)
    regfile.write_text(regstr)
    rfile = Region(str(regfile), True)

    assert str(rfile) == 'Circle(50,50,30)|Ellipse(40,75,30,20,320)&!RotBox(30,30,10,5,45)'
    check_complex_region_1245_no_exclude(rfile)


@pytest.mark.xfail
def test_complex_region_manual_explicit():
    """Based on #1245.

    Here we create the region like test_complex_region_file_explicit
    """

    rmanual = Region(COMPLEX_REGION[0] + COMPLEX_REGION[2])
    rmanual = rmanual.combine(Region(COMPLEX_REGION[1] + COMPLEX_REGION[2]))

    assert str(rmanual) == 'Circle(50,50,30)&!RotBox(30,30,10,5,45)&Ellipse(40,75,30,20,320)&!RotBox(30,30,10,5,45)'
    check_complex_region_1245(rmanual)


@pytest.mark.xfail
def test_complex_region_manual_implicit():
    """Based on #1245.

    Here we create the region like test_complex_region_file_implicit
    """

    rmanual = Region(COMPLEX_REGION[0])
    rmanual = rmanual.combine(Region(COMPLEX_REGION[1]))
    rmanual = rmanual.combine(Region(COMPLEX_REGION[2][1:]), exclude=True)

    assert str(rmanual) == 'Circle(50,50,30)&Ellipse(40,75,30,20,320)&!RotBox(30,30,10,5,45)'
    check_complex_region_1245(rmanual)
