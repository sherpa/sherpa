#
#  Copyright (C) 2025
#  Smithsonian Astrophysical Observatory
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

from sherpa.astro.data import DataIMG, DataIMGInt
from sherpa.models import Gauss2D


@pytest.mark.parametrize("cls", [DataIMG, DataIMGInt])
def test_create_dataIMG_direct(cls):
    """Create an array and flatten it. Then check that we can get the
    image back without wrong ordering from row-first vs column-first confusion."""
    # Note that this say x1, x0 (in that order!)
    x1, x0 = np.mgrid[20:30, 5:20]
    datashape = x0.shape
    y = np.sqrt((x0 - 10)**2 + (x1 - 31)**2)
    x0 = x0.flatten()
    x1 = x1.flatten()
    y = y.flatten()
    if cls is DataIMG:
        image = cls("bimage", x0=x0, x1=x1, y=y, shape=datashape)
    elif cls is DataIMGInt:
        image = cls("bimage",
                    x0lo=x0-0.5, x1lo=x1-0.5, x0hi=x0+0.5, x1hi=x1+0.5,
                    y=y, shape=datashape)
    else:
        raise ValueError("Unexpected class")

    out = image.get_img()
    assert out.shape == datashape
    # Check it's not totally disordered.
    assert np.unravel_index(np.argmin(out), out.shape) == (9, 5)

    # Could express this in array notation, but this is readable and simple.
    # If array ordering (n,m) was reversed to (m,n) it would scramble
    # values and show here.
    for i in range(9):
        np.all(out[i, :] > out[i+1, :])
    # Unlikely that the contents of an array changes, but check just in case.
    assert out.min() == pytest.approx(2.0)

    indep = image.get_indep()
    if cls is DataIMG:
        # Testing the order here.
        assert indep[0][:15] == pytest.approx(np.arange(5, 20))
        assert indep[1][:15] == pytest.approx(20)
    elif cls is DataIMGInt:
        # Testing the order here.
        assert indep[0][:15] == pytest.approx(np.arange(5, 20) - 0.5)
        assert indep[1][:15] == pytest.approx(20 - 0.5)
        assert indep[2][:15] == pytest.approx(np.arange(5, 20) + 0.5)
        assert indep[3][:15] == pytest.approx(20 + 0.5)


@pytest.mark.parametrize("cls", [DataIMG, DataIMGInt])
@pytest.mark.parametrize("use_2d_array", [True, False])
def test_create_dataIMG_from_2darray(cls, use_2d_array):
    """Similar to test_create_dataIMG_direct but using a 2D array.

    In this case, we all define a model using the 2D array, so that we
    can check the evaluated values as well.
    """
    rng = np.random.default_rng(12345)

    x0 = rng.normal(size=1000, loc=21.6, scale=0.03)
    x1 = rng.normal(size=1000, loc=0.2, scale=0.03)
    x0range = np.arange(21.4, 21.7, 0.01)
    x1range = np.arange(-.25, .250001, 0.01)
    hist, x0edges, x1edges = np.histogram2d(x0, x1, bins=(x0range, x1range))
    x0lo, x1lo = np.meshgrid(x0edges[:-1], x1edges[:-1])
    x0hi, x1hi = np.meshgrid(x0edges[1:], x1edges[1:])

    # Make a model
    gline = Gauss2D("gline")
    gline.fwhm = 0.05
    gline.fwhm.freeze()
    gline.ampl = 115
    gline.xpos = 21.6
    gline.ypos.frozen = False
    gline.ypos = 0.2
    gline.xpos.frozen = False

    if cls == DataIMG:
        x01 = {'x0': x0lo, 'x1': x1lo}
    elif cls == DataIMGInt:
        x01 = {'x0lo': x0lo, 'x1lo': x1lo,
                'x0hi': x0hi, 'x1hi': x1hi}
    else:
        raise ValueError("Unexpected class")

    if use_2d_array:
        if cls == DataIMG:
            x01 = {'x0': x0edges[:-1], 'x1': x1edges[:-1]}
        elif cls == DataIMGInt:
            x01 = {'x0_bounds': x0edges, 'x1_bounds': x1edges}
        else:
            raise ValueError("Unexpected class")

        image = cls.from_2d_array("image", y=hist, **x01)
    else:
        if cls == DataIMG:
            x01 = {'x0': x0lo, 'x1': x1lo}
        elif cls == DataIMGInt:
            x01 = {'x0lo': x0lo, 'x1lo': x1lo,
                    'x0hi': x0hi, 'x1hi': x1hi}
        else:
            raise ValueError("Unexpected class")

        image = cls("image",
                    **{k: v.flatten() for k, v in x01.items()},
                    y=hist.T.flatten(), shape=hist.T.shape,
                    staterror=np.sqrt(hist).T.flatten())
    out = image.get_img(gline)
    assert out[0].shape == (50, 30)
    assert out[1].shape == (50, 30)
    indep = image.get_indep()
    # Check the first 30 elements.
    # That's enough to see a row-first vs. column-first ordering.
    assert indep[0][:30] == pytest.approx(x0range[:-1])
    assert indep[1][:30] == pytest.approx(-0.25)
    if cls == DataIMGInt:
        assert indep[2][:30] == pytest.approx(x0range[1:])
        assert indep[3][:30] == pytest.approx(-0.24)


@pytest.mark.parametrize("cls", [DataIMG, DataIMGInt])
def test_from2D_wrong_dimension(cls):
    """Check that we catch wrong dimension in from_2d_array"""
    y = np.zeros((10, 20, 30))
    x0 = np.arange(11)
    x1 = np.arange(21)
    with pytest.raises(ValueError) as exc:
        if cls == DataIMG:
            cls.from_2d_array("image", y=y, x0=x0[:-1], x1=x1[:-1])
        elif cls == DataIMGInt:
            cls.from_2d_array("image", y=y, x0_bounds=x0, x1_bounds=x1)
        else:
            raise ValueError("Unexpected class")
    assert "Expected 2D array for y, got 3D array instead." in str(exc.value)

@pytest.mark.parametrize("cls", [DataIMG, DataIMGInt])
def test_wrong_order_x0x1(cls):
    """Check that we catch wrong dimension in from_2d_array"""
    y = np.zeros((10, 20))
    x0 = np.arange(21)
    x1 = np.arange(11)
    with pytest.raises(ValueError,
                       match=r"Length of x0[\s\S]+ \(10\)"):
        if cls == DataIMG:
            cls.from_2d_array("image", y=y, x0=x0[:-1], x1=x1[:-1])
        elif cls == DataIMGInt:
            cls.from_2d_array("image", y=y, x0_bounds=x0, x1_bounds=x1)
        else:
            raise ValueError("Unexpected class")