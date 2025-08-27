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

