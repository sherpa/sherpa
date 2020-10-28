#
#  Copyright (C) 2020  Smithsonian Astrophysical Observatory
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

import numpy

import pytest

from sherpa.utils import akima


def test_akima():
    x = numpy.sort(numpy.random.random(10) * 10)
    y = numpy.random.normal(0.0, 0.1, size=len(x))
    myakima = akima.Akima(x, y)

    num = 3
    for ii in numpy.linspace(0.01, 0.05, num):
        xx = numpy.arange(x[0], x[-1], 0.05)
        z = akima.interpolate(x, y, xx)
        akima_z = myakima(xx)
        assert akima_z == pytest.approx(z)
