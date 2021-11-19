#
#  Copyright (C) 2021
#  MIT
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

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.astro import io
from sherpa.astro import ui


@requires_fits
@requires_data
def test_backend_switch(make_data_path):
    '''Test that we can switch backends.

    Of course, switching only makes sense if there is more than one backend
    to try, so this is listed as "requires_fits".
    '''
    pha = make_data_path("3c273.pi")
    ui.clean()
    ui.load_pha(pha)
    assert ui.get_data().name.endswith('3c273.pi')

    orig_backend = io.backend

    import sherpa.astro.io.dummy_backend
    io.backend = sherpa.astro.io.dummy_backend

    with pytest.raises(NotImplementedError,
                       match="No usable I/O backend was imported."):
         ui.load_pha(pha)

    io.backend = orig_backend
    infile = make_data_path("9774.pi")
    ui.load_pha(infile)
    assert ui.get_data().name.endswith('9774.pi')
