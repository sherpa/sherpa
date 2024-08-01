#
#  Copyright (C) 2013, 2015, 2018, 2021  Smithsonian Astrophysical Observatory
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

from sherpa.utils.testing import requires_data, requires_fits
from sherpa.astro import ui


# bug 12784
@requires_data
@requires_fits
def test_filter_nan(run_thread, clean_astro_ui):
    run_thread('filternan')
    assert not numpy.isnan(ui.get_fit_results().statval)
