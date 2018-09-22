#
#  Copyright (C) 2017, 2018  Smithsonian Astrophysical Observatory
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

from sherpa.astro import ui
from sherpa.astro.utils import filter_resp
from sherpa.utils.testing import requires_data, requires_fits


# See https://github.com/sherpa/sherpa/issues/405
def test_filter_resp_nochans():
    """What happens if no channels?

    It could be an error, or empty arrays could be returned.
    """

    # fake up an RMF
    ngrp = [1, 1]
    f_chan = [1, 2]
    n_chan = [1, 1]
    matrix = [1.0, 1.0]
    offset = 1

    with pytest.raises(ValueError) as excinfo:
        filter_resp([], ngrp, f_chan, n_chan, matrix, offset)

    emsg = "There are no noticed channels"
    assert str(excinfo.value) == emsg


# See https://github.com/sherpa/sherpa/issues/405
@requires_data
@requires_fits
def test_rmf_filter_no_chans(make_data_path):
    """This is not really a unit test but placed here as it
    is related to test_filter_resp_nochans.
    """

    rmffile = make_data_path('3c273.rmf')
    rmf = ui.unpack_rmf(rmffile)
    with pytest.raises(ValueError) as excinfo:
        rmf.notice([])

    emsg = "There are no noticed channels"
    assert str(excinfo.value) == emsg
