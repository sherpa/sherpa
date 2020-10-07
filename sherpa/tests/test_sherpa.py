#
#  Copyright (C) 2007, 2015, 2016, 2018, 2020
#          Smithsonian Astrophysical Observatory
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

import os.path

import pytest

import sherpa
from sherpa import ui
from sherpa.utils.testing import requires_data


def test_include_dir():
    incdir = os.path.join(sherpa.get_include(), 'sherpa')
    assert os.path.isdir(incdir)


@requires_data
def test_not_reading_header_without_comment(make_data_path):
    with pytest.raises(ValueError):
        ui.load_data(make_data_path('agn2'))


@requires_data
def test_require_float(make_data_path):
    with pytest.raises(ValueError):
        ui.load_data(make_data_path('agn2'))


@requires_data
def test_reading_floats(make_data_path):
    ui.load_data(make_data_path('agn2_fixed'))


@requires_data
def test_reading_strings(make_data_path):
    ui.load_data(make_data_path('table.txt'), require_floats=False)
