#
#  Copyright (C) 2016, 2018, 2020  Smithsonian Astrophysical Observatory
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

# At present this only contains a very basic test for fix 309.

import os

from sherpa.ui.utils import Session
from sherpa.testing import requires_data


@requires_data
def test_309(make_data_path):

    idval = 'bug309'

    # have values near unity for the data
    ynorm = 1e9

    session = Session()

    dname = make_data_path('load_template_with_interpolation-bb_data.dat')

    session.load_data(idval, dname)
    session.get_data(idval).y *= ynorm

    indexname = 'bb_index.dat'
    datadir = make_data_path('')

    # Need to load the data from the same directory as the index
    basedir = os.getcwd()
    os.chdir(datadir)
    try:
        session.load_template_model('bbtemp', indexname)
    finally:
        os.chdir(basedir)

    bbtemp = session.get_model_component('bbtemp')
    session.set_source(idval, bbtemp * ynorm)

    session.set_method('gridsearch')
    session.set_method_opt('sequence', None)
    session.fit(idval)
