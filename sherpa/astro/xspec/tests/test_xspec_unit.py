#
#  Copyright (C) 2016  Smithsonian Astrophysical Observatory
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

# Supplement test_xspec.py with py.test tests
#

import pytest

from numpy.testing import assert_almost_equal

from sherpa.utils import requires_data, requires_fits, requires_xspec


@requires_data
@requires_fits
@requires_xspec
def test_read_xstable_model(make_data_path):
    """Limited test (only one file).

    Evaluation tests using this model are in
    sherpa.astro.xspec.tests.test_xspec.
    """

    from sherpa.astro import xspec

    path = make_data_path('xspec-tablemodel-RCS.mod')
    tbl = xspec.read_xstable_model('bar', path)

    assert tbl.name == 'bar'
    assert isinstance(tbl, xspec.XSTableModel)
    assert tbl.addmodel

    assert len(tbl.pars) == 4
    assert tbl.pars[0].name == 'tau'
    assert tbl.pars[1].name == 'beta'
    assert tbl.pars[2].name == 't'
    assert tbl.pars[3].name == 'norm'

    assert_almost_equal(tbl.tau.val, 1)
    assert_almost_equal(tbl.tau.min, 1)
    assert_almost_equal(tbl.tau.max, 10)

    assert_almost_equal(tbl.beta.val, 0.1)
    assert_almost_equal(tbl.beta.min, 0.1)
    assert_almost_equal(tbl.beta.max, 0.5)

    assert_almost_equal(tbl.t.val, 0.1)
    assert_almost_equal(tbl.t.min, 0.1)
    assert_almost_equal(tbl.t.max, 1.3)

    assert_almost_equal(tbl.norm.val, 1)
    assert_almost_equal(tbl.norm.min, 0)
    assert_almost_equal(tbl.norm.max, 1e24)

    for p in tbl.pars:
        assert not(p.frozen)
