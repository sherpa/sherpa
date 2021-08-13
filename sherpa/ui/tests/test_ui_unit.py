#
#  Copyright (C) 2017, 2018, 2020, 2021
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

"""
Should these tests be moved to test_session.py?

Note that this test is almost duplicated in
sherpa/astro/ui/tests/test_astro_ui_unit.py
"""

import numpy as np

import pytest

from sherpa import ui
from sherpa.utils.err import ArgumentTypeErr, IdentifierErr


# This is part of #397
#
def test_list_samplers():
    """Ensure list_samplers returns a list."""

    samplers = ui.list_samplers()

    assert isinstance(samplers, list)
    assert len(samplers) > 0


def test_list_samplers_contents():
    """Are the expected values included"""

    # Test that the expected values exist in this list,
    # but do not enforce these are the only values.
    #
    samplers = ui.list_samplers()
    for expected in ['mh', 'metropolismh']:
        assert expected in samplers


def test_all_has_no_repeated_elements():
    """Ensure __all__ does not contain repeated elements.

    It is not actually a problem if this happens, but it does
    indicate the possibility of confusion over what functionality
    the repeated symbol has (originally noticed with erf, which
    could be either sherpa.utils.erf or the ModelWrapper version
    of sherpa.models.basic.Erf). See
    https://github.com/sherpa/sherpa/issues/502
    """

    n1 = len(ui.__all__)
    n2 = len(set(ui.__all__))
    assert n1 == n2


@pytest.mark.parametrize("func", [ui.notice_id, ui.ignore_id])
@pytest.mark.parametrize("lo", [None, 1, "1:5"])
def test_check_ids_not_none(func, lo):
    """Check they error out when id is None

    There used to be the potential for different behavior depending
    if the lo argument was a string or not,hence the check.
    """

    with pytest.raises(ArgumentTypeErr) as exc:
        func(None, lo)

    assert str(exc.value) == "'ids' must be an identifier or list of identifiers"


@pytest.mark.parametrize("func", [ui.notice, ui.ignore])
@pytest.mark.parametrize("lo,hi", [(1, 5), (1, None), (None, 5), (None, None),
                                   ("1:5", None)])
def test_filter_no_data_is_an_error(func, lo, hi, clean_ui):
    """Does applying a filter lead to an error?

    This test was added because it was noted that an update for Python 3
    had lead to an error state not being reached.
    """

    with pytest.raises(IdentifierErr) as ie:
        func(lo, hi)

    assert str(ie.value) == 'No data sets found'


@pytest.mark.xfail
def test_save_filter_data1d(tmp_path, clean_ui):
    """Check save_filter [Data1D]"""

    x = np.arange(1, 11, dtype=np.int16)
    ui.load_arrays(1, x, x)

    ui.notice(2, 4)
    ui.notice(6, 8)

    outfile = tmp_path / "filter.dat"
    ui.save_filter(str(outfile))

    expected = [0, 1, 1, 1, 0, 1, 1, 1, 0, 0]

    d = ui.unpack_data(str(outfile), colkeys=['X', 'FILTER'])
    assert isinstance(d, ui.Data1D)
    assert d.x == pytest.approx(x)
    assert d.y == pytest.approx(expected)
    assert d.staterror is None
    assert d.syserror is None
