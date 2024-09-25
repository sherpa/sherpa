#
#  Copyright (C) 2007, 2015, 2016, 2018, 2020, 2021, 2022
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

import os.path

import pytest

import sherpa
from sherpa import ui
from sherpa.utils.testing import requires_data


def test_has_version():
    """Just check we have a version setting.

    It's not obvious how we should check the version string - it
    can be 0+untagged.blah or a.b.c+n.... so it doesn't seem
    worth checking.
    """
    assert isinstance(sherpa.__version__, str)
    assert sherpa.__version__ != ''


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


def test_citation_hardcoded(tmp_path):
    """Check citation works for hardcoded queries."""

    citefile = tmp_path / 'cite.txt'
    with citefile.open(mode='w') as fh:
        sherpa.citation('4.8.0', filename=fh, clobber=True)

    cts = citefile.read_text('ascii').split('\n')
    assert cts[0].startswith('You are using Sherpa ')
    assert cts[1] == ''
    assert cts[2] == 'Sherpa 4.8.0 was released on January 27, 2016.'
    assert cts[3] == ''
    assert cts[4] == '@software{sherpa_2016_45243,'
    assert len(cts) == 86


@pytest.mark.zenodo
def test_todo_latest_success(tmp_path):
    """Check Zenodo knows about the current version.

    Since we don't know what the actual text should be, this
    is not a full test, but enough to know that the
    information has been extracted.
    """

    citefile = tmp_path / 'cite.txt'
    with citefile.open(mode='w') as fh:
        sherpa.citation('latest', filename=fh)

    cts = citefile.read_text('utf-8').split('\n')
    assert cts[0].startswith('You are using Sherpa ')
    assert cts[1] == ''
    assert cts[2].startswith('The latest release of Sherpa is ')
    assert cts[3].startswith('released on ')
    assert cts[4] == ''
    assert cts[5].startswith('@software{sherpa_')
    assert 'Please review the Zenodo Sherpa page at' in cts
    assert '@INPROCEEDINGS{2001SPIE.4477...76F,' in cts
    assert '@INPROCEEDINGS{2007ASPC..376..543D,' in cts


@pytest.mark.zenodo
def test_todo_version_fail(tmp_path):
    """Check Zenodo can not found an invalid version"""

    citefile = tmp_path / 'cite.txt'
    with citefile.open(mode='w') as fh:
        sherpa.citation('4.7.0', filename=fh)

    cts = citefile.read_text('utf-8').split('\n')
    assert cts[0].startswith('You are using Sherpa ')
    assert cts[1] == ''
    assert cts[2] == 'There was a problem retrieving the data from Zenodo:'
    assert cts[3] == 'Zenodo has no information for version 4.7.0.'
    assert cts[4] == ''
    assert cts[5] == 'Please review the Zenodo Sherpa page at'
    assert '@INPROCEEDINGS{2001SPIE.4477...76F,' in cts
    assert '@INPROCEEDINGS{2007ASPC..376..543D,' in cts
