#
#  Copyright (C) 2021  Smithsonian Astrophysical Observatory
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
Split out testing of parts of sherpa.utils.testing to here.
"""

import os

import pytest

from sherpa.utils.testing import get_datadir, set_datadir


def test_set_datadir_does_not_exist(tmp_path):
    """We error out when the directory does not exist."""

    d = tmp_path / "does-not-exist"

    with pytest.raises(OSError) as exc:
        set_datadir(d)

    assert str(exc.value).startswith('datadir=')
    assert str(exc.value).endswith('/does-not-exist is empty or not a directory')


def test_set_datadir_is_empty(tmp_path):
    """We error out when the directory is empty."""

    d = tmp_path

    with pytest.raises(OSError) as exc:
        set_datadir(d)

    assert str(exc.value).startswith('datadir=')
    assert str(exc.value).endswith(' is empty or not a directory')


def test_set_datadir_is_not_a_directory(tmp_path):
    """We error out when the directory is actually a file."""

    p = tmp_path / "is-a-file"
    p.write_text('something')

    with pytest.raises(OSError) as exc:
        set_datadir(p)

    assert str(exc.value).startswith('datadir=')
    assert str(exc.value).endswith('/is-a-file is empty or not a directory')


def test_set_datadir():
    """Check we can set it, then restore the previous setting."""

    # Since set_datadir does not accept None, and I don't want to deal with
    # the low-level internals of the testing module, skip this set if
    # the directory is not set.
    #
    opath = get_datadir()
    if opath is None:
        pytest.skip("datadir is not set")

    # Ideally we'd wrap this up so that a problem would not bleed into
    # other tests, but this would be a fairly serious problem so do not
    # bother at this time.
    #
    # Assume the current directory is a reasonable directory to use.
    #
    npath = os.getcwd()
    assert opath != npath  # just a sanity check
    set_datadir(npath)

    assert get_datadir() == npath
    set_datadir(opath)
    assert get_datadir() != npath
