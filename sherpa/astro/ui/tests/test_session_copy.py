#
#  Copyright (C) 2016, 2017, 2021  Smithsonian Astrophysical Observatory
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
This is a copy of sherpa/ui/tests/test_session changed to use
the astro Session object.

If we want to do this then some way of factoring out the common
code and tests is likely needed. I have started this, but not
in a principled manner.
"""

from numpy.testing import assert_array_equal

import pytest

from sherpa.astro.ui.utils import Session as AstroSession
from sherpa.ui.utils import Session
from sherpa.utils.err import ClobberErr, IOErr
from sherpa.utils.testing import requires_plotting


TEST = [1, 2, 3]
TEST2 = [4, 5, 6]


# bug #303
#
# Note: this test could be run even if plotting is not available,
# since in that case the initial get_data_plot_prefs() call will
# return an empty dictionary rather than a filled one.
#
@requires_plotting
def test_set_log():
    session = AstroSession()
    assert not session.get_data_plot_prefs()['xlog']
    assert not session.get_data_plot_prefs()['ylog']
    session.set_xlog()
    assert session.get_data_plot_prefs()['xlog']
    session.set_ylog()
    assert session.get_data_plot_prefs()['ylog']
    session.set_xlinear()
    assert not session.get_data_plot_prefs()['xlog']
    session.set_ylinear()
    assert not session.get_data_plot_prefs()['ylog']


# bug #262
def test_list_ids():
    session = AstroSession()
    session.load_arrays(1, TEST, TEST)
    session.load_arrays("1", TEST, TEST2)

    # order of 1 and "1" is not determined
    assert {1, "1"} == set(session.list_data_ids())
    assert_array_equal(TEST2, session.get_data('1').get_dep())
    assert_array_equal(TEST, session.get_data(1).get_dep())


# bug #297
def test_save_restore(tmpdir):
    outfile = tmpdir.join("sherpa.save")
    session = AstroSession()
    session.load_arrays(1, TEST, TEST2)
    session.save(str(outfile), clobber=True)
    session.clean()
    assert set() == set(session.list_data_ids())

    session.restore(str(outfile))
    assert {1, } == set(session.list_data_ids())
    assert_array_equal(TEST, session.get_data(1).get_indep()[0])
    assert_array_equal(TEST2, session.get_data(1).get_dep())


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_save_clobber_check(session, tmp_path):
    """save does not clobber"""

    out = tmp_path / 'save.file'
    out.write_text('x')

    s = session()
    with pytest.raises(ClobberErr):
        s.save(str(out))

    assert out.read_text() == 'x'


# Apparently sherpa.ui.utils does not have save_all
@pytest.mark.parametrize("session", [AstroSession])
def test_save_all_clobber_check(session, tmp_path):
    """save_all does not clobber"""

    out = tmp_path / 'save.file'
    out.write_text('x')

    s = session()
    with pytest.raises(ClobberErr):
        s.save_all(str(out))

    assert out.read_text() == 'x'


@pytest.mark.parametrize("session", [Session, AstroSession])
def test_load_template_not_binary(session, tmp_path):
    """load_template_model doesn't like binary files"""

    # Creating a file which Sherpa thinks is a binary file is
    # an art.
    #
    out = tmp_path / 'save.file'
    out.write_bytes(b'x\0y')

    s = session()
    with pytest.raises(IOErr) as exc:
        s.load_template_model('xmodel', str(out))

    assert str(exc.value).startswith("file '")
    assert str(exc.value).endswith("' does not appear to be ASCII")
