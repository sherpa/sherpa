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

"""
Should these tests be moved to test_astro_session.py?

This is almost a copy of sherpa/ui/tests/test_ui_unit.py, but
some of the answers are different. Is there a better way than
this duplication?
"""

from sherpa.astro import ui


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
    # but do not enforce these are the only values. This is
    # a slightly-different return list to the non-astro version.
    #
    samplers = ui.list_samplers()
    for expected in ['mh', 'metropolismh', 'pragbayes', 'fullbayes']:
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
