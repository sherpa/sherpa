#
#  Copyright (C) 2023, 2024
#  MIT
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
'''Helper functions for testing the plots

'''

from sherpa import plot
# This not the most elegant solution, but it makes sense to
# have an independent check here and not rely on what is
# done in the sherpa.plot.__init__ module, so that the tests
# stay independent of that particular implementation.
try:
    from sherpa.plot.pylab_backend import PylabBackend
    HAS_PYLAB = True
except ModuleNotFoundError:
    HAS_PYLAB = False

try:
    from sherpa.plot.bokeh_backend import BokehBackend
    HAS_BOKEH = True
except ModuleNotFoundError:
    HAS_BOKEH = False


__all__ = ('check_empty', 'check_full')


def check_empty(r, summary, nsummary=0):
    """Is this an 'empty' response?

    Parameters
    ----------
    r : str
        The HTML response.
    summary : str
        The summary string to look for.
    nsummary : int, optional
        The number of lines in the detailed table.
    """

    if (HAS_PYLAB and isinstance(plot.backend, PylabBackend)) or \
        (HAS_BOKEH and isinstance(plot.backend, BokehBackend)):
        assert r is None
        return

    assert r is not None
    assert f"<summary>{summary} ({nsummary})</summary>" in r


def check_full(r, summary, label='', title='', nsummary=0, test_other=None):
    """Is this a 'full' response?

    This test runs different checks for functional (pylab and bokeh)
    backends and dummy backends that do no produce graphics.

    Parameters
    ----------
    r : str
        The HTML response.
    summary : str
        The summary string to look for.
    label : str
        The label string to look for.
    title : str
        The title string to look for.
    nsummary : int, optional
        The number of lines in the detailed table.
    test_other : list of str, optional
        Usually dummy backends get tested for a list of
        strings including the label and title.
        However, in some cases, that is not appropriate
        and a list of strings can be passed to check
        one the existence of all those literal strings
        in the HTML response.
    """

    assert r is not None

    if (HAS_PYLAB and isinstance(plot.backend, PylabBackend)) or \
        (HAS_BOKEH and isinstance(plot.backend, BokehBackend)):
        assert f"<summary>{summary}</summary>" in r

        if (HAS_PYLAB and isinstance(plot.backend, PylabBackend)):
            assert "<svg " in r
        elif (HAS_BOKEH and isinstance(plot.backend, BokehBackend)):
            assert "BokehJS library" in r
        return

    assert '<svg ' not in r
    assert "BokehJS library" not in r

    if test_other is None:
        assert f"<summary>{summary} ({nsummary})</summary>" in r
        assert f'<div class="dataval">{label}</div>' in r
        assert f'<div class="dataval">{title}</div>' in r
    else:
        for line in test_other:
            assert line in r
