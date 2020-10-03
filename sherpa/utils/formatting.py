#
# Copyright (C) 2020  Smithsonian Astrophysical Observatory
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


"""Display representations of Sherpa objects.

This is aimed at IPython/Jupiter support but may be useful elsewhere.

"""

import contextlib
import html
import pkg_resources

import numpy as np


# Used by printoptions
DISPLAY_WIDTH = 80

# The CSS file for the Notebook HTML code
CSS_FILE_PATH = "/".join(("static", "css", "style.css"))
CSS_STYLE = pkg_resources.resource_string("sherpa", CSS_FILE_PATH).decode("utf8")


@contextlib.contextmanager
def printoptions(*args, **kwargs):
    """Temporarily over-ride the NumPy print options."""

    orig = np.get_printoptions()
    np.set_printoptions(*args, **kwargs)
    try:
        yield
    finally:
        np.set_printoptions(**orig)


def show_array(x):
    """Show array, possibly eliding elements.

    The idea is to show a "useful" representation of the data, and is
    aimed at supporting 1D and 2D arrays. At the moment this uses the
    str, rather than repr, representation (i.e. no surrounding
    'array(...)'.

    """

    x = np.asarray(x)

    # Take advantage of work by dask/xarray to come up with
    # useful values.
    #
    with printoptions(precision=6, threshold=200,
                      linewidth=DISPLAY_WIDTH):
        return str(x)


def clean_bracket(s):
    """Remove the outermost brackets for a model string

    Parameter
    ---------
    s : str
        A model expression (but can be any text).

    Returns
    -------
    cleaned : s
         The outermost () is removed (if there is one).

    Examples
    --------

    >>> clean_bracket('(xsphabs.gal * powlaw1d.pl)')
    'xsphabs.gal * powlaw1d.pl'

    """
    if s.startswith('(') and s.endswith(')'):
        return s[1:-1]

    return s


def html_from_sections(obj, ls):
    """Convert the list of sections to HTML.

    This adds the stylesheet, since the sections are assumed to
    already be converted to HTML. It also creates a fall-through
    version for untrusted input, using the object's repr.

    Parameters
    ----------
    obj
        The object
    ls : sequence of str
        The HTML blocks used to combine.

    Returns
    -------
    html : str
        The HTML, including style blocks.

    Notes
    -----
    This is low-level, and so may change significantly.

    """

    out = '<style>'
    out += CSS_STYLE
    out += '</style>'

    # fall through for non-CSS-respecting displays
    #
    out += '<div class="sherpa-text-fallback">'
    out += html.escape(repr(obj))
    out += '</div>'

    out += '<div hidden class="sherpa">'
    out += ''.join(ls)
    out += '</div>'
    return out


def html_section(rows, summary=None, open_block=False):
    """Create a section (a details block) for HTML.

    Parameters
    ----------
    rows : sequence of (string, value)
        The data to display, identified by a label and then the value
        to display.
    summary : str or None, optional
        If not None, the summary for the details block
    open_block : bool, optional
        If True then the details block is opened (it is closed by
        default).

    Notes
    -----

    There is no attempt to catch values that could lead to surprising
    output (e.g. injection of HTML).

    """

    out = '<details'
    if open_block:
        out += ' open'
    out += '>'

    if summary is not None:
        out += '<summary>{} ({})</summary>'.format(summary,
                                                   len(rows))

    out += '<div class="datavals">'
    for l, r in rows:
        out += '<div class="dataname">{}</div>'.format(l)

        if isinstance(r, np.ndarray):
            r = show_array(r)

        out += '<div class="dataval">{}</div>'.format(r)

    out += '</div></details>'
    return out


def html_table(header, rows, summary=None, caption=None,
               rowcount=True, classname=None, open_block=True):
    """Create a section (a details block) for a table.

    Parameters
    ----------
    header : list of str
        The column headers.
    rows : list of list of (string, value)
        The data to display, identified by a label and then the value
        to display. The number of columns should match the header.
    summary : str or None, optional
        If not None, the summary for the details block
    caption : str or None, optional
        The caption, if set.
    rowcount : bool, optional
        If set then the number of rows is added to the summary
        element.
    classname : str or None, optional
        The class name for the table.
    open_block : bool, optional
        If True then the details block is opened (it is closed by
        default).

    Notes
    -----

    There is no attempt to catch values that could lead to surprising
    output (e.g. injection of HTML).

    """

    out = '<details'
    if open_block:
        out += ' open'
    out += '>'

    if summary is not None:
        out += '<summary>{}'.format(summary)
        if rowcount:
            out += ' ({})'.format(len(rows))

        out += '</summary>'

    out += '<table'
    if classname is not None:
        out += ' class="{}"'.format(classname)

    out += '>'

    if caption is not None:
        out += '<caption>{}</caption>'.format(caption)

    out += '<thead><tr>'
    for h in header:
        out += '<th>{}</th>'.format(h)

    out += '</tr></thead><tbody>'

    for row in rows:
        out += '<tr>'
        for r in row:
            out += '<td>{}</td>'.format(r)

        out += '</tr>'

    out += '</tbody></table></details>'
    return out


def html_svg(svg, summary, open_block=True):
    """Create a section (a details block) for HTML from a SVG.

    Parameters
    ----------
    svg : str
        The plot to display.
    summary : str
        The the summary for the details block
    open_block : bool, optional
        If True then the details block is opened (it is open by
        default).

    """

    out = '<details'
    if open_block:
        out += ' open'
    out += '>'

    out += '<summary>{}</summary>'.format(summary)
    out += svg
    out += '</details>'
    return out
