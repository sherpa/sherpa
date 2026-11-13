#
# Copyright (C) 2020, 2023, 2026
# Smithsonian Astrophysical Observatory
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

This includes HTML and LaTeX formatting helpers.
"""

import contextlib
import html
from importlib.resources import files

import numpy as np



__all__ = ["DISPLAY_WIDTH", "printoptions",
           "show_array", "clean_bracket", "html_from_sections",
           "html_section", "html_table", "html_svg",
           "smallest_significant_order", "latex_val", "latex_val_errors",
           ]


DISPLAY_WIDTH = 80
"Width of the display in characters; used by `show_array`."

# The CSS file for the Notebook HTML code
#
CSS_FILE_PATH = files("sherpa").joinpath("static", "css", "style.css")
CSS_STYLE = CSS_FILE_PATH.read_text(encoding="utf-8")


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


def smallest_significant_order(values: np.ndarray) -> int:
    """Return the smallest significant order of magnitude in the array.

    Parameters
    ----------
    values : array
        The array of values to check.

    Returns
    -------
    int
        The smallest significant order of magnitude in the array.

    Examples
    --------
    >>> smallest_significant_order(np.array([0.01, 0.03, 0.06, 0.1, 0.3, 0.6, 1]))
    -2
    """
    signif_mag = np.log10(np.abs(values))
    signif_mag = np.floor(signif_mag)
    return int(np.min(signif_mag))


def latex_val(val : float) -> str:
    r"""Format a value for LaTeX output.

    A small number is just passed through, but a number that is written in scientific
    notation is converted to a string with the \times 10^{...} notation.

    Parameters
    ----------
    val : float
        The value to format.

    Returns
    -------
    str
        The LaTeX formatted string.

    Examples
    --------
    >>> latex_val(12345)
    '$12345$'
    >>> latex_val(0.0000123)
    '$1.23 \\times 10^{-5}$'
    """
    fval = f'{val:g}'
    if not 'e' in fval:
        return f'${fval}$'
    return f'${fval.split("e")[0]} \\times 10^{{{int(fval.split("e")[1])}}}$'


def latex_val_errors(val : float, up : float | None, down : float | None) -> str:
    r"""Format a value with errors for LaTeX output.

    The precision of the printed value is determined by
    the magnitude of the smallest uncertainties.

    Parameters
    ----------
    val : float
        The value to format.
    up : float | None
        The upper error. If the error cannot be determined
        (e.g. it hits a hard boundary) then this can be None.
    down : float | None
        The lower error. This can be a positive or negative number to allows for two
        conventions: Either "the size of the down error" or " the size and direction
        of the down error". The output is independent of the sign of down.

    Returns
    -------
    str
        The LaTeX formatted string with the value and errors.

    Examples
    --------
    >>> latex_val_errors(12275, 100, 200)
    '$12300^{+100}_{-200}$'

    Note that 12275 has been rounded to 12300 to match the scale of the
    uncertainties (hundreds). The function will also handle scientific notation:

    >>> latex_val_errors(0.000123, 0.0000134, 0.00002123)
    '$12^{+1}_{-2} \\times 10^{-5}$'

    The function will also handle None values for the uncertainties:
    >>> latex_val_errors(0.000123, None, 0.00002123)
    '$12^{+\\infty}_{-2} \\times 10^{-5}$'

    The use of +/-infinity might not always be physically correct,
    e.g. a temperature cannot be negative, but this requires knowledge
    of the model that is not encoded in Sherpa. Values without constraints
    usually indicate fit problems and need to be investigated further or
    discussed by the user.
    """
    if down is not None:
        down = - abs(down)  # make sure down is negative
    valid = [v for v in [val, up, down] if v is not None]
    if len(valid) == 0:
        return '$?^{+\\infty}_{-\\infty}$'
    precision = smallest_significant_order(valid)
    fval = f'{np.round(val, -precision):g}'
    fup = '+\\infty' if up is None else f'{np.round(up, -precision):+g}'
    fdown = '-\\infty' if down is None else f'{np.round(down, -precision):+g}'
    if (not 'e' in fval) and (not 'e' in fup) and (not 'e' in fdown):
        return f'${fval}^{{{fup}}}_{{{fdown}}}$'

    fval = f'{val/10**precision:.0f}'
    fup = '+\\infty' if up is None else f'{up/10**precision:+.0f}'
    fdown = '-\\infty' if down is None else f'{down/10**precision:+.0f}'
    return f'${fval}^{{{fup}}}_{{{fdown}}} \\times 10^{{{precision}}}$'
