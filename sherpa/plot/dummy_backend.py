#
#  Copyright (C) 2007, 2015, 2020, 2021  Smithsonian Astrophysical Observatory
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
'''A dummy backend for plotting.

This backend implements only minimal functionality (some formatting of
strings as HTML or LaTeX which are usually used as axis labels), but no real
plotting capabilities. It is here to ensure that the `sherpa.plot` module can
be imported, even if no plotting backend is installed.
'''

import logging

from sherpa.utils import get_keyword_defaults
from sherpa.utils import formatting


__all__ = ('clear_window', 'plot', 'contour', 'point', 'set_subplot',
           'get_split_plot_defaults', 'get_confid_point_defaults',
           'get_plot_defaults', 'get_point_defaults', 'begin', 'end',
           'get_data_plot_defaults', 'get_model_plot_defaults',
           'get_fit_plot_defaults', 'get_resid_plot_defaults',
           'get_ratio_plot_defaults', 'get_contour_defaults',
           'get_data_contour_defaults', 'get_model_contour_defaults',
           'get_fit_contour_defaults', 'get_resid_contour_defaults',
           'get_ratio_contour_defaults', 'get_confid_plot_defaults',
           'get_confid_contour_defaults', 'set_window_redraw', 'set_jointplot',
           'get_histo_defaults', 'get_model_histo_defaults',
           'get_component_plot_defaults', 'get_component_histo_defaults',
           'vline', 'hline', 'get_scatter_plot_defaults',
           'get_cdf_plot_defaults', 'get_latex_for_string', 'name')


lgr = logging.getLogger(__name__)

# Identify the backend
name = 'dummy'

warning = lgr.warning
warning('Failed to import usable sherpa.plotting backends.' +
        ' Plotting routines will not be available')


def point(*args, **kwargs):
    """A do-nothing operation"""
    pass


clear_window = point
set_window_redraw = point
end = point
begin = point
exceptions = point

plot = point
histo = point
contour = point
set_subplot = point
set_jointplot = point

vline = point
hline = point


def get_split_plot_defaults():
    return get_keyword_defaults(set_subplot, 3)


def get_plot_defaults():
    return get_keyword_defaults(plot, 7)


def get_point_defaults():
    return get_keyword_defaults(point, 2)


def get_contour_defaults():
    return get_keyword_defaults(contour, 6)


def get_histo_defaults():
    return get_keyword_defaults(histo, 6)


def get_dummy_defaults():
    return {}


get_data_plot_defaults = get_dummy_defaults
get_model_plot_defaults = get_dummy_defaults
get_fit_plot_defaults = get_dummy_defaults
get_resid_plot_defaults = get_dummy_defaults
get_ratio_plot_defaults = get_dummy_defaults

get_data_contour_defaults = get_dummy_defaults
get_model_contour_defaults = get_dummy_defaults
get_fit_contour_defaults = get_dummy_defaults
get_resid_contour_defaults = get_dummy_defaults
get_ratio_contour_defaults = get_dummy_defaults

get_confid_point_defaults = get_dummy_defaults
get_confid_plot_defaults = get_dummy_defaults
get_confid_contour_defaults = get_dummy_defaults
get_model_histo_defaults = get_dummy_defaults
get_component_plot_defaults = get_dummy_defaults
get_component_histo_defaults = get_dummy_defaults
get_scatter_plot_defaults = get_dummy_defaults
get_cdf_plot_defaults = get_dummy_defaults


def get_latex_for_string(txt):
    """Convert to LaTeX form for the dummy back end.

    Parameters
    ----------
    txt : str
        The text component in LaTeX form (e.g. r'\alpha^2'). It
        should not contain any non-LaTeX content.

    Returns
    -------
    latex : str
        The input text (i.e. no change).

    """

    return txt


# HTML representation as tabular data
#
def as_html(data, fields):
    """Create HTML representation of a plot

    Parameters
    ----------
    data : Plot instance
        The plot object to display.
    fields : sequence of strings
        The fields of data to use.

    """

    # Would like a nicer way to set the summary label, but without
    # adding a per-class field for this it is safest just to use
    # the object name.

    meta = []
    for name in fields:
        # skip records which we don't know about. This indicates
        # an error in the calling code, but we don't want it to
        # stop the generation of the HTML.
        #
        try:
            val = getattr(data, name)
        except Exception as e:
            lgr.debug("Skipping field {}: {}".format(name, e))
            continue

        meta.append((name, val))

    ls = [formatting.html_section(meta, open_block=True,
                                  summary=type(data).__name__)]
    return formatting.html_from_sections(data, ls)


def as_html_histogram(plot):
    return as_html(plot,
                   ['xlo', 'xhi', 'y', 'title', 'xlabel', 'ylabel'])


def as_html_pdf(plot):
    return as_html(plot,
                   ['points', 'xlo', 'xhi', 'y', 'title', 'xlabel', 'ylabel'])


def as_html_cdf(plot):
    return as_html(plot,
                   ['points', 'x', 'y',
                    'median', 'lower', 'upper',
                    'title', 'xlabel', 'ylabel'])


def as_html_lr(plot):
    return as_html(plot,
                   ['ratios', 'lr', 'xlo', 'xhi', 'y',
                    'title', 'xlabel', 'ylabel'])


def as_html_data(plot):
    return as_html(plot,
                   ['x', 'xerr', 'y', 'yerr',
                    'title', 'xlabel', 'ylabel'])


def as_html_datacontour(plot):
    return as_html(plot,
                   ['x0', 'x1', 'y', 'levels',
                    'title', 'xlabel', 'ylabel'])


def as_html_model(plot):
    return as_html(plot,
                   ['x', 'xerr', 'y', 'yerr',
                    'title', 'xlabel', 'ylabel'])


def as_html_modelcontour(plot):
    return as_html(plot,
                   ['x0', 'x1', 'y', 'levels',
                    'title', 'xlabel', 'ylabel'])


def get_html(attr):
    if attr is None:
        return ''
    return attr._repr_html_()


def as_html_fit(plot):
    # Would like to do a better combination than this
    dplot = get_html(plot.dataplot)
    mplot = get_html(plot.modelplot)

    if dplot == '' and mplot == '':
        return None

    return dplot + mplot


def as_html_fitcontour(plot):
    # Would like to do a better combination than this
    dplot = get_html(plot.datacontour)
    mplot = get_html(plot.modelcontour)

    if dplot == '' and mplot == '':
        return None

    return dplot + mplot


def as_html_contour1d(plot):
    return as_html(plot,
                   ['x', 'y', 'min', 'max', 'nloop',
                    'delv', 'fac', 'log'])


def as_html_contour2d(plot):
    return as_html(plot,
                   ['parval0', 'parval1', 'sigma',
                    'x0', 'x1', 'y', 'levels',
                    'min', 'max', 'nloop',
                    'delv', 'fac', 'log'])
